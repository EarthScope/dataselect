/***************************************************************************
 * dataselect.c - miniSEED data selection.
 *
 * Opens one or more user specified files, applys filtering criteria
 * and outputs any matched data while time-ordering the data and
 * optionally pruning any overlap (at record or sample level) and
 * splitting records on day, hour or minute boundaries.
 *
 * In general critical error messages are prefixed with "ERROR:" and
 * the return code will be 1.  On successfull operation the return
 * code will be 0.
 *
 * Written by Chad Trabant, EarthScope Data Services.
 ***************************************************************************/

/***************************************************************************
 *
 * Data structures and operational overview
 *
 * The data map (using actual structure names):
 *
 * MS3TraceList
 *   |-MS3TraceID
 *   |   |-MS3TraceSeg
 *   |        |-RecordMap
 *   |            |-Record
 *   |            |-Record
 *   |            |-...
 *   |
 *   |-MS3TraceID
 *   |   |-MS3TraceSeg
 *   |        |-RecordMap
 *   |            |-Record
 *   |            |-Record
 *   |            |-...
 *   |
 *   |-...
 *
 * The program goes through the following stages:
 *
 * 1) Read all input files constructing a data map of contiguous trace
 * segments and the data records that comprise them.  This operation
 * is done by reading each file record-by-record and as each record is
 * read, search the MS3TraceList for an MS3TraceSeg that the record "fits"
 * with (same channel and adjacent in time).  If the record is found
 * to fit with a specific MS3TraceSeg its coverage will be added to the
 * MS3TraceSeg information otherwise a new MS3TraceSeg is created and
 * added to the MS3TraceList.
 *
 * Each MS3TraceSeg has an associated RecordMap which includes a list of
 * Records.  A Record structure is meta information about the record
 * coverage and location (file and offset).  The list of Records is
 * always in time order.
 *
 * When splitting on a time boundary and an input record crosses the
 * specified boundary, two Record structures will be created for the
 * two sides of the boundary and the new start and end time will be
 * set accordingly.  When writing this data out the data record will
 * effectively be split in two.
 *
 * There is no relationship between the location of input records in
 * specific files or offsets into files.  In other words, the program
 * will reconstruct the most contiguous, time-ordered data segments
 * possible from all the input records regardless of how they are
 * organized in the input files.  The resulting time-ordering of the
 * data records and contiguous segments is a characteristic of the
 * internal data structures and cannot be turned off.
 *
 * As each record is read the input data selection criteria are
 * applied.  For example, start/end time selections and regular
 * expression matching/rejecting selections are applied.
 *
 * 2) If data pruning (removing overlap data) has been selected the
 * data map will be processed to identify all overlapping data and to
 * mark individual Record structures either for complete removal or
 * for partial record trimming (when pruning at the sample level).
 * When a complete record is pruned from the ouput its record length
 * member will be set to 0 indicating that it is no longer
 * contributing to the segment, a special case understood by
 * downstream functions.  Note that no actual data records are changed
 * in this operation, modification of the data records occurs when the
 * data is written to the new output files.
 *
 * 3) Write all contributing data records in the data map out to the
 * output files.  After each record is read into memory it's associated
 * Record structure is checked to see if the record needs to be
 * trimmed due either to sample level pruning or time boundary
 * splitting.  Trimming a data record involves unpacking, sample
 * removal and repacking.  After trimming or if no trimming is
 * required the data record is written to the appropriate output file.
 * In this way only the minimal number of records needing modification
 * (trimming) are repacked.
 *
 ***************************************************************************/

/* _ISOC9X_SOURCE needed to get a declaration for llabs on some archs */
#define _ISOC9X_SOURCE

#define __STDC_FORMAT_MACROS
#include <errno.h>
#include <inttypes.h>
#include <math.h>
#include <regex.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/resource.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <time.h>

#include <libmseed.h>

#include "dsarchive.h"

#define VERSION "4.0DEV"
#define PACKAGE "dataselect"

/* Input/output file selection information containers */
typedef struct Filelink_s
{
  char *infilename_raw;   /* Input file name with potential annotation (byte range) */
  char *infilename;       /* Input file name without annotation (byte range) */
  FILE *infp;             /* Input file descriptor */
  uint64_t reordercount;  /* Number of records re-ordered */
  uint64_t recsplitcount; /* Number of records split */
  uint64_t recrmcount;    /* Number of records removed */
  uint64_t rectrimcount;  /* Number of records trimmed */
  nstime_t earliest;      /* Earliest data time in this file selection */
  nstime_t latest;        /* Latest data time in this file selection */
  uint64_t byteswritten;  /* Number of bytes written out */
  struct Filelink_s *next;
} Filelink;

/* Archive output structure definition containers */
typedef struct Archive_s
{
  DataStream datastream;
  struct Archive_s *next;
} Archive;

/* miniSEED record information structures */
typedef struct Record_s
{
  struct Filelink_s *flp;
  off_t offset;
  int reclen;
  nstime_t starttime;
  nstime_t endtime;
  uint8_t pubversion;
  nstime_t selectstart;
  nstime_t selectend;
  nstime_t newstart;
  nstime_t newend;
  struct Record_s *prev;
  struct Record_s *next;
} Record;

/* Record map, holds Record structures for a given MS3Trace */
typedef struct RecordMap_s
{
  long long int recordcnt;
  struct Record_s *first;
  struct Record_s *last;
} RecordMap;

/* Container for coverage entries used to prune data */
typedef struct Coverage_s
{
  nstime_t starttime;
  nstime_t endtime;
  uint8_t pubversion;
  double samprate;
  struct Coverage_s *next;
} Coverage;

/* Holder for data passed to the record writer */
typedef struct WriterData_s
{
  FILE *ofp;
  Record *rec;
  MS3TraceID *id;
  long *suffixp;
  int8_t *errflagp;
} WriterData;

static int readfiles (MS3TraceList **ppmstl);
static int processtraces (MS3TraceList *mstl);
static int writetraces (MS3TraceList *mstl);
static int trimrecord (Record *rec, char *recbuf, WriterData *writerdata);
static void writerecord (char *record, int reclen, void *handlerdata);

static int prunetraces (MS3TraceList *mstl);
static int findcoverage (MS3TraceList *mstl, MS3TraceID *targetid,
                         MS3TraceSeg *targetseg, Coverage **ppcoverage);
static int trimtrace (MS3TraceSeg *targetseg, const char *targetsourceid,
                      Coverage *coverage);
static int reconcile_tracetimes (MS3TraceList *mstl);

static void printmodsummary (flag nomods);
static void printtracemap (MS3TraceList *mstl);
static void printrecordmap (RecordMap *recmap, flag details);
static void printwritten (MS3TraceList *mstl);

static int findselectlimits (const MS3Selections *select, const char *sourceid,
                             nstime_t starttime, nstime_t endtime, Record *rec);
static int sortrecmap (RecordMap *recmap);
static int recordcmp (Record *rec1, Record *rec2);

static int processparam (int argcount, char **argvec);
static char *getoptval (int argcount, char **argvec, int argopt);
static int setofilelimit (int limit);
static int addfile (char *filename);
static int addlistfile (char *filename);
static int addarchive (const char *path, const char *layout);
static int readregexfile (char *regexfile, char **pppattern);
static void usage (int level);

static flag verbose = 0;
static flag basicsum = 0;             /* Controls printing of basic summary */
static flag bestversion = 1;          /* Use publication version to retain the "best" data when pruning */
static flag prunedata = 0;            /* Prune data: 'r= record level, 's' = sample level, 'e' = edges only */
static char restampqind = 0;          /* Re-stamp data record/quality indicator */
static flag modsummary = 0;           /* Print modification summary after all processing */
static nstime_t starttime = NSTUNSET; /* Limit to records containing or after starttime */
static nstime_t endtime = NSTUNSET;   /* Limit to records containing or before endtime */
static char splitboundary = 0;        /* Split records on day 'd', hour 'h' or minute 'm' boundaries */
static char splitreclen = 0;          /* Split output files on record length changes */
static double timetol = -1.0;         /* Time tolerance for continuous traces */
static double sampratetol = -1.0;     /* Sample rate tolerance for continuous traces */
static MS3Tolerance tolerance = {.time = NULL, .samprate = NULL};

/* Trivial callback functions for fixed time and sample rate tolerances */
double
timetol_callback (const MS3Record *msr)
{
  (void)msr;
  return timetol;
}
double
samprate_callback (const MS3Record *msr)
{
  (void)msr;
  return sampratetol;
}

static regex_t *match = NULL;  /* Compiled match regex */
static regex_t *reject = NULL; /* Compiled reject regex */

static char *outputfile = NULL;  /* Single output file */
static flag outputmode = 0;      /* Mode for single output file: 0=overwrite, 1=append */
static Archive *archiveroot = 0; /* Output file structures */

static char recordbuf[16384]; /* Global record buffer */

static Filelink *filelist = NULL;        /* List of input files */
static Filelink *filelisttail = NULL;    /* Tail of list of input files */
static MS3Selections *selections = NULL; /* List of data selections */

static char *writtenfile = NULL;       /* File to write summary of output records */
static char *writtenprefix = NULL;     /* Prefix for summary of output records */
static MS3TraceList *writtentl = NULL; /* TraceList of output records */

int
main (int argc, char **argv)
{
  MS3TraceList *mstl = NULL;
  char *leapsecondfile = NULL;

  /* Set default error message prefix */
  ms_loginit (NULL, NULL, NULL, "ERROR: ");

  /* Process input parameters */
  if (processparam (argc, argv) < 0)
    return 1;

  /* Read leap second list file if env. var. LIBMSEED_LEAPSECOND_FILE is set */
  if ((leapsecondfile = getenv ("LIBMSEED_LEAPSECOND_FILE")))
  {
    if (strcmp (leapsecondfile, "NONE"))
      ms_readleapsecondfile (leapsecondfile);
  }
  else if (verbose >= 1)
  {
    ms_log (1, "Warning: No leap second file specified with LIBMSEED_LEAPSECOND_FILE\n");
    ms_log (1, "  This is highly recommended, see man page for details.\n");
  }

  /* Data stream archiving maximum concurrent open files */
  if (archiveroot)
    ds_maxopenfiles = 50;

  /* Init written MS3TraceList */
  if (writtenfile)
    if ((writtentl = mstl3_init (writtentl)) == NULL)
      return 1;

  if (verbose > 2)
    ms_log (1, "Processing input files\n");

  /* Read and process all files specified on the command line */
  if (readfiles (&mstl))
    return 1;

  /* Processes traces */
  if (mstl->numtraceids > 0 && processtraces (mstl))
    return 1;

  if (modsummary)
    printmodsummary (verbose);

  if (writtenfile)
  {
    printwritten (writtentl);
    mstl3_free (&writtentl, 1);
  }

  /* The main MS3TraceList (mstl) is not freed on purpose: the structure has a
   * potentially huge number of sub-structures (Records in the RecordMap of
   * each ID) which would take a long time to iterate through.  This would be a
   * waste of time given that the program is now done.
   * Mind that this may show up as a memory leak for some profilers. */

  return 0;
} /* End of main() */

/***************************************************************************
 * readfiles:
 *
 * Read input files specified as a Filelink list and populate an
 * MS3TraceList and record maps for each trace.  All input files are
 * renamed with a ".orig" suffix before being read.
 *
 * Returns 0 on success and -1 otherwise.
 ***************************************************************************/
static int
readfiles (MS3TraceList **ppmstl)
{
  MS3FileParam *msfp = NULL;
  Filelink *flp;
  MS3Record *msr = NULL;
  MS3TraceSeg *seg = NULL;

  int totalrecs = 0;
  int totalsamps = 0;
  int totalfiles = 0;

  const MS3Selections *matchsp = NULL;
  const MS3SelectTime *matchstp = NULL;

  RecordMap *recmap = NULL;
  Record *rec = NULL;

  RecordMap newrecmap;
  Record *newrec = NULL;

  nstime_t recstarttime;
  nstime_t recendtime;

  char stime[30] = {0};

  uint32_t flags = 0;
  int retcode;
  flag whence;

  if (!ppmstl)
    return -1;

  /* Initialize MS3TraceList */
  *ppmstl = mstl3_init (*ppmstl);

  if (!*ppmstl)
  {
    ms_log (2, "readfiles(): cannot (re)initialize MS3TraceList\n");
    return -1;
  }

  /* Read all input files and construct continuous traces, using the
   * libmseed MS3TraceList.  For each trace maintain a list of each
   * data record that contributed to the trace, implemented as a
   * RecordMap struct (MS3TraceSeg->prvtptr) where a linked list of
   * Record structs is maintained.  The records are always listed in
   * time order.
   */

  flp = filelist;

  while (flp)
  {
    if (verbose)
    {
      if (strcmp (flp->infilename, flp->infilename_raw) == 0)
        ms_log (1, "Reading: %s\n", flp->infilename);
      else
        ms_log (1, "Reading: %s (specified as %s)\n", flp->infilename, flp->infilename_raw);
    }

    /* Set bit flags to validate CRC and extrace start-stop range from file names */
    flags |= MSF_VALIDATECRC;
    flags |= MSF_PNAMERANGE;

    /* Loop over the input file */
    while ((retcode = ms3_readmsr_selection (&msfp, &msr, flp->infilename_raw, flags,
                                             selections, verbose - 2)) == MS_NOERROR)
    {
      recstarttime = msr->starttime;
      recendtime = msr3_endtime (msr);

      /* Check if record matches start time criteria: starts after or contains starttime */
      if ((starttime != NSTUNSET) && (recstarttime < starttime && !(recstarttime <= starttime && recendtime >= starttime)))
      {
        if (verbose >= 3)
        {
          ms_nstime2timestr (recstarttime, stime, ISOMONTHDAY_Z, NANO_MICRO);
          ms_log (1, "Skipping (starttime) %s, %s\n", msr->sid, stime);
        }
        continue;
      }

      /* Check if record matches end time criteria: ends after or contains endtime */
      if ((endtime != NSTUNSET) && (recendtime > endtime && !(recstarttime <= endtime && recendtime >= endtime)))
      {
        if (verbose >= 3)
        {
          ms_nstime2timestr (recstarttime, stime, ISOMONTHDAY_Z, NANO_MICRO);
          ms_log (1, "Skipping (endtime) %s, %s\n", msr->sid, stime);
        }
        continue;
      }

      /* Check if record is matched by the match regex */
      if (match)
      {
        if (regexec (match, msr->sid, 0, 0, 0) != 0)
        {
          if (verbose >= 3)
          {
            ms_nstime2timestr (recstarttime, stime, ISOMONTHDAY_Z, NANO_MICRO);
            ms_log (1, "Skipping (match) %s, %s\n", msr->sid, stime);
          }
          continue;
        }
      }

      /* Check if record is rejected by the reject regex */
      if (reject)
      {
        if (regexec (reject, msr->sid, 0, 0, 0) == 0)
        {
          if (verbose >= 3)
          {
            ms_nstime2timestr (recstarttime, stime, ISOMONTHDAY_Z, NANO_MICRO);
            ms_log (1, "Skipping (reject) %s, %s\n", msr->sid, stime);
          }
          continue;
        }
      }

      /* Check if record is matched by selection */
      if (selections)
      {
        if (!(matchsp = ms3_matchselect (selections, msr->sid, recstarttime,
                                         recendtime, msr->pubversion, &matchstp)))
        {
          if (verbose >= 3)
          {
            ms_nstime2timestr (recstarttime, stime, ISOMONTHDAY_Z, NANO_MICRO);
            ms_log (1, "Skipping (selection) %s, %s\n", msr->sid, stime);
          }
          continue;
        }
      }

      if (verbose > 2)
        msr3_print (msr, verbose - 3);

      /* Add record to the MS3TraceList */
      if (!(seg = mstl3_addmsr (*ppmstl, msr, bestversion, 0, 0, &tolerance)))
      {
        ms_nstime2timestr (recstarttime, stime, ISOMONTHDAY_Z, NANO_MICRO);
        ms_log (2, "Cannot add record to trace list, %s, %s\n", msr->sid, stime);
        continue;
      }

      /* Determine where the record fit with this MS3TraceSeg
       * whence:
       *   0 = New MS3Trace
       *   1 = End of MS3Trace
       *   2 = Beginning of MS3Trace
       */
      whence = 0;
      if (seg->prvtptr)
      {
        if (seg->endtime == recendtime)
          whence = 1;
        else if (seg->starttime == recstarttime)
          whence = 2;
        else if (recendtime == recstarttime)
        {
          /* Determine best fit for records with no span */
          if (llabs (recstarttime - seg->endtime) < llabs (recstarttime - seg->starttime))
            whence = 1;
          else
            whence = 2;
        }
        else
        {
          ms_log (2, "Cannot determine where record fit relative to trace segment\n");
          msr3_print (msr, 1);
          continue;
        }
      }

      /* Create and populate new Record structure */
      if (!(rec = (Record *)malloc (sizeof (Record))))
      {
        ms_log (2, "Cannot allocate memory for Record entry\n");
        return -1;
      }

      rec->flp = flp;
      rec->offset = (msfp->streampos - msr->reclen);
      rec->reclen = msr->reclen;
      rec->starttime = recstarttime;
      rec->endtime = recendtime;
      rec->pubversion = msr->pubversion;
      rec->selectstart = NSTUNSET;
      rec->selectend = NSTUNSET;
      rec->newstart = NSTUNSET;
      rec->newend = NSTUNSET;
      rec->prev = NULL;
      rec->next = NULL;

      /* Populate a new record map */
      newrecmap.recordcnt = 1;
      newrecmap.first = rec;
      newrecmap.last = rec;

      /* If record is not completely selected search for joint selection limits */
      if (matchstp && !(matchstp->starttime <= recstarttime && matchstp->endtime >= recendtime))
      {
        if (findselectlimits (matchsp, msr->sid, recstarttime, recendtime, rec))
        {
          ms_log (2, "Problem in findselectlimits(), please report\n");
        }
      }

      /* If pruning at the sample level trim right at the start/end times */
      if (prunedata == 's' || prunedata == 'e')
      {
        nstime_t seltime;

        /* Determine strictest start time (selection time or global start time) */
        if (starttime != NSTUNSET && rec->selectstart != NSTUNSET)
          seltime = (starttime > rec->selectstart) ? starttime : rec->selectstart;
        else if (rec->selectstart != NSTUNSET)
          seltime = rec->selectstart;
        else
          seltime = starttime;

        /* If the Record crosses the start time */
        if (seltime != NSTUNSET && (seltime > recstarttime) && (seltime <= recendtime))
        {
          rec->newstart = seltime;
        }

        /* Determine strictest end time (selection time or global end time) */
        if (endtime != NSTUNSET && rec->selectend != NSTUNSET)
          seltime = (endtime < rec->selectend) ? endtime : rec->selectend;
        else if (rec->selectend != NSTUNSET)
          seltime = rec->selectend;
        else
          seltime = endtime;

        /* If the Record crosses the end time */
        if (seltime != NSTUNSET && (seltime >= recstarttime) && (seltime < recendtime))
        {
          rec->newend = seltime;
        }
      }

      /* Create extra Record structures if splitting on a time boundary */
      if (splitboundary)
      {
        uint16_t year;
        uint16_t yday;
        uint8_t hour;
        uint8_t min;
        uint8_t sec;
        uint32_t nsec;
        nstime_t boundary = NSTUNSET;
        nstime_t effstarttime;

        for (;;)
        {
          effstarttime = (rec->newstart != NSTUNSET) ? rec->newstart : rec->starttime;
          ms_nstime2time (effstarttime, &year, &yday, &hour, &min, &sec, &nsec);

          /* Determine next split boundary */
          if (splitboundary == 'd') /* Days */
          {
            yday += 1;
            hour = min = sec = nsec = 0;
            boundary = ms_time2nstime (year, yday, hour, min, sec, nsec);
          }
          else if (splitboundary == 'h') /* Hours */
          {
            hour += 1;
            min = sec = nsec = 0;
            boundary = ms_time2nstime (year, yday, hour, min, sec, nsec);
          }
          else if (splitboundary == 'm') /* Minutes */
          {
            min += 1;
            sec = nsec = 0;
            boundary = ms_time2nstime (year, yday, hour, min, sec, nsec);
          }
          else
          {
            ms_log (2, "Split boundary code unrecognized: '%c'\n", splitboundary);
            break;
          }

          /* If end time is beyond the boundary create a new Record */
          if (rec->endtime >= boundary)
          {
            if (!(newrec = (Record *)malloc (sizeof (Record))))
            {
              ms_log (2, "Cannot allocate memory for Record entry");
              return -1;
            }

            memcpy (newrec, rec, sizeof (Record));

            /* Set current Record and next Record new boundary times */
            rec->newend = boundary - 1;
            newrec->newstart = boundary;

            /* Update new record map */
            newrecmap.recordcnt++;
            newrecmap.last = newrec;

            /* Insert the new Record in chain and set as current */
            rec->next = newrec;
            newrec->prev = rec;
            rec = newrec;

            flp->recsplitcount++;
          }
          /* Otherwise we are done */
          else
          {
            break;
          }
        }
      } /* Done splitting on time boundary */

      /* Add the new Record(s) to the RecordMap associated with the MS3TraceSeg */

      /* Add new Record(s) to end of the RecordMap */
      if (whence == 1)
      {
        recmap = (RecordMap *)seg->prvtptr;

        recmap->last->next = newrecmap.first;
        newrecmap.first->prev = recmap->last;

        recmap->last = newrecmap.last;

        recmap->recordcnt += newrecmap.recordcnt;
      }
      /* Add new Record(s) to beginning of the RecordMap */
      else if (whence == 2)
      {
        recmap = (RecordMap *)seg->prvtptr;

        recmap->first->prev = newrecmap.last;
        newrecmap.last->next = recmap->first;

        recmap->first = newrecmap.first;

        recmap->recordcnt += newrecmap.recordcnt;

        /* Increment reordered count */
        flp->reordercount++;
      }
      /* First Record(s) for this MS3TraceSeg, allocate RecordMap */
      else
      {
        if (seg->prvtptr)
          ms_log (2, "Supposedly first record, but RecordMap not empty, report this\n");

        if (!(recmap = (RecordMap *)malloc (sizeof (RecordMap))))
        {
          ms_log (2, "Cannot allocate memory for new RecordMap\n");
          return -1;
        }
        recmap->first = newrecmap.first;
        recmap->last = newrecmap.last;
        recmap->recordcnt = newrecmap.recordcnt;

        seg->prvtptr = recmap;
      }

      totalrecs++;
      totalsamps += msr->samplecnt;
    } /* End of looping through records in file */

    /* Critical error if file was not read properly */
    if (retcode != MS_ENDOFFILE)
    {
      ms_log (2, "Cannot read %s: %s\n", flp->infilename, ms_errorstr (retcode));
      ms3_readmsr_selection (&msfp, &msr, NULL, flags, NULL, 0);
      return -1;
    }

    /* Make sure everything is cleaned up */
    ms3_readmsr_selection (&msfp, &msr, NULL, flags, NULL, 0);

    totalfiles++;
    flp = flp->next;
  } /* End of looping over file list */

  /* Increase open file limit if necessary, in general we need the
   * filecount + ds_maxopenfiles and some wiggle room. */
  setofilelimit (totalfiles + ds_maxopenfiles + 20);

  if (basicsum)
    ms_log (0, "Files: %d, Records: %d, Samples: %d\n", totalfiles, totalrecs, totalsamps);

  return 0;
} /* End of readfiles() */

/***************************************************************************
 * processtraces:
 *
 * Process all data represented by the MS3TraceList by first pruning
 * them and then writing out the remaining data.
 *
 * Returns 0 on success and -1 on error.
 ***************************************************************************/
static int
processtraces (MS3TraceList *mstl)
{
  if (verbose > 2)
    printtracemap (mstl);

  /* Prune data */
  if (prunedata)
  {
    /* Prune overlaps */
    if (prunedata == 'r' || prunedata == 's')
      if (prunetraces (mstl))
        return -1;

    /* Reconcile MS3TraceID times with associated record maps */
    if (reconcile_tracetimes (mstl))
      return -1;
  }

  /* Write all MS3TraceSeg associated records to output file(s) */
  if (writetraces (mstl))
    return -1;

  return 0;
} /* End of processtraces() */

/***************************************************************************
 * writetraces():
 *
 * Write all MS3TraceSeg associated records to output file(s).  If an
 * output file is specified all records will be written to it,
 * otherwise records will be written to the original files and
 * (optionally) backups of the original files will remain.
 *
 * This routine will also call trimrecord() to trim a record when data
 * suturing is requested.  Record trimming is triggered when
 * Record.newstart or Record.newend are set for any output records.
 *
 * TODO allow setting quality somewhere before here?  Not sure how
 * The quality flag is optionally set for all output records.
 *
 * Returns 0 on success and 1 on error.
 ***************************************************************************/
static int
writetraces (MS3TraceList *mstl)
{
  static uint64_t totalrecsout = 0;
  static uint64_t totalbytesout = 0;
  char *recordptr = NULL;
  char *wb = "wb";
  char *ab = "ab";
  char *mode;
  int8_t errflag = 0;
  long suffix = 0;
  int rv;

  MS3TraceID *id;
  MS3TraceID *groupid;
  MS3TraceSeg *seg;

  RecordMap *groupmap = NULL;
  RecordMap *recmap = NULL;
  Record *rec;
  Record *recnext;
  Filelink *flp;

  FILE *ofp = NULL;
  WriterData writerdata;

  writerdata.errflagp = &errflag;
  writerdata.suffixp = &suffix;

  if (!mstl)
    return 1;

  if (!mstl->traces.next[0])
    return 1;

  if (verbose)
    ms_log (1, "Writing output data\n");

  /* Open the output file if specified */
  if (outputfile)
  {
    /* Decide if we are appending or overwriting */
    mode = (totalbytesout || outputmode) ? ab : wb;

    if (strcmp (outputfile, "-") == 0)
    {
      ofp = stdout;
    }
    else if ((ofp = fopen (outputfile, mode)) == NULL)
    {
      ms_log (2, "Cannot open output file: %s (%s)\n",
              outputfile, strerror (errno));
      return 1;
    }
  }

  /* Re-link records into write lists, from per-segment lists to per ID lists.
   * This allows (later) sorting of data records as logical groups regardless
   * from which segment the record was originally associated. */
  id = mstl->traces.next[0];
  groupid = id;
  while (id)
  {
    /* Determine if this is a new group */
    if (groupid != id)
    {
      /* If data was pruned, group by source ID.
       * When pruning, data is grouped by source ID and publication version.
       * This groups the write lists by source ID without verison considered. */
      if (prunedata)
      {
        if (strcmp (id->sid, groupid->sid))
        {
          groupid = id;
        }
      }
      else
      {
        groupid = id;
      }
    }

    seg = id->first;
    if (seg)
    {
      /* Allocate group ID RecordMap if needed */
      if (!groupid->prvtptr)
      {
        if (!(id->prvtptr = calloc (1, sizeof (RecordMap))))
        {
          ms_log (2, "writetraces(): Cannot allocate memory\n");
          return 1;
        }
        groupmap = (RecordMap *)id->prvtptr;
      }

      while (seg)
      {
        recmap = (RecordMap *)seg->prvtptr;
        rec = recmap->first;
        while (rec && !errflag)
        {
          recnext = rec->next;

          /* Free and skip marked (pre-identified as non-contributing) records */
          if (rec->reclen == 0)
          {
            free (rec);
            rec = recnext;
            continue;
          }

          /* Add record to group ID write list */
          if (!groupmap->first)
          {
            groupmap->first = rec;
            rec->prev = NULL;
            rec->next = NULL;
            groupmap->last = rec;
            groupmap->recordcnt = 1;
          }
          else
          {
            groupmap->last->next = rec;
            rec->prev = groupmap->last;
            rec->next = NULL;
            groupmap->last = rec;
            groupmap->recordcnt++;
          }

          rec = recnext;
        } /* Done looping through Records in the RecordMap */

        /* Free segment RecordMap and cauterize */
        if (seg->prvtptr)
          free (seg->prvtptr);
        seg->prvtptr = NULL;

        seg = seg->next;
      } /* Done looping through MS3TraceSegs in the MS3TraceID */
    }

    id = id->next[0];
  } /* Done looping through MS3TraceIDs in the MS3TraceList */

  /* Loop through MS3TraceList and write records in write lists */
  id = mstl->traces.next[0];
  while (id && errflag != 1)
  {
    /* Skip when no write list is present */
    if (id->prvtptr == NULL)
    {
      id = id->next[0];
      continue;
    }

    /* Reset error flag for continuation errors */
    if (errflag == 2)
      errflag = 0;

    /* Set suffix if splitting on record length changes */
    if (splitreclen)
      suffix = 1;

    recmap = (RecordMap *)id->prvtptr;

    /* Sort record list if overlaps have been pruned, if the data has not been
     * pruned it is already in time order. */
    if (prunedata == 'r' || prunedata == 's')
    {
      sortrecmap (recmap);
    }

    /* Loop through each Record in the write list RecordMap.
     * After records are read from the input files, perform any
     * pre-identified pruning before writing data back out */
    rec = recmap->first;
    while (rec && !errflag)
    {
      /* Read the record from the input file */

      if (rec->reclen > sizeof (recordbuf))
      {
        ms_log (2, "Record length (%d bytes) larger than buffer (%llu bytes)\n",
                rec->reclen, (long long unsigned int)sizeof (recordbuf));
        errflag = 1;
        break;
      }

      /* Open file for reading if not already done */
      if (!rec->flp->infp)
        if (!(rec->flp->infp = fopen (rec->flp->infilename, "rb")))
        {
          ms_log (2, "Cannot open '%s' for reading: %s\n",
                  rec->flp->infilename, strerror (errno));
          errflag = 1;
          break;
        }

      /* Seek to record offset */
      if (lmp_fseek64 (rec->flp->infp, rec->offset, SEEK_SET) == -1)
      {
        ms_log (2, "Cannot seek in '%s': %s\n",
                rec->flp->infilename, strerror (errno));
        errflag = 1;
        break;
      }

      /* Read record into buffer */
      if (fread (recordbuf, rec->reclen, 1, rec->flp->infp) != 1)
      {
        ms_log (2, "Cannot read %d bytes at offset %llu from '%s'\n",
                rec->reclen, (long long unsigned)rec->offset,
                rec->flp->infilename);
        errflag = 1;
        break;
      }

      recordptr = recordbuf;

      /* Re-stamp quality indicator if specified */
      // TODO fix for two versions
      if (restampqind)
      {
        if (verbose > 3)
          ms_log (1, "Re-stamping data quality indicator to '%c'\n", restampqind);

        *(recordptr + 6) = restampqind;
      }

      /* Setup writer data */
      writerdata.ofp = ofp;
      writerdata.rec = rec;
      writerdata.id = id;

      /* Write out the data, either the record needs to be trimmed (and will be
       * send to the record writer) or we send it directly to the record writer.
       */

      /* Trim data from the record if new start or end times are specifed */
      if (rec->newstart != NSTUNSET || rec->newend != NSTUNSET)
      {
        rv = trimrecord (rec, recordptr, &writerdata);

        if (rv == -1)
        {
          rec = rec->next;
          continue;
        }
        if (rv == -2)
        {
          ms_log (2, "Cannot unpack miniSEED from byte offset %lld in %s\n",
                  rec->offset, rec->flp->infilename);
          ms_log (2, "  Expecting %s, skipping the rest of this channel\n", id->sid);
          errflag = 2;
          break;
        }
      }
      else
      {
        writerecord (recordptr, rec->reclen, &writerdata);
      }

      if (errflag)
        break;

      /* Update file entry time stamps and counts */
      if (!rec->flp->earliest || (rec->flp->earliest > rec->starttime))
      {
        rec->flp->earliest = rec->starttime;
      }
      if (!rec->flp->latest || (rec->flp->latest < rec->endtime))
      {
        rec->flp->latest = rec->endtime;
      }

      rec->flp->byteswritten += rec->reclen;

      totalrecsout++;
      totalbytesout += rec->reclen;

      /* Increment suffix if splitting and record length changes */
      if (splitreclen && rec->next && rec->next->reclen != rec->reclen)
      {
        suffix++;
      }

      rec = rec->next;
    } /* Done looping through Records in the RecordMap */

    id = id->next[0];
  } /* Done looping through MS3TraceIDs the MS3TraceList */

  /* Close all open input & output files and remove backups if requested */
  flp = filelist;
  while (flp)
  {
    if (!ofp && verbose)
    {
      ms_log (1, "Wrote %" PRId64 " bytes from file %s\n",
              flp->byteswritten, flp->infilename);
    }

    if (flp->infp)
    {
      fclose (flp->infp);
      flp->infp = NULL;
    }

    flp = flp->next;
  }

  /* Close output file if used */
  if (ofp)
  {
    fclose (ofp);
    ofp = NULL;
  }

  if (verbose)
  {
    ms_log (1, "Wrote %llu bytes of %llu records to output file(s)\n",
            totalbytesout, totalrecsout);
  }

  return (errflag) ? 1 : 0;
} /* End of writetraces() */

/***************************************************************************
 * trimrecord():
 *
 * Unpack a data record and trim samples, either from the beginning or
 * the end, to fit the Record.newstart and/or Record.newend boundary
 * times and pack the record.  Record.newstart and Record.newend are
 * treated as arbitrary boundaries, not as explicit new start/end
 * times, this routine calculates which samples fit within the new
 * boundaries.
 *
 * Return 0 on success, -1 on failure or skip and -2 on unpacking errors.
 ***************************************************************************/
static int
trimrecord (Record *rec, char *recordbuf, WriterData *writerdata)
{
  MS3Record *msr = NULL;
  nstime_t nsdelta;
  nstime_t ostarttime = rec->starttime;

  char stime[30] = {0};
  char etime[30] = {0};

  int trimsamples;
  uint8_t samplesize;
  char sampletype;
  int64_t packedsamples;
  int packedrecords;
  int retcode;

  if (!rec || !recordbuf)
    return -1;

  /* Sanity check for new start/end times */
  if ((rec->newstart != NSTUNSET && rec->newend != NSTUNSET && rec->newstart > rec->newend) ||
      (rec->newstart != NSTUNSET && (rec->newstart < rec->starttime || rec->newstart > rec->endtime)) ||
      (rec->newend != NSTUNSET && (rec->newend > rec->endtime || rec->newend < rec->starttime)))
  {
    ms_log (2, "Problem with new start/end record bound times.\n");
    // TODO, add SourceID to log message if restructured and information is available.
    ms_log (2, "  Original record %s from %s (byte offset: %llu)\n",
            "SourceID", rec->flp->infilename, (unsigned long long)rec->offset);
    ms_nstime2timestr (rec->starttime, stime, ISOMONTHDAY_Z, NANO_MICRO);
    ms_nstime2timestr (rec->endtime, etime, ISOMONTHDAY_Z, NANO_MICRO);
    ms_log (2, "       Start: %s       End: %s\n", stime, etime);
    if (rec->newstart == NSTUNSET)
      strcpy (stime, "NONE");
    else
      ms_nstime2timestr (rec->newstart, stime, ISOMONTHDAY_Z, NANO_MICRO);
    if (rec->newend == NSTUNSET)
      strcpy (etime, "NONE");
    else
      ms_nstime2timestr (rec->newend, etime, ISOMONTHDAY_Z, NANO_MICRO);
    ms_log (2, " Start bound: %-24s End bound: %-24s\n", stime, etime);

    return -1;
  }

  /* Parse data record header without decoding samples */
  if ((retcode = msr3_parse (recordbuf, rec->reclen, &msr, 0, 0)))
  {
    ms_log (2, "Cannot unpack miniSEED record: %s\n", ms_errorstr (retcode));
    return -2;
  }

  if (ms_encoding_sizetype (msr->encoding, &samplesize, &sampletype))
  {
    ms_log (2, "Cannot determine sample size and type for encoding %d\n", msr->encoding);
    msr3_free (&msr);
    return -1;
  }

  /* Check for supported samples types, can only trim what can be packed */
  if (sampletype != 'i' && sampletype != 'f' && sampletype != 'd')
  {
    if (verbose)
    {
      ms_nstime2timestr (rec->starttime, stime, ISOMONTHDAY_Z, NANO_MICRO);
      ms_log (1, "Skipping trim of %s (%s), unsupported encoding (%d: %s)\n",
              msr->sid, stime, msr->encoding, ms_encodingstr (msr->encoding));
    }

    msr3_free (&msr);
    return 0;
  }

  /* Decode data samples */
  if ((retcode = msr3_unpack_data (msr, 0)) < 0)
  {
    ms_log (2, "Cannot unpack miniSEED record: %s\n", ms_errorstr (retcode));
    msr3_free (&msr);
    return -2;
  }

  if (verbose > 1)
  {
    ms_log (1, "Triming record: %s (%u)\n", msr->sid, msr->pubversion);
    ms_nstime2timestr (rec->starttime, stime, ISOMONTHDAY_Z, NANO_MICRO);
    ms_nstime2timestr (rec->endtime, etime, ISOMONTHDAY_Z, NANO_MICRO);
    ms_log (1, "       Start: %s        End: %s\n", stime, etime);
    if (rec->newstart == NSTUNSET)
      strcpy (stime, "NONE");
    else
      ms_nstime2timestr (rec->newstart, stime, ISOMONTHDAY_Z, NANO_MICRO);
    if (rec->newend == NSTUNSET)
      strcpy (etime, "NONE");
    else
      ms_nstime2timestr (rec->newend, etime, ISOMONTHDAY_Z, NANO_MICRO);
    ms_log (1, " Start bound: %-24s  End bound: %-24s\n", stime, etime);
  }

  /* Determine sample period in nanosecond time ticks */
  nsdelta = (msr->samprate) ? (nstime_t)(NSTMODULUS / msr->samprate) : 0;

  /* Remove samples from the beginning of the record */
  if (rec->newstart != NSTUNSET && nsdelta)
  {
    nstime_t newstarttime;

    /* Determine new start time and the number of samples to trim */
    trimsamples = 0;
    newstarttime = rec->starttime;

    while (newstarttime < rec->newstart && trimsamples < msr->samplecnt)
    {
      newstarttime += nsdelta;
      trimsamples++;
    }

    if (trimsamples >= msr->samplecnt)
    {
      if (verbose > 1)
        ms_log (1, "All samples would be trimmed from record, skipping\n");

      msr3_free (&msr);
      return -1;
    }

    if (verbose > 2)
    {
      ms_nstime2timestr (newstarttime, stime, ISOMONTHDAY_Z, NANO_MICRO);
      ms_log (1, "Removing %d samples from the start, new start time: %s\n", trimsamples, stime);
    }

    memmove (msr->datasamples,
             (char *)msr->datasamples + (samplesize * trimsamples),
             samplesize * (msr->numsamples - trimsamples));

    msr->numsamples -= trimsamples;
    msr->samplecnt -= trimsamples;
    msr->starttime = newstarttime;
    rec->starttime = newstarttime;
  }

  /* Remove samples from the end of the record */
  if (rec->newend != NSTUNSET && nsdelta)
  {
    nstime_t newendtime;

    /* Determine new end time and the number of samples to trim */
    trimsamples = 0;
    newendtime = rec->endtime;

    while (newendtime > rec->newend && trimsamples < msr->samplecnt)
    {
      newendtime -= nsdelta;
      trimsamples++;
    }

    if (trimsamples >= msr->samplecnt)
    {
      if (verbose > 1)
        ms_log (1, "All samples would be trimmed from record, skipping\n");

      msr3_free (&msr);
      return -1;
    }

    if (verbose > 2)
    {
      ms_nstime2timestr (newendtime, etime, ISOMONTHDAY_Z, NANO_MICRO);
      ms_log (1, "Removing %d samples from the end, new end time: %s\n", trimsamples, etime);
    }

    msr->numsamples -= trimsamples;
    msr->samplecnt -= trimsamples;
    rec->endtime = newendtime;
  }

  /* Pack the data record into the global record buffer used by writetraces() */
  packedrecords = msr3_pack (msr, &writerecord, writerdata,
                             &packedsamples, MSF_FLUSHDATA, verbose - 1);

  // TODO fix below logic, wasted time when packed records > 1
  if (packedrecords != 1)
  {
    ms_nstime2timestr (ostarttime, stime, ISOMONTHDAY_Z, NANO_MICRO);

    if (packedrecords <= 0)
    {
      ms_log (2, "%s(): Cannot pack miniSEED record for %s %s\n",
              __func__, msr->sid, stime);
      msr3_free (&msr);
      return -2;
    }
  }

  /* Clean up MS3Record */
  msr3_free (&msr);

  return 0;
} /* End of trimrecord() */

/***************************************************************************
 * writerecord():
 *
 * Used by writetraces() directly, and trimrecord() when called, to save
 * repacked miniSEED to global record buffer.
 ***************************************************************************/
static void
writerecord (char *record, int reclen, void *handlerdata)
{
  WriterData *writerdata = handlerdata;
  MS3Record *msr = NULL;
  Archive *arch;
  int retcode;

  if (!record || reclen <= 0 || !handlerdata)
    return;

  /* Write to a single output file if specified */
  if (writerdata->ofp)
  {
    if (fwrite (record, reclen, 1, writerdata->ofp) != 1)
    {
      ms_log (2, "Cannot write to '%s'\n", outputfile);
      *writerdata->errflagp = 1;
    }
  }

  /* Write to Archive(s) if specified and/or add to written list */
  if (archiveroot || writtenfile)
  {
    /* Parse data record header without decoding samples */
    if ((retcode = msr3_parse (record, reclen, &msr, 0, 0)))
    {
      ms_log (2, "Cannot unpack miniSEED record: %s\n", ms_errorstr (retcode));
      ms_log (2, "  From byte offset %lld in %s, file changed?\n",
              writerdata->rec->offset, writerdata->rec->flp->infilename);
      ms_log (2, "  Expecting %s, skipping the rest of this channel\n",
              writerdata->id->sid);
      *writerdata->errflagp = 2;
    }

    if (archiveroot)
    {
      arch = archiveroot;
      while (arch)
      {
        if (ds_streamproc (&arch->datastream, msr, *writerdata->suffixp, verbose - 1))
        {
          *writerdata->errflagp = 1;
        }

        arch = arch->next;
      }
    }

    if (writtenfile)
    {
      MS3TraceSeg *seg;

      if ((seg = mstl3_addmsr (writtentl, msr, 0, 0, 0, NULL)) == NULL)
      {
        ms_log (2, "Error adding MS3Record to MS3TraceList, bah humbug.\n");
      }
      else
      {
        if (!seg->prvtptr)
        {
          if ((seg->prvtptr = malloc (sizeof (int64_t))) == NULL)
          {
            ms_log (2, "Error allocating memory for written count, bah humbug.\n");
            *writerdata->errflagp = 1;
          }
          else
          {
            *((int64_t *)seg->prvtptr) = 0;
          }
        }

        *((int64_t *)seg->prvtptr) += msr->reclen;
      }
    }

    msr3_free (&msr);
  }
} /* End of writerecord() */

/***************************************************************************
 * prunetraces():
 *
 * Prune all redundant data from the RecordMap entries associated with
 * the specified MS3TraceSegs.
 *
 * For each MS3TraceSeg determine the coverage of the RecordMap associated
 * with each overlapping, higher-priority MS3TraceSeg using findcoverage().
 * If some higher-priority overlap was determined to exist modify the
 * RecordMap of the MS3TraceSeg in question to mark the overlapping data
 * using trimtrace().
 *
 * Return 0 on success and -1 on failure.
 ***************************************************************************/
static int
prunetraces (MS3TraceList *mstl)
{
  MS3TraceID *id = NULL;
  MS3TraceSeg *seg = NULL;
  Coverage *coverage = NULL;
  int retval;

  if (!mstl)
    return -1;

  if (!mstl->traces.next[0])
    return -1;

  if (verbose)
    ms_log (1, "Pruning trace data\n");

  /* For each MS3TraceSeg determine the coverage of the overlapping
     Records from the other traces with a higher priority and prune
     the overlap. */
  id = mstl->traces.next[0];
  while (id)
  {
    seg = id->first;
    while (seg)
    {
      /* Determine overlapping trace coverage */
      retval = findcoverage (mstl, id, seg, &coverage);

      if (retval)
      {
        ms_log (2, "cannot findcoverage()\n");
        return -1;
      }
      else if (coverage)
      {
        if (trimtrace (seg, id->sid, coverage) < 0)
        {
          ms_log (2, "cannot trimtraces()\n");
          return -1;
        }
      }

      /* Free the coverage */
      while (coverage)
      {
        Coverage *next = coverage->next;
        free (coverage);
        coverage = next;
      }

      seg = seg->next;
    }

    id = id->next[0];
  }

  return 0;
} /* End of prunetraces() */

/***************************************************************************
 * findcoverage():
 *
 * Search an MS3TraceList for entries that overlap the target MS3TraceSeg
 * and, from the Record entries of the overlapping MS3TraceSegs, build a
 * coverage list.
 *
 * Only data with a higher priority than the target MS3TraceSeg will be
 * added to the overlap coverage.  Priority is determined using the
 * publication versions and if the versions are equal the
 * longest time-series will be given priority.
 *
 * On success a new Coverage will be allocated and returned, it is
 * up to the caller to properly free this memory.
 *
 * When no overlap coverage is found *ppcoverage will be 0, otherwise
 * it will contain a list of MS3Traces representing the overlap
 * coverage.
 *
 * Returns 0 on success and -1 on error.
 ***************************************************************************/
static int
findcoverage (MS3TraceList *mstl, MS3TraceID *targetid, MS3TraceSeg *targetseg,
              Coverage **ppcoverage)
{
  MS3TraceID *id = NULL;
  MS3TraceSeg *seg = NULL;
  Coverage *coverage = NULL;
  Coverage *prevcoverage = NULL;
  RecordMap *recmap;
  Record *rec;
  nstime_t nsdelta, nstimetol;
  nstime_t effstarttime, effendtime;
  int priority;
  int newsegment;

  if (!mstl || !targetseg || !ppcoverage)
    return -1;

  /* Determine sample period in high precision time ticks */
  nsdelta = (targetseg->samprate) ? (nstime_t)(NSTMODULUS / targetseg->samprate) : 0;

  /* Determine time tolerance in high precision time ticks */
  nstimetol = (timetol == -1) ? (nsdelta / 2) : (nstime_t)(NSTMODULUS * timetol);

  /* Loop through each MS3TraceID in the list */
  id = mstl->traces.next[0];
  while (id)
  {
    /* Continue with next if source ID is different */
    if (targetid != id)
    {
      if (strcmp (id->sid, targetid->sid))
      {
        id = id->next[0];
        continue;
      }
    }

    prevcoverage = NULL;
    seg = id->first;
    while (seg)
    {
      /* Skip target segment */
      if (seg == targetseg)
      {
        seg = seg->next;
        continue;
      }

      /* Stop searching if target segment is before segment start time,
       * assuming the segments are in time order nothing later will overlap. */
      if ((targetseg->endtime + nstimetol) < seg->starttime)
      {
        break;
      }

      /* Skip out-of-band (0 samprate) trace */
      if (seg->samprate == 0.0)
      {
        seg = seg->next;
        continue;
      }

      /* Continue with next if sample rate are different */
      if (!MS_ISRATETOLERABLE (seg->samprate, targetseg->samprate))
      {
        seg = seg->next;
        continue;
      }

      /* Check for duplicate or overlap NSLCQ last coverage entry */
      if (prevcoverage)
      {
        /* At this point the NSLCQ and rate are the same, check if the
         * segment is completly contained by the previous coverage entry. */
        if (seg->starttime >= prevcoverage->starttime &&
            seg->endtime <= prevcoverage->endtime)
        {
          seg = seg->next;
          continue;
        }
      }

      /* Test for overlap with targetseg */
      if ((targetseg->endtime + nstimetol) >= seg->starttime &&
          (targetseg->starttime - nstimetol) <= seg->endtime)
      {
        /* Determine priority:
         *  -1 : seg > targetseg
         *   0 : seg == targetseg
         *   1 : seg < targetseg */
        priority = 0;

        /* If best version is requested compare the qualities to determine priority */
        if (bestversion)
        {
          if (id->pubversion > targetid->pubversion)
            priority = -1;
          else if (id->pubversion < targetid->pubversion)
            priority = 1;
        }

        /* If priorities are equal (pubversions are equal or no checking)
         * give priority to the longest segment */
        if (priority == 0)
        {
          if ((seg->endtime - seg->starttime) >= (targetseg->endtime - targetseg->starttime))
            priority = -1;
          else
            priority = 1;
        }

        /* If overlapping trace is a higher priority than targetseg add to coverage */
        if (priority == -1)
        {
          /* Loop through list of associated Records, determine
           * contiguous coverage and store in an MSTraceGroup */
          recmap = (RecordMap *)seg->prvtptr;
          rec = recmap->first;
          newsegment = 1;
          while (rec)
          {
            /* Check if record has been marked as non-contributing */
            if (rec->reclen == 0)
            {
              rec = rec->next;
              continue;
            }

            /* Determine effective record start and end times */
            effstarttime = (rec->newstart != NSTUNSET) ? rec->newstart : rec->starttime;
            effendtime = (rec->newend != NSTUNSET) ? rec->newend : rec->endtime;

            /* Create a new segment if a break in the time-series is detected */
            if (coverage)
              if (llabs ((coverage->endtime + nsdelta) - effstarttime) > nstimetol)
                newsegment = 1;

            if (newsegment)
            {
              newsegment = 0;

              if ((coverage = (Coverage *)malloc (sizeof (Coverage))) == NULL)
              {
                ms_log (2, "Cannot allocate memory for coverage, bah humbug.\n");
                return -1;
              }

              prevcoverage = coverage;
              coverage->pubversion = id->pubversion;
              coverage->samprate = seg->samprate;
              coverage->starttime = effstarttime;
              *ppcoverage = coverage;
            }

            if (coverage)
              coverage->endtime = effendtime;
            else
              ms_log (2, "ACK! covergage is not allocated!?  PLEASE REPORT\n");

            rec = rec->next;
          }
        }
      }

      seg = seg->next;
    }

    id = id->next[0];
  }

  return 0;
} /* End of findcoverage() */

/***************************************************************************
 * trimtrace():
 *
 * Adjust Record entries associated with the target MS3TraceSeg that
 * are overlapping the time represented by the Coverage
 * in two different ways: 1) mark records that are completely
 * overlapped and 2) determine partial record trim boundaries (new
 * record times) if sample level pruning is requested.
 *
 * Completely overlapping Record entries are marked for omission by
 * setting Record.reclen = 0.  Partial Record overlaps are noted by
 * setting Record.newstart and Record.newend when sample level pruning
 * is requested.  The actual trimming of the data records, complete or
 * partial, is performed during the output sequence, not in this
 * routine.
 *
 * Returns the number of Record modifications on success and -1 on error.
 ***************************************************************************/
static int
trimtrace (MS3TraceSeg *targetseg, const char *targetsourceid, Coverage *coverage)
{
  RecordMap *recmap;
  Record *rec;
  Coverage *cov;
  nstime_t effstarttime, effendtime;
  nstime_t nsdelta, nstimetol;
  char stime[30] = {0};
  char etime[30] = {0};
  int modcount = 0;

  if (!targetseg || !coverage)
    return -1;

  /* Determine sample period in high precision time ticks */
  nsdelta = (targetseg->samprate) ? (nstime_t)(NSTMODULUS / targetseg->samprate) : 0;

  /* Determine time tolerance in high precision time ticks */
  nstimetol = (timetol == -1) ? (nsdelta / 2) : (nstime_t)(NSTMODULUS * timetol);

  /* Traverse the Record chain for the target MS3Trace and mark Records
   * that are completely overlapped by the MSTraceGroup coverage */
  recmap = (RecordMap *)targetseg->prvtptr;
  rec = recmap->first;
  while (rec)
  {
    cov = coverage;
    while (cov)
    {
      if (!rec->reclen) /* Skip if marked non-contributing */
        break;

      /* Determine effective record start and end times for comparison */
      effstarttime = (rec->newstart != NSTUNSET) ? rec->newstart : rec->starttime;
      effendtime = (rec->newend != NSTUNSET) ? rec->newend : rec->endtime;

      /* Use selection start and end for pruning if they are stricter */
      if (rec->selectstart != NSTUNSET && rec->selectstart > effstarttime)
        effstarttime = rec->selectstart;
      if (rec->selectend != NSTUNSET && rec->selectend < effendtime)
        effendtime = rec->selectend;

      /* Mark Record if it is completely overlapped by the coverage including tolerance */
      if (effstarttime >= (coverage->starttime - nstimetol) &&
          effendtime <= (coverage->endtime + nstimetol))
      {
        if (verbose > 1)
        {
          ms_nstime2timestr (rec->starttime, stime, ISOMONTHDAY_Z, NANO_MICRO);
          ms_nstime2timestr (rec->endtime, etime, ISOMONTHDAY_Z, NANO_MICRO);
          ms_log (1, "Removing Record %s (%u) :: %s  %s  offset: %lld, reclen: %d\n",
                  targetsourceid, rec->pubversion, stime, etime,
                  (long long int)rec->offset, rec->reclen);
        }

        rec->flp->recrmcount++;
        rec->reclen = 0;
        modcount++;
      }

      /* Determine the new start/end times if pruning at the sample level */
      if (prunedata == 's' && rec->reclen != 0)
      {
        /* Record overlaps beginning of coverage */
        if (effstarttime < coverage->starttime &&
            (effendtime + nstimetol) >= coverage->starttime)
        {
          /* Set Record new end time boundary including specified time tolerance */
          rec->newend = coverage->starttime - nsdelta + nstimetol;

          if ((starttime != NSTUNSET) && (rec->newend < starttime))
          {
            if (verbose > 1)
            {
              ms_nstime2timestr (rec->starttime, stime, ISOMONTHDAY_Z, NANO_MICRO);
              ms_nstime2timestr (rec->endtime, etime, ISOMONTHDAY_Z, NANO_MICRO);
              ms_log (1, "Removing Record %s (%u) :: %s  %s\n",
                      targetsourceid, rec->pubversion, stime, etime);
            }

            rec->flp->recrmcount++;
            rec->reclen = 0;
            modcount++;
          }
          else
          {
            effendtime = rec->newend;
            rec->flp->rectrimcount++;
            modcount++;
          }
        }

        /* Record overlaps end of coverage */
        if ((effstarttime - nstimetol) <= coverage->endtime &&
            effendtime > coverage->endtime)
        {
          /* Set Record new start time boundary including specified time tolerance */
          rec->newstart = coverage->endtime + nsdelta - nstimetol;

          if ((endtime != NSTUNSET) && (rec->newstart > endtime))
          {
            if (verbose > 1)
            {
              ms_nstime2timestr (rec->starttime, stime, ISOMONTHDAY_Z, NANO_MICRO);
              ms_nstime2timestr (rec->endtime, etime, ISOMONTHDAY_Z, NANO_MICRO);
              ms_log (1, "Removing Record %s (%u) :: %s  %s\n",
                      targetsourceid, rec->pubversion, stime, etime);
            }

            rec->flp->recrmcount++;
            rec->reclen = 0;
            modcount++;
          }
          else
          {
            effstarttime = rec->newstart;
            rec->flp->rectrimcount++;
            modcount++;
          }
        }

        /* Remove record if all samples have been pruned within tolerance,
         * test for special cases of:
         * a) no time coverage (single sample) and no pruning
         * b) no time coverage (single last sample) and split boundary usage */
        if (effstarttime >= (effendtime - nstimetol) &&
            !(rec->starttime == rec->endtime &&
              rec->starttime == effstarttime &&
              rec->endtime == effendtime) &&
            !(splitboundary &&
              (effstarttime == effendtime && effendtime == rec->endtime)))
        {
          if (verbose > 1)
          {
            ms_nstime2timestr (rec->starttime, stime, ISOMONTHDAY_Z, NANO_MICRO);
            ms_nstime2timestr (rec->endtime, etime, ISOMONTHDAY_Z, NANO_MICRO);
            ms_log (1, "Removing Record %s (%u) :: %s  %s\n",
                    targetsourceid, rec->pubversion, stime, etime);
          }

          rec->flp->recrmcount++;
          rec->reclen = 0;
          modcount++;
        }

      } /* Done pruning at sample level */

      cov = cov->next;
    }

    rec = rec->next;
  }

  return modcount;
} /* End of trimtrace() */

/***************************************************************************
 * reconcile_tracetimes():
 *
 * Reconcile the start and end times of the traces in a specified
 * trace group with the list of records in an associated record map.
 * In other words, set the start and end times of each MS3TraceSeg in
 * the MS3TraceList according to the start time of the first and end
 * time of the last contributing records in the associated record map;
 * this should be performed after the pruning process which could mark
 * complete records as pruned (non-contributing).
 *
 * Returns 0 on success and -1 on error.
 ***************************************************************************/
static int
reconcile_tracetimes (MS3TraceList *mstl)
{
  MS3TraceID *id;
  MS3TraceSeg *seg;
  RecordMap *recmap;
  Record *rec;
  Record *first = NULL;
  Record *last = NULL;

  if (!mstl)
    return -1;

  if (!mstl->traces.next[0])
    return -1;

  id = mstl->traces.next[0];
  while (id)
  {
    seg = id->first;
    while (seg)
    {
      recmap = (RecordMap *)seg->prvtptr;

      /* Find first contributing record (reclen != 0) */
      rec = recmap->first;
      while (rec)
      {
        if (rec->reclen > 0)
        {
          first = rec;
          break;
        }

        rec = rec->next;
      }

      /* Find last contributing record (reclen != 0) */
      rec = recmap->last;
      while (rec)
      {
        if (rec->reclen > 0)
        {
          last = rec;
          break;
        }

        rec = rec->prev;
      }

      /* Set a new MS3TraceSeg start time */
      if (first)
      {
        /* Use the new boundary start time if set and sane */
        if (first->newstart != NSTUNSET && first->newstart > first->starttime)
          seg->starttime = first->newstart;
        /* Otherwise use the record start time */
        else
          seg->starttime = first->starttime;
      }

      /* Set a new MS3TraceSeg end time */
      if (last)
      {
        /* Use the new boundary end time if set and sane */
        if (last->newend != NSTUNSET && last->newend < last->endtime)
          seg->endtime = last->newend;
        /* Otherwise use the record end time */
        else
          seg->endtime = last->endtime;
      }

      first = NULL;
      last = NULL;
      seg = seg->next;
    }

    id = id->next[0];
  }

  return 0;
} /* End of reconcile_tracetimes() */

/***************************************************************************
 * printmodsummary():
 *
 * Print a summary of modifications to stdout.  If 'nomods' is true
 * include files that were not modified.
 ***************************************************************************/
static void
printmodsummary (flag nomods)
{
  Filelink *flp;

  ms_log (0, "File modification summary:\n");

  flp = filelist;

  while (flp)
  {
    if (!nomods && !flp->reordercount && !flp->recrmcount && !flp->rectrimcount)
    {
      flp = flp->next;
      continue;
    }

    ms_log (0, " Records split: %" PRId64 " trimmed: %" PRId64 " removed: %" PRId64 ", Segments reordered: %" PRId64 " :: %s\n",
            flp->recsplitcount, flp->rectrimcount, flp->recrmcount, flp->reordercount, flp->infilename);

    flp = flp->next;
  }

  return;
} /* End of printmodsummary() */

/***************************************************************************
 * printtracemap():
 *
 * Print record map for each MS3TraceSeg to stdout.
 ***************************************************************************/
static void
printtracemap (MS3TraceList *mstl)
{
  MS3TraceID *id = NULL;
  MS3TraceSeg *seg = NULL;
  char stime[30] = {0};
  char etime[30] = {0};
  int segcnt = 0;

  if (!mstl)
    return;

  id = mstl->traces.next[0];

  /* Print out the appropriate header */
  ms_log (0, "\nTrace Map:\n");
  ms_log (0, "   Source              Start sample             End sample        Hz   Samples\n");

  while (id)
  {
    seg = id->first;

    while (seg)
    {
      /* Create formatted time strings */
      if (ms_nstime2timestr (seg->starttime, stime, ISOMONTHDAY_Z, NANO_MICRO) == NULL)
        ms_log (2, "Cannot convert trace start time for %s\n", id->sid);

      if (ms_nstime2timestr (seg->endtime, etime, ISOMONTHDAY_Z, NANO_MICRO) == NULL)
        ms_log (2, "Cannot convert trace end time for %s\n", id->sid);

      /* Print MS3TraceSeg header */
      ms_log (0, "%-15s %-24s %-24s %-4.4g %-" PRId64 "\n",
              id->sid, stime, etime, seg->samprate, seg->samplecnt);

      if (!seg->prvtptr)
      {
        ms_log (2, "No record map associated with this MS3TraceSeg.\n");
      }
      else
      {
        printrecordmap ((RecordMap *)seg->prvtptr, 0);
      }

      segcnt++;
      seg = seg->next;
    }

    id = id->next[0];
  }

  ms_log (0, "End of trace map: %d trace segment(s)\n\n", segcnt);

} /* End of printtracemap() */

/***************************************************************************
 * printrecordmap():
 *
 * Print record map to stdout.
 ***************************************************************************/
static void
printrecordmap (RecordMap *recmap, flag details)
{
  char stime[30] = {0};
  char etime[30] = {0};
  Record *rec;

  if (!recmap)
    return;

  rec = recmap->first;

  ms_log (0, "Record map contains %lld records:\n", recmap->recordcnt);

  while (rec)
  {
    ms_log (0, "  Filename: %s  Offset: %llu  RecLen: %d  PubVersion: %u\n",
            rec->flp->infilename, (long long unsigned)rec->offset, rec->reclen, rec->pubversion);

    ms_nstime2timestr (rec->starttime, stime, ISOMONTHDAY_Z, NANO_MICRO);
    ms_nstime2timestr (rec->endtime, etime, ISOMONTHDAY_Z, NANO_MICRO);
    ms_log (0, "        Start: %s       End: %s\n", stime, etime);

    if (details)
    {
      if (rec->newstart == NSTUNSET)
        strcpy (stime, "NONE");
      else
        ms_nstime2timestr (rec->newstart, stime, ISOMONTHDAY_Z, NANO_MICRO);
      if (rec->newend == NSTUNSET)
        strcpy (etime, "NONE");
      else
        ms_nstime2timestr (rec->newend, etime, ISOMONTHDAY_Z, NANO_MICRO);
      ms_log (0, "  Start bound: %-24s  End bound: %-24s\n", stime, etime);

      if (rec->selectstart == NSTUNSET)
        strcpy (stime, "NONE");
      else
        ms_nstime2timestr (rec->selectstart, stime, ISOMONTHDAY_Z, NANO_MICRO);
      if (rec->selectend == NSTUNSET)
        strcpy (etime, "NONE");
      else
        ms_nstime2timestr (rec->selectend, etime, ISOMONTHDAY_Z, NANO_MICRO);
      ms_log (0, " Select start: %-24s Select end: %-24s\n", stime, etime);
    }

    rec = rec->next;
  }
} /* End of printrecordmap() */

/***************************************************************************
 * printwritten():
 *
 * Print summary of output records.
 ***************************************************************************/
static void
printwritten (MS3TraceList *mstl)
{
  MS3TraceID *id = NULL;
  MS3TraceSeg *seg = NULL;
  char stime[30] = {0};
  char etime[30] = {0};
  FILE *ofp;

  if (!mstl)
    return;

  if (strcmp (writtenfile, "-") == 0)
  {
    ofp = stdout;
  }
  else if (strcmp (writtenfile, "--") == 0)
  {
    ofp = stderr;
  }
  else if ((ofp = fopen (writtenfile, "ab")) == NULL)
  {
    ms_log (2, "Cannot open output file: %s (%s)\n",
            writtenfile, strerror (errno));
    return;
  }

  /* Loop through trace list */
  id = mstl->traces.next[0];
  while (id)
  {
    /* Loop through segment list */
    seg = id->first;
    while (seg)
    {
      if (ms_nstime2timestr (seg->starttime, stime, ISOMONTHDAY_Z, NANO_MICRO) == NULL)
        ms_log (2, "Cannot convert trace start time for %s\n", id->sid);

      if (ms_nstime2timestr (seg->endtime, etime, ISOMONTHDAY_Z, NANO_MICRO) == NULL)
        ms_log (2, "Cannot convert trace end time for %s\n", id->sid);

      fprintf (ofp, "%s%s|%u|%-30s|%-30s|%lld|%lld\n",
               (writtenprefix) ? writtenprefix : "",
               id->sid, id->pubversion, stime, etime,
               (long long int)*((int64_t *)seg->prvtptr),
               (long long int)seg->samplecnt);

      seg = seg->next;
    }

    id = id->next[0];
  }

  if (ofp != stdout && fclose (ofp))
    ms_log (2, "Cannot close output file: %s (%s)\n",
            writtenfile, strerror (errno));

} /* End of printwritten() */

/***************************************************************************
 * findselectlimits():
 *
 * Determine selection limits for the given record based on all
 * matching selection entries.
 *
 * Return 0 on success and -1 on error.
 ***************************************************************************/
static int
findselectlimits (const MS3Selections *select, const char *sid, nstime_t starttime,
                  nstime_t endtime, Record *rec)
{
  const MS3SelectTime *selecttime = NULL;
  char timestring[100];

  if (!rec || !sid || !select)
    return -1;

  while ((select = ms3_matchselect (select, sid, starttime, endtime, 0, &selecttime)))
  {
    while (selecttime)
    {
      /* Continue if selection edge time does not intersect with record coverage */
      if ((starttime < selecttime->starttime && !(starttime <= selecttime->starttime && endtime >= selecttime->starttime)))
      {
        selecttime = selecttime->next;
        continue;
      }
      else if ((endtime > selecttime->endtime && !(starttime <= selecttime->endtime && endtime >= selecttime->endtime)))
      {
        selecttime = selecttime->next;
        continue;
      }

      /* Check that the selection intersects previous selection range if set,
       * otherwise the combined selection is not possible. */
      if (rec->selectstart != NSTUNSET && rec->selectend != NSTUNSET &&
          !(rec->selectstart <= selecttime->endtime && rec->selectend >= selecttime->starttime))
      {
        ms_nstime2timestr (starttime, timestring, ISOMONTHDAY_Z, NANO_MICRO);
        ms_log (1, "Warning: impossible combination of selections for record (%s, %s), not pruning.\n",
                sid, timestring);
        rec->selectstart = NSTUNSET;
        rec->selectend = NSTUNSET;
        return 0;
      }

      if (rec->selectstart == NSTUNSET || rec->selectstart > selecttime->starttime)
      {
        rec->selectstart = selecttime->starttime;
      }

      if (rec->selectend == NSTUNSET || rec->selectend < selecttime->endtime)
      {
        rec->selectend = selecttime->endtime;
      }

      /* Shortcut if the entire record is already selected */
      if (rec->starttime >= rec->selectstart && rec->endtime <= rec->selectend)
        return 0;

      selecttime = selecttime->next;
    }

    select = select->next;
  }

  return 0;
} /* End of findselectlimits() */

/***************************************************************************
 * sortrecmap():
 *
 * Sort a RecordMap so that records are in time order using a
 * mergesort algorithm.
 *
 * The mergesort implementation was inspired by the listsort function
 * published and copyright 2001 by Simon Tatham.
 *
 * Return 0 on success and -1 on error.
 ***************************************************************************/
static int
sortrecmap (RecordMap *recmap)
{
  Record *p, *q, *e, *top, *tail;
  int nmerges;
  int insize, psize, qsize, i;

  if (!recmap)
    return -1;

  if (!recmap->first || !recmap->last) /* Done if empty */
    return 0;

  top = recmap->first;
  insize = 1;

  for (;;)
  {
    p = top;
    top = NULL;
    tail = NULL;

    nmerges = 0; /* count number of merges we do in this pass */

    while (p)
    {
      nmerges++; /* there exists a merge to be done */

      /* step `insize' places along from p */
      q = p;
      psize = 0;
      for (i = 0; i < insize; i++)
      {
        psize++;
        q = q->next;
        if (!q)
          break;
      }

      /* if q hasn't fallen off end, we have two lists to merge */
      qsize = insize;

      /* now we have two lists; merge them */
      while (psize > 0 || (qsize > 0 && q))
      {
        /* decide whether next element of merge comes from p or q */
        if (psize == 0)
        { /* p is empty; e must come from q. */
          e = q;
          q = q->next;
          qsize--;
        }
        else if (qsize == 0 || !q)
        { /* q is empty; e must come from p. */
          e = p;
          p = p->next;
          psize--;
        }
        else if (recordcmp (p, q) <= 0)
        { /* First element of p is lower (or same), e must come from p. */
          e = p;
          p = p->next;
          psize--;
        }
        else
        { /* First element of q is lower; e must come from q. */
          e = q;
          q = q->next;
          qsize--;
        }

        /* add the next element to the merged list */
        if (tail)
          tail->next = e;
        else
          top = e;

        tail = e;
      }

      /* now p has stepped `insize' places along, and q has too */
      p = q;
    }

    tail->next = NULL;

    /* If we have done only one merge, we're finished. */
    if (nmerges <= 1) /* allow for nmerges==0, the empty list case */
    {
      recmap->first = top;
      recmap->last = tail;

      return 0;
    }

    /* Otherwise repeat, merging lists twice the size */
    insize *= 2;
  }
} /* End of sortrecmap() */

/***************************************************************************
 * recordcmp():
 *
 * Compare the start times of each Record for the purposes of sorting
 * a RecordMap.
 *
 * Return 1 if rec1 is "greater" than rec2, otherwise return 0.
 ***************************************************************************/
static int
recordcmp (Record *rec1, Record *rec2)
{
  nstime_t *start1;
  nstime_t *start2;

  if (!rec1 || !rec2)
    return -1;

  /* Determine effective start times */
  start1 = (rec1->newstart != NSTUNSET) ? &(rec1->newstart) : &(rec1->starttime);
  start2 = (rec2->newstart != NSTUNSET) ? &(rec2->newstart) : &(rec2->starttime);

  if (*start1 > *start2)
  {
    return 1;
  }

  return 0;
} /* End of recordcmp() */

/***************************************************************************
 * processparam():
 * Process the command line parameters.
 *
 * Returns 0 on success, and -1 on failure
 ***************************************************************************/
static int
processparam (int argcount, char **argvec)
{
  int optind;
  char *selectfile = NULL;
  char *matchpattern = NULL;
  char *rejectpattern = NULL;
  char *tptr;

  /* Process all command line arguments */
  for (optind = 1; optind < argcount; optind++)
  {
    if (strcmp (argvec[optind], "-V") == 0)
    {
      ms_log (1, "%s version: %s\n", PACKAGE, VERSION);
      exit (0);
    }
    else if (strcmp (argvec[optind], "-h") == 0)
    {
      usage (0);
      exit (0);
    }
    else if (strcmp (argvec[optind], "-H") == 0)
    {
      usage (1);
      exit (0);
    }
    else if (strncmp (argvec[optind], "-v", 2) == 0)
    {
      verbose += strspn (&argvec[optind][1], "v");
    }
    else if (strcmp (argvec[optind], "-tt") == 0)
    {
      timetol = strtod (getoptval (argcount, argvec, optind++), NULL);
      tolerance.time = timetol_callback;
    }
    else if (strcmp (argvec[optind], "-rt") == 0)
    {
      sampratetol = strtod (getoptval (argcount, argvec, optind++), NULL);
      tolerance.samprate = samprate_callback;
    }
    else if (strcmp (argvec[optind], "-E") == 0)
    {
      bestversion = 0;
    }
    else if (strcmp (argvec[optind], "-s") == 0)
    {
      selectfile = getoptval (argcount, argvec, optind++);
    }
    else if (strcmp (argvec[optind], "-ts") == 0)
    {
      starttime = ms_timestr2nstime (getoptval (argcount, argvec, optind++));
      if (starttime == NSTERROR)
        return -1;
    }
    else if (strcmp (argvec[optind], "-te") == 0)
    {
      endtime = ms_timestr2nstime (getoptval (argcount, argvec, optind++));
      if (endtime == NSTERROR)
        return -1;
    }
    else if (strcmp (argvec[optind], "-M") == 0)
    {
      matchpattern = strdup (getoptval (argcount, argvec, optind++));
    }
    else if (strcmp (argvec[optind], "-R") == 0)
    {
      rejectpattern = strdup (getoptval (argcount, argvec, optind++));
    }
    else if (strcmp (argvec[optind], "-m") == 0)
    {
      tptr = getoptval (argcount, argvec, optind++);

      if (ms3_addselect (&selections, tptr, NSTUNSET, NSTUNSET, 0) < 0)
      {
        ms_log (2, "Unable to add selection: '%s'\n", tptr);
        return -1;
      }
    }
    else if (strcmp (argvec[optind], "-o") == 0)
    {
      outputfile = getoptval (argcount, argvec, optind++);
      outputmode = 0;
    }
    else if (strcmp (argvec[optind], "+o") == 0)
    {
      outputfile = getoptval (argcount, argvec, optind++);
      outputmode = 1;
    }
    else if (strcmp (argvec[optind], "-A") == 0)
    {
      if (addarchive (getoptval (argcount, argvec, optind++), NULL) == -1)
        return -1;
    }
    else if (strcmp (argvec[optind], "-Pr") == 0)
    {
      prunedata = 'r';
    }
    else if (strcmp (argvec[optind], "-Ps") == 0 || strcmp (argvec[optind], "-P") == 0)
    {
      prunedata = 's';
    }
    else if (strcmp (argvec[optind], "-Pe") == 0)
    {
      prunedata = 'e';
    }
    else if (strcmp (argvec[optind], "-Sd") == 0)
    {
      splitboundary = 'd';
    }
    else if (strcmp (argvec[optind], "-Sh") == 0)
    {
      splitboundary = 'h';
    }
    else if (strcmp (argvec[optind], "-Sm") == 0)
    {
      splitboundary = 'm';
    }
    else if (strcmp (argvec[optind], "-rls") == 0)
    {
      splitreclen = 1;
    }
    else if (strcmp (argvec[optind], "-Q") == 0)
    {
      tptr = getoptval (argcount, argvec, optind++);
      restampqind = *tptr;

      if (restampqind != 'D' && restampqind != 'R' &&
          restampqind != 'Q' && restampqind != 'M')
      {
        ms_log (2, "Invalid data indicator: '%c'\n", restampqind);
        exit (1);
      }
    }
    else if (strcmp (argvec[optind], "-sum") == 0)
    {
      basicsum = 1;
    }
    else if (strcmp (argvec[optind], "-mod") == 0)
    {
      modsummary = 1;
    }
    else if (strcmp (argvec[optind], "-out") == 0)
    {
      writtenfile = getoptval (argcount, argvec, optind++);
    }
    else if (strcmp (argvec[optind], "-outprefix") == 0)
    {
      writtenprefix = getoptval (argcount, argvec, optind++);
    }
    else if (strcmp (argvec[optind], "-CHAN") == 0)
    {
      if (addarchive (getoptval (argcount, argvec, optind++), CHANLAYOUT) == -1)
        return -1;
    }
    else if (strcmp (argvec[optind], "-QCHAN") == 0)
    {
      if (addarchive (getoptval (argcount, argvec, optind++), QCHANLAYOUT) == -1)
        return -1;
    }
    else if (strcmp (argvec[optind], "-CDAY") == 0)
    {
      if (addarchive (getoptval (argcount, argvec, optind++), CDAYLAYOUT) == -1)
        return -1;
    }
    else if (strcmp (argvec[optind], "-SDAY") == 0)
    {
      if (addarchive (getoptval (argcount, argvec, optind++), SDAYLAYOUT) == -1)
        return -1;
    }
    else if (strcmp (argvec[optind], "-BUD") == 0)
    {
      if (addarchive (getoptval (argcount, argvec, optind++), BUDLAYOUT) == -1)
        return -1;
    }
    else if (strcmp (argvec[optind], "-SDS") == 0)
    {
      if (addarchive (getoptval (argcount, argvec, optind++), SDSLAYOUT) == -1)
        return -1;
    }
    else if (strcmp (argvec[optind], "-CSS") == 0)
    {
      if (addarchive (getoptval (argcount, argvec, optind++), CSSLAYOUT) == -1)
        return -1;
    }
    else if (strncmp (argvec[optind], "-", 1) == 0 &&
             strlen (argvec[optind]) > 1)
    {
      ms_log (2, "Unknown option: %s\n", argvec[optind]);
      exit (1);
    }
    else
    {
      tptr = argvec[optind];

      /* Check for an input file list */
      if (tptr[0] == '@')
      {
        if (addlistfile (tptr + 1) < 0)
        {
          ms_log (2, "Error adding list file %s", tptr + 1);
          exit (1);
        }
      }
      /* Otherwise this is an input file */
      else
      {
        /* Add file to global file list */
        if (addfile (tptr))
        {
          ms_log (2, "Error adding file to input list %s", tptr);
          exit (1);
        }
      }
    }
  }

  /* Make sure input file(s) were specified */
  if (!filelist)
  {
    ms_log (2, "No input files were specified\n\n");
    ms_log (1, "%s version %s\n\n", PACKAGE, VERSION);
    ms_log (1, "Try %s -h for usage\n", PACKAGE);
    exit (0);
  }

  /* Make sure output file(s) were specified or replacing originals */
  if (!archiveroot && !outputfile)
  {
    ms_log (2, "No output files were specified\n\n");
    ms_log (1, "%s version %s\n\n", PACKAGE, VERSION);
    ms_log (1, "Try %s -h for usage\n", PACKAGE);
    exit (0);
  }

  /* Read data selection file */
  if (selectfile)
  {
    if (ms3_readselectionsfile (&selections, selectfile) < 0)
    {
      ms_log (2, "Cannot read data selection file\n");
      exit (1);
    }
  }

  /* Expand match pattern from a file if prefixed by '@' */
  if (matchpattern)
  {
    if (*matchpattern == '@')
    {
      tptr = strdup (matchpattern + 1); /* Skip the @ sign */
      free (matchpattern);
      matchpattern = NULL;

      if (readregexfile (tptr, &matchpattern) <= 0)
      {
        ms_log (2, "Cannot read match pattern regex file\n");
        exit (1);
      }

      free (tptr);
    }
  }

  /* Expand reject pattern from a file if prefixed by '@' */
  if (rejectpattern)
  {
    if (*rejectpattern == '@')
    {
      tptr = strdup (rejectpattern + 1); /* Skip the @ sign */
      free (rejectpattern);
      rejectpattern = NULL;

      if (readregexfile (tptr, &rejectpattern) <= 0)
      {
        ms_log (2, "Cannot read reject pattern regex file\n");
        exit (1);
      }

      free (tptr);
    }
  }

  /* Compile match and reject patterns */
  if (matchpattern)
  {
    if (!(match = (regex_t *)malloc (sizeof (regex_t))))
    {
      ms_log (2, "Cannot allocate memory for match expression\n");
      exit (1);
    }

    if (regcomp (match, matchpattern, REG_EXTENDED) != 0)
    {
      ms_log (2, "Cannot compile match regex: '%s'\n", matchpattern);
    }

    free (matchpattern);
  }

  if (rejectpattern)
  {
    if (!(reject = (regex_t *)malloc (sizeof (regex_t))))
    {
      ms_log (2, "Cannot allocate memory for reject expression\n");
      exit (1);
    }

    if (regcomp (reject, rejectpattern, REG_EXTENDED) != 0)
    {
      ms_log (2, "Cannot compile reject regex: '%s'\n", rejectpattern);
    }

    free (rejectpattern);
  }

  /* Report the program version */
  if (verbose)
    ms_log (1, "%s version: %s\n", PACKAGE, VERSION);

  return 0;
} /* End of processparam() */

/***************************************************************************
 * getoptval:
 * Return the value to a command line option; checking that the value is
 * itself not an option (starting with '-') and is not past the end of
 * the argument list.
 *
 * argcount: total arguments in argvec
 * argvec: argument list
 * argopt: index of option to process, value is expected to be at argopt+1
 *
 * Returns value on success and exits with error message on failure
 ***************************************************************************/
static char *
getoptval (int argcount, char **argvec, int argopt)
{
  if (argvec == NULL || argvec[argopt] == NULL)
  {
    ms_log (2, "getoptval(): NULL option requested\n");
    exit (1);
    return 0;
  }

  /* Special case of '-o -' usage */
  if ((argopt + 1) < argcount && strcmp (argvec[argopt], "-o") == 0)
    if (strcmp (argvec[argopt + 1], "-") == 0)
      return argvec[argopt + 1];

  /* Special case of '+o -' usage */
  if ((argopt + 1) < argcount && strcmp (argvec[argopt], "+o") == 0)
    if (strcmp (argvec[argopt + 1], "-") == 0)
      return argvec[argopt + 1];

  /* Special case of '-s -' usage */
  if ((argopt + 1) < argcount && strcmp (argvec[argopt], "-s") == 0)
    if (strcmp (argvec[argopt + 1], "-") == 0)
      return argvec[argopt + 1];

  /* Special case of '-out -' or '-out --' usage */
  if ((argopt + 1) < argcount && strcmp (argvec[argopt], "-out") == 0)
    if (strcmp (argvec[argopt + 1], "-") == 0 ||
        strcmp (argvec[argopt + 1], "--") == 0)
      return argvec[argopt + 1];

  if ((argopt + 1) < argcount && *argvec[argopt + 1] != '-')
    return argvec[argopt + 1];

  ms_log (2, "Option %s requires a value, try -h for usage\n", argvec[argopt]);
  exit (1);
  return 0;
} /* End of getoptval() */

/***************************************************************************
 * setofilelimit:
 *
 * Check the current open file limit and if it is not >= 'limit' try
 * to increase it to 'limit'.
 *
 * Returns the open file limit on success and -1 on error.
 ***************************************************************************/
static int
setofilelimit (int limit)
{
  struct rlimit rlim;
  rlim_t oldlimit;

  /* Get the current soft open file limit */
  if (getrlimit (RLIMIT_NOFILE, &rlim) == -1)
  {
    ms_log (2, "getrlimit() failed to get open file limit\n");
    return -1;
  }

  if (rlim.rlim_cur < limit)
  {
    oldlimit = rlim.rlim_cur;
    rlim.rlim_cur = (rlim_t)limit;

    if (verbose > 1)
      ms_log (1, "Setting open file limit to %d\n",
              (int)rlim.rlim_cur);

    if (setrlimit (RLIMIT_NOFILE, &rlim) == -1)
    {
      ms_log (2, "setrlimit failed to raise open file limit from %d to %d (max: %d)\n",
              (int)oldlimit, limit, (int)rlim.rlim_max);
      return -1;
    }
  }

  return (int)rlim.rlim_cur;
} /* End of setofilelimit() */

/***************************************************************************
 * addfile:
 *
 * Add file to end of the specified file list.
 *
 * Check for and parse start and end byte offsets (a read range)
 * embedded in the file name.  The form for specifying a read range is:
 *  filename@startoffset:endoffset
 * where both start and end offsets are optional.

 * Returns 0 on success and -1 on error.
 ***************************************************************************/
static int
addfile (char *filename)
{
  Filelink *newlp;
  char *at;
  char *colon;

  if (!filename)
  {
    ms_log (2, "%s(): No file name specified\n", __func__);
    return -1;
  }

  if (!(newlp = (Filelink *)calloc (1, sizeof (Filelink))))
  {
    ms_log (2, "%s(): Cannot allocate memory, out of memory?\n", __func__);
    return -1;
  }

  if (!(newlp->infilename_raw = strdup (filename)))
  {
    ms_log (2, "%s(): Cannot allocate memory, out of memory?\n", __func__);
    return -1;
  }

  /* Check for optional read byte range specifiers
   * Convert legacy byte separator of ":" to "-" as used by libmseed
   * Legacy form: "filename@startoffset:endoffset"
   * Needed form: "filename@startoffset-endoffset"
   */
  if ((at = strrchr (newlp->infilename_raw, '@')))
  {
    if ((colon = strrchr (at, ':')))
    {
      *colon = '-';
    }
  }

  if (!(newlp->infilename = strndup (filename, strcspn (filename, "@"))))
  {
    ms_log (2, "%s(): Cannot duplicate string, out of memory?\n", __func__);
    return -1;
  }

  /* Add new file to the end of the list */
  if (filelisttail == NULL)
  {
    filelist = newlp;
    filelisttail = newlp;
  }
  else
  {
    filelisttail->next = newlp;
    filelisttail = newlp;
  }

  return 0;
} /* End of addfile() */

/***************************************************************************
 * addlistfile:
 *
 * Add files listed in the specified file to the global input file list.
 *
 * Returns count of files added on success and -1 on error.
 ***************************************************************************/
static int
addlistfile (char *filename)
{
  FILE *fp;
  char filelistent[1024];
  int filecount = 0;

  if (verbose >= 1)
    ms_log (1, "Reading list file '%s'\n", filename);

  if (!(fp = fopen (filename, "rb")))
  {
    ms_log (2, "Cannot open list file %s: %s\n", filename, strerror (errno));
    return -1;
  }

  while (fgets (filelistent, sizeof (filelistent), fp))
  {
    char *cp;

    /* End string at first newline character */
    if ((cp = strchr (filelistent, '\n')))
      *cp = '\0';

    /* Skip empty lines */
    if (!strlen (filelistent))
      continue;

    /* Skip comment lines */
    if (*filelistent == '#')
      continue;

    if (verbose > 1)
      ms_log (1, "Adding '%s' from list file\n", filelistent);

    if (addfile (filelistent))
      return -1;

    filecount++;
  }

  fclose (fp);

  return filecount;
} /* End of addlistfile() */

/***************************************************************************
 * addarchive:
 * Add entry to the data stream archive chain.  'layout' if defined
 * will be appended to 'path'.
 *
 * Returns 0 on success, and -1 on failure
 ***************************************************************************/
static int
addarchive (const char *path, const char *layout)
{
  Archive *newarch;
  int pathlayout;

  if (!path)
  {
    ms_log (2, "addarchive(): cannot add archive with empty path\n");
    return -1;
  }

  if (!(newarch = (Archive *)malloc (sizeof (Archive))))
  {
    ms_log (2, "addarchive(): cannot allocate memory for new archive definition\n");
    return -1;
  }

  /* Setup new entry and add it to the front of the chain */
  pathlayout = strlen (path) + 2;
  if (layout)
    pathlayout += strlen (layout);

  if (!(newarch->datastream.path = (char *)malloc (pathlayout)))
  {
    ms_log (2, "addarchive(): cannot allocate memory for new archive path\n");
    if (newarch)
      free (newarch);
    return -1;
  }

  if (layout)
    snprintf (newarch->datastream.path, pathlayout, "%s/%s", path, layout);
  else
    snprintf (newarch->datastream.path, pathlayout, "%s", path);

  newarch->datastream.idletimeout = 60;
  newarch->datastream.grouproot = NULL;

  newarch->next = archiveroot;
  archiveroot = newarch;

  return 0;
} /* End of addarchive() */

/***************************************************************************
 * readregexfile:
 *
 * Read a list of regular expressions from a file and combine them
 * into a single, compound expression which is returned in *pppattern.
 * The return buffer is reallocated as need to hold the growing
 * pattern.  When called *pppattern should not point to any associated
 * memory.
 *
 * Returns the number of regexes parsed from the file or -1 on error.
 ***************************************************************************/
static int
readregexfile (char *regexfile, char **pppattern)
{
  FILE *fp;
  char line[1024];
  char linepattern[1024];
  int regexcnt = 0;
  int lengthbase;
  int lengthadd;

  if (!regexfile)
  {
    ms_log (2, "readregexfile: regex file not supplied\n");
    return -1;
  }

  if (!pppattern)
  {
    ms_log (2, "readregexfile: pattern string buffer not supplied\n");
    return -1;
  }

  /* Open the regex list file */
  if ((fp = fopen (regexfile, "rb")) == NULL)
  {
    ms_log (2, "Cannot open regex list file %s: %s\n",
            regexfile, strerror (errno));
    return -1;
  }

  if (verbose)
    ms_log (1, "Reading regex list from %s\n", regexfile);

  *pppattern = NULL;

  while ((fgets (line, sizeof (line), fp)) != NULL)
  {
    /* Trim spaces and skip if empty lines */
    if (sscanf (line, " %s ", linepattern) != 1)
      continue;

    /* Skip comment lines */
    if (*linepattern == '#')
      continue;

    regexcnt++;

    /* Add regex to compound regex */
    if (*pppattern)
    {
      lengthbase = strlen (*pppattern);
      lengthadd = strlen (linepattern) + 4; /* Length of addition plus 4 characters: |()\0 */

      *pppattern = realloc (*pppattern, lengthbase + lengthadd);

      if (*pppattern)
      {
        snprintf ((*pppattern) + lengthbase, lengthadd, "|(%s)", linepattern);
      }
      else
      {
        ms_log (2, "Cannot allocate memory for regex string\n");
        return -1;
      }
    }
    else
    {
      lengthadd = strlen (linepattern) + 3; /* Length of addition plus 3 characters: ()\0 */

      *pppattern = malloc (lengthadd);

      if (*pppattern)
      {
        snprintf (*pppattern, lengthadd, "(%s)", linepattern);
      }
      else
      {
        ms_log (2, "Cannot allocate memory for regex string\n");
        return -1;
      }
    }
  }

  fclose (fp);

  return regexcnt;
} /* End of readregexfile() */

/***************************************************************************
 * usage():
 * Print the usage message.
 ***************************************************************************/
static void
usage (int level)
{
  fprintf (stderr, "%s - select, sort and prune miniSEED: %s\n\n", PACKAGE, VERSION);
  fprintf (stderr, "Usage: %s [options] file1 [file2] [file3] ...\n\n", PACKAGE);
  fprintf (stderr,
           " ## Options ##\n"
           " -V           Report program version\n"
           " -h           Show this usage message\n"
           " -H           Show usage message with 'format' details (see -A option)\n"
           " -v           Be more verbose, multiple flags can be used\n"
           " -tt secs     Specify a time tolerance for continuous traces\n"
           " -rt diff     Specify a sample rate tolerance for continuous traces\n"
           " -E           Consider all qualities equal instead of 'best' prioritization\n"
           "\n"
           " ## Data selection options ##\n"
           " -s file      Specify a file containing selection criteria\n"
           " -ts time     Limit to records that contain or start after time\n"
           " -te time     Limit to records that contain or end before time\n"
           "                time format: 'YYYY[,DDD,HH,MM,SS,FFFFFF]' delimiters: [,:.]\n"
           " -M match     Limit to records matching the specified regular expression\n"
           " -R reject    Limit to records not matching the specfied regular expression\n"
           "                Regular expressions are applied to: 'NET_STA_LOC_CHAN_QUAL'\n"
           "\n"
           " ## Output options ##\n"
           " -o file      Specify a single output file, use +o file to append\n"
           " -A format    Write all records in a custom directory/file layout (try -H)\n"
           " -Pr          Prune data at the record level using 'best' version priority\n"
           " -Ps          Prune data at the sample level using 'best' version priority\n"
           " -Pe          Prune traces at user specified edges only, leave overlaps\n"
           " -S[dhm]      Split records on day, hour or minute boundaries\n"
           " -rls         Add suffixes to output files to split on record length changes\n"
           " -Q DRQM      Re-stamp output data records with quality code: D, R, Q or M\n"
           "\n"
           " ## Diagnostic output ##\n"
           " -sum         Print a basic summary after reading all input files\n"
           " -mod         Print summary of file modifications after processing\n"
           " -out file    Write a summary of output records to specified file\n"
           " -outprefix X Include prefix on summary output lines for identification\n"
           "\n"
           " ## Input data ##\n"
           " file#        Files(s) of miniSEED records\n"
           "\n");

  if (level)
  {
    fprintf (stderr,
             "\n"
             "  # Preset format layouts #\n"
             " -CHAN dir    Write records into separate Net.Sta.Loc.Chan files\n"
             " -QCHAN dir   Write records into separate Net.Sta.Loc.Chan.Quality files\n"
             " -CDAY dir    Write records into separate Net.Sta.Loc.Chan.Year:Yday:<time> files\n"
             " -SDAY dir    Write records into separate Net.Sta.Year:Yday files\n"
             " -BUD BUDdir  Write records in a BUD file layout\n"
             " -SDS SDSdir  Write records in a SDS file layout\n"
             " -CSS CSSdir  Write records in a CSS-like file layout\n"
             "\n"
             "The archive 'format' argument is expanded for each record using the\n"
             "following flags:\n"
             "\n"
             "  n : network code, white space removed\n"
             "  s : station code, white space removed\n"
             "  l : location code, white space removed\n"
             "  c : channel code, white space removed\n"
             "  Y : year, 4 digits\n"
             "  y : year, 2 digits zero padded\n"
             "  j : day of year, 3 digits zero padded\n"
             "  H : hour, 2 digits zero padded\n"
             "  M : minute, 2 digits zero padded\n"
             "  S : second, 2 digits zero padded\n"
             "  F : fractional seconds, 4 digits zero padded\n"
             "  q : single character record quality indicator (D, R, Q, M)\n"
             "  L : data record length in bytes\n"
             "  r : Sample rate (Hz) as a rounded integer\n"
             "  R : Sample rate (Hz) as a float with 6 digit precision\n"
             "  %% : the percent (%%) character\n"
             "  # : the number (#) character\n"
             "\n"
             "The flags are prefaced with either the %% or # modifier.  The %% modifier\n"
             "indicates a defining flag while the # indicates a non-defining flag.\n"
             "All records with the same set of defining flags will be written to the\n"
             "same file. Non-defining flags will be expanded using the values in the\n"
             "first record for the resulting file name.\n"
             "\n");
  }
} /* End of usage() */
