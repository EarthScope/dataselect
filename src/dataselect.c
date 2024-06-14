/***************************************************************************
 * dataselect.c - miniSEED data selection.
 *
 * Opens one or more user specified files, applys filtering criteria
 * and outputs any matched data while time-ordering the data and
 * optionally pruning any overlap (at record or sample level).
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
 * The data structure (using actual structure names):
 *
 * MS3TraceList
 *   |-MS3TraceID
 *   |   |-MS3TraceSeg
 *   |        |-MS3RecordList
 *   |            |-MS3RecordPtr
 *   |            |-MS3RecordPtr
 *   |            |-...
 *   |
 *   |-MS3TraceID
 *   |   |-MS3TraceSeg
 *   |        |-MS3RecordList
 *   |            |-MS3RecordPtr
 *   |            |-MS3RecordPtr
 *   |            |-...
 *   |
 *   |-...
 *
 * The program goes through the following stages:
 *
 * 1) Read all input files constructing a view of contiguous trace
 * segments and the data records that comprise them.
 *
 * There is no relationship between the location of input records in
 * specific files or offsets into files.  In other words, the program
 * will reconstruct the most contiguous, time-ordered data segments
 * possible from all the input records regardless of how they are
 * organized in the input files.  The resulting time-ordering of the
 * data records and contiguous segments is a characteristic of the
 * internal data structures and cannot be disabled.
 *
 * 2) If data pruning (removing overlap data) has been selected the
 * data view will be processed to identify all overlapping data and to
 * mark individual Record structures either for complete removal or
 * for partial record trimming (when pruning at the sample level).
 * When a complete record is pruned from the ouput its record length
 * member will be set to 0 indicating that it is no longer
 * contributing to the segment, a special case understood by
 * downstream functions.  Note that no actual data records are changed
 * in this operation, modification of the data records occurs when the
 * data is written to the new output files.
 *
 * 3) Write all contributing data records in the data list out to the
 * output files.  After each record is read into memory it's associated
 * structure is checked to see if the record needs to be
 * trimmed due to sample level pruning.  Trimming a data record involves
 * unpacking, sample removal and repacking.  After trimming or if no
 * trimming is required the data record is written to the appropriate
 * output file. In this way only the minimal number of records needing
 * modification (trimming) are repacked.
 *
 ***************************************************************************/

/* _ISOC9X_SOURCE needed to get a declaration for llabs on some archs */
#define _ISOC9X_SOURCE

#define __STDC_FORMAT_MACROS
#include <errno.h>
#include <inttypes.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/resource.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <time.h>

#include <libmseed.h>
#include <mseedformat.h>

#include "dsarchive.h"

#define VERSION "4.0.1"
#define PACKAGE "dataselect"

/* Input/output file selection information containers */
typedef struct Filelink_s
{
  char *infilename_raw;   /* Input file name with potential annotation (byte range) */
  char *infilename;       /* Input file name without annotation (byte range) */
  FILE *infp;             /* Input file descriptor */
  struct Filelink_s *next;
} Filelink;

/* Archive output structure definition containers */
typedef struct Archive_s
{
  DataStream datastream;
  struct Archive_s *next;
} Archive;

/* Container for coverage entries used to prune data */
typedef struct TimeRange_s
{
  nstime_t starttime;
  nstime_t endtime;
} TimeRange;

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
  MS3RecordPtr *recptr;
  Filelink *flp;
  int8_t *errflagp;
} WriterData;

static int setselectionlimits (MS3TraceList *mstl);

static int writetraces (MS3TraceList *mstl);
static int trimrecord (MS3RecordPtr *rec, char *recbuf, WriterData *writerdata);
static void writerecord (char *record, int reclen, void *handlerdata);

static int prunetraces (MS3TraceList *mstl);
static int findcoverage (MS3TraceList *mstl, MS3TraceID *targetid,
                         MS3TraceSeg *targetseg, Coverage **ppcoverage);
static int trimtrace (MS3TraceSeg *targetseg, const char *targetsourceid,
                      Coverage *coverage);
static int reconcile_tracetimes (MS3TraceList *mstl);

static void printtracelist (MS3TraceList *mstl, uint8_t details);
static void printwritten (MS3TraceList *mstl);

static int sortrecordlist (MS3RecordList *reclist);
static int recordcmp (MS3RecordPtr *rec1, MS3RecordPtr *rec2);

static int processparam (int argcount, char **argvec);
static char *getoptval (int argcount, char **argvec, int argopt);
static int setofilelimit (int limit);
static int addfile (char *filename);
static int addlistfile (char *filename);
static int addarchive (const char *path, const char *layout);
static void usage (int level);

static int8_t verbose = 0;
static int8_t skipnotdata = 0;    /* Controls skipping of non-miniSEED data */
static int8_t bestversion = 1;    /* Use publication version to retain the "best" data when pruning */
static int8_t prunedata = 0;      /* Prune data: 'r= record level, 's' = sample level, 'e' = edges only */
static uint8_t setpubver = 0;     /* Set publication version/quality indicator on output records */
static double timetol = -1.0;     /* Time tolerance for continuous traces */
static double sampratetol = -1.0; /* Sample rate tolerance for continuous traces */
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

static char *outputfile = NULL;  /* Single output file */
static int8_t outputmode = 0;    /* Mode for single output file: 0=overwrite, 1=append */
static Archive *archiveroot = 0; /* Output file structures */

static char recordbuf[MAXRECLEN]; /* Global record buffer */

static Filelink *filelist = NULL;        /* List of input files */
static Filelink *filelisttail = NULL;    /* Tail of list of input files */
static MS3Selections *selections = NULL; /* Data selection criteria, SIDs and time ranges */

static char *writtenfile = NULL;       /* File to write summary of output records */
static char *writtenprefix = NULL;     /* Prefix for summary of output records */
static MS3TraceList *writtentl = NULL; /* TraceList of output records */

int
main (int argc, char **argv)
{
  Filelink *flp;
  MS3TraceList *mstl = NULL;

  uint32_t flags = 0;
  int totalfiles = 0;
  int retcode;

  /* Set default error message prefix */
  ms_loginit (NULL, NULL, NULL, "ERROR: ");

  /* Process input parameters */
  if (processparam (argc, argv) < 0)
    return 1;

  /* Read leap second list file if env. var. LIBMSEED_LEAPSECOND_FILE is set */
  ms_readleapseconds ("LIBMSEED_LEAPSECOND_FILE");

  /* Data stream archiving maximum concurrent open files */
  if (archiveroot)
    ds_maxopenfiles = 50;

  /* Initialize written MS3TraceList */
  if (writtenfile)
    if ((writtentl = mstl3_init (writtentl)) == NULL)
      return 1;

  /* Set flags to:
   * - validate CRCs (if present)
   * - extract start-stop range from file names
   * - construct a record-list for each segment */
  flags |= MSF_VALIDATECRC;
  flags |= MSF_PNAMERANGE;
  flags |= MSF_RECORDLIST;

  if (skipnotdata)
    flags |= MSF_SKIPNOTDATA;

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

    /* Read all miniSEED into a trace list, limiting to selections */
    retcode = ms3_readtracelist_selection (&mstl, flp->infilename_raw, &tolerance,
                                           selections, bestversion, flags, verbose);

    /* Critical error if file was not read properly */
    if (retcode != MS_NOERROR)
    {
      ms_log (2, "Cannot read %s: %s\n", flp->infilename, ms_errorstr (retcode));
      return -1;
    }

    totalfiles++;
    flp = flp->next;
  } /* End of looping over file list */

  /* Increase open file limit if necessary, in general we need the
   * filecount + ds_maxopenfiles and some wiggle room. */
  setofilelimit (totalfiles + ds_maxopenfiles + 20);

  /* Set time limits based on selections when pruning to specific time limits */
  if ((prunedata == 's' || prunedata == 'e') &&
      selections && setselectionlimits (mstl))
    return 1;

  if (verbose > 2)
  {
    ms_log (1, "== Input data ==\n");
    printtracelist (mstl, 1);
  }

  if (mstl->numtraceids == 0)
  {
    if (verbose)
      ms_log (1, "No data selected\n");

    return 0;
  }

  /* Prune data */
  if (prunedata)
  {
    /* Prune overlaps */
    if (prunedata == 'r' || prunedata == 's')
      if (prunetraces (mstl))
        return 1;

    /* Reconcile MS3TraceID times with associated record times */
    if (reconcile_tracetimes (mstl))
      return 1;
  }

  if (verbose > 2)
  {
    ms_log (1, "== Pruned data ==\n");
    printtracelist (mstl, 1);
  }

  /* Write all MS3TraceSeg associated records to output file(s) */
  if (writetraces (mstl))
    return 1;

  if (writtenfile)
  {
    printwritten (writtentl);
    mstl3_free (&writtentl, 1);
  }

  /* The main MS3TraceList (mstl) is not freed on purpose: the structure has a
   * potentially huge number of sub-structures which would take a long time to
   * iterate through.  This would be a waste of time given the program is now done.
   *
   * This may show up as a memory leak for some profilers. */

  return 0;
} /* End of main() */

/***************************************************************************
 * Determine selection limits for each record based on all
 * matching selection entries.
 *
 * At this point data selection has already been performed at the record
 * level by the libmseed logic.  This routine will set new record start
 * and end times when they intersect the record coverage.
 *
 * Return 0 on success and -1 on error.
 ***************************************************************************/
static int
setselectionlimits (MS3TraceList *mstl)
{
  const MS3Selections *select = NULL;
  const MS3SelectTime *selecttime = NULL;
  MS3TraceID *id = NULL;
  MS3TraceSeg *seg = NULL;
  MS3RecordPtr *recptr = NULL;
  TimeRange *timerange = NULL;
  nstime_t newstart;
  nstime_t newend;

  if (!mstl)
    return -1;

  /* Set new record times based on selection times */
  id = mstl->traces.next[0];
  while (id)
  {
    seg = id->first;
    while (seg)
    {
      recptr = seg->recordlist->first;
      while (recptr)
      {
        select = selections;
        while ((select = ms3_matchselect (select,
                                          recptr->msr->sid,
                                          recptr->msr->starttime,
                                          recptr->endtime,
                                          0,
                                          &selecttime)))
        {
          while (selecttime)
          {
            /* Records are either completely or partially selected by time limits */
            newstart = NSTUNSET;
            newend = NSTUNSET;

            if (selecttime->starttime != NSTUNSET &&
                selecttime->starttime > recptr->msr->starttime &&
                selecttime->starttime < recptr->endtime)
            {
              newstart = selecttime->starttime;
            }

            if (selecttime->endtime != NSTUNSET &&
                selecttime->endtime > recptr->msr->starttime &&
                selecttime->endtime < recptr->endtime)
            {
              newend = selecttime->endtime;
            }

            if (newstart == NSTUNSET && newend == NSTUNSET)
            {
              selecttime = selecttime->next;
              continue;
            }

            /* Allocate TimeRange for new time boundaries */
            if (recptr->prvtptr == NULL)
            {
              if ((recptr->prvtptr = (TimeRange *)malloc (sizeof (TimeRange))) == NULL)
              {
                ms_log (2, "%s(): Cannot allocate memory\n", __func__);
                return -1;
              }

              ((TimeRange *)recptr->prvtptr)->starttime = NSTUNSET;
              ((TimeRange *)recptr->prvtptr)->endtime = NSTUNSET;
            }

            timerange = (TimeRange *)recptr->prvtptr;

            if (newstart != NSTUNSET &&
                (timerange->starttime == NSTUNSET || newstart < timerange->starttime))
              timerange->starttime = newstart;

            if (newend != NSTUNSET &&
                (timerange->endtime == NSTUNSET || newend > timerange->endtime))
              timerange->endtime = newend;

            selecttime = selecttime->next;
          }
          select = select->next;
        }
        recptr = recptr->next;
      }
      seg = seg->next;
    }
    id = id->next[0];
  }

  return 0;
} /* End of setselectionlimits() */

/***************************************************************************
 * Write all MS3TraceSeg associated records to output file(s).  If an
 * output file is specified all records will be written to it,
 * otherwise records will be written to specified archive layouts.
 *
 * This routine will also call trimrecord() to trim a record when new
 * start and end times have been identified in earlier processing.
 * Record trimming is triggered when the RecordPtr.prvtptr has new
 * TimeRange.starttime or TimeRange.endtime values.
 *
 * Returns 0 on success and 1 on error.
 ***************************************************************************/
static int
writetraces (MS3TraceList *mstl)
{
  static uint64_t totalrecsout = 0;
  static uint64_t totalbytesout = 0;
  char *wb = "wb";
  char *ab = "ab";
  char *mode;
  int8_t errflag = 0;
  int rv;

  MS3TraceID *id;
  MS3TraceID *groupid;
  MS3TraceSeg *seg;
  MS3RecordPtr *recptr;
  MS3RecordPtr *recptrprev;
  MS3RecordPtr *recptrnext;

  MS3RecordList *groupreclist = NULL;

  TimeRange *newrange;
  Filelink *flpsearch;
  Filelink *flp;

  FILE *ofp = NULL;
  WriterData writerdata;

  writerdata.errflagp = &errflag;

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

  /* Re-link records into write lists, from per-segment lists to per-ID lists.
   * This allows (later) sorting of data records as logical groups regardless
   * from which segment the record was originally associated. */
  id = mstl->traces.next[0];
  groupid = id;
  while (id)
  {
    /* Check if new group ID is needed */
    if (groupid != id && strcmp (groupid->sid, id->sid) != 0)
    {
      groupid = id;
    }

    if (groupid->prvtptr == NULL)
    {
      /* Allocate MS3RecordList for ID-level list */
      if ((id->prvtptr = (MS3RecordList *)malloc (sizeof (MS3RecordList))) == NULL)
      {
        ms_log (2, "%s(): Cannot allocate memory\n", __func__);
        return 1;
      }

      groupreclist = (MS3RecordList *)id->prvtptr;
      groupreclist->first = NULL;
      groupreclist->last = NULL;
      groupreclist->recordcnt = 0;
    }

    seg = id->first;
    while (seg)
    {
      /* Remove non-contributing records from list denoted with 0 reclen */
      if (prunedata)
      {
        recptr = seg->recordlist->first;
        recptrprev = NULL;
        while (recptr)
        {
          recptrnext = recptr->next;

          /* Re-link list to remove recptr, maintaining first and last */
          if (recptr->msr->reclen == 0)
          {
            if (recptr == seg->recordlist->first)
              seg->recordlist->first = recptr->next;
            else if (recptrprev)
              recptrprev->next = recptr->next;

            if (recptr == seg->recordlist->last)
              seg->recordlist->last = recptrprev;

            msr3_free (&recptr->msr);
            free (recptr);
            seg->recordlist->recordcnt--;
          }
          else
          {
            recptrprev = recptr;
          }

          recptr = recptrnext;
        }
      }

      /* Append record list to ID-level list */
      if (seg->recordlist->first != NULL)
      {
        if (groupreclist && groupreclist->first == NULL)
        {
          groupreclist->first = seg->recordlist->first;
          groupreclist->last = seg->recordlist->last;
          groupreclist->recordcnt = seg->recordlist->recordcnt;
        }
        else
        {
          groupreclist->last->next = seg->recordlist->first;
          groupreclist->last = seg->recordlist->last;
          groupreclist->recordcnt += seg->recordlist->recordcnt;
        }
      }

      seg->recordlist->first = NULL;
      seg->recordlist->last = NULL;
      seg->recordlist->recordcnt = 0;

      seg = seg->next;
    }

    id = id->next[0];
  } /* Done combining pruned records into SourceID groups */

  /* Loop through MS3TraceList and write records */
  id = mstl->traces.next[0];
  while (id && errflag == 0)
  {
    groupreclist = (MS3RecordList *)id->prvtptr;

    if (groupreclist && groupreclist->recordcnt > 0)
    {
      /* Sort record list if overlaps have been pruned, if the data has not been
       * pruned it is already in time order. */
      if (prunedata == 'r' || prunedata == 's')
      {
        sortrecordlist (groupreclist);
      }

      /* Write each record.
       * After records are read from the input files, perform any
       * pre-identified pruning before writing data. */
      recptr = groupreclist->first;
      while (recptr && errflag == 0)
      {
        if ((size_t)recptr->msr->reclen > sizeof (recordbuf))
        {
          ms_log (2, "Record length (%d bytes) larger than buffer (%llu bytes)\n",
                  recptr->msr->reclen, (long long unsigned int)sizeof (recordbuf));
          errflag = 1;
          break;
        }

        /* Find the matching input file entry */
        flp = NULL;
        flpsearch = filelist;
        while (flpsearch)
        {
          if (flpsearch->infilename_raw == recptr->filename)
          {
            flp = flpsearch;
            break;
          }

          flpsearch = flpsearch->next;
        }

        if (flp == NULL)
        {
          ms_log (2, "Cannot find input file entry for %s\n", recptr->filename);
          errflag = 1;
          break;
        }

        /* Open file for reading if not already done */
        if (!flp->infp)
          if (!(flp->infp = fopen (flp->infilename, "rb")))
          {
            ms_log (2, "Cannot open '%s' for reading: %s\n",
                    flp->infilename, strerror (errno));
            errflag = 1;
            break;
          }

        /* Seek to record offset */
        if (lmp_fseek64 (flp->infp, recptr->fileoffset, SEEK_SET) == -1)
        {
          ms_log (2, "Cannot seek in '%s': %s\n",
                  flp->infilename, strerror (errno));
          errflag = 1;
          break;
        }

        /* Read record into buffer */
        if (fread (recordbuf, recptr->msr->reclen, 1, flp->infp) != 1)
        {
          ms_log (2, "Cannot read %d bytes at offset %llu from '%s'\n",
                  recptr->msr->reclen, (long long unsigned)recptr->fileoffset,
                  flp->infilename);
          errflag = 1;
          break;
        }

        /* Setup writer data */
        writerdata.ofp = ofp;
        writerdata.recptr = recptr;
        writerdata.flp = flp;

        /* Write out the data, either the record needs to be trimmed (and will be
         * send to the record writer) or we send it directly to the record writer. */
        newrange = (TimeRange *)(recptr->prvtptr);

        /* Trim data from the record if new start or end times are specifed */
        if (newrange && (newrange->starttime != NSTUNSET || newrange->endtime != NSTUNSET))
        {
          rv = trimrecord (recptr, recordbuf, &writerdata);

          if (rv == -1)
          {
            recptr = recptr->next;
            continue;
          }
          if (rv == -2)
          {
            ms_log (1, "Cannot unpack miniSEED from byte offset %" PRId64 " in %s\n",
                    recptr->fileoffset, flp->infilename);
            ms_log (1, "  Writing %s record without trimming\n", id->sid);

            writerecord (recordbuf, recptr->msr->reclen, &writerdata);
          }
        }
        else
        {
          writerecord (recordbuf, recptr->msr->reclen, &writerdata);
        }

        if (errflag)
          break;

        totalrecsout++;
        totalbytesout += recptr->msr->reclen;

        recptr = recptr->next;
      } /* Done looping through record list */
    }

    id = id->next[0];
  } /* Done looping through MS3TraceIDs */

  /* Close all open input & output files and remove backups if requested */
  flp = filelist;
  while (flp)
  {
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
    ms_log (1, "Wrote %" PRIu64 " bytes of %" PRIu64 " records to output file(s)\n",
            totalbytesout, totalrecsout);
  }

  return (errflag) ? 1 : 0;
} /* End of writetraces() */

/***************************************************************************
 * Unpack a data record and trim samples, either from the beginning or
 * the end, to fit the TimeRange.starttime and TimeRange.endtime boundary
 * times and pack the record.
 *
 * Data samples times are not modified.  The new start and end times
 * are treated as arbitrary boundaries, not as explicit new start/end
 * times, this routine calculates which samples fit within the new
 * boundaries.
 *
 * Return 0 on success, -1 on failure or skip and -2 on unpacking errors.
 ***************************************************************************/
static int
trimrecord (MS3RecordPtr *recptr, char *recordbuf, WriterData *writerdata)
{
  nstime_t nsperiod;
  nstime_t ostarttime;
  TimeRange *newrange;

  char stime[32] = {0};
  char etime[32] = {0};

  int trimsamples;
  uint8_t samplesize;
  char sampletype;
  int64_t packedsamples;
  int packedrecords;
  int retcode;

  if (!recptr || !recordbuf)
    return -1;

  ostarttime = recptr->msr->starttime;
  newrange = (TimeRange *)(recptr->prvtptr);

  /* Sanity check for new start/end times */
  if ((newrange->starttime != NSTUNSET && newrange->endtime != NSTUNSET && newrange->starttime > newrange->endtime) ||
      (newrange->starttime != NSTUNSET && (newrange->starttime < recptr->msr->starttime || newrange->starttime > recptr->endtime)) ||
      (newrange->endtime != NSTUNSET && (newrange->endtime > recptr->endtime || newrange->endtime < recptr->msr->starttime)))
  {
    ms_log (2, "Problem with new start/end record bound times.\n");
    ms_log (2, "  Original record %s from %s (byte offset: %llu)\n",
            "SourceID", writerdata->flp->infilename, (unsigned long long)recptr->fileoffset);
    ms_nstime2timestr (recptr->msr->starttime, stime, ISOMONTHDAY_Z, NANO_MICRO);
    ms_nstime2timestr (recptr->endtime, etime, ISOMONTHDAY_Z, NANO_MICRO);
    ms_log (2, "       Start: %s       End: %s\n", stime, etime);
    if (newrange->starttime == NSTUNSET)
      strcpy (stime, "NONE");
    else
      ms_nstime2timestr (newrange->starttime, stime, ISOMONTHDAY_Z, NANO_MICRO);
    if (newrange->endtime == NSTUNSET)
      strcpy (etime, "NONE");
    else
      ms_nstime2timestr (newrange->endtime, etime, ISOMONTHDAY_Z, NANO_MICRO);
    ms_log (2, " Start bound: %-24s End bound: %-24s\n", stime, etime);

    return -1;
  }

  if (ms_encoding_sizetype (recptr->msr->encoding, &samplesize, &sampletype))
  {
    ms_log (2, "Cannot determine sample size and type for encoding %d\n", recptr->msr->encoding);
    return -1;
  }

  /* Check for supported samples types, can only trim what can be packed */
  if (sampletype != 'i' && sampletype != 'f' && sampletype != 'd')
  {
    if (verbose)
    {
      ms_nstime2timestr (recptr->msr->starttime, stime, ISOMONTHDAY_Z, NANO_MICRO);
      ms_log (1, "Skipping trim of %s (%s), unsupported encoding (%d: %s)\n",
              recptr->msr->sid, stime, recptr->msr->encoding, ms_encodingstr (recptr->msr->encoding));
    }

    return 0;
  }

  /* Decode data samples */
  recptr->msr->record = recordbuf;
  if ((retcode = msr3_unpack_data (recptr->msr, 0)) < 0)
  {
    ms_log (2, "Cannot unpack miniSEED record: %s\n", ms_errorstr (retcode));

    return -2;
  }

  if (verbose > 1)
  {
    ms_log (1, "Triming record: %s (%u)\n", recptr->msr->sid, recptr->msr->pubversion);
    ms_nstime2timestr (recptr->msr->starttime, stime, ISOMONTHDAY_Z, NANO_MICRO);
    ms_nstime2timestr (recptr->endtime, etime, ISOMONTHDAY_Z, NANO_MICRO);
    ms_log (1, "       Start: %s        End: %s\n", stime, etime);
    if (newrange->starttime == NSTUNSET)
      strcpy (stime, "NONE");
    else
      ms_nstime2timestr (newrange->starttime, stime, ISOMONTHDAY_Z, NANO_MICRO);
    if (newrange->endtime == NSTUNSET)
      strcpy (etime, "NONE");
    else
      ms_nstime2timestr (newrange->endtime, etime, ISOMONTHDAY_Z, NANO_MICRO);
    ms_log (1, " Start bound: %-24s  End bound: %-24s\n", stime, etime);
  }

  /* Determine sample period in nanosecond time ticks */
  nsperiod = msr3_nsperiod(recptr->msr);

  /* Remove samples from the beginning of the record */
  if (newrange->starttime != NSTUNSET && nsperiod)
  {
    nstime_t newstarttime;

    /* Determine new start time and the number of samples to trim */
    trimsamples = 0;
    newstarttime = recptr->msr->starttime;

    while (newstarttime < newrange->starttime && trimsamples < recptr->msr->samplecnt)
    {
      newstarttime += nsperiod;
      trimsamples++;
    }

    if (trimsamples >= recptr->msr->samplecnt)
    {
      if (verbose > 1)
        ms_log (1, "All samples would be trimmed from record, skipping\n");

      return -1;
    }

    if (verbose > 2)
    {
      ms_nstime2timestr (newstarttime, stime, ISOMONTHDAY_Z, NANO_MICRO);
      ms_log (1, "Removing %d samples from the start, new start time: %s\n", trimsamples, stime);
    }

    memmove (recptr->msr->datasamples,
             (char *)recptr->msr->datasamples + (samplesize * trimsamples),
             samplesize * (recptr->msr->numsamples - trimsamples));

    recptr->msr->numsamples -= trimsamples;
    recptr->msr->samplecnt -= trimsamples;
    recptr->msr->starttime = newstarttime;
    newrange->starttime = newstarttime;
  }

  /* Remove samples from the end of the record */
  if (newrange->endtime != NSTUNSET && nsperiod)
  {
    nstime_t newendtime;

    /* Determine new end time and the number of samples to trim */
    trimsamples = 0;
    newendtime = recptr->endtime;

    while (newendtime > newrange->endtime && trimsamples < recptr->msr->samplecnt)
    {
      newendtime -= nsperiod;
      trimsamples++;
    }

    if (trimsamples >= recptr->msr->samplecnt)
    {
      if (verbose > 1)
        ms_log (1, "All samples would be trimmed from record, skipping\n");

      return -1;
    }

    if (verbose > 2)
    {
      ms_nstime2timestr (newendtime, etime, ISOMONTHDAY_Z, NANO_MICRO);
      ms_log (1, "Removing %d samples from the end, new end time: %s\n", trimsamples, etime);
    }

    recptr->msr->numsamples -= trimsamples;
    recptr->msr->samplecnt -= trimsamples;
    newrange->endtime = newendtime;
  }

  /* Add the v2 "sequence number" to extra headers so it is included in output */
  if (recptr->msr->formatversion == 2)
  {
    int64_t seqnum = 0;
    char seqstr[7];
    char *endptr;

    memcpy (seqstr, recordbuf, 6);
    seqstr[6] = '\0';

    seqnum = (int64_t)strtoll (seqstr, &endptr, 10);

    if (endptr != seqstr)
    {
      if (mseh_set (recptr->msr, "/FDSN/Sequence", &seqnum, 'i'))
      {
        ms_log (2, "Cannot set sequence number in extra headers\n");
      }
    }
  }

  /* Pack the data record into the global record buffer used by writetraces() */
  packedrecords = msr3_pack (recptr->msr, &writerecord, writerdata,
                             &packedsamples, MSF_FLUSHDATA, verbose - 1);

  if (packedrecords <= 0)
  {
    ms_nstime2timestr (ostarttime, stime, ISOMONTHDAY_Z, NANO_MICRO);
    ms_log (2, "%s(): Cannot pack miniSEED record for %s %s\n",
            __func__, recptr->msr->sid, stime);

    return -2;
  }

  /* Free allocated samples */
  libmseed_memory.free(recptr->msr->datasamples);
  recptr->msr->datasamples = NULL;
  recptr->msr->numsamples = 0;

  return 0;
} /* End of trimrecord() */

/***************************************************************************
 * Used by writetraces() directly, and trimrecord() when called, to save
 * repacked miniSEED to global record buffer.
 ***************************************************************************/
static void
writerecord (char *record, int reclen, void *handlerdata)
{
  WriterData *writerdata = handlerdata;
  Archive *arch;

  if (!record || reclen <= 0 || !handlerdata)
    return;

  /* Set v3 publication version or v2 data quality indicator */
  if (setpubver)
  {
    if (writerdata->recptr->msr->formatversion == 2)
    {
      unsigned char dataquality;

      if (setpubver == 1)
        dataquality = 'R';
      else if (setpubver == 2)
        dataquality = 'D';
      else if (setpubver == 3)
        dataquality = 'Q';
      else
        dataquality = 'M';

      if (verbose > 2)
        ms_log (1, "Setting v2 data quality indicator to '%c'\n", setpubver);

      *pMS2FSDH_DATAQUALITY (recordbuf) = dataquality;
    }
    else if (writerdata->recptr->msr->formatversion == 3)
    {
      if (verbose > 2)
        ms_log (1, "Setting publication version to %u\n", setpubver);

      *pMS3FSDH_PUBVERSION (record) = setpubver;

      /* Recalculate CRC */
      *pMS3FSDH_CRC (record) = 0;
      uint32_t crc = ms_crc32c ((uint8_t *)record, reclen, 0);
      *pMS3FSDH_CRC (record) = HO4u (crc, ms_bigendianhost ());
    }
    else
    {
      ms_log (2, "Cannot set publication version for format version %d\n",
              writerdata->recptr->msr->formatversion);
    }
  }

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
    if (archiveroot)
    {
      arch = archiveroot;
      while (arch)
      {
        if (ds_streamproc (&arch->datastream,
                           writerdata->recptr->msr,
                           verbose - 1, NULL))
        {
          *writerdata->errflagp = 1;
        }

        arch = arch->next;
      }
    }

    if (writtenfile)
    {
      MS3TraceSeg *seg;

      if ((seg = mstl3_addmsr (writtentl, writerdata->recptr->msr, 0, 0, 0, NULL)) == NULL)
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

        *((int64_t *)seg->prvtptr) += reclen;
      }
    }
  }
} /* End of writerecord() */

/***************************************************************************
 * Prune all redundant data from the records list entries associated with
 * the specified MS3TraceSegs.
 *
 * For each MS3TraceSeg determine the coverage of the record list associated
 * with each overlapping, higher-priority MS3TraceSeg using findcoverage().
 * If some higher-priority overlap was determined to exist modify the
 * record list of the MS3TraceSeg in question to mark the overlapping data
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
   * Records from the other traces with a higher priority and prune
   * the overlap. */
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
 * Search an MS3TraceList for entries that overlap the target MS3TraceSeg
 * and, from the record entries of the overlapping MS3TraceSegs, build a
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
 * When no overlap coverage is found *ppcoverage will be NULL, otherwise
 * it will contain a list of representing the higher-priority overlap
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
  MS3RecordPtr *recptr;
  Coverage *coverage = NULL;
  Coverage *prevcoverage = NULL;
  TimeRange *newrange;
  nstime_t nsperiod, nstimetol;
  nstime_t effstarttime, effendtime;
  int priority;
  int newsegment;

  if (!mstl || !targetid || !targetseg || !ppcoverage)
    return -1;

  *ppcoverage = NULL;

  /* Determine sample period in high precision time ticks */
  nsperiod = (targetseg->samprate) ? (nstime_t)(NSTMODULUS / targetseg->samprate + 0.5) : 0;

  /* Determine time tolerance in high precision time ticks */
  nstimetol = (timetol == -1.0) ? (nsperiod / 2) : (nstime_t)(NSTMODULUS * timetol);

  /* Loop through each MS3TraceID in the list */
  id = mstl->traces.next[0];
  while (id)
  {
    /* Continue with next if SourceID is different */
    if (targetid != id)
    {
      if (strcmp (id->sid, targetid->sid))
      {
        id = id->next[0];
        continue;
      }
    }

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

      /* Skip segments with no time coverage (0 samprate) */
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

      /* Check for duplicate or overlap SourceIDs last coverage entry */
      if (coverage)
      {
        /* At this point the SourceID and rate are the same, check if the
         * segment is completly contained by the previous coverage entry. */
        if (seg->starttime >= coverage->starttime &&
            seg->endtime <= coverage->endtime)
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
          /* Loop through list of records, and determine contiguous coverage */
          recptr = seg->recordlist->first;
          newsegment = 1;
          while (recptr)
          {
            /* Check if record has been marked as non-contributing */
            if (recptr->msr->reclen == 0)
            {
              recptr = recptr->next;
              continue;
            }

            newrange = (TimeRange *)recptr->prvtptr;

            /* Determine effective record start and end times */
            effstarttime = (newrange && newrange->starttime != NSTUNSET) ? newrange->starttime : recptr->msr->starttime;
            effendtime = (newrange && newrange->endtime != NSTUNSET) ? newrange->endtime : recptr->endtime;

            /* Create a new segment if a break in the time-series is detected */
            if (coverage)
              if (llabs ((coverage->endtime + nsperiod) - effstarttime) > nstimetol)
                newsegment = 1;

            if (newsegment)
            {
              newsegment = 0;

              prevcoverage = coverage;

              if ((coverage = (Coverage *)malloc (sizeof (Coverage))) == NULL)
              {
                ms_log (2, "Cannot allocate memory for coverage, bah humbug.\n");
                return -1;
              }

              if (*ppcoverage == NULL)
                *ppcoverage = coverage;
              else
                prevcoverage->next = coverage;

              coverage->pubversion = id->pubversion;
              coverage->samprate = seg->samprate;
              coverage->starttime = effstarttime;
              coverage->next = NULL;
            }

            if (coverage)
              coverage->endtime = effendtime;
            else
              ms_log (2, "ACK! covergage is not allocated!?  PLEASE REPORT\n");

            recptr = recptr->next;
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
 * Adjust Record entries associated with the target MS3TraceSeg that
 * are overlapping the time represented by the Coverage
 * in two different ways: 1) mark records that are completely
 * overlapped and 2) determine partial record trim boundaries (new
 * record times) if sample level pruning is requested.
 *
 * Completely overlapping record entries are marked for omission by
 * setting reclen = 0.  Partial Record overlaps are noted by setting
 * newrange.starttime and newrange.newend when sample level pruning
 * is requested.  The actual trimming of the data records, complete or
 * partial, is performed during the output sequence, not in this
 * routine.
 *
 * Returns the number of Record modifications on success and -1 on error.
 ***************************************************************************/
static int
trimtrace (MS3TraceSeg *targetseg, const char *targetsourceid, Coverage *coverage)
{
  MS3RecordPtr *recptr;
  TimeRange *newrange;
  Coverage *cov;
  nstime_t effstarttime, effendtime;
  nstime_t nsperiod, nstimetol;
  char stime[32] = {0};
  char etime[32] = {0};
  int modcount = 0;

  if (!targetseg || !coverage)
    return -1;

  /* Determine sample period in high precision time ticks */
  nsperiod = (targetseg->samprate) ? (nstime_t)(NSTMODULUS / targetseg->samprate + 0.5) : 0;

  /* Determine time tolerance in high precision time ticks */
  nstimetol = (timetol == -1.0) ? (nsperiod / 2) : (nstime_t)(NSTMODULUS * timetol);

  /* Traverse the record list for the target segment and mark records
   * that overlap or intersect with the coverage */
  recptr = targetseg->recordlist->first;
  while (recptr)
  {
    cov = coverage;
    while (cov)
    {
      if (!recptr->msr->reclen) /* Skip if marked non-contributing */
        break;

      newrange = (TimeRange *)recptr->prvtptr;

      /* Determine effective record start and end times for comparison */
      effstarttime = (newrange && newrange->starttime != NSTUNSET) ? newrange->starttime : recptr->msr->starttime;
      effendtime = (newrange && newrange->endtime != NSTUNSET) ? newrange->endtime : recptr->endtime;

      /* Mark record if it is completely overlapped by the coverage including tolerance */
      if (effstarttime >= (cov->starttime - nstimetol) &&
          effendtime <= (cov->endtime + nstimetol))
      {
        if (verbose > 1)
        {
          ms_nstime2timestr (recptr->msr->starttime, stime, ISOMONTHDAY_Z, NANO_MICRO);
          ms_nstime2timestr (recptr->endtime, etime, ISOMONTHDAY_Z, NANO_MICRO);
          ms_log (1, "Removing Record [complete overlap] %s (%u) :: %s  %s  offset: %" PRId64 ", reclen: %d\n",
                  targetsourceid, recptr->msr->pubversion, stime, etime,
                  recptr->fileoffset, recptr->msr->reclen);
        }

        recptr->msr->reclen = 0;
        modcount++;
      }

      /* Determine the new start/end times if pruning at the sample level */
      if (prunedata == 's' && recptr->msr->reclen != 0)
      {
        /* Record intersects beginning of coverage */
        if (effstarttime < cov->starttime &&
            (effendtime + nstimetol) >= cov->starttime)
        {
          if (recptr->prvtptr == NULL)
          {
            if ((recptr->prvtptr = (TimeRange *)malloc (sizeof (TimeRange))) == NULL)
            {
              ms_log (2, "Cannot allocate memory for TimeRange, bah humbug.\n");
              return -1;
            }

            ((TimeRange *)recptr->prvtptr)->starttime = NSTUNSET;
            ((TimeRange *)recptr->prvtptr)->endtime = NSTUNSET;
          }

          newrange = (TimeRange *)recptr->prvtptr;

          /* Set new end time boundary including specified time tolerance */
          newrange->endtime = cov->starttime - nsperiod + nstimetol;

          if (newrange->starttime != NSTUNSET && newrange->endtime < newrange->starttime)
          {
            if (verbose > 1)
            {
              ms_nstime2timestr (recptr->msr->starttime, stime, ISOMONTHDAY_Z, NANO_MICRO);
              ms_nstime2timestr (recptr->endtime, etime, ISOMONTHDAY_Z, NANO_MICRO);
              ms_log (1, "Removing record [start intersect] %s (%u) :: %s  %s\n",
                      targetsourceid, recptr->msr->pubversion, stime, etime);
            }

            recptr->msr->reclen = 0;
            modcount++;
          }
          else
          {
            effendtime = newrange->endtime;
            modcount++;
          }
        }

        /* Record intersects end of coverage */
        if ((effstarttime - nstimetol) <= cov->endtime &&
            effendtime > cov->endtime)
        {
          if (recptr->prvtptr == NULL)
          {
            if ((recptr->prvtptr = (TimeRange *)malloc (sizeof (TimeRange))) == NULL)
            {
              ms_log (2, "Cannot allocate memory for TimeRange, bah humbug.\n");
              return -1;
            }

            ((TimeRange *)recptr->prvtptr)->starttime = NSTUNSET;
            ((TimeRange *)recptr->prvtptr)->endtime = NSTUNSET;
          }

          newrange = (TimeRange *)recptr->prvtptr;

          /* Set Record new start time boundary including specified time tolerance */
          newrange->starttime = cov->endtime + nsperiod - nstimetol;

          if (newrange->endtime != NSTUNSET && newrange->starttime > newrange->endtime)
          {
            if (verbose > 1)
            {
              ms_nstime2timestr (recptr->msr->starttime, stime, ISOMONTHDAY_Z, NANO_MICRO);
              ms_nstime2timestr (recptr->endtime, etime, ISOMONTHDAY_Z, NANO_MICRO);
              ms_log (1, "Removing record [end intersect] %s (%u) :: %s  %s\n",
                      targetsourceid, recptr->msr->pubversion, stime, etime);
            }

            recptr->msr->reclen = 0;
            modcount++;
          }
          else
          {
            effstarttime = newrange->starttime;
            modcount++;
          }
        }

        /* Remove record if all samples have been pruned within tolerance,
         * test for special cases of:
         * a) no time coverage (single sample) and no pruning
         * b) no time coverage (single last sample) and split boundary usage */
        if (effstarttime >= (effendtime - nstimetol) &&
            !(recptr->msr->starttime == recptr->endtime &&
              recptr->msr->starttime == effstarttime &&
              recptr->endtime == effendtime))
        {
          if (verbose > 1)
          {
            ms_nstime2timestr (recptr->msr->starttime, stime, ISOMONTHDAY_Z, NANO_MICRO);
            ms_nstime2timestr (recptr->endtime, etime, ISOMONTHDAY_Z, NANO_MICRO);
            ms_log (1, "Removing record [all pruned] %s (%u) :: %s  %s\n",
                    targetsourceid, recptr->msr->pubversion, stime, etime);
          }

          recptr->msr->reclen = 0;
          modcount++;
        }

      } /* Done pruning at sample level */

      cov = cov->next;
    }

    recptr = recptr->next;
  }

  return modcount;
} /* End of trimtrace() */

/***************************************************************************
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
  MS3RecordPtr *recptr;
  MS3RecordPtr *first = NULL;
  MS3RecordPtr *last = NULL;
  TimeRange *newrange;

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
      /* Find first contributing record (reclen != 0) */
      recptr = seg->recordlist->first;
      while (recptr)
      {
        if (recptr->msr->reclen > 0)
        {
          if (!first)
            first = recptr;

          last = recptr;
        }

        recptr = recptr->next;
      }

      /* Set a new MS3TraceSeg start time */
      if (first)
      {
        newrange = (TimeRange *)first->prvtptr;

        /* Use the new boundary start time if set and sane */
        if (newrange && newrange->starttime != NSTUNSET &&
            newrange->starttime > first->msr->starttime)
          seg->starttime = newrange->starttime;
        /* Otherwise use the record start time */
        else
          seg->starttime = first->msr->starttime;
      }

      /* Set a new MS3TraceSeg end time */
      if (last)
      {
        newrange = (TimeRange *)last->prvtptr;

        /* Use the new boundary end time if set and sane */
        if (newrange && newrange->endtime != NSTUNSET &&
            newrange->endtime < last->endtime)
          seg->endtime = newrange->endtime;
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
 * Print record list for each MS3TraceSeg to stdout.
 ***************************************************************************/
static void
printtracelist (MS3TraceList *mstl, uint8_t details)
{
  MS3TraceID *id = NULL;
  MS3TraceSeg *seg = NULL;
  MS3RecordPtr *recptr = NULL;
  TimeRange *newrange = NULL;
  char stime[32] = {0};
  char etime[32] = {0};
  int segcnt = 0;

  if (!mstl)
    return;

  /* Print out the appropriate header */
  ms_log (0, "   Source              Start sample             End sample        Hz   Samples\n");

  id = mstl->traces.next[0];
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

      if (!seg->recordlist)
      {
        ms_log (2, "No record list associated with this MS3TraceSeg.\n");
      }
      else
      {
        recptr = seg->recordlist->first;
        while (recptr)
        {
          ms_log (0, "  Filename: %s  Offset: %" PRId64 "  RecLen: %d  PubVersion: %u\n",
                  (recptr->filename) ? recptr->filename : "NONE", recptr->fileoffset,
                  recptr->msr->reclen, recptr->msr->pubversion);

          ms_nstime2timestr (recptr->msr->starttime, stime, ISOMONTHDAY_Z, NANO_MICRO);
          ms_nstime2timestr (recptr->endtime, etime, ISOMONTHDAY_Z, NANO_MICRO);
          ms_log (0, "        Start: %s        End: %s\n", stime, etime);

          if (details && recptr->prvtptr != NULL)
          {
            newrange = (TimeRange *)recptr->prvtptr;

            if (newrange->starttime == NSTUNSET)
              strcpy (stime, "NONE");
            else
              ms_nstime2timestr (newrange->starttime, stime, ISOMONTHDAY_Z, NANO_MICRO);
            if (newrange->endtime == NSTUNSET)
              strcpy (etime, "NONE");
            else
              ms_nstime2timestr (newrange->endtime, etime, ISOMONTHDAY_Z, NANO_MICRO);

            ms_log (0, " Select start: %-24s Select end: %-24s\n", stime, etime);
          }

          recptr = recptr->next;
        }
      }

      segcnt++;
      seg = seg->next;
    }

    id = id->next[0];
  }

  ms_log (0, "End of trace list: %d trace segment(s)\n\n", segcnt);

} /* End of printtracelist() */

/***************************************************************************
 * Print summary of output records.
 ***************************************************************************/
static void
printwritten (MS3TraceList *mstl)
{
  MS3TraceID *id = NULL;
  MS3TraceSeg *seg = NULL;
  char stime[32] = {0};
  char etime[32] = {0};
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

      fprintf (ofp, "%s%s|%u|%s|%s|%" PRId64 "|%" PRId64 "\n",
               (writtenprefix) ? writtenprefix : "",
               id->sid, id->pubversion, stime, etime,
               *((int64_t *)seg->prvtptr),
               seg->samplecnt);

      seg = seg->next;
    }

    id = id->next[0];
  }

  if (ofp != stdout && fclose (ofp))
    ms_log (2, "Cannot close output file: %s (%s)\n",
            writtenfile, strerror (errno));

} /* End of printwritten() */

/***************************************************************************
 * Sort a record list so that records are in time order using a
 * mergesort algorithm.
 *
 * The mergesort implementation was inspired by the listsort function
 * published and copyright 2001 by Simon Tatham.
 *
 * Return 0 on success and -1 on error.
 ***************************************************************************/
static int
sortrecordlist (MS3RecordList *reclist)
{
  MS3RecordPtr *p, *q, *e, *top, *tail;
  int nmerges;
  int insize, psize, qsize, i;

  if (reclist == NULL)
    return -1;

  /* Done if no records in list */
  if (reclist->recordcnt == 0)
    return 0;

  top = reclist->first;
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
      reclist->first = top;
      reclist->last = tail;

      return 0;
    }

    /* Otherwise repeat, merging lists twice the size */
    insize *= 2;
  }
} /* End of sortrecordlist() */

/***************************************************************************
 * Compare the start times of each Record for the purposes of sorting
 * a record list.
 *
 * Return 1 if rec1 is "greater" than rec2, otherwise return 0.
 ***************************************************************************/
static int
recordcmp (MS3RecordPtr *rec1, MS3RecordPtr *rec2)
{
  TimeRange *newrange1;
  TimeRange *newrange2;
  nstime_t start1;
  nstime_t start2;

  if (!rec1 || !rec2)
    return -1;

  /* Determine effective start times */
  newrange1 = (TimeRange *)rec1->prvtptr;
  start1 = (newrange1 && newrange1->starttime != NSTUNSET) ? newrange1->starttime : rec1->msr->starttime;

  newrange2 = (TimeRange *)rec2->prvtptr;
  start2 = (newrange2 && newrange2->starttime != NSTUNSET) ? newrange2->starttime : rec2->msr->starttime;

  if (start1 > start2)
  {
    return 1;
  }

  return 0;
} /* End of recordcmp() */

/***************************************************************************
 * Process the command line parameters.
 *
 * Returns 0 on success, and -1 on failure
 ***************************************************************************/
static int
processparam (int argcount, char **argvec)
{
  nstime_t timestart = NSTUNSET;
  nstime_t timeend = NSTUNSET;
  char matchpattern[100] = {0};
  char *selectfile = NULL;
  char *tptr = NULL;
  char *endptr = NULL;
  unsigned long ulong;
  int optind;

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
    else if (strcmp (argvec[optind], "-snd") == 0)
    {
      skipnotdata = 1;
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
      timestart = ms_timestr2nstime (getoptval (argcount, argvec, optind++));
      if (timestart == NSTERROR)
        return -1;
    }
    else if (strcmp (argvec[optind], "-te") == 0)
    {
      timeend = ms_timestr2nstime (getoptval (argcount, argvec, optind++));
      if (timeend == NSTERROR)
        return -1;
    }
    else if (strcmp (argvec[optind], "-M") == 0)
    {
      /* Accept value if it is valid globbbing characters an FDSN SourceID */
      tptr = getoptval (argcount, argvec, optind++);
      if (strspn (tptr, "-[]*?:_ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrtsuvwxyz0123456789") == strlen (tptr))
      {
        strncpy (matchpattern, tptr, sizeof (matchpattern) - 1);
      }
      else
      {
        ms_log (2, "Invalid globbing pattern: %s\n", tptr);
        ms_log (2, "Regular expressions are no longer supported, see the -m option\n");
        return -1;
      }
    }
    else if (strcmp (argvec[optind], "-m") == 0)
    {
      strncpy (matchpattern, getoptval (argcount, argvec, optind++), sizeof (matchpattern) - 1);
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
    else if (strcmp (argvec[optind], "-Q") == 0)
    {
      tptr = getoptval (argcount, argvec, optind++);

      if (tptr[0] == 'R' && tptr[1] == '\0')
        setpubver = 1;
      else if (tptr[0] == 'D' && tptr[1] == '\0')
        setpubver = 2;
      else if (tptr[0] == 'Q' && tptr[1] == '\0')
        setpubver = 3;
      else if (tptr[0] == 'M' && tptr[1] == '\0')
        setpubver = 4;
      else
      {
        ulong = strtoul (tptr, &endptr, 10);

        if (*endptr == '\0' && ulong > 0 && ulong <= UINT8_MAX)
        {
          setpubver = ulong;
        }
        else
        {
          ms_log (2, "Invalid publication version/quality indicator: %s\n", tptr);
          return -1;
        }
      }
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
    else if (strcmp (argvec[optind], "-VCHAN") == 0)
    {
      if (addarchive (getoptval (argcount, argvec, optind++), VCHANLAYOUT) == -1)
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

  /* Combine SourceID match pattern, time start and end into a selection entry */
  if (matchpattern[0] || timestart != NSTUNSET || timeend != NSTUNSET)
  {
    size_t mplength = strlen (matchpattern);

    /* Add wildcards to match pattern for logical "contains" */
    if (matchpattern[0] && mplength < (sizeof (matchpattern) - 3))
    {
      memmove (matchpattern + 1, matchpattern, mplength);
      matchpattern[0] = '*';
      matchpattern[mplength + 1] = '*';
      matchpattern[mplength + 2] = '\0';
    }
    else if (matchpattern[0] == 0)
    {
      matchpattern[0] = '*';
      matchpattern[1] = '\0';
    }

    if (ms3_addselect (&selections, matchpattern, timestart, timeend, 0))
    {
      ms_log (2, "Unable to add selection: '%s'\n", tptr);
      return -1;
    }
  }

  /* Report the program version */
  if (verbose)
    ms_log (1, "%s version: %s\n", PACKAGE, VERSION);

  return 0;
} /* End of processparam() */

/***************************************************************************
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

  if (rlim.rlim_cur < (rlim_t)limit)
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
           " -snd         Skip non-miniSEED data, otherwise quit on unrecognized input\n"
           " -E           Consider all qualities equal instead of 'best' prioritization\n"
           "\n"
           " ## Data selection options ##\n"
           " -s file      Specify a file containing selection criteria\n"
           " -ts time     Limit to records that contain or start after time\n"
           " -te time     Limit to records that contain or end before time\n"
           "                time format: 'YYYY-MM-DD[THH:MM:SS.FFFFFFFFF]'\n"
           " -m match     Limit to records containing the specified pattern\n"
           "                Patterns are applied to: 'FDSN:NET_STA_LOC_BAND_SOURCE_SS'\n"
           "\n"
           " ## Output options ##\n"
           " -o file      Specify a single output file, use +o file to append\n"
           " -A format    Write all records in a custom directory/file layout (try -H)\n"
           " -Pr          Prune data at the record level using 'best' version priority\n"
           " -Ps          Prune data at the sample level using 'best' version priority\n"
           " -Pe          Prune traces at user specified edges only, leave overlaps\n"
           " -Q #DRQM     Specify publication version of all output records\n"
           "\n"
           " ## Logging ##\n"
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
             " -VCHAN dir   Write records into separate Net.Sta.Loc.Chan.PubVersion files\n"
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
             "  v : publication version, 1-255\n"
             "  q : data quality if possible, otherwise pub version (D, R, Q, M, or #)\n"
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
