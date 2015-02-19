/***************************************************************************
 * dataselect.c - Mini-SEED data selection.
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
 * Written by Chad Trabant, IRIS Data Management Center.
 *
 * modified 2015.050
 ***************************************************************************/

/***************************************************************************
 *
 * Data structures and operational overview
 *
 * The data map (using actual structure names):
 *
 * MSTraceList
 *   |-MSTraceID
 *   |   |-MSTraceSeg
 *   |        |-RecordMap
 *   |            |-Record
 *   |            |-Record
 *   |            |-...
 *   |
 *   |-MSTraceID
 *   |   |-MSTraceSeg
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
 * read, search the MSTraceList for an MSTraceSeg that the record "fits"
 * with (same channel and adjacent in time).  If the record is found
 * to fit with a specific MSTraceSeg its coverage will be added to the
 * MSTraceSeg information otherwise a new MSTraceSeg is created and
 * added to the MSTraceList.
 *
 * Each MSTraceSeg has an associated RecordMap which includes a list of
 * Records.  A Record structure is not the actual data record itself
 * but meta information about the record coverage and location (file
 * and offset).  The list of Records is always in time order.
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
 * output files.  If the original files are to be replaced they will
 * be renamed by adding a ".orig" prefix before the writing process
 * begins.  After each record is read into memory it's associated
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
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/resource.h>
#include <errno.h>
#include <time.h>
#include <regex.h>

#include <libmseed.h>

#include "dsarchive.h"

#define VERSION "3.17"
#define PACKAGE "dataselect"

/* Input/output file selection information containers */
typedef struct Filelink_s {
  char    *infilename;     /* Input file name */
  FILE    *infp;           /* Input file descriptor */
  char    *outfilename;    /* Output file name */
  FILE    *outfp;          /* Output file descriptor */
  uint64_t startoffset;    /* Byte offset to start reading, 0 = unused */
  uint64_t endoffset;      /* Byte offset to end reading, 0 = unused */
  int      reordercount;   /* Number of records re-ordered */
  int      recsplitcount;  /* Number of records split */
  int      recrmcount;     /* Number of records removed */
  int      rectrimcount;   /* Number of records trimmed */
  hptime_t earliest;       /* Earliest data time in this file selection */
  hptime_t latest;         /* Latest data time in this file selection */
  int      byteswritten;   /* Number of bytes written out */
  struct Filelink_s *next;
} Filelink;

/* Archive output structure definition containers */
typedef struct Archive_s {
  DataStream  datastream;
  struct Archive_s *next;
} Archive;

/* Mini-SEED record information structures */
typedef struct Record_s {
  struct Filelink_s *flp;
  off_t     offset;
  off_t     stageoffset;
  int       reclen;
  hptime_t  starttime;
  hptime_t  endtime;
  char      quality;
  hptime_t  selectstart;
  hptime_t  selectend;
  hptime_t  newstart;
  hptime_t  newend;
  struct Record_s *prev;
  struct Record_s *next;
} Record;

/* Record map, holds Record structures for a given MSTrace */
typedef struct RecordMap_s {
  long long int    recordcnt;
  struct Record_s *first;
  struct Record_s *last;
} RecordMap;


static int readfiles (MSTraceList **ppmstl);
static int processtraces (MSTraceList *mstl);
static int writetraces (MSTraceList *mstl);
static int trimrecord (Record *rec, char *recbuf);
static void record_handler (char *record, int reclen, void *handlerdata);

static int prunetraces (MSTraceList *mstl);
static int findcoverage (MSTraceList *mstl, MSTraceID *targetid,
			 MSTraceSeg *targetseg, MSTraceGroup **ppcoverage);
static int trimtrace (MSTraceSeg *targetseg, char *targetsrcname,
		      MSTraceGroup *coverage);
static int reconcile_tracetimes (MSTraceList *mstl);
static int qcompare (const char quality1, const char quality2);

static int minsegmentlength (MSTraceList *mstl, double minseconds);
static int longestsegmentonly (MSTraceList *mstl);

static void printmodsummary (flag nomods);
static void printtracemap (MSTraceList *mstl);
static void printrecordmap (RecordMap *recmap, flag details);
static void printwritten (MSTraceList *mstl);

static int findselectlimits (Selections *select, char *srcname,
			     hptime_t starttime, hptime_t endtime, Record *rec);
static int sortrecmap (RecordMap *recmap);
static int recordcmp (Record *rec1, Record *rec2);

static int  processparam (int argcount, char **argvec);
static char *getoptval (int argcount, char **argvec, int argopt);
static int  setofilelimit (int limit);
static int  addfile (char *filename);
static int  addlistfile (char *filename);
static int  addarchive (const char *path, const char *layout);
static int  readregexfile (char *regexfile, char **pppattern);
static off_t calcsize (char *sizestr);
static void usage (int level);


static flag     verbose       = 0;
static flag     basicsum      = 0;    /* Controls printing of basic summary */
static flag     bestquality   = 1;    /* Use D, R, Q, M quality to retain the "best" data when pruning */
static flag     prunedata     = 0;    /* Prune data: 'r= record level, 's' = sample level, 'e' = edges only */
static double   timetol       = -1.0; /* Time tolerance for continuous traces */
static double   sampratetol   = -1.0; /* Sample rate tolerance for continuous traces */
static char     restampqind   = 0;    /* Re-stamp data record/quality indicator */
static int      reclen        = -1;   /* Input data record length, autodetected in most cases */
static flag     modsummary    = 0;    /* Print modification summary after all processing */
static hptime_t starttime     = HPTERROR;  /* Limit to records containing or after starttime */
static hptime_t endtime       = HPTERROR;  /* Limit to records containing or before endtime */
static char     splitboundary = 0;    /* Split records on day 'd', hour 'h' or minute 'm' boundaries */
static char     splitreclen   = 0;    /* Split output files on record length changes */

static regex_t *match         = 0;    /* Compiled match regex */
static regex_t *reject        = 0;    /* Compiled reject regex */

static flag     skipzerosamps = 0;    /* Controls skipping of records with zero samples */
static flag     longestonly   = 0;    /* Controls output of longest segment/channel only */
static double   minseglength  = 0.0;  /* Minimum segment length in seconds */

static flag     replaceinput  = 0;    /* Replace input files */
static flag     nobackups     = 0;    /* Remove re-named original files when done with them */
static char    *outputfile    = 0;    /* Single output file */
static flag     outputmode    = 0;    /* Mode for single output file: 0=overwrite, 1=append */
static Archive *archiveroot   = 0;    /* Output file structures */

static char     recordbuf[16384];     /* Global record buffer */
static char    *stagebuffer   = 0;    /* Global record staging buffer */
static off_t    stagelength   = 0;    /* Global record staging buffer length */
static off_t    stageoffset   = 0;    /* Global record staging buffer end offset */
static flag     stagefull     = 0;    /* Global record staging buffer full flag */

static Filelink *filelist     = 0;    /* List of input files */
static Filelink *filelisttail = 0;    /* Tail of list of input files */
static Selections *selections = 0;    /* List of data selections */

static char *writtenfile      = 0;    /* File to write summary of output records */
static char *writtenprefix    = 0;    /* Prefix for summary of output records */
static MSTraceList *writtentl = 0;    /* TraceList of output records */

int
main ( int argc, char **argv )
{
  MSTraceList *mstl = 0;
  
  /* Set default error message prefix */
  ms_loginit (NULL, NULL, NULL, "ERROR: ");
  
  /* Process input parameters */
  if ( processparam (argc, argv) < 0 )
    return 1;
  
  /* Data stream archiving maximum concurrent open files */
  if ( archiveroot )
    ds_maxopenfiles = 50;
  
  /* Init written MSTraceList */
  if ( writtenfile )
    if ( (writtentl = mstl_init (writtentl)) == NULL )
      return 1;
  
  if ( verbose > 2 )
    ms_log (1, "Processing input files\n");
  
  /* Read and process all files specified on the command line */
  if ( readfiles (&mstl) )
    return 1;

  /* Processes traces */
  if ( mstl->numtraces > 0 && processtraces (mstl) )
    return 1;
  
  if ( modsummary )
    printmodsummary (verbose);
  
  if ( writtenfile )
    {
      printwritten (writtentl);
      mstl_free (&writtentl, 1);
    }
  
  /* The main MSTraceList (mstl) is not freed on purpose: the structure has a
   * potentially huge number of sub-structures (Records in the RecordMap of
   * each ID) which would take a long time iterate through.  This would be a
   * waste of time given that the program is now done.
   * Mind that this may show up as a memory leak for some profilers. */

  return 0;
}  /* End of main() */


/***************************************************************************
 * readfiles:
 *
 * Read input files specified as a Filelink list and populate an
 * MSTraceList and record maps for each trace.  All input files are
 * renamed with a ".orig" suffix before being read.
 *
 * Returns 0 on success and -1 otherwise.
 ***************************************************************************/
static int
readfiles (MSTraceList **ppmstl)
{
  MSFileParam *msfp = NULL;
  Filelink *flp;
  MSRecord *msr = 0;
  MSTraceSeg *seg = 0;
  
  int totalrecs  = 0;
  int totalsamps = 0;
  int totalfiles = 0;
  
  Selections *matchsp = 0;
  SelectTime *matchstp = 0;
  
  RecordMap *recmap = 0;
  Record *rec = 0;
  
  RecordMap newrecmap;
  Record *newrec = 0;
  
  off_t fpos = 0;
  hptime_t recstarttime;
  hptime_t recendtime;
  
  char srcname[50];
  char stime[30];
  
  int infilenamelen = 0;
  int retcode;
  flag whence;

  if ( ! ppmstl )
    return -1;
  
  /* Initialize MSTraceList */
  *ppmstl = mstl_init (*ppmstl);
  
  if ( ! *ppmstl )
    {
      ms_log (2, "readfiles(): cannot (re)initialize MSTraceList\n");
      return -1;
    }
  
  /* Read all input files and construct continuous traces, using the
   * libmseed MSTraceList.  For each trace maintain a list of each
   * data record that contributed to the trace, implemented as a
   * RecordMap struct (MSTraceSeg->prvtptr) where a linked list of
   * Record structs is maintained.  The records are always listed in
   * time order.
   */
  
  flp = filelist;
  
  while ( flp != 0 )
    {
      /* Add '.orig' suffix to input file if it will be replaced */
      if ( replaceinput )
	{
	  /* The output file name is the original input file name */
	  flp->outfilename = flp->infilename;
	  
	  infilenamelen = strlen(flp->outfilename) + 6;
	  
	  if ( ! (flp->infilename = (char *) malloc (infilenamelen)) )
	    {
	      ms_log (2, "Cannot allocate memory for input file name\n");
	      return -1;
	    }
	  
	  snprintf (flp->infilename, infilenamelen, "%s.orig", flp->outfilename);
	  
	  if ( rename (flp->outfilename, flp->infilename) )
	    {
	      ms_log (2, "Cannot rename %s -> %s : '%s'\n",
		      flp->outfilename, flp->infilename, strerror(errno));
	      return -1;
	    }
	}
      
      if ( verbose )
	{
	  if ( replaceinput )
            {
              if ( flp->startoffset || flp->endoffset )
                ms_log (1, "Reading: %s (was %s) [range %"PRIu64":%"PRIu64"]\n",
                        flp->infilename, flp->outfilename, flp->startoffset, flp->endoffset);
              else
                ms_log (1, "Reading: %s (was %s)\n",
                        flp->infilename, flp->outfilename);
            }
	  else
            {
              if ( flp->startoffset || flp->endoffset )
                ms_log (1, "Reading: %s [range %"PRIu64":%"PRIu64"]\n",
                        flp->infilename, flp->startoffset, flp->endoffset);
              else
                ms_log (1, "Reading: %s\n", flp->infilename);
            }
	}
      
      /* Instruct libmseed to start at specified offset by setting a negative file position */
      fpos = - flp->startoffset;  /* Unset value is a 0, making this a non-operation */
      
      /* Loop over the input file */
      while ( (retcode = ms_readmsr_main (&msfp, &msr, flp->infilename, reclen, &fpos, NULL, 1, 0, selections, verbose-2))
	      == MS_NOERROR )
	{
          /* Break out as EOF if we have read past end offset */
          if ( flp->endoffset > 0 && fpos >= flp->endoffset )
            {
              retcode = MS_ENDOFFILE; 
              break;
            }
          
          recstarttime = msr->starttime;
	  recendtime = msr_endtime (msr);
	  
	  /* Generate the srcname with the quality code */
	  msr_srcname (msr, srcname, 1);
	  
	  /* Check if record should be skipped due to zero samples */
	  if ( skipzerosamps && msr->samplecnt == 0 )
	    {
	      if ( verbose >= 3 )
		{
		  ms_hptime2seedtimestr (recstarttime, stime, 1);
		  ms_log (1, "Skipping (zero samples) %s, %s\n", srcname, stime);
		}
	      continue;
	    }
	  
	  /* Check if record matches start time criteria: starts after or contains starttime */
	  if ( (starttime != HPTERROR) && (recstarttime < starttime && ! (recstarttime <= starttime && recendtime >= starttime)) )
	    {
	      if ( verbose >= 3 )
		{
		  ms_hptime2seedtimestr (recstarttime, stime, 1);
		  ms_log (1, "Skipping (starttime) %s, %s\n", srcname, stime);
		}
	      continue;
	    }
	  
	  /* Check if record matches end time criteria: ends after or contains endtime */
	  if ( (endtime != HPTERROR) && (recendtime > endtime && ! (recstarttime <= endtime && recendtime >= endtime)) )
	    {
	      if ( verbose >= 3 )
		{
		  ms_hptime2seedtimestr (recstarttime, stime, 1);
		  ms_log (1, "Skipping (endtime) %s, %s\n", srcname, stime);
		}
	      continue;
	    }
	  
	  /* Check if record is matched by the match regex */
	  if ( match )
	    {
	      if ( regexec ( match, srcname, 0, 0, 0) != 0 )
		{
		  if ( verbose >= 3 )
		    {
		      ms_hptime2seedtimestr (recstarttime, stime, 1);
		      ms_log (1, "Skipping (match) %s, %s\n", srcname, stime);
		    }
		  continue;
		}
	    }
	  
	  /* Check if record is rejected by the reject regex */
	  if ( reject )
	    {
	      if ( regexec ( reject, srcname, 0, 0, 0) == 0 )
		{
		  if ( verbose >= 3 )
		    {
		      ms_hptime2seedtimestr (recstarttime, stime, 1);
		      ms_log (1, "Skipping (reject) %s, %s\n", srcname, stime);
		    }
		  continue;
		}
	    }
	  
	  /* Check if record is matched by selection */
	  if ( selections )
	    {
	      if ( ! (matchsp = ms_matchselect (selections, srcname, recstarttime, recendtime, &matchstp)) )
		{
		  if ( verbose >= 3 )
		    {
		      ms_hptime2seedtimestr (recstarttime, stime, 1);
		      ms_log (1, "Skipping (selection) %s, %s\n", srcname, stime);
		    }
		  continue;
		}
	    }
	  
	  if ( verbose > 2 )
	    msr_print (msr, verbose - 3);
	  
	  /* Add record to the MSTraceList */
	  if ( ! (seg = mstl_addmsr (*ppmstl, msr, bestquality, 0, timetol, sampratetol)) )
	    {
	      ms_hptime2seedtimestr (recstarttime, stime, 1);
	      ms_log (2, "Cannot add record to trace list, %s, %s\n", srcname, stime);
	      continue;
	    }
	  
	  /* Determine where the record fit with this MSTraceSeg
	   * whence:
	   *   0 = New MSTrace
	   *   1 = End of MSTrace
	   *   2 = Beginning of MSTrace
	   */
	  whence = 0;
	  if ( seg->prvtptr )
	    {
	      if ( seg->endtime == recendtime )
		whence = 1;
	      else if ( seg->starttime == recstarttime )
		whence = 2;
	      else if ( recendtime == recstarttime )
		{
		  /* Determine best fit for records with no span */
		  if ( llabs (recstarttime - seg->endtime) < llabs (recstarttime - seg->starttime) )
		    whence = 1;
		  else
		    whence = 2;
		}
	      else
		{
		  ms_log (2, "Cannot determine where record fit relative to trace segment\n");
		  msr_print (msr, 1);
		  continue;
		}
	    }
	  
	  /* Create and populate new Record structure */
	  if ( ! (rec = (Record *) malloc (sizeof(Record))) )
	    {
	      ms_log (2, "Cannot allocate memory for Record entry\n");
	      return -1;
	    }

	  rec->flp = flp;
	  rec->offset = fpos;
	  rec->stageoffset = -1;
	  rec->reclen = msr->reclen;
	  rec->starttime = recstarttime;
	  rec->endtime = recendtime;
	  rec->quality = msr->dataquality;
	  rec->selectstart = HPTERROR;
	  rec->selectend = HPTERROR;
	  rec->newstart = HPTERROR;
	  rec->newend = HPTERROR;
	  rec->prev = 0;
	  rec->next = 0;
	  
	  /* Populate a new record map */
	  newrecmap.recordcnt = 1;
	  newrecmap.first = rec;
	  newrecmap.last = rec;
	  
	  /* If record is not completely selected search for joint selection limits */
	  if ( matchstp && ! (matchstp->starttime <= recstarttime && matchstp->endtime >= recendtime) )
	    {
	      if ( findselectlimits (matchsp, srcname, recstarttime, recendtime, rec) )
		{
		  ms_log (2, "Problem in findselectlimits(), please report\n");
		}
	    }
	  
	  /* If pruning at the sample level trim right at the start/end times */
	  if ( prunedata == 's' || prunedata == 'e' )
	    {
	      hptime_t seltime;
	      
	      /* Determine strictest start time (selection time or global start time) */
	      if ( starttime != HPTERROR && rec->selectstart != HPTERROR )
		seltime = ( starttime > rec->selectstart ) ? starttime : rec->selectstart;
	      else if ( rec->selectstart != HPTERROR )
		seltime = rec->selectstart;
	      else
		seltime = starttime;
	      
	      /* If the Record crosses the start time */
	      if ( seltime != HPTERROR && (seltime > recstarttime) && (seltime <= recendtime) )
		{
		  rec->newstart = seltime;
		}
	      
	      /* Determine strictest end time (selection time or global end time) */
	      if ( endtime != HPTERROR && rec->selectend != HPTERROR )
		seltime = ( endtime < rec->selectend ) ? endtime : rec->selectend;
	      else if ( rec->selectend != HPTERROR )
		seltime = rec->selectend;
	      else
		seltime = endtime;
	      
	      /* If the Record crosses the end time */
	      if ( seltime != HPTERROR && (seltime >= recstarttime) && (seltime < recendtime) )
		{
		  rec->newend = seltime;
		}
	    }
	  
	  /* Create extra Record structures if splitting on a time boundary */
	  if ( splitboundary )
	    {
	      BTime startbtime;
	      hptime_t boundary = HPTERROR;
	      hptime_t effstarttime;
	      
	      for (;;)
		{
		  effstarttime = (rec->newstart != HPTERROR) ? rec->newstart : rec->starttime;
		  ms_hptime2btime (effstarttime, &startbtime);
		  
		  /* Determine next split boundary */
		  if ( splitboundary == 'd' ) /* Days */
		    {
		      startbtime.day += 1;
		      startbtime.hour = startbtime.min = startbtime.sec = startbtime.fract = 0;
		      boundary = ms_btime2hptime (&startbtime);
		    }
		  else if ( splitboundary == 'h' ) /* Hours */
		    {
		      startbtime.hour += 1;
		      startbtime.min = startbtime.sec = startbtime.fract = 0;
		      boundary = ms_btime2hptime (&startbtime);
		    }
		  else if ( splitboundary == 'm' ) /* Minutes */
		    {
		      startbtime.min += 1;
		      startbtime.sec = startbtime.fract = 0;
		      boundary = ms_btime2hptime (&startbtime);
		    }
		  else
		    {
		      ms_log (2, "Split boundary code unrecognized: '%c'\n", splitboundary);
		      break;
		    }
		  
		  /* If end time is beyond the boundary create a new Record */
		  if ( rec->endtime >= boundary )
		    {
		      if ( ! (newrec = (Record *) malloc (sizeof(Record))) )
			{
			  ms_log (2, "Cannot allocate memory for Record entry");
			  return -1;
			}
		      
		      memcpy (newrec, rec, sizeof(Record));
		      
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
	  
	  /* Add the new Record(s) to the RecordMap associated with the MSTraceSeg */
	  
	  /* Add new Record(s) to end of the RecordMap */
	  if ( whence == 1 )
	    {
	      recmap = (RecordMap *) seg->prvtptr;
	      
	      recmap->last->next = newrecmap.first;
	      newrecmap.first->prev = recmap->last;
	      
	      recmap->last = newrecmap.last;
	      
	      recmap->recordcnt += newrecmap.recordcnt;
	    }
	  /* Add new Record(s) to beginning of the RecordMap */
	  else if ( whence == 2 )
	    {
	      recmap = (RecordMap *) seg->prvtptr;
	      
	      recmap->first->prev = newrecmap.last;
	      newrecmap.last->next = recmap->first;
	      
	      recmap->first = newrecmap.first;
	      
	      recmap->recordcnt += newrecmap.recordcnt;
	      
	      /* Increment reordered count */
	      flp->reordercount++;
	    }
	  /* First Record(s) for this MSTraceSeg, allocate RecordMap */
	  else
	    {
	      if ( seg->prvtptr )
		ms_log (2, "Supposedly first record, but RecordMap not empty, report this\n");
	      
	      if ( ! (recmap = (RecordMap *) malloc (sizeof(RecordMap))) )
		{
		  ms_log (2, "Cannot allocate memory for new RecordMap\n");
		  return -1;
		}
	      recmap->first = newrecmap.first;
	      recmap->last = newrecmap.last;
	      recmap->recordcnt = newrecmap.recordcnt;

	      seg->prvtptr = recmap;
	    }
	  
	  /* Copy record to staging buffer if space is available */
	  if ( stagebuffer && ! stagefull )
	    {
	      if ( (stagelength - stageoffset) >= msr->reclen )
		{
		  /* Copy record to end of staging buffer contents */
		  memcpy ((stagebuffer + stageoffset), msr->record, msr->reclen);
		  
		  /* Set offset to record in staging buffer */
		  rec->stageoffset = stageoffset;
		  
		  /* Update offset to end of staging buffer */
		  stageoffset += msr->reclen;
		}
	      else if ( ! stagefull )
		{
		  ms_log (1, "Warning: record staging buffer full, using input files for remaining records\n");
		  
		  stagefull = 1;
		}
	    }
	  
	  totalrecs++;
	  totalsamps += msr->samplecnt;
          
          /* Break out as EOF if record is at or beyond end offset */
          if ( flp->endoffset > 0 && (fpos + msr->reclen) >= flp->endoffset )
            {
              retcode = MS_ENDOFFILE; 
              break;
            }
        } /* End of looping through records in file */
      
      /* Critical error if file was not read properly */
      if ( retcode != MS_ENDOFFILE )
	{
	  ms_log (2, "Cannot read %s: %s\n", flp->infilename, ms_errorstr(retcode));
	  ms_readmsr_main (&msfp, &msr, NULL, 0, NULL, NULL, 0, 0, NULL, 0);
	  return -1;
	}
      
      /* Make sure everything is cleaned up */
      ms_readmsr_main (&msfp, &msr, NULL, 0, NULL, NULL, 0, 0, NULL, 0);
      
      totalfiles++;
      flp = flp->next;
    } /* End of looping over file list */
  
  /* Increase open file limit if necessary, in general we need the
   * filecount + ds_maxopenfiles and some wiggle room. */
  setofilelimit (totalfiles + ds_maxopenfiles + 20);
  
  if ( basicsum )
    ms_log (0, "Files: %d, Records: %d, Samples: %d\n", totalfiles, totalrecs, totalsamps);
  
  return 0;
}  /* End of readfiles() */


/***************************************************************************
 * processtraces:
 *
 * Process all data represented by the MSTraceList by first pruning
 * them and then writing out the remaining data.
 *
 * Returns 0 on success and -1 on error.
 ***************************************************************************/
static int
processtraces (MSTraceList *mstl)
{
  if ( verbose > 2 )
    printtracemap (mstl);
  
  /* Prune data */
  if ( prunedata )
    {
      /* Prune overlaps */
      if ( prunedata == 'r' || prunedata == 's' )
	if ( prunetraces (mstl) )
	  return -1;
      
      /* Reconcile MSTraceID times with associated record maps */
      if ( reconcile_tracetimes (mstl) )
	return -1;
    }
  
  /* Enforce minimum segment length by removing shorter */
  if ( minseglength )
    {
      if ( minsegmentlength (mstl, minseglength) )
	return -1;
    }
  
  /* Trim segment list to longest segment */
  if ( longestonly )
    {
      if ( longestsegmentonly (mstl) )
	return -1;
    }
  
  /* Write all MSTraceSeg associated records to output file(s) */
  if ( writetraces (mstl) )
    return -1;
  
  return 0;
}  /* End of processtraces() */


/***************************************************************************
 * writetraces():
 *
 * Write all MSTraceSeg associated records to output file(s).  If an
 * output file is specified all records will be written to it,
 * otherwise records will be written to the original files and
 * (optionally) backups of the original files will remain.
 * 
 * This routine will also call trimrecord() to trim a record when data
 * suturing is requested.  Record trimming is triggered when
 * Record.newstart or Record.newend are set for any output records.
 *
 * The quality flag is optionally set for all output records.
 *
 * Returns 0 on success and 1 on error.
 ***************************************************************************/
static int
writetraces (MSTraceList *mstl)
{
  static uint64_t totalrecsout = 0;
  static uint64_t totalbytesout = 0;
  char *recordptr = NULL;
  char *wb = "wb";
  char *ab = "ab";
  char *mode;
  char errflag = 0;
  long suffix = 0;
  int rv;
  
  MSTraceID *id;
  MSTraceID *groupid;
  MSTraceSeg *seg;
  
  RecordMap *groupmap = 0;
  RecordMap *recmap = 0;
  Record *rec;
  Record *recnext;
  Filelink *flp;
  
  FILE *ofp = 0;
  
  if ( ! mstl )
    return 1;
  
  if ( ! mstl->traces )
    return 1;
  
  if ( verbose )
    ms_log (1, "Writing output data\n");
  
  /* Open the output file if specified */
  if ( outputfile )
    {
      /* Decide if we are appending or overwriting */
      mode = ( totalbytesout || outputmode ) ? ab : wb;
      
      if ( strcmp (outputfile, "-") == 0 )
        {
          ofp = stdout;
        }
      else if ( (ofp = fopen (outputfile, mode)) == NULL )
        {
          ms_log (2, "Cannot open output file: %s (%s)\n",
		  outputfile, strerror(errno));
          return 1;
        }
    }
  
  /* Re-link records into write lists, from per-segment lists to per ID lists.
   * This allows (later) sorting of data records as logical groups regardless
   * from which segment the record was originally associated. */
  id = mstl->traces;
  groupid = id;
  while ( id )
    {
      /* Determine if this is a new group */
      if ( groupid != id )
	{
	  /* If data was pruned, group by network, station, location and channel.
	   * When pruning, data is grouped by NSLC and quality.  This groups the write 
	   * lists by NSLC without quality considered. */
	  if ( prunedata )
	    {
	      if ( strcmp (id->network, groupid->network) || strcmp (id->station, groupid->station) ||
		   strcmp (id->location, groupid->location) || strcmp (id->channel, groupid->channel) )
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
      if ( seg )
        {
          /* Allocate group ID RecordMap if needed */
          if ( ! groupid->prvtptr )
	    {
	      if ( ! (id->prvtptr = calloc (1, sizeof(RecordMap))) )
	        {
	          ms_log (2, "writetraces(): Cannot allocate memory\n");
	          return 1;
	        }
	      groupmap = (RecordMap *) id->prvtptr;
	    }
	  
          while ( seg )
	    {
	      recmap = (RecordMap *) seg->prvtptr;
	      rec = recmap->first;
	      while ( rec && ! errflag )
	        {
	          recnext = rec->next;
		  
	          /* Free and skip marked (pre-identified as non-contributing) records */
	          if ( rec->reclen == 0 )
		    {
		      free (rec);
		      rec = recnext;
		      continue;
		    }
		  
	          /* Add record to group ID write list */
	          if ( ! groupmap->first )
		    {
		      groupmap->first = rec;
		      rec->prev = 0;
		      rec->next = 0;
		      groupmap->last = rec;
		      groupmap->recordcnt = 1;
		    }
	          else
		    {
		      groupmap->last->next = rec;
		      rec->prev = groupmap->last;
		      rec->next = 0;
		      groupmap->last = rec;
		      groupmap->recordcnt++;
		    }
		  
	          rec = recnext;
	        } /* Done looping through Records in the RecordMap */
	      
	      /* Free segment RecordMap and cauterize */
	      if ( seg->prvtptr )
	        free (seg->prvtptr);
	      seg->prvtptr = 0;
	      
	      seg = seg->next;
	    } /* Done looping through MSTraceSegs in the MSTraceID */
        }
      
      id = id->next;
    } /* Done looping through MSTraceIDs in the MSTraceList */  
  
  /* Loop through MSTraceList and write records in write lists */
  id = mstl->traces;
  while ( id && errflag != 1 )
    {
      /* Skip when no write list is present */
      if ( id->prvtptr == 0 )
	{
	  id = id->next;
	  continue;
	}
      
      /* Reset error flag for continuation errors */
      if ( errflag == 2 )
	errflag = 0;
      
      /* Set suffix if splitting on record length changes */
      if ( splitreclen )
	suffix = 1;
      
      recmap = (RecordMap *) id->prvtptr;
      
      /* Sort record list if overlaps have been pruned, if the data has not been
       * pruned it is already in time order. */
      if ( prunedata == 'r' || prunedata == 's' )
	{
          sortrecmap (recmap);
	}
      
      /* Loop through each Record in the write list RecordMap.
       * After records are read from the input files, perform any
       * pre-identified pruning before writing data back out */
      rec = recmap->first;
      while ( rec && ! errflag )
	{
	  /* Set record pointer to location in the staging buffer */
	  if ( rec->stageoffset >= 0 )
	    {
	      recordptr = stagebuffer + rec->stageoffset;
	    }
	  /* Otherwise read the record from the input file */
	  else
	    {
	      /* Make sure the record buffer is large enough */
	      if ( rec->reclen > sizeof(recordbuf) )
		{
		  ms_log (2, "Record length (%d bytes) larger than buffer (%llu bytes)\n",
			  rec->reclen, (long long unsigned int) sizeof(recordbuf));
		  errflag = 1;
		  break;
		}
	      
	      /* Open file for reading if not already done */
	      if ( ! rec->flp->infp )
		if ( ! (rec->flp->infp = fopen (rec->flp->infilename, "rb")) )
		  {
		    ms_log (2, "Cannot open '%s' for reading: %s\n",
			    rec->flp->infilename, strerror(errno));
		    errflag = 1;
		    break;
		  }
	      
	      /* Seek to record offset */
	      if ( lmp_fseeko (rec->flp->infp, rec->offset, SEEK_SET) == -1 )
		{
		  ms_log (2, "Cannot seek in '%s': %s\n",
			  rec->flp->infilename, strerror(errno));
		  errflag = 1;
		  break;
		}
	      
	      /* Read record into buffer */
	      if ( fread (recordbuf, rec->reclen, 1, rec->flp->infp) != 1 )
		{
		  ms_log (2, "Cannot read %d bytes at offset %llu from '%s'\n",
			  rec->reclen, (long long unsigned)rec->offset,
			  rec->flp->infilename);
		  errflag = 1;
		  break;
		}
	      
	      recordptr = recordbuf;
	    }
	  
	  /* Trim data from the record if new start or end times are specifed */
	  if ( rec->newstart != HPTERROR || rec->newend != HPTERROR )
	    {
	      rv = trimrecord (rec, recordptr);
	      
	      if ( rv == -1 )
		{
		  rec = rec->next;
		  continue;
		}
	      if ( rv == -2 )
		{
		  ms_log (2, "Cannot unpack Mini-SEED from byte offset %lld in %s\n",
			  rec->offset, rec->flp->infilename);
		  ms_log (2, "  Expecting %s, skipping the rest of this channel\n", id->srcname);
		  errflag = 2;
		  break;
		}
	      
	      /* Trimmed records are re-packed in the global record buffer */
	      recordptr = recordbuf;
	    }
	  
	  /* Re-stamp quality indicator if specified */
	  if ( restampqind )
	    {
	      if ( verbose > 3 )
		ms_log (1, "Re-stamping data quality indicator to '%c'\n", restampqind);
	      
	      *(recordptr + 6) = restampqind;
	    }
	  
	  /* Write to a single output file if specified */
	  if ( ofp )
	    {
	      if ( fwrite (recordptr, rec->reclen, 1, ofp) != 1 )
		{
		  ms_log (2, "Cannot write to '%s'\n", outputfile);
		  errflag = 1;
		  break;
		}
	    }
	  
	  /* Write to Archive(s) if specified and/or add to written list */
	  if ( archiveroot || writtenfile )
	    {
	      MSRecord *msr = 0;
	      Archive *arch;
	      
	      if ( msr_unpack (recordptr, rec->reclen, &msr, 0, verbose) != MS_NOERROR )
		{
		  ms_log (2, "Cannot unpack Mini-SEED from byte offset %lld in %s, file changed?\n",
			  rec->offset, rec->flp->infilename);
		  ms_log (2, "  Expecting %s, skipping the rest of this channel\n", id->srcname);
		  errflag = 2;
		  break;
		}
	      
	      if ( archiveroot )
		{
		  arch = archiveroot;
		  while ( arch )
		    {
		      ds_streamproc (&arch->datastream, msr, suffix, verbose-1);
		      arch = arch->next;
		    }
		}
	      
	      if ( writtenfile )
		{
		  MSTraceSeg *seg;
		  
		  if ( (seg = mstl_addmsr (writtentl, msr, 1, 1, -1.0, -1.0)) == NULL )
		    {
		      ms_log (2, "Error adding MSRecord to MSTraceList, bah humbug.\n");
		    }
		  else
		    {
		      if ( ! seg->prvtptr )
			{
			  if ( (seg->prvtptr = malloc (sizeof(int64_t))) == NULL )
			    {
			      ms_log (2, "Error allocating memory for written count, bah humbug.\n");
			      errflag = 1;
			      break;
			    }
			  else
			    {
			      *((int64_t *)seg->prvtptr) = 0;
			    }
			}
		      
		      *((int64_t *)seg->prvtptr) += msr->reclen;
		    }
		}
	      
	      msr_free (&msr);
	    }
	  
	  /* Open original file for output if replacing input and write */
	  if ( replaceinput )
	    {
	      if ( ! rec->flp->outfp )
		if ( ! (rec->flp->outfp = fopen (rec->flp->outfilename, "wb")) )
		  {
		    ms_log (2, "Cannot open '%s' for writing: %s\n",
			    rec->flp->outfilename, strerror(errno));
		    errflag = 1;
		    break;
		  }
	      
	      if ( fwrite (recordptr, rec->reclen, 1, rec->flp->outfp) != 1 )
		{
		  ms_log (2, "Cannot write to '%s'\n", rec->flp->outfilename);
		  errflag = 1;
		  break;
		}
	    }
	  
	  /* Update file entry time stamps and counts */
	  if ( ! rec->flp->earliest || (rec->flp->earliest > rec->starttime) )
	    {
	      rec->flp->earliest = rec->starttime;
	    }
	  if ( ! rec->flp->latest || (rec->flp->latest < rec->endtime) )
	    {
	      rec->flp->latest = rec->endtime;
	    }
	  
	  rec->flp->byteswritten += rec->reclen;
	  
	  totalrecsout++;
	  totalbytesout += rec->reclen;
	  
	  /* Increment suffix if splitting and record length changes */
	  if ( splitreclen && rec->next && rec->next->reclen != rec->reclen )
	    {
	      suffix++;
	    }
	  
	  rec = rec->next;
	} /* Done looping through Records in the RecordMap */
      
      id = id->next;
    } /* Done looping through MSTraceIDs the MSTraceList */
  
  /* Close all open input & output files and remove backups if requested */
  flp = filelist;
  while ( flp )
    {
      if ( ! ofp && verbose )
	{
	  if ( replaceinput )
	    ms_log (1, "Wrote %d bytes from file %s (was %s)\n",
		    flp->byteswritten, flp->infilename, flp->outfilename);
	  else
	    ms_log (1, "Wrote %d bytes from file %s\n",
		    flp->byteswritten, flp->infilename);
	}
      
      if ( flp->infp )
	{
	  fclose (flp->infp);
	  flp->infp = 0;
	}
      
      if ( flp->outfp )
	{
	  fclose (flp->outfp);
	  flp->outfp = 0;
	}
      
      if ( nobackups && ! ofp )
	if ( unlink (flp->infilename) )
	  {
	    ms_log (2, "Cannot remove '%s': %s\n",
		    flp->infilename, strerror(errno));
	    errflag = 1;
	  }
      
      flp = flp->next;
    }
  
  /* Close output file if used */
  if ( ofp )
    {
      fclose (ofp);
      ofp = 0;
    }
  
  if ( verbose )
    {
      ms_log (1, "Wrote %llu bytes of %llu records to output file(s)\n",
	      totalbytesout, totalrecsout);
    }
  
  return (errflag) ? 1 : 0;
}  /* End of writetraces() */


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
trimrecord (Record *rec, char *recordbuf)
{
  MSRecord *msr = 0;
  hptime_t hpdelta;

  char srcname[50];
  char stime[30];
  char etime[30];
  
  int trimsamples;
  int samplesize;
  int64_t packedsamples;
  int packedrecords;
  int retcode;
  
  if ( ! rec || ! recordbuf )
    return -1;

  srcname[0] = '\0';
  stime[0] = '\0';
  etime[0] = '\0';
  
  /* Sanity check for new start/end times */
  if ( (rec->newstart != HPTERROR && rec->newend != HPTERROR && rec->newstart > rec->newend) ||
       (rec->newstart != HPTERROR && (rec->newstart < rec->starttime || rec->newstart > rec->endtime)) ||
       (rec->newend != HPTERROR && (rec->newend > rec->endtime || rec->newend < rec->starttime)) )
    {
      ms_log (2, "Problem with new start/end record bound times.\n");
      ms_recsrcname (recordbuf, srcname, 1);
      ms_log (2, "  Original record %s from %s (byte offset: %llu)\n",
	      srcname, rec->flp->infilename, (unsigned long long)rec->offset);
      ms_hptime2seedtimestr (rec->starttime, stime, 1);
      ms_hptime2seedtimestr (rec->endtime, etime, 1);
      ms_log (2, "       Start: %s       End: %s\n", stime, etime);
      if ( rec->newstart == HPTERROR ) strcpy (stime, "NONE");
      else ms_hptime2seedtimestr (rec->newstart, stime, 1);
      if ( rec->newend == HPTERROR ) strcpy (etime, "NONE");
      else ms_hptime2seedtimestr (rec->newend, etime, 1);
      ms_log (2, " Start bound: %-24s End bound: %-24s\n", stime, etime);
      
      return -1;
    }
  
  /* Unpack data record header without data samples */
  if ( (retcode = msr_unpack(recordbuf, rec->reclen, &msr, 0, 0)) != MS_NOERROR )
    {
      ms_log (2, "Cannot unpack Mini-SEED record: %s\n", ms_errorstr(retcode));
      return -2;
    }
  
  /* Check for missing B1000 or unsupported data encoding, can only trim what can be packed */
  if ( ! msr->Blkt1000 || (msr->encoding != DE_INT16 && msr->encoding != DE_INT32 &&
			   msr->encoding != DE_FLOAT32 && msr->encoding != DE_FLOAT64 &&
			   msr->encoding != DE_STEIM1 && msr->encoding != DE_STEIM2) )
    {
      if ( verbose )
	{
	  msr_srcname (msr, srcname, 0);
	  ms_hptime2seedtimestr (rec->starttime, stime, 1);
	  if ( ! msr->Blkt1000 )
	    ms_log (1, "Skipping trim of %s (%s), missing blockette 1000\n",
		    srcname, stime, msr->encoding);
	  else if ( msr->encoding == DE_ASCII )
	    ms_log (1, "Skipping trim of %s (%s), ASCII encoded data\n",
		    srcname, stime, msr->encoding);
	  else
	    ms_log (1, "Skipping trim of %s (%s), unsupported encoding (%d: %s)\n",
		    srcname, stime, msr->encoding, ms_encodingstr(msr->encoding));
	}
      
      msr_free (&msr);
      return 0;
    }
  
  msr_free (&msr);
  
  /* Unpack data record header including data samples */
  if ( (retcode = msr_unpack(recordbuf, rec->reclen, &msr, 1, verbose-1)) != MS_NOERROR )
    {
      ms_log (2, "Cannot unpack Mini-SEED record: %s\n", ms_errorstr(retcode));
      return -2;
    }
  
  if ( verbose > 1 )
    {
      msr_srcname (msr, srcname, 0);
      ms_log (1, "Triming record: %s (%c)\n", srcname, msr->dataquality);
      ms_hptime2seedtimestr (rec->starttime, stime, 1);
      ms_hptime2seedtimestr (rec->endtime, etime, 1);
      ms_log (1, "       Start: %s        End: %s\n", stime, etime);
      if ( rec->newstart == HPTERROR ) strcpy (stime, "NONE");
      else ms_hptime2seedtimestr (rec->newstart, stime, 1);
      if ( rec->newend == HPTERROR ) strcpy (etime, "NONE");
      else ms_hptime2seedtimestr (rec->newend, etime, 1);
      ms_log (1, " Start bound: %-24s  End bound: %-24s\n", stime, etime);
    }
  
  /* Determine sample period in high precision time ticks */
  hpdelta = ( msr->samprate ) ? (hptime_t) (HPTMODULUS / msr->samprate) : 0;
  
  /* Remove samples from the beginning of the record */
  if ( rec->newstart != HPTERROR && hpdelta )
    {
      hptime_t newstarttime;
      
      /* Determine new start time and the number of samples to trim */
      trimsamples = 0;
      newstarttime = rec->starttime;
      
      while ( newstarttime < rec->newstart && trimsamples < msr->samplecnt )
	{
	  newstarttime += hpdelta;
	  trimsamples++;
	}
      
      if ( trimsamples >= msr->samplecnt )
	{
	  if ( verbose > 1 )
	    ms_log (1, "All samples would be trimmed from record, skipping\n");
	  
	  msr_free (&msr);
	  return -1;
	}
      
      if ( verbose > 2 )
	{
	  ms_hptime2seedtimestr (newstarttime, stime, 1);
	  ms_log (1, "Removing %d samples from the start, new start time: %s\n", trimsamples, stime);
	}
      
      samplesize = ms_samplesize (msr->sampletype);
      
      memmove (msr->datasamples,
	       (char *)msr->datasamples + (samplesize * trimsamples),
	       samplesize * (msr->numsamples - trimsamples));
      
      msr->numsamples -= trimsamples;
      msr->samplecnt -= trimsamples;
      msr->starttime = newstarttime;
      rec->starttime = newstarttime;
    }
  
  /* Remove samples from the end of the record */
  if ( rec->newend != HPTERROR && hpdelta )
    {
      hptime_t newendtime;
      
      /* Determine new end time and the number of samples to trim */
      trimsamples = 0;
      newendtime = rec->endtime;
      
      while ( newendtime > rec->newend && trimsamples < msr->samplecnt )
	{
	  newendtime -= hpdelta;
	  trimsamples++;
	}
      
      if ( trimsamples >= msr->samplecnt )
	{
	  if ( verbose > 1 )
	    ms_log (1, "All samples would be trimmed from record, skipping\n");
	  
	  msr_free (&msr);
	  return -1;
	}
      
      if ( verbose > 2 )
	{
	  ms_hptime2seedtimestr (newendtime, etime, 1);
	  ms_log (1, "Removing %d samples from the end, new end time: %s\n", trimsamples, etime);
	}
      
      msr->numsamples -= trimsamples;
      msr->samplecnt -= trimsamples;
      rec->endtime = newendtime;
    }

  /* Repacking the record will apply any unapplied time corrections to the start time,
   * make sure the flag is set to indicate that the correction has been applied. */
  if ( msr->fsdh && msr->fsdh->time_correct != 0 && ! (msr->fsdh->act_flags & 0x02) )
    {
      msr->fsdh->act_flags |= (1<<1);
    }
  
  /* Pack the data record into the global record buffer used by writetraces() */
  packedrecords = msr_pack (msr, &record_handler, NULL, &packedsamples, 1, verbose-1);
  
  /* Clean up MSRecord */
  msr_free (&msr);
  
  return 0;
}  /* End of trimrecord() */


/***************************************************************************
 * record_handler():
 *
 * Used by trimrecord() to save repacked Mini-SEED to global record
 * buffer.
 ***************************************************************************/
static void
record_handler (char *record, int reclen, void *handlerdata)
{
  /* Copy record to global record buffer */
  memcpy (recordbuf, record, reclen);
  
}  /* End of record_handler() */


/***************************************************************************
 * prunetraces():
 *
 * Prune all redundant data from the RecordMap entries associated with
 * the specified MSTraceSegs.
 *
 * For each MSTraceSeg determine the coverage of the RecordMap associated
 * with each overlapping, higher-priority MSTraceSeg using findcoverage().
 * If some higher-priority overlap was determined to exist modify the
 * RecordMap of the MSTraceSeg in question to mark the overlapping data
 * using trimtrace().
 *
 * Return 0 on success and -1 on failure.
 ***************************************************************************/
static int
prunetraces (MSTraceList *mstl)
{
  MSTraceID *id = 0;
  MSTraceSeg *seg = 0;
  MSTraceGroup *coverage = 0;
  int retval;
  
  if ( ! mstl )
    return -1;
  
  if ( ! mstl->traces )
    return -1;
  
  if ( verbose )
    ms_log (1, "Pruning trace data\n");
  
  /* For each MSTraceSeg determine the coverage of the overlapping
     Records from the other traces with a higher priority and prune
     the overlap. */
  id = mstl->traces;
  while ( id )
    {
      seg = id->first;
      while ( seg )
	{
	  /* Determine overlapping trace coverage */
	  retval = findcoverage (mstl, id, seg, &coverage);
	  
	  if ( retval )
	    {
	      ms_log (2, "cannot findcoverage()\n");
	      return -1;
	    }
	  else if ( coverage )
	    {
	      if ( trimtrace (seg, id->srcname, coverage) < 0 )
		{
		  ms_log (2, "cannot trimtraces()\n");
		  return -1;
		}
	    }
	  
	  /* Free the coverage MSTraceGroup for reuse */
	  mst_freegroup (&coverage);
	  
	  seg = seg->next;
	}
      
      id = id->next;
    }
  
  return 0;
}  /* End of prunetraces() */


/***************************************************************************
 * findcoverage():
 *
 * Search an MSTraceList for entries that overlap the target MSTraceSeg
 * and, from the Record entries of the overlapping MSTraceSegs, build a
 * coverage map in the form of a sparsely-populated MSTraceGroup.
 *
 * Only data with a higher priority than the target MSTraceSeg will be
 * added to the overlap coverage.  Priority is determined using the
 * quality codes via qcompare() and if the qualities are equal the
 * longest time-series will be given priority.
 *
 * On success a new MSTraceGroup will be allocated and returned, it is
 * up to the caller to properly free this memory.
 *
 * When no overlap coverage is found *ppcoverage will be 0, otherwise
 * it will contain a list of MSTraces representing the overlap
 * coverage.
 *
 * Returns 0 on success and -1 on error.
 ***************************************************************************/
static int
findcoverage (MSTraceList *mstl, MSTraceID *targetid, MSTraceSeg *targetseg,
	      MSTraceGroup **ppcoverage)
{
  MSTraceID *id = 0;
  MSTraceSeg *seg = 0;
  MSTrace *cmst = 0;
  MSTrace *prevcmst = 0;
  RecordMap *recmap;
  Record *rec;
  hptime_t hpdelta, hptimetol;
  hptime_t effstarttime, effendtime;
  int priority;
  int newsegment;
  
  if ( ! mstl || ! targetseg || ! ppcoverage )
    return -1;
  
  /* Reset coverage MSTraceGroup if needed */
  if ( *ppcoverage )
    mst_freegroup (ppcoverage);
  
  /* Determine sample period in high precision time ticks */
  hpdelta = ( targetseg->samprate ) ? (hptime_t) (HPTMODULUS / targetseg->samprate) : 0;
  
  /* Determine time tolerance in high precision time ticks */
  hptimetol = ( timetol == -1 ) ? (hpdelta / 2) : (hptime_t) (HPTMODULUS * timetol);
  
  /* Loop through each MSTraceID in the list */
  id = mstl->traces;
  while ( id )
    {
      /* Continue with next if network, station, location or channel are different */
      if ( targetid != id )
	{
	  if ( strcmp (id->network, targetid->network) || strcmp (id->station, targetid->station) ||
	       strcmp (id->location, targetid->location) || strcmp (id->channel, targetid->channel) )
	    {
	      id = id->next;
	      continue;
	    }
	}
      
      prevcmst = 0;
      seg = id->first;
      while ( seg )
	{
	  /* Skip target segment */
	  if ( seg == targetseg )
	    {
	      seg = seg->next;
	      continue;
	    }
	  
	  /* Stop searching if target segment is before segment start time,
	   * assuming the segments are in time order nothing later will overlap. */
	  if ( (targetseg->endtime + hptimetol) < seg->starttime )
	    {
	      break;
	    }
	  
	  /* Skip out-of-band (0 samprate) trace */
	  if ( seg->samprate == 0.0 )
	    {
	      seg = seg->next;
	      continue;
	    }
	  
	  /* Continue with next if sample rate are different */
	  if ( ! MS_ISRATETOLERABLE (seg->samprate, targetseg->samprate) )
	    {
	      seg = seg->next;
	      continue;
	    }

          /* Check for duplicate or overlap NSLCQ last coverage entry */
          if ( prevcmst )
            {
              /* At this point the NSLCQ and rate are the same, check if the
               * segment is completly contained by the previous coverage entry. */
              if ( seg->starttime >= prevcmst->starttime &&
                   seg->endtime <= prevcmst->endtime )
                {
                  seg = seg->next;
                  continue;
                }
            }
          
	  /* Test for overlap with targetseg */
	  if ( (targetseg->endtime + hptimetol) >= seg->starttime &&
	       (targetseg->starttime - hptimetol) <= seg->endtime )
	    {
	      /* Determine priority:
	       *  -1 : seg > targetseg
	       *   0 : seg == targetseg
	       *   1 : seg < targetseg */
	      priority = 0;
	      
	      /* If best quality is requested compare the qualities to determine priority */
	      if ( bestquality )
		priority = qcompare (id->dataquality, targetid->dataquality);
	      
	      /* If priorities are equal (qualities are equal or no checking)
	       * give priority to the longest segment */
	      if ( priority == 0 )
		{
		  if ( (seg->endtime - seg->starttime) >= (targetseg->endtime - targetseg->starttime) )
		    priority = -1;
		  else
		    priority = 1;
		}

	      /* If overlapping trace is a higher priority than targetseg add to coverage */
	      if ( priority == -1 )
		{
		  /* Loop through list of associated Records, determine
		     contiguous coverage and store in an MSTraceGroup */
		  recmap = (RecordMap *) seg->prvtptr;
		  rec = recmap->first;
		  newsegment = 1;
		  while ( rec )
		    {
		      /* Check if record has been marked as non-contributing */
		      if ( rec->reclen == 0 )
			{
			  rec = rec->next;
			  continue;
			}
		      
		      /* Determine effective record start and end times */
		      effstarttime = ( rec->newstart != HPTERROR ) ? rec->newstart : rec->starttime;
		      effendtime = ( rec->newend != HPTERROR ) ? rec->newend : rec->endtime;
		      
		      /* Create a new segment if a break in the time-series is detected */
		      if ( cmst )
			if ( llabs((cmst->endtime + hpdelta) - effstarttime) > hptimetol )
			  newsegment = 1;
		      
		      if ( newsegment )
			{
			  newsegment = 0;
			  
			  cmst = mst_init (NULL);
			  
			  if ( cmst )
			    {
			      if ( ! *ppcoverage )
				*ppcoverage = mst_initgroup (NULL);
			      
			      mst_addtracetogroup (*ppcoverage, cmst);

			      prevcmst = cmst;
			      cmst->dataquality = id->dataquality;
			      cmst->samprate = seg->samprate;
			      cmst->starttime = effstarttime;
			    }
			}
		      
		      if ( cmst )
			cmst->endtime = effendtime;
		      else
			ms_log (2, "ACK! covergage MSTrace is not allocated!?  PLEASE REPORT\n");
		      
		      rec = rec->next;
		    }
		}
	    }
	  
	  seg = seg->next;
	}
      
      id = id->next;
    }
  
  return 0;
}  /* End of findcoverage() */


/***************************************************************************
 * trimtrace():
 *
 * Adjust Record entries associated with the target MSTraceSeg that
 * are overlapping the time represented by the coverage MSTraceGroup
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
trimtrace (MSTraceSeg *targetseg, char *targetsrcname, MSTraceGroup *coverage)
{
  RecordMap *recmap;
  Record *rec;
  MSTrace *cmst;
  hptime_t effstarttime, effendtime;
  hptime_t hpdelta, hptimetol;
  char stime[30];
  char etime[30];
  int modcount = 0;
  
  if ( ! targetseg || ! coverage )
    return -1;
  
  /* Determine sample period in high precision time ticks */
  hpdelta = ( targetseg->samprate ) ? (hptime_t) (HPTMODULUS / targetseg->samprate) : 0;
  
  /* Determine time tolerance in high precision time ticks */
  hptimetol = ( timetol == -1 ) ? (hpdelta / 2) : (hptime_t) (HPTMODULUS * timetol);
  
  /* Traverse the Record chain for the target MSTrace and mark Records
   * that are completely overlapped by the MSTraceGroup coverage */
  recmap = (RecordMap *) targetseg->prvtptr;
  rec = recmap->first;
  while ( rec )
    {
      cmst = coverage->traces;
      while ( cmst )
	{
          if ( ! rec->reclen ) /* Skip if marked non-contributing */
            break;
          
	  /* Determine effective record start and end times for comparison */
	  effstarttime = ( rec->newstart != HPTERROR ) ? rec->newstart : rec->starttime;
	  effendtime = ( rec->newend != HPTERROR ) ? rec->newend : rec->endtime;
	  
	  /* Use selection start and end for pruning if they are stricter */
	  if ( rec->selectstart != HPTERROR && rec->selectstart > effstarttime )
	    effstarttime = rec->selectstart;
	  if ( rec->selectend != HPTERROR && rec->selectend < effendtime )
	    effendtime = rec->selectend;
	  
	  /* Mark Record if it is completely overlapped by the coverage including tolerance */
	  if ( effstarttime >= (cmst->starttime - hptimetol) &&
	       effendtime <= (cmst->endtime + hptimetol) )
	    {
	      if ( verbose > 1 )
		{
		  ms_hptime2seedtimestr (rec->starttime, stime, 1);
		  ms_hptime2seedtimestr (rec->endtime, etime, 1);
		  ms_log (1, "Removing Record %s (%c) :: %s  %s  offset: %lld, reclen: %d\n",
			  targetsrcname, rec->quality, stime, etime,
                          (long long int)rec->offset, rec->reclen);
		}
	      
	      rec->flp->recrmcount++;
	      rec->reclen = 0;
	      modcount++;
	    }
	  
	  /* Determine the new start/end times if pruning at the sample level */
	  if ( prunedata == 's' && rec->reclen != 0 )
	    {
	      /* Record overlaps beginning of HP coverage */
	      if ( effstarttime < cmst->starttime &&
		   (effendtime + hptimetol) >= cmst->starttime )
		{
		  /* Set Record new end time boundary including specified time tolerance */
		  rec->newend = cmst->starttime - hpdelta + hptimetol;
		  
		  if ( (starttime != HPTERROR) && (rec->newend < starttime) )
		    {
		      if ( verbose > 1 )
			{
			  ms_hptime2seedtimestr (rec->starttime, stime, 1);
			  ms_hptime2seedtimestr (rec->endtime, etime, 1);
			  ms_log (1, "Removing Record %s (%c) :: %s  %s\n",
				  targetsrcname, rec->quality, stime, etime);
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
	      
	      /* Record overlaps end of HP coverage */
	      if ( (effstarttime - hptimetol) <= cmst->endtime &&
		   effendtime > cmst->endtime )
		{
		  /* Set Record new start time boundary including specified time tolerance */
		  rec->newstart = cmst->endtime + hpdelta - hptimetol;
		  
		  if ( (endtime != HPTERROR) && (rec->newstart > endtime) )
		    {
		      if ( verbose > 1 )
			{
			  ms_hptime2seedtimestr (rec->starttime, stime, 1);
			  ms_hptime2seedtimestr (rec->endtime, etime, 1);
			  ms_log (1, "Removing Record %s (%c) :: %s  %s\n",
				  targetsrcname, rec->quality, stime, etime);
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
	       * test for special case of no time coverage (single sample) and no pruning */
	      if ( effstarttime >= (effendtime - hptimetol) &&
		   ! (rec->starttime == rec->endtime &&
		      rec->starttime == effstarttime &&
		      rec->endtime == effendtime) )
		{
		  if ( verbose > 1 )
		    {
		      ms_hptime2seedtimestr (rec->starttime, stime, 1);
		      ms_hptime2seedtimestr (rec->endtime, etime, 1);
		      ms_log (1, "Removing Record %s (%c) :: %s  %s\n",
			      targetsrcname, rec->quality, stime, etime);
		    }
		  
		  rec->flp->recrmcount++;
		  rec->reclen = 0;
		  modcount++;
		}
	      
	    } /* Done pruning at sample level */
	  
	  cmst = cmst->next;
	}
      
      rec = rec->next;
    }
  
  return modcount;
}  /* End of trimtrace() */


/***************************************************************************
 * reconcile_tracetimes():
 *
 * Reconcile the start and end times of the traces in a specified
 * trace group with the list of records in an associated record map.
 * In other words, set the start and end times of each MSTraceSeg in
 * the MSTraceList according to the start time of the first and end
 * time of the last contributing records in the associated record map;
 * this should be performed after the pruning process which could mark
 * complete records as pruned (non-contributing).
 *
 * Returns 0 on success and -1 on error.
 ***************************************************************************/
static int
reconcile_tracetimes (MSTraceList *mstl)
{
  MSTraceID *id;
  MSTraceSeg *seg;
  RecordMap *recmap;
  Record *rec;
  Record *first = 0;
  Record *last = 0;
  
  if ( ! mstl )
    return -1;
  
  if ( ! mstl->traces )
    return -1;
  
  id = mstl->traces;
  while ( id )
    {
      seg = id->first;
      while ( seg )
	{
	  recmap = (RecordMap *) seg->prvtptr;
	  
	  /* Find first contributing record (reclen != 0) */
	  rec = recmap->first;
	  while ( rec )
	    {
	      if ( rec->reclen > 0 )
		{
		  first = rec;
		  break;
		}
	      
	      rec = rec->next;
	    }
	  
	  /* Find last contributing record (reclen != 0) */
	  rec = recmap->last;
	  while ( rec )
	    {
	      if ( rec->reclen > 0 )
		{
		  last = rec;
		  break;
		}
	      
	      rec = rec->prev;
	    }
	  
	  /* Set a new MSTraceSeg start time */
	  if ( first )
	    {
	      /* Use the new boundary start time if set and sane */
	      if ( first->newstart != HPTERROR && first->newstart > first->starttime )
		seg->starttime = first->newstart;
	      /* Otherwise use the record start time */
	      else
		seg->starttime = first->starttime;
	    }
	  
	  /* Set a new MSTraceSeg end time */
	  if ( last )
	    {
	      /* Use the new boundary end time if set and sane */
	      if ( last->newend != HPTERROR && last->newend < last->endtime )
		seg->endtime = last->newend;
	      /* Otherwise use the record end time */
	      else
		seg->endtime = last->endtime;
	    }
	  
	  first = 0;
	  last = 0;
	  seg = seg->next;
	}
      
      id = id->next;
    }
  
  return 0;
}  /* End of reconcile_tracetimes() */


/***************************************************************************
 * qcompare:
 *
 * Compare two different quality codes.
 *
 * Returns:
 * -1 = quality1 is greater than quality2
 *  0 = qualities are equal
 *  1 = quality2 is greater than quality1
 ***************************************************************************/
static int
qcompare (const char quality1, const char quality2)
{
  if ( quality1 == quality2 )
    return 0;
  
  if ( quality1 == 'R' && (quality2 == 'D' || quality2 == 'Q' || quality2 == 'M') )
    return 1;
  
  if ( quality1 == 'D' && (quality2 == 'Q' || quality2 == 'M') )
    return 1;
  
  if ( quality1 == 'Q' && quality2 == 'M' )
    return 1;
  
  return -1;
} /* End of qcompare() */


/***************************************************************************
 * minsegmentlength:
 *
 * Enforce a minimum segment length by removing segments shorter than
 * minseconds from the collection.
 *
 * Returns 0 on success and -1 otherwise.
 ***************************************************************************/
static int
minsegmentlength (MSTraceList *mstl, double minseconds)
{
  MSTraceID *id = 0;
  MSTraceSeg *seg = 0;
  MSTraceSeg *freeseg = 0;
  hptime_t hpminimum;
  hptime_t segmentlength;
  char timestr[50];
  RecordMap *recmap = 0;
  Record *rec;
  Record *recnext;
  
  if ( ! mstl )
    return -1;
  
  hpminimum = minseconds * HPTMODULUS;
  
  /* Loop through trace list */
  id = mstl->traces;  
  while ( id )
    {
      /* Loop through segment list */
      seg = id->first;
      while ( seg )
        {
	  segmentlength = seg->endtime - seg->starttime;
	  
	  if ( segmentlength < hpminimum )
	    {
	      if ( verbose > 2 )
		{
		  ms_hptime2seedtimestr (seg->starttime, timestr, 1);
		  ms_log (1, "Removing segment of %g seconds for %s starting at %s\n",
			  (double)segmentlength / HPTMODULUS, id->srcname, timestr);
		}
	      
	      /* Relink segment list to remove MSTraceSeg */
	      if ( id->first == seg )
		id->first = seg->next;
	      
	      if ( id->last == seg )
		id->last = seg->prev;
	      
	      if ( seg->prev )
		seg->prev->next = seg->next;
	      
	      if ( seg->next )
		seg->next->prev = seg->prev;
	      
	      id->numsegments--;
	      
	      if ( seg->prvtptr )
		{
		  recmap = (RecordMap *) seg->prvtptr;
		  rec = recmap->first;
		  while ( rec )
		    {
		      recnext = rec->next;
		      free (rec);
		      rec = recnext;
		    }
		  free (recmap);
		}
	      
	      freeseg = seg;
	      seg = seg->next;
	      free (freeseg);
	    }
	  else
	    {
	      seg = seg->next;
	    }
	}
      
      id = id->next;
    }
  
  return 0;
} /* End of minsegmentlength() */


/***************************************************************************
 * longestsegmentonly:
 *
 * Determine the longest trace segment for each channel and remove all
 * others.
 *
 * Returns 0 on success and -1 otherwise.
 ***************************************************************************/
static int
longestsegmentonly (MSTraceList *mstl)
{
  MSTraceID *id = 0;
  MSTraceSeg *seg = 0;
  MSTraceSeg *longestseg = 0;
  hptime_t longestsegment;
  hptime_t segmentlength;
  char timestr[50];
  RecordMap *recmap = 0;
  Record *rec;
  Record *recnext;
  
  if ( ! mstl )
    return -1;
  
  /* Loop through trace list */
  id = mstl->traces;
  while ( id )
    {
      /* Shortcut when not multiple segments */
      if ( id->numsegments <= 1 )
	{
	  id = id->next;
	  continue;
	}
      
      longestsegment = 0;
      longestseg = id->first;
      
      /* 1st pass: identify longest segment */
      seg = id->first;
      while ( seg )
        {
	  segmentlength = seg->endtime - seg->starttime;
	  
	  if ( segmentlength > longestsegment )
	    {
	      longestseg = seg;
	      longestsegment = segmentlength;
	    }
	  
          seg = seg->next;
	}
      
      if ( verbose > 2 )
	{
	  ms_hptime2seedtimestr (longestseg->starttime, timestr, 1);
	  ms_log (1, "Removing all but segment of %g seconds for %s starting at %s\n",
		  (double)longestsegment / HPTMODULUS, id->srcname, timestr);
	}
      
      /* 2nd pass: remove all but the longest segment */
      seg = id->first;
      while ( seg )
        {
	  if ( seg != longestseg )
	    {
	      if ( seg->prvtptr )
		{
		  recmap = (RecordMap *) seg->prvtptr;
		  rec = recmap->first;
		  while ( rec )
		    {
		      recnext = rec->next;
		      free (rec);
		      rec = recnext;
		    }
		  free (recmap);
		}
	      
	      free (seg);
	    }
	  
          seg = seg->next;
	}
      
      /* Longest segment is the only segment */
      id->first = longestseg;
      id->last = longestseg;
      id->numsegments = 1;
      longestseg->prev = NULL;
      longestseg->next = NULL;
      
      id = id->next;
    }
  
  return 0;
} /* End of longestsegmentonly() */


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
  
  while ( flp != 0 )
    {
      if ( ! nomods && ! flp->reordercount && ! flp->recrmcount && ! flp->rectrimcount )
	{
	  flp = flp->next;
	  continue;
	}
      
      if ( replaceinput )
	ms_log (0, " Records split: %3d trimmed: %3d removed: %3d, Segments reordered: %3d :: %s\n",
		flp->recsplitcount, flp->rectrimcount, flp->recrmcount, flp->reordercount, flp->outfilename);
      else
	ms_log (0, " Records split: %3d trimmed: %3d removed: %3d, Segments reordered: %3d :: %s\n",
		flp->recsplitcount, flp->rectrimcount, flp->recrmcount, flp->reordercount, flp->infilename);
      
      flp = flp->next;
    }
  
  return;
}  /* End of printmodsummary() */


/***************************************************************************
 * printtracemap():
 *
 * Print record map for each MSTraceSeg to stdout.
 ***************************************************************************/
static void
printtracemap (MSTraceList *mstl)
{
  MSTraceID *id = 0;
  MSTraceSeg *seg = 0;
  char stime[30];
  char etime[30];
  int segcnt = 0;
  
  if ( ! mstl )
    return;
  
  id = mstl->traces;
  
  /* Print out the appropriate header */
  ms_log (0, "\nTrace Map:\n");
  ms_log (0, "   Source              Start sample             End sample        Hz   Samples\n");
  
  while ( id )
    {
      seg = id->first;
      
      while ( seg )
	{
	  /* Create formatted time strings */
	  if ( ms_hptime2seedtimestr (seg->starttime, stime, 1) == NULL )
	    ms_log (2, "Cannot convert trace start time for %s\n", id->srcname);
	  
	  if ( ms_hptime2seedtimestr (seg->endtime, etime, 1) == NULL )
	    ms_log (2, "Cannot convert trace end time for %s\n", id->srcname);
	  
	  /* Print MSTraceSeg header */
	  ms_log (0, "%-15s %-24s %-24s %-4.4g %-d\n",
		  id->srcname, stime, etime, seg->samprate, seg->samplecnt);
	  
	  if ( ! seg->prvtptr )
	    {
	      ms_log (2, "No record map associated with this MSTraceSeg.\n");
	    }
	  else
	    {
	      printrecordmap ((RecordMap *) seg->prvtptr, 0);
	    }
	  
	  segcnt++;
	  seg = seg->next;
	}
      
      id = id->next;
    }
  
  ms_log (0, "End of trace map: %d trace segment(s)\n\n", segcnt);
  
}  /* End of printtracemap() */


/***************************************************************************
 * printrecordmap():
 *
 * Print record map to stdout.
 ***************************************************************************/
static void
printrecordmap (RecordMap *recmap, flag details)
{
  char stime[30];
  char etime[30];
  Record *rec;
  
  if ( ! recmap )
    return;
  
  rec = recmap->first;
  
  ms_log (0, "Record map contains %lld records:\n", recmap->recordcnt);
  
  while ( rec )
    {
      ms_log (0, "  Filename: %s  Offset: %llu  RecLen: %d  Quality: %c\n",
	      rec->flp->infilename, (long long unsigned)rec->offset, rec->reclen, rec->quality);
      
      ms_hptime2seedtimestr (rec->starttime, stime, 1);
      ms_hptime2seedtimestr (rec->endtime, etime, 1);
      ms_log (0, "        Start: %s       End: %s\n", stime, etime);

      if ( details )
	{
	  if ( rec->newstart == HPTERROR ) strcpy (stime, "NONE");
	  else ms_hptime2seedtimestr (rec->newstart, stime, 1);
	  if ( rec->newend == HPTERROR ) strcpy (etime, "NONE" );
	  else ms_hptime2seedtimestr (rec->newend, etime, 1);
	  ms_log (0, "  Start bound: %-24s  End bound: %-24s\n", stime, etime);
	  
	  if ( rec->selectstart == HPTERROR ) strcpy (stime, "NONE");
	  else ms_hptime2seedtimestr (rec->selectstart, stime, 1);
	  if ( rec->selectend == HPTERROR ) strcpy (etime, "NONE" );
	  else ms_hptime2seedtimestr (rec->selectend, etime, 1);
	  ms_log (0, " Select start: %-24s Select end: %-24s\n", stime, etime);
	}
      
      rec = rec->next;
    }
}  /* End of printrecordmap() */


/***************************************************************************
 * printwritten():
 *
 * Print summary of output records.
 ***************************************************************************/
static void
printwritten (MSTraceList *mstl)
{
  MSTraceID *id = 0;
  MSTraceSeg *seg = 0;
  char stime[30];
  char etime[30];
  FILE *ofp;
  
  if ( ! mstl )
    return;
  
  if ( strcmp (writtenfile, "-") == 0 )
    {
      ofp = stdout;
    }
  else if ( strcmp (writtenfile, "--") == 0 )
    {
      ofp = stderr;
    }
  else if ( (ofp = fopen (writtenfile, "ab")) == NULL )
    {
      ms_log (2, "Cannot open output file: %s (%s)\n",
	      writtenfile, strerror(errno));
      return;
    }
  
  /* Loop through trace list */
  id = mstl->traces;  
  while ( id )
    {
      /* Loop through segment list */
      seg = id->first;
      while ( seg )
        {
	  if ( ms_hptime2seedtimestr (seg->starttime, stime, 1) == NULL )
	    ms_log (2, "Cannot convert trace start time for %s\n", id->srcname);
	  
	  if ( ms_hptime2seedtimestr (seg->endtime, etime, 1) == NULL )
	    ms_log (2, "Cannot convert trace end time for %s\n", id->srcname);
	  
	  fprintf (ofp, "%s%s|%s|%s|%s|%c|%-24s|%-24s|%lld|%lld\n",
		   ( writtenprefix )? writtenprefix : "",
		   id->network, id->station, id->location, id->channel, id->dataquality,
		   stime, etime, (long long int) *((int64_t *)seg->prvtptr),
		   (long long int)seg->samplecnt);
	  
	  seg = seg->next;
	}
      
      id = id->next;
    }
  
  if ( ofp != stdout && fclose (ofp) )
    ms_log (2, "Cannot close output file: %s (%s)\n",
	    writtenfile, strerror(errno));
  
}  /* End of printwritten() */


/***************************************************************************
 * findselectlimits():
 *
 * Determine selection limits for the given record based on all
 * matching selection entries.
 *
 * Return 0 on success and -1 on error.
 ***************************************************************************/
static int
findselectlimits (Selections *select, char *srcname, hptime_t starttime,
		  hptime_t endtime, Record *rec)
{
  SelectTime *selecttime;
  
  if ( ! rec || ! srcname || ! select )
    return -1;
  
  while ( (select = ms_matchselect(select, srcname, starttime, endtime, &selecttime)) )
    {
      while ( selecttime )
	{
	  if ( (starttime < selecttime->starttime && ! (starttime <= selecttime->starttime && endtime >= selecttime->starttime)) )
	    { selecttime = selecttime->next; continue; }
	  else if ( (endtime > selecttime->endtime && ! (starttime <= selecttime->endtime && endtime >= selecttime->endtime)) )
	    { selecttime = selecttime->next; continue; }
	  
	  if ( rec->selectstart == HPTERROR || rec->selectstart > selecttime->starttime )
	    rec->selectstart = selecttime->starttime;
	  
	  if ( rec->selectend == HPTERROR || rec->selectend < selecttime->endtime )
	    rec->selectend = selecttime->endtime;
	  
	  /* Shortcut if the entire record is already selected */
	  if ( rec->starttime >= rec->selectstart && rec->endtime <= rec->selectend )
	    return 0;
	  
	  selecttime = selecttime->next;
	}
      
      select = select->next;
    }
  
  return 0;
}  /* End of findselectlimits() */


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
  
  if ( ! recmap )
    return -1;
  
  if ( ! recmap->first || ! recmap->last ) /* Done if empty */
    return 0;
  
  top = recmap->first;
  insize = 1;
  
  for (;;)
    {
      p = top;
      top = NULL;
      tail = NULL;
      
      nmerges = 0;  /* count number of merges we do in this pass */
      
      while ( p )
        {
          nmerges++;  /* there exists a merge to be done */
          
          /* step `insize' places along from p */
          q = p;
          psize = 0;
          for (i = 0; i < insize; i++)
            {
              psize++;
              q = q->next;
              if ( ! q )
                break;
            }
          
          /* if q hasn't fallen off end, we have two lists to merge */
          qsize = insize;
          
          /* now we have two lists; merge them */
          while ( psize > 0 || (qsize > 0 && q) )
            {
              /* decide whether next element of merge comes from p or q */
              if ( psize == 0 )
                {  /* p is empty; e must come from q. */
                  e = q; q = q->next; qsize--;
                }
              else if ( qsize == 0 || ! q )
                {  /* q is empty; e must come from p. */
                  e = p; p = p->next; psize--;
                }
              else if ( recordcmp (p, q) <= 0 )
                {  /* First element of p is lower (or same), e must come from p. */
                  e = p; p = p->next; psize--;
                }
              else
                {  /* First element of q is lower; e must come from q. */
                  e = q; q = q->next; qsize--;
                }
              
              /* add the next element to the merged list */
              if ( tail )
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
      if ( nmerges <= 1 )   /* allow for nmerges==0, the empty list case */
        {
          recmap->first = top;
          recmap->last = tail;
	  
          return 0;
        }
      
      /* Otherwise repeat, merging lists twice the size */
      insize *= 2;
    }
}  /* End of sortrecmap() */


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
  hptime_t *start1;
  hptime_t *start2;
  
  if ( ! rec1 || ! rec2 )
    return -1;
  
  /* Determine effective start times */
  start1 = ( rec1->newstart != HPTERROR ) ? &(rec1->newstart) : &(rec1->starttime);
  start2 = ( rec2->newstart != HPTERROR ) ? &(rec2->newstart) : &(rec2->starttime);
  
  if ( *start1 > *start2 )
    {
      return 1;
    }
  
  return 0;
}  /* End of recordcmp() */


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
  char *stagelengthstr = 0;
  char *selectfile = 0;
  char *matchpattern = 0;
  char *rejectpattern = 0;
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
	  timetol = strtod (getoptval(argcount, argvec, optind++), NULL);
	}
      else if (strcmp (argvec[optind], "-rt") == 0)
	{
	  sampratetol = strtod (getoptval(argcount, argvec, optind++), NULL);
	}
      else if (strcmp (argvec[optind], "-E") == 0)
	{
	  bestquality = 0;
	}
      else if (strcmp (argvec[optind], "-sb") == 0)
	{
	  stagelengthstr = getoptval(argcount, argvec, optind++);
	}
      else if (strcmp (argvec[optind], "-s") == 0)
	{
	  selectfile = getoptval(argcount, argvec, optind++);
	}
      else if (strcmp (argvec[optind], "-ts") == 0)
	{
	  starttime = ms_seedtimestr2hptime (getoptval(argcount, argvec, optind++));
	  if ( starttime == HPTERROR )
	    return -1;
	}
      else if (strcmp (argvec[optind], "-te") == 0)
	{
	  endtime = ms_seedtimestr2hptime (getoptval(argcount, argvec, optind++));
	  if ( endtime == HPTERROR )
	    return -1;
	}
      else if (strcmp (argvec[optind], "-M") == 0)
	{
	  matchpattern = strdup(getoptval(argcount, argvec, optind++));
	}
      else if (strcmp (argvec[optind], "-R") == 0)
	{
	  rejectpattern = strdup(getoptval(argcount, argvec, optind++));
	}
      else if (strcmp (argvec[optind], "-szs") == 0)
	{
	  skipzerosamps = 1;
	}
      else if (strcmp (argvec[optind], "-lso") == 0)
	{
	  longestonly = 1;
	}
      else if (strcmp (argvec[optind], "-msl") == 0)
	{
	  minseglength = strtod (getoptval(argcount, argvec, optind++), NULL);
	}
      else if (strcmp (argvec[optind], "-m") == 0)
	{
          tptr = getoptval(argcount, argvec, optind++);
	  
	  if ( ms_addselect (&selections, tptr, HPTERROR, HPTERROR) < 0 )
	    {
	      ms_log (2, "Unable to add selection: '%s'\n", tptr);
	      return -1;
	    }
	}
      else if (strcmp (argvec[optind], "-rep") == 0)
        {
          replaceinput = 1;
        }
      else if (strcmp (argvec[optind], "-nb") == 0)
        {
          nobackups = 1;
        }
      else if (strcmp (argvec[optind], "-o") == 0)
        {
          outputfile = getoptval(argcount, argvec, optind++);
	  outputmode = 0;
        }
      else if (strcmp (argvec[optind], "+o") == 0)
        {
          outputfile = getoptval(argcount, argvec, optind++);
	  outputmode = 1;
        }
      else if (strcmp (argvec[optind], "-A") == 0)
        {
          if ( addarchive(getoptval(argcount, argvec, optind++), NULL) == -1 )
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
          tptr = getoptval(argcount, argvec, optind++);
          restampqind = *tptr;
          
          if ( ! MS_ISDATAINDICATOR (restampqind) )
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
	  writtenfile = getoptval(argcount, argvec, optind++);
        }
      else if (strcmp (argvec[optind], "-outprefix") == 0)
        {
	  writtenprefix = getoptval(argcount, argvec, optind++);
        }
      else if (strcmp (argvec[optind], "-CHAN") == 0)
        {
          if ( addarchive(getoptval(argcount, argvec, optind++), CHANLAYOUT) == -1 )
            return -1;
        }
      else if (strcmp (argvec[optind], "-QCHAN") == 0)
        {
          if ( addarchive(getoptval(argcount, argvec, optind++), QCHANLAYOUT) == -1 )
            return -1;
        }
      else if (strcmp (argvec[optind], "-CDAY") == 0)
        {
          if ( addarchive(getoptval(argcount, argvec, optind++), CDAYLAYOUT) == -1 )
            return -1;
        }
      else if (strcmp (argvec[optind], "-SDAY") == 0)
        {
          if ( addarchive(getoptval(argcount, argvec, optind++), SDAYLAYOUT) == -1 )
            return -1;
        }
      else if (strcmp (argvec[optind], "-BUD") == 0)
        {
          if ( addarchive(getoptval(argcount, argvec, optind++), BUDLAYOUT) == -1 )
            return -1;
        }
      else if (strcmp (argvec[optind], "-CSS") == 0)
        {
          if ( addarchive(getoptval(argcount, argvec, optind++), CSSLAYOUT) == -1 )
            return -1;
        }
      else if (strncmp (argvec[optind], "-", 1) == 0 &&
	       strlen (argvec[optind]) > 1 )
	{
	  ms_log (2, "Unknown option: %s\n", argvec[optind]);
	  exit (1);
	}
      else
	{
	  tptr = argvec[optind];
	  
          /* Check for an input file list */
          if ( tptr[0] == '@' )
            {
              if ( addlistfile (tptr+1) < 0 )
                {
                  ms_log (2, "Error adding list file %s", tptr+1);
                  exit (1);
                }
            }
          /* Otherwise this is an input file */
          else
            {
	      /* Add file to global file list */
	      if ( addfile (tptr) )
		{
                  ms_log (2, "Error adding file to input list %s", tptr);
                  exit (1);
		}
	    }
	}
    }
  
  /* Make sure input file(s) were specified */
  if ( filelist == 0 )
    {
      ms_log (2, "No input files were specified\n\n");
      ms_log (1, "%s version %s\n\n", PACKAGE, VERSION);
      ms_log (1, "Try %s -h for usage\n", PACKAGE);
      exit (0);
    }

  /* Make sure output file(s) were specified or replacing originals */
  if ( archiveroot == 0 && outputfile == 0 && replaceinput == 0 )
    {
      ms_log (2, "No output files were specified\n\n");
      ms_log (1, "%s version %s\n\n", PACKAGE, VERSION);
      ms_log (1, "Try %s -h for usage\n", PACKAGE);
      exit (0);
    }
  
  /* The no backups option cannot be used when not replacing the originals */
  if ( nobackups && ! replaceinput )
    nobackups = 0;
  
  /* Read data selection file */
  if ( selectfile )
    {
      if ( ms_readselectionsfile (&selections, selectfile) < 0 )
	{
	  ms_log (2, "Cannot read data selection file\n");
	  exit (1);
	}
    }
  
  /* Expand match pattern from a file if prefixed by '@' */
  if ( matchpattern )
    {
      if ( *matchpattern == '@' )
	{
	  tptr = strdup(matchpattern + 1); /* Skip the @ sign */
	  free (matchpattern);
	  matchpattern = 0;
	  
	  if ( readregexfile (tptr, &matchpattern) <= 0 )
	    {
	      ms_log (2, "Cannot read match pattern regex file\n");
	      exit (1);
	    }
	  
	  free (tptr);
	}
    }
  
  /* Expand reject pattern from a file if prefixed by '@' */
  if ( rejectpattern )
    {
      if ( *rejectpattern == '@' )
	{
	  tptr = strdup(rejectpattern + 1); /* Skip the @ sign */
	  free (rejectpattern);
	  rejectpattern = 0;
	  
	  if ( readregexfile (tptr, &rejectpattern) <= 0 )
	    {
	      ms_log (2, "Cannot read reject pattern regex file\n");
	      exit (1);
	    }
	  
	  free (tptr);
	}
    }
  
  /* Compile match and reject patterns */
  if ( matchpattern )
    {
      if ( ! (match = (regex_t *) malloc (sizeof(regex_t))) )
	{
	  ms_log (2, "Cannot allocate memory for match expression\n");	  
	  exit (1);
	}
      
      if ( regcomp (match, matchpattern, REG_EXTENDED) != 0)
	{
	  ms_log (2, "Cannot compile match regex: '%s'\n", matchpattern);
	}
      
      free (matchpattern);
    }
  
  if ( rejectpattern )
    {
      if ( ! (reject = (regex_t *) malloc (sizeof(regex_t))) )
	{
	  ms_log (2, "Cannot allocate memory for reject expression\n");
	  exit (1);
	}
      
      if ( regcomp (reject, rejectpattern, REG_EXTENDED) != 0)
	{
	  ms_log (2, "Cannot compile reject regex: '%s'\n", rejectpattern);
	}
      
      free (rejectpattern);
    }
  
  /* Allocate record staging buffer */
  if ( stagelengthstr )
    {
      /* Parse size from argument */
      if ( ! (stagelength = calcsize(stagelengthstr)) )
	{
	  ms_log (2, "Cannot parse staging buffer length: '%s'\n", stagelengthstr);
	  exit (1);
	}
      
      /* Parse buffer */
      if ( ! (stagebuffer = (char *) malloc (stagelength)) )
	{
	  ms_log (2, "Cannot allocate memory for staging buffer (%lld bytes requested)\n",
		  (long long int) stagelength);
	  exit (1);
	}
      
      if ( verbose > 1 )
	ms_log (1, "Allocated staging buffer (%lld bytes)\n",
		(long long int) stagelength);
    }
  
  /* Report the program version */
  if ( verbose )
    ms_log (1, "%s version: %s\n", PACKAGE, VERSION);
  
  return 0;
}  /* End of processparam() */


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
  if ( argvec == NULL || argvec[argopt] == NULL ) {
    ms_log (2, "getoptval(): NULL option requested\n");
    exit (1);
    return 0;
  }
  
  /* Special case of '-o -' usage */
  if ( (argopt+1) < argcount && strcmp (argvec[argopt], "-o") == 0 )
    if ( strcmp (argvec[argopt+1], "-") == 0 )
      return argvec[argopt+1];

  /* Special case of '+o -' usage */
  if ( (argopt+1) < argcount && strcmp (argvec[argopt], "+o") == 0 )
    if ( strcmp (argvec[argopt+1], "-") == 0 )
      return argvec[argopt+1];
  
  /* Special case of '-s -' usage */
  if ( (argopt+1) < argcount && strcmp (argvec[argopt], "-s") == 0 )
    if ( strcmp (argvec[argopt+1], "-") == 0 )
      return argvec[argopt+1];

  /* Special case of '-out -' or '-out --' usage */
  if ( (argopt+1) < argcount && strcmp (argvec[argopt], "-out") == 0 )
    if ( strcmp (argvec[argopt+1], "-") == 0 ||
	 strcmp (argvec[argopt+1], "--") == 0)
      return argvec[argopt+1];
  
  if ( (argopt+1) < argcount && *argvec[argopt+1] != '-' )
    return argvec[argopt+1];
  
  ms_log (2, "Option %s requires a value, try -h for usage\n", argvec[argopt]);
  exit (1);
  return 0;
}  /* End of getoptval() */


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
  if ( getrlimit (RLIMIT_NOFILE, &rlim) == -1 )
    {
      ms_log (2, "getrlimit() failed to get open file limit\n");
      return -1;
    }
  
  if ( rlim.rlim_cur < limit )
    {
      oldlimit = rlim.rlim_cur;
      rlim.rlim_cur = (rlim_t) limit;
      
      if ( verbose > 1 )
	ms_log (1, "Setting open file limit to %d\n",
		(int) rlim.rlim_cur);
      
      if ( setrlimit (RLIMIT_NOFILE, &rlim) == -1 )
	{
	  ms_log (2, "setrlimit failed to raise open file limit from %d to %d (max: %d)\n",
		  (int) oldlimit, (int) limit, rlim.rlim_max);
	  return -1;
	}
    }
  
  return (int) rlim.rlim_cur;
}  /* End of setofilelimit() */


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
  
  if ( ! filename )
    {
      ms_log (2, "addfile(): No file name specified\n");
      return -1;
    }
  
  if ( ! (newlp = (Filelink *) calloc (1, sizeof (Filelink))) )
    {
      ms_log (2, "addfile(): Cannot allocate memory\n");
      return -1;
    }
  
  /* Check for optional read byte range specifiers
   * Expected form: "filename@startoffset:endoffset"
   * Both start are optional */
  if ( (at = strrchr (filename, '@')) )
    {
      *at++ = '\0';
      
      if ( (colon = strrchr (at, ':')) )
        {
          *colon++ = '\0';
          newlp->endoffset = strtoull (colon, NULL, 10);
        }
      
      newlp->startoffset = strtoull (at, NULL, 10);
    }
  
  if ( ! (newlp->infilename = strdup(filename)) )
    {
      ms_log (2, "addfile(): Cannot duplicate string\n");
      return -1;
    }
  
  /* Add new file to the end of the list */
  if ( filelisttail == 0 )
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
}  /* End of addfile() */


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
  
  if ( verbose >= 1 )
    ms_log (1, "Reading list file '%s'\n", filename);
  
  if ( ! (fp = fopen(filename, "rb")) )
    {
      ms_log (2, "Cannot open list file %s: %s\n", filename, strerror(errno));
      return -1;
    }
  
  while ( fgets (filelistent, sizeof(filelistent), fp) )
    {
      char *cp;
      
      /* End string at first newline character */
      if ( (cp = strchr(filelistent, '\n')) )
        *cp = '\0';
      
      /* Skip empty lines */
      if ( ! strlen (filelistent) )
        continue;
      
      /* Skip comment lines */
      if ( *filelistent == '#' )
        continue;
      
      if ( verbose > 1 )
	ms_log (1, "Adding '%s' from list file\n", filelistent);
      
      if ( addfile (filelistent) )
	return -1;
      
      filecount++;
    }
  
  fclose (fp);
  
  return filecount;
}  /* End of addlistfile() */


/***************************************************************************
 * addarchive:
 * Add entry to the data stream archive chain.  'layout' if defined
 * will be appended to 'path'.
 *
 * Returns 0 on success, and -1 on failure
 ***************************************************************************/
static int
addarchive ( const char *path, const char *layout )
{
  Archive *newarch;
  int pathlayout;
  
  if ( ! path )
    {
      ms_log (2, "addarchive(): cannot add archive with empty path\n");
      return -1;
    }
  
  if ( ! (newarch = (Archive *) malloc (sizeof (Archive))) )
    {
      ms_log (2, "addarchive(): cannot allocate memory for new archive definition\n");
      return -1;
    }
  
  /* Setup new entry and add it to the front of the chain */
  pathlayout = strlen (path) + 2;
  if ( layout )
    pathlayout += strlen (layout);
  
  if ( ! (newarch->datastream.path = (char *) malloc (pathlayout)) )
    {
      ms_log (2, "addarchive(): cannot allocate memory for new archive path\n");
      if ( newarch )
        free (newarch);
      return -1;
    }
  
  if ( layout )
    snprintf (newarch->datastream.path, pathlayout, "%s/%s", path, layout);
  else
    snprintf (newarch->datastream.path, pathlayout, "%s", path);
  
  newarch->datastream.idletimeout = 60;
  newarch->datastream.grouproot = NULL;
  
  newarch->next = archiveroot;
  archiveroot = newarch;
  
  return 0;
}  /* End of addarchive() */


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
  char  line[1024];
  char  linepattern[1024];
  int   regexcnt = 0;
  int   lengthbase;
  int   lengthadd;
  
  if ( ! regexfile )
    {
      ms_log (2, "readregexfile: regex file not supplied\n");
      return -1;
    }
  
  if ( ! pppattern )
    {
      ms_log (2, "readregexfile: pattern string buffer not supplied\n");
      return -1;
    }
  
  /* Open the regex list file */
  if ( (fp = fopen (regexfile, "rb")) == NULL )
    {
      ms_log (2, "Cannot open regex list file %s: %s\n",
	      regexfile, strerror (errno));
      return -1;
    }
  
  if ( verbose )
    ms_log (1, "Reading regex list from %s\n", regexfile);
  
  *pppattern = NULL;
  
  while ( (fgets (line, sizeof(line), fp)) !=  NULL)
    {
      /* Trim spaces and skip if empty lines */
      if ( sscanf (line, " %s ", linepattern) != 1 )
	continue;
      
      /* Skip comment lines */
      if ( *linepattern == '#' )
	continue;
      
      regexcnt++;
      
      /* Add regex to compound regex */
      if ( *pppattern )
	{
	  lengthbase = strlen(*pppattern);
	  lengthadd = strlen(linepattern) + 4; /* Length of addition plus 4 characters: |()\0 */
	  
	  *pppattern = realloc (*pppattern, lengthbase + lengthadd);
	  
	  if ( *pppattern )
	    {
	      snprintf ((*pppattern)+lengthbase, lengthadd, "|(%s)", linepattern);
	    }
	  else
	    {
	      ms_log (2, "Cannot allocate memory for regex string\n");
	      return -1;
	    }
	}
      else
	{
	  lengthadd = strlen(linepattern) + 3; /* Length of addition plus 3 characters: ()\0 */
	  
	  *pppattern = malloc (lengthadd);
	  
	  if ( *pppattern )
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
}  /* End of readregexfile() */


/***************************************************************************
 * calcsize:
 * 
 * Calculate a size in bytes for the specified size string.  If the
 * string is terminated with the following suffixes the specified
 * scaling will be applied:
 *
 * 'K' or 'k' : kilobytes - value * 1024
 * 'M' or 'm' : megabytes - value * 1024*1024
 * 'G' or 'g' : gigabytes - value * 1024*1024*1024
 *
 * Returns a size in bytes on success and 0 on error.
 ***************************************************************************/
static off_t
calcsize (char *sizestr)
{
  off_t size = 0;
  char *parsestr;
  int termchar;
  
  if ( ! sizestr )
    return 0;

  if ( ! (parsestr = strdup (sizestr)) )
    return 0;
  
  termchar = strlen (parsestr) - 1;
  
  if ( termchar <= 0 )
    return 0;
  
  /* For kilobytes */
  if ( parsestr[termchar] == 'K' || parsestr[termchar] == 'k' )
    {
      parsestr[termchar] = '\0';
      size = strtoull (parsestr, NULL, 10);
      if ( ! size )
        {
	  ms_log (1, "calcsize(): Error converting %s to integer", parsestr);
          return 0;
        }
      size *= 1024;
    }  
  /* For megabytes */
  else if ( parsestr[termchar] == 'M' || parsestr[termchar] == 'm' )
    {
      parsestr[termchar] = '\0';
      size = strtoull (parsestr, NULL, 10);
      if ( ! size )
        {
	  ms_log (1, "calcsize(): Error converting %s to integer", parsestr);
          return 0;
        }
      size *= 1024*1024;
    }
  /* For gigabytes */
  else if ( parsestr[termchar] == 'G' || parsestr[termchar] == 'g' )
    {
      parsestr[termchar] = '\0';
      size = strtoull (parsestr, NULL, 10);
      if ( ! size )
        {
	  ms_log (1, "calcsize(): Error converting %s to integer", parsestr);
          return 0;
        }
      size *= 1024*1024*1024;
    }
  else
    {
      size = strtoull(parsestr, NULL, 10);
      if ( ! size )
        {
	  ms_log (1, "calcsize(): Error converting %s to integer", parsestr);
	  return 0;
        }
    }
  
  if ( parsestr )
    free (parsestr);
  
  return size;
}  /* End of calcsize() */


/***************************************************************************
 * usage():
 * Print the usage message.
 ***************************************************************************/
static void
usage (int level)
{
  fprintf (stderr, "%s - select, sort and prune Mini-SEED: %s\n\n", PACKAGE, VERSION);
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
	   " -sb size     Use a record staging buffer of specified size for speed\n"
	   "\n"
	   " ## Data selection options ##\n"
	   " -s file      Specify a file containing selection criteria\n"
	   " -ts time     Limit to records that contain or start after time\n"
	   " -te time     Limit to records that contain or end before time\n"
	   "                time format: 'YYYY[,DDD,HH,MM,SS,FFFFFF]' delimiters: [,:.]\n"
	   " -M match     Limit to records matching the specified regular expression\n"
	   " -R reject    Limit to records not matching the specfied regular expression\n"
	   "                Regular expressions are applied to: 'NET_STA_LOC_CHAN_QUAL'\n"
	   " -szs         Skip input records that contain zero samples\n"
	   " -lso         Longest segment only, output only longest segment per channel\n"
	   " -msl secs    Minimum segment length, skip any shorter segments\n"
	   "\n"
	   " ## Output options ##\n"
	   " -rep         Replace input files, creating backup (.orig) files\n"
	   " -nb          Do not keep backups of original input files if replacing them\n"
	   " -o file      Specify a single output file, use +o file to append\n"
	   " -A format    Write all records in a custom directory/file layout (try -H)\n"
	   " -Pr          Prune data at the record level using 'best' quality priority\n"
	   " -Ps          Prune data at the sample level using 'best' quality priority\n"
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
	   " file#        Files(s) of Mini-SEED records\n"
	   "\n");
  
  if  ( level )
    {
      fprintf (stderr,
               "\n"
	       "  # Preset format layouts #\n"
	       " -CHAN dir    Write records into separate Net.Sta.Loc.Chan files\n"
	       " -QCHAN dir   Write records into separate Net.Sta.Loc.Chan.Quality files\n"
	       " -CDAY dir    Write records into separate Net.Sta.Loc.Chan.Year:Yday:<time> files\n"
	       " -SDAY dir    Write records into separate Net.Sta.Year:Yday files\n"
	       " -BUD BUDdir  Write records in a BUD file layout\n"
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
}  /* End of usage() */
