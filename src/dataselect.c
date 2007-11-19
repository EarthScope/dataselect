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
 * modified 2007.323
 ***************************************************************************/

/***************************************************************************
 *
 * Data structures and operational overview
 *
 * The data map (using actual structure names):
 *
 * MSTraceGroup
 *   |-MSTrace
 *   |   |-RecordMap
 *   |        |-Record
 *   |        |-Record
 *   |        |-...
 *   |
 *   |-MSTrace
 *   |   |-RecordMap
 *   |        |-Record
 *   |        |-Record
 *   |        |-...
 *   |-...
 *
 * The program goes through the following stages:
 * 
 * 1) Read all input files constructing a data map of contiguous trace
 * segments and the data records that comprise them.  This operation
 * is done by reading each file record-by-record and as each record is
 * read, search the MSTraceGroup for an MSTrace that the record "fits"
 * with (same channel and adjacent in time).  If the record is found
 * to fit with a specific MSTrace its coverage will be added to the
 * MSTrace information otherwise a new MSTrace is created and added to
 * the MSTraceGroup.
 *
 * Each MSTrace has an associated RecordMap which includes a list of
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

#define VERSION "0.9.3"
#define PACKAGE "dataselect"

/* For a linked list of strings, as filled by strparse() */
typedef struct StrList_s {
  char             *element;
  struct StrList_s *next;
} StrList;

/* POD request record containers */
typedef struct ReqRec_s {
  struct StrList_s *strlist;
  char  *network;
  char  *station;
  char  *location;
  char  *channel;
  time_t datastart;
  time_t dataend;
  char  *filename;
  char  *headerdir;
  time_t reqstart;
  time_t reqend;
  char   pruned;
  struct ReqRec_s *prev;
  struct ReqRec_s *next;
} ReqRec;

/* Input/output file information containers */
typedef struct Filelink_s {
  struct ReqRec_s *reqrec;
  char    *infilename;     /* Input file name */
  FILE    *infp;           /* Input file descriptor */
  char    *outfilename;    /* Output file name */
  FILE    *outfp;          /* Output file descriptor */
  int      reordercount;   /* Number of records re-ordered */
  int      recsplitcount;  /* Number of records split */
  int      recrmcount;     /* Number of records removed */
  int      rectrimcount;   /* Number of records trimed */
  hptime_t earliest;       /* Earliest data time in this file */
  hptime_t latest;         /* Latest data time in this file */
  int      byteswritten;   /* Number of bytes written out */
  struct Filelink_s *next;
} Filelink;

/* Archive (output structure) definition containers */
typedef struct Archive_s {
  DataStream  datastream;
  struct Archive_s *next;
} Archive;

/* Mini-SEED record information structures */
typedef struct Record_s {
  struct Filelink_s *flp;
  off_t     offset;
  int       reclen;
  hptime_t  starttime;
  hptime_t  endtime;
  char      quality;
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


static int processpod (char *requestfile, char *datadir);
static int setofilelimit (int limit);
static ReqRec *readreqfile (char *requestfile);
static int writereqfile (char *requestfile, ReqRec *rrlist);

static int processtraces (MSTraceGroup *mstg, Filelink *filelist);
static int writetraces (MSTraceGroup *mstg, Filelink *filelist);
static int trimrecord (Record *rec, char *recbuf);
static void record_handler (char *record, int reclen, void *handlerdata);

static int prunetraces (MSTraceGroup *mstg);
static int getcoverage (MSTraceGroup *mstg, MSTrace *targetmst,
			MSTraceGroup **ppcoverage);
static int trimtrace (MSTrace *targetmst, MSTraceGroup *coverage);
static int reconcile_tracetimes (MSTraceGroup *mstg);
static int qcompare (const char quality1, const char quality2);

static int readfiles (Filelink *flp, MSTraceGroup **ppmstg);
static MSTraceGroup *reinitgroup (MSTraceGroup *mstg);
static void printmodsummary (Filelink *filelist, flag nomods);
static void printtracemap (MSTraceGroup *mstg);
static void printrecordmap (RecordMap *recmap, flag details);

static int processparam (int argcount, char **argvec);
static char *getoptval (int argcount, char **argvec, int argopt);
static int strparse (const char *string, const char *delim, StrList **list);
static void addfile (Filelink **ppfilelist, char *filename, ReqRec *reqrec);
static int  addarchive (const char *path, const char *layout);
static int readregexfile (char *regexfile, char **pppattern);
static void freereqrec (ReqRec *reqrec);
static void freefilelist (Filelink **ppfilelist);
static void usage (int level);

static flag     verbose       = 0;
static flag     basicsum      = 0;    /* Controls printing of basic summary */
static flag     bestquality   = 1;    /* Use D, R, Q, M quality to retain the "best" data when pruning */
static flag     prunedata     = 0;    /* Prune data: 'r= record level, 's' = sample level */
static double   timetol       = -1.0; /* Time tolerance for continuous traces */
static double   sampratetol   = -1.0; /* Sample rate tolerance for continuous traces */
static char     restampqind   = 0;    /* Re-stamp data record/quality indicator */
static int      reclen        = -1;   /* Input data record length, autodetected in most cases */
static flag     modsummary    = 0;    /* Print modification summary after all processing */
static hptime_t starttime     = HPTERROR;  /* Limit to records after starttime */
static hptime_t endtime       = HPTERROR;  /* Limit to records before endtime */
static char     splitboundary = 0;    /* Split records on day(d), hour(h) or minute(m) boundaries */

static regex_t *match         = 0;    /* Compiled match regex */
static regex_t *reject        = 0;    /* Compiled reject regex */

static char    *podreqfile    = 0;    /* POD h. request file */
static char    *poddatadir    = 0;    /* POD data directory */

static flag     replaceinput  = 0;    /* Replace input files */
static flag     nobackups     = 0;    /* Remove re-named original files when done with them */
static char    *outputfile    = 0;    /* Single output file */
static Archive *archiveroot   = 0;    /* Output file structures */

static char     recordbuf[16384];     /* Global record buffer */

static Filelink *gfilelist    = 0;    /* List of input files */


int
main ( int argc, char **argv )
{
  MSTraceGroup *mstg = 0;
  
  /* Set default error message prefix */
  ms_loginit (NULL, NULL, NULL, "ERROR: ");
  
  /* Process input parameters */
  if ( processparam (argc, argv) < 0 )
    return 1;
  
  if ( podreqfile && poddatadir )
    {
      if ( verbose > 2 )
	ms_log (1, "Procesing POD structure:\nrequest file: %s, data dir: %s\n",
		podreqfile, poddatadir);
      
      if ( processpod (podreqfile, poddatadir) )
	{
	  ms_log (2, "Cannot process POD structure\n");
	  return 1;
	}
    }
  else
    {
      if ( verbose > 2 )
	ms_log (1, "Processing input files\n");
      
      /* Read and process all files specified on the command line */
      if ( readfiles (gfilelist, &mstg) )
	return 1;
      
      if ( processtraces (mstg, gfilelist) )
	return 1;
      
      if ( modsummary )
	printmodsummary (gfilelist, verbose);
      
      /* Clean up MSTraceGroup and file list */      
      mst_freegroup (&mstg);      
      freefilelist (&gfilelist);
    }
  
  return 0;
}  /* End of main() */


/***************************************************************************
 * processpod:
 *
 * Process data from a POD structure.  Group input data files by
 * channel and prune each group separately using a M > Q > D > R priority.
 *
 * Returns 0 on success and -1 on error.
 ***************************************************************************/
static int
processpod (char *requestfile, char *datadir)
{
  Filelink *filelist = 0;
  MSTraceGroup *mstg = 0;
  ReqRec *reqrecs;
  ReqRec *fox, *hound;
  ReqRec *reciter, *recnext;
  Filelink *flp;
  int filecount;

  char tmpfilename[1024];
  
  /* Read the request file */
  if ( ! (reqrecs = readreqfile(requestfile)) )
    return -1;
  
  /* Group the input files by channel (complete matching NSLC) and
   * prune as a group.
   *
   * All files in the request file will be pruned, grouping by channel
   * is a minor memory footprint optimization.
   */
  
  hound = reqrecs;
  while ( hound )
    {
      if ( hound->pruned )
	{
	  hound = hound->next;
	  continue;
	}
      
      /* Build appropriate data file name using base data dir and record file name */
      snprintf (tmpfilename, sizeof(tmpfilename), "%s/%s/%s",
		poddatadir, hound->station, hound->filename);
      
      /* Add file to list to be pruned and mark it as pruned */
      addfile (&filelist, tmpfilename, hound);
      hound->pruned = 1;
      filecount = 1;
      
      /* Find all the other matching NSLC files and add them to the list */
      fox = reqrecs;
      while ( fox )
	{
	  if ( fox == hound || fox->pruned )
	    {
	      fox = fox->next;
	      continue;
	    }
	  
	  if ( ! strcmp (hound->network, fox->network) &&
	       ! strcmp (hound->station, fox->station) &&
	       ! strcmp (hound->location, fox->location) &&
	       ! strcmp (hound->channel, fox->channel) )
	    {
	      /* Build appropriate data file name using base data dir and record file name */
	      snprintf (tmpfilename, sizeof(tmpfilename), "%s/%s/%s",
			poddatadir, fox->station, fox->filename);
	      
	      /* Add file to list to be pruned and mark it as pruned */
	      addfile (&filelist, tmpfilename, fox);
	      fox->pruned = 1;
	      filecount++;
	    }
	  
	  fox = fox->next;
	}
      
      /* Increase open file limit if necessary, in general
       * we need 2 X the filecount and some wiggle room. */
      if ( setofilelimit ((filecount * 2) + 20) == -1 )
	{
	  freefilelist (&filelist);
	  hound = hound->next;
	  continue;
	}
      
      /* Read data files & prune time coverage */
      if ( readfiles (filelist, &mstg) )
	return -1;
      
      if ( processtraces (mstg, filelist) )
	return -1;
      
      /* Update the time values in the request records */
      flp = filelist;
      while ( flp )
	{
	  if ( flp->byteswritten == 0 )
	    {
	      /* This signals no coverage, request record will be removed below */
	      flp->reqrec->datastart = flp->reqrec->dataend = 0;
	    }
	  else
	    {
	      flp->reqrec->datastart = (time_t) (MS_HPTIME2EPOCH(flp->earliest));
	      flp->reqrec->dataend = (time_t) (MS_HPTIME2EPOCH(flp->latest));
	    }
	  
	  flp = flp->next;
	}
      
      if ( modsummary )
	printmodsummary (filelist, verbose);
      
      /* Clean up MSTraceGroup and file list */
      mst_freegroup (&mstg);
      freefilelist (&filelist);
      
      hound = hound->next;
    }
  
  if ( verbose > 1 )
    ms_log (1, "Renaming and rewriting request file (h.)\n");
  
  /* Remove request records that have no time coverage, i.e. all data pruned */
  reciter = reqrecs;
  while ( reciter )
    {
      recnext = reciter->next;
      
      if ( reciter->datastart == 0 && reciter->dataend == 0 )
	{
	  if ( reciter == reqrecs )  /* First in chain */
	    reqrecs = reciter->next;
	  else if ( reciter->next == 0 )  /* Last in chain */
	    reciter->prev->next = 0;
	  else  /* In middle of chain */
	    {
	      reciter->prev->next = reciter->next;
	      reciter->next->prev = reciter->prev;
	    }
	  
	  freereqrec (reciter);
	}
      
      reciter = recnext;
    }
  
  /* Rename original request file (add '.orig') */
  snprintf (tmpfilename, sizeof(tmpfilename), "%s.orig", requestfile);
  
  if ( rename (requestfile, tmpfilename) )
    {
      ms_log (2, "cannot rename original request file: %s -> %s : '%s'\n",
	      tmpfilename, requestfile, strerror(errno));
      return -1;
    }
  
  /* Write the request file */
  if ( writereqfile (requestfile, reqrecs) )
    return -1;
  
  /* Free memory associated with the request records */
  reciter = reqrecs;
  while ( reciter )
    {
      recnext = reciter->next;
      
      freereqrec (reciter);
      
      reciter = recnext;
    }
  
  return 0;
}  /* End of processpod() */


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
  
  /* Get the current soft open file limit */
  if ( getrlimit (RLIMIT_NOFILE, &rlim) == -1 )
    {
      ms_log (2, "getrlimit() failed to get open file limit\n");
      return -1;
    }
  
  if ( rlim.rlim_cur < limit )
    {
      rlim.rlim_cur = (rlim_t) limit;
      
      if ( verbose > 1 )
	ms_log (1, "Setting open file limit to %d\n",
		(int) rlim.rlim_cur);
      
      if ( setrlimit (RLIMIT_NOFILE, &rlim) == -1 )
	{
	  ms_log (2, "ERROR setrlimit failed to set open file limit\n");
	  return -1;
	}
    }
  
  return (int) rlim.rlim_cur;
}  /* End of setofilelimit() */


/***************************************************************************
 * readreqfile:
 *
 * Parse the specified request file and return a chain of ReqRec structs.
 *
 * The request file (h. file) is can reference the same file multiple
 * times, but these request lines will be grouped together for
 * processing and the output request file will only have one line per
 * input file.  The data and request start/end times will be the
 * outermost times of any grouping.
 *
 * The request file is assumed to NEVER reference a file for more than
 * one channel, in otherwords there is never two or more unique
 * channels in a given file.
 *
 * Returns a pointer on success and NULL on error.
 ***************************************************************************/
static ReqRec *
readreqfile (char *requestfile)
{
  FILE *rf;
  ReqRec *rr = 0;
  ReqRec *lastrr = 0;
  ReqRec *newrr;
  ReqRec *looprr;
  char reqline[200];
  StrList *list = 0;
  StrList *listptr;
  int reqreccnt = 0;

  char tmpfilename[1024];

  if ( verbose )
    ms_log (1, "Reading request file: %s\n", requestfile);
  
  if ( ! (rf = fopen (requestfile, "rb")) )
    {
      ms_log (2, "Cannot open request file '%s': %s\n",
	      requestfile, strerror(errno));
      return NULL;
    }
  
  while ( fgets(reqline, sizeof(reqline), rf) )
    {
      if ( strparse(reqline, "\t", &list) != 10 )
	{
	  if ( verbose > 2 )
	    ms_log (1, "skipping request line: '%s'\n", reqline);
	  continue;
	}
      
      reqreccnt++;
      
      newrr = (ReqRec *) malloc (sizeof(ReqRec));
      newrr->strlist = list;
      
      listptr = list;
      
      newrr->station = strdup(listptr->element);   /* 1 */
      listptr = listptr->next;
      newrr->network = strdup(listptr->element);   /* 2 */
      listptr = listptr->next;
      newrr->channel = strdup(listptr->element);   /* 3 */
      listptr = listptr->next;
      newrr->location = strdup(listptr->element);  /* 4 */
      listptr = listptr->next;
      newrr->datastart = (time_t) (MS_HPTIME2EPOCH (ms_seedtimestr2hptime(listptr->element))); /* 5 */
      listptr = listptr->next;
      newrr->dataend = (time_t) (MS_HPTIME2EPOCH (ms_seedtimestr2hptime(listptr->element)));   /* 6 */
      listptr = listptr->next;
      newrr->filename = strdup(listptr->element);  /* 7 */
      listptr = listptr->next;
      newrr->headerdir = strdup(listptr->element); /* 8 */
      listptr = listptr->next;
      newrr->reqstart = (time_t) (MS_HPTIME2EPOCH (ms_seedtimestr2hptime(listptr->element)));  /* 9 */
      listptr = listptr->next;
      newrr->reqend = (time_t) (MS_HPTIME2EPOCH (ms_seedtimestr2hptime(listptr->element)));   /* 10 */
      listptr = listptr->next;
      
      newrr->pruned = 0;
      newrr->prev = 0;
      newrr->next = 0;
      
      /* Free the parsed string list */
      strparse (0, 0, &list);

      /* Build appropriate data file name using base data dir and record file name */
      snprintf (tmpfilename, sizeof(tmpfilename), "%s/%s/%s",
		poddatadir, newrr->station, newrr->filename);
      
      /* Check for file access */
      if ( access(tmpfilename, F_OK) )
	{
	  ms_log (2, "Cannot find file '%s', keeping a placeholder\n", tmpfilename);
	  newrr->pruned = 1;
	}
      
      /* If this is the first ReqRec */
      if ( ! rr )
	{
	  rr = newrr;
	  lastrr = newrr;
	}
      /* Other wise add to or update the ReqRec list */
      else
	{
	  /* Search the list for matching files */
	  looprr = rr;
	  while ( looprr )
	    {
	      /* If a matching file name was found update the time stamps and
	       * free this new one.  One buried assumption with this approach
	       * is that there will never be more than one channel in a given
	       * data file (unique Net, Sta, Loc, Chan & Quality) */
	      if ( ! strcmp (looprr->filename, newrr->filename) )
		{
		  if ( looprr->datastart > newrr->datastart )
		    looprr->datastart = newrr->datastart;
		  if ( looprr->dataend < newrr->dataend )
		    looprr->dataend = newrr->dataend;
		  
		  if ( looprr->reqstart > newrr->reqstart )
		    looprr->reqstart = newrr->reqstart;
		  if ( looprr->reqend < newrr->reqend )
		    looprr->reqend = newrr->reqend;
		  
		  free ( newrr );
		  break;
		}
	      
	      looprr = looprr->next;
	    }
	  
	  /* If no match add new ReqRec to the list */
	  if ( ! looprr )
	    {
	      newrr->prev = lastrr;
	      lastrr->next = newrr;
	      lastrr = newrr;
	    }
	}
    }
  
  if ( verbose )
    ms_log (1, "Read %d request records (lines)\n", reqreccnt);
  
  fclose (rf);
  
  return rr;
} /* End of readreqfile() */


/***************************************************************************
 * writereqfile:
 *
 * Write a request file for a chained ReqRec list.
 *
 * Returns 0 on success and -1 on error.
 ***************************************************************************/
static int
writereqfile (char *requestfile, ReqRec *rrlist)
{
  FILE *rf;
  ReqRec *rr = 0;
  char datastart[30];
  char dataend[30];
  char reqstart[30];
  char reqend[30];
  
  int reqreccnt = 0;

  if ( verbose )
    ms_log (1, "Writing request file: %s\n", requestfile);
  
  if ( ! (rf = fopen (requestfile, "wb")) )
    {
      ms_log (2, "Cannot open request file '%s': %s\n",
	      requestfile, strerror(errno));
      return -1;
    }
  
  rr = rrlist;
  
  while ( rr )
    {
      reqreccnt++;
      
      strftime (datastart, sizeof(datastart), "%Y,%j,%H:%M:%S", gmtime(&rr->datastart));
      strftime (dataend, sizeof(dataend), "%Y,%j,%H:%M:%S", gmtime(&rr->dataend));
      strftime (reqstart, sizeof(reqstart), "%Y,%j,%H:%M:%S", gmtime(&rr->reqstart));
      strftime (reqend, sizeof(reqend), "%Y,%j,%H:%M:%S", gmtime(&rr->reqend));
      
      fprintf (rf,"%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
	       rr->station,rr->network,rr->channel,rr->location,
	       datastart, dataend,
	       rr->filename, rr->headerdir,
	       reqstart, reqend);
      
      rr = rr->next;
    }
  
  fclose (rf);
  
  if ( verbose )
    ms_log (1, "Wrote %d request records (lines)\n", reqreccnt);

  return 0;
} /* End of writereqfile() */


/***************************************************************************
 * processtraces:
 *
 * Process all MSTrace entries in the global MSTraceGroups by first
 * pruning them and then writing out the remaining data.
 *
 * Returns 0 on success and -1 on error.
 ***************************************************************************/
static int
processtraces (MSTraceGroup *mstg, Filelink *filelist)
{
  
  /* Sort global trace group by srcname, sample rate, starttime and descending end time */
  mst_groupsort (mstg, 0);
  
  if ( verbose > 2 )
    printtracemap (mstg);
  
  /* Prune data */
  if ( prunedata )
    {
      /* Perform pre-identified pruning actions */
      if ( prunetraces (mstg) )
	return -1;
      
      /* Reconcile MSTrace times with associated record maps */
      if ( reconcile_tracetimes (mstg) )
	return -1;
      
      /* Re-sort trace group by srcname, sample rate, starttime and descending end time */
      mst_groupsort (mstg, 0);
    }
  
  /* Write all MSTrace associated records to output file(s) */
  if ( writetraces (mstg, filelist) )
    return -1;
  
  return 0;
}  /* End of processtraces() */


/***************************************************************************
 * writetraces():
 *
 * Write all MSTrace associated records to output file(s).  If an output
 * file is specified all records will be written to it, otherwise
 * records will be written to the original files and (optionally)
 * backups of the original files will remain.
 * 
 * This routine will also call trimrecord() to trim a record when data
 * suturing is requested.  Record trimming is triggered when
 * MSRecord.newstart or MSRecord.newend are set for any output
 * records.
 *
 * The quality flag is optionally set for all output records.
 *
 * Returns 0 on success and 1 on error.
 ***************************************************************************/
static int
writetraces (MSTraceGroup *mstg, Filelink *filelist)
{
  static int totalrecsout = 0;
  static int totalbytesout = 0;
  char *wb = "wb";
  char *ab = "ab";
  char *mode;
  char errflag = 0;
  
  hptime_t hpdelta;
  
  MSTrace *mst;
  RecordMap *recmap;
  Record *rec;
  Filelink *flp;
  Archive *arch;
  
  FILE *ofp = 0;
  
  if ( ! mstg )
    return 1;
  
  if ( ! mstg->traces )
    return 1;
  
  /* Open the output file if specified */
  if ( outputfile )
    {
      /* Decide if we are appending or overwriting */
      mode = ( totalbytesout ) ? ab : wb;
      
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
  
  mst = mstg->traces;
  
  /* Loop through each MSTrace in the MSTraceGroup */
  while ( mst && ! errflag )
    {
      recmap = (RecordMap *) mst->prvtptr;
      rec = recmap->first;
      
      /* Loop through each Record in the MSTrace's RecordMap.
       * After records are read from the input files, perform any
       * pre-identified pruning before writing data back out */
      while ( rec && ! errflag )
	{
	  /* Skip marked (pre-identified as non-contributing) records */
	  if ( rec->reclen == 0 )
	    {
	      rec = rec->next;
	      continue;
	    }
	  
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
	  
	  /* Trim data from the record if new start or end times are specifed */
	  if ( rec->newstart || rec->newend )
	    {
	      if ( trimrecord (rec, recordbuf) )
		{
		  rec = rec->next;
		  continue;
		}
	    }
	  
	  /* Re-stamp quality indicator if specified */
          if ( restampqind )
            {
              if ( verbose > 1 )
                ms_log (1, "Re-stamping data quality indicator to '%c'\n", restampqind);
              
              *(recordbuf + 6) = restampqind;
            }
	  
	  /* Write to a single output file if specified */
	  if ( ofp )
	    {
	      if ( fwrite (recordbuf, rec->reclen, 1, ofp) != 1 )
		{
		  ms_log (2, "Cannot write to '%s'\n", outputfile);
		  errflag = 1;
		  break;
		}
	    }
	  
	  /* Write to Archive(s) if specified */
	  if ( archiveroot )
	    {
	      MSRecord *msr = 0;
	      
	      if ( msr_unpack (recordbuf, rec->reclen, &msr, 0, verbose) != MS_NOERROR )
		{
		  ms_log (2, "Cannot unpack Mini-SEED, cannot write to archive\n");
		}
	      else
		{
		  arch = archiveroot;
		  while ( arch )
		    {
		      ds_streamproc (&arch->datastream, msr, 0, verbose-1);
		      arch = arch->next;
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
	      
	      if ( fwrite (recordbuf, rec->reclen, 1, rec->flp->outfp) != 1 )
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
	      hpdelta = ( mst->samprate ) ? (hptime_t) (HPTMODULUS / mst->samprate) : 0;
	      rec->flp->latest = rec->endtime + hpdelta;
	    }
	  
	  rec->flp->byteswritten += rec->reclen;
	  
	  totalrecsout++;
	  totalbytesout += rec->reclen;
	  
	  rec = rec->next;
	} /* Done looping through Records in the RecordMap */
      
      mst = mst->next;
    } /* Done looping through MStraces in the MSTraceGroup */
  
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
      ms_log (1, "Wrote %d bytes of %d records to output file(s)\n",
	      totalbytesout, totalrecsout);
    }
  
  return errflag;
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
 * Return 0 on success and -1 on failure.
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
  int packedsamples;
  int packedrecords;
  int retcode;
  
  if ( ! rec || ! recordbuf )
    return -1;
  
  /* Sanity check for new start/end times */
  if ( (rec->newstart && rec->newend && rec->newstart >= rec->newend) ||
       (rec->newstart && (rec->newstart < rec->starttime || rec->newstart >= rec->endtime)) ||
       (rec->newend && (rec->newend > rec->endtime || rec->newend <= rec->starttime)) )
    {
      ms_log (2, "Problem with new start/end record bound times, skipping.\n");
      ms_log (2, "  Original record from %s\n", rec->flp->infilename);
      ms_hptime2seedtimestr (rec->starttime, stime, 1);
      ms_hptime2seedtimestr (rec->endtime, etime, 1);
      ms_log (2, "       Start: %s       End: %s\n", stime, etime);
      if ( rec->newstart == 0 ) strcpy (stime, "NONE");
      else ms_hptime2seedtimestr (rec->newstart, stime, 1);
      if ( rec->newend == 0 ) strcpy (etime, "NONE");
      else ms_hptime2seedtimestr (rec->newend, etime, 1);
      ms_log (2, " Start bound: %-24s End bound: %-24s\n", stime, etime);
    }
  
  /* Unpack data record */
  if ( (retcode = msr_unpack(recordbuf, rec->reclen, &msr, 1, verbose-1)) != MS_NOERROR )
    {
      ms_log (2, "Cannot unpack Mini-SEED record: %s\n", ms_errorstr(retcode));
      return -1;
    }
  
  if ( verbose > 1 )
    {
      msr_srcname (msr, srcname, 0);
      ms_log (1, "Triming record: %s (%c)\n", srcname, msr->dataquality);
      ms_hptime2seedtimestr (rec->starttime, stime, 1);
      ms_hptime2seedtimestr (rec->endtime, etime, 1);
      ms_log (1, "       Start: %s        End: %s\n", stime, etime);
      if ( rec->newstart == 0 ) strcpy (stime, "NONE");
      else ms_hptime2seedtimestr (rec->newstart, stime, 1);
      if ( rec->newend == 0 ) strcpy (etime, "NONE");
      else ms_hptime2seedtimestr (rec->newend, etime, 1);
      ms_log (1, " Start bound: %-24s  End bound: %-24s\n", stime, etime);
    }
  
  /* Determine sample period in high precision time ticks */
  hpdelta = ( msr->samprate ) ? (hptime_t) (HPTMODULUS / msr->samprate) : 0;
  
  /* Remove samples from the beginning of the record */
  if ( rec->newstart && hpdelta )
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
	ms_log (2, "All %d samples trimmed from record ??\n", trimsamples);
      
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
  if ( rec->newend && hpdelta )
    {
      hptime_t newendtime;
      
      /* Determine new start time and the number of samples to trim */
      trimsamples = 0;
      newendtime = rec->endtime;
      
      while ( newendtime > rec->newend && trimsamples < msr->samplecnt )
	{
	  newendtime -= hpdelta;
	  trimsamples++;
	}
      
      if ( trimsamples >= msr->samplecnt )
	ms_log (2, "All %d samples trimmed from record ??\n", trimsamples);
      
      if ( verbose > 2 )
	{
	  ms_hptime2seedtimestr (newendtime, etime, 1);
	  ms_log (1, "Removing %d samples from the end, new end time: %s\n", trimsamples, etime);
	}
      
      msr->numsamples -= trimsamples;
      msr->samplecnt -= trimsamples;
      rec->endtime = newendtime;
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
 * the specified MSTraces.
 *
 * For each MSTrace determine the coverage of the RecordMap associated
 * with each overlapping, higher-priority MSTrace using getcoverage().
 * If some higher-priority overlap was determined to exist modify the
 * RecordMap of the MSTrace in question to mark the overlapping data
 * using trimtrace().
 *
 *
 * Return 0 on success and -1 on failure.
 ***************************************************************************/
static int
prunetraces (MSTraceGroup *mstg)
{
  MSTrace *mst = 0;
  MSTraceGroup *coverage = 0;
  int retval;
  
  if ( ! mstg )
    return -1;
  
  if ( ! mstg->traces )
    return -1;
  
  if ( verbose )
    ms_log (1, "Pruning MSTrace data\n");
  
  /* For each MSTrace determine the coverage of the overlapping
     Records from the other traces with a higher priority and prune
     the overlap. */
  mst = mstg->traces;
  while ( mst )
    {
      /* Determine overlapping MSTrace coverage */
      retval = getcoverage (mstg, mst, &coverage);
      
      if ( retval )
	{
	  ms_log (2, "cannot getcoverage()\n");
	  return -1;
	}
      else if ( coverage )
	{
	  if ( trimtrace (mst, coverage) < 0 )
	    {
	      ms_log (2, "cannot trimtraces()\n");
	      return -1;
	    }
	}
      
      /* Free the coverage MSTraceGroup for reuse */
      mst_freegroup (&coverage);
      
      mst = mst->next;
    }
  
  return 0;
}  /* End of prunetraces() */


/***************************************************************************
 * getcoverage():
 *
 * Search an MSTraceGroup for entries that overlap the target MSTrace
 * and, from the Record entries of the overlapping MSTraces, build a
 * coverage map in the form of a sparsely-populated MSTraceGroup.
 *
 * Only data with a higher priority than the target MSTrace will be
 * added to the overlap coverage.  Priority is determined using the
 * quality codes via qcompare() and if the qualities are equal the
 * longest time-series will be given priority.
 *
 * On success a new MSTraceGroup will be allocated and returned, it is
 * up to the caller to properly free this memory.
 *
 * To avoid numerous string operations this routine will compare
 * MSTrace srcnames by comparing 44 bytes starting from
 * MSTrace->network which should be the strings for network, station,
 * location and channel.  As long as the MSTrace struct does not
 * change this shortcut will be valid.
 *
 * When no overlap coverage is found *ppcoverage will be 0, otherwise
 * it will contain a list of MSTraces representing the overlap
 * coverage.
 *
 * Returns 0 on success and -1 on error.
 ***************************************************************************/
static int
getcoverage (MSTraceGroup *mstg, MSTrace *targetmst, MSTraceGroup **ppcoverage)
{
  MSTrace *cmst = 0;
  MSTrace *mst;
  RecordMap *recmap;
  Record *rec;
  hptime_t hpdelta, hptimetol;
  hptime_t effstarttime, effendtime;
  int priority;
  int newsegment;
  
  if ( ! mstg || ! targetmst || ! ppcoverage )
    return -1;
  
  /* Reset coverage MSTraceGroup if needed */
  if ( *ppcoverage )
    mst_freegroup (ppcoverage);
  
  /* Determine sample period in high precision time ticks */
  hpdelta = ( targetmst->samprate ) ? (hptime_t) (HPTMODULUS / targetmst->samprate) : 0;
  
  /* Determine time tolerance in high precision time ticks */
  hptimetol = ( timetol == -1 ) ? (hpdelta / 2) : (hptime_t) (HPTMODULUS * timetol);
  
  /* Loop through each MSTrace in the group searching for target coverage */
  mst = mstg->traces;
  while ( mst )
    {
      /* Skip target trace */
      if ( mst == targetmst )
	{
	  mst = mst->next;
	  continue;
	}
      
      /* Continue with next if srcname or sample rate are different */
      if ( memcmp (mst->network, targetmst->network, 44) ||
	   ! MS_ISRATETOLERABLE (mst->samprate, targetmst->samprate) )
	{
	  mst = mst->next;
	  continue;
	}
      
      /* Test for overlap with targetmst */
      if ( targetmst->endtime > mst->starttime &&
	   targetmst->starttime < mst->endtime )
	{
	  /* Determine priority:
	   *  -1 : mst > targetmst
	   *   0 : mst == targetmst
	   *   1 : mst < targetmst */
	  priority = 0;
	  
	  /* If best quality is requested compare the qualities to determine priority */
	  if ( bestquality )
	    priority = qcompare (mst->dataquality, targetmst->dataquality);
	  
	  /* If priorities are equal (qualities are equal or no checking) 
	   * give priority to the longest MSTrace */
	  if ( priority == 0 )
	    {
	      if ( (mst->endtime - mst->starttime) > (targetmst->endtime - targetmst->starttime) )
		priority = -1;
	      else
		priority = 1;
	    }
	  
	  /* If overlapping trace is a higher priority than targetmst add to coverage */
	  if ( priority == -1 )
	    {
	      /* Loop through list of associated Records, determine
		 contiguous coverage and store in an MSTraceGroup */
	      recmap = (RecordMap *) mst->prvtptr;
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
		  effstarttime = ( rec->newstart ) ? rec->newstart : rec->starttime;
		  effendtime = ( rec->newend ) ? rec->newend : rec->endtime;
		  
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
			  
			  cmst->dataquality = mst->dataquality;
			  cmst->samprate = mst->samprate;
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
      
      mst = mst->next;
    }
  
  /* Heal the coverage MSTraceGroup */
  if ( *ppcoverage )
    if ( mst_groupheal (*ppcoverage, timetol, sampratetol) < 0 )
      {
	ms_log (2, "cannot heal coverage trace group\n");
	return -1;
      }
  
  return 0;
}  /* End of getcoverage() */


/***************************************************************************
 * trimtrace():
 *
 * Adjust Record entries associated with the target MSTrace that are
 * overlapping the time represented by the coverage MSTraceGroup in
 * two different ways: 1) mark records that are completely overlapped
 * and 2) determine partial record trim boundaries (new record times)
 * if sample level pruning is requested.
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
trimtrace (MSTrace *targetmst, MSTraceGroup *coverage)
{
  RecordMap *recmap;
  Record *rec;
  MSTrace *cmst;
  hptime_t effstarttime, effendtime;
  hptime_t hpdelta, hptimetol;
  char srcname[50];
  char stime[30];
  char etime[30];
  int modcount = 0;
  
  if ( ! targetmst || ! coverage )
    return -1;
  
  /* Determine sample period in high precision time ticks */
  hpdelta = ( targetmst->samprate ) ? (hptime_t) (HPTMODULUS / targetmst->samprate) : 0;
  
  /* Determine time tolerance in high precision time ticks */
  hptimetol = ( timetol == -1 ) ? (hpdelta / 2) : (hptime_t) (HPTMODULUS * timetol);

  /* Traverse the Record chain for the target MSTrace and mark Records
   * that are completely overlapped by the HP MSTraceGroup coverage */
  recmap = (RecordMap *) targetmst->prvtptr;
  rec = recmap->first;
  while ( rec )
    {
      cmst = coverage->traces;
      while ( cmst )
	{
	  /* Determine effective record start and end times for comparison */
	  effstarttime = ( rec->newstart ) ? rec->newstart : rec->starttime;
	  effendtime = ( rec->newend ) ? rec->newend : rec->endtime;
	  
	  /* Mark Record if it is completely overlaped by the coverage */
	  if ( effstarttime >= cmst->starttime &&
	       effendtime <= cmst->endtime )
	    {
	      if ( verbose > 1 )
		{
		  mst_srcname (targetmst, srcname, 1);
		  ms_hptime2seedtimestr (rec->starttime, stime, 1);
		  ms_hptime2seedtimestr (rec->endtime, etime, 1);
		  ms_log (1, "Removing Record %s (%c) :: %s  %s\n",
			  srcname, rec->quality, stime, etime);
		}
	      
	      rec->flp->recrmcount++;
	      rec->reclen = 0;
	      modcount++;
	    }
	  
	  /* Determine the new start/end times if pruning at the sample level */
	  if ( prunedata == 's' && rec->reclen != 0 )
	    {
	      /* Record overlaps beginning of HP coverage */
	      if ( effstarttime <= cmst->starttime &&
		   effendtime >= cmst->starttime )
		{
		  /* Set Record new end time boundary including specified time tolerance */
		  rec->newend = cmst->starttime - hpdelta + hptimetol;
		  rec->flp->rectrimcount++;
		  modcount++;
		}
	      
	      /* Record overlaps end of HP coverage */
	      if ( effstarttime <= cmst->endtime &&
		   effendtime >= cmst->endtime )
		{
		  /* Set Record new start time boundary including specified time tolerance */
		  rec->newstart = cmst->endtime + hpdelta - hptimetol;
		  rec->flp->rectrimcount++;
		  modcount++;
		}
	    }
	  
	  cmst = cmst->next;
	}
      
      rec = rec->next;
    }
  
  return modcount;
}  /* End of trimtraces() */


/***************************************************************************
 * reconcile_tracetimes():
 *
 * Reconcile the start and end times of the traces in a specified
 * trace group with the list of records in an associated record map.
 * In other words, set the start and end times of each MSTrace in the
 * MSTraceGroup according to the start time of the first and end time
 * of the last contributing records in the associated record map; this
 * should be preformed after the pruning process which could mark
 * complete records as pruned (non-contributing).
 *
 * Returns 0 on success and -1 on error.
 ***************************************************************************/
static int
reconcile_tracetimes (MSTraceGroup *mstg)
{
  MSTrace *mst;
  RecordMap *recmap;
  Record *rec;
  Record *first = 0;
  Record *last = 0;
  
  if ( ! mstg )
    return -1;
  
  if ( ! mstg->traces )
    return -1;
  
  mst = mstg->traces;
  
  while ( mst )
    {
      recmap = (RecordMap *) mst->prvtptr;
      
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
      
      /* Set a new MSTrace start time */
      if ( first )
	{
	  /* Use the new boundary start time if set and sane */
	  if ( first->newstart && first->newstart > first->starttime )
	    mst->starttime = first->newstart;
	  /* Otherwise use the record start time */
	  else
	    mst->starttime = first->starttime;
	}
      
      /* Set a new MSTrace end time */
      if ( last )
	{
	  /* Use the new boundary end time if set and sane */
	  if ( last->newend && last->newend < last->endtime )
	    mst->endtime = last->newend;
	  /* Otherwise use the record end time */
	  else
	    mst->endtime = last->endtime;
	}
      
      first = 0;
      last = 0;
      mst = mst->next;
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
 * readfiles:
 *
 * Read input files specified as a Filelink list and populate an
 * MSTraceGroup and record maps for each trace.  All input files are
 * renamed with a ".orig" suffix before being read.
 *
 * Returns 0 on success and -1 otherwise.
 ***************************************************************************/
static int
readfiles (Filelink *filelist, MSTraceGroup **ppmstg)
{
  Filelink *flp;
  MSRecord *msr = 0;
  MSTrace *mst = 0;
  int retcode;
  flag whence;
  
  int totalrecs  = 0;
  int totalsamps = 0;
  int totalfiles = 0;
  
  RecordMap *recmap;
  Record *rec;
  
  RecordMap newrecmap;
  Record *newrec;
  
  off_t fpos;
  hptime_t recstarttime;
  hptime_t recendtime;
  
  char srcname[50];
  char stime[30];
  
  int infilenamelen;
  
  if ( ! filelist || ! ppmstg )
    return -1;
  
  /* (Re)Initialize MSTraceGroup */
  *ppmstg = reinitgroup (*ppmstg);
  
  if ( ! *ppmstg )
    {
      ms_log (2, "readfiles(): cannot initialize MSTraceGroup\n");
      return -1;
    }
  
  /* Read all input files and construct continuous traces, using the
   * libmseed MSTrace and MSTraceGroup functionality.  For each trace
   * maintain a list of each data record that contributed to the
   * trace, implemented as a RecordMap struct (MSTrace->prvtptr) where
   * a linked list of Record structs is maintained.  The records are
   * always listed in time order.
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
	  flp->infilename = (char *) malloc (infilenamelen);
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
	    ms_log (1, "Processing: %s (was %s)\n", flp->infilename, flp->outfilename);
	  else
	    ms_log (1, "Processing: %s\n", flp->infilename);
	}
      
      /* Loop over the input file */
      while ( (retcode = ms_readmsr (&msr, flp->infilename, reclen, &fpos, NULL, 1, 0, verbose-2))
	      == MS_NOERROR )
	{
	  recstarttime = msr_starttime (msr);
	  recendtime = msr_endtime (msr);
	  
	  /* Generate the srcname with the quality code */
	  msr_srcname (msr, srcname, 1);
	  
	  /* Generate an ASCII start time string */
	  ms_hptime2seedtimestr (recstarttime, stime, 1);
	  
	  /* Check if record matches start time criteria */
	  if ( (starttime != HPTERROR) && (recstarttime < starttime) )
	    {
	      if ( verbose >= 3 )
		ms_log (1, "Skipping (starttime) %s, %s\n", srcname, stime);
	      continue;
	    }
	  
	  /* Check if record matches end time criteria */
	  if ( (endtime != HPTERROR) && (recendtime > endtime) )
	    {
	      if ( verbose >= 3 )
		ms_log (1, "Skipping (endtime) %s, %s\n", srcname, stime);
	      continue;
	    }
	  
	  /* Check if record is matched by the match regex */
	  if ( match )
	    {
	      if ( regexec ( match, srcname, 0, 0, 0) != 0 )
		{
		  if ( verbose >= 3 )
		    ms_log (1, "Skipping (match) %s, %s\n", srcname, stime);
		  continue;
		}
	    }
	  
	  /* Check if record is rejected by the reject regex */
	  if ( reject )
	    {
	      if ( regexec ( reject, srcname, 0, 0, 0) == 0 )
		{
		  if ( verbose >= 3 )
		    ms_log (1, "Skipping (reject) %s, %s\n", srcname, stime);
		  continue;
		}
	    }
	  
	  if ( verbose > 2 )
	    msr_print (msr, verbose - 3);
	  
	  /* Add record to the MSTraceGroup */
	  if ( ! (mst = mst_addmsrtogroup (*ppmstg, msr, bestquality, timetol, sampratetol)) )
	    {
	      ms_log (2, "Cannot add record to trace group, %s, %s\n", srcname, stime);
	    }
	  
	  /* Determine where the record fit this MSTrace
	   * whence:
	   * 0 = New MSTrace
	   * 1 = End of MSTrace
	   * 2 = Beginning of MSTrace
	   */
	  whence = 0;
	  if ( mst->prvtptr )
	    {
	      if ( mst->endtime == recendtime )
		whence = 1;
	      else if ( mst->starttime == recstarttime )
		whence = 2;
	      else if ( recendtime == recstarttime )
		{
		  /* Determine best fit for records with no span (not added to the MSTrace coverage) */
		  if ( llabs (recstarttime - mst->endtime) < llabs (recstarttime - mst->starttime) )
		    whence = 1;
		  else
		    whence = 2;
		}
	      else
		{
		  ms_log (2, "Cannot determine where record fit relative to trace\n");
		  msr_print (msr, 1);
		  continue;
		}
	    }
	  
	  /* Create and populate new Record structure */
	  rec = (Record *) malloc (sizeof(Record));
	  rec->flp = flp;
	  rec->offset = fpos;
	  rec->reclen = msr->reclen;
	  rec->starttime = recstarttime;
	  rec->endtime = recendtime;
	  rec->quality = msr->dataquality;
	  rec->newstart = 0;
	  rec->newend = 0;
	  rec->prev = 0;
	  rec->next = 0;
	  
	  /* Populate a new record map */
	  newrecmap.recordcnt = 1;
	  newrecmap.first = rec;
	  newrecmap.last = rec;
	  
	  /* If pruning at the sample level trim right at the start/end times */
	  if ( prunedata == 's' )
	    {
	      /* If the Record crosses the start time */
	      if ( starttime != HPTERROR && (starttime > recstarttime) && (starttime < recendtime) )
		{
		  rec->newstart = starttime;
		}
	      
	      /* If the Record crosses the end time */
	      if ( endtime != HPTERROR && (endtime > recstarttime) && (endtime < recendtime) )
		{
		  rec->newend = endtime;
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
		  effstarttime = (rec->newstart) ? rec->newstart : rec->starttime;
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
		  if ( rec->endtime > boundary )
		    {
		      newrec = (Record *) malloc (sizeof(Record));
		      memcpy (newrec, rec, sizeof(Record));
		      
		      /* Set current Record and next Record new boundary times */
		      rec->newend = boundary;
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
	  
	  /* Add the new Record(s) to the RecordMap associated with the MSTrace */

	  /* Add new Record(s) to end of the RecordMap */
	  if ( whence == 1 )
	    {
	      recmap = (RecordMap *) mst->prvtptr;
	      
	      recmap->last->next = newrecmap.first;
	      newrecmap.first->prev = recmap->last;

	      recmap->last = newrecmap.last;
	      
	      recmap->recordcnt += newrecmap.recordcnt;
	    }
	  /* Add new Record(s) to beginning of the RecordMap */
	  else if ( whence == 2 )
	    {
	      recmap = (RecordMap *) mst->prvtptr;
	      
	      recmap->first->prev = newrecmap.last;
	      newrecmap.last->next = recmap->first;
	      
	      recmap->first = newrecmap.first;
	      
	      recmap->recordcnt += newrecmap.recordcnt;
	      
	      /* Increment reordered count */
	      flp->reordercount++;
	    }
	  /* First Record(s) for this MSTrace, allocate RecordMap */
	  else
	    {
	      if ( mst->prvtptr )
		ms_log (2, "Supposedly first record, but RecordMap not empty, report this\n");
	      
	      recmap = (RecordMap *) malloc (sizeof(RecordMap));
	      recmap->first = newrecmap.first;
	      recmap->last = newrecmap.last;
	      recmap->recordcnt = newrecmap.recordcnt;
	      
	      mst->prvtptr = recmap;
	    }
	  
	  totalrecs++;
	  totalsamps += msr->samplecnt;
	}
      
      /* Critical error if file was not read properly */
      if ( retcode != MS_ENDOFFILE )
	{
	  ms_log (2, "Cannot read %s: %s\n", flp->infilename, ms_errorstr(retcode));
	  ms_readmsr (&msr, NULL, 0, NULL, NULL, 0, 0, 0);
	  return -1;
	}
      
      /* Make sure everything is cleaned up */
      ms_readmsr (&msr, NULL, 0, NULL, NULL, 0, 0, 0);
      
      totalfiles++;
      flp = flp->next;
    } /* End of looping over file list */
  
  if ( basicsum )
    ms_log (0, "Files: %d, Records: %d, Samples: %d\n", totalfiles, totalrecs, totalsamps);
  
  return 0;
}  /* End of readfiles() */


/***************************************************************************
 * reinitgroup():
 *
 * (Re)Initialize a MSTraceGroup, freeing all associated memory.
 *
 * Return pointer to (re)initialized MSTraceGroup.
 *********************.******************************************************/
static MSTraceGroup *
reinitgroup (MSTraceGroup *mstg)
{
  MSTrace *mst;
  RecordMap *recmap;
  Record *rec, *nextrec;
  
  if ( ! mstg )
    {
      return mst_initgroup (mstg);
    }
  
  /* Free all Records from each MSTrace in group */
  mst = mstg->traces;
  
  while ( mst )
    {
      recmap = (RecordMap *) mst->prvtptr;
      rec = recmap->first;
      
      while (rec)
	{
	  nextrec = rec->next;
	  free (rec);
	  rec = nextrec;
	}
      
      free (mst->prvtptr);
      mst->prvtptr = 0;
      
      mst = mst->next;
    }
  
  mst_initgroup (mstg);
  
  return mstg;
}  /* End of reinitgroup() */


/***************************************************************************
 * printmodsummary():
 *
 * Print a summary of modifications to stdout.  If 'nomods' is true
 * include files that were not modified.
 ***************************************************************************/
static void
printmodsummary (Filelink *filelist, flag nomods)
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
 * Print record map for each MSTrace to stdout.
 ***************************************************************************/
static void
printtracemap (MSTraceGroup *mstg)
{
  MSTrace *mst = 0;
  char srcname[50];
  char stime[30];
  char etime[30];
  int tracecnt = 0;
  
  if ( ! mstg )
    return;
  
  mst = mstg->traces;
  
  /* Print out the appropriate header */
  ms_log (0, "\nMSTrace Map:\n");
  ms_log (0, "   Source              Start sample             End sample        Hz   Samples\n");
  
  while ( mst )
    {
      mst_srcname (mst, srcname, 1);
      
      /* Create formatted time strings */
      if ( ms_hptime2seedtimestr (mst->starttime, stime, 1) == NULL )
	ms_log (2, "Cannot convert trace start time for %s\n", srcname);
      
      if ( ms_hptime2seedtimestr (mst->endtime, etime, 1) == NULL )
	ms_log (2, "Cannot convert trace end time for %s\n", srcname);
      
      /* Print MSTrace header */
      ms_log (0, "%-15s %-24s %-24s %-4.4g %-d\n",
	      srcname, stime, etime, mst->samprate, mst->samplecnt);
      
      if ( ! mst->prvtptr )
	{
	  ms_log (2, "No record map associated with this MSTrace.\n");
	}
      else
	{
	  printrecordmap ((RecordMap *) mst->prvtptr, 0);
	}
      
      tracecnt++;
      mst = mst->next;
    }
  
  if ( tracecnt != mstg->numtraces )
    ms_log (2, "printtracemap(): number of traces in trace group is inconsistent\n");
  
  ms_log (0, "End of MSTrace Map: %d trace(s)\n\n", tracecnt);
  
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
      ms_log (0, "       Start: %s       End: %s\n", stime, etime);

      if ( details )
	{
	  if ( rec->newstart == 0 ) strcpy (stime, "NONE");
	  else ms_hptime2seedtimestr (rec->newstart, stime, 1);
	  if ( rec->newend == 0 ) strcpy (etime, "NONE" );
	  else ms_hptime2seedtimestr (rec->newend, etime, 1);
	  ms_log (0, " Start bound: %-24s End bound: %-24s\n", stime, etime);
	}
      
      rec = rec->next;
    }
}  /* End of printrecordmap() */


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
	  matchpattern = getoptval(argcount, argvec, optind++);
	}
      else if (strcmp (argvec[optind], "-R") == 0)
	{
	  rejectpattern = getoptval(argcount, argvec, optind++);
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
      else if (strcmp (argvec[optind], "-POD") == 0)
	{
	  if ( argvec[optind+1] )
	    if ( (optind+1) < argcount && *argvec[optind+1] != '-' )
	      podreqfile = argvec[optind+1];
	  optind++;
	  
	  if ( argvec[optind+1] )
	    if ( (optind+1) < argcount && *argvec[optind+1] != '-' )
	      poddatadir = argvec[optind+1];
	  optind++;
	  
	  if ( ! podreqfile || ! poddatadir )
	    {
	      ms_log (2, "Option -POD requires two values, try -h for usage\n");
	      exit (1);
	    }
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
	  /* Add file to global file list */
	  addfile (&gfilelist, argvec[optind], NULL);
	}
    }
  
  /* Cannot specify both input files and POD */
  if ( gfilelist && (podreqfile && poddatadir) )
    {
      ms_log (2, "Cannot specify both input files and POD structure\n");
      exit (1);
    }
  
  /* Make sure input file(s) or POD were specified */
  if ( gfilelist == 0 && ! (podreqfile && poddatadir) )
    {
      ms_log (2, "No input files were specified\n\n");
      ms_log (1, "%s version %s\n\n", PACKAGE, VERSION);
      ms_log (1, "Try %s -h for usage\n", PACKAGE);
      exit (0);
    }
  
  /* Expand match pattern from a file if prefixed by '@' */
  if ( matchpattern )
    {
      if ( *matchpattern == '@' )
	{
	  tptr = matchpattern + 1; /* Skip the @ sign */
	  matchpattern = 0;
	  
	  if ( readregexfile (tptr, &matchpattern) <= 0 )
	    {
	      ms_log (2, "Cannot read match pattern regex file\n");
	      exit (1);
	    }
	}
    }
  
  /* Expand reject pattern from a file if prefixed by '@' */
  if ( rejectpattern )
    {
      if ( *rejectpattern == '@' )
	{
	  tptr = rejectpattern + 1; /* Skip the @ sign */
	  rejectpattern = 0;
	  
	  if ( readregexfile (tptr, &rejectpattern) <= 0 )
	    {
	      ms_log (2, "Cannot read reject pattern regex file\n");
	      exit (1);
	    }
	}
    }
  
  /* Compile match and reject patterns */
  if ( matchpattern )
    {
      match = (regex_t *) malloc (sizeof(regex_t));
      
      if ( regcomp (match, matchpattern, REG_EXTENDED) != 0)
	{
	  ms_log (2, "Cannot compile match regex: '%s'\n", matchpattern);
	}
    }
  
  if ( rejectpattern )
    {
      reject = (regex_t *) malloc (sizeof(regex_t));
      
      if ( regcomp (reject, rejectpattern, REG_EXTENDED) != 0)
	{
	  ms_log (2, "Cannot compile reject regex: '%s'\n", rejectpattern);
	}
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
  
  if ( (argopt+1) < argcount && *argvec[argopt+1] != '-' )
    return argvec[argopt+1];
  
  ms_log (2, "Option %s requires a value, try -h for usage\n", argvec[argopt]);
  exit (1);
  return 0;
}  /* End of getoptval() */


/***************************************************************************
 * addfile:
 *
 * Add file to end of the specified file list.
 ***************************************************************************/
static void
addfile (Filelink **ppfilelist, char *filename, ReqRec *reqrec)
{
  Filelink *lastlp, *newlp;
  
  if ( filename == NULL )
    {
      ms_log (2, "addfile(): No file name specified\n");
      return;
    }
  
  lastlp = *ppfilelist;
  while ( lastlp != 0 )
    {
      if ( lastlp->next == 0 )
	break;
      
      lastlp = lastlp->next;
    }
  
  newlp = (Filelink *) malloc (sizeof (Filelink));
  newlp->reqrec = reqrec;
  newlp->infilename = strdup(filename);
  newlp->infp = 0;
  newlp->outfilename = 0;
  newlp->outfp = 0;
  newlp->reordercount = 0;
  newlp->recsplitcount = 0;
  newlp->recrmcount = 0;
  newlp->rectrimcount = 0;
  newlp->earliest = 0;
  newlp->latest = 0;
  newlp->byteswritten = 0;
  newlp->next = 0;
  
  if ( lastlp == 0 )
    *ppfilelist = newlp;
  else
    lastlp->next = newlp;
  
}  /* End of addfile() */


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

  newarch = (Archive *) malloc (sizeof (Archive));
  
  if ( newarch == NULL )
    {
      ms_log (2, "addarchive(): cannot allocate memory for new archive definition\n");
      return -1;
    }
  
  /* Setup new entry and add it to the front of the chain */
  pathlayout = strlen (path) + 2;
  if ( layout )
    pathlayout += strlen (layout);
  
  newarch->datastream.path = (char *) malloc (pathlayout);
  
  if ( layout )
    snprintf (newarch->datastream.path, pathlayout, "%s/%s", path, layout);
  else
    snprintf (newarch->datastream.path, pathlayout, "%s", path);
  
  newarch->datastream.grouproot = NULL;
  
  if ( newarch->datastream.path == NULL )
    {
      ms_log (2, "addarchive(): cannot allocate memory for new archive path\n");
      if ( newarch )
        free (newarch);
      return -1;
    }
  
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
  int   newpatternsize;
  
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
	  newpatternsize = strlen(*pppattern) + strlen(linepattern) + 4;
	  *pppattern = realloc (*pppattern, newpatternsize);	  
	  snprintf (*pppattern, newpatternsize, "%s|(%s)", *pppattern, linepattern);
	}
      else
	{
	  newpatternsize = strlen(linepattern) + 3;
	  *pppattern = realloc (*pppattern, newpatternsize);
	  snprintf (*pppattern, newpatternsize, "(%s)", linepattern);
	}
    }
  
  fclose (fp);
  
  return regexcnt;
}  /* End of readregexfile() */


/***************************************************************************
 * freereqrec:
 *
 * Free all memory assocated with a ReqRec.
 ***************************************************************************/
static void
freereqrec (ReqRec *reqrec)
{
  if ( ! reqrec )
    return;
  
  if ( reqrec->station )
    free (reqrec->station);
  if ( reqrec->network )
    free (reqrec->network);
  if ( reqrec->channel )
    free (reqrec->channel);
  if ( reqrec->location )
    free (reqrec->location);
  if ( reqrec->filename )
    free (reqrec->filename);
  if ( reqrec->headerdir )
    free (reqrec->headerdir);
  
  free (reqrec);
  
  return;
}  /* End of freereqrec() */


/***************************************************************************
 * freefilelist:
 *
 * Free all memory assocated with a specified file list.
 ***************************************************************************/
static void
freefilelist (Filelink **ppfilelist)
{
  Filelink *flp, *nextflp;
  
  if ( ! ppfilelist )
    return;
  
  flp = *ppfilelist;
  
  while ( flp )
    {
      nextflp = flp->next;
      
      if ( flp->infilename )
	free (flp->infilename);
      if ( flp->outfilename )
	free (flp->outfilename);
      
      free (flp);
      
      flp = nextflp;
    }
  
  *ppfilelist = 0;
  
  return;
}  /* End of freefilelist() */


/***************************************************************************
 * strparse:
 *
 * splits a 'string' on 'delim' and puts each part into a linked list
 * pointed to by 'list' (a pointer to a pointer).  The last entry has
 * it's 'next' set to 0.  All elements are NULL terminated strings.
 * If both 'string' and 'delim' are 0 then the linked list is
 * traversed and the memory used is free'd and the list pointer is set
 * to 0.
 *
 * Returns the number of elements added to the list, or 0 when freeing
 * the linked list.
 ***************************************************************************/
static int
strparse (const char *string, const char *delim, StrList **list)
{
  const char *beg;			/* beginning of element */
  const char *del;			/* delimiter */
  int stop = 0;
  int count = 0;
  int total;

  StrList *curlist = 0;
  StrList *tmplist = 0;

  if ( string && delim )
    {
      total = strlen (string);
      beg = string;
      
      while ( ! stop )
	{
	  /* Find delimiter */
	  del = strstr (beg, delim);
	  
	  /* Delimiter not found or empty */
	  if (del == NULL || strlen (delim) == 0)
	    {
	      del = string + strlen (string);
	      stop = 1;
	    }

	  tmplist = (StrList *) malloc (sizeof (StrList));
	  tmplist->next = 0;

	  tmplist->element = (char *) malloc (del - beg + 1);
	  strncpy (tmplist->element, beg, (del - beg));
	  tmplist->element[(del - beg)] = '\0';

	  /* Add this to the list */
	  if (count++ == 0)
	    {
	      curlist = tmplist;
	      *list = curlist;
	    }
	  else
	    {
	      curlist->next = tmplist;
	      curlist = curlist->next;
	    }

	  /* Update 'beg' */
	  beg = (del + strlen (delim));
	  if ((beg - string) > total)
	    break;
	}

      return count;
    }
  else
    {
      curlist = *list;
      while ( curlist )
	{
	  tmplist = curlist->next;
	  free (curlist->element);
	  free (curlist);
	  curlist = tmplist;
	}
      *list = 0;
      
      return 0;
    }
}  /* End of strparse() */


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
	   "\n"
	   " ## Data selection options ##\n"
	   " -ts time     Limit to records that start after time\n"
	   " -te time     Limit to records that end before time\n"
	   "                time format: 'YYYY[,DDD,HH,MM,SS,FFFFFF]' delimiters: [,:.]\n"
	   " -M match     Limit to records matching the specified regular expression\n"
	   " -R reject    Limit to records not matchint the specfied regular expression\n"
	   "                Regular expressions are applied to: 'NET_STA_LOC_CHAN_QUAL'\n"
	   "\n"
	   " ## Output options ##\n"
	   " -rep         Replace input files, default leaves .orig files\n"
	   " -nb          Do not keep backups of original input files if replacing them\n"
	   " -o file      Specify a single output file\n"
	   " -A format    Write all records in a custom directory/file layout (try -H)\n"
	   " -Pr          Prune data at the record level using 'best' quality priority\n"
	   " -Ps          Prune data at the sample level using 'best' quality priority\n"
	   " -S[dhm]      Split records on day, hour or minute boundaries\n"
	   " -Q DRQM      Re-stamp output data records with quality code: D, R, Q or M\n"
           "\n"
	   " ## Diagnostic output ##\n"
	   " -sum         Print a basic summary after reading all input files\n"
	   " -mod         Print summary of file modifications after processing\n"
	   "\n"
	   " ## Input data ##\n"
	   " -POD reqfile datadir\n"
	   "              Prune data from a POD structure\n"
	   " file#        Files(s) of Mini-SEED records\n"
	   "\n");
  
  if  ( level )
    {
      fprintf (stderr,
               "\n"
	       "  # Preset format layouts #\n"
	       " -CHAN dir    Write records into separate Net.Sta.Loc.Chan files\n"
	       " -QCHAN dir   Write records into separate Net.Sta.Loc.Chan.Quality files\n"
	       " -CDAY dir    Write records into separate Net.Sta.Loc.Chan-day files\n"
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
