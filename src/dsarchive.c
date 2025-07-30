/***************************************************************************
 * dsarchive.c
 * Routines to archive miniSEED data records.
 *
 * Written by Chad Trabant
 *   ORFEUS/EC-Project MEREDIAN
 *   IRIS Data Management Center
 *   EarthScope Data Services
 *
 * The philosophy: a "DataStream" describes an archive that miniSEED
 * records will be saved to.  Each archive can be separated into
 * "DataStreamGroup"s, each unique group will be saved into a unique
 * file.  The definition of the groups is implied by the format of the
 * archive.
 *
 * Version 2024.5.16
 ***************************************************************************/

#include <errno.h>
#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/resource.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <time.h>
#include <unistd.h>

#include <libmseed.h>

#include "dsarchive.h"

/* Maximum number of open files */
int ds_maxopenfiles = 0;
int ds_openfilecount = 0;

/* For a linked list of strings, as filled by strparse() */
typedef struct strlist_s
{
  char *element;
  struct strlist_s *next;
} strlist;

/* Functions internal to this source file */
static DataStreamGroup *ds_getstream (DataStream *datastream, const char *defkey,
                                      const char *filename);
static int ds_openfile (DataStream *datastream, const char *filename);
static int ds_closeidle (DataStream *datastream, int idletimeout);
static void ds_shutdown (DataStream *datastream);
static int strparse (const char *string, const char *delim, strlist **list);

static int dsverbose;

/***************************************************************************
 * ds_streamproc:
 *
 * Save miniSEED records in a custom directory/file structure.  The
 * appropriate directories and files are created if necessary.  If
 * files already exist they are appended to.  If 'msr' is NULL then
 * ds_shutdown() will be called to close all open files and free all
 * associated memory.
 *
 * NOTE: the expand_code() callback function is not yet implemented.
 *
 * Returns 0 on success, -1 on error.
 ***************************************************************************/
int
ds_streamproc (DataStream *datastream, MS3Record *msr, int verbose,
               int (expand_code) (const char *code, MS3Record *msr,
                                  char *expanded, int expandedlen))
{
  (void) expand_code; /* Suppress warnings while unused */
  DataStreamGroup *foundgroup = NULL;
  strlist *fnlist, *fnptr;
  char network[11] = {0};
  char station[11] = {0};
  char location[11] = {0};
  char channel[11] = {0};
  uint16_t year;
  uint16_t yday;
  uint8_t hour;
  uint8_t min;
  uint8_t sec;
  uint32_t nsec;

  char filename[400] = {0};
  char definition[400] = {0};
  char pathformat[600] = {0};
  char tstr[20] = {0};
  int fnlen = 0;

  /* Set Verbosity for ds_ functions */
  dsverbose = verbose;

  /* Special case for stream shutdown */
  if (!msr)
  {
    if (dsverbose >= 1)
      fprintf (stderr, "Closing archiving for: %s\n", datastream->path);

    ds_shutdown (datastream);
    return 0;
  }

  /* Build file path and name from datastream->path */
  filename[0] = '\0';
  definition[0] = '\0';
  snprintf (pathformat, sizeof (pathformat), "%s", datastream->path);

  if (strparse (pathformat, "/", &fnlist) < 0)
  {
    fprintf (stderr, "%s(): error parsing path format: '%s'\n",
             __func__, pathformat);
    return -1;
  }

  fnptr = fnlist;

  /* Special case of an absolute path (first entry is empty) */
  if (*fnptr->element == '\0')
  {
    if (fnptr->next != 0)
    {
      strncat (filename, "/", sizeof (filename) - 1);
      fnptr = fnptr->next;
    }
    else
    {
      fprintf (stderr, "%s(): empty path format\n", __func__);
      strparse (NULL, NULL, &fnlist);
      return -1;
    }
  }

  /* Decompose start time to date-time values */
  if (ms_nstime2time (msr->starttime, &year, &yday, &hour, &min, &sec, &nsec))
  {
    fprintf (stderr, "%s(): cannot convert start time to separate fields\n", __func__);
    strparse (NULL, NULL, &fnlist);
    return -1;
  }

  /* Decompose SourceID to network station, location and channel codes */
  if (ms_sid2nslc_n (msr->sid, network, sizeof (network), station, sizeof (station),
                     location, sizeof (location), channel, sizeof (channel)))
  {
    fprintf (stderr, "%s(): cannot convert SourceID to network, station, location and channel\n", __func__);
    strparse (NULL, NULL, &fnlist);
    return -1;
  }

  while (fnptr != 0)
  {
    uint16_t tdy;
    char quality;
    char *w, *p, def;

    p = fnptr->element;

    /* Special case of no file given */
    if (*p == '\0' && fnptr->next == 0)
    {
      fprintf (stderr, "%s(): no file name specified, only %s\n",
               __func__, filename);
      strparse (NULL, NULL, &fnlist);
      return -1;
    }

    while ((w = strpbrk (p, "%#")) != NULL)
    {
      def = (*w == '%');
      *w = '\0';
      strncat (filename, p, (sizeof (filename) - fnlen));
      fnlen = strlen (filename);

      w += 1;

      switch (*w)
      {
      case 'n':
        strncat (filename, network, (sizeof (filename) - fnlen));
        if (def)
          strncat (definition, network, (sizeof (definition) - fnlen));
        fnlen = strlen (filename);
        p = w + 1;
        break;
      case 's':
        strncat (filename, station, (sizeof (filename) - fnlen));
        if (def)
          strncat (definition, station, (sizeof (definition) - fnlen));
        fnlen = strlen (filename);
        p = w + 1;
        break;
      case 'l':
        strncat (filename, location, (sizeof (filename) - fnlen));
        if (def)
          strncat (definition, location, (sizeof (definition) - fnlen));
        fnlen = strlen (filename);
        p = w + 1;
        break;
      case 'c':
        strncat (filename, channel, (sizeof (filename) - fnlen));
        if (def)
          strncat (definition, channel, (sizeof (definition) - fnlen));
        fnlen = strlen (filename);
        p = w + 1;
        break;
      case 'Y':
        snprintf (tstr, sizeof (tstr), "%04u", year);
        strncat (filename, tstr, (sizeof (filename) - fnlen));
        if (def)
          strncat (definition, tstr, (sizeof (definition) - fnlen));
        fnlen = strlen (filename);
        p = w + 1;
        break;
      case 'y':
        tdy = year;
        while (tdy > 100)
        {
          tdy -= 100;
        }
        snprintf (tstr, sizeof (tstr), "%02d", tdy);
        strncat (filename, tstr, (sizeof (filename) - fnlen));
        if (def)
          strncat (definition, tstr, (sizeof (definition) - fnlen));
        fnlen = strlen (filename);
        p = w + 1;
        break;
      case 'j':
        snprintf (tstr, sizeof (tstr), "%03u", yday);
        strncat (filename, tstr, (sizeof (filename) - fnlen));
        if (def)
          strncat (definition, tstr, (sizeof (definition) - fnlen));
        fnlen = strlen (filename);
        p = w + 1;
        break;
      case 'H':
        snprintf (tstr, sizeof (tstr), "%02u", hour);
        strncat (filename, tstr, (sizeof (filename) - fnlen));
        if (def)
          strncat (definition, tstr, (sizeof (definition) - fnlen));
        fnlen = strlen (filename);
        p = w + 1;
        break;
      case 'M':
        snprintf (tstr, sizeof (tstr), "%02d", (int)min);
        strncat (filename, tstr, (sizeof (filename) - fnlen));
        if (def)
          strncat (definition, tstr, (sizeof (definition) - fnlen));
        fnlen = strlen (filename);
        p = w + 1;
        break;
      case 'S':
        snprintf (tstr, sizeof (tstr), "%02d", (int)sec);
        strncat (filename, tstr, (sizeof (filename) - fnlen));
        if (def)
          strncat (definition, tstr, (sizeof (definition) - fnlen));
        fnlen = strlen (filename);
        p = w + 1;
        break;
      case 'N':
        snprintf (tstr, sizeof (tstr), "%09u", (int)nsec);
        strncat (filename, tstr, (sizeof (filename) - fnlen));
        if (def)
          strncat (definition, tstr, (sizeof (definition) - fnlen));
        fnlen = strlen (filename);
        p = w + 1;
        break;
      case 'q':
        /* Map publication version into a legacy data quality indicator */
        if (msr->pubversion == 1)
          quality = 'R';
        else if (msr->pubversion == 2)
          quality = 'D';
        else if (msr->pubversion == 3)
          quality = 'Q';
        else if (msr->pubversion == 4)
          quality = 'M';
        else
          quality = '?';

        /* Use quality code if mapping possible, otherwise publication version */
        if (quality != '?')
          snprintf (tstr, sizeof (tstr), "%c", quality);
        else
          snprintf (tstr, sizeof (tstr), "%u", msr->pubversion);

        strncat (filename, tstr, (sizeof (filename) - fnlen));
        if (def)
          strncat (definition, tstr, (sizeof (definition) - fnlen));
        fnlen = strlen (filename);
        p = w + 1;
        break;
      case 'v':
        snprintf (tstr, sizeof (tstr), "%u", msr->pubversion);
        strncat (filename, tstr, (sizeof (filename) - fnlen));
        if (def)
          strncat (definition, tstr, (sizeof (definition) - fnlen));
        fnlen = strlen (filename);
        p = w + 1;
        break;
      case 'L':
        snprintf (tstr, sizeof (tstr), "%d", msr->reclen);
        strncat (filename, tstr, (sizeof (filename) - fnlen));
        if (def)
          strncat (definition, tstr, (sizeof (definition) - fnlen));
        fnlen = strlen (filename);
        p = w + 1;
        break;
      case 'r':
        snprintf (tstr, sizeof (tstr), "%ld", (long int)(msr->samprate + 0.5));
        strncat (filename, tstr, (sizeof (filename) - fnlen));
        if (def)
          strncat (definition, tstr, (sizeof (definition) - fnlen));
        fnlen = strlen (filename);
        p = w + 1;
        break;
      case 'R':
        snprintf (tstr, sizeof (tstr), "%.6f", msr->samprate);
        strncat (filename, tstr, (sizeof (filename) - fnlen));
        if (def)
          strncat (definition, tstr, (sizeof (definition) - fnlen));
        fnlen = strlen (filename);
        p = w + 1;
        break;
      case '%':
        strncat (filename, "%", (sizeof (filename) - fnlen));
        fnlen = strlen (filename);
        p = w + 1;
        break;
      case '#':
        strncat (filename, "#", (sizeof (filename) - fnlen));
        fnlen = strlen (filename);
        p = w + 1;
        break;
      default:
        fprintf (stderr, "Unknown layout format code: '%c'\n", *w);
        p = w;
        break;
      }
    }

    strncat (filename, p, (sizeof (filename) - fnlen));
    fnlen = strlen (filename);

    /* If not the last entry then it should be a directory */
    if (fnptr->next != 0)
    {
      if (access (filename, F_OK))
      {
        if (errno == ENOENT)
        {
          if (dsverbose >= 1)
            fprintf (stderr, "Creating directory: %s\n", filename);

          if (mkdir (filename, S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH))
          {
            fprintf (stderr, "%s: mkdir(%s) %s\n", __func__, filename, strerror (errno));
            strparse (NULL, NULL, &fnlist);
            return -1;
          }
        }
        else
        {
          fprintf (stderr, "%s: access denied, %s\n", filename, strerror (errno));
          strparse (NULL, NULL, &fnlist);
          return -1;
        }
      }

      strncat (filename, "/", (sizeof (filename) - fnlen));
      fnlen++;
    }

    fnptr = fnptr->next;
  }

  strparse (NULL, NULL, &fnlist);

  /* Make sure the filename and definition are NULL terminated */
  *(filename + sizeof (filename) - 1) = '\0';
  *(definition + sizeof (definition) - 1) = '\0';

  /* Check for previously used stream entry, otherwise create it */
  foundgroup = ds_getstream (datastream, definition, filename);

  if (foundgroup != NULL)
  {
    if (dsverbose >= 3)
      fprintf (stderr, "Writing data record to data stream file %s\n", filename);

    if (!write (foundgroup->filed, msr->record, msr->reclen))
    {
      fprintf (stderr, "%s: failed to write data record\n", __func__);
      return -1;
    }
    else
    {
      foundgroup->modtime = time (NULL);
    }

    return 0;
  }

  return -1;
} /* End of ds_streamproc() */

/***************************************************************************
 * ds_getstream:
 *
 * Find the DataStreamGroup entry that matches the definition key, if
 * no matching entries are found allocate a new entry and open the
 * given file.
 *
 * Resource maintenance is performed here: the modification time of
 * each stream, modtime, is compared to the current time.  If the
 * stream entry has been idle for 'DataStream.idletimeout' seconds
 * then the stream will be closed (file closed and memory freed).
 *
 * Returns a pointer to a DataStreamGroup on success or NULL on error.
 ***************************************************************************/
static DataStreamGroup *
ds_getstream (DataStream *datastream, const char *defkey, const char *filename)
{
  DataStreamGroup *foundgroup = NULL;
  DataStreamGroup *searchgroup = NULL;
  DataStreamGroup *prevgroup = NULL;
  time_t curtime;

  if (!datastream)
    return NULL;

  searchgroup = datastream->grouproot;
  curtime = time (NULL);

  /* Traverse the stream chain looking for matching streams */
  while (searchgroup != NULL)
  {
    DataStreamGroup *nextgroup = (DataStreamGroup *)searchgroup->next;

    if (!strcmp (searchgroup->defkey, defkey))
    {
      if (dsverbose >= 3)
        fprintf (stderr, "Found data stream entry for key %s\n", defkey);

      foundgroup = searchgroup;

      /* Keep ds_closeidle from closing this stream */
      if (foundgroup->modtime > 0)
      {
        foundgroup->modtime *= -1;
      }

      break;
    }

    prevgroup = searchgroup;
    searchgroup = nextgroup;
  }

  /* If not found, create a stream entry */
  if (foundgroup == NULL)
  {
    if (dsverbose >= 2)
      fprintf (stderr, "Creating data stream entry for key %s\n", defkey);

    if (!(foundgroup = (DataStreamGroup *)malloc (sizeof (DataStreamGroup))))
    {
      fprintf (stderr, "%s(): ERROR, Cannot allocate memory for DataStreamGroup\n",
               __func__);
      return NULL;
    }

    foundgroup->defkey = strdup (defkey);
    foundgroup->filed = 0;
    foundgroup->modtime = -curtime;
    foundgroup->next = NULL;

    /* Set the stream root if this is the first entry */
    if (datastream->grouproot == NULL)
    {
      datastream->grouproot = foundgroup;
    }
    /* Otherwise add to the end of the chain */
    else if (prevgroup != NULL)
    {
      prevgroup->next = foundgroup;
    }
    else
    {
      fprintf (stderr, "%s(): ERROR, stream chain is broken!\n", __func__);
      return NULL;
    }
  }

  /* Close idle stream files */
  ds_closeidle (datastream, datastream->idletimeout);

  /* If no file is open, well, open it */
  if (foundgroup->filed == 0)
  {
    off_t filepos;

    if (dsverbose >= 1)
      fprintf (stderr, "Opening data stream file %s\n", filename);

    if ((foundgroup->filed = ds_openfile (datastream, filename)) == -1)
    {
      fprintf (stderr, "%s(): ERROR, cannot open data stream file, %s\n",
               __func__, strerror (errno));
      return NULL;
    }

    if ((filepos = (off_t)lseek (foundgroup->filed, (off_t)0, SEEK_END)) < 0)
    {
      fprintf (stderr, "%s(): ERROR, cannot seek in data stream file, %s\n",
               __func__, strerror (errno));
      return NULL;
    }
  }

  return foundgroup;
} /* End of ds_getstream() */

/***************************************************************************
 * ds_openfile:
 *
 * Open a specified file, if the open file limit has been reached try
 * once to increase the limit, if that fails or has already been done
 * start closing idle files with decreasing idle timeouts until a file
 * can be opened.
 *
 * Return the result of open(2), normally this a the file descriptor
 * on success and -1 on error.
 ***************************************************************************/
static int
ds_openfile (DataStream *datastream, const char *filename)
{
  static char rlimit = 0;
  struct rlimit rlim;
  int idletimeout = datastream->idletimeout;
  int oret = 0;
  int flags = (O_RDWR | O_CREAT | O_APPEND);
  mode_t mode = (S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH); /* Mode 0644 */

  /* Lookup process open file limit and change ds_maxopenfiles if needed */
  if (!rlimit)
  {
    rlimit = 1;

    if (getrlimit (RLIMIT_NOFILE, &rlim) == -1)
    {
      fprintf (stderr, "%s(): getrlimit() failed to get open file limit\n", __func__);
    }
    else
    {
      /* Increase process open file limit to ds_maxopenfiles or hard limit */
      if (ds_maxopenfiles && (rlim_t)ds_maxopenfiles > rlim.rlim_cur)
      {
        if ((rlim_t)ds_maxopenfiles > rlim.rlim_max)
          rlim.rlim_cur = rlim.rlim_max;
        else
          rlim.rlim_cur = (rlim_t)ds_maxopenfiles;

        if (dsverbose >= 2)
          fprintf (stderr, "Setting open file limit to %lld\n", (long long int)rlim.rlim_cur);

        if (setrlimit (RLIMIT_NOFILE, &rlim) == -1)
        {
          fprintf (stderr, "%s(): setrlimit() failed to set open file limit\n", __func__);
        }

        ds_maxopenfiles = rlim.rlim_cur;
      }
      /* Set max to current soft limit if not already specified */
      else if (!ds_maxopenfiles)
      {
        ds_maxopenfiles = rlim.rlim_cur;
      }
    }
  }

  /* Close open files from the DataStream if already at the limit of (ds_maxopenfiles - 10) */
  if ((ds_openfilecount + 10) > ds_maxopenfiles)
  {
    if (dsverbose >= 1)
      fprintf (stderr, "Maximum open archive files reached (%d), closing idle stream files\n",
               (ds_maxopenfiles - 10));

    /* Close idle streams until we have free descriptors */
    while (ds_closeidle (datastream, idletimeout) == 0 && idletimeout >= 0)
    {
      idletimeout = (idletimeout / 2) - 1;
    }
  }

  /* Open file */
  if ((oret = open (filename, flags, mode)) != -1)
  {
    ds_openfilecount++;
  }

  return oret;
} /* End of ds_openfile() */

/***************************************************************************
 * ds_closeidle:
 *
 * Close all stream files that have not been active for the specified
 * idletimeout.
 *
 * Return the number of files closed.
 ***************************************************************************/
static int
ds_closeidle (DataStream *datastream, int idletimeout)
{
  int count = 0;
  DataStreamGroup *searchgroup = NULL;
  DataStreamGroup *prevgroup = NULL;
  DataStreamGroup *nextgroup = NULL;
  time_t curtime;

  searchgroup = datastream->grouproot;
  curtime = time (NULL);

  /* Traverse the stream chain */
  while (searchgroup != NULL)
  {
    nextgroup = searchgroup->next;

    if (searchgroup->modtime > 0 && (curtime - searchgroup->modtime) > idletimeout)
    {
      if (dsverbose >= 2)
        fprintf (stderr, "Closing idle stream with key %s\n", searchgroup->defkey);

      /* Re-link the stream chain */
      if (prevgroup != NULL)
      {
        if (searchgroup->next != NULL)
          prevgroup->next = searchgroup->next;
        else
          prevgroup->next = NULL;
      }
      else
      {
        if (searchgroup->next != NULL)
          datastream->grouproot = searchgroup->next;
        else
          datastream->grouproot = NULL;
      }

      /* Close the associated file */
      if (close (searchgroup->filed))
        fprintf (stderr, "%s(), closing data stream file, %s\n",
                 __func__, strerror (errno));
      else
        count++;

      free (searchgroup->defkey);
      free (searchgroup);
    }
    else
    {
      prevgroup = searchgroup;
    }

    /* Avoid closing all files if we are under the limit and idle timeout is less
     * than possible  due to progresive decreases. */
    if (idletimeout < 0 &&
        (ds_openfilecount - 10 - count) < ds_maxopenfiles)
      break;

    searchgroup = nextgroup;
  }

  ds_openfilecount -= count;

  return count;
} /* End of ds_closeidle() */

/***************************************************************************
 * ds_shutdown:
 *
 * Close all stream files and release all of the DataStreamGroup memory
 * structures.
 ***************************************************************************/
static void
ds_shutdown (DataStream *datastream)
{
  DataStreamGroup *curgroup = NULL;
  DataStreamGroup *prevgroup = NULL;

  curgroup = datastream->grouproot;

  while (curgroup != NULL)
  {
    prevgroup = curgroup;
    curgroup = curgroup->next;

    if (dsverbose >= 2)
      fprintf (stderr, "Shutting down stream with key: %s\n", prevgroup->defkey);

    if (close (prevgroup->filed))
      fprintf (stderr, "%s(), closing data stream file, %s\n",
               __func__, strerror (errno));

    free (prevgroup->defkey);
    free (prevgroup);
  }
} /* End of ds_shutdown() */

/***************************************************************************
 * strparse:
 *
 * splits a 'string' on 'delim' and puts each part into a linked list
 * pointed to by 'list' (a pointer to a pointer).  The last entry has
 * it's 'next' set to 0.  All elements are NULL terminated strings.
 * If both 'string' and 'delim' are NULL then the linked list is
 * traversed and the memory used is free'd and the list pointer is
 * set to NULL.
 *
 * Returns the number of elements added to the list, 0 when freeing
 * the linked list and -1 on error.
 ***************************************************************************/
static int
strparse (const char *string, const char *delim, strlist **list)
{
  const char *beg; /* beginning of element */
  const char *del; /* delimiter */
  int stop = 0;
  int count = 0;
  int total;

  strlist *curlist = 0;
  strlist *tmplist = 0;

  if (string != NULL && delim != NULL)
  {
    total = strlen (string);
    beg = string;

    while (!stop)
    {
      /* Find delimiter */
      del = strstr (beg, delim);

      /* Delimiter not found or empty */
      if (del == NULL || strlen (delim) == 0)
      {
        del = string + strlen (string);
        stop = 1;
      }

      if (!(tmplist = (strlist *)malloc (sizeof (strlist))))
      {
        fprintf (stderr, "%s(): ERROR, Cannot allocate memory for string parsing\n", __func__);
        return -1;
      }

      tmplist->next = 0;

      if (!(tmplist->element = (char *)malloc (del - beg + 1)))
      {
        fprintf (stderr, "%s(): ERROR, Cannot allocate memory for string parsing\n", __func__);
        return -1;
      }

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
    while (curlist != NULL)
    {
      tmplist = curlist->next;
      free (curlist->element);
      free (curlist);
      curlist = tmplist;
    }
    *list = NULL;

    return 0;
  }
} /* End of strparse() */
