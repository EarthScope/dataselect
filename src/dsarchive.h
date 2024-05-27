
#ifndef DSARCHIVE_H
#define DSARCHIVE_H

#include <time.h>

#include <libmseed.h>

/* Define pre-formatted archive layouts */
#define CHANLAYOUT  "%n.%s.%l.%c"
#define VCHANLAYOUT "%n.%s.%l.%c.%v"
#define QCHANLAYOUT "%n.%s.%l.%c.%q"
#define CDAYLAYOUT  "%n.%s.%l.%c.%Y:%j:#H:#M:#S"
#define SDAYLAYOUT  "%n.%s.%Y:%j"
#define BUDLAYOUT   "%n/%s/%s.%n.%l.%c.%Y.%j"
#define CSSLAYOUT   "%Y/%j/%s.%c.%Y:%j:#H:#M:#S"
#define SDSLAYOUT   "%Y/%n/%s/%c.D/%n.%s.%l.%c.D.%Y.%j"

typedef struct DataStreamGroup_s
{
  char   *defkey;
  int     filed;
  time_t  modtime;
  struct  DataStreamGroup_s *next;
}
DataStreamGroup;

typedef struct DataStream_s
{
  char   *path;
  int     idletimeout;
  struct  DataStreamGroup_s *grouproot;
}
DataStream;

/* Maximum number of open files for all DataStreams */
extern int ds_maxopenfiles;

extern int ds_streamproc (DataStream *datastream, MS3Record *msr, int verbose,
                          int (expand_code) (const char *code, MS3Record *msr,
                                             char *expanded, int expandedlen));

#endif /* DSARCHIVE_H */
