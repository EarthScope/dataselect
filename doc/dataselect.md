# dataselect — miniSEED data selection, sorting and pruning

1. [Synopsis](#synopsis)
1. [Description](#description)
1. [Options](#options)
1. [The Pruning Process](#the-pruning-process)
1. [Selection File](#selection-file)
1. [Input List File](#input-list-file)
1. [Input File Range](#input-file-range)
1. [Archive Format](#archive-format)
1. [Archive Format Examples](#archive-format-examples)
1. [Leap Second List File](#leap-second-list-file)
1. [Error Handling And Return Codes](#error-handling-and-return-codes)
1. [Caveats And Limitations](#caveats-and-limitations)
1. [Author](#author)

## <a id="synopsis">Synopsis</a>

```
dataselect [options] file1 [file2 file3 ...]
```

## <a id="description">Description</a>

<b>dataselect</b> selects, sorts and prunes miniSEED data.  Any output data will always be sorted by ascending time-series segments.  Various data selection operations are possible including time criteria, the removal of duplicate data and the splicing of partially overlapping data records to form continuous time-series.

Data pruning, removal of overlap, can be performed at either the record or sample level.  Pruning at the record level guarantees that records are never unpacked/repacked, but this could potentially leave small amounts of overlap in the data.  Pruning at the sample level will remove any overlap and splice data to within the time-series tolerance, this requires the unpacking and repacking of data records and is done in a modification-minimizing way.

When removing overlapping data records or samples the concept of priority is used to determine from which time-series data should be removed if overlaps are detected.  By default the priority is given to the highest publication version (or v2 quality data).  When the qualities are the same priority is given to the longer segment.

An output destination must be specified, either a single file with <b>-o</b> or an archive layout with <b>-A</b> (or a pre-defined layout option).  Both may be used together, and <b>-A</b> may be used more than once, to write the same output to multiple destinations.  The input files are never modified.

Multiple input files will be read in the order specified and processed all together as if all the data records were from the same file.  The program must read input data from regular, seekable files; input from pipes, standard input and URLs is not possible.

Records are checked against their CRC when one is present (miniSEED v3), and a mismatch is a fatal error unless <b>-snd</b> is specified.

Files on the command line prefixed with a '@' character are input list files and are expected to contain a simple list of input files, see <b>INPUT LIST FILE</b> for more details.

Each input file may be specified with an explicit byte range to read. The program will begin reading at the specified start offset and stop reading at the specified end range.  See <b>INPUT FILE RANGE</b> for more details.

## <a id="options">Options</a>

- <b>-V</b>
  Print program version and exit.

- <b>-h</b>
  Print program usage and exit.

- <b>-H</b>
  Print verbose program usage including details of archive format specification and exit.

- <b>-v</b>
  Be more verbose.  This flag can be used multiple times ("-v -v" or "-vv") for more verbosity.

- -tt <i>secs</i>
  Specify a time tolerance for constructing continuous trace segments. The tolerance is specified in seconds.  The default tolerance is 1/2 of the sample period.

- -rt <i>diff</i>
  Specify a sample rate tolerance for constructing continuous trace segments. The tolerance is specified as the difference between two sampling rates.  The default tolerance is tested as: (abs(1-sr1/sr2) &lt; 0.0001).

- <b>-snd</b>
  Skip non-miniSEED records.  By default the program will stop when it encounters data that cannot be identified as a miniSEED record. This option can be useful with full SEED volumes or files with bad data.  Records that fail CRC validation are also skipped.

- <b>-E</b>
  Consider all publication versions (or v2 qualities) equal when determining priority for pruning.  By default priority is given to the data with the highest publication version.

- <b>-F</b>
  Use input file order for 'best' prioritization, lowest to highest, for pruning regardless of publication version (or v2 quality).  Data in input files specified later will be considered higher priority than data in input files specified earlier.

- -s <i>selectfile</i>
  Limit processing to miniSEED records that match a selection in the specified file.  The selection file contains parameters to match the SourceID (network, station, location, channel), publication version (or v2 quality), and time range for input records. As a special case, specifying "-" will result in selection lines being read from stdin.  For more details see the \fBSELECTION FILE\fP section below.

- -ts <i>time</i>
  Limit processing to miniSEED records that start after or contain <i>time</i>.  The preferred format of the <i>time</i> argument is: 'YYYY-MM-DD[THH:MM:SS.FFFFFFFFF]'.

- -te <i>time</i>
  Limit processing to miniSEED records that end before or contain <i>time</i>.  The preferred format of the <i>time</i> argument is: 'YYYY-MM-DD[THH:MM:SS.FFFFFFFFF]'.

- -m <i>match</i>
  Limit input to records that match this globbing pattern, the <i>match</i> is tested against the full FDSN Source ID \&amp;'FDSN:NET_STA_LOC_B_S_SS'.  The pattern is applied as a logical "contains", i.e. it does not need to match the entire Source ID.

  The <b>-m</b>, <b>-ts</b> and <b>-te</b> options are combined into a single selection.  When they are used together with <b>-s</b> this selection is <i>additional</i> to those in the selection file, records matching any selection are processed.

- -r <i>reject</i>
  Limit input by excluding records that match this globbing pattern, the inverse of <b>-m</b>.  The <i>reject</i> is tested against the full FDSN Source ID \&amp;'FDSN:NET_STA_LOC_B_S_SS' as a logical "contains".  This option may be used multiple times; a Source ID matching any <b>-r</b> pattern is rejected.  Rejection is applied after selection, taking precedence over <b>-s</b>, <b>-m</b>, <b>-ts</b> and <b>-te</b>.

- -o <i>file</i>
  Write all output data to output <i>file</i>.  If '-' is specified as the output file all output data will be written to standard out.  By default the output file will be overwritten, changing the option to <i>+o file</i> appends to the output file.  This option may be combined with <b>-A</b> to write the same output to more than one destination.

- -A <i>format</i>
  All output records will be written to a directory/file layout defined by <i>format</i>.  All directories implied in the <i>format</i> string will be created if necessary.  The option may be used multiple times to write input records to multiple archives.  See the \fBARCHIVE FORMAT\fP section below for more details including pre-defined archive layouts.

- -CHAN <i>directory</i>

- -VCHAN <i>directory</i>

- -QCHAN <i>directory</i>

- -CDAY <i>directory</i>

- -SDAY <i>directory</i>

- -BUD <i>directory</i>

- -SDS <i>directory</i>

- -CSS <i>directory</i>
  Pre-defined output archive formats, see the <b>ARCHIVE FORMAT</b> section below for more details.

- <b>-Pr</b>
  Prune, remove overlap data, at the record level.  This will result in removal of all the completely overlapped data records.  Small, partial-record overlaps might remain in the data.

- <b>-Ps</b>
  Prune, remove overlap data, at the sample level.  This will result in removal of all the completely overlapped data samples.  When data records partially overlap the lowest priority record is unpacked, trimmed and repacked.  Record trimming requires a supported data encoding, if unsupported (primarily older encodings) the record will be in the output untrimmed.

- <b>-Pe</b>
  Prune (trim) returned traces to user specified edges (start and end times) at the sample level. This option will not remove overlap data within specified start and end time window.  Caveats the same as for <b>-Ps</b>.  This option has no effect unless time limits are specified with <b>-ts</b>, <b>-te</b> or a selection file.

- <b>-Q pubversion</b>
  Change the data publication version or quality indicator for all output records to the specified value.  If this value is one of the letters: R, D, Q or M it will be translated to the appropriate publication of 1, 2, 3, 4 respectively.  If the value is not one of these letters it must be a number between 1 and 255.  Note that miniSEED v2 data quality indicators only support values 1-4, and all higher values will result in a publication version of 4 (aka data quality 'M').

- <b>-out file</b>
  Print a summary of output records to the specified file.  Any existing file will be appended to.  Specify the file as '-' to print to stdout or '--' to print to stderr.  One line is printed per output trace segment, containing the following '|' separated fields:

  ```
  SourceID|PubVersion|Starttime|Endtime|Bytes|Samples
  ```

- <b>-outprefix prefix</b>
  Include the specified prefix string at the beginning of each line of summary output when using the <i>-out</i> option.  This is useful to identify the summary output in a stream that is potentially mixed with other output.

## <a id="the-pruning-process">The Pruning Process</a>

The pruning algorithm used is independent of the file structure and organization.  Data from all input files are parsed and a map created for every data record and their relationship in continuous time series segments.

Each data record time coverage in each continuous time-series is compared to the time coverage of every other continuous time-series. When overlap is detected, data is optionally removed from the lower priority time-series until the overlap is minimized or completely removed depending on the pruning option specified.

## <a id="selection-file">Selection File</a>

A selection file is used to match input data records based on their FDSN Source ID, which contains the network, station, location and channel information.  Optionally a publication version (or v2 quality) and time range may also be specified for more refined selection.  The non-time fields may use the '*' wildcard to match multiple characters and the '?' wildcard to match single characters.  Character sets may also be used, for example '[ENZ]' will match either E, N or Z. Empty lines and lines with a '#' as their first non-whitespace character are ignored, trailing comments are not recognized.  A record is selected if it matches any entry in the file.

The Source ID pattern is matched against the entire Source ID, which includes the leading 'FDSN:' namespace, so wildcards are needed to match a prefix such as a network code.  The time fields, when present, must be times; there is no wildcard placeholder for them, and the publication version cannot be specified without them.

Example selection file entries (only the Source ID is required)

```
#SourceID                  Starttime             Endtime               Pubversion
FDSN:IU_ANMO_*_B_H_?
FDSN:IU_COLA_00_L_H_[ENZ]  2008-04-09T10:00:00Z
FDSN:IU_COLA_00_L_H_Z      2008-04-09T10:00:00Z  2008-04-09T10:30:00Z
FDSN:II_*                  2008-04-09T00:00:00Z  2008-04-10T00:00:00Z  3
```

For compatibility, entries may alternatively be specified as separate network, station, location and channel fields, where all four fields are required:

```
#Net  Sta   Loc  Chan  Pubversion  Starttime             Endtime
IU    ANMO  *    BH?
IU    COLA  00   LHZ   1           2008-04-09T10:00:00Z  2008-04-09T10:30:00Z
```

Entries that match neither form are skipped with a warning.

<b>Warning:</b> with a selection file it is possible to specify multiple, arbitrary selections.  Some combinations of these selects are not possible.  See <b>CAVEATS AND LIMITATIONS</b> for more details.

## <a id="input-list-file">Input List File</a>

A list file can be used to specify input files, one file per line. The initial '@' character indicating a list file is not considered part of the file name.  As an example, if the following command line option was used:

```
@files.list
```

The 'files.list' file might look like this:

```
data/day1.mseed
data/day2.mseed
data/day3.mseed
```

Empty lines and lines beginning with a '#' character are ignored.  Each entry may include a byte range as described in <b>INPUT FILE RANGE</b>. List files cannot be nested.

## <a id="input-file-range">Input File Range</a>

Each input file may be specified with an associated byte range to read.  The program will begin reading at the specified start offset and finish reading when at or beyond the end offset.  The range is specified by appending an '@' character to the filename with the start and end offsets separated by a dash:

```
filename.mseed@[startoffset][-][endoffset]
```

For example: "filename.mseed@4096-8192".  Both the start and end offsets are optional.  The dash separator is optional if no end offset is specified.  A ':' is also accepted as the separator for compatibility with previous versions.

The start offset must be the beginning of a data record.

## <a id="archive-format">Archive Format</a>

The pre-defined archive layouts are as follows:

```
-CHAN dir   :: dir/%n.%s.%l.%c
-VCHAN dir  :: dir/%n.%s.%l.%c.%v
-QCHAN dir  :: dir/%n.%s.%l.%c.%q
-CDAY dir   :: dir/%n.%s.%l.%c.%Y:%j:#H:#M:#S
-SDAY dir   :: dir/%n.%s.%Y:%j
-BUD dir    :: dir/%n/%s/%s.%n.%l.%c.%Y.%j
-SDS dir    :: dir/%Y/%n/%s/%c.D/%n.%s.%l.%c.D.%Y.%j
-CSS dir    :: dir/%Y/%j/%s.%c.%Y:%j:#H:#M:#S
```

An archive format is expanded for each record using the following substitution flags:

```
  n : network code, white space removed
  s : station code, white space removed
  l : location code, white space removed
  c : channel code, white space removed
  Y : year, 4 digits
  y : year, 2 digits zero padded
  j : day of year, 3 digits zero padded
  H : hour, 2 digits zero padded
  M : minute, 2 digits zero padded
  S : second, 2 digits zero padded
  N : nanoseconds, 9 digits zero padded
  v : publication version, 1-255
  q : data quality if possible, otherwise pub version (D, R, Q, M, or #)
  L : data record length in bytes
  r : sample rate (Hz) as a rounded integer
  R : sample rate (Hz) as a float with 6 digit precision
  % : the percent (%) character
  # : the number (#) character
```

The flags are prefaced with either the <b>%</b> or <b>#</b> modifier. The <b>%</b> modifier indicates a defining flag while the <b>#</b> indicates a non-defining flag.  All records with the same set of defining flags will be written to the same file.  Non-defining flags will be expanded using the values in the first record for the resulting file name.

Time flags are based on the start time of the given record.

## <a id="archive-format-examples">Archive Format Examples</a>

The format string for the predefined <i>BUD</i> layout:

<b>/archive/%n/%s/%s.%n.%l.%c.%Y.%j</b>

would expand to day length files named something like:

<b>/archive/NL/HGN/HGN.NL..BHE.2003.055</b>

As an example of using non-defining flags the format string for the predefined <i>CSS</i> layout:

<b>/data/%Y/%j/%s.%c.%Y:%j:#H:#M:#S</b>

would expand to:

<b>/data/2003/055/HGN.BHE.2003:055:14:17:54</b>

resulting in day length files because the hour, minute and second are specified with the non-defining modifier.  The hour, minute and second fields are from the first record in the file.

## <a id="leap-second-list-file">Leap Second List File</a>

NOTE: A list of leap seconds is included in the program and no external list should be needed unless a leap second is added after year 2023.

If the environment variable LIBMSEED_LEAPSECOND_FILE is set it is expected to indicate a file containing a list of leap seconds in NTP leap second list format. Some locations where this file can be obtained are indicated in RFC 8633 section 3.7: https://www.rfc-editor.org/rfc/rfc8633.html#section-3.7

If present, the leap seconds listed in this file will be used to adjust the time coverage for records that contain a leap second. Also, leap second indicators in the miniSEED headers will be ignored.

## <a id="error-handling-and-return-codes">Error Handling And Return Codes</a>

Any significant error message will be pre-pended with "ERROR" which can be parsed to determine run-time errors.  Additionally the program will return an exit code of 0 on successful operation and 1 when any errors were encountered.

## <a id="caveats-and-limitations">Caveats And Limitations</a>

With the ability to specify multiple, arbitrary data selections it is possible to specify very complex and pathological compound selections. When pruning samples from records in order to fit the requested selections, this program is limited to trimming samples from the beginning and/or end of the record.  This means it is not possible to select two or more non-intersecting time ranges from a single record. Put another way, one cannot select data from the beginning and end, but not the middle of a record.  The work-around for this limitation is to run the program once for each selection.

## <a id="author">Author</a>

```
Chad Trabant
EarthScope Data Services
```

---

*Generated from man page dated 2026/07/29.*
