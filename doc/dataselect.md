# <p >miniSEED data selection, sorting and pruning</p>

1. [Name](#)
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

## <a id='synopsis'>Synopsis</a>

<pre >
dataselect [options] file1 [file2 file3 ...]
</pre>

## <a id='description'>Description</a>

<p ><b>dataselect</b> selects, sorts and prunes miniSEED data.  Any output data will always be sorted by ascending time-series segments.  Various data selection operations are possible including time criteria, the removal of duplicate data and the splicing of partially overlapping data records to form continuous time-series.</p>

<p >Data pruning, removal of overlap, can be performed at either the record or sample level.  Pruning at the record level guarantees that records are never unpacked/repacked, but this could potentially leave small amounts of overlap in the data.  Pruning at the sample level will remove any overlap and splice data to within the time-series tolerance, this requires the unpacking and repacking of data records and is done in a modification-minimizing way.</p>

<p >When removing overlapping data records or samples the concept of priority is used to determine from which time-series data should be removed if overlaps are detected.  By default the priority is given to the highest publication version (or v2 quality data).  When the qualities are the same priority is given to the longer segment.</p>

<p >Multiple input files will be read in the order specified and processed all together as if all the data records were from the same file.  The program must read input data from files; input from pipes, etc. is not possible.</p>

<p >Files on the command line prefixed with a '@' character are input list files and are expected to contain a simple list of input files, see <b>INPUT LIST FILE</b> for more details.</p>

<p >Each input file may be specified with an explicit byte range to read. The program will begin reading at the specified start offset and stop reading at the specified end range.  See <b>INPUT FILE RANGE</b> for more details.</p>

## <a id='options'>Options</a>

<b>-V</b>

<p style="padding-left: 30px;">Print program version and exit.</p>

<b>-h</b>

<p style="padding-left: 30px;">Print program usage and exit.</p>

<b>-H</b>

<p style="padding-left: 30px;">Print verbose program usage including details of archive format specification and exit.</p>

<b>-v</b>

<p style="padding-left: 30px;">Be more verbose.  This flag can be used multiple times ("-v -v" or "-vv") for more verbosity.</p>

<b>-tt </b><i>secs</i>

<p style="padding-left: 30px;">Specify a time tolerance for constructing continuous trace segments. The tolerance is specified in seconds.  The default tolerance is 1/2 of the sample period.</p>

<b>-rt </b><i>diff</i>

<p style="padding-left: 30px;">Specify a sample rate tolerance for constructing continuous trace segments. The tolerance is specified as the difference between two sampling rates.  The default tolerance is tested as: (abs(1-sr1/sr2) < 0.0001).</p>

<b>-snd</b>

<p style="padding-left: 30px;">Skip non-miniSEED records.  By default the program will stop when it encounters data that cannot be identified as a miniSEED record. This option can be useful with full SEED volumes or files with bad data.</p>

<b>-E</b>\fP

<p style="padding-left: 30px;">Consider all publication versions (or v2 qualities) equal when determining priority for pruning.  By default priority is given to the data with the highest publication version.</p>

<b>-s </b><i>selectfile</i>

<p style="padding-left: 30px;">Limit processing to miniSEED records that match a selection in the specified file.  The selection file contains parameters to match the SourceID (network, station, location, channel), publication version (or v2 quality), and time range for input records. As a special case, specifying "-" will result in selection lines being read from stdin.  For more details see the <b>SELECTION FILE</b> section below.</p>

<b>-ts </b><i>time</i>

<p style="padding-left: 30px;">Limit processing to miniSEED records that start after or contain <i>time</i>.  The preferred format of the <i>time</i> argument is: 'YYYY-MM-DD[THH:MM:SS.FFFFFFFFF]'.</p>

<b>-te </b><i>time</i>

<p style="padding-left: 30px;">Limit processing to miniSEED records that end before or contain <i>time</i>.  The preferred format of the <i>time</i> argument is: 'YYYY-MM-DD[THH:MM:SS.FFFFFFFFF]'.</p>

<b>-m </b><i>match</i>

<p style="padding-left: 30px;">Limit input to records that match this globbing pattern, the <i>match</i> is tested against the full FDSN Source ID: 'FDSN:NET_STA_LOC_B_S_SS'.</p>

<b>-o </b><i>file</i>

<p style="padding-left: 30px;">Write all output data to output <i>file</i> instead of replacing the original files.  If '-' is specified as the output file all output data will be written to standard out.  By default the output file will be overwritten, changing the option to <i>+o file</i> appends to the output file.</p>

<b>-A </b><i>format</i>

<p style="padding-left: 30px;">All output records will be written to a directory/file layout defined by <i>format</i>.  All directories implied in the <i>format</i> string will be created if necessary.  The option may be used multiple times to write input records to multiple archives.  See the <b>ARCHIVE FORMAT</b> section below for more details including pre-defined archive layouts.</p>

<b>-CHAN </b><i>directory</i>

<b>-CHAN </b>\fVCHANLAYOUT\fP

<b>-QCHAN </b><i>directory</i>

<b>-CDAY </b><i>directory</i>

<b>-SDAY </b><i>directory</i>

<b>-BUD </b><i>directory</i>

<b>-SDS </b><i>directory</i>

<b>-CSS </b><i>directory</i>

<p style="padding-left: 30px;">Pre-defined output archive formats, see the <b>ARCHIVE FORMAT</b> section below for more details.</p>

<b>-Pr</b>

<p style="padding-left: 30px;">Prune, remove overlap data, at the record level.  This will result in removal of all the completely overlapped data records.  Small, partial-record overlaps might remain in the data.</p>

<b>-Ps</b>

<p style="padding-left: 30px;">Prune, remove overlap data, at the sample level.  This will result in removal of all the completely overlapped data samples.  When data records partially overlap the lowest priority record is unpacked, trimmed and repacked.  Record trimming requires a supported data encoding, if unsupported (primarily older encodings) the record will be in the output untrimmed.</p>

<b>-Pe</b>

<p style="padding-left: 30px;">Prune (trim) returned traces to user specified edges (start and end times) at the sample level. This option will not remove overlap data within specified start and end time window.  Caveats the same as for <b>-Ps</b>.</p>

<b>-Q pubversion</b>

<p style="padding-left: 30px;">Change the data publication version or quality indicator for all output records to the specified value.  If this value is one of the letters: R, D, Q or M it will be translated to the appropriate publication of 1, 2, 3, 4 respectively.  If the value is not one of these letters it must be a number between 1 and 255.  Note that miniSEED v2 data quality indicators only support values 1-4, and all higher values will result in a publication version of 4 (aka data quality 'M').</p>

<b>-out file</b>

<p style="padding-left: 30px;">Print a summary of output records to the specified file.  Any existing file will be appended to.  Specify the file as '-' to print to stdout or '--' to print to stderr.  Each summary line contains FDSN Source ID, publication version, start time, end time, byte count, and sample count for each output trace segment.</p>

<b>-outprefix prefix</b>

<p style="padding-left: 30px;">Include the specified prefix string at the beginning of each line of summary output when using the <i>-out</i> option.  This is useful to identify the summary output in a stream that is potentially mixed with other output.</p>

## <a id='the-pruning-process'>The Pruning Process</a>

<p >The pruning algorithm used is independent of the file structure and organization.  Data from all input files are parsed and a map created for every data record and their relationship in continuous time series segments.</p>

<p >Each data record time coverage in each continuous time-series is compared to the time coverage of every other continuous time-series. When overlap is detected, data is optionally removed from the lower priority time-series until the overlap is minimized or completely removed depending on the pruning option specified.</p>

## <a id='selection-file'>Selection File</a>

<p >A selection file is used to match input data records based on network, station, location and channel information.  Optionally a quality and time range may also be specified for more refined selection.  The non-time fields may use the '*' wildcard to match multiple characters and the '?' wildcard to match single characters.  Character sets may also be used, for example '[ENZ]' will match either E, N or Z. The '#' character indicates the remaining portion of the line will be ignored.</p>

<p >Example selection file entries (the first field is required)</p>
<pre >
#SourceID                  Starttime              Endtime             Pubversion
FDSN:IU_ANMO_*_B_H_?
FDSN:II                    *                      *                   3
FDSN:IU_COLA_00_L_H_[ENZ]  *                      *                   1
FDSN:IU_COLA_00_L_H_Z      2008-4-9T10:00:00Z    2008-4-9T10:30:00Z
</pre>

<p ><b>Warning:</b> with a selection file it is possible to specify multiple, arbitrary selections.  Some combinations of these selects are not possible.  See <b>CAVEATS AND LIMITATIONS</b> for more details.</p>

## <a id='input-list-file'>Input List File</a>

<p >A list file can be used to specify input files, one file per line. The initial '@' character indicating a list file is not considered part of the file name.  As an example, if the following command line option was used:</p>

<pre >
<b>@files.list</b>
</pre>

<p >The 'files.list' file might look like this:</p>

<pre >
data/day1.mseed
data/day2.mseed
data/day3.mseed
</pre>

## <a id='input-file-range'>Input File Range</a>

<p >Each input file may be specified with an associated byte range to read.  The program will begin reading at the specified start offset and finish reading when at or beyond the end offset.  The range is specified by appending an '@' character to the filename with the start and end offsets separated by a dash:</p>

<pre >
filename.mseed@[startoffset][-][endoffset]
</pre>

<p >For example: "filename.mseed@4096-8192".  Both the start and end offsets are optional.  The dash separator is optional if no end offset is specified.</p>

## <a id='archive-format'>Archive Format</a>

<p >The pre-defined archive layouts are as follows:</p>

<pre >
-CHAN dir   :: dir/%n.%s.%l.%c
-VCHAN dir  :: dir/%n.%s.%l.%c.%v
-QCHAN dir  :: dir/%n.%s.%l.%c.%q
-CDAY dir   :: dir/%n.%s.%l.%c.%Y:%j:#H:#M:#S
-SDAY dir   :: dir/%n.%s.%Y:%j
-BUD dir    :: dir/%n/%s/%s.%n.%l.%c.%Y.%j
-SDS dir    :: dir/%Y/%n/%s/%c.D/%n.%s.%l.%c.D.%Y.%j
-CSS dir    :: dir/%Y/%j/%s.%c.%Y:%j:#H:#M:#S
</pre>

<p >An archive format is expanded for each record using the following substitution flags:</p>

<pre >
  <b>n</b> : network code, white space removed
  <b>s</b> : station code, white space removed
  <b>l</b> : location code, white space removed
  <b>c</b> : channel code, white space removed
  <b>Y</b> : year, 4 digits
  <b>y</b> : year, 2 digits zero padded
  <b>j</b> : day of year, 3 digits zero padded
  <b>H</b> : hour, 2 digits zero padded
  <b>M</b> : minute, 2 digits zero padded
  <b>S</b> : second, 2 digits zero padded
  <b>F</b> : fractional seconds, 4 digits zero padded
  <b>v</b> : publication version, 1-255
  <b>q</b> : data quality if possible, otherwise pub version (D, R, Q, M, or #)
  <b>L</b> : data record length in bytes
  <b>r</b> : sample rate (Hz) as a rounded integer
  <b>R</b> : sample rate (Hz) as a float with 6 digit precision
  <b>%</b> : the percent (%) character
  <b>#</b> : the number (#) character
</pre>

<p >The flags are prefaced with either the <b>%</b> or <b>#</b> modifier. The <b>%</b> modifier indicates a defining flag while the <b>#</b> indicates a non-defining flag.  All records with the same set of defining flags will be written to the same file.  Non-defining flags will be expanded using the values in the first record for the resulting file name.</p>

<p >Time flags are based on the start time of the given record.</p>

## <a id='archive-format-examples'>Archive Format Examples</a>

<p >The format string for the predefined <i>BUD</i> layout:</p>

<p ><b>/archive/%n/%s/%s.%n.%l.%c.%Y.%j</b></p>

<p >would expand to day length files named something like:</p>

<p ><b>/archive/NL/HGN/HGN.NL..BHE.2003.055</b></p>

<p >As an example of using non-defining flags the format string for the predefined <i>CSS</i> layout:</p>

<p ><b>/data/%Y/%j/%s.%c.%Y:%j:#H:#M:#S</b></p>

<p >would expand to:</p>

<p ><b>/data/2003/055/HGN.BHE.2003:055:14:17:54</b></p>

<p >resulting in day length files because the hour, minute and second are specified with the non-defining modifier.  The hour, minute and second fields are from the first record in the file.</p>

## <a id='leap-second-list-file'>Leap Second List File</a>

<p >NOTE: A list of leap seconds is included in the program and no external list should be needed unless a leap second is added after year 2023.</p>

<p >If the environment variable LIBMSEED_LEAPSECOND_FILE is set it is expected to indicate a file containing a list of leap seconds in NTP leap second list format. Some locations where this file can be obtained are indicated in RFC 8633 section 3.7: https://www.rfc-editor.org/rfc/rfc8633.html#section-3.7</p>

<p >If present, the leap seconds listed in this file will be used to adjust the time coverage for records that contain a leap second. Also, leap second indicators in the miniSEED headers will be ignored.</p>

## <a id='error-handling-and-return-codes'>Error Handling And Return Codes</a>

<p >Any significant error message will be pre-pended with "ERROR" which can be parsed to determine run-time errors.  Additionally the program will return an exit code of 0 on successful operation and 1 when any errors were encountered.</p>

## <a id='caveats-and-limitations'>Caveats And Limitations</a>

<p >With the ability to specify multiple, arbitrary data selections it is possible to specify very complex and pathological compound selections. When pruning samples from records into order to fit the requested selections, this program is limited to trimming samples from the beginning and/or end of the record.  This means it is not possible to select two or more non-intersecting time ranges from a single record. Put another way, one cannot select select data from the beginning and end, but not the middle of a record.  The work-around for this limitation is to run the program once for each selection.</p>

## <a id='author'>Author</a>

<pre >
Chad Trabant
EarthScope Data Services
</pre>


(man page 2024/10/23)
