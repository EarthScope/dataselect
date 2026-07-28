## dataselect test suite

`run-tests.py` exercises the documented dataselect command line usage and checks
the resulting miniSEED.  Run it after building:

    make test

or directly:

    python3 test/run-tests.py

Requires Python 3, no other packages.  A single test or group can be selected by
name, for example:

    python3 test/run-tests.py Archive
    python3 test/run-tests.py Output.test_quality_indicator_v2

### Test data

No test data is stored here, all input comes from `../libmseed/test/data`.  Cases
that need data which does not exist there are created during the run: overlapping
coverage by reading the same file twice, and records with a specific year or
encoding by copying a file and patching the header (recomputing the CRC).

Output is written to `test-tmp/`, which is removed when all tests pass and left in
place when one does not.

### What is covered

Program invocation (`-V`, `-h`, `-H`), input forms (plain files, byte ranges, list
files, miniSEED 2 and 3), selection (`-ts`, `-te`, `-s`, `-m`), pruning (`-Pr`,
`-Ps`, `-Pe`, `-E`, `-F`, `-tt`, `-rt`), output (`-o`, `+o`, `-Q`), archive writing
(`-A`, every layout code, all eight presets), summary logging (`-out`,
`-outprefix`), and exit status.

Beyond option behavior the suite checks the written miniSEED: records are
byte-identical when copied, the archive and `-o` paths produce the same bytes,
every miniSEED 3 record has a valid CRC and a self consistent length, and extra
headers survive the re-packing done when records are trimmed.

Not covered: paths that require injecting a failure, such as a full filesystem, a
failed allocation, or an archive file that cannot be opened.
