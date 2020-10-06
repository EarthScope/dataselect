# dataselect - Selection and sorting for data in miniSEED format.

This program is a general use tool for extracting a subset and sorting data in
miniSEED format.

For usage infromation see the [dataselect manual](doc/dataselect.md) in the
'doc' directory.

## Downloading

The [releases](https://github.com/iris-edu/dataselect/releases) area contains release versions.

## Building and installing

Building the program requires a C compiler (C99-compatible) and the GNU make program. On a mac system, this can be accomplished with the installation of the Xcode Developer Tools.

In most environments a simple `make`, in the source directory, will build the program.
The CC and CFLAGS environment variables can be used to configure the build parameters.

After successfully compiling the program, the `dataselect` binary may be copied to
any desired location, normally in your PATH.  The man page, in the `doc` directory, may
be copied to somewhere in your MANPATH for use with the `man` program.

## Licensing 

GNU GPL version 3.  See included LICENSE file for details.

Copyright (c) 2017 Chad Trabant
