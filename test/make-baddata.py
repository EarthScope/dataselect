#!/usr/bin/env python3
"""Generate miniSEED 2 test data with controlled record overlaps.

Used to benchmark dataselect pruning, which degrades on data where records
partially overlap each other.  Such records cannot be merged into a single
trace segment, so the segment count approaches the record count.

A record from ../libmseed/test/data is replicated with patched source
identifier codes and start times; the encoded samples are left untouched.

    test/make-baddata.py scattered -r 20000 -c 4 -o bad.mseed

Shapes:
    clean      contiguous records, healed into one segment per channel
    scattered  overlapping pairs separated by gaps, spread over time
    dense      a rolling pile where each record overlaps many neighbors
    layered    the same span delivered repeatedly at slightly offset record
               boundaries, giving a few record-rich overlapping segments

The overlapping shapes place record starts so that no record is contiguous
with any other within the time tolerance, which is what stops libmseed from
healing them and leaves one trace segment per record.
"""

import argparse
import os
import struct
import sys

TESTDIR = os.path.dirname(os.path.abspath(__file__))
TEMPLATE = os.path.join(
    TESTDIR, os.pardir, "libmseed", "test", "data", "reference-testdata-steim2.mseed2"
)

RECLEN = 512

# miniSEED 2 fixed section field offsets, from libmseed/mseedformat.h
MS2_SEQNUM = 0
MS2_STATION = 8
MS2_LOCATION = 13
MS2_CHANNEL = 15
MS2_NETWORK = 18
MS2_YEAR = 20
MS2_NUMSAMPLES = 30

# Ticks of the miniSEED 2 start time fraction field, 0.0001 seconds
TICKSPERSEC = 10000
SECPERDAY = 86400


def daysinyear(year):
    leap = (year % 4 == 0 and year % 100 != 0) or year % 400 == 0
    return 366 if leap else 365


def btime(year, tick):
    """Return the packed BTIME for 'tick' 0.0001s ticks into 'year'."""
    while tick >= daysinyear(year) * SECPERDAY * TICKSPERSEC:
        tick -= daysinyear(year) * SECPERDAY * TICKSPERSEC
        year += 1

    fract = tick % TICKSPERSEC
    second = tick // TICKSPERSEC
    yday = second // SECPERDAY + 1
    second %= SECPERDAY

    return struct.pack(
        ">HHBBBBH",
        year,
        yday,
        second // 3600,
        (second // 60) % 60,
        second % 60,
        0,
        fract,
    )


def channelcodes(index):
    """Return (network, station, location, channel) for channel 'index'."""
    return (
        b"XX",
        b"T%03d " % (index % 1000),
        b"  ",
        b"BH" + bytes([ord("A") + (index % 26)]),
    )


def makerecord(template, seqnum, codes, year, tick):
    record = bytearray(template)

    record[MS2_SEQNUM : MS2_SEQNUM + 6] = b"%06d" % (seqnum % 1000000)
    record[MS2_NETWORK : MS2_NETWORK + 2] = codes[0]
    record[MS2_STATION : MS2_STATION + 5] = codes[1]
    record[MS2_LOCATION : MS2_LOCATION + 2] = codes[2]
    record[MS2_CHANNEL : MS2_CHANNEL + 3] = codes[3]
    record[MS2_YEAR : MS2_YEAR + 10] = btime(year, tick)

    return bytes(record)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("shape", choices=["clean", "scattered", "dense", "layered"])
    parser.add_argument("-l", "--layers", type=int, default=8,
                        help="deliveries for the layered shape (default: 8)")
    parser.add_argument("-r", "--records", type=int, default=20000,
                        help="records per channel (default: 20000)")
    parser.add_argument("-c", "--channels", type=int, default=4,
                        help="number of channels (default: 4)")
    parser.add_argument("-o", "--output", required=True, help="output file")
    parser.add_argument("-y", "--year", type=int, default=2012, help="start year")
    args = parser.parse_args()

    with open(TEMPLATE, "rb") as fp:
        template = fp.read(RECLEN)

    if len(template) != RECLEN:
        sys.exit("%s: short read of template record" % TEMPLATE)

    samples = struct.unpack(">H", template[MS2_NUMSAMPLES : MS2_NUMSAMPLES + 2])[0]

    # Ticks from one record start to the next when contiguous, 247 samples at 40 Hz
    span = samples * TICKSPERSEC // 40

    def layertick(index, layer):
        """Contiguous records within a layer, layers offset so they cannot heal."""
        return index * span + layer * (span // (args.layers + 1) + 13)

    def starttick(index):
        if args.shape == "clean":
            return index * span

        if args.shape == "scattered":
            # Pairs offset by half a record, pairs separated by a full record,
            # so the two records of a pair overlap and nothing else does
            return (index // 2) * 2 * span + (index % 2) * (span // 2)

        # A rolling pile, each record overlapping the surrounding ~20
        return index * (span // 20 + 13)

    layers = args.layers if args.shape == "layered" else 1

    with open(args.output, "wb") as out:
        for index in range(args.records):
            for layer in range(layers):
                for channel in range(args.channels):
                    tick = (layertick (index, layer) if args.shape == "layered"
                            else starttick (index))
                    out.write (makerecord (template, index, channelcodes (channel),
                                           args.year, tick))

    total = args.records * args.channels * layers
    print("Wrote %d records (%d channels x %d x %d layers) of %d bytes to %s"
          % (total, args.channels, args.records, layers, RECLEN, args.output))
    print("Shape '%s': record span %.4f s"
          % (args.shape, span / TICKSPERSEC))


if __name__ == "__main__":
    main()
