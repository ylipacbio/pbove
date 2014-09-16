#!/usr/bin/env python
"""Define Class `FastaSplitter` which splits a fasta file into
smaller files each containing `reads_per_split` reads."""

import os
import os.path as op
from pbcore.io import FastaReader, FastaWriter
from pbove.utils.Utils import mkdir

class FastaSplitter(object):
    """An object of `FastaSplitter` splits a fasta file into
    smaller chunks with a given prefix."""
    def __init__(self, input_fasta, reads_per_split, out_dir, out_prefix):
        self.input_fasta = input_fasta
        self.out_dir = out_dir
        self.reads_per_split = reads_per_split  # Number of reads per split
        self.out_prefix = out_prefix
        self.out_fns = None
        mkdir(self.out_dir)

    def __str__(self):
        if self.out_fns is None or len(self.out_fns) == 0:
            return "{input_fasta} ".format(input_fasta=self.input_fasta) + \
                "will be splitted into files each has " + \
                "{n} reads.".format(n=self.reads_per_split)
        else:
            return "{input_fasta} has been splitted into ".\
                   format(input_fasta=self.input_fasta) + \
                   "{m} files each has {n} reads:\n".\
                   format(m=len(self.out_fns),
                n=self.reads_per_split) + ";".join(self.out_fns)

    def _out_fn(self, split_index):
        """Return name of the `split_index`-th splitted file."""
        if split_index > 99:
            raise ValueError("Too many splitted files to generate: number " +
                "of splitted files exceed 100.")
        name = "{prefix}_{idx:02d}.fa".format(prefix=self.out_prefix,
                                              idx=split_index)
        return op.join(self.out_dir, name)

    def split(self):
        """Split `input_fasta` into smaller files each containing
        `reads_per_split` reads. Return splitted fasta."""
        split_index = 0
        self.out_fns = []
        writer = FastaWriter(self._out_fn(split_index))
        self.out_fns.append(self._out_fn(split_index))
        with FastaReader(self.input_fasta) as reader:
            for ridx, r in enumerate(reader):
                if ridx % self.reads_per_split == 0 and ridx != 0:
                    split_index += 1
                    writer.close()
                    writer = FastaWriter(self._out_fn(split_index))
                    self.out_fns.append(self._out_fn(split_index))
                writer.writeRecord(r.name, r.sequence)

        writer.close()
        return list(self.out_fns)

    def rmOutFNs(self):
        """Remove splitted files."""
        for f in self.out_fns:
            os.remove(f)
        self.out_fns = []


def splitFasta(input_fasta, reads_per_split, out_dir, out_prefix):
    """
    Split input_fasta into small fasta files each containing at most
    reads_per_split reads. All splitted fasta files will be placed under
    out_dir with out_prefix. Return paths to splitted files in a list.
    """
    obj = FastaSplitter(input_fasta=input_fasta,
                        reads_per_split=reads_per_split,
                        out_dir=out_dir, out_prefix=out_prefix)
    return obj.split()


def get_args():
    """Get arguments."""
    import argparse
    parser = argparse.ArgumentParser(
        description="Split a fasta file into smaller chunks.")
    parser.add_argument("input_fasta",
                        type=str,
                        help="Input fasta to be splitted.")
    parser.add_argument("reads_per_split",
                        type=int,
                        help="Reads per split.")
    parser.add_argument("out_dir",
                        type=str,
                        help="Output directory.")
    parser.add_argument("out_prefix",
                        type=str,
                        help="Output files prefix.")
    this_args = parser.parse_args()
    return this_args

def main():
    """Main function, split a fasta into smaller chunks."""
    import logging
    from pbtools.pbtranscript.__init__ import get_version
    log = logging.getLogger(__name__)
    args = get_args()
    from pbtools.pbtranscript.Utils import setup_log
    setup_log(alog=log, level=logging.DEBUG)
    log.info("Running {f} v{v}.".format(f=op.basename(__file__),
                                        v=get_version()))

    splitFasta(input_fasta=args.input_fasta,
               reads_per_split=args.reads_per_split,
               out_dir=args.out_dir,
               out_prefix=args.out_prefix)

if __name__ == "__main__":
    import sys
    sys.exit(main())
