"""Define util functions."""
import os.path as op
import os
import shutil
from pbcore.io import FastaReader, FastaWriter


def realpath(f):
    """Return absolute, user expanded path."""
    if f is None:
        return None
    return op.abspath(op.expanduser(f))


def mkdir(path):
    """Create a directory if it does not pre-exist,
    otherwise, pass."""
    if not op.exists(path):
        os.makedirs(path)


def mknewdir(path):
    """Create a new directory if it does not pre-exist,
    otherwise, delete it and then re-create it."""
    if op.exists(path):
        shutil.rmtree(path)
    os.makedirs(path)


def revcmp(seq):
    """Given a sequence return its reverse complement sequence."""
    NTMAP = {'a':'t', 'c':'g', 't':'a', 'g':'c',
             'A':'T', 'C':'G', 'T':'A', 'G':'C'}
    return "".join([NTMAP[x] for x in seq])[::-1]


def revcmp_fasta(in_fasta, out_fasta):
    """Reverse compelement every reads in in_fasta and
    output it to out_fasta.
    """
    if realpath(in_fasta) == realpath(out_fasta):
        raise ValueError("revcmp_fasta input and output fasta files " +
                         "are identical.")

    num_reads = 0
    with FastaReader(in_fasta) as reader, \
         FastaWriter(out_fasta) as writer:
        for r in reader:
            num_reads += 1
            writer.writeRecord(r.name, revcmp(r.sequence))
    return num_reads


def cat_files(src, dst):
    """Concatenate files in src and save to dst.
       src --- source file names in a list
       dst --- destinate file name
    """
    if src is None or len(src) == 0:
        raise ValueError("src should contain at least one file.")
    if dst in src:
        raise IOError("Unable to cat a file and save to itself.")

    with open (dst, 'w') as writer:
        for src_f in src:
            with open(src_f, 'r') as reader:
                for line in reader:
                    writer.write(line.rstrip() + '\n')

