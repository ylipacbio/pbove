#!/usr/bin/env python
"""Given an input.fofn and a reference sequence,
* call blasr to align reads in fofn to reference sequence
  in order to get ground truth mapping.
* call pls2fasta to get all reads and seed reads which add
  up to 30X of reference and no less than 6000 bp in length
* call blasr to align all reads to seed reads and get
  predicted overlaps
* call pbove_eval module to evaluate sensitivity & specificity ...
"""

import os.path as op
import sys
import logging
from pbcore.io import FastaReader
from pbcore.util.ToolRunner import PBToolRunner
from pbove.__init__ import get_version
from pbove.utils.Utils import realpath, mkdir
from pbove.pbove_filter_subreads import FilterSubreads
from pbove.pbove_reseq import DoReseq
from pbove.pbove_preassembly import DoPreassembly
from pbove.pbove_eval import DoEval
from pbalign.utils.fileutil import checkReferencePath

class Pbove(object):
    """Class of pbove."""
    def __init__(self, input_fofn, ref, out_tb, out_dir,
                 reseq_blasr_opts, preassembly_blasr_opts,
                 force_redo, min_seed_len, split_palindrome,
                 ovl_cut_off, palindrome_score_cutoff,
                 gt_overlaps_file):
        print gt_overlaps_file
        self.input_fofn = realpath(input_fofn)

        self.ref = ref
        self.ref_path, self.ref_fasta, self.ref_sa, self.in_refrepo, \
        _gff = checkReferencePath(self.ref)

        self.out_tb = realpath(out_tb)
        self.out_dir = realpath(out_dir)
        mkdir(self.out_dir)
        self.reseq_blasr_opts = reseq_blasr_opts
        self.preassembly_blasr_opts = preassembly_blasr_opts
        self.force_redo = force_redo
        self.min_seed_len = min_seed_len
        self.split_palindrome = split_palindrome
        self.ovl_cut_off = int(ovl_cut_off)
        self.palindrome_score_cutoff = int(palindrome_score_cutoff)
        self.gt_overlaps_file = realpath(gt_overlaps_file)

    @property
    def reseq_m4(self):
        """Return resequencing m4 output."""
        return op.join(self.out_dir, "reseq_out.m4")

    @property
    def preassembly_m4(self):
        """Return preassembly m4 output."""
        return op.join(self.out_dir, "preassembly_out.m4")

    @property
    def all_reads_fasta(self):
        """Return all reads extracted from input.fofn."""
        return op.join(self.out_dir, "all_reads.fasta")

    @property
    def seed_reads_fasta(self):
        """Return seed reads extracted from all_reads_fasta."""
        return op.join(self.out_dir, "seed_reads.fasta")

    @property
    def ref_sz(self):
        """Return number of bases in reference."""
        with FastaReader(self.ref_fasta) as reader:
            sz = 0
            for r in reader:
                sz += len(r.sequence)
        return sz

    def run(self):
        """Run"""

        logging.info("pbove started.")
        logging.info("Filter subreads from movies.")

        fsubreads = FilterSubreads(input_fofn=self.input_fofn,
                                   all_reads_fasta=self.all_reads_fasta,
                                   out_dir=self.out_dir,
                                   split_palindrome=self.split_palindrome,
                                   nproc=12,
                                   palindrome_score_cutoff=
                                   self.palindrome_score_cutoff)
        fsubreads.run()

        logging.info("resequencing started.")
        reseq = DoReseq(input_reads=self.all_reads_fasta,
                        ref=self.ref,
                        out_m4=self.reseq_m4,
                        blasr_opts=self.reseq_blasr_opts,
                        force_redo=self.force_redo)
        reseq.run()

        logging.info("preassembly started.")
        preassembly = DoPreassembly(all_reads_fasta=self.all_reads_fasta,
                                    seed_reads_fasta=self.seed_reads_fasta,
                                    out_m4=self.preassembly_m4,
                                    ref_sz=self.ref_sz,
                                    out_dir=self.out_dir,
                                    min_seed_len=self.min_seed_len,
                                    blasr_opts=self.preassembly_blasr_opts,
                                    force_redo=self.force_redo)

        preassembly.run()

        logging.info("eval started.")
        evl = DoEval(query_fasta=self.all_reads_fasta,
                     target_fasta=self.seed_reads_fasta,
                     reseq_m4=self.reseq_m4,
                     preassembly_m4=self.preassembly_m4,
                     out_dir=self.out_dir,
                     out_tb=self.out_tb,
                     ovl_cut_off=self.ovl_cut_off,
                     gt_overlaps_file=self.gt_overlaps_file)

        evl.run()

        logging.info("pbove completed.")

def add_params_to_parser(parser):
    """Add resequencing params, preassembly params and
       filter params to parser."""
    defaultstr = "-minMatch 12 -minPctIdentity 70 -minSubreadLength 200 " + \
                 "-bestn 10"
    helpstr = "Advanced blasr parameters for resequencing, not including " + \
              "fixed parameters fixed such as -m, -out, -nproc, " + \
              "-placeRepeatsRandomly, and -useQuality." + \
              "Default: {df}".format(df=defaultstr)
    parser.add_argument("--reseq_blasr_opts", type=str,
                        default=defaultstr, help=helpstr)

    defaultstr = "-bestn 20 -nCandidates 20 -minMatch 10 " + \
                 "-noSplitSubreads -minReadLength 200 -maxLCPLength 16"
    helpstr = "Advanced blasr parameters for preassembly, not including " + \
              "fixed parameters, such as -m, -out, -nproc. " + \
              "Default: {df}".format(df=defaultstr)
    parser.add_argument("--preassembly_blasr_opts", type=str,
                        default=defaultstr, help=helpstr)

    helpstr = "Force to recompute even if outupt files exist."
    parser.add_argument("--force_redo", action="store_true", help=helpstr)

    helpstr = "Minimum seed read length"
    parser.add_argument("--min_seed_len", type=int, default=6000, help=helpstr)

    parser.add_argument("--ovl_cut_off", type=int, default=200,
                        help="Minimum number of overlapping base pairs to " +
                             "consider two reads as positive overlap.")

    parser.add_argument("--split_palindrome", default=False,
                        action='store_true',
                        help="Split palindrome reads (e.g., reads " +
                             "with missing adapters) into two subreads.")

    parser.add_argument("--palindrome_score_cutoff", default=-9000,
                        help="Score cut off to consider a subread palindrome.")
    return parser


def set_parser(parser):
    """Set parser."""
    helpstr = "Input FOFN of bax.h5 files."
    parser.add_argument("input_fofn", type=str, help=helpstr)

    helpstr = "Reference sequence or reference repository."
    parser.add_argument("ref", type=str, help=helpstr)

    helpstr = "Output table."
    parser.add_argument("out_tb", type=str, help=helpstr)

    parser.add_argument('-d', "--out_dir", dest="out_dir",
                        type=str, default="pbove_out", help="Output directory.")

    parser.add_argument('-g', "--ground_truth_overlaps_file",
                        dest="gt_overlaps_file", type=str, default=None,
                        help="Print out ground truth overlap pairs to file.")

    return add_params_to_parser(parser)


class PboveRunner(PBToolRunner):
    """pbove runner"""
    def __init__(self):
        desc = "Pbove to evaluate overlap sensitivity & sepcificity & so on."
        PBToolRunner.__init__(self, desc)
        set_parser(self.parser)

    def getVersion(self):
        """Get version string"""
        return get_version()

    def run(self):
        """Run"""
        logging.info("Running {f} v{v}.".format(f=op.basename(__file__),
                                                v=self.getVersion()))
        args = self.args
        try:
            obj = Pbove(input_fofn=args.input_fofn,
                        ref=args.ref,
                        out_tb=args.out_tb,
                        out_dir=args.out_dir,
                        reseq_blasr_opts=args.reseq_blasr_opts,
                        preassembly_blasr_opts=args.preassembly_blasr_opts,
                        force_redo=args.force_redo,
                        min_seed_len=args.min_seed_len,
                        split_palindrome=args.split_palindrome,
                        ovl_cut_off=args.ovl_cut_off,
                        palindrome_score_cutoff=args.palindrome_score_cutoff,
                        gt_overlaps_file=args.gt_overlaps_file)
            obj.run()
        except ValueError as e:
            logging.error(str(e))
            return 1
        return 0


def main():
    """Main function."""
    runner = PboveRunner()
    return runner.start()

if __name__ == "__main__":
    sys.exit(main())

