#!/usr/bin/env python

"""Given a FOFN of bax.h5 files, filter subreads from movies."""

import os
import os.path as op
import sys
import logging
import shutil
from multiprocessing import Pool
from pbove.__init__ import get_version
from pbove.utils.Utils import realpath, mkdir, revcmp_fasta, cat_files
from pbove.io.SDPReader import SDPReader
from pbove.io.FastaSplitter import FastaSplitter
from pbcore.util.ToolRunner import PBToolRunner
from pbcore.io import FastaReader, FastaWriter
from pbcore.util.Process import backticks


class FilterSubreads(object):
    """pbove filter subreads."""
    def __init__(self, input_fofn, all_reads_fasta, out_dir,
                 force_redo=False, split_palindrome=False,
                 nproc=12, palindrome_score_cutoff=-9000):
        self.input_fofn = realpath(input_fofn)
        self.all_reads_fasta = realpath(all_reads_fasta)
        self.force_redo = force_redo
        self.split_palindrome = split_palindrome
        self.out_dir = realpath(out_dir)
        if not op.exists(self.out_dir):
            mkdir(self.out_dir)

        self.filtered_region_dir = op.join(self.out_dir, "filtered_region")
        self.region_fofn = op.join(self.filtered_region_dir,
                                   "filtered_regions.fofn")
        self.nproc = int(nproc)
        self.palindrome_score_cutoff = int(palindrome_score_cutoff)

    def _filter_subreads(self):
        """Filter subreads from input_fofn using pls2fasta, and create
        all_reads_fasta."""
        logging.info("Start to filter subreads in fofn.")
        if op.exists(self.ori_all_reads_fasta) and self.force_redo is not True:
            msg = "{fa} already exists, skip pls2fasta".format(fa=self.ori_all_reads_fasta)
            logging.warn(msg)
        else:
            logging.debug("{f} does not exist, call pls2fasta".
                          format(f=self.ori_all_reads_fasta))
            filter_summary = op.join(self.filtered_region_dir,
                                     "filtered_summary.csv")
            cmd = "filter_plsh5.py --debug " + \
                  "--filter='MinReadScore=0.80,MinSRL=500,MinRL=100' " + \
                  "--trim='True' --outputDir={fr} ".format(
                          fr=self.filtered_region_dir) + \
                  "--outputSummary={sm} ".format(sm=filter_summary) + \
                  "--outputFofn={rgn} ".format(rgn=self.region_fofn) + \
                  "{in_fofn}".format(in_fofn=self.input_fofn)
            logging.info("CMD: {cmd}".format(cmd=cmd))
            _o, _c, _m = backticks(cmd)
            if _c != 0:
                raise RuntimeError("CMD failed. " + str(_o) + ' ' + str(_m))

            cmd = "pls2fasta -trimByRegion " + \
                  "-regionTable {rgn} ".format(rgn=self.region_fofn) + \
                  "{fofn} {fa} ".format(fofn=self.input_fofn,
                                        fa=self.ori_all_reads_fasta)
            logging.info("CMD: {cmd}".format(cmd=cmd))
            _o, _c, _m = backticks(cmd)
            if _c != 0:
                raise RuntimeError("CMD failed. " + str(_o) + ' ' + str(_m))
            logging.info("{f} created.".format(f=self.ori_all_reads_fasta))

        logging.debug("Copying {ori_f} to {f}.".format(
            ori_f=self.ori_all_reads_fasta, f=self.all_reads_fasta))
        shutil.copyfile(self.ori_all_reads_fasta, self.all_reads_fasta)

    @property
    def rc_all_reads_fasta(self):
        """Return a fasta file containing the reverse complement sequence of
        each read in self.all_reads_fasta."""
        return op.join(self.out_dir, "rc_all_reads.fasta")

    @property
    def sdp_out_file(self):
        """Return output file of sdpMatcher aligning all_reads_fasta to
        rc_all_reads_fasta."""
        return op.join(self.out_dir, "read_to_rc_read.sdp")

    @property
    def tmp_all_reads_fasta(self):
        """Return tmp file for all_reads_fasta."""
        return op.join(self.out_dir, "all_reads.fasta.tmp")

    @property
    def palindrome_reads_fasta(self):
        """Return palindrome reads to palindrome_subreads.fasta."""
        return op.join(self.out_dir, "palindrome_subreads.fasta")

    @property
    def ori_all_reads_fasta(self):
        """Return original file for all_reads_fasta."""
        return op.join(self.out_dir, "all_reads.fasta.ori")

    def _self_align(self):
        """Call SDPMatcher to align every read to its reverse
        complementary reads."""
        logging.info("Splitting palindrome.")
        logging.debug("Making reverse complement sequences of reads in " +
                      "{i} to {o}".format(i=self.ori_all_reads_fasta,
                                          o=self.rc_all_reads_fasta))
        num_reads = revcmp_fasta(self.ori_all_reads_fasta,
                                 self.rc_all_reads_fasta)

        reads_per_split = max(1, int(num_reads/self.nproc) + 1)
        logging.debug("Splitting {f} to small files each containing {n} reads.".
                      format(f=self.ori_all_reads_fasta, n=reads_per_split))
        fs = FastaSplitter(input_fasta=self.ori_all_reads_fasta,
                           reads_per_split=reads_per_split,
                           out_dir=self.out_dir,
                           out_prefix="reads.split.")
        fs.split()
        sp_fasta_files = fs.out_fns

        logging.debug("Splitting {f} to smaller files.".
                      format(f=self.rc_all_reads_fasta))
        rc_fs = FastaSplitter(input_fasta=self.rc_all_reads_fasta,
                              reads_per_split=reads_per_split,
                              out_dir=self.out_dir,
                              out_prefix="rc_reads.split.")
        rc_fs.split()
        rc_sp_fasta_files = rc_fs.out_fns

        logging.debug("Aligning each read in {i} to its revese compelement " +
                      "read using sdpMatcher.".format(i=self.ori_all_reads_fasta))

        sdps = ["{f}.sdp".format(f=f) for f in sp_fasta_files]
        jobs = []
        for f, rc_f, sdp in zip(sp_fasta_files, rc_sp_fasta_files, sdps):
            cmd = "sdpMatcher {f} {rc_f} ".format(f=f, rc_f=rc_f) + \
                  "10 -local > {sdp} ".format(sdp=sdp)
            logging.debug("CMD: {cmd}".format(cmd=cmd))
            jobs.append(cmd)

        pool = Pool(processes=self.nproc)
        rets = pool.map(backticks, jobs)
        pool.close()
        pool.join()

        for i, job in enumerate(jobs):
            if rets[i][1] != 0:
                errMsg = "Job {j} failed.".format(j=job) + str(rets[i][2])
                raise RuntimeError(errMsg)

        logging.debug("Concatenating all sdp outputs to {f}".
                      format(f=self.sdp_out_file))
        cat_files(src=sdps, dst=self.sdp_out_file)

        logging.debug("Cleaning intermediate fasta & sdp files.")
        fs.rmOutFNs()
        rc_fs.rmOutFNs()

        for f in sdps:
            os.remove(f)

    def _split_palindrome(self):
        """There exist some chimeric reads in which adapters are either missing
        or not recognizable. These are called palindrome reads conform to the
        following template:
               read NNN rev_comp(read)
        In order to remove chimeras, we will align each read to its reverse
        compelementary sequence using sdpMatcher. If both forward and backward
        alignments are found for a read (i.e., alignments with tstrand 0
        and 1 both exist), then we will split this read from middle of the query
        range.
        The side effect of this process is that true plindrome reads will be
        cut short.
        """
        if not op.exists(self.sdp_out_file) or self.force_redo is True:
            self._self_align()

        logging.debug("Parsing sdp and detect plindrome reads")
        split_table = {}
        with SDPReader(self.sdp_out_file) as reader:
            for sdp in reader:
                if sdp.score <= self.palindrome_score_cutoff:
                    split_table[str(sdp.qID)] = sdp

        logging.debug("Splitting palindrom reads.")
        with FastaReader(self.ori_all_reads_fasta) as reader, \
             FastaWriter(self.tmp_all_reads_fasta) as writer, \
             FastaWriter(self.palindrome_reads_fasta) as palindrome_writer:
            for r in reader:
                if r.name in split_table:
                    # found a palindrome
                    sdp = split_table[r.name]
                    # Write palindrome subreads to palindrome_subreads.fasta
                    palindrome_writer.writeRecord(r.name, r.sequence)
#
#                    # split this read in the middle
#                    split_point = int(sdp.qstart +
#                                      (sdp.alnqstart + sdp.alnqend)/2)
#                    # Write the first half
#                    rname_1 = "{movie}/{zmw}/{s}_{e}".format(
#                              movie=sdp.movie, zmw=sdp.zmw, s=sdp.qstart,
#                              e=split_point)
#                    writer.writeRecord(rname_1,
#                                       r.sequence[0:(split_point-sdp.qstart)])
#
#                    # Write the second half
#                    rname_2 = "{movie}/{zmw}/{s}_{e}".format(
#                              movie=sdp.movie, zmw=sdp.zmw,
#                              s=(split_point+1), e=sdp.qend)
#                    writer.writeRecord(rname_2,
#                                       r.sequence[(split_point-sdp.qstart):])
                else:
                    writer.writeRecord(r.name, r.sequence)

        logging.debug("Moving {i} to {o}.".format(i=self.tmp_all_reads_fasta,
                                                  o=self.all_reads_fasta))
        shutil.move(self.tmp_all_reads_fasta, self.all_reads_fasta)

    def run(self):
        """Run"""
        self._filter_subreads()

        if self.split_palindrome:
            self._split_palindrome()


def set_parser(parser):
    """Set parser."""
    helpstr = "Input FOFN of bax.h5."
    parser.add_argument("input_fofn", type=str, help=helpstr)

    parser.add_argument('-d', "--out_dir", dest="out_dir",
                        type=str, default="pbove_out", help="Output directory.")

    parser.add_argument("all_reads_fasta", type=str,
                        help="Filtered subreads in FASTA.")

    helpstr = "Force to recompute even if outupt files exist."
    parser.add_argument("--force_redo", action="store_true", help=helpstr)

    parser.add_argument("--split_palindrome", default=False,
                        action='store_true',
                        help="Split palindrome reads (e.g., reads with " +
                        "missing adapters) into two subreads.")
    parser.add_argument("--palindrome_score_cutoff", default=-9000,
                        help="Score cut off to consider a subread palindrome.")
    parser.add_argument("--nproc", default=12,
                        help="Number of threads to call SDPMatcher.")
    return parser


class FilterSubreadsRunner(PBToolRunner):
    """pbove filter subreads runner"""
    def __init__(self):
        desc = "Filter subreads from movies."
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
            obj = FilterSubreads(input_fofn=args.input_fofn,
                                 out_dir=args.out_dir,
                                 all_reads_fasta=args.all_reads_fasta,
                                 force_redo=args.force_redo,
                                 split_palindrome=args.split_palindrome,
                                 nproc=args.nproc,
                                 palindrome_score_cutoff=
                                 args.palindrome_score_cutoff)
            obj.run()
        except ValueError as e:
            logging.error(str(e))
            return 1
        return 0


def main():
    """Main function."""
    runner = FilterSubreadsRunner()
    return runner.start()

if __name__ == "__main__":
    sys.exit(main())

