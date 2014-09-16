#!/usr/bin/env python
"""Given a fasta file of subreads and length of reference sequence,
do:
    (1) obtain the longest subreads which add up to 30X coverage of reference
    (2) align all subreads to the 30X longeset subreads usinb blasr
"""


import os.path as op
import sys
import logging
from pbove.__init__ import get_version
from pbove.utils.Utils import realpath, mkdir
from pbcore.util.ToolRunner import PBToolRunner
from pbcore.util.Process import backticks

class DoPreassembly(object):
    """pbove do preassembly."""
    def __init__(self, all_reads_fasta, seed_reads_fasta,
                 out_m4, ref_sz, out_dir, min_seed_len,
                 blasr_opts, force_redo=False):
        self.all_reads_fasta = realpath(all_reads_fasta)
        self.ref_sz = int(ref_sz)
        self.out_dir = realpath(out_dir)
        self.min_seed_len = min_seed_len
        self.blasr_opts = blasr_opts
        self.force_redo = force_redo

        if not op.exists(self.out_dir):
            mkdir(self.out_dir)

        self.out_m4 = realpath(out_m4) if out_m4 is not None \
                      else op.join(self.out_dir, "preassembly_out.m4")

        self.seed_reads_fasta = realpath(seed_reads_fasta) \
                               if seed_reads_fasta is not None \
                               else op.join(self.out_dir, "seed_reads.fasta")
        self.seed_reads_sa = self.seed_reads_fasta + ".sa"

        self._validate_blasr_opts(self.blasr_opts)

    def _validate_blasr_opts(self, blasr_opts):
        """Validate additional blasr options."""
        not_supported_opts = ["-m ", "-out ", '-nproc ',
                              '-placeRepeatsRandomly']
        for opt in not_supported_opts:
            if opt in blasr_opts:
                raise ValueError("-blasr_opts should not contain {opt}".
                        format(opt=opt))

    def create_seed_reads_fasta(self):
        """Create seed_reads_fasta from all_reads_fasta."""
        logging.info("Start to create seed reads.")
        # First get read lengths of reads in all_reads_fasta
        cmd = "fastalength {all_reads} ".\
              format(all_reads=self.all_reads_fasta) + \
              "| cut -d ' ' -f 1 | sort -nr | " + \
              "awk '{{t+=$1;if(t>={ref_sz}*30){{print $1;exit;}}}}'".\
              format(ref_sz=self.ref_sz)
        logging.info("CMD: {cmd}".format(cmd=cmd))
        _o, _c, _m = backticks(cmd)
        if _c != 0:
            raise RuntimeError("CMD failed. " + str(_o) + ' ' + str(_m))

        min_seed_len = max(int(_o[0] if len(_o) > 0 else 0), self.min_seed_len)

        reads_to_rm = op.join(self.out_dir, "reads_to_rm.txt")
        cmd = "fastalength {all_reads} ".\
              format(all_reads=self.all_reads_fasta) + \
              "| awk '($1 < {len}){{print $2}} ' > ".\
              format(len=min_seed_len) + \
              "{reads_to_rm}".format(reads_to_rm=reads_to_rm)
        logging.info("CMD: {cmd}".format(cmd=cmd))
        _o, _c, _m = backticks(cmd)
        if _c != 0:
            raise RuntimeError("CMD failed. " + str(_o) + ' ' + str(_m))

        cmd = "fastaremove {all_reads} ".\
              format(all_reads=self.all_reads_fasta) + \
              "{reads_to_rm} ".format(reads_to_rm=reads_to_rm) + \
              "> {seed_reads} ".format(seed_reads=self.seed_reads_fasta)
        logging.info("CMD: {cmd}".format(cmd=cmd))
        _o, _c, _m = backticks(cmd)
        if _c != 0:
            raise RuntimeError("CMD failed. " + str(_o) + ' ' + str(_m))

    def align(self):
        """Align all_reads_fasta to seed_reads_fasta"""
        logging.info("Start to align all reads to seed reads")
        if op.exists(self.seed_reads_sa) and self.force_redo is not True:
            msg = "sa file {sa} already exist, skip sawriter.".\
                    format(sa=self.seed_reads_fasta)
            logging.warn(msg)
        else:
            cmd = "sawriter {sa} {fa} -blt 10".format(sa=self.seed_reads_sa,
                                                      fa=self.seed_reads_fasta)
            logging.info("CMD: {cmd}".format(cmd=cmd))
            _o, _c, _m = backticks(cmd)
            if _c != 0:
                raise RuntimeError("CMD failed. " + str(_o) + ' ' + str(_m))

        if op.exists(self.out_m4) and self.force_redo is not True:
            msg = "preasembly output {m4} already exists, skip blasr.".\
                    format(m4=self.out_m4)
            logging.warn(msg)
        else:
            cmd = 'blasr ' + \
                  self.all_reads_fasta + ' ' + \
                  self.seed_reads_fasta + ' ' + \
                  '-m 4 -nproc 12 ' + \
                  '-sa {sa} '.format(sa=self.seed_reads_sa) + \
                  '-out ' + self.out_m4 + ' ' + \
                  self.blasr_opts
            logging.info("CMD: {cmd}".format(cmd=cmd))
            _o, _c, _m = backticks(cmd)
            if _c != 0:
                raise RuntimeError("CMD failed. " + str(_o) + ' ' + str(_m))
        logging.info("Preassembly m4 output done.")

    def run(self):
        """Run"""
        self.create_seed_reads_fasta()

        self.align()


def set_parser(parser):
    """Set parser."""
    helpstr = "Reference sequence size."
    parser.add_argument("ref_sz", type=int, help=helpstr)

    parser.add_argument("--out_dir", type=str, default="pbove_out",
                        help="Output directory")

    helpstr = "Output preassembly results in m4 format."
    parser.add_argument("--out_m4", type=str, default=None, help=helpstr)

    helpstr = "Fasta file containing all subreads from FOFN."
    parser.add_argument("--all_reads_fasta", type=str, default=None,
                        help=helpstr)

    helpstr = "Fasta file containing the seed subreads from FOFN. " + \
              "Obtain the longest subreads which add upt to 30X of " + \
              "reference sequence as seed reads candidates, and then " + \
              "select reads whose length which are longer than the " + \
              "minimum seed length as seeds."
    parser.add_argument("--seed_reads_fasta", type=str, default=None,
                        help=helpstr)

    helpstr = "Minimum seed read length"
    parser.add_argument("--min_seed_len", type=int, default=6000, help=helpstr)

    defaultstr = "-bestn 24 -nCandidates 24 " + \
                 "-noSplitSubreads -minReadLength 200 -maxLCPLength 16"
    helpstr = "Advanced blasr parameters for preassembly, not including " + \
              "fixed parameters, such as -m, -out, -nproc. " + \
              "Default: {df}".format(df=defaultstr)
    parser.add_argument("--blasr_opts", type=str, default=defaultstr,
                        help=helpstr)

    helpstr = "Force to recompute even if outupt files exist."
    parser.add_argument("--force_redo", action="store_true", help=helpstr)

    return parser


class DoPreassemblyRunner(PBToolRunner):
    """pbove do_preassembly runner"""
    def __init__(self):
        desc = "Mimic HGAP Preassembly."
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
            obj = DoPreassembly(all_reads_fasta=args.all_reads_fasta,
                                seed_reads_fasta=args.seed_reads_fasta,
                                out_m4=args.out_m4,
                                ref_sz=args.ref_sz,
                                out_dir=args.out_dir,
                                min_seed_len=args.min_seed_len,
                                blasr_opts=args.blasr_opts,
                                force_redo=args.force_redo)
            obj.run()
        except ValueError as e:
            logging.error(str(e))
            return 1
        return 0


def main():
    """Main function."""
    runner = DoPreassemblyRunner()
    return runner.start()

if __name__ == "__main__":
    sys.exit(main())

