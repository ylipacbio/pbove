#!/usr/bin/env python

"""Given a FOFN of bax.h5 files and a reference, do resequencing
(i.e.,align all subreads in FOFN to reference in order to get mapping
'ground truth')."""


import os.path as op
import sys
import logging
from pbove.__init__ import get_version
from pbove.utils.Utils import realpath
from pbcore.util.ToolRunner import PBToolRunner
from pbcore.util.Process import backticks
from pbalign.utils.fileutil import checkReferencePath


class DoReseq(object):
    """pbove do resequencing."""
    def __init__(self, input_reads, ref, out_m4,
                 blasr_opts="", force_redo=False):
        self.input_reads = realpath(input_reads)
        self.ref = realpath(ref)
        self.out_m4 = realpath(out_m4)
        self.ref_path, self.ref_fasta, self.ref_sa, self.in_refrepo, \
        _gff = checkReferencePath(self.ref)
        self.blasr_opts = blasr_opts
        self.force_redo = force_redo
        self._validate_blasr_opts(self.blasr_opts)

    def _validate_blasr_opts(self, blasr_opts):
        """Validate additional blasr options."""
        not_supported_opts = ["-m ", "-out ", '-nproc ',
                              '-placeRepeatsRandomly', '-useQuality']
        for opt in not_supported_opts:
            if opt in blasr_opts:
                raise ValueError("-blasr_opts should not contain {opt}".
                        format(opt=opt))

    def run(self):
        """Run"""
        cmd = 'blasr ' + \
               self.input_reads + ' ' + \
               self.ref_fasta + ' ' + \
               '-m 4 -nproc 12 -placeRepeatsRandomly ' + \
               '-out ' + self.out_m4 + ' ' + \
               ('' if self.ref_sa is None else \
                '-sa {sa} '.format(sa=self.ref_sa)) + \
               self.blasr_opts

        if op.exists(self.out_m4) and self.force_redo is False:
            msg = "Output m4 file {out} exists! Skip blasr ... ".\
                    format(out=self.out_m4)
            logging.warn(msg)
        else:
            logging.info("CMD: {cmd}".format(cmd=cmd))
            _o, _c, _m = backticks(cmd)
            if _c != 0:
                raise RuntimeError("CMD failed: " + str(_o) + ' ' + str(_m))


def set_parser(parser):
    """Set parser."""
    helpstr = "Input reads in FASTA, bax.h5 or FOFN of bax.h5 formats."
    parser.add_argument("input_reads", type=str, help=helpstr)

    helpstr = "Reference sequence or reference repository."
    parser.add_argument("ref", type=str, help=helpstr)

    helpstr = "Output resequencing results in m4 format."
    parser.add_argument("out_m4", type=str, help=helpstr)

    defaultstr = "-minMatch 12 -minPctIdentity 70 -minSubreadLength 200 " + \
                 "-bestn 10"
    helpstr = "Advanced blasr parameters, not including parameters fixed " + \
              "for resequencing, such as -m, -out, -nproc, " + \
              "-placeRepeatsRandomly, and -useQuality." + \
              "Default: {df}".format(df=defaultstr)
    parser.add_argument("--blasr_opts", type=str, default=defaultstr,
                        help=helpstr)

    helpstr = "Force to recompute even if outupt files exist."
    parser.add_argument("--force_redo", action="store_true", help=helpstr)

    return parser


class DoReseqRunner(PBToolRunner):
    """pbove do_reseq runner"""
    def __init__(self):
        desc = "Mimic RS_Resequencing protocol."
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
            obj = DoReseq(input_reads=args.input_reads,
                          ref=args.ref,
                          out_m4=args.out_m4,
                          blasr_opts=args.blasr_opts,
                          force_redo=args.force_redo)
            obj.run()
        except ValueError as e:
            logging.error(str(e))
            return 1
        return 0


def main():
    """Main function."""
    runner = DoReseqRunner()
    return runner.start()

if __name__ == "__main__":
    sys.exit(main())

