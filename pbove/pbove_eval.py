#!/usr/bin/env python

"""
Evaluate overlap detection of preassembly.
"""
import os.path as op
from pbcore.util.ToolRunner import PBToolRunner
from pbove.__init__ import get_version
from pbove.Reseq import ReseqGroundTruth
from pbove.utils.Utils import mkdir
from pbove.Preassembly import PreassemblyPrediction
from pbove.io.PBIReadFastaHeadIO import PBIReadFastaHeadReader
from pbove.utils.compute import ComputeAllPosNegNumbers, \
        write_gt_overlaps
import pbove.QTSO as QTSO
import sys
import logging

class Summary(object):
    """Brief summary"""
    def __init__(self):
        self.numQ, self.numT, \
        self.numGTPos, self.numGTNeg, \
        self.numGTWeak, self.numUnmappableAlns, \
        self.numMappableAlns, self.numAlns = 0, 0, 0, 0, 0, 0, 0, 0

    def __str__(self):
        ret = "Number of query reads: " + str(self.numQ) + "\n" + \
            "Number of target reads: " + str(self.numT) + "\n" + \
            "Number of ground truth positive overlaps: " + \
            str(self.numGTPos) + "\n" + \
            "Number of ground truth negative overlaps: " + \
            str(self.numGTNeg) + "\n" + \
            "Number of ground truth weak overlaps: " + \
            str(self.numGTWeak) + "\n" + \
            "Number of ground truth unmappable pairs of reads: " + \
            str(self.numUnmappableAlns) + "\n" + \
            "Number of ground truth mappable pairs of reads: " + \
            str(self.numMappableAlns) + "\n" + \
            "Total number of pairs of reads: " + str(self.numAlns)
        return ret


class DoEval(object):
    """
    Evaluate the sensitivity and sepcificity of overlap detection
    in preassembly.
    """
    def __init__(self, query_fasta, target_fasta, reseq_m4, preassembly_m4,
                 out_dir, out_tb, out_qtso=None, out_dtb=None,
                 ovl_cut_off=200, gt_overlaps_file=None):
        self.query_fasta = query_fasta
        self.target_fasta = target_fasta
        self.reseq_m4 = reseq_m4
        self.preassembly_m4 = preassembly_m4
        self.out_dir = out_dir
        self.out_tb = out_tb
        self.out_qtso = out_qtso if out_qtso is not None else \
                        (op.join(self.out_dir, "out.qtso"))
        self.out_dtb = out_dtb if out_dtb is not None else \
                        (op.join(self.out_dir, "out.delta"))
        self.ovl_cut_off = int(ovl_cut_off)
        self.gt_overlaps_file = gt_overlaps_file

        mkdir(self.out_dir)
        self.summary = Summary()
        self.summary_f = op.join(out_dir, "summary.txt")

    def run(self):
        """Run"""
        logging.info("Read resequencing M4 file: {f}".format(f=self.reseq_m4))
        # Get ground truth from reseq_m4.
        gt = ReseqGroundTruth(self.reseq_m4)

        logging.info("Get query reads from {f}".format(f=self.query_fasta))
        # Get query reads from query_fasta.
        queryReads = PBIReadFastaHeadReader(self.query_fasta)

        logging.info("Find positions of query reads in coordinate of " +
                     "reference genome.")
        q = gt.MapPBISubreadsToReference(queryReads.reads)

        # Get target reads from target_fasta.
        logging.info("Get target reads from {f}".format(f=self.target_fasta))
        targetReads = PBIReadFastaHeadReader(self.target_fasta)

        logging.info("Find positions of target reads in coordinate of " +
                     "reference genome.")
        t = gt.MapPBISubreadsToReference(targetReads.reads)

        # Compute number of
        #     ground truth positive overlap,
        #     ground truth negative overlap,
        #     ground truth weak overlap
        #     un-mappable query reads
        #     mappable query reads
        #     alignments
        # positive: query and target overlap by at least
        #           'OverlapLengthCutoff' bases in coordinate of reference.
        # negative: query and target do not overlap in coordinate of reference.
        # weak: query and target overlap length >= 0, < OverlapLengthCutoff.
        #
        logging.info("Computing numbers of ground truth posivitive, negative.")
        (numGTPos, numGTNeg, numGTWeak, numUnmappableAlns, numMappableAlns,
         numAlns) = ComputeAllPosNegNumbers(q, t, self.ovl_cut_off)

        self.summary.numQ = len(queryReads.reads)
        self.summary.numT = len(targetReads.reads)
        self.summary.numGTPos = numGTPos
        self.summary.numGTNeg = numGTNeg
        self.summary.numGTWeak = numGTWeak
        self.summary.numUnmappableAlns = numUnmappableAlns
        self.summary.numMappableAlns = numMappableAlns
        self.summary.numAlns = numAlns

        with open(self.summary_f, 'w') as writer:
            writer.write(str(self.summary) + "\n")

        if (self.gt_overlaps_file is not None):
            logging.info("Writing ground truth overlap pairs to {f}".
                         format(f=self.gt_overlaps_file))
            write_gt_overlaps(query=q, target=t,
                              out_file=self.gt_overlaps_file)

        # Query-target overlap relationships read from preassembly_m4.
        logging.info("Reading overlap relations from {f}".
                     format(f=self.preassembly_m4))
        pred = PreassemblyPrediction(self.preassembly_m4)

        logging.info("Retrieve overlap lengths from ground truth.")
        pred.OverlapLengthsInReference(gt, infer=True)

        logging.info("Write QTSO info to {f}.".format(f=self.out_qtso))
        pred.ToQTSO(self.out_qtso)

        qtso = QTSO.QTSO(self.out_qtso)

        logging.info("Write delta table to {f}.".format(f=self.out_dtb))
        qtso.getDeltaTable(stepSize=100,
                           outfile=self.out_dtb,
                           overlapLengthCutoff=self.ovl_cut_off)

        logging.info("Write output to {f}.".format(f=self.out_tb))
        qtso.getTable(numGTPos=numGTPos, numGTNeg=numGTNeg,
                numGTWeak=numGTWeak, numUnmappableAlns=numUnmappableAlns,
                numMappableAlns=numMappableAlns,
                numAlns=numAlns, outfile=self.out_tb)


def set_parser(parser):
    """Set parser arguments."""
    parser.add_argument("query_fasta", type=str, help="Query reads in Fasta.")
    parser.add_argument("target_fasta", type=str, help="Target reads in Fasta.")
    parser.add_argument("reseq_m4", type=str,
        help="Resequencing results in M4 format, produced by aligning " +
             "query_fasta to reference genome using BLASR.")

    parser.add_argument("preassembly_m4", type=str,
        help="Overlapping results in M4 format, produced by aligning " +
             "query_fasta to target_fasta using BLASR.")

    parser.add_argument("out_tb", type=str, help="Output table results.")

    parser.add_argument("--ovl_cut_off", type=int, default=200,
        help="Minimum number of overlapping base pairs to consider two " +
             "reads as positive overlap.")

    parser.add_argument('-d', "--out_dir", dest="out_dir",
                        type=str, default="pbove_out", help="Output directory")

    parser.add_argument("--out_qtso", type=str, default=None,
        help="Each line is a tab-delimited record of a query read, a " +
             "target read, blasr score and ground truth overlap length.")

    parser.add_argument("--out_dtb", type=str, default=None,
        help="Delta results in a table.")
    return parser


class DoEvalRunner(PBToolRunner):
    """pbove eval runner"""
    def __init__(self):
        desc = "Given resequencing m4 output and preassembly m4 output, " + \
               "evaluate sensitivity & sepcificity of overlap detection."
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
            obj = DoEval(query_fasta=args.query_fasta,
                         target_fasta=args.target_fasta,
                         reseq_m4=args.reseq_m4,
                         preassembly_m4=args.preassembly_m4,
                         out_dir=args.out_dir,
                         out_tb=args.out_tb,
                         out_qtso=args.out_qtso,
                         out_dtb=args.out_dtb,
                         ovl_cut_off=args.ovl_cut_off)
            obj.run()
        except ValueError as e:
            logging.error(str(e))
            return 1
        return 0


def main():
    """Main function."""
    runner = DoEvalRunner()
    return runner.start()

if __name__ == "__main__":
    sys.exit(main())

