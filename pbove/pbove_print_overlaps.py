#!/usr/bin/env python

"""
Given a query fasta file, a target fasta file, and a resequencing m4
file produced by aligning the uion of query and target reads to a
reference, print ground truth overlap pairs.
"""
import os.path as op
from pbcore.util.ToolRunner import PBToolRunner
from pbove.__init__ import get_version
from pbove.Reseq import ReseqGroundTruth
from pbove.io.PBIReadFastaHeadIO import PBIReadFastaHeadReader
from pbove.utils.compute import write_gt_overlaps
import sys
import logging


class PrintGTOverlaps(object):
    """Print ground truth overlap pairs, given
       a query fasta file, a target fasta file
       and a resequencing m4 file produced by aligning
       the uion of query and target reads to a reference.
    """
    def __init__(self, query_fasta, target_fasta,
                 qt_ref_reseq_m4, gt_overlaps_file):
        self.query_fasta = query_fasta
        self.target_fasta = target_fasta
        self.qt_ref_reseq_m4 = qt_ref_reseq_m4
        self.gt_overlaps_file = gt_overlaps_file

    def run(self):
        """Run"""
        logging.info("Read resequencing M4 file: {f}".format(f=self.qt_ref_reseq_m4))
        # Get ground truth from qt_ref_reseq_m4.
        gt = ReseqGroundTruth(self.qt_ref_reseq_m4)

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

        logging.info("Writing ground truth overlaps (including weak overlaps) " +
                     "to {f} ...".format(f=self.gt_overlaps_file))
        write_gt_overlaps(query=q, target=t,
                          out_file=self.gt_overlaps_file)
        logging.info("pbove_print_overlaps.py completed.")

def set_parser(parser):
    """Set parser arguments."""
    parser.add_argument("query_fasta", type=str, help="Query reads in Fasta.")
    parser.add_argument("target_fasta", type=str, help="Target reads in Fasta.")
    parser.add_argument("qt_ref_reseq_m4", type=str,
        help="Resequencing results in M4 format, produced by aligning " +
             "the union of query and target reads to a reference genome.")
    parser.add_argument("out_file", type=str,
        help="Output file to save ground truth overlap pairs.")
    return parser


class PrintGTOverlapsRunner(PBToolRunner):
    """runner"""
    def __init__(self):
        desc = "Given a query fasta, a target fasta and a resequencing m4 " + \
               "file produced by aligning the uion of query and target " + \
               "reads to a reference, print ground truth overlap pairs " + \
               "(including weak overlaps) to out_file."
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
            obj = PrintGTOverlaps(query_fasta=args.query_fasta,
                                  target_fasta=args.target_fasta,
                                  qt_ref_reseq_m4=args.qt_ref_reseq_m4,
                                  gt_overlaps_file=args.out_file)
            obj.run()
        except ValueError as e:
            logging.error(str(e))
            return 1
        return 0


def main():
    """Main function."""
    runner = PrintGTOverlapsRunner()
    return runner.start()

if __name__ == "__main__":
    sys.exit(main())

