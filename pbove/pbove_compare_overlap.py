#!/usr/bin/env python

"""Compare multiple aligners (e.g., different versions of blasr,
or aligner with different parameters)."""

import os.path as op
import sys
import logging
from collections import namedtuple, defaultdict

from pbcore.util.ToolRunner import PBToolRunner
from pbcore.io import FastaReader
from pbalign.utils.fileutil import checkReferencePath

from pbove.io.M4IO import M4Reader
from pbove.utils.Utils import realpath #, mkdir
from pbove.__init__ import get_version
from pbove.utils.Interval import Interval, RefIntervals

RefInfo = namedtuple('RefInfo', ['name', 'len', 'index'])


def div(x, y):
    """return str(x/y)."""
    if x == 0:
        return "0"
    if y == 0:
        return "NA"
    else:
        return "{ret:.5f}".format(ret=(float(x))/y)

def sensitivity(intersect_len, sub_1_2_len, sub_2_1_len, take_1_as_gold=True):
    """To compute sensitivity.
       intersect_len: length of intersected interval
       sub_1_2_len: length of interval in 1 not in 2
       sub_2_1_len: length of interval in 2 not in 1
       take_1_as_gold: True, taking 1 as the gold standard, and compute
                       sensitivity of 2 compared with 1.
                       False: taking 2 as the gold standard, and compute
                       sensitivity of 1 compared with 2.
    """
    sub_len = sub_1_2_len if take_1_as_gold else sub_2_1_len
    return div(intersect_len, intersect_len + sub_len)


def fdr(intersect_len, sub_1_2_len, sub_2_1_len, take_1_as_gold=True):
    """To compute false discovery rate."""
    sub_len = sub_2_1_len if take_1_as_gold else sub_1_2_len
    return div(sub_len, sub_len + intersect_len)


class PboveCompareOverlap(object):
    """Class of pbove_compare_overlap."""
    def __init__(self, query_reads, ref, m4_1, m4_2, out_file):
        self.query_reads = query_reads
        self.ref = ref
        _a, self.ref_fasta, _b, _c, _d = checkReferencePath(self.ref)
        self.m4_1 = realpath(m4_1)
        self.m4_2 = realpath(m4_2)
        self.out_file = realpath(out_file)
        #self.plot_dir = plot_dir
        #mkdir(self.plot_dir)

        self.ref_infos = {}
        self.alns_1 = []
        self.alns_2 = []

    def get_ref_infos(self):
        """Get reference infos, including names, lengths, and
           indices in ref repo.
           ref_infos: dictionary,
                      ref_name -> RefInfo(ref_name, length, index_of_this_ref)
        """
        with FastaReader(self.ref_fasta) as reader:
            self.ref_infos = {r.name.split()[0]:
                              RefInfo(r.name.split()[0], len(r.sequence), idx)
                              for (idx, r) in enumerate(reader)}

    def get_aln_infos(self, m4):
        """From m4, for each query, save its target mapping intervals
        (RefIntervals) in aln_infos and return aln_infos.
        aln_infos: dictionary of RefIntervals,
                   qname -> {tindex -> intervals}
        """
        aln_infos = defaultdict(RefIntervals)
        reader = M4Reader(m4)
        for entry in reader:
            aln_infos[entry.qname] \
                     [self.ref_infos[entry.tname].index] \
                     .add(Interval(entry.abs_tstart, entry.abs_tend))
        reader.close()
        return aln_infos

    def cmp_aln_infos(self, ais_1, ais_2):
        """Taking aln_infos_1 as gold standard, compare aln_infos_2 with
        aln_info_1.
        Write (query name,
               total alignment length of this query in aln_infos_1,
               total alignment length of this query in aln_infos_2,
               total alignment length of this query in aln_infos_1 intersect with 2,
               total alignment length of this query in aln_infos_1, but not in 2,
               total alignment length of this query in aln_infos_2, but not in 1)

        Return (total alignment length of all queries in aln_infos_1,
                total alignment length of all queries in aln_infos_2,
                total alignment length of all queries in aln_infos_1 intersect with 2,
                total alignment length of all queries in aln_infos_1 but not in 2,
                total alignment length of all queries in aln_infos_2 but not in 1)
        """
        total_intersect_len, total_len_1, total_len_2 = 0, 0, 0
        sub_2_1_len, sub_1_2_len = 0, 0

        qnames = set(ais_1.keys()).union(ais_2.keys())

        with open(self.out_file, 'w') as f:
            f.write("#subreads:{sr}\n".format(sr=self.query_reads))
            f.write("#reference:{rf}\n".format(rf=self.ref))
            f.write("#m4_1:{f1}\n".format(f1=self.m4_1))
            f.write("#m4_2:{f2}\n".format(f2=self.m4_2))
            f.write("\t".join(["#qname", "q_aln_len_in_1", "q_aln_len_in_2",
                               "q_aln_intersect_len", "q_aln_not_in_1_len",
                               "q_aln_not_in_2_len", "sensitivity_1_as_gold",
                               "fdr_1_as_gold", "sensitivity_2_as_gold",
                               "fdr_2_as_gold"]) + "\n")
            for qname in qnames:
                logging.debug("Processing {qname}".format(qname=qname))
                ris_1 = ais_1[qname]
                ris_2 = ais_2[qname]

                q_intersect_len_1, q_intersect_len_2 = 0, 0
                q_total_len_1, q_total_len_2 = 0, 0
                q_sub_2_1_len, q_sub_1_2_len = 0, 0

                for tindex, intvs_1 in ris_1.items():
                    # Taking aln_infos_1[qname][tindex] as gold standard,
                    # compare with aln_infos_2[qname][tindex]
                    intersect_intvs = intvs_1.intersect(ris_2[tindex])
                    q_intersect_len_1 += intersect_intvs.length
                    sub_2_1_intvs = ris_2[tindex] - intvs_1
                    q_sub_2_1_len += sub_2_1_intvs.length
                    q_total_len_1 += intvs_1.length

                for tindex, intvs_2 in ris_2.items():
                    # Taking aln_infos_2[qname][tindex] as gold standard,
                    # compare with aln_infos_1[qname][tindex]
                    intersect_intvs = intvs_2.intersect(ris_1[tindex])
                    q_intersect_len_2 += intersect_intvs.length
                    sub_1_2_intvs = ris_1[tindex] - intvs_2
                    q_sub_1_2_len += sub_1_2_intvs.length
                    q_total_len_2 += intvs_2.length

                assert(q_intersect_len_1 == q_intersect_len_2)

                total_intersect_len += q_intersect_len_1
                sub_2_1_len += q_sub_2_1_len
                sub_1_2_len += q_sub_1_2_len
                total_len_1 += q_total_len_1
                total_len_2 += q_total_len_2

                sensitivity_1 = sensitivity(q_intersect_len_1, q_sub_1_2_len,
                                            q_sub_2_1_len, True)
                fdr_1 = fdr(q_intersect_len_1, q_sub_1_2_len,
                            q_sub_2_1_len, True)
                sensitivity_2 = sensitivity(q_intersect_len_1, q_sub_1_2_len,
                                            q_sub_2_1_len, False)
                fdr_2 = fdr(q_intersect_len_1, q_sub_1_2_len,
                            q_sub_2_1_len, False)

                fields = [qname, q_total_len_1, q_total_len_2,
                          q_intersect_len_1, q_sub_2_1_len, q_sub_1_2_len,
                          sensitivity_1, fdr_1, sensitivity_2, fdr_2]
                f.write("\t".join([str(x) for x in fields]) + "\n")

        return (total_len_1, total_len_2,
                total_intersect_len, sub_2_1_len, sub_1_2_len)

    def run(self):
        """Run"""
        logging.info("pbove_compare_overlap started.")

        logging.info("Get reference infos, including names, lengths " +
                     "and indices in refrepo.")
        self.get_ref_infos()

        logging.info("pbove_compare_overlap will compare alignment " +
                     "out files {f1} with {f2}.".
                     format(f1=self.m4_1, f2=self.m4_2))

        ais_1 = self.get_aln_infos(self.m4_1)
        ais_2 = self.get_aln_infos(self.m4_2)

        logging.info("Comparing {f2} and {f1}".
                     format(f1=self.m4_1, f2=self.m4_2))
        (total_len_1, total_len_2, total_intersect_len,
         sub_2_1_len, sub_1_2_len) = self.cmp_aln_infos(ais_1, ais_2)

        sensitivity_1 = sensitivity(total_intersect_len, sub_1_2_len,
                                    sub_2_1_len, True)
        fdr_1 = fdr(total_intersect_len, sub_1_2_len, sub_2_1_len, True)
        sensitivity_2 = sensitivity(total_intersect_len, sub_1_2_len,
                                    sub_2_1_len, False)
        fdr_2 = fdr(total_intersect_len, sub_1_2_len, sub_2_1_len, False)

        output_1 = "Taking {f1} as gold standard: \n".format(f1=self.m4_1) + \
                   "Sensitivity: {sen}\n".format(sen=sensitivity_1) + \
                   "False discovery rate: {fdr}".format(fdr=fdr_1)
        logging.info(output_1)

        output_2 = "Taking {f2} as gold standard: \n".format(f2=self.m4_2) + \
                   "Sensitivity: {sen}\n".format(sen=sensitivity_2) + \
                   "False discovery rate: {fdr}".format(fdr=fdr_2)
        logging.info(output_2)

        output_3 = "Comparing total length of intervals in reference:\n" + \
                   "{t1}\t{t2}\t{ratio}".format(t1=total_len_1,
                   t2=total_len_2, ratio=div(total_len_2, total_len_1))
        logging.info(output_3)

        with open(self.out_file, 'a') as f:
            f.write("#" + output_1 + "\n#" + output_2 + "\n#" + output_3 + "\n")

        logging.info("pbove_compare_overlap completed.")


def set_parser(parser):
    """Set parser."""
    helpstr = "Input query reads in fasta format."
    parser.add_argument("query_fasta", type=str, help=helpstr)

    helpstr = "Reference sequence or reference repository."
    parser.add_argument("ref", type=str, help=helpstr)

    helpstr = "Alignment results in m4 format produced " + \
              "by an aligner (e.g., blasr) aligning reads " + \
              "to the reference sequence."
    parser.add_argument("m4_1", type=str, help=helpstr)

    helpstr = "Alignment results in m4 format produced " + \
              "by a different aligner (e.g., blasr of a " + \
              "different version)."
    parser.add_argument("m4_2", type=str, help=helpstr)

    helpstr = "Output."
    parser.add_argument("out_file", type=str, help=helpstr)

    return parser


class PboveCompareOverlapRunner(PBToolRunner):
    """pbove_compare_runs runner"""
    def __init__(self):
        desc = "Pbove_compare_overlap to compare alignment results of " + \
               "two different aligners."
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
            obj = PboveCompareOverlap(
                query_reads=args.query_fasta,
                ref=args.ref,
                m4_1=args.m4_1,
                m4_2=args.m4_2,
                out_file=args.out_file)
            obj.run()
        except ValueError as e:
            logging.error(str(e))
            return 1
        return 0


def main():
    """Main function."""
    runner = PboveCompareOverlapRunner()
    return runner.start()

if __name__ == "__main__":
    sys.exit(main())


