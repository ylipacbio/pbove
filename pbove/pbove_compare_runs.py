#!/usr/bin/env python
"""Compare multiple runs (sets of movies) in the context of preassembly and
plot a ROC-like figure."""

import os.path as op
import os, sys
import logging
from pbcore.util.ToolRunner import PBToolRunner
from pbcore.util.Process import backticks
from pbove.utils.Utils import realpath, mkdir, cat_files
from pbove.io.RunInfoReader import RunInfo, RunInfoReader
from pbove.pbove_main import add_params_to_parser
from pbove.__init__ import get_version, get_dir


class PboveComparePreassembly(object):
    """Class of pbove_compare_preassembly."""
    def __init__(self, runinfos_fn, ref, plot_dir,
                 reseq_blasr_opts, preassembly_blasr_opts,
                 force_redo, min_seed_len, split_palindrome,
                 ovl_cut_off, palindrome_score_cutoff):
        self.runinfos_fn = realpath(runinfos_fn)
        self.ref = ref
        self.plot_dir = plot_dir
        mkdir(self.plot_dir)

        self.reseq_blasr_opts = reseq_blasr_opts
        self.preassembly_blasr_opts = preassembly_blasr_opts
        self.force_redo = force_redo
        self.min_seed_len = min_seed_len
        self.split_palindrome = split_palindrome
        self.ovl_cut_off = int(ovl_cut_off)
        self.palindrome_score_cutoff = int(palindrome_score_cutoff)

        self.runinfos = []
        with RunInfoReader(self.runinfos_fn) as reader:
            for runinfo in reader:
                runinfo.fofn = realpath(runinfo.fofn)
                runinfo.out_dir = realpath(runinfo.out_dir)
                self.runinfos.append(runinfo)
        self.job_fns = []

    def pbove_job(self, runinfo):
        """Return pbove_job.sh for the given run."""
        return op.join(runinfo.out_dir, "pbove_job.sh")

    def pbove_log(self, runinfo):
        """Return log for pbove_job.sh for the given run."""
        return op.join(runinfo.out_dir, "pbove_job.log")

    def pbove_out_csv(self, runinfo):
        """Return out table, e.g., "out.csv", for the given run."""
        return op.join(runinfo.out_dir, "out.csv")

    def create_pbove_job(self, runinfo):
        """Create a job script, self.pbove_job(runinfo), for a give run,
        chmod +x, and return path to the job script.
        """
        mkdir(runinfo.out_dir)

        cmd = "#!/bin/bash \n"
        cmd += "pbove -vv {fofn} {ref} {out} --out_dir {od} ".\
               format(fofn=runinfo.fofn,
                      ref=self.ref,
                      out=self.pbove_out_csv(runinfo),
                      od=runinfo.out_dir) + \
               "--reseq_blasr_opts=\"{rbo}\" ".\
               format(rbo=self.reseq_blasr_opts) + \
               "--preassembly_blasr_opts=\"{pbo}\" ".\
               format(pbo=self.preassembly_blasr_opts) + \
               ("--force_redo " if self.force_redo else "") + \
               "--min_seed_len={msl} ".format(msl=self.min_seed_len) + \
               "--ovl_cut_off={oco} ".format(oco=self.ovl_cut_off)

        if self.split_palindrome:
            cmd += "--split_palindrome " + \
                   "--palindrome_score_cutoff={psc} ".\
                   format(psc=self.palindrome_score_cutoff)

        cmd += "2>{log}".format(log=self.pbove_log(runinfo))

        job_fn = self.pbove_job(runinfo)
        with open(job_fn, 'w') as writer:
            writer.write(cmd + "\n")

        os.chmod(job_fn, 0744)
        return job_fn

    @property
    def all_pbove_jobs(self):
        """Return all_pbove_jobs.txt"""
        return realpath("all_pbove_jobs.txt")

    def create_pbove_jobs(self):
        """Create a job script for every run folder in run_folders,
        and save job scripts to file: self.all_pbove_jobs.
        """
        self.job_fns = []
        with open(self.all_pbove_jobs, 'w') as writer:
            for runinfo in self.runinfos:
                logging.debug("Creating a job script at {f}".
                              format(f=self.pbove_job(runinfo)))
                job_fn = self.create_pbove_job(runinfo)
                writer.write(job_fn + "\n")
                self.job_fns.append(job_fn)

        os.chmod(self.all_pbove_jobs, 0744)
        logging.info("Writing all scripts to {f}.".
                     format(f=self.all_pbove_jobs))

    def execute_pbove_jobs(self):
        """Execute pbove jobs."""
        for job_fn in self.job_fns:
            cmd = "/bin/bash {f}".format(f=job_fn)
            logging.debug("CMD: {cmd}".format(cmd=cmd))
            _out, _code, _msg = backticks(cmd)
            if _code != 0:
                raise RuntimeError("Failed to run {cmd}. ".
                                   format(cmd=cmd) +
                                   str(_out) + " " + str(_msg))

    @property
    def R_input(self):
        """return R input."""
        return realpath(op.join(self.plot_dir, "pbove_R_input.txt"))

    def create_R_input(self):
        """Create R input."""
        with open(self.R_input, 'w') as writer:
            writer.write("\t".join(["out_csv", "name", "group"]) + "\n")
            for runinfo in self.runinfos:
                writer.write("\t".join([self.pbove_out_csv(runinfo),
                                       runinfo.name, runinfo.group]) + "\n")

    @property
    def R_script(self):
        """Return R script."""
        return realpath(op.join(self.plot_dir, "pbove.R"))

    def create_R_script(self):
        """Create and execute R scripts"""
        tmp_R = op.join(self.plot_dir, "tmp.R")
        #R_out_dir = op.join(os.getcwd(), self.plot_dir)
        with open(tmp_R, 'w') as writer:
            writer.write("#!/usr/bin/env Rscript\n\n")
            writer.write("rootDir = \"{rd}\"\n".format(rd=os.getcwd()))
            writer.write("outDir  = \"{rod}\"\n".format(rod=self.plot_dir))
            writer.write("all_out_fn = \"{f}\"\n".format(f=self.R_input))

        saved_R = op.join(get_dir(), "R/pbove_compare_runs.R")
        cat_files([tmp_R, saved_R], self.R_script)

    def plot_figures(self):
        """Plot figures using R."""
        logging.debug("Creating R input {f}.".format(f=self.R_input))
        self.create_R_input()

        logging.debug("Creating R script {f}.".format(f=self.R_script))
        self.create_R_script()

        logging.debug("Executing R script.")
        cmd = "Rscript {f}".format(f=self.R_script)
        logging.debug("CMD:{cmd}".format(cmd=cmd))
        _out, _code, _msg = backticks(cmd)
        if _code != 0:
            raise RuntimeError("Failed to run R script {f}. ".
                               format(f=self.R_script) +
                               str(_out) + " " + str(_msg))

    def run(self):
        """Run"""
        logging.info("pbove_compare_runs started.")

        logging.info("pbove_compare_runs will compare {n} runs.".
                     format(n=len(self.runinfos)))
        for i, runinfo in enumerate(self.runinfos):
            logging.debug("[" + str(i) + "]: " + runinfo.fofn +
                          ", " + runinfo.name + ", " + runinfo.group)

        logging.info("Creating scripts for runs.")
        self.create_pbove_jobs()

        logging.info("Executing scripts for runs.")
        self.execute_pbove_jobs()

        logging.info("Plotting figures.")
        self.plot_figures()

        logging.info("pbove_compare_runs completed.")


def set_parser(parser):
    """Set parser."""
    parser.add_argument("runinfos", type=str,
                        help="A file containing info of jobs/movies " +
                             "to compare. Each line has three fields, " +
                             "fofn of a job, job name and an output " +
                             "directory for analyzing this job.")

    helpstr = "Reference sequence or reference repository."
    parser.add_argument("ref", type=str, help=helpstr)

    helpstr = "Save generated plots to this directory."
    parser.add_argument("plot_dir", type=str, help=helpstr)

    return add_params_to_parser(parser)


class PboveComparePreassemblyRunner(PBToolRunner):
    """pbove_compare_runs runner"""
    def __init__(self):
        desc = "Pbove_compare_runs to compare performance of multiple " + \
               "jobs/runs in the context of preassembly overlap prediction."
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
            obj = PboveComparePreassembly(
                runinfos_fn=args.runinfos,
                ref=args.ref,
                plot_dir=args.plot_dir,
                reseq_blasr_opts=args.reseq_blasr_opts,
                preassembly_blasr_opts=args.preassembly_blasr_opts,
                force_redo=args.force_redo,
                min_seed_len=args.min_seed_len,
                split_palindrome=args.split_palindrome,
                ovl_cut_off=args.ovl_cut_off,
                palindrome_score_cutoff=args.palindrome_score_cutoff)
            obj.run()
        except ValueError as e:
            logging.error(str(e))
            return 1
        return 0


def main():
    """Main function."""
    runner = PboveComparePreassemblyRunner()
    return runner.start()

if __name__ == "__main__":
    sys.exit(main())
