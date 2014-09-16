#!/usr/bin/env python
"""Plot ROC-like figures for comparison."""

import os.path as op
import os, sys
import logging
from pbcore.util.ToolRunner import PBToolRunner
from pbcore.util.Process import backticks
from pbove.utils.Utils import realpath, mkdir, cat_files
from pbove.io.RunInfoReader import RunInfo, RunInfoReader
from pbove.pbove_main import add_params_to_parser
from pbove.__init__ import get_version, get_dir

#run_folders = [
#"11063P_4pctFMP/",
#"11063P_8pctFMP/",
#"P5C3/",
#"P4C2/"]

class PbovePlotRuns(object):
    """Class of pbove_plot_runs."""
    def __init__(self, runinfos_fn, plot_dir):

        self.runinfos_fn = realpath(runinfos_fn)
        self.plot_dir = plot_dir
        mkdir(self.plot_dir)

        self.runinfos = []
        with RunInfoReader(self.runinfos_fn) as reader:
            for runinfo in reader:
                runinfo.fofn = realpath(runinfo.fofn)
                runinfo.out_dir = realpath(runinfo.out_dir)
                self.runinfos.append(runinfo)
        self.job_fns = []

    def pbove_out_csv(self, runinfo):
        """Return out table, e.g., "out.csv", for the given run."""
        return op.join(runinfo.out_dir, "out.csv")

    @property
    def R_input(self):
        """return R input."""
        return realpath(op.join(self.plot_dir, "pbove_R_input.txt"))

    def create_R_input(self):
        """Create R input."""
        with open(self.R_input, 'w') as writer:
            writer.write("\t".join(["out_csv", "name", "group"]) + "\n")
            for runinfo in self.runinfos:
                out_csv = self.pbove_out_csv(runinfo)
                if not op.exists(out_csv):
                    raise IOError("{f} does not exists.".format(f=out_csv))
                writer.write("\t".join([out_csv, runinfo.name,
                                        runinfo.group]) + "\n")

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
        logging.info("pbove_plot_runs started.")

        logging.info("pbove_plot_runs will plot figures for {n} runs.".
                     format(n=len(self.runinfos)))
        for i, runinfo in enumerate(self.runinfos):
            logging.debug("[" + str(i) + "]: " + runinfo.fofn +
                          ", " + runinfo.name + ", " + runinfo.group)

        logging.info("Plotting figures.")
        self.plot_figures()

        logging.info("pbove_plot_runs completed.")


def set_parser(parser):
    """Set parser."""
    parser.add_argument("runinfos", type=str,
                        help="A file containing info of jobs/movies " +
                             "to compare. Each line has three fields, " +
                             "fofn of a job, job name and an output " +
                             "directory for analyzing this job.")

    helpstr = "Save generated plots to this directory."
    parser.add_argument("plot_dir", type=str, help=helpstr)

    return parser

class PbovePlotRunner(PBToolRunner):
    """pbove plot runner"""
    def __init__(self):
        desc = "To plot performance of multiple " + \
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
            obj = PbovePlotRuns(
                runinfos_fn=args.runinfos,
                plot_dir=args.plot_dir)
            obj.run()
        except ValueError as e:
            logging.error(str(e))
            return 1
        return 0


def main():
    """Main function."""
    runner = PbovePlotRunner()
    return runner.start()

if __name__ == "__main__":
    sys.exit(main())
