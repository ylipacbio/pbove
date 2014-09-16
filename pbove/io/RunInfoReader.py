#!/usr/bin/env python

"""Define class RunInfo and RunInfoReader."""


class RunInfo:
    """Info of a run/job, including fofn of movies, run name,
    and an output directory for analyzing this run."""
    def __init__(self, fofn, name, group, out_dir):
        self.fofn = fofn
        self.name = name
        self.group = group
        self.out_dir = out_dir

    def __str__(self):
        return "fofn: {fofn}\n".format(fofn=self.fofn) + \
               "name: {name}\n".format(name=self.name) + \
               "group: {group}\n".format(group=self.group) + \
               "dir : {out_dir}\n".format(out_dir=self.out_dir)

    @classmethod
    def from_string(cls, line):
        try:
            fds = line.split()
            if len(fds) != 4:
                raise AssertionError("{l} has more than three fields.".
                                     format(l=line))
            return RunInfo(fofn=fds[0], name=fds[1],
                           group=fds[2], out_dir=fds[3])
        except (AssertionError, ValueError) as e:
            msg = "String not recognized as a valid RunInfo record. " + \
                    str(e)
            raise ValueError(msg)


class RunInfoReader:
    """Read all RunInfo objects from a file."""
    def __init__(self, fn):
        self.fn = fn
        try:
            self.infile = open(self.fn, 'r')
        except IOError as e:
            msg = "RunInfoReader: could not read file: " + \
                    self.fn + "\n" + str(e)
            raise IOError(msg)

    def __iter__(self):
        try:
            for line in self.infile:
                line = line.strip()
                if line.startswith("#") or len(line) == 0:
                    continue
                yield RunInfo.from_string(line)
        except (ValueError, IOError) as e:
            raise IOError("Failed to read {fn}. ".format(fn=self.fn) + str(e))

    def __enter__(self):
        return self

    def __close__(self):
        self.infile.close()

    def __exit__(self, exc_type, exc_valuue, traceback):
        self.__close__()



