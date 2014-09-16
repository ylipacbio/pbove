"""
Module for parsing simple blasr m4 output.

M4 format:
Each line contains the following 12 fields.
      Field
0     qname
1     tname
2     score
3     pctSimilarity
4     qstrand
5     qstart
6     qend
7     qseqlength
8     tstrand
9     tstart
10    tend
11    tseqlength
12    mapqv
"""
import os
import sys
from operator import  attrgetter
from pbove.utils.PBIReadNameUtils import PBISubsubreadName, PBISubreadName

M4DELIMITER = " "
M4FIELDS = ["qname", "tname", "score", "pctsimilarity",
            "qstrand", "qstart", "qend", "qseqlength",
            "tstrand", "tstart", "tend", "tseqlength",
            "mapqv"]
M4HEADER = M4DELIMITER.join (M4FIELDS)

def parseStrand(strand):
    """Return 0 if positive strand, 1 if negative."""
    if (strand == "+" or strand == "0"):
        return '+'
    elif (strand == "-" or strand == "1"):
        return '-'
    else:
        assert False, "Failed to parse strand {0}.\n".format(strand)

class M4Entry:
    """Storage class for alignment hit records in a M4 file.
    """
    def __init__( self, line, delimiter=" "):
        (self.qname,   self.tname,  self.score, self.pctSimilarity,
         self.qstrand, self.qstart, self.qend,  self.qseqlength,
         self.tstrand, self.tstart, self.tend,  self.tseqlength,
         self.mapqv) = self.__parseLine( line, delimiter )
        # For an M4 entry whose target strand is "-", always map
        # target start and end in '-' strand to end and start in
        # the '+' strand.
        if (self.tstrand == "-"):
            self.abs_tend, self.abs_tstart = self.tseqlength - self.tstart, \
                                     self.tseqlength - self.tend
        else:
            self.abs_tstart, self.abs_tend = self.tstart, self.tend
        # Assume qstrand is always "+"
        assert(self.qstrand == "+")

        try:
        # By default, an M4 entry's query is a PBI sub-subread
            if len(self.qname.split("/")) == 3:
                self.qpbi = PBISubreadName(self.qname)
                self.abs_qstart = self.qstart
                self.abs_qend = self.qend
            elif len(self.qname.split("/")) == 4:
                self.qpbi = PBISubsubreadName(self.qname)
                self.abs_qstart = self.qstart + self.qpbi.start
                self.abs_qend = self.qend + self.qpbi.start
            else:
                raise Exception, "Could not recognize {0} as a PacBio read."\
                        .format(self.qname)
        except Exception as e:
            raise Exception, "Could not parse read name {0}.\n"\
                   .format(self.qname) + str(e)

        try:
        # For an M4 entry whose target is also a PBI subread,
        # map target start and end's coordinate relative to the
        # full length zmw read.
            self.tpbi = PBISubreadName(self.tname)
            self.abs_tstart = self.abs_tstart + self.tpbi.start
            self.abs_tend = self.abs_tend + self.tpbi.start
        except ValueError:
            pass
        except IndexError:
            pass

    def __str__( self ):
        ret = "qname = {0}, tname = {1}, score = {2}," \
              "pctSimilarity = {3}\n".format(self.qname, self.tname,
              self.score, self.pctSimilarity)
        ret += "qstrand = {0}, qstart = {1}, qend = {2}, qseqlength = {3}\n" \
               .format(self.qstrand, self.qstart, self.qend, self.qseqlength)
        ret += "tstrand = {0}, tstart = {1}, tend = {2}, tseqlength = {3}\n" \
               .format(self.tstrand, self.tstart, self.tend, self.tseqlength)
        ret += "abs_qstart (in full length read) = {0}, abs_qend = {1}\n" \
               .format(self.abs_qstart, self.abs_qend)
        ret += "abs_tstart (in + strand target) = {0}, abs_tend = {1}\n\n" \
               .format(self.abs_tstart, self.abs_tend)
        return ret

    def __parseLine(self, line, delimiter):
        """Parse a line in a M4 file."""
        try:
            (qname,   tname,  score, pctSimilarity,
             qstrand, qstart, qend,  qseqlength,
             tstrand, tstart, tend,  tseqlength,
             mapqv) = line.rstrip("\n").split(delimiter)
            assert(len(mapqv)>0)
            score = int(score)
            pctSimilarity = float(pctSimilarity)
            qstrand = parseStrand(qstrand)
            qstart, qend, qseqlength = int(qstart), int(qend), int(qseqlength)
            assert(qstart >= 0 and qend >= 0 and qseqlength > 0 and
                   qstart <= qseqlength and qend <= qseqlength)
            tstrand = parseStrand(tstrand)
            tstart, tend, tseqlength = int(tstart), int(tend), int(tseqlength)
            assert(tstart >= 0 and tend >= 0 and tseqlength > 0 and
                   tstart <= tseqlength and tend <= tseqlength)
            mapqv = int(mapqv)
            return (qname,   tname,  score, pctSimilarity,
                    qstrand, qstart, qend,  qseqlength,
                    tstrand, tstart, tend,  tseqlength,
                    mapqv)
        except (AssertionError, ValueError):
            raise AssertionError("{l} is an incorrect M4 record.".format(l=line))

    def InferAPointFromQToT(self, point):
        """If a point is out of range [abs_qstart, abs_qend], infer its
           mapped point in range [abs_tstart, abs_tend]. Otherwise, return
           MapAPointFromQToT(point).
        """
        if (self.tstrand == "+"):
            return int ((point - self.abs_qstart) *
                     float(self.abs_tend - self.abs_tstart) /
                     float(self.abs_qend - self.abs_qstart)  +
                     self.abs_tstart)
        elif (self.tstrand == "-"):
            return int (self.abs_tend -
                     (self.abs_tend - self.abs_tstart) *
                     float(point - self.abs_qstart) /
                     float(self.abs_qend - self.abs_qstart))
        else :
            raise ValueError, "Unknown tstrand {0}.".format(self.tstrand)
        assert(self.qstrand == "+")

    def MapAPointFromQToT(self, point, infer=False):
        """If a point is in range [abs_qstart, abs_qend], return its mapped
           point in range [abs_tstart, abs_tend]. Otherwise, return -1."""
        if (not infer):
            if not (point >= self.abs_qstart and point <= self.abs_qend):
                return -1
        return self.InferAPointFromQToT(point)

    def MapAnIntervalFromQToT(self, spoint, epoint, infer=False):
        """Given an interval [spoint, epoint] in the coordinate of query, return
           an interval in the target coordinate where [spoint, epoint) maps to.
           If infer is False, the returned interval should be within the
           target range of this read (i.e. [self.abs_tstart, self.abs_tend)).
           If infer is True, the returned interval can be outside of the
           target range.
        """
        mappedSpoint = self.MapAPointFromQToT(spoint, infer)
        mappedEpoint = self.MapAPointFromQToT(epoint, infer)

        if (not infer):
            if (mappedSpoint == -1 and mappedEpoint != -1):
                if (self.tstrand == "+"):
                    mappedSpoint = self.abs_tstart
                else:
                    mappedSpoint = mappedEpoint
                    mappedEpoint = self.abs_tend
            elif (mappedSpoint != -1 and mappedEpoint == -1):
                if (self.tstrand == "+"):
                    mappedEpoint = self.abs_tend
                else:
                    mappedEpoint = mappedSpoint
                    mappedSpoint = self.abs_tstart
            elif (mappedSpoint == -1 and mappedEpoint == -1):
                if (spoint < self.abs_qstart and epoint > self.abs_qend):
                    mappedSpoint, mappedEpoint = self.abs_tstart, self.abs_tend
            else:
                if (self.tstrand == "-"):
                    tmp = mappedSpoint
                    mappedSpoint = mappedEpoint
                    mappedEpoint = tmp
        else:
            if (self.tstrand == "-"):
                tmp = mappedSpoint
                mappedSpoint = mappedEpoint
                mappedEpoint = tmp

        assert(mappedSpoint <= mappedEpoint)
        return (mappedSpoint, mappedEpoint)


class M4StreamReader:
    """Useful for parsing M4 streams as opposed to files"""
    def __init__(self, iterator):
        self.iterator = iterator
        self.delimiter = M4DELIMITER

    def __iter__( self ):
        for line in self.iterator:
            line = line.rstrip()
            if line[0] == '#' or len(line) == 0 or line == M4HEADER:
                continue

            yield M4Entry(line, self.delimiter)

    def setDelimiter(self, delimiter):
        """Set M4 delimiter."""
        self.delimiter = delimiter


class M4Reader(object):
    """M4 reader."""
    def __init__( self, fileName ):
        self.fileName = fileName
        if not os.path.exists( self.fileName ):
            sys.stderr.write( "Can't find file %s\n" % fileName )
            raise IOError, "M4Reader: can't find file %s" % fileName
        self.infile = open( self.fileName, 'r' )
        self.streamReader = M4StreamReader(self.infile)

    def __iter__(self):
        return self.streamReader.__iter__()

    def setDelimiter(self, delimiter):
        """Set M4 delimiter separating fields."""
        self.streamReader.setDelimiter(delimiter)

    def close(self):
        """Close the M4 file."""
        self.infile.close()

def by_tname(myBuffer):
    """Sort by target name."""
    return sorted(myBuffer, key=attrgetter('tname'))

def by_tmovie(myBuffer):
    """Sort by target movie."""
    return sorted(myBuffer, key=lambda anEntry:anEntry.tpbi.movie)

def by_tholeNumber(myBuffer):
    """Sort by target movie and then holeNumber."""
    newBuffer = sorted(myBuffer, key=lambda anEntry:anEntry.tpbi.holeNumber)
    return by_tmovie(newBuffer)

def by_abststart(myBuffer, isTargetPBIRead = False):
    """Sort by target movie, holeNumber and then abs_tstart."""
    newBuffer = sorted(myBuffer, key=attrgetter('abs_tstart'))
    # sorting by abs_tstart does not make sense unless alignments are
    # grouped by target name.
    if (not isTargetPBIRead):
        return by_tname(newBuffer)
    else:
        return by_tholeNumber(newBuffer)

def by_abststart_abstend(myBuffer):
    """Sort by target movie, holeNumnber, abs_tstart and then abs_tend."""
    # sorting by abs_tstart/end does not make sense unless alignments are
    # grouped by target name.
    newBuffer = sorted(myBuffer, key=attrgetter('abs_tend'), reverse=True)
    return sorted(newBuffer, key=attrgetter('abs_tstart'), reverse=False)

def by_qmovie(myBuffer):
    """Sort by query movie."""
    return sorted(myBuffer, key=lambda anEntry:anEntry.qpbi.movie)

def by_qholeNumber(myBuffer):
    """Sort by query movie and then hole number."""
    newBuffer = sorted(myBuffer, key=lambda anEntry:anEntry.qpbi.holeNumber)
    # sorting by holeNumber only makes sense when holeNumbers are grouped by
    # movie name.
    return by_qmovie(newBuffer)

def by_absqstart(myBuffer):
    """Sort by query movie, hole number, and abs_qstart."""
    newBuffer = sorted(myBuffer, key=attrgetter('abs_qstart'))
    # sorting by abs_qstart only makes sense when alignments are grouped
    # by hole numbers.
    return by_qholeNumber(newBuffer)

def by_score(myBuffer, isNegativeBetter=True):
    """Sort by alignment score."""
    return sorted(myBuffer, key=attrgetter('score'),
                  reverse=not isNegativeBetter)

