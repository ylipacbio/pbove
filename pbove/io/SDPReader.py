"""Define SDPReader."""
from pbove.utils.PBIReadNameUtils import PBISubreadName

# An example of SDP file
#qid,tid,qstart,qend,qlen,tstart,tend,tlen,score
#m130812_185809_42141_c100533960310000001823079711101380_s1_p0/405/0_7884,m130812_185809_42141_c100533960310000001823079711101380_s1_p0/405/0_7884,0,7884,7884,0,7884,7884,-39420


class SDPRecord(object):
    """SDP Alignment."""
    def __init__(self, qID, tID, qStart, qEnd, qLength,
                 tStart, tEnd, tLength, score):
        self.qID = qID
        self.tID = tID
        self.qStart = int(qStart)
        self.qEnd = int(qEnd)
        self.qLength = int(qLength)
        self.tStart = int(tStart)
        self.tEnd = int(tEnd)
        self.tLength = int(tLength)
        self.score = float(score)
        try:
            self.qID = PBISubreadName(self.qID)
        except ValueError:
            pass

    def __eq__(self, another):
        return self.__dict__ == another.__dict__

    def __str__(self):
        msg = """
        qID: {qID}
        tID: {tID}
        qStart: {qStart}
        qEnd: {qEnd}
        qLength: {qLength}
        tStart: {tStart}
        tEnd: {tEnd}
        tLength: {tLength}
        """.format(qID=self.qID, tID=self.tID, qStart=self.qStart,
                   qEnd=self.qEnd, qLength=self.qLength,
                   tStart=self.tStart, tEnd=self.tEnd,
                   tLength=self.tLength, score=self.score)
        return msg

    @classmethod
    def fromString(cls, line, delimiter=','):
        """Interpret a string as a SDP record."""
        try:
            fields = line.rstrip().split(delimiter)
            assert(len(fields) == 9)
            qID = fields[0]
            tID = fields[1]
            qStart = fields[2]
            qEnd = fields[3]
            qLength = fields[4]
            tStart = fields[5]
            tEnd = fields[6]
            tLength = fields[7]
            score = fields[8]
            return SDPRecord(qID=qID, tID=tID, qStart=qStart,
                             qEnd=qEnd, qLength=qLength,
                             tStart=tStart, tEnd=tEnd,
                             tLength=tLength, score=score)
        except (AssertionError, ValueError):
            errMsg = "String not recognized as a valid SDP record."
            raise ValueError(errMsg)


class SDPReader(object):
    """SDP Reader."""
    def __init__(self, fileName):
        self.fileName = fileName
        try:
            self.infile = open(self.fileName, 'r')
        except IOError as e:
            errMsg = "SDPReader: could not read file " + \
                     fileName + "\n" + str(e)
            raise IOError(errMsg)

    def __iter__(self):
        try:
            for line in self.infile:
                line = line.strip()
                if len(line) == 0 or line[0] == "#" or \
                    line.endswith("score"):
                    continue
                yield SDPRecord.fromString(line=line)
        except ValueError as e:
            raise ValueError("SDPReader failed to parse " +
                             self.fileName + "\n" + str(e))

    def __enter__(self):
        return self

    def close(self):
        """Close the file."""
        self.infile.close()

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()

