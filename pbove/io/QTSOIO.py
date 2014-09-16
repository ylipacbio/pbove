""" Define class QTSOIO which reads a QTSO file.
 Each line of a QTSO file has the following four fields:
  [0] qname - query read name
  [1] tname - target read name
  [2] score - score
  [3] overlap - overlap length
 Comments in a QTSO file starts with '#':
 '#numGroundTruthOVLPos=100" means that the number of
 ground truth overlapping query-target pairs is 100.
 '#numGroundTruthOVLNeg=200" means that the number of
 ground truth non-overlapping query-target pairs is 200.
"""


import os
import sys

QTSO_DELIMITER = "\t"
QTSO_FIELDS = ["qname", "tname", "score", "overlap"]
QTSO_HEADER = QTSO_DELIMITER.join(QTSO_FIELDS)


class QTSOEntry(object):
    """An QTSOEntry class object represents a record line of
    a QTSO file.
    """
    def __init__(self, line , delimiter="\t"):
        try:
            fields = line.rstrip().split(delimiter)
            self.qname, self.tname, self.score, self.overlap = \
            str(fields[0]), str(fields[1]), int(fields[2]), int(fields[3])
        except ValueError as e:
            raise ValueError("Could not recognize {l} as a QTSO entry.\n".
                             format(l=line)  + str(e))

    def __str__( self ):
        ret = "{0}, {1}, {2}, {3}".format(self.qname, self.tname,
              self.score, self.overlap)
        return ret


class QTSOReader(object):
    """A class for reading a QTSO file."""
    numP = 0
    numN = 0

    def __init__(self, fileName):
        self.fileName = fileName
        if not os.path.exists( self.fileName ):
            sys.stderr.write( "Can't find file %s\n" % fileName )
            raise IOError, "QTSOReader: can't find file %s" % fileName
        self.infile = open( self.fileName, 'r' )
        self.streamReader = QTSOStreamReader(self.infile)

    def __iter__(self):
        return self.streamReader.__iter__()

    def setDelimiter(self, delimiter):
        """Set QTSO delimiter separating fields."""
        self.streamReader.setDelimiter(delimiter)

    def close(self):
        """Close the QTSO file."""
        self.infile.close()


class QTSOStreamReader:
    """Useful for parsing QTSO streams as opposed to files"""
    def __init__(self, iterator):
        self.iterator = iterator
        self.delimiter = QTSO_DELIMITER

    def __iter__( self ):
        for line in self.iterator:
            line = line.rstrip()
            if line[0] == '#' or len(line) == 0 or line == QTSO_HEADER:
                continue

            yield QTSOEntry(line, self.delimiter)

    def setDelimiter(self, delimiter):
        """Set QTSO delimiter."""
        self.delimiter = delimiter

