""" Define class QTSO, which is short for Query_Target_Score_Overlap."""

from sys import maxint
from pbove.io.QTSOIO import QTSOReader


class QTSO(object):
    """ Define QTSO. """
    def __init__(self, fileName):
        self.fileName = fileName
        self.reader = QTSOReader(fileName)
        self.records = []
        for i in self.reader:
            self.records.append(i)
        self.sortByScore()
        self.deltaTable = []

    def sortByScore(self, reverse=False):
        """Sort QTSO entries by score."""
        self.records = sorted(self.records, key=lambda a:a.score,
                              reverse=reverse)

    def __str__(self):
        ret = "An QTSO object with {0} records.\n".format(len(self.records))
        for i in self.records[0:min(10, len(self.records))]:
            ret += str(i) + "\n"

        return ret

    def getMinMaxScores(self):
        """Return the minimum and maximum scores in all records."""
        minScore, maxScore = maxint, -maxint-1
        for i in self.records:
            if i.score < minScore:
                minScore = i.score
            if i.score > maxScore:
                maxScore = i.score
        return (minScore, maxScore)

    def getDeltaTable(self, stepSize, outfile="", overlapLengthCutoff=200):
        """ Generate a table each line of which has five fields:
        score_lower_bound (inclusive)
        score_upper_bound (exclusive),
        numdeltaTP:
            number of true positive read-read alignments whose blasr score is
            in [score_lower_bound, score_upper_bound) and overlap length is
            greater than or equal to overlapLengthCutoff.
        numdeltaFP:
            number of false positive read-read alignments whose blasr score is
            in [score_lower_bound, socre_upper_bound), and overlap length is
            less than overlapLengthCutoff.
        numdeltaPW:
            number of PW alignments whose score is in [score_lower_bound,
            score_upper_bound) and overlap length is in [1, overLaplengthCutoff).
        """
        (minScore, maxScore) = self.getMinMaxScores()
        overlapLengthCutoff = int(overlapLengthCutoff)
        indx = 0
        self.deltaTable = []
        s = (minScore/100) * 100 - 100

        headers = ("scoreLowerBound", "scoreUpperBound",
                  "numDeltaTruePositive", "numDeltaFalsePositive",
                  "numDeltaWeakPositive")

        if outfile != "" and outfile is not None:
            of = open (outfile, 'w')
            of.write("#" + "\t".join(headers) + "\n")

        while (s <= maxScore):
            # read-read overlap length l,
            # if l >=overlapLengthCutoff:
            #     True Positive
            #     predicted=positive, ground truth = positive
            # if l <= 0 :
            #     False Positive
            #     predicted=positive, ground truth = negative
            # otherwise:
            #     Positive Weak,
            #     predicted=Positive, ground truth = weak positive
            numdeltaTP, numdeltaFP, numdeltaPW = 0, 0, 0
            # Range [s, s+stepSize)
            while (indx < len(self.records)):
                r = self.records[indx]
                if r.score < s:
                    e = "ERROR computing delta table r.score = {0} < s = {1}".\
                        format(r.score, s)
                    raise Exception (e)
                if r.score >= s + stepSize:
                    break
                if r.score >= s and r.score < s + stepSize:
                    # print r
                    if (r.overlap >= overlapLengthCutoff):
                        numdeltaTP += 1
                    #    print " >= {0}, numdeltaTP = {1}".format(
                    #        overlapLengthCutoff, numdeltaTP)
                    elif (r.overlap <= 0):
                        numdeltaFP += 1
                    #    print " <= {0}, numdeltaFP = {1}".format(
                    #        overlapLengthCutoff, numdeltaFP)
                    elif (r.overlap < overlapLengthCutoff and
                        r.overlap > 0):
                        numdeltaPW += 1
                    #    print "0 < r.overlap < {0}, numdeltaPW = {1}".format(
                    #        overlapLengthCutoff, numdeltaPW)

                indx += 1
            if (numdeltaTP != 0 or numdeltaFP != 0 or numdeltaPW != 0):
                res = (s, s+stepSize, numdeltaTP, numdeltaFP, numdeltaPW)
                if outfile == "":
                    print "\t".join([str(item) for item in res])
                else:
                    of.write("\t".join([str(item) for item in res]) + "\n")
                self.deltaTable.append(res)
            s += stepSize

        if outfile != "":
            of.close()

    def getTable(self, numGTPos, numGTNeg, numGTWeak,
                 numUnmappableAlns, numMappableAlns, numAlns,
                 outfile=""):
        """
            Return a table each row of which has the following fields:
            "ScoreCutoff",
            "numTP", "numFP", "numFN", "numTN", "numPW", "numNW",
            "numGTPos", "numGTNeg", "numGTWeak",
            "numPredPos", "numPredNeg",
            "numUnmappableAlns", "numMappalbeAlns", "numAlns"

            'GT' means ground truth. The mapped locations of reads to reference genome
            using resequencing protocols (parameters) are considered as ground truth.

            numGTPos: number of ground truth positive (i.e., the number of <readi, readj>
            pairs that overlap by more than 200 bases according to ground truth)
            numGTNeg: number of ground truth negative (i.e., the number of <readi, readj>
            pairs that donot overlap at all according to ground truth)
            numGTWeak: number of <readi, readj> pairs that overlap by less than
            200 bases and more than 0 bases, according to the ground truth

            numPredPos: number of readi readj pairs whose blasr score is less than the
            score cutoff.
            numPredNeg: number of readi readj pairs which can not align to each other
            at all whose blasr score is greater than the score cutoff.

            numUnmappableAlns: number of <readi, readj> pairs of which either readi
            or readj can not map to the refernece genome.
            numMappableAln: number of <readi, readj> pairs of which both readi and readj
            can map to the reference genome.
            numAlns total number of alignments = |# of query reads| * |# of target reads|
        """
        numTP, numFP, numFN, numTN, numPW, numNW = 0, 0, 0, 0, 0, 0
        numGTPos  = long(numGTPos)
        numGTNeg  = long(numGTNeg)
        numGTWeak = long(numGTWeak)
        numPredPos, numPredNeg = 0, 0
        numUnmappableAlns, numMappableAlns, numAlns = long(numUnmappableAlns), \
        long(numMappableAlns), long(numAlns)
        header = ("ScoreCutoff",
            "numTP", "numFP", "numFN", "numTN", "numPW", "numNW",
            "numGTPos", "numGTNeg", "numGTWeak",
            "numPredPos", "numPredNeg",
            "numUnmappableAlns", "numMappableAlns", "numAlns")

        if outfile == "":
            print "\t".join(header)
        else:
            of = open(outfile, 'w')
            of.write("\t".join(header) + "\n")

        for deltaItem in self.deltaTable:
            # each deltaItem has five fields:
            # score lower bound, score upper bound, numdeltaTP, numdeltaFP,
            # numdeltaPW
            _lb, ub, numdeltaTP, numdeltaFP, numdeltaPW = deltaItem
            numTP += numdeltaTP
            numFP += numdeltaFP
            numPW += numdeltaPW
            numFN = numGTPos - numTP
            numTN = numGTNeg - numFP
            numNW = numGTWeak - numPW
            numPredPos = numTP + numFP + numPW
            numPredNeg = numFN + numTN + numNW
            assert(numPredPos + numPredNeg == numMappableAlns)
            assert(numGTPos + numGTNeg + numGTWeak == numMappableAlns)
            res = (ub,
                   numTP, numFP, numFN, numTN, numPW, numNW,
                   numGTPos, numGTNeg, numGTWeak,
                   numPredPos, numPredNeg,
                   numUnmappableAlns, numMappableAlns, numAlns)
            if outfile == "":
                print "\t".join([str(item) for item in res])
            else:
                of.write("\t".join([str(item) for item in res]) + "\n")

        if outfile != "":
            of.close()












