"""Define class PreassmblyPrediction."""
from pbove.io.M4IO import M4Reader
from datetime import datetime

class PreassemblyPrediction(object):
    """ Read an M4 file which contains read-read alignments as
        preassembly predictions into memory."""
    def __init__(self, fileName):
        reader = M4Reader(fileName)
        self.readToRead = []
        for i in reader:
            self.readToRead.append(i)

    def __str__(self):
        return "{0} read to read prediction.\n" \
               .format(len(self.readToRead))

    def QInReference(self, groundTruth, index, infer=False):
        """Map query read of self.readToRead[index] to reference."""
        i = self.readToRead[index]
        return groundTruth.MapPBISubreadToReference(
                i.qpbi.movie, i.qpbi.holeNumber,
                i.abs_qstart, i.abs_qend, infer)

    def TInReference(self, groundTruth, index, infer=False):
        """Map target read of self.readToRead[index] to reference."""
        i = self.readToRead[index]
        return groundTruth.MapPBISubreadToReference(
                i.tpbi.movie, i.tpbi.holeNumber,
                i.abs_tstart, i.abs_tend, infer)

    def OverlapLengthOfQandTInReference(self, groundTruth, index, infer=False):
        """Map a read-read M4 entry, compute the overlap length between
           query and Target."""
        return groundTruth.OverlapLengthOfQandTInReference(
                self.readToRead[index], infer)


    def OverlapLengthsInReference(self, groundTruth, infer=False):
        """For each read to read alignment in the M4 file, compute the
        overlap length of their mapped intervals in the reference genome. """
        for index, i in enumerate(self.readToRead):
            if index % 1000 == 0:
                print "Processing {0} / {1} read-read alignment at {2}" \
                      .format(index, len(self.readToRead), str(datetime.now()))
            i.QMappable = groundTruth.IsMappable(i.qpbi)
            i.TMappable = groundTruth.IsMappable(i.tpbi)
            i.overlapLength = groundTruth.OverlapLengthOfQandTInReference(
                i, infer)

    def ToQTSO(self, outfile=""):
        """For each read to read alignment, print in QTSO format."""
        if outfile != "":
            of = open (outfile, 'w')
        for i in self.readToRead:
            if not i.QMappable or not i.TMappable:
                continue
            record = "{0}\t{1}\t{2}\t{3}\n"\
                    .format(i.qpbi, i.tpbi, i.score, i.overlapLength)
            if outfile == "":
                print record
            else:
                of.write(record)

        if outfile != "":
            of.close()


    def ComputeTPFP(self, groundTruth, overlapLengthCutoff=200, infer=False):
        """Compute the number of true positive and false positive overlapping
        read-read alignment."""
        numTP, numFP, numBoundary = 0, 0, 0
        # TP: true positive. An alignment R is true positive iff map(R.query,
        # reference) and map(R.target, reference) overlap by more than 200 bps.
        # FP: false positive. An alignment R is false positive iff map(R.query,
        # reference) and map(R.target, reference) overlap by 0 bps.
        # Boundry) An alignment R is boundary iff map(R.query, reference)
        # and map(R.target, reference) overlap by > 0 and < 200 bps
        for index, i in enumerate(self.readToRead):
            if index % 1000 == 0:
                print "Processing {0} / {1} read-read alignment at {2}" \
                      .format(index, len(self.readToRead), str(datetime.now()))
            i.overlapLength = groundTruth.OverlapLengthOfQandTInReference(
                i, infer)
            if i.overlapLength >= overlapLengthCutoff:
                numTP += 1
            elif i.overlapLength <= 0:
                numFP += 1
            else:
                numBoundary += 1
        return (numTP, numFP, numBoundary)

