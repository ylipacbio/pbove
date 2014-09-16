"""Define class ReseqGroundTruth."""
from pbove.io.M4IO import M4Reader, by_absqstart
from pbove.io.PBIReadFastaHeadIO import PBIReadFastaHeadReader
from pbove.utils.compute import *
from operator import itemgetter
from sys import maxint

class ReseqGroundTruth(object):
    """Resequencing ground truth."""
    def __init__(self, fileName):
        reader = M4Reader(fileName)
        self.readToReference = []
        for i in reader:
            self.readToReference.append(i)
        # Sort readToReference by query qstart
        self.readToReference = by_absqstart(self.readToReference)

    def __str__(self):
        return "{0} reads in ground truth.\n" \
               .format(len(self.readToReference))

    def searchRead(self, movie, holeNumber):
        """Search reads which have the specified movie and
        hole number."""
        searchStart, searchEnd = 0, len(self.readToReference) - 1
        while(searchStart <= searchEnd):
            searchMid = int ((searchStart + searchEnd) / 2)
            midEntry = self.readToReference[searchMid]
            if (midEntry.qpbi.movie < movie or
                (midEntry.qpbi.movie == movie and
                 midEntry.qpbi.holeNumber < holeNumber)):
                searchStart = searchMid + 1
            elif (midEntry.qpbi.movie > movie or
                  (midEntry.qpbi.movie == movie and
                   midEntry.qpbi.holeNumber > holeNumber)):
                searchEnd = searchMid - 1
            else:
                break

        subBuffer = []
        if (searchStart <= searchEnd):
            searchMid = int ((searchStart + searchEnd) / 2)
            while (searchMid >= 0 and searchMid < len(self.readToReference)):
                midEntry = self.readToReference[searchMid]
                if (midEntry.qpbi.movie == movie and
                    midEntry.qpbi.holeNumber ==  holeNumber):
                    subBuffer.append(midEntry)
                    searchMid -= 1
                else:
                    break

            searchMid = int (searchStart + searchEnd) / 2 + 1
            while (searchMid >= 0 and searchMid < len(self.readToReference)):
                midEntry = self.readToReference[searchMid]
                if (midEntry.qpbi.movie == movie and
                    midEntry.qpbi.holeNumber ==  holeNumber):
                    subBuffer.append(midEntry)
                    searchMid += 1
                else:
                    break
        return subBuffer

    def MapPBISubreadToReference(self, movie, holeNumber, start, end,
            infer=False):
        """Map a PBI subread (given its movie, holeNumber, start and end) to
        a reference and return (reference_name, mapped_start, mapped_end).

        If no reads which have the specified movie and holeNumber map to
        the reference, return ("", -1, -1).
        If one or more than one read, has the specified movie and holeNumber,
        and can map to the reference, we denote one of them as the 'best',
        which should overlap with the subread in the unrolled zmw sequence
        by more bases than the others.
        If infer is True, we infer where the subread [start, end) has
        mapped to based on where the 'best' read has mapped to.
        If infer is False, return ("", -1, -1) if either endpoint of the
        input subread is not included by the best read.

        Note that coordinate of start and end is relative to the full length
        unrolled read.
        """
        # Get reads which have the specified movie and holeNumber.
        subBuffer = self.searchRead(movie, holeNumber)
        # Get the best read to query which overlap with interval
        # [start, end) by most bases.
        bestRead = GetBestReadToQueryForInterval(subBuffer, start, end)
        if bestRead is not None:
            mappedStart, mappedEnd = bestRead.MapAnIntervalFromQToT(
                start, end, infer)
            assert(mappedStart <= mappedEnd)
            if (mappedStart != -1 and mappedEnd != -1):
                return (bestRead.tname, mappedStart, mappedEnd)
            # Could not map this PBI subread to anywhere in the reference
        return ("", -1, -1)

    def IsMappable(self, read, infer=False):
        """Can this read be mapped to the reference?"""
        ref, _s, _e = self.MapPBISubreadToReference(read.movie,
                read.holeNumber, read.start, read.end, infer)
        if ref == "":
            return False
        return True

    def OverlapLengthOfQandTInReference(self, rr, infer=False):
        """Given a read-read M4 entry, compute the overlap length between
           its query and target."""
        qmovie, qholeNumber, qstart, qend = (
            rr.qpbi.movie, rr.qpbi.holeNumber,
            rr.abs_qstart, rr.abs_qend)
        tmovie, tholeNumber, tstart, tend  = (
            rr.tpbi.movie, rr.tpbi.holeNumber,
            rr.abs_tstart, rr.abs_tend)
        return self.OverlapLengthInReference(
                qmovie, qholeNumber, qstart, qend,
                tmovie, tholeNumber, tstart, tend, infer)

    def OverlapLengthInReference(self,
            movie1, holeNumber1, start1, end1,
            movie2, holeNumber2, start2, end2,
            infer=False):
        """Map two PBI subreads (given movie, holeNumber, start, end)
        to the reference and return overlap length of intervals
        where the two reads are mapped to. If the two PBI subreads
        can not map to the same reference, return (-sys.maxint).
        If two PBI subreads can map to the same reference but are
        apart, return (-distance) between them."""
        refName1, refStart1, refEnd1 = self.MapPBISubreadToReference(
            movie1, holeNumber1, start1, end1, infer)
        if (refName1 == ""):
            return -maxint

        refName2, refStart2, refEnd2 = self.MapPBISubreadToReference(
           movie2, holeNumber2, start2, end2, infer)

        if refName1 != refName2 or refName2 == "":
            return -maxint
        else:
            return GetOverlapLengthOfTwoIntervals(refStart1, refEnd1,
                                                  refStart2, refEnd2)

    def MapPBISubreadsToReference(self, reads, infer=False):
        """Map a list of PBI Subreads to references and produce a list of
        [(readobj, ref, start, end), ..., ()]. Sort this list by reference
        name, start, and end and return the sorted list.
        """
        retList = []
        (refcol, startcol, endcol) = (1, 2, 3)
        for read in reads:
            ref, refstart, refend = self.MapPBISubreadToReference(
                read.movie, read.holeNumber, read.start, read.end, infer)
            retList.append((read, ref, refstart, refend))
        retList = sorted(retList, key=itemgetter(endcol), reverse=True)
        retList = sorted(retList, key=itemgetter(startcol))
        retList = sorted(retList, key=itemgetter(refcol))
        return retList


    def ComputeAllPosNegNumbersFromFiles(self, queryFile, targetFile,
            overlapLengthCutoff=200):
        """Given a query read file, and a target reads file,
        compute the number of ground truth positive (a query read and
        a target read overlap by at least overlapLengthCutoff bases)
        and ground truth negative (a query read and a target read
        overlap by less than overlapLengthCutoff bases).
        Return (num_positive, num_negative)."""
        queryReads = PBIReadFastaHeadReader(queryFile)
        print "There are {0} query reads in total." \
              .format(len(queryReads.reads))

        targetReads = PBIReadFastaHeadReader(targetFile)
        print "There are {0} target reads in total." \
              .format(len(targetReads.reads))

        q = self.MapPBISubreadsToReference(queryReads.reads)
        t = self.MapPBISubreadsToReference(targetReads.reads)

        return ComputeAllPosNegNumbers(q, t, overlapLengthCutoff)



