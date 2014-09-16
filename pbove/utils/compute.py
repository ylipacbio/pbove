"""Provide util functions for computing mappable/unmapble reads,
number of good/bad overlaps."""
from sys import maxint

def GetOverlapLengthOfTwoIntervals(start1, end1, start2, end2):
    """Given two intervals [start1, end1) and [start2, end2),
    if they overlap, return length of the overlapping region.
    Otherwise, return (0 - distance), where distance between
    two intervals equals the minimum of distances between any
    two points each in one interval.
    """
    ovlLen = min(end1, end2) - max(start1, start2)
    if ovlLen >= 0:
        return ovlLen
    else:
        return - min(abs(min(start1, end1 - 1) - max(start2, end2 - 1) - 1),
                     abs(max(start1, end1 - 1) - min(start2, end2 - 1) - 1))

def GetDistanceBetweenTwoIntervals(start1, end1, start2, end2):
    """Given two intervals [start1, end1) and [start2, end2),
    if they overlap, return 0 - length of the overlapping region.
    Otherwise, return the minimum distance between any two
    points each in one interval."""
    return - GetOverlapLengthOfTwoIntervals(
             start1, end1, start2, end2)

def GetBestReadToQueryForInterval(reads, start, end):
    """Given an array of reads and an interval [start, end),
    return the read whose query interval [abs_qstart, abs_qend)
    has the longest overlapping regions with the interval;
    otherwise, return the read whose query interval has the
    shortest distance to the interval."""
    bestIndex = -1
    minDist = maxint
    for index, read in enumerate(reads):
        dist = GetDistanceBetweenTwoIntervals(
                read.abs_qstart, read.abs_qend, start, end)
        if dist < minDist:
            minDist = dist
            bestIndex = index
    if bestIndex != -1:
        return reads[bestIndex]
    return None

def GetNumOfUnmappable(q):
    """Given a sorted list of [(readobj, reference, refStart, refEnd), ...],
    return the number of unmappable read objects in q, where q is
    sorted by reference, refStart, and then -refEnd.
    """
    numUnmappable = 0
    for index, qitem in enumerate(q):
        _qread, qref, qrefstart, qrefend = qitem
        if qref == "":
            numUnmappable += 1
        if index == 0:
            continue
        else:
            _qqread, qqref, qqrefstart, qqrefend = q[index-1]
            assert(qqref < qref or
                   (qqref == qref and qqrefstart < qrefstart) or
                   (qqref == qref and qqrefstart == qrefstart and
                    qqrefend >= qrefend)
                  )
    return numUnmappable


def ComputeAllPosNegNumbers(q, t, overlapLengthCutoff=200):
    """Given a list of sorted query reads q and a list of sorted
    target reads t, for each pair of (q_i, t_j) where q_i in q,
    and t_j in t, identify whether they overlap or not. Return
    numGTPos:  number of alns that indeed overlap by more than
               overlapLengthCutoff base pairs.
    numGTNeg:  number of alns that do not overlap (negative),
    numGTWeak: number of alns that overlap by a short reigon
               with less than overlaplengthCutoff base pairs
    numUnmappableAlns: number of alns that either query or target
               read is not mappable to the reference.
    numMappableAlns: number of alns that both query and target reads
               can map to the reference.
    numAlns: |query reads| * |target reads|

    q: query, each item has four columns, sorted by reference, then
    refStart, then -refEnd:
       (query_read, reference, refStart, refEnd).
    t: target, each item has four columns, sorted by reference, then
    refStart, then -refEnd:
       (target_read, reference, refStart, refEnd).
    overlapLengthCutoff: the minumum number of overlapping bases to be
    considered positive overlap.
    """
    numQMappable = len(q) - GetNumOfUnmappable(q)
    numTMappable = len(t) - GetNumOfUnmappable(t)
    numQTMappable = numQMappable * numTMappable
    numQTUnmappable  = len(q) * len(t) - numQTMappable
    numAlns = len(q) * len(t)

    # Get the maximum length of intervals in the reference which the target
    # reads can map to.
    (_objcol, _refcol, startcol, endcol) = (0, 1, 2, 3)
    maxtlen = 0
    for titem in t:
        maxtlen = max(maxtlen, titem[endcol] - titem[startcol])
    print "maxtlen = {0}".format(maxtlen)

    # if q & t overlap length >= OverlapLengthCutoff,
    #    it is a ground truth positive overlap
    # if q & t overlap length > 0 && < OverlapLengthCutoff,
    #    it is a gt weak overlap
    # otherwise, negative overlap, which is computed as
    #     numGTNeg = numQTMappable - numGTPos - numGTWeak
    numGTPos, numGTWeak = 0, 0
    for qindex, qitem in enumerate(q):
        # for each query read, first identify a range of target
        # reads which the query read may overlap with.
        _qread, qref, qrefstart, qrefend = qitem
        if (qindex % 100 == 0):
            print "Processing query {0} / {1} ".format(
                    qindex, len(q))
        if (qref == ""):
            # this query read is not mappable to the reference,
            continue

        searchStart, searchEnd = 0, len(t) - 1
        searchMid = -1
        while (searchStart <= searchEnd and
               searchStart >= 0 and searchEnd >= 0 and
               searchStart < len(t) and searchEnd < len(t)):
            searchMid = int((searchStart + searchEnd) / 2)
            _midread, midref, midrefstart, _midrefend = t[searchMid]
            if (qref < midref or
                (qref == midref and qrefend < midrefstart)):
                searchEnd = searchMid - 1
                continue
            elif (qref > midref):
                searchStart = searchMid + 1
                continue
            elif (qref == midref and qrefstart > midrefstart + maxtlen):
                searchStart = searchMid + 1
                continue
            else:
                break

        if (searchStart > searchEnd):
            # This query read is mappable to the reference, but
            # do not overlap with any target reads.
            pass
        else:
            # This query read may overlap with some target
            # reads in the search range.
            for titem in t[searchStart: searchEnd + 1]:
                _tread, tref, trefstart, trefend = titem
                if (tref == "" or qref != tref ):
                    continue
                elif (trefend < qrefstart or qrefend < trefstart):
                    continue

                overlapLength = int(GetOverlapLengthOfTwoIntervals(
                        qrefstart, qrefend,
                        trefstart, trefend))
                if int(overlapLength) >= int(overlapLengthCutoff):
                    numGTPos += 1
                elif overlapLength > 0:
                    numGTWeak += 1

    numGTNeg = numQTMappable - numGTPos - numGTWeak
    print ("numGTPos={0}, numGTNeg = {1}, numGTWeak={2}, ".\
            format(numGTPos, numGTNeg, numGTWeak))

    print ("numUnmappable={0}, numQTMappable = {1}, numAlns ={2}, ".\
            format(numQTUnmappable, numQTMappable, numAlns))
    return (numGTPos, numGTNeg, numGTWeak, numQTUnmappable,
            numQTMappable, numAlns)


def write_gt_overlaps(query, target, out_file):
    """Given a list of sorted query reads q and a list of sorted
    target reads t, for each pair of (q_i, t_j) where q_i in q,
    and t_j in t, identify whether they (ground truth) overlap or
    not and save tuples (qread, tread, overlap_len) to out_file.

    query      : query, each item has four columns, sorted by reference, then
                 refStart, then -refEnd:
                 (query_read, reference, refStart, refEnd).
    target     : target, each item has four columns, sorted by reference, then
                 refStart, then -refEnd:
                 (target_read, reference, refStart, refEnd).
    #ovl_cut_off: the minumum number of overlapping bases to be
    #             considered positive overlap.
    """
    # if q & t overlap length >= OverlapLengthCutoff,
    #    it is a ground truth positive overlap
    # if q & t overlap length > 0 && < OverlapLengthCutoff,
    #    it is a gt weak overlap
    # Here, we report all gt overlaps including weak overlaps
    of = open(out_file, 'w')
    of.write("#query\ttarget\toverlap_len\n")

    for qindex, qitem in enumerate(query):
        qread, qref, qrefstart, qrefend = qitem
        #if (qindex % 1000 == 0):
        #    print "Processing query {0} / {1} ".format(qindex, len(query))
        if (qref == ""):
            # this query read is not mappable to the reference,
            continue
        for _tindex, titem in enumerate(target):
            tread, tref, trefstart, trefend = titem
            if (qref == tref):
                ol = max(0, min(trefend, qrefend) - max(trefstart, qrefstart))
                if ol > 0:
                    of.write("{qread}\t{tread}\t{ovl_len}\n".format(
                             qread=qread, tread=tread, ovl_len=ol))
    of.close()
