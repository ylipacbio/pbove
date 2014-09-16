#!/usr/bin/env python

"""Compare multiple aligners (e.g., different versions of blasr)
in the context of resequencing."""

from collections import defaultdict
from copy import deepcopy

class Interval(object):
    """Interval [start, end) contains integers [start, start+1, ..., end-1].
    Discrete, not continuous.
    """
    def __init__(self, start, end):
        self.start = int(start)
        self.end = int(end)
        if end < start:
            #(1) does not allow end < start,
            #(2) end == start means there is no element in the interval.
            raise AssertionError ("Value of interval [s, e) is invalid".
                                   format(s=start, e=end))

    def __eq__(self, right):
        return self.start == right.start and self.end == right.end

    def isempty(self):
        """Return True if this interval contains no element, otherwise
        False."""
        return (self.start == self.end)

    def intersect(self, right):
        """Return intersection of self and right."""
        if self.isintersect(right):
            return Interval(max(self.start, right.start),
                            min(self.end, right.end))
        else:
            return Interval(0, 0)

    def __str__(self):
        if self.isempty():
            return "[,)"
        return "[{start}, {end})".format(start=self.start, end=self.end)

    def isintersect(self, another):
        """Return whether or not this interval object intersect with
        another."""
        if self.isempty() or another.isempty():
            return False
        return not (self.start >= another.end or self.end <= another.start)

    def issubinterval(self, another):
        """Return whether or not this interval is a sub-interval of
        another."""
        if self.isempty():
            return True
        return (self.start >= another.start and self.end <= another.end)

    def issuperinterval(self, another):
        """Return whether or not this interval is a super-interval of
        another, which contains all integers in another."""
        if another.isempty():
            return True
        return (self.start <= another.start and self.end >= another.end)


class Intervals(object):
    """A list of intervals."""
    def __init__(self, in_intvs=None):
        self.intvs = []
        if in_intvs is not None:
            try:
                for (s, e) in in_intvs:
                    self.add(Interval(s, e))
            except (AssertionError, ValueError) as e:
                raise AssertionError (
                    "Failed to cpontruct an object of Intervals from a list.")

    def to_list(self):
        """Convert Intervals to a list of (start, end) tuples."""
        return [(intv.start, intv.end) for intv in self.intvs]

    def __str__(self):
        return ",".join([str(intv) for intv in self.intvs])

    def __eq__(self, right):
        if type(right) is Intervals:
            right = right.to_list()
        return self.to_list() == right

    @property
    def length(self):
        """Return total number of integers contained by this Intervals object."""
        return sum([intv.end-intv.start for intv in self.intvs])

    def add(self, new_intv):
        """Add a new interval to this object."""
        if type(new_intv) is not Interval:
            raise TypeError ("new_intv is not of type Interval but "
                             % type(new_intv))

        new_s, new_e = new_intv.start, new_intv.end
        if len(self.intvs) == 0:
            self.intvs = [new_intv]
        else:
            done = False
            for (idx, intv) in enumerate(self.intvs):
                end_updated = False
                if new_e < intv.start:
                    self.intvs.insert(idx, new_intv)
                    done = True
                elif new_s <= intv.start and new_e >= intv.start:
                    self.intvs[idx].start = new_s
                    if new_e >= intv.end:
                        intv.end = new_e
                        end_updated = True
                    done = True
                elif new_s >= intv.start and new_s <= intv.end:
                    if new_e > intv.end:
                        intv.end = new_e
                        end_updated = True
                    done = True
                elif new_s > intv.end:
                    pass
                else:
                    raise AssertionError ("Bad logic, please check " +
                        "[{new_s}, {new_e}) and [{s}, {e})".
                        format(new_s=new_s, new_e=new_e,
                        s=intv.start, e=intv.end))
                if done is True:
                    if end_updated is True:
                        stop = idx + 1
                        while stop < len(self.intvs):
                            if self.intvs[stop].end <= new_e:
                                del self.intvs[stop]
                            elif self.intvs[stop].end > new_e and \
                                 self.intvs[stop].start <= new_e:
                                intv.end = self.intvs[stop].end
                                del self.intvs[stop]
                                break
                            elif self.intvs[stop].start > new_e:
                                break
                            else:
                                raise AssertionError (
                                    "Bad logic, please check add")
                    break

            if done is not True:
                self.intvs.append(new_intv)

    def __add__(self, right):
        ret = deepcopy(self)
        if type(right) is list:
            right = Intervals(right)
        elif type(right) is not Intervals:
            raise AssertionError(
                "The right operand of 'Intervals +' should be of type " +
                "'Intervals' or list.")

        for intv in right.intvs:
            ret.add(intv)
        return ret

    def remove(self, rm_intv):
        """Remove an interval from this object."""
        if type(rm_intv) is not Interval:
            raise TypeError ("rm_intv is not of type " % Interval)

        if len(self.intvs) == 0 :
            pass
        else:
            rm_s, rm_e = rm_intv.start, rm_intv.end
            idx = 0
            while (idx < len(self.intvs)):
                intv = self.intvs[idx]
                if not intv.isintersect(rm_intv):
                    pass
                else:
                    if rm_intv.issuperinterval(intv):
                        del self.intvs[idx]
                        idx -= 1
                    elif rm_intv.issubinterval(intv):
                        if rm_s == intv.start:
                            intv.start = rm_e
                        elif rm_e == intv.end:
                            intv.end = rm_s
                        else:
                            self.intvs.insert(idx,
                                              Interval(intv.start, rm_s))
                            self.intvs[idx+1].start = rm_e
                            break
                    else: #intersect, not a super/sub-interval
                        if (rm_s > intv.start and rm_e >= intv.end):
                            intv.end = rm_intv.start
                        elif (rm_e > intv.start and rm_e <= intv.end):
                            intv.start = rm_intv.end
                        else:
                            raise AssertionError(
                                "Bad logic, please check remove")
                idx += 1

    def intersect(self, right):
        """Return intersection of this object and right."""
        if type(right) is list:
            right = Intervals(right)
        elif type(right) is not Intervals:
            raise AssertionError(
                "Intervals.intersect(x), where x should be of type " +
                "'Intervals' or list.")

        ret = Intervals()
        i, j = 0, 0
        while (i < len(self.intvs) and j < len(right.intvs)):
            intv_i = self.intvs[i]
            intv_j = right.intvs[j]
            if not intv_i.isintersect(intv_j):
                if intv_i.end <= intv_j.start:
                    i += 1
                elif intv_j.end <= intv_i.start:
                    j += 1
                else:
                    raise AssertionError("{i} not intersect with {j}".
                                         format(i=intv_i, j=intv_j))
            else: # i, j intersect
                ret.add(intv_i.intersect(intv_j))
                if intv_i.end < intv_j.end:
                    i += 1
                elif intv_j.end < intv_i.end:
                    j += 1
                else: # intv_i.end == intv_j.end
                    i += 1
                    j += 1
        return ret

    def __sub__(self, right):
        ret = deepcopy(self)
        if type(right) is list:
            right = Intervals(right)
        elif type(right) is not Intervals:
            raise AssertionError(
                "The right operand of 'Intervals -' should be of type " +
                "'Intervals' or list.")
        for intvs in right.intvs:
            ret.remove(intvs)
        return ret


class RefIntervals(object):
    """Intervals associated with a ref_id
    ref_id --> an object of class `Intervals` associated this ref id
    """
    def __init__(self):
        self.refintervals = defaultdict(Intervals)

    def __str__(self):
        ret = "RefIntervals: \n"
        for ref_id, intvs in self.refintervals.items():
            ret += "{ref_id} -> [".format(ref_id=ref_id)
            ret += "{intvs}]\n".format(intvs=intvs)
        ret += "\n"
        return ret

    def __len__(self):
        return len(self.refintervals)

    def __setitem__(self, key, item):
        if type(item) is list:
            item = Intervals(item)
        self.refintervals[key] = item

    def __getitem__(self, key):
        return self.refintervals[key]

    def __delitem__(self, key):
        del self.refintervals[key]

    def __iter__(self):
        return iter(self.refintervals)

    def items(self):
        return self.refintervals.items()

#
# a = Intervals([(1,100),(200,300), (500,900), (1, 10000)])
# b = Intervals([(1,10000), (200,300), (1,100), (500,900)])
# a == b
# a - [(500, 900)] == [(1, 500), (900, 10000)]
# a - [(500, 900)] - [(-1000, 300)] == [(300, 500), (900, 10000)]
# a - [(500, 900), (-1000, 300), (450, 1000)] == [(300, 450), (1000, 10000)]

# ris = RefIntervals()
# ris["m1"] = [(1, 100), (2000, 3000), (200, 300), (1100, 5500), (8000, 9000)]
# ris["m2"] = [(2000, 3000)]
# print ris

# its["m1"].add(Interval(0, 10000))
#Intervals:
#m1 -> [[0, 10000)]
#m2 -> [[2000, 3000)]

