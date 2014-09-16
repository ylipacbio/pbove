"""PacBio SMRT Read title utils."""

class PBIReadName(object):
    """Parse a PacBio SMRT read title string in format:
       'movie/holeNumber'.
    """
    def __init__(self, readTitle):
        try:
            self.movie, self.holeNumber = readTitle.rstrip().split("/")[0:2]
            self.holeNumber = int(self.holeNumber)
        except ValueError as e:
            raise ValueError, "Could not parse PBI read name {0} as " \
                "movie/holeNumber.\n".format(readTitle) + str(e)

    def __str__(self):
        return "{0}/{1}".format(self.movie, self.holeNumber)

class PBISubreadName(PBIReadName):
    """Parse a PacBio subread title string in format:
       'movie/holeNumber/start1_end1'.
    """
    def __init__(self, subreadTitle):
        super(PBISubreadName, self).__init__(subreadTitle)
        try:
            start1_end1 = subreadTitle.rstrip().split("/")[2]
            self.start1, self.end1 = start1_end1.split("_")
            self.start1, self.end1 = int(self.start1), int(self.end1)
            self.start, self.end = self.start1, self.end1
        except ValueError as e:
            raise ValueError, "Could not parse PBI subread name {0} as " \
                "movie/holeNumber/start1_end1\n".format(subreadTitle) + str(e)

    def __str__(self):
        return "{0}/{1}/{2}_{3}".format(self.movie, self.holeNumber,
                                        self.start1, self.end1)

class PBISubsubreadName(PBISubreadName):
    """Parse a PacBio subsubread title string in format:
       'movie/holeNumber/start1_end1/start2_end2'.
       Note that the coordinate of start1 and end1 is relative to the
       full length PacBio reads, while the coordinate of start2 and end2
       is relative to the subread (within start1 and end1).
    """
    def __init__(self, subsubreadTitle):
        super(PBISubsubreadName, self).__init__(subsubreadTitle)
        try:
            fields = subsubreadTitle.rstrip().split("/")
            if (len(fields) == 4):
                start2_end2 = fields[3]
                self.start2, self.end2 = start2_end2.split("_")
                self.start2, self.end2 = int(self.start2), int(self.end2)
                self.start = self.start1 + self.start2
                self.end = self.start1 + self.end2
            else:
                raise ValueError
        except ValueError as e:
            raise ValueError, "Could not parse PBI subsubread name " \
                " {0} as movie/holeNumber/start1_end1/start2_end2.\n"\
                .format(subsubreadTitle) + str(e)
        if (self.start2 > self.end - self.start or
            self.end2 > self.end - self.start or
            self.start2 < 0 or self.end2 < 0):
            raise ValueError, "In PBI subsubread name {0}, " \
                " either start2 or end2 out of boundary." \
                .format(subsubreadTitle)

    def __str__(self):
        return "{0}/{1}/{2}_{3}/{4}_{5}".format(self.movie, self.holeNumber,
                                                self.start1, self.end1,
                                                self.start2, self.end2)


