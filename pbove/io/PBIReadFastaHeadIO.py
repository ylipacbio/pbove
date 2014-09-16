"""PBI reads fasta header IO."""
from pbove.utils.PBIReadNameUtils import PBISubreadName
from operator import attrgetter

class PBIReadFastaHeadReader:
    """PBI read Fasta head."""
    def __init__(self, fileName):
        self.reads = []
        for line in open(fileName):
            if line[0] == ">":
                self.reads.append (PBISubreadName(line[1:]))

    def sortby_movie(self):
        """Sort reads by movie name."""
        self.reads = sorted(self.reads, key=attrgetter('movie'))

    def sortby_holeNumber(self):
        """Sort reads by movie and hole number."""
        self.reads = sorted(self.reads, key=attrgetter('holeNumber'))
        self.sortby_movie()

    def sortby_start(self):
        """Sort reads by movie, hole number and start."""
        self.reads = sorted(self.reads, key=attrgetter('start'))
        self.sortby_holeNumber()

    def sortby_start_end(self):
        """Sort reads by movie, hole number, start and -end."""
        self.reads = sorted(self.reads, key=attrgetter('end'),
                            reverse=True)
        self.sortby_start()
