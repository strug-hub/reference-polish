class SimpleRegion:
    def __init__(self, chrom, start, end, length=None):     
        self.chrom = chrom
        self.start = start
        self.end = end
        self.length = length
            
    def contains(self, chrom, pos):
        if chrom != self.chrom: return False
        if pos < self.start: return False
        if pos > self.end: return False
        return True
    
    def contains_pos(self, pos):
        if pos < self.start: return False
        if pos > self.end: return False
        return True

    def overlaps(self, otherRegion):
        if self.chrom != otherRegion.chrom: return False
        return max(self.start, otherRegion.start) <= min(self.end, otherRegion.end)
        
    def str_base1(self):
        self.start += 1 ; self.end += 1
        string = str(self)
        self.start -= 1 ; self.end -= 1
        return string

    def extend_left(self, size):
        newStart = max(0, self.start - size)
        return SimpleRegion(self.chrom, newStart, self.end)
        
    def extend_right(self, size):
        if self.length is None:
            newEnd = self.end + size
        else:
            newEnd = min(self.length-1, self.end + size)
        return SimpleRegion(self.chrom, self.start, newEnd)

    def extend(self, size):
        newStart = max(0, self.start - size)
        if self.length is None:
            newEnd = self.end + size
        else:
            newEnd = min(self.length-1, self.end + size)
        return SimpleRegion(self.chrom, newStart, newEnd)
        
    def __sub__(self, other):
        if other.chrom != self.chrom: return [self]
     
        if self.start < other.start:
            if self.end < other.start:
                return [self]
            elif self.end > other.end:
                sr1 = SimpleRegion(self.chrom, self.start, other.start-1)
                sr2 = SimpleRegion(self.chrom, other.end+1, self.end)
                return [sr1, sr2]
            elif self.end >= other.start:
                return [SimpleRegion(self.chrom, self.start, other.start-1)]
            
        if self.start >= other.start:
            if other.end < self.start:
                return [self]
            if self.end <= other.end:
                return []
            if other.end >= self.start:
                return [SimpleRegion(self.chrom, other.end+1, self.end)]
    
    def __and__(self, other):
        if other.chrom != self.chrom: return []
        
        if (other.start > self.end) or (self.start > other.end):
            return []
        else:
            s = max(self.start, other.start)
            e = min(self.end, other.end)
            return [SimpleRegion(self.chrom, s, e)]

    def __eq__(self, other):
        return self.chrom == other.chrom and  \
               self.start == other.start and self.end == other.end
        
    def __len__(self):
        return self.end - self.start + 1
    def __repr__(self):
        return self.chrom + ":" + str(self.start) + "-" + str(self.end) + \
                " (" + str(len(self)) + "bp)"
            
    def clean_str(self):
        return self.chrom + "_" + str(self.start) + "_" + str(self.end)
    def __str__(self):
        return self.chrom + ":" + str(self.start) + "-" + str(self.end)

def region_from_string(string, lengthData=None):
    chrom = string.split(":")[0]
    start = int(string.split(":")[1].split("-")[0])
    end = int(string.split(":")[1].split("-")[1])
    return SimpleRegion(chrom, start, end, lengthData)     

def region_difference(region, otherRegions, minSize=2):
    
    if minSize is not None and len(region) < minSize: return []

    split = [region]
    for splitter in otherRegions:
        newSplit = []
        for s in split: newSplit.extend(s - splitter)
        split = newSplit
        if minSize is not None:
            split = [s for s in split if len(s) >= minSize]

    return split


def region_overlap(region, otherRegions, minOverlap=0):

    ovelaps = []
    for otherRegion in otherRegions:
        intersection = region & otherRegion
        
        if len(intersection) > 0:
            if len(intersection[0]) > minOverlap:
                ovelaps.append(otherRegion)
            
    return ovelaps