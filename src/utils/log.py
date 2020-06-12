import numpy as np

class Logger(object):
    __instance = None
    level = 0
    wait = True

    def __new__(cls, level=0, wait=True):
        if Logger.__instance is None:
            Logger.__instance = object.__new__(cls)
            Logger.level = level
            Logger.wait = wait

        return Logger.__instance
    
    def out(self, text, verboseLevel=-1, wait=False):
        if self.level >= verboseLevel:
            print(text)
            if wait and self.wait:
                input()


class FileLogger(object):
    __instance = None
    lines = {}
    flushWhen = 255
    outdir = ""

    def __new__(cls, clean=False, outdir="."):
        if FileLogger.__instance is None or clean:
            FileLogger.__instance = object.__new__(cls)
            FileLogger.__instance.outdir = outdir
            FileLogger.__instance.lines = {}
        return FileLogger.__instance
    
    def flush(self, file):
        if file not in self.lines:
            return
        
        with open(self.outdir + "/" + file, "a") as writer:
            for line in self.lines[file]:
                writer.write(line + '\n')
            self.lines[file] = []

    def flush_all(self):
        for file in self.lines:
            if len(self.lines[file]) > 0:
                self.flush(file)

    def write(self, file, line):
        if file not in self.lines:
            self.lines[file] = [line]
        else:
           self.lines[file].append(line)
           if len(self.lines[file]) >= self.flushWhen:
               self.flush(file)
           
    def write_cols(self, file, cols, sep='\t'):
        self.write(file, sep.join([str(c) for c in cols]))
        
def rbg_generator(string): 
    if string is None: return [255, 255, 255]
    np.random.seed((hash(string) % 2**30))
    rands = [np.random.rand(), np.random.rand(), np.random.rand()]
    rounded = [round(1e7*r) % 256 for r in rands]
    rgb = ",".join([str(c) for c in rounded])
    return rgb
