import re
import log
import numpy as np
import matplotlib.pyplot as plt
import random
import matplotlib.patheffects as pe
#import path_helper as ph

def plot_identity(contig, outputPath=None):
    
    if len(contig.mblocks) < 1:
        return
    
    qid = contig.mblocks[0].qid
    identities=[]
    avg=[]

    blockidx=[]
    fig = plt.figure()

    for mblock in contig.mblocks:
        c=int(255*random.random())
        for block in mblock:
            ids = block.percent_identities()
            identities.extend(ids)
            avg.extend([np.mean(ids)] * len(ids))
            blockidx.extend([c] * len(ids))

    y = np.array(identities)
    g = np.array(blockidx)
    n=len(identities)
    x = np.array(range(n))
    
    minval=94
    #wsize=n/500
    #w=np.array([ np.mean(y[max(0,i-int(wsize/2)):min(n-1, i+int(wsize/2))]) for i in range(n)])
    #w=100-w+minval
    w=np.array(avg)
    plt.plot(x,w, linewidth=1, alpha=0.75, color="black", \
             path_effects=[pe.Stroke(linewidth=5, foreground='white', alpha=0.7), pe.Normal()])

    plt.scatter(x, y, c=g, alpha=0.4, s=1)
    plt.ylim((minval,100))
    
    if outputPath is not None:
        fig.savefig(outputPath + str(qid) + ".png")
    
    plt.close('all')
 
def analyze_paths(paths, lengthData, param):

    lens = []
    ridCount = []
    
    for path in paths:

        rids = set()
        for fork in path:
            rids.add(fork.qid)
        ridCount.append(len(rids))
        lens.append(ph.path_length(path))
    
    ridCount=[x for _,x in sorted(zip(lens,ridCount))]
    lens=[x for x in sorted(lens)]

    plt.scatter(range(len(lens)), np.log10(lens), s=np.array(ridCount)*np.array(ridCount), alpha=0.5, c=np.array(ridCount)*100)



def analyze_scaffolds(scaffolds, leftovers, lengthData, param):

    def normalize(fork, q=True):
        pos = fork.qpos if q else fork.rpos
        strand = (fork.qstrand if q else fork.rstrand)
        tigId = str(fork.rid if q else fork.qid)
        
        if strand == -1: 
            return lengthData[tigId] - pos
        return pos

    rtails = 0
    tigs = dict()
    rtigs = set()
    rtigsEnds = set()

    for scaffold in scaffolds + leftovers:

    
        rtigsEnds.add(scaffold[0].rid)
        rtigsEnds.add(scaffold[-1].rid)
        printed=False
        for fork in scaffold:
            rtigs.add(fork.rid)
            
            if not printed and fork.rid == "tig00983269_pilon_pilon":
                print(scaffold.__repr__())
                printed=True
    
    rtigsDiff = rtigs.difference(rtigsEnds)
        
    
'''
 
   for scaffold in scaffolds:

        firstFork = scaffold[0]
        rlen = lengthData[str(firstFork.rid)]
        #rstrand =  fork.rstrand
        print( firstFork )

        rtail = firstFork.rpos
        rtails = rtails + rtail
        print(rtail)
        if firstFork.rid in tigs: tigs[firstFork.rid] = tigs[firstFork.rid] +1
        else: tigs[firstFork.rid] = 1
        
        
        
        normalize(firstFork, q=False)
        
        id = fork.before_id()
        pos = fork.before_pos()
        seq = seqData[str(id)]
        if fork.before_strand() == -1:
            pos = len(seq) - pos

        forkSeq = seq[pos-nbases:pos+nbases]

        if fork.before_strand() == -1:
            forkSeq = reverse_complement(forkSeq)
        
        print(forkSeq[:nbases] + nbases*"-" + source)
        print(forkSeq + source)
        
        #---------------------
        
        source = " CANU" if fork.switch == 'r' else " NOVA"

        id = fork.after_id()
        pos = fork.after_pos()
        seq = seqData[str(id)]
        if fork.after_strand() == -1:
            pos = len(seq) - pos

        forkSeq = seq[pos-nbases:pos+nbases]

        if fork.after_strand() == -1:
            forkSeq = reverse_complement(forkSeq)
        
        print(forkSeq + source)       
        print(nbases*"-" + forkSeq[nbases:] + source)

def sequence_validation(sequences):
    f = open("validate.fasta", "w+")
    i = 0
    for sequence in sequences:
        
        if len(sequence) < 3001:
            f.write(">segment" + "_" + str(i) + "\n")
            f.write(re.sub("(.{64})", "\\1\n", sequence, 0, re.DOTALL) + "\n")
        else:
            f.write(">segment" + "_" + str(i) + "start\n")
            f.write(re.sub("(.{64})", "\\1\n", sequence[:1000], 0, re.DOTALL) + "\n")
            f.write(">segment" + "_" + str(i) + "end\n")
            f.write(re.sub("(.{64})", "\\1\n", sequence[-1000:], 0, re.DOTALL) + "\n")

        i = i + 1
    f.close()      
'''

def plot_scaffold_hist(scaffoldList):
    if len(scaffoldList) > 3:
        scaffoldList = [scaffoldList]
    alpha=1
    for scaffolds in scaffoldList:
        x=np.array([ph.path_length(path) for path in scaffolds])
        plt.hist(np.log10(x), bins=50, alpha=alpha)
        alpha=0.5
    plt.show()
    plt.close("all")
    
def plot_scaffold_line(scaffoldList):
    if len(scaffoldList) > 3:
        scaffoldList = [scaffoldList]
    alpha=1
    for scaffolds in scaffoldList:
        y = np.sort(np.array([ph.path_length(path) for path in scaffolds]))[::-1]
        x = np.array(range(len(y)))
        z = [np.sum(y[:n]) for n in x]
        
        plt.plot(np.log10(x), z, alpha=alpha)
        alpha=0.5
    plt.show()
    plt.close("all") 
        
def plot_scaffold_line2(scaffoldList):
    if len(scaffoldList) > 3:
        scaffoldList = [scaffoldList]
    alpha=1
    for scaffolds in scaffoldList:
        y = np.sort(np.array([ph.path_length(path) for path in scaffolds]))[::-1]
        x = np.array(range(len(y)))
        z = [np.sum(y[:n]) for n in x]
        
        plt.plot(x, z, alpha=alpha)
        alpha=0.5
    plt.show()
    plt.close("all") 

def validate_forks(path, seqData, nbases=25):

    for fork in path:
        
        print( "=========================" )
        
        source = " NOVA" if fork.switch == 'r' else " CANU"

        id = fork.before_id()
        pos = fork.before_pos()
        seq = seqData[str(id)]
        if fork.before_strand() == -1:
            pos = len(seq) - pos

        forkSeq = seq[pos-nbases:pos+nbases]

        if fork.before_strand() == -1:
            forkSeq = reverse_complement(forkSeq)
        
        print(forkSeq[:nbases] + nbases*"-" + source)
        print(forkSeq + source)
        
        #---------------------
        
        source = " CANU" if fork.switch == 'r' else " NOVA"

        id = fork.after_id()
        pos = fork.after_pos()
        seq = seqData[str(id)]
        if fork.after_strand() == -1:
            pos = len(seq) - pos

        forkSeq = seq[pos-nbases:pos+nbases]

        if fork.after_strand() == -1:
            forkSeq = reverse_complement(forkSeq)
        
        print(forkSeq + source)       
        print(nbases*"-" + forkSeq[nbases:] + source)

'''
def sequence_validation(sequences):
    f = open("validate.fasta", "w+")
    i = 0
    for sequence in sequences:
        
        if len(sequence) < 3001:
            f.write(">segment" + "_" + str(i) + "\n")
            f.write(re.sub("(.{64})", "\\1\n", sequence, 0, re.DOTALL) + "\n")
        else:
            f.write(">segment" + "_" + str(i) + "start\n")
            f.write(re.sub("(.{64})", "\\1\n", sequence[:1000], 0, re.DOTALL) + "\n")
            f.write(">segment" + "_" + str(i) + "end\n")
            f.write(re.sub("(.{64})", "\\1\n", sequence[-1000:], 0, re.DOTALL) + "\n")

        i = i + 1
    f.close()      
'''

def get_edges(sequence, nbases=1000, toPrint=False):
    edges = []
    for i in range(len(sequence)-1):
        edge = sequence[i][-int(nbases/2):] +\
            sequence[i+1][:int(nbases/2)]
        edges.append(edge)
        if toPrint:
            print(">" + str(i))
            print(edge)
        
    return edges


def path_to_sequence(path, seqData, file=None, invert=False):
    sequence = []
    source = []
    
    def add_seq(startFork, endFork):
        
        tigId = startFork.after_id()
        s = startFork.after_pos()
        e = endFork.before_pos()
        strand = startFork.after_strand()
        src = startFork.after_switch()        

        if invert:
            tigId = startFork.before_id()
            s = startFork.before_pos()
            e = endFork.after_pos()
            strand = startFork.before_strand()
            src = startFork.before_switch()  
            
            if tigId is None or e is None:
                sequence.append("-")
                source.append("-")
                return
            
            if(startFork.before_strand() != endFork.after_strand()):
                e = len(seqData[str(endFork.after_id())]) - endFork.after_pos()
                
            if strand == 1 and s > e: (e,s) = (s,e)
            if strand == -1 and e < s: (e,s) = (s,e)

        seq = seqData[str(tigId)]
        
        if strand == 1:
            start, end = s, e
        elif strand == -1:
            start = None if e is None else (len(seq) - e)
            end = None if s is None else (len(seq) - s)
            
        segment = seq[start:end]
        if strand == -1:
            segment = reverse_complement(segment)

        sequence.append(segment)
        source.append(src*len(segment))
    
    startFork = None
    endFork = path.head()       
    
    for fork in path:
        startFork = endFork
        endFork = fork
        add_seq(startFork, endFork)

    add_seq(endFork, path.tail())
    
    if file is not None:
        
        f = open(file, "w+")
        f.write(">hybrid" + "\n")
        f.write(re.sub("(.{64})", "\\1\n", "".join(sequence), 0, re.DOTALL) + "\n")
        f.close()   
        return 
    
    return (sequence, source)


def path_to_sequence2(path, seqData):
    sequence = []
    source = []
    
    def add_seq(startFork, endFork):
        
        if startFork.is_Nfork() or endFork.is_Nfork():
            sequence.append("N"*32)
            source.append("N"*32)
            return
        
        tigId = startFork.after_id()
        start = startFork.after_pos()
        end = endFork.before_pos()
        strand = startFork.after_strand()
        src = startFork.after_switch()        
        
        if(startFork.after_strand() != endFork.before_strand()):
            print(startFork)
            print(endFork)

            print("strand issue")
            #input()
            
        seq = seqData[str(tigId)]
        
        if strand == -1:
            start = len(seq) - start
            end = len(seq) - end
            start, end = end, start
            
        segment = seq[start:end]
        if strand == -1:
            segment = reverse_complement(segment)

        sequence.append(segment)
        source.append(src*len(segment))
    
    startFork = None
    endFork = path[0]     
    
    for fork in path[1:]:
        startFork = endFork
        endFork = fork
        add_seq(startFork, endFork)

    #add_seq(endFork, path.tail())
    
    return (sequence, source)


def output_contigs(tigList, seqData, file):

    f = open(file, "w+")
    
    for tigId in tigList:
        tigId = str(tigId)
        f.write(">" + tigId + "\n")
        f.write(re.sub("(.{64})", "\\1\n", seqData[tigId], 0, re.DOTALL) + "\n")
    f.close()   
    return



def nfix_report_output(param):
    successes=0
    fails=0
    reasons=dict()
    left=0
    right=0

    for reportSet in param.reports:
        
        if reportSet.type == log.NFIX_ATTEMPT:
            print(reportSet.reports)
            if reportSet.success():
                successes = successes + 1
            else:
                fails = fails + 1
                for report in reportSet:
                    if not report.success:
                        reason = report.details[log.REASON]
                        if reason in reasons: reasons[reason] = reasons[reason] + 1
                        else: reasons[reason] = 1
                        
                        if report.details[log.SIDE] == log.RIGHT:
                            right = right + 1
                        elif report.details[log.SIDE] == log.LEFT:
                            left = left + 1
            
    print("--N FIXING RESULTS--")
    print("successes = "  + str(successes))
    print("fails = "  + str(fails))
    print("left fails = "  + str(left) + "    right fails = " + str(right))
    for key in reasons.keys():
        print(str(key) + " = "  + str(reasons[key]))



def scaffold_report_output(param):
    successes=0
    fails=0
    reasons=dict()

    for reportSet in param.reports:
        if reportSet.type == log.SCAFFOLD_ATTEMPT:
            if reportSet.has_success():
                successes = successes + 1
            else:
                fails = fails + 1
                report = reportSet[-1]
                reason = report.details[log.REASON]
                if reason in reasons: reasons[reason] = reasons[reason] + 1
                else: reasons[reason] = 1
            
                if reason == log.OVERLAP:
                    print (reportSet.details)
                    
            
    print("--SCAFFOLDING RESULTS--")
    print("successes = "  + str(successes))
    print("fails = "  + str(fails))
    for key in reasons.keys():
        print(str(key) + " = "  + str(reasons[key]))
        
    
    
    
    
    


