import sys
sys.path.append("..")

import helper
import external_tools as tools
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

def get_tag(reads, tag):
    return np.array([read.get_tag(tag) for read in reads])

def get_tag_from_dict(dictionary, tag, average=True):
    tagDict = dict()
    for name in dictionary:
        tagDict[name] = np.mean(get_tag(dictionary[name], tag))
    return tagDict
    
def join_dicts(dictionary1, dictionary2):
    list1, list2 = [], []
    list1Unmapped, list2Unmapped = [], []

    for name in dictionary1:
        if name in dictionary2:
            list1.append(dictionary1[name])
            list2.append(dictionary2[name])
        else:
            list1Unmapped.append(dictionary1[name])
    for name in dictionary2:
        if name not in dictionary1:
            list2Unmapped.append(dictionary2[name])
    
    return (np.array(list1), np.array(list2), np.array(list1Unmapped), np.array(list2Unmapped))
     
def to_array(dictionary, withNames=False):
    arrayList = []
    nameList = []
    for name in dictionary:
        arrayList.append(dictionary[name])
        nameList.append(name)

    return np.array(arrayList) if not withNames else (np.array(arrayList), np.array(nameList))
               
def to_dict(reads):
    dictionary = dict()
    for read in reads:
        name = read.qname
        if name in dictionary:
            dictionary[name].append(read)
        else:
            dictionary[name] = [read]
    return dictionary

def plot_identity(x, y):
    low, high = min(min(x), min(y)), max(max(x), max(y))
    x=np.linspace(low,high,10) 
    plt.plot(x,x,'r-')
    
def plot_scatter(x, y, alpha=0.5, colour="c"):
    plt.scatter(x, y, alpha=alpha, color=colour)
sns.set(style="whitegrid")

def plot_hist(x, bins=40, alpha=0.5, colour=None):
    if colour is not None:
        plt.hist(x, bins=bins, alpha=0.5, color=colour)
    else:
        plt.hist(x, bins=bins, alpha=0.5)

def plot_mc(originalReads, newReads, simple=True):
    originalReadDict = to_dict(originalReads)
    newReadDict = to_dict(newReads)

    originalMCDict = get_tag_from_dict(originalReadDict, "mc")
    newMCDict = get_tag_from_dict(newReadDict, "mc")

    mc1, mc2, mc1_, mc2_ = join_dicts(originalMCDict, newMCDict)

    if not simple:
        plot_identity(mc1, mc2)
        plot_scatter(mc1, mc2)
        plt.show()
    
    plot_hist(mc1, bins=30, colour="grey")
    plot_hist(mc2, bins=30, colour="forestgreen")
    plt.show()
    print("mc")
    print(np.mean(mc1), np.mean(mc2))


def plot_cigar(originalReads, newReads, simple=True):
    originalReadDict = to_dict(originalReads)
    newReadDict = to_dict(newReads)
    
    originalCigarDict = dict()
    newCigarDict = dict()

    for name in originalReadDict:
        originalCigarDict[name] = helper.cigar_summary(originalReadDict[name], block=False)
    for name in newReadDict:
        newCigarDict[name] = helper.cigar_summary(newReadDict[name], block=False)

    def get_cigar_dict(dictionary, key):
        return { name : dictionary[name][key] for name in dictionary }
    
    '''   
    indel1, indel2, indel1Unmapped, indel2Unmapped = join_dicts(get_cigar_dict(originalCigarDict, "indel"), 
                                                get_cigar_dict(newCigarDict, "indel"))
    plot_identity(indel1, indel2)
    plot_scatter(indel1, indel2)
    '''
    total1, total2, total1Unmapped, total2Unmapped = join_dicts(get_cigar_dict(originalCigarDict, "total"), 
                                            get_cigar_dict(newCigarDict, "total"))
    
    clip1, clip2, clip1Unmapped, clip2Unmapped = join_dicts(get_cigar_dict(originalCigarDict, "clip"), 
                                            get_cigar_dict(newCigarDict, "clip"))

    norm1, norm2 = (total1 - clip1), (total2 - clip2)

        
    I1, I2, I1Unmapped, I2Unmapped = join_dicts(get_cigar_dict(originalCigarDict, "I"), 
                                                get_cigar_dict(newCigarDict, "I"))
    
    D1, D2, D1Unmapped, D2Unmapped = join_dicts(get_cigar_dict(originalCigarDict, "D"), 
                                                get_cigar_dict(newCigarDict, "D"))
    
    M1, M2, M1Unmapped, M2Unmapped = join_dicts(get_cigar_dict(originalCigarDict, "match"), 
                                                get_cigar_dict(newCigarDict, "match"))

    X1, X2, X1Unmapped, X2Unmapped = join_dicts(get_cigar_dict(originalCigarDict, "X"), 
                                                get_cigar_dict(newCigarDict, "X"))

    if not simple:

        normI1, normI2 = I1/norm1, I2/norm2
        plot_identity(normI1, normI2)
        plot_scatter(normI1, normI2)
        plt.show()
        plot_hist(normI1, colour="grey")
        plot_hist(normI2, colour="forestgreen")
        plt.show()
        print("insertion")
        print(round(np.mean(normI1),4), round(np.mean(normI2),4))
    
        normD1, normD2 = D1/norm1, D2/norm2
        plot_identity(normD1, normD2)
        plot_scatter(normD1, normD2)
        plt.show()
        plot_hist(normD1, colour="grey")
        plot_hist(normD2, colour="forestgreen")
        plt.show()
        print("deletion")
        print(round(np.mean(normD1),4), round(np.mean(normD2),4))
    
        '''
        print("")
        normX1, normX2 = X1/norm1, X2/norm2
        plot_identity(normX1, normX2)
        plot_scatter(normX1, normX2)
        plt.show()
        '''
    
    discord1, discord2 = (I1+D1+X1)/norm1, (I2+D2+X2)/norm2
    plot_identity(discord1, discord2)
    plot_scatter(discord1, discord2)
    plt.show()
    plot_hist(discord1, colour="grey")
    plot_hist(discord2, colour="forestgreen")
    plt.show()
    print("discordance")
    print(round(np.mean(discord1),4), round(np.mean(discord2),4))

    concord1, concord2 = M1/norm1, M2/norm2
    if not simple:
        plot_identity(concord1, concord2)
        plot_scatter(concord1, concord2)
        plt.show()
    plot_hist(concord1, colour="grey")
    plot_hist(concord2, colour="forestgreen")
    plt.show()
    print("concordance")
    print(round(np.mean(concord1),4), round(np.mean(concord2),4))

def calculate_difference(originalReads, newReads):
    plot_mc(originalReads, newReads)
    plot_cigar(originalReads, newReads)
    

def assess_difference(segment, seqData, param):
    
    headerBam = segment.bamA if segment.bamA is not None else segment.bamB
    headerBam = segment.bamU if headerBam is None else headerBam

    if segment.regionA is not None:
        originalFa = helper.isolate_region_fasta(segment.regionA, seqData, param)
        newFaDict = {"polished" : segment.seqA}
        prefix = helper.file_prefix(segment.regionA, param)
        newFa = tools.dict2fasta(newFaDict, prefix + "_polished_TEMP_")
        
        if segment.bamA is not None:
            originalBamA = None if segment.bamA is None else tools.align_pacbio(originalFa, segment.bamA, prefix + "_original_TEMP_A")
            newBamA = None if segment.bamA is None else tools.align_pacbio(newFa, segment.bamA, prefix + "_new_TEMP_A")
            originalReads = tools.samtools_fetch(originalBamA)
            newReads = tools.samtools_fetch(newBamA)
            calculate_difference(originalReads, newReads)
            
        if segment.bamU is not None:
            originalBamU = None if segment.bamU is None else tools.align_pacbio(originalFa, segment.bamU, prefix + "_original_TEMP_Au")
            newBamU = None if segment.bamU is None else tools.align_pacbio(newFa, segment.bamU, prefix + "_new_TEMP_Au")
            originalReads = tools.samtools_fetch(originalBamU)
            newReads = tools.samtools_fetch(newBamU)
            calculate_difference(originalReads, newReads)

    if segment.regionB is not None:
        originalFa = helper.isolate_region_fasta(segment.regionB, seqData, param)
        newFaDict = {"polished" : segment.seqB}
        prefix = helper.file_prefix(segment.regionB, param)
        newFa = tools.dict2fasta(newFaDict, prefix + "_polished_TEMP_B")
        
        if segment.bamB is not None:
            originalBamB = None if segment.bamB is None else tools.align_pacbio(originalFa, segment.bamB, prefix + "_original_TEMP_B")
            newBamB = None if segment.bamB is None else tools.align_pacbio(newFa, segment.bamB, prefix + "_new_TEMP_B")
            originalReads = tools.samtools_fetch(originalBamB)
            newReads = tools.samtools_fetch(newBamB)
            calculate_difference(originalReads, newReads)

        if segment.bamU is not None:
            originalBamU = None if segment.bamU is None else tools.align_pacbio(originalFa, segment.bamU, prefix + "_original_TEMP_Bu")
            newBamU = None if segment.bamU is None else tools.align_pacbio(newFa, segment.bamU, prefix + "_new_TEMP_Bu")
            originalReads = tools.samtools_fetch(originalBamU)
            newReads = tools.samtools_fetch(newBamU)       
            calculate_difference(originalReads, newReads)


def plot_alignments(faDict, bamDict, prefix, faDictOrder=None, image="svg"):
    
    def get_reads(fasta, bam):
        if bam is None:
            return []
        file = tools.align_pacbio(fasta, bam, prefix + "_TEMP_realigned_plot")
        reads = tools.samtools_fetch(file)
        helper.delete_file(file)
        return reads

    mc, bamName, faName, readName = [], [], [], []
    
    for faFileName in faDict:
        faFile = faDict[faFileName]
    
        for bamFileName in bamDict:
            bamFile = bamDict[bamFileName]
            reads = get_reads(faFile, bamFile)
            readDict = to_dict(reads)
            readMC,names = to_array(get_tag_from_dict(readDict, "mc"), withNames=True)
            readMC,names = list(readMC),list(names)

            mc = mc + readMC
            bamName = bamName + [bamFileName]*len(readMC)
            faName = faName + [faFileName]*len(readMC)
            readName = readName +  names
            
    plotDf = pd.DataFrame({"mc" : mc, "alignment" : bamName, "reference" : faName, "read" : readName})
    
    palette = "muted"
    sns.set(rc={'figure.figsize':(11.7,8.27)})
    sns.set_style("whitegrid")
    
    data = plotDf.groupby(['alignment', 'reference']).mean().reset_index()
    #sd = plotDf.groupby(['alignment', 'reference']).std().reset_index()
    #data["mc_sd"] = sd["mc"]
    mean = np.mean(plotDf["mc"])
    
    plt.figure(figsize=(12, 8))
    ax = sns.pointplot(x="reference", y="mc",data=plotDf, palette=palette,
                         hue="alignment", size=10, linewidth=1, edgecolor='black',
                         order=faDictOrder, dodge=True, join=False, capsize=.1,
                         errwidth=2)
    ax.axhline(y = mean, color='red', linewidth=2, alpha=.4)
    ax.set(xlabel="", ylabel="Average Mapped Concordance Percentage")
    plt.legend(title='Reads')
    figure = ax.get_figure()
    figure.savefig(prefix + "_average_mapped_concordance." + image)
    
    for faFileName in faDict:

        data = plotDf[plotDf["reference"] == faFileName]
        plt.figure(figsize=(12, 8))

        sns.swarmplot("alignment", y="mc", data=data,
                 color="black", edgecolor="black", size=5, alpha=0.3)

        ax = sns.violinplot(x="alignment", y="mc", 
                             data=data, palette=palette, inner="quartile")
        ax.set_title("Reads Aligned Against " + faFileName)
        ax.set(xlabel="Reads", ylabel="Mapped Concordance Percentage")

        figure = ax.get_figure()
        figure.savefig(prefix + "_" + faFileName + "_violin." + image)

    
        plt.figure(figsize=(12, 8))
        ax = sns.boxplot(x="alignment", y="mc", data=data, palette=palette)
        ax.set_title("Reads Aligned Against " + faFileName)
        ax.set(xlabel="Reads", ylabel="Mapped Concordance Percentage")
        figure = ax.get_figure()
        figure.savefig(prefix + "_" + faFileName + "_boxplot." + image)



    
