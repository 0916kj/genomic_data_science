### CSE181: Week 1 ###
### Objective: Locating the ori in an organism's genome
### Heavily repeated sequences in the same region of the genome
### indicate the point of origin of replication in some organisms,
### such as vibrio cholerae
### write code to identify frequently repeated sequences in
### the vibrio cholerae and E Coli genomes

import numpy as np
import pandas as pd

def count_substring(string, substring):
    ##count number of occurrences of a sub-sequence in 
    ## a nucleotide sequence
    length=len(string)
    sub_length=len(substring)
    count=0
    for i in range(0,length-sub_length+1):
        if string[i:i+sub_length]==substring:
            count+=1
    return count

def FrequentWords(Text, k):
    frequent_patterns=[]
    counts=[]
    length = len(Text)
    for i in range(0,length-k+1):
        pattern = Text[i:i+k]
        counts.append(count_substring(Text,pattern))
    max_count=max(counts)
    for i in range(0,length-k+1):
        pattern = Text[i:i+k]
        if count_substring(Text,pattern) == max_count:
            if pattern not in frequent_patterns:
                frequent_patterns.append(pattern)
    return frequent_patterns
## The brute force FrequentWords search is too slow for practical use
    ## a faster approach: convert k-mers to number, where each number
    ## in range 0 -> 4**k represents one of all possible k-mers
    ## creating range 0 -> 4**k also generates all possible k-mers

def SymbolToNumber(symbol):
    symbol_dict={'A':0,'C':1,'G':2,'T':3}
    return symbol_dict[symbol]

def NumberToSymbol(number):
    number_dict = {0:'A',1:'C',2:'G',3:'T'}
    return number_dict[number]

def PatternToNumber(Pattern):
    is_pattern=True
    for letter in {'A','C','G','T'}:
        if letter in Pattern:
            is_pattern=True
            break
        else:
            is_pattern=False
    if is_pattern == False:
        return 0            
    symbol = Pattern[-1]
    Prefix = Pattern[0:-1]
    return 4*PatternToNumber(Prefix) + SymbolToNumber(symbol)

def NumberToPattern(index,k):
    pattern = ''
    while len(pattern)<k:
        quotient = index//4
        remainder = index%4
        symbol = NumberToSymbol(remainder)
        pattern = symbol + pattern
        index = quotient
    return pattern

def ComputingFrequencies(Text, k):
    ## generate numpy array representing all k-mers in Text
    ## and their frequency
    FrequencyArray={}
    num_permutations=4**k
    for i in range(0,num_permutations):
        FrequencyArray[i] = 0
    text_length = len(Text)-k+1
    for i in range(0,text_length):
        Pattern = Text[i:i+k]
        j = PatternToNumber(Pattern)
        if j in FrequencyArray.keys():
            FrequencyArray[j] += 1    
        else:
            FrequencyArray[j] = 1
    return np.array(list(FrequencyArray.items()),dtype=object) 

def FasterFrequentWords(Text, k):
    ## return set of most frequent k-mers based on frequency array
    FrequentPatterns = []
    FrequencyArray = ComputingFrequencies(Text, k)
    maxCount = np.amax(FrequencyArray,axis=0)[1]
    num_permutations=4**k
    for i in range(0,num_permutations):
        if FrequencyArray[i][1] == maxCount:
            Pattern = NumberToPattern(i,k)
            FrequentPatterns.append(Pattern)
    return set(FrequentPatterns)

def frequency_df(text,k):
    ## return Pandas dataframe of frequency array
    array = ComputingFrequencies(text,k)
    pattern_list=[]
    for number in array[:,0]:
        pattern_list.append(NumberToPattern(number,k))
    num_permutations = 4**k
    pattern_array=np.array(pattern_list,dtype=object).reshape(num_permutations,1)
    new_array = np.append(array, pattern_array, 1)
    frequency_df = pd.DataFrame(data=new_array[:,1:],index=new_array[:,0],columns=['Frequency','Pattern'])
    frequency_df = frequency_df[["Pattern","Frequency"]]
    return frequency_df
### Use DataFrame.sort_values(["Frequency"],ascending=False) to sort
    #### from greatest to lowest frequency

def ClumpFinding(Genome, k, t, L):
    ## window search of genome: return frequent k-mers
    ## t is threshold for k-mer frequency
    ## L is window size for genome search
    FrequentPatterns = []
    Clump = []
    for i in range(0,4**k-1):
        Clump.append(0)
    genome_length = len(Genome)
    for i in range(0, genome_length-L):
        Text = Genome[i:i+L]
        FrequencyArray = ComputingFrequencies(Text, k)
        ##create frequency array for each "window" of genome
        for index in range(0,4**k-1):
            if FrequencyArray[index][1] >= t:
                Clump[index] = 1
        for i in range(0,4**k-1):
            if Clump[i] == 1:
                Pattern = NumberToPattern(i, k)
                FrequentPatterns.append(Pattern)
        return FrequentPatterns

def BetterClumpFinding(Genome, k, t, L):
    FrequentPatterns = []
    Clump = []
    genome_length = len(Genome)
    Text = Genome[0:L]
    FrequencyArray = ComputingFrequencies(Text, k)
    for i in range(0,4**k-1):
        ##generate clump array as you search for frequent k-mers
        if FrequencyArray[i][1] >= t:
            Clump.append(1)
        else:
            Clump.append(0)
    for i in range(1,genome_length - L):
        FirstPattern = Genome[i-1:i+k-1]
        index = PatternToNumber(FirstPattern)
        FrequencyArray[index][1] = FrequencyArray[index][1] - 1
        LastPattern = Genome[i+L-k:i+L]
        index = PatternToNumber(LastPattern)
        FrequencyArray[index][1] = FrequencyArray[index][1]+1
        if FrequencyArray[index][1] >= t:
            Clump[index] = 1
    for i in range(0,4**k-1):
        if Clump[i] == 1:
            Pattern = NumberToPattern(i, k)
            FrequentPatterns.append(Pattern)
    return set(FrequentPatterns)

def ReverseComplement(Pattern):
    complement_dict={'A':'T','C':'G','G':'C','T':'A'}
    new_sequence=''
    for base in Pattern[::-1]:
        new_sequence += complement_dict[base]
    return new_sequence

def FindRCs(Genome, k, t, L):
    ### return reverse complement
    frequent_kmers = BetterClumpFinding(Genome, k, t, L)
    RCs = []
    for kmer in frequent_kmers:
        if ReverseComplement(kmer) not in RCs:
            for other_kmer in frequent_kmers:
                if ReverseComplement(kmer) == other_kmer:
                    RCs.append(kmer)
    return set(RCs)

def PatternMatching(Pattern, Genome):
    ### return starting indices of pattern within genome/longer sequences
    indices=[]
    k = len(Pattern)
    for i in range(0,len(Genome)-k+1):
        if Genome[i:i+k] == Pattern:
            indices.append(i)
    return indices

## Search real genomes for 9-mer "clumps"            
with open("vibrio_cholerae_genome.txt") as f:
    genome=f.read()
    print(BetterClumpFinding(genome,9,3,1000))

with open("e_coli.txt") as f:
     genome=f.read()
     clumps = BetterClumpFinding(genome,9,3,500)
     print(len(clumps))
     print(clumps)
     