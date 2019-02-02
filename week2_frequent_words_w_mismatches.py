### Intro to Genomic Data Science Week 2
### Identify frequent words as well as close mismatches and 
### their reverse complements to locate ori
#### Asymmetrical rates of DNA replication - use to find ori

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import  week1_frequent_words as week_1

#### Skew: measure and plot difference between occurences of G and C
##### up to base i in genome
#### C decreases on the forward/lagging strand in the 5'-3' direction
##### due to lagging strand's vulnerability to deamination mutation
#### The ori is located where G-C goes from decreasing to increasing (minimum)

def skew_i(genome,i):
    ###return list of numbers representing G-C up to point i
    skew_list=[0] ##skew_0 always = 0
    for k in range(0,i):
        if genome[k]=="C":
            skew_list.append(skew_list[k]-1)
        elif genome[k]=="G":
            skew_list.append(skew_list[k]+1)
        else:
            skew_list.append(skew_list[k])
    return skew_list
    
def plot_skew(skew_list):
    ##visualize skew vs location in genome
    plt.plot(skew_list)
    plt.ylabel('Skew')
    plt.xlabel('Locn in Genome')
    plt.show()
    
def minimum_skew(genome):
    ###Return list of indices where skew is minimized
    skew_list=skew_i(genome,len(genome))
    minimum = min(skew_list)
    indices=[]
    for i in range(0,len(skew_list)):
        if skew_list[i] == minimum:
            indices.append(i)
    return indices

def hamming_distance(p,q):
    ###calculate # of mismatches between two strings p and q
    if len(p) != len(q):
        print("Strings are different lengths")
    else:
        mismatch_count = 0
        length = len(p)
        for i in range(0,length):
            if p[i] != q[i]:
                mismatch_count += 1
            else:
                continue
        return mismatch_count
    
def approx_pattern_matching(text, pattern, d):
    ##return starting indices of subpatterns <= d mismatches
    ###from pattern in text
    length = len(pattern)
    indices=[]
    for i in range(0,len(text)-length+1):
        pattern2 = text[i:i+length]
        if hamming_distance(pattern,pattern2) <= d:
            indices.append(i)
    return indices
 
def count_d(text, pattern, d):
    ###count total appearances of pattern with <=d mismatches in text
    list_of_approxs = approx_pattern_matching(text,pattern,d)
    return len(list_of_approxs)

"""
Addressing the frequent words problem again:
    A most frequent k-mer with up to d mismatches in Text is simply
    a string Pattern maximizing Countd(Text, Pattern) among all k-mers. 
    Note that Pattern does not need to actually appear as a substring of Text
"""

def immediate_neighbors(pattern):
    ### return set of patterns with 1 mismatch compared to input
    neighborhood = [pattern]
    for i in range(0,len(pattern)):
        symbol = pattern[i]
        pattern_list=list(pattern)
        for nucleotide in "ACGT":
            if not nucleotide == symbol:
                neighbor = pattern_list
                neighbor[i] = nucleotide
                neighborhood.append(''.join(neighbor))
    return set(neighborhood)

def neighbors(pattern,d):
    ### return set of all patterns with <= d mismatches
    if d == 0:
        return {pattern}
    if len(pattern) == 1:
        return {'A','C','G','T'}
    neighborhood = [pattern]
    suffix = pattern[1:]
    first_symbol = pattern[0]
    ###create set of neighbors to suffix of pattern
    ####these neighbors have <=d mismatches with pattern
    #####when we put the first symbol in front
    suffix_neighbors = neighbors(suffix, d)
    for neighbor in suffix_neighbors:
        if hamming_distance(neighbor,suffix) == d:
            neighborhood.append(first_symbol + neighbor)
        elif hamming_distance(neighbor,suffix) < d:
            ###for suffix neighbors with <d mismatches, we can put
            ####any symbol in front and have <=d mismatches total
            for symbol in ['A','C','G','T']:
                 neighborhood.append(symbol + neighbor)
    return set(neighborhood)

###Use neighbors to create array of frquent words and mismatches

def computing_frequencies_with_mismatches(text, k, d):
    ### return numpy array containing frequencies of k-mers
    #### count mismatches <=d in frequency count
    FrequencyArray={}
    num_permutations=4**k
    for i in range(0,num_permutations):
        FrequencyArray[i] = 0
    for i in range(0,len(text)-k+1):
        pattern = text[i:i+k]
        neighborhood = neighbors(pattern,d)
        for neighbor in neighborhood:    
            j = week_1.PatternToNumber(neighbor)
            FrequencyArray[j] += 1    
    return np.array(list(FrequencyArray.items()),dtype=object)

def mismatch_frequency_df(text,k,d):
    ##return dataframe containing frequency array with
    ###string pattern as a column
    array = computing_frequencies_with_mismatches(text,k,d)
    pattern_list=[]
    for number in array[:,0]:
        pattern_list.append(week_1.NumberToPattern(number,k))
    num_permutations = 4**k
    pattern_array=np.array(pattern_list,dtype=object).reshape(num_permutations,1)
    new_array = np.append(array, pattern_array, 1)
    frequency_df = pd.DataFrame(data=new_array[:,1:],index=new_array[:,0],columns=['Frequency','Pattern'])
    frequency_df = frequency_df[["Pattern","Frequency"]]
    return frequency_df

def frequent_words_with_mismatches(text, k, d):
    ##return most frequent words including mismatches <= d
    FrequentPatterns = []
    FrequencyArray = computing_frequencies_with_mismatches(text, k, d)
    maxCount = np.amax(FrequencyArray,axis=0)[1]
    num_permutations=4**k
    for i in range(0,num_permutations):
        if FrequencyArray[i][1] == maxCount:
            Pattern = week_1.NumberToPattern(i,k)
            FrequentPatterns.append(Pattern)
    return set(FrequentPatterns)

def frequent_words_w_mismatches_sorting(text, k, d):
    ##return frequent patterns, sorted by "number" (ie PatterntoNumber)
    ## of pattern
    frequent_patterns = []
    neighborhoods = []
    for i in range(0,len(text)-k+1):
        for neighbor in neighbors(text[i:i+k],d):
            neighborhoods.append(neighbor)
    neighborhood_dict = {"Pattern":[],"Index":[],"Count":[]}
    for pattern in neighborhoods:
        neighborhood_dict["Pattern"].append(pattern)
        neighborhood_dict["Index"].append(week_1.PatternToNumber(pattern))
        neighborhood_dict["Count"].append(1)
    sorted_index = sorted(neighborhood_dict["Index"])
    for i in range(0,len(neighborhoods)-1):
        if sorted_index[i] == sorted_index[i+1]:
            neighborhood_dict["Count"][i+1] = neighborhood_dict["Count"][i] + 1
    max_count = max(neighborhood_dict["Count"])
    for i in range(0,len(neighborhoods)):
        if neighborhood_dict["Count"][i] == max_count:
            pattern = week_1.NumberToPattern(sorted_index[i], k)
            frequent_patterns.append(pattern)
    return frequent_patterns

### Now find frequent words including mismatches AND reverse complements
def computing_freqs_w_mismatches_and_rcs(text, k, d):
    ### return numpy array containing frequencies of k-mers
    #### count mismatches <=d in frequency count
    count=[]
    count_rc=[]
    num_permutations=4**k
    for i in range(0,num_permutations):
        count.append(0)
        count_rc.append(0)
    for i in range(0,len(text)-k+1):
        pattern = text[i:i+k]
        neighborhood = neighbors(pattern,d)
        for neighbor in neighborhood:
            if neighbor != pattern:
                j = week_1.PatternToNumber(neighbor)
                count[j] += 1
        complement = week_1.ReverseComplement(pattern)
        if complement != pattern and complement not in neighborhood:
            m = week_1.PatternToNumber(complement)
            count_rc[m] += 1
    total_count=[]
    for i in range(0,len(count)):
        total_count.append(count[i]+count_rc[i])
    return total_count

def freq_words_w_mismatches_and_rcs(text, k, d):
    FrequentPatterns = []
    count = computing_freqs_w_mismatches_and_rcs(text, k, d)
    max_count = max(count)
    num_permutations=4**k
    for i in range(0,num_permutations):
        if count[i] == max_count:
            Pattern = week_1.NumberToPattern(i,k)
            FrequentPatterns.append(Pattern)
    return set(FrequentPatterns)