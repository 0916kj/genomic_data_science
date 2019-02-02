### CSE181 Week 3: Identifying Regulatory Motifs & Analysis of
### Motif Matrices
### Searching for Regulatory Motifs
    #### motifs appear in upstream regions of genome
    #### however if these motifs are not well conserved, they may
    ###have multiple mutations, and the length of the motifs and number
    ### of mismatches may make the frequent words with mismatches
    ### algorithm too slow, so we need a new approach
    #### also, motifs do not clump together, but appear once at several
    #### different regions of genome
import week1_frequent_words as week1
import week2_frequent_words_w_mismatches as week2
import numpy as np
import math

def probability_of_kmer_in_seq(sequence,k):
    ## return expected value for no. of occurences of k-mer in sequence
    ### assume random probability of 0.25 for each base
    prob = 0.25**k
    num_of_kmers = len(sequence) - k + 1
    return prob*num_of_kmers
            
### Brute force algorithm for finding k-mer motifs
    ### Given a collection of strings Dna and an integer d, 
    ### a k-mer is a (k,d)-motif if it appears in every string 
    ### from Dna with at most d mismatches

def motif_enumeration(dna, k, d):
    ### dna is a list of sequences
    ### return k-mer motifs found in each seq with <= d mismatches
    patterns = []
    for i in range(0,len(dna[0])-k+1):
        ##initialize with each k-mer found in dna
        pattern = dna[0][i:i+k]
        neighborhood = week2.neighbors(pattern,d)
        ##get list of all d-neighbors of k-mer
        for neighbor in neighborhood:
            neighbor_in_seqs = [0]*len(dna)
            neighborhood2 = week2.neighbors(neighbor,d)
            ## test if neighbor in each sequence with <=d mismatches
            for neighbor2 in neighborhood2:
                for i in range(0,len(dna)):
                    if neighbor2 in dna[i]:
                        neighbor_in_seqs[i] = 1
            #print(neighbor_in_seqs)
            #print("\n")
            #if all in neighbor_in_seqs == 1:
             #   patterns.append(neighbor)
            if 0 not in neighbor_in_seqs:
                patterns.append(neighbor)
    return set(patterns)

def motifs_matrix(list_of_motifs):
    ## return numpy array where rows are motifs being compared
    ## list of motifs can be an iterator of motifs stored as strings
    nested_list = []
    for motif in list_of_motifs:
        nested_list.append(list(motif))
    #t = len(nested_list)
    same_length = True
    for motif in nested_list:
        for element in nested_list:
            if len(motif) != len(element):
                same_length = False
        if same_length == False:        
            print("Motifs are not of equal length")
            break
        else:
            #length = len(motif)
            motif_array = np.array(nested_list)
    ### create t x length 2D array from nested list
    #print(motif_array.shape)
    return(motif_array)

def score_motifs(motifs_matrix):
    ## return numeric score equal to the sum mismatches in each column 
    ## store columns as nested list
    score = 0
    t = len(motifs_matrix)
    length = len(motifs_matrix[0])
    columns = []
    for element in motifs_matrix[0]:
        columns.append([])
        ## yield columns list of length length
    for i in range(0,t):
        for n in range(0,length):
            columns[n].append(motifs_matrix[i][n])
    for column in columns:
        count_list=[column.count("A"),column.count("C"),column.count("G"),
                    column.count("T")]
        max_count = max(count_list)
        for i in range(0,len(count_list)):
            if count_list[i] == max_count:
                if i == 0:
                    most_freq_base = "A"
                elif i == 1:
                    most_freq_base = "C"
                elif i == 2:
                    most_freq_base = "G"
                elif i == 3:
                    most_freq_base = "T"
                break
        for base in column:
            if base != most_freq_base:
                 score += 1
    return score
                
## We can construct the 4 × Length count matrix Count(Motifs) 
### counting the number of occurrences of each nucleotide in each column 
### of the motif matrix; the (i, j)-th element of Count(Motifs) stores 
### the number of times that nucleotide i appears in column j of Motifs
    
def count_motifs(motifs_matrix):
    ## return 4 x Length array of # of times base A, C, G, or T
    ### appeared at motif[j]
    ### row 0 is A, 1 is C, 2 is G, 3 is T
    t = len(motifs_matrix)
    length = len(motifs_matrix[0])
    count_array = np.zeros((4, length))
    for i in range(0,t):
        for n in range(0,length):
            base = motifs_matrix[i][n]
            if base == "A":
                count_array[0][n] += 1
            elif base == "C":
                count_array[1][n] += 1
            elif base == "G":
                count_array[2][n] += 1
            elif base == "T":
                count_array[3][n] += 1

    return count_array
            
def profile_motifs(motifs_matrix):
####return array of frequency of base i at position j
### based on count_motifs 
    t = len(motifs_matrix)
    count_array = count_motifs(motifs_matrix)
    profile_array = count_array / t
    return profile_array

#### Finally, we form a consensus string, denoted Consensus(Motifs), 
#### from the most popular letters in each column of the motif matrix.
### If we select Motifs correctly from the collection of upstream regions, 
#### then Consensus(Motifs) provides an ideal candidate regulatory motif 
#### for these regions.    
    
def consensus_motifs(motifs_matrix):
    ### return string comprised of most frequent base at each position
    consensus_motif = ""
    t = len(motifs_matrix)
    length = len(motifs_matrix[0])
    columns = []
    for element in motifs_matrix[0]:
        columns.append([])
        ## yield columns list of length length
    for i in range(0,t):
        for n in range(0,length):
            columns[n].append(motifs_matrix[i][n])
    from statistics import mode
    for column in columns:
        most_freq_base = mode(column)
        consensus_motif += (most_freq_base)
    return consensus_motif

### However, it doesn't make sense to score columns purely by mismatches
    ### For many biological motifs, certain positions feature two nucleotides
    ### with roughly the same ability to bind to a transcription factor
    ### a more appropriate representation of the consensus string 
    ### should include viable alternatives to the most popular nucleotides 
    ### in each column
    
### To help determine the "conservedness" of a position, we can
    ### calculate the entropy of the probability distribution of each column
    ### (the prob distr of each column ~= 1)
    ### Entropy is a measure of the uncertainty of a probability distribution
    ### and is defined as follows:
        #### S = -sum(Pi*log_2(Pi))

### entropy offers an improved method of scoring motif matrices: 
#### the entropy of a motif matrix is defined as the sum of the entropies
#### of its columns.

def entropy_of_prob_distr(prob_distr):
    # input prob distribution as list of floats
    ## use one column of profile matrix
    entropy = 0
    for prob in prob_distr:
        if prob > 0:
            entropy_i = math.log2(prob)
            entropy_i *= prob
        else:
            entropy_i = 0
        entropy -= entropy_i
    return entropy

def entropy_of_motif_matrix(motif_matrix):
    profile = profile_motifs(motif_matrix)
    length = len(profile[0])
    ## length = number of columns
    total_entropy = 0
    for i in range(0,length):
        column = profile[:,i]
        column_entropy = entropy_of_prob_distr(column)
        total_entropy += column_entropy
    return total_entropy

NF_kB_motifs=["TCGGGGGTTTTT","CCGGTGACTTAC",
              "ACGGGGATTTTC","TTGGGGACTTTT",
              "AAGGGGACTTCC","TTGGGGACTTCC",
              "TCGGGGATTCAT","TCGGGGATTCCT",
              "TAGGGGAACTAC","TCGGGTATAACC"]

nf_kb_matrix = motifs_matrix(NF_kB_motifs)
profile_matrix = profile_motifs(nf_kb_matrix)
#print(profile_matrix[:,1])

#print(entropy_of_prob_distr(profile_matrix[:,0]))
#print(entropy_of_motif_matrix(nf_kb_matrix))
#print(consensus_motifs(nf_kb_matrix))

### Motif Finding Problem: Given a collection of strings, 
### find a set of k-mers, one from each string, that minimizes
### the score of the resulting motif.
     ### Input: A collection of strings Dna and an integer k.
     ### Output: A collection Motifs of k-mers, one from each string in Dna, 
     ##### minimizing Score(Motifs) among all possible choices of k-mers.
     
### Brute force search takes too long to solve for minimum score(motifs)
     ### better method:  explore all potential k-mer consensus strings first 
     ### and then find the best possible collection Motifs 
     ### for each consensus string: consensus -> motifs
     
### instead of computing score column-by-column, we can compute it row-by-row
     ### the score(i) of each row = the hamming distance between the 
     ### consensus string and the motif in row i of the matrix
     ### if d(Pattern, Motifs) = sum[i=0,t](hamming_distance(Pattern,Motif_i))
     ### search for Pattern that minimizes d (Pattern is now consensus)
###  Input: A collection of strings Dna and an integer k.
     ### Output: A k-mer Pattern and a collection of k-mers, 
     #### one from each string in Dna, minimizing d(Pattern, Motifs) 
     ### among all possible choices of Pattern and Motifs.﻿
     
def d_matrix(pattern,motifs):
    ## return d score for pattern and collection of motifs
    ## pattern is a potential consensus string
    ## dna is an iterable containing strings that contain motifs
    d_total = 0
    for motif in motifs:
        d = week2.hamming_distance(pattern,motif)
        d_total += d
    return d_total

def motifs(pattern,dna):
    ## return collection of motifs that minimize d_matrix
    ## find k-mer in each string of dna individually that minimizes 
        ## hamming_distance(pattern,motif)
    k = len(pattern)
    length = len(dna[0])
    motifs = []
    for string in dna:
        hamming_distances = []
        for i in range(0,length-k+1):
            k_mer = string[i:i+k]
            hamming_distances.append(week2.hamming_distance(pattern,k_mer))
        minimum = min(hamming_distances)
        for i in range(0,len(hamming_distances)+1):
            if hamming_distances[i] == minimum:
                motifs.append(string[i:i+k])
                break
    return motifs
                
dna_1 = ["TTACCTTAAC","GATATCTGTC","ACGGCGTTCG",
         "CCCTAAAGAG","CGTCAGAGGT"]
dna_2 = ["AAATTGACGCAT", "GACGACCACGTT", "CGTCAGCGCCTG",
"GCTGAGCACCGG", "AGTACGGGACAG"]

#print(motifs("AAA",dna_1))
### Our goal is to find a k-mer Pattern that minimizes d(Pattern, Dna)
### over all k-mers Pattern, the same task that the 
### Equivalent Motif Finding Problem is trying to achieve. 
### We call such a k-mer a median string for Dna
### generate patterns to search via motifs(pattern,dna), then find
    ### the set motifs and pattern that minimizes d_matrix
    
### using brute force approach:
    


def median_string(dna,k):
    distance = 1000
    for i in range(0,4**k):
        pattern = week1.NumberToPattern(i,k)
        motifs_list = motifs(pattern,dna)
        if distance > d_matrix(pattern,motifs_list):
            distance = d_matrix(pattern,motifs_list)
            median = pattern
    return median

#print(median_string(dna_2,3))    
from itertools import product
bases = 'A', 'C', 'G', 'T'

def MedianString(dna, k):
    distance = 1000
    kmers = [''.join(i) for i in product(bases, repeat = k)]
    for pattern in kmers:
        motifs_list = motifs(pattern,dna)
        if distance > d_matrix(pattern,motifs_list):
            distance = d_matrix(pattern,motifs_list)
            median = pattern
    return median

### Unfortunately, since MedianString has to consider 4k k-mers, 
###     it becomes too slow for the Subtle Motif Problem, for which k = 15
### We have thus far assumed that the value of k is known in advance, 
    ### which is not the case in practice. 
    ### As a result, we are forced to run our motif finding algorithms 
    ### for different values of k and then try to deduce the correct motif 
    ### length. Since some regulatory motifs are rather long — later in the 
    ### chapter, we will search for a biologically important motif 
    ### of length 20 — MedianString may be too slow to find them.

### THE GREEDY MOTIF SEARCH ###
    ### Many algorithms are iterative procedures that must choose 
    ### among many alternatives at each iteration. 
    ### Some of these alternatives may lead to correct solutions, 
    ### whereas others may not. 
    ### Greedy algorithms select the “most attractive”
    ### alternative at each iteration.
    
### Probability(k-mer,motifs) = product of probability of each base
    ### in k-mer appearing in its position according to profile(motifs)
    ### A k-mer tends to have a higher probability when it is more similar 
    ### to the consensus string of a profile.
    
def probability(k_mer,profile):
    ### motifs is a motifs matrix
    prob = 1
    for i in range(0,len(k_mer)):
        if k_mer[i] == "A":
            prob *= profile[0][i]
        elif k_mer[i] == "C":
            prob *= profile[1][i]
        elif k_mer[i] == "G":
            prob *= profile[2][i]
        elif k_mer[i] == "T":
            prob *= profile[3][i]
    return prob

#motifs1 = motifs_matrix(dna_1)
#print(probability("TTACCTTAAG",motifs1))

### Given a profile matrix Profile, we can evaluate the probability 
### of every k-mer in a string Text and find a Profile-most probable k-mer 
### in Text, i.e., a k-mer that was most likely to have been generated 
### by Profile among all k-mers in Text.

def probability2(k_mer,profile):
    ## assume profile is dictionary where keys = bases
    prob = 1
    for i in range(0,len(k_mer)):
        base = k_mer[i]
        prob *= profile[base][i]
    return prob

def profile_most_probable_kmer(text, k, profile):
    prob_list = []
    for i in range(0,len(text)-k+1):
        k_mer = text[i:i+k]
        prob = probability2(k_mer,profile)
        prob_list.append(prob)
    max_prob = max(prob_list)
    for n in range(0,len(prob_list)):
        if prob_list[n] == max_prob:
            return text[n:n+k]

def profile_most_probable_kmer2(text, k, profile):
    ## use profile array not dict
    prob_list = []
    for i in range(0,len(text)-k+1):
        k_mer = text[i:i+k]
        prob = probability(k_mer,profile)
        prob_list.append(prob)
    max_prob = max(prob_list)
    for n in range(0,len(prob_list)):
        if prob_list[n] == max_prob:
            return text[n:n+k]

def greedy_motif_search(dna, k, t):
    motifs_list = []
    for string in dna:
        motifs_list.append(string[0:k])
        ## initialize motif array as first k-mer of each string in dna
    best_motifs = motifs_matrix(motifs_list)
    motifs_list2 = []
    for i in range(0, len(dna[0])-k+1):
        motifs_list2 = []
        motifs_list2.append(dna[0][i:i+k])
        ## create motifs array, make first motif in each array each k-mer
        ## in first dna string
        for n in range(1,t):
            ## add motifs to array
            ## motif added is profile most probable k-mer for profile
            ## created from motif matrix up to motif_matrix[n-1]
            motifs_profile = profile_motifs(motifs_matrix(motifs_list2))
            motifs_list2.append(profile_most_probable_kmer2(dna[n],k,motifs_profile))
        if score_motifs(motifs_matrix(motifs_list2)) < score_motifs(best_motifs):
            ## use score to determine if each array is better than current
            ## "best" array
            best_motifs = motifs_matrix(motifs_list2)
    return best_motifs
            
### Because Profile(motifs) can contain zeroes, which makes probability = 0,
### greedy motif search is often inaccurate
    ### to combat this, bioinformaticists often replace 0's with
    ### pseudocounts
    ### In the case of motifs, pseudocounts often amount to adding 1 
    ### (or some other small number) to each element of Count(Motifs)
    
def pseudocount_motifs(motifs_matrix,pseudocount):
    ## return 4 x Length array of # of times base A, C, G, or T
    ### appeared at motif[j]
    ### row 0 is A, 1 is C, 2 is G, 3 is T
    t = len(motifs_matrix)
    length = len(motifs_matrix[0])
    count_array = np.full((4, length),pseudocount)
    for i in range(0,t):
        for n in range(0,length):
            base = motifs_matrix[i][n]
            if base == "A":
                count_array[0][n] += 1
            elif base == "C":
                count_array[1][n] += 1
            elif base == "G":
                count_array[2][n] += 1
            elif base == "T":
                count_array[3][n] += 1
    return count_array

def pseudocount_profile_motifs(motifs_matrix,pseudocount):
####return array of frequency of base i at position j
### based on count_motifs 
    t = len(motifs_matrix)
    count_array = pseudocount_motifs(motifs_matrix,pseudocount)
    profile_array = count_array / t
    return profile_array
           
def greedy_motif_search_w_pseudocounts(dna, k, t, pseudocount):
    motifs_list = []
    for string in dna:
        motifs_list.append(string[0:k])
        ## initialize motif array as first k-mer of each string in dna
    best_motifs = motifs_list
    motifs_list2 = []
    for i in range(0, len(dna[0])-k+1):
        motifs_list2 = []
        motifs_list2.append(dna[0][i:i+k])
        ## create motifs array, make first motif in each array each k-mer
        ## in first dna string
        for n in range(1,t):
            ## add motifs to array
            ## motif added is profile most probable k-mer for profile
            ## created from motif matrix up to motif_matrix[n-1]
            motifs_profile = pseudocount_profile_motifs(motifs_matrix(motifs_list2),pseudocount)
            motifs_list2.append(profile_most_probable_kmer2(dna[n],k,motifs_profile))
        if score_motifs(motifs_matrix(motifs_list2)) < score_motifs(motifs_matrix(best_motifs)):
            ## use score to determine if each array is better than current
            ## "best" array
            best_motifs = motifs_list2
    return best_motifs
