### Genomic Data Science Week 4 ###
### random motifs search: randomly select motifs from set dna
    ### construct profile(motifs), then find motifs(profile,dna)
    ### where each motif is the profile-most-probable k-mer
    ### for each string dna
    ### compare motif scores, then iterate until score no longer increases

import week3_motif_matrices as week3
import random

def random_motifs(dna,k):
    ## choose a random k-mer from each string dna
    length = len(dna[0])
    motifs_list = []
    for string in dna:
        i = random.randint(0,length-k)
        motifs_list.append(string[i:i+k])
    return week3.motifs_matrix(motifs_list)

def motifs_from_profile(dna,k,t,profile):
    ## construct motifs matrix from given profile
    ## the motifs are the profile-most-probable (see week3) k-mer
    ## for each sequence in set dna
    motifs_list=[]
    for i in range(0,t):
        motifs_list.append(week3.profile_most_probable_kmer2(dna[i], k, profile))
    return week3.motifs_matrix(motifs_list)

def random_motifs_search(dna,k,t):
    ## iterate to find probable motifs of set dna:
    ## start with random set of motifs, then iteratively
    ## produce new motif matrix based on produced profile
    ##until highest score is achieved
    motifs = random_motifs(dna,k)
    while True:
        score = week3.score_motifs(motifs)
        ## use pseudocounts, in this case, 1
        profile = week3.pseudocount_profile_motifs(motifs,1)
        new_motifs = motifs_from_profile(dna,k,t,profile)
        if week3.score_motifs(new_motifs) < score:
            motifs = new_motifs
        else:
            return motifs

## Then, write a function here called RepeatedRandomizedMotifSearch() 
## that takes Dna, k, t, and a parameter
## N that returns the best collection of motifs 
## resulting from running RandomizedMotifSearch()
## N times.  It should return a list of strings.

def repeated_random_motifs_search(dna,k,t,n):
    ## bioinformaticians run random searches thousands
    ## of times to eliminate error
    motifs = random_motifs_search(dna,k,t)
    count = 0
    while count < n:
        score = week3.score_motifs(motifs)
        new_motifs = random_motifs_search(dna,k,t)
        if week3.score_motifs(new_motifs) < score:
            motifs = new_motifs
        count +=1
    return motifs

dna_1 = ["CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA",
"GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG",
"TAGTACCGAGACCGAAAGAAGTATACAGGCGT",
"TAGATCAAGTTTCAGGTGCACGTCGGTGAACC",
"AATCCACCAGCTCCACGTGCAATGTTGGCCTA"]

dna_2 = ["AATTGGCACATCATTATCGATAACGATTCGCCGCATTGCC",
"GGTTAACATCGAATAACTGACACCTGCTCTGGCACCGCTC",
"AATTGGCGGCGGTATAGCCAGATAGTGCCAATAATTTCCT",
"GGTTAATGGTGAAGTGTGGGTTATGGGGAAAGGCAGACTG",
"AATTGGACGGCAACTACGGTTACAACGCAGCAAGAATATT",
"GGTTAACTGTTGTTGCTAACACCGTTAAGCGACGGCAACT",
"AATTGGCCAACGTAGGCGCGGCTTGGCATCTCGGTGTGTG",
"GGTTAAAAGGCGCATCTTACTCTTTTCGCTTTCAAAAAAA"]

## find motifs of dna sets from 1000 random searches
print(repeated_random_motifs_search(dna_1,8,5,1000))
print(repeated_random_motifs_search(dna_2,6,8,1000).tolist())