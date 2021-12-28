# -*- coding: utf-8 -*-
"""
THE GENE FINDER PROGRAM
@author Brenda Muthoni Karumbo
"""

import random
from amino_acids import aa, codons, aa_table   # you may find these useful
from load import load_seq


def shuffle_string(s):
    """Shuffles the characters in the input string
        NOTE: this is a helper function, you do not
        have to modify this in any way """
    return ''.join(random.sample(s, len(s)))

# YOU WILL START YOUR IMPLEMENTATION FROM HERE DOWN ##
def get_complement(nucleotide):
    """ Returns the complementary nucleotide
        nucleotide: a nucleotide (A, C, G, or T) represented as a string
        returns: the complementary nucleotide
    >>> get_complement('A')
    'T'
    >>> get_complement('C')
    'G'
    """
    complement = nucleotide.replace("A","t")
    complement = complement.replace("T","a")
    complement = complement.replace("C","g")
    complement = complement.replace("G","c")
    complement = complement.upper()
    return complement

    # TODO: implement this
    pass


def get_reverse_complement(dna):
    complement = dna.replace("A","t")
    complement = complement.replace("T","a")
    complement = complement.replace("C","g")
    complement = complement.replace("G","c")
    complement = complement.upper()
    reverse = complement[::-1]
    return reverse
    """ Computes the reverse complementary sequence of DNA for the specfied DNA
        sequence
        dna: a DNA sequence represented as a string
        returns: the reverse complementary DNA sequence represented as a string
    >>> get_reverse_complement("ATGCCCGCTTT")
    'AAAGCGGGCAT'
    >>> get_reverse_complement("CCGCGTTCA")
    'TGAACGCGG'
    """
    # TODO: implement this
    pass

def rest_of_ORF(dna):
   
    """ Takes a DNA sequence that is assumed to begin with a start
        codon and returns the sequence up to but not including the
        first in frame stop codon.  If there is no in frame stop codon,
        returns the whole string.
        dna: a DNA sequence
        returns: the open reading frame represented as a string
    >>> rest_of_ORF("ATGTGAA")
    'ATG'
    >>> rest_of_ORF("ATGAGATAGG")
    'ATGAGA'
    """
   
    Orf1 = []
    for i in range(0, len(dna), 3 ):
        Orf1.append( dna[0+int(i):3+int(i)] )
    #print("ORF IS ---->", Orf1, "\n")
    lis = []
    stops = ['TAG', 'TAA', 'TGA']
    for i in stops:
        if i in Orf1:
            lis.append(i)
            #print("All stop codons in this orf ---->", lis, "\n")
    indexs_stops = []
    for i in lis:
        indexs_stops.append(int(Orf1.index(i)))
    sorted_indexs_stops = sorted(indexs_stops)
    start = (Orf1.index("ATG")) #start index
    #print(start)
    if len(sorted_indexs_stops) >= 1:
        stop = sorted_indexs_stops[0]#stop index
        return ("".join(Orf1[start:stop]))
       
    else:
        return ("".join(Orf1[start:]))

 

    # TODO: implement this
    pass


def find_all_ORFs_oneframe(dna):
    """ Finds all non-nested open reading frames in the given DNA
        sequence and returns them as a list.  This function should
        only find ORFs that are in the default frame of the sequence
        (i.e. they start on indices that are multiples of 3).
        By non-nested we mean that if an ORF occurs entirely within
        another ORF, it should not be included in the returned list of ORFs.
        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGCCC")
    ['ATGCATGAATGTAGA', 'ATGTGCCC']
    """
    listOfOrf = list()
    
    frames = [] # storing the default reading frame
        # create the  frame
        # split the frames into codons for better performance
    frames.append([dna[i:i + 3] for i in range(0, len(dna), 3)])
    #print(frames)
    for i in range(0,len(frames),1): #looping  the dna frame
        start=0
        while start <len(frames[i]): #looping the frame for start and stop codons
            if frames[i][start]=="ATG":
                for stop in range(start+1,len(frames[i]),1):
                    if frames[i][stop]=="TAA" or  frames[i][stop]=="TAG" or  frames[i][stop]=="TGA" :
                        listOfOrf.append(' '.join(frames[i][start:stop])) # retrieve the orf
                        break
                else:
                     listOfOrf.append(' '.join(frames[i][start:]))
            start+=1
    one_f_orf =(",".join(listOfOrf).replace(" ",""))
    one_f_orf = one_f_orf.split(",")
    return one_f_orf

    pass


def find_all_ORFs(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence in
        all 3 possible frames and returns them as a list.  By non-nested we
        mean that if an ORF occurs entirely within another ORF and they are
        both in the same frame, it should not be included in the returned list
        of ORFs.
        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs("ATGCATGAATGTAG")
    ['ATGCATGAATGTAG', 'ATGAATGTAG', 'ATG']
    """
    listOfOrf = list()
    frames = [] # storing the three reading frames 
        # create the positive frames
        # split the frames into codons for better performance
    frames.append([dna[i:i + 3] for i in range(0, len(dna), 3)])
    frames.append([dna[i:i + 3] for i in range(1, len(dna), 3)])
    frames.append([dna[i:i + 3] for i in range(2, len(dna), 3)])
    #print(frames)
    for i in range(0,len(frames),1): #looping all the frames
        start=0
        while start <len(frames[i]): #looping each frame for start and stop codons
            if frames[i][start]=="ATG":
                for stop in range(start+1,len(frames[i]),1):
                    if frames[i][stop]=="TAA" or  frames[i][stop]=="TAG" or  frames[i][stop]=="TGA" :
                        listOfOrf.append(' '.join(frames[i][start:stop])) # retrieve the orf
                        break
                else:
                     listOfOrf.append(' '.join(frames[i][start:]))
            start+=1
    all_orf= ",".join(listOfOrf).replace(" ","")
    all_orf = all_orf.split(",")
    return all_orf

    # TODO: implement this
    pass

def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string
    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    """
   
    The_orf = find_all_ORFs_both_strands(dna)
    thelength=[]
    for i in The_orf:
        thelength.append(len(i))
    seqlen=dict((j,i) for j,i in zip(thelength,The_orf))
    orderedlength=sorted(thelength)
    return seqlen[thelength[-1]]

    # TODO: implement this
    pass



def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.
ATGAGGCTCAGGGATGATCTTGGGTTTTGTAATGGTCGCTGTACGATTATGATCG
​
        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATG', 'ATGCTACATTCGCAT']
    """
    
    listOfOrf = list()
    frames = [] # storing the six frames that would be extacted from the fragments
    reverseCdna = [] # storing the reverse compliments
        # create the foward strand frames
        # split the frames into codons for better performance
    frames.append([dna[i:i + 3] for i in range(0, len(dna), 3)])
    frames.append([dna[i:i + 3] for i in range(1, len(dna), 3)])
    frames.append([dna[i:i + 3] for i in range(2, len(dna), 3)])
    # reverse compliment of the fragment
    reverse = {"A": "T", "C": "G", "T": "A", "G": "C"}
    for i in range(len(dna)):
        reverseCdna.append(reverse[dna[-i - 1]]) if dna[-i - 1] in reverse.keys() else reverseCdna.append(dna[-i - 1])  # if any  contamination found we keep it for further more check
    reverseCdna = ''.join(reverseCdna) # joining
    # create the reverse complement  frames
    frames.append([reverseCdna[i:i + 3] for i in range(0, len(reverseCdna), 3)])
    frames.append([reverseCdna[i:i + 3] for i in range(1, len(reverseCdna), 3)])
    frames.append([reverseCdna[i:i + 3] for i in range(2, len(reverseCdna), 3)])
    #print(frames)
    #print(reverseCdna)
    for i in range(0,len(frames),1): #looping all the frames
        start=0
        while start <len(frames[i]): #looping each frame for start and stop codons
            if frames[i][start]=="ATG":
                for stop in range(start+1,len(frames[i]),1):
                    if frames[i][stop]=="TAA" or  frames[i][stop]=="TAG" or  frames[i][stop]=="TGA" :
                        listOfOrf.append(' '.join(frames[i][start:stop])) # retrieve the orf
                        break
                else:
                     listOfOrf.append(' '.join(frames[i][start:]))
            start+=1
    my_string =",".join(listOfOrf).replace(" ","")
    my_list = my_string.split(",")
    return my_list
    # TODO: implement this
    pass


def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence
        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF 
    """
    import random
    i=0
    frames3 =[]
    listofshuffled=[]
    while i <num_trials:
        listofshuffled.append(shuffle_string(dna))
        i+=1
    for i in listofshuffled:
        frames3.append(longest_ORF(i))
    return (max(frames3,key=len))
    
    # TODO: implement this
    pass


def coding_strand_to_AA(dna):
   
    """ Computes the Protein encoded by a sequence of DNA.  This function
        does not check for start and stop codons (it assumes that the input
        DNA sequence represents an protein coding region).
        dna: a DNA sequence represented as a string
        returns: a string containing the sequence of amino acids encoded by the
                 the input DNA fragment
        >>> coding_strand_to_AA("ATGCGA")
        'MR'
        >>> coding_strand_to_AA("ATGCCCGCTTT")
        'MPA'
    """
    aa = ['F', 'L', 'I', 'M', 'V', 'S', 'P', 'T', 'A', 'Y',
      '|', 'H', 'Q', 'N', 'K', 'D', 'E', 'C', 'W', 'R',
      'G']

    codons = [['TTT', 'TTC'],
              ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'],
              ['ATT', 'ATC', 'ATA'],
              ['ATG'],
              ['GTT', 'GTC', 'GTA', 'GTG'],
              ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'],
              ['CCT', 'CCC', 'CCA', 'CCG'],
              ['ACT', 'ACC', 'ACA', 'ACG'],
              ['GCT', 'GCC', 'GCA', 'GCG'],
              ['TAT', 'TAC'],
              ['TAA', 'TAG', 'TGA'],
              ['CAT', 'CAC'],
              ['CAA', 'CAG'],
              ['AAT', 'AAC'],
              ['AAA', 'AAG'],
              ['GAT', 'GAC'],
              ['GAA', 'GAG'],
              ['TGT', 'TGC'],
              ['TGG'],
              ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
              ['GGT', 'GGC', 'GGA', 'GGG']]

    # create a dictionary lookup table for mapping codons into amino acids
    aa_table = {}
    for i in range(len(aa)):
        for codon in codons[i]:
            aa_table[codon] = aa[i]
    init_pos = 0
    return''.join([aa_table[dna[pos:pos + 3]] for pos in range (init_pos, len(dna) -2, 3)])
    # TODO: implement this
    
    pass


def gene_finder(dna):
    """ Returns the amino acid sequences that are likely coded by the specified dna
        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.
    """
    my_genes=[]
    for i in find_all_ORFs_both_strands(dna):
        my_genes.append(coding_strand_to_AA(i))
    return my_genes

    # TODO: implement this
    pass

if __name__ == "__main__":
    import doctest
    doctest.testmod()
