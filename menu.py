'''the  menu program that contains the main menu function for calling functions in gene_finder.py'''

import sys
import gene_finder
from gene_finder import get_complement
from load import load_seq
from gene_finder import get_reverse_complement
from gene_finder import rest_of_ORF
from gene_finder import find_all_ORFs_oneframe
from gene_finder import find_all_ORFs
from gene_finder import find_all_ORFs
from gene_finder import find_all_ORFs_both_strands
from gene_finder import longest_ORF
from gene_finder import longest_ORF_noncoding
from gene_finder import coding_strand_to_AA
from gene_finder import gene_finder
def MainMenu():
    option = True
    while option:
        print('''
                         a. Get the complementary nucleotide
                         b. Get the reverse complementary sequence of DNA
                         c. Find the rest_of_ORF
                         d. Find all non-nested open reading frames
                         e. Find all non-nested open reading frames in all 3 possible frames
                         f. Find all non-nested open reading frames on both strands
                         g. Find the longest ORF on both strands
                         h. Compute the maximum length of the longest non-coding ORF
                         i. Compute the Protein encoded by a sequence of DNA
                         j. Find the amino acid sequences that are likely coded by the specified dna
                         Q. Exit:
                         ''')

        dna = load_seq('./data/X73525.fa')
        
        nucleotide = load_seq('./data/X73525.fa')
        num_trials = 15
        option = input("Please enter an option from the above menu ")
        if option == "A" or option == "a":
            print(get_complement(nucleotide))
        if option == "B" or option == "b":
            print(get_reverse_complement(dna))
        elif option == "C" or option == "c":
            print(rest_of_ORF(dna))
        elif option == "D" or option == "d":
            print(find_all_ORFs_oneframe(dna))
        elif option == "E" or option == "e":
            print(find_all_ORFs(dna))
        elif option == "F" or option == "f":
            print(find_all_ORFs_both_strands(dna))
        elif option == "G" or option == "g":
            print(longest_ORF(dna))
        elif option == "H" or option == "h":  
            print(longest_ORF_noncoding(dna,num_trials))
        elif option == "i" or option == "":
            print(coding_strand_to_AA(dna))
        elif option == "J" or option == "j":
                    print(gene_finder(dna))  
        elif option == "Q" or option =="q":
            print('Goodbye')
            option = None
            sys.exit
        else:
            print("You must input an option")
            print("Please try again")
    return MainMenu
MainMenu()
