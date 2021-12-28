 
# The Gene_finder program

### This program is designed to manipulate data from DNA and find the genes that are encoded by the DNA sequence. It contains 11 independent functions.These functions are :

* get_complement(nucleotide)
    * Returns the complementary nucleotide
        nucleotide: a nucleotide (A, C, G, or T) represented as a string
        returns: the complementary nucleotide
    `get_complement('A')
    'T'`
* get_reverse_complement(dna)
    * Computes the reverse complementary sequence of DNA for the specfied DNA
        sequence
* rest_of_ORF(dna)
    * Takes a DNA sequence that is assumed to begin with a start
        codon and returns the sequence up to but not including the
        first in frame stop codon.  If there is no in frame stop codon,
        returns the whole string.
* find_all_ORFs_oneframe(dna)
    * Finds all non-nested open reading frames in the given DNA
        sequence and returns them as a list.  This function should
        only find ORFs that are in the default frame of the sequence
        (i.e. they start on indices that are multiples of 3).
        By non-nested we mean that if an ORF occurs entirely within
        another ORF, it should not be included in the returned list of ORFs.
* find_all_ORFs(dna)
    * Finds all non-nested open reading frames in the given DNA sequence in
        all 3 possible frames and returns them as a list.  By non-nested we
        mean that if an ORF occurs entirely within another ORF and they are
        both in the same frame, it should not be included in the returned list
        of ORFs.
* longest_ORF(dna)
    * Finds the longest ORF on both strands of the specified DNA and returns it
        as a string
*  find_all_ORFs_both_strands(dna)
    * Finds all non-nested open reading frames in the given DNA sequence on both
        strands
* shuffle_string(s)
    * Shuffles the characters in the input string

* longest_ORF_noncoding(dna, num_trials)
    * Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence
* coding_strand_to_AA(dna)
    * Computes the Protein encoded by a sequence of DNA.  This function
        does not check for start and stop codons (it assumes that the input
        DNA sequence represents an protein coding region)
* gene_finder(dna)
    *  Returns the amino acid sequences that are likely coded by the specified dna


To start the program import the `MainMenu()` function from the menu module:
### disclaimer!!

in the menu module, MainMenu function ... you need to ensure the path matches your path. The path given here matches the data give
