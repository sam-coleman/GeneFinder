# -*- coding: utf-8 -*-
"""
Gene Finder Mini Project 1

@author: scoleman
"""

import random
from amino_acids import aa, codons, aa_table   # you may find these useful
from load import load_seq
dna = load_seq("./data/X73525.fa")

def shuffle_string(s):
    """Shuffles the characters in the input string
        NOTE: this is a helper function, you do not
        have to modify this in any way """
    return ''.join(random.sample(s, len(s)))

# YOU WILL START YOUR IMPLEMENTATION FROM HERE DOWN ###


def get_complement(nucleotide):
    """ Returns the complementary nucleotide

        nucleotide: a nucleotide (A, C, G, or T) represented as a string
        returns: the complementary nucleotide
    Added unit test for if there is no A T C G entered
    >>> get_complement('A')
    'T'
    >>> get_complement('C')
    'G'
    >>> get_complement('')
    ''
    """
    # TODO: implement this
    if nucleotide == 'A':
        return 'T'
    elif nucleotide == 'T':
        return 'A'
    elif nucleotide == 'C':
        return 'G'
    elif nucleotide == 'G':
        return 'C'
    else:
        return ''

def get_reverse_complement(dna):
    """ Computes the reverse complementary sequence of DNA for the specfied DNA
        sequence

        dna: a DNA sequence represented as a string
        returns: the reverse complementary DNA sequence represented as a string

        The unit tests covers the cases
    >>> get_reverse_complement("ATGCCCGCTTT")
    'AAAGCGGGCAT'
    >>> get_reverse_complement("CCGCGTTCA")
    'TGAACGCGG'
    """
    # TODO: implement this
    reverse_comp = ''
    reverse = ''
    index = len(dna) - 1

    #reverse order of dna sequence
    while index >= 0:
        reverse += dna[index]
        index -= 1

    #get completment for each nucleotide in reverse
    for letter in reverse:
        reverse_comp += get_complement(letter)

    return reverse_comp

def rest_of_ORF(dna):
    """ Takes a DNA sequence that is assumed to begin with a start
        codon and returns the sequence up to but not including the
        first in frame stop codon.  If there is no in frame stop codon,
        returns the whole string.

        dna: a DNA sequence
        returns: the open reading frame represented as a string
        Added unit test for when the dna sequence does not end with a stop codon (it returns the whole string)

    >>> rest_of_ORF("ATGTGAA")
    'ATG'
    >>> rest_of_ORF("ATGAGATAGG")
    'ATGAGA'
    >>> rest_of_ORF('ATGAGATCCC')
    'ATGAGATCCC'
    """
    # TODO: implement this
    orf = ''
    i = 0
    while i < len(dna):
        current_letter = dna[i]
        if dna[i:i+3] == 'TAG' or dna[i:i+3] == 'TAA' or dna[i:i+3] == 'TGA':
            return orf
        orf += dna[i:i+3]
        i+=3

    return orf

def find_all_ORFs_oneframe(dna):
    """ Finds all non-nested open reading frames in the given DNA
        sequence and returns them as a list.  This function should
        only find ORFs that are in the default frame of the sequence
        (i.e. they start on indices that are multiples of 3).
        By non-nested we mean that if an ORF occurs entirely within
        another ORF, it should not be included in the returned list of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs

        added unit test for if dna doesn't begin with start codon
    >>> find_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGCCC")
    ['ATGCATGAATGTAGA', 'ATGTGCCC']
    >>> find_all_ORFs_oneframe("ATGCGTAATGCCT")
    ['ATGCGTAATGCCT']
    >>> find_all_ORFs_oneframe("CTTATGCGTAATGCCT")
    ['ATGCGTAATGCCT']
    """
    # TODO: implement this

    orfs = []
    i = 0
    while i <  len(dna):
        if dna[i:i+3] == 'ATG':
            orf = rest_of_ORF(dna[i:])
            orfs.append(orf)
            i += len(orf) + 3
        else:
            i += 3
    return orfs

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
    # TODO: implement this
    all_orfs = []
    for i in range(0,3):
        all_orfs += find_all_ORFs_oneframe(dna[i:])
        i += 1

    return all_orfs

def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCAT']
    """
    # TODO: implement this
    reverse_dna = get_reverse_complement(dna)
    orfs_both_strands = []
    orfs_both_strands += find_all_ORFs(dna)
    orfs_both_strands += find_all_ORFs(reverse_dna)
    return orfs_both_strands

def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string
        Added unit test is for when there are multiple orfs to iterate through
    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    >>> longest_ORF("ATGCATGAATGTAG")
    'ATGCATGAATGTAG'
    """
    # TODO: implement this
    pass
    all_orfs = find_all_ORFs_both_strands(dna)
    longest_orf = all_orfs[0]
    for item in all_orfs:
        if len(item) > len(longest_orf):
            longest_orf = item

    return longest_orf

def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence

        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF
        To test, I entered a short string of dna and a small number of trials. I added the lengths of each orf into a list and printed the list.
        I then ensured the longest length is what is returned. Once I tested this a few times, I commented out the lines I used for testing.
        """
    # TODO: implement this

    maxOrfLength = 0
    #allLengths = [] #for testing
    for i in range(num_trials):
        sequence = shuffle_string(dna)
        orf = longest_ORF(sequence)
        #allLengths.append(len(orf)) #for testing
        if len(orf) > maxOrfLength:
            maxOrfLength = len(orf)

    #print(allLengths) #for testing
    #print(maxOrfLength)
    return maxOrfLength

def coding_strand_to_AA(dna):
    """ Computes the Protein encoded by a sequence of DNA.  This function
        does not check for start and stop codons (it assumes that the input
        DNA sequence represents an protein coding region).

        dna: a DNA sequence represented as a string
        returns: a string containing the sequence of amino acids encoded by the
                 the input DNA fragment
        These unit tests cover what is needed

        >>> coding_strand_to_AA("ATGCGA")
        'MR'
        >>> coding_strand_to_AA("ATGCCCGCTTT")
        'MPA'
    """
    # TODO: implement this
    amino_acid = ''
    i = 0
    while i <=  len(dna) - 3:
        amino_acid += aa_table[dna[i:i+3]]
        i += 3

    return amino_acid

def gene_finder(dna):
    """ Returns the amino acid sequences that are likely coded by the specified dna

        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.
    """
    # TODO: implement this
    threshold = longest_ORF_noncoding(dna, 1500)
    all_orf_both = find_all_ORFs_both_strands(dna)
    i = 0
    aa = []

    #have only orfs longer than threshold
    while i < len(all_orf_both):
        if len(all_orf_both[i]) < threshold:
            all_orf_both.pop(i)
        else:
            i += 1

    #get list of amino amino_acids
    for item in all_orf_both:
        aa.append(coding_strand_to_AA(item))

    print(aa)
    #make it easier to copy and paste aa into Blast
    for item in aa:
        print(item)

gene_finder(dna)

if __name__ == "__main__":
    import doctest
    #doctest.testmod()
    doctest.run_docstring_examples(find_all_ORFs_oneframe, globals(), verbose = False)
