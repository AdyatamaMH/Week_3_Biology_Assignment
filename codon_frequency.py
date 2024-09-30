from Bio.Data import CodonTable
from itertools import product

# Codon
codon_table = CodonTable.unambiguous_rna_by_name["Standard"]

def codon_frequency(amino_acids):
    codon_frequencies = []

    for amino_acid in amino_acids:
        possible_codons = [codon for codon, amino in codon_table.forward_table.items() if amino == amino_acid]
        
        if possible_codons:
            codon_frequencies.append(possible_codons)

    all_combinations = product(*codon_frequencies)
    
    mrna_sequences = [''.join(combination) for combination in all_combinations]
    
    return mrna_sequences

def count_codons_per_sequence(mrna_sequence, amino_acids):
    codon_count = {}
    for amino_acid in amino_acids:
        codons = [codon for codon, amino in codon_table.forward_table.items() if amino == amino_acid]
        for codon in codons:
            codon_count[codon] = mrna_sequence.count(codon)
    return codon_count

amino_acids = input("Enter amino acid sequence (e.g., W Y N, Maximum 3 amino acids): ").split()

mrna_sequences = codon_frequency(amino_acids)

# Display results
print(f"\nConstructed mRNA Sequences (9 nucleotides each) and their codon counts:")


for sequence in mrna_sequences:
    codon_counts = count_codons_per_sequence(sequence, amino_acids)
    

    print(f"\n{sequence}")
    

    for codon, count in codon_counts.items():
        print(f"{codon}: {count}")
