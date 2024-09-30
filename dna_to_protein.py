from Bio.Seq import Seq
from Bio.Data import CodonTable

# Amino Dictionary
amino_acid_names = {
    'A': 'Ala',  
    'R': 'Arg', 
    'N': 'Asn', 
    'D': 'Asp',  
    'C': 'Cys',  
    'E': 'Glu', 
    'Q': 'Gln',  
    'G': 'Gly',  
    'H': 'His',  
    'I': 'Ile',  
    'L': 'Leu',  
    'K': 'Lys',  
    'M': 'Met',  
    'F': 'Phe', 
    'P': 'Pro', 
    'S': 'Ser',  
    'T': 'Thr', 
    'W': 'Trp',  
    'Y': 'Tyr',  
    'V': 'Val',  
    '*': 'Stop' 
}

# DNA to Protein Converter
def dna_to_protein(dna_sequence):
    dna_sequence = Seq(dna_sequence)
    complement_dna = dna_sequence.complement()

    mrna_sequence = complement_dna.transcribe()
   
    protein_sequence = mrna_sequence.translate(to_stop=True)

    formatted_protein = ' - '.join(f"{amino_acid_names[aa]} ({aa})" for aa in protein_sequence)
    
    return complement_dna, mrna_sequence, formatted_protein
dna_sequence = input("Enter DNA sequence: ")

complement_dna, mrna, formatted_protein = dna_to_protein(dna_sequence)

print(f"Complementary DNA: {complement_dna}")
print(f"mRNA: {mrna}")
print(f"Amino acids: {formatted_protein}")
