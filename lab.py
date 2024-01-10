from scipy.spatial.distance import euclidean
import numpy as np

def parse_fasta(filename):
    with open(filename, "r") as f:
        sequences = {}
        header = None
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                header = line[1:].split()[0]  # Take only the first word of the header
                sequences[header] = ""
            else:
                sequences[header] += line
        return sequences


def find_start_stop_pairs(sequence):
    start_codon = "ATG"
    stop_codons = {"TAA", "TAG", "TGA"}

    # Translating to reverse complement
    trans = str.maketrans("ATCG", "TAGC")
    rev_comp = sequence[::-1].translate(trans)

    pairs = []

    for strand_seq in [sequence, rev_comp]:
        start_indices = []
        stop_indices = []

        # First pass: Record positions of all start and stop codons
        for i in range(0, len(strand_seq) - 2, 3):
            codon = strand_seq[i:i+3]
            if codon == start_codon:
                start_indices.append(i)
            elif codon in stop_codons:
                stop_indices.append(i)

        # Second pass: Pair start and nearest stop codons
        for start in start_indices:
            for stop in stop_indices:
                if stop > start:
                    pairs.append((start, stop + 3))  # Adding 3 to include the stop codon
                    break

    return pairs


def filter_pairs_by_length(sequence, pairs):
    return [(start, end) for start, end in pairs if end - start+3>= 100]


def translate_dna_to_protein(dna_sequence):
    genetic_code = {
        "ATA":"I", "ATC":"I", "ATT":"I", "ATG":"M",
        "ACA":"T", "ACC":"T", "ACG":"T", "ACT":"T",
        "AAC":"N", "AAT":"N", "AAA":"K", "AAG":"K",
        "AGC":"S", "AGT":"S", "AGA":"R", "AGG":"R",                 
        "CTA":"L", "CTC":"L", "CTG":"L", "CTT":"L",
        "CCA":"P", "CCC":"P", "CCG":"P", "CCT":"P",
        "CAC":"H", "CAT":"H", "CAA":"Q", "CAG":"Q",
        "CGA":"R", "CGC":"R", "CGG":"R", "CGT":"R",
        "GTA":"V", "GTC":"V", "GTG":"V", "GTT":"V",
        "GCA":"A", "GCC":"A", "GCG":"A", "GCT":"A",
        "GAC":"D", "GAT":"D", "GAA":"E", "GAG":"E",
        "GGA":"G", "GGC":"G", "GGG":"G", "GGT":"G",
        "TCA":"S", "TCC":"S", "TCG":"S", "TCT":"S",
        "TTC":"F", "TTT":"F", "TTA":"L", "TTG":"L",
        "TAC":"Y", "TAT":"Y", "TAA":"_", "TAG":"_",
        "TGC":"C", "TGT":"C", "TGA":"_", "TGG":"W",
    }

    protein = ""
    for i in range(0, len(dna_sequence) - 2, 3):
        codon = dna_sequence[i:i+3]
        protein += genetic_code.get(codon, "?")

    return protein


def update_amino_acid_frequencies(protein_sequence, amino_acid_freq):
    total_amino_acids = len(protein_sequence)
    for amino_acid in protein_sequence:
        if amino_acid in amino_acid_freq:
            amino_acid_freq[amino_acid] += 1 / total_amino_acids

def update_dicodon_frequencies(dna_sequence, dicodon_freq):
    total_dicodons = len(dna_sequence) - 5
    for i in range(total_dicodons):
        dicodon = dna_sequence[i:i+6]
        if dicodon in dicodon_freq:
            dicodon_freq[dicodon] += 1 / total_dicodons


def initialize_amino_acid_frequencies():
    amino_acids = ['A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
    amino_acid_freq = {aa: 0 for aa in amino_acids}
    return amino_acid_freq

def initialize_diamino_acid_frequencies():
    amino_acids = ['A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
    diamino_acid_freq = {aa1+aa2: 0 for aa1 in amino_acids for aa2 in amino_acids}
    return diamino_acid_freq


def calculate_euclidean_distance(vec1, vec2):
    return np.sqrt(np.sum((vec1 - vec2) ** 2))

def frequency_dict_to_vector(frequency_dict):
    if not isinstance(frequency_dict, dict):
        raise TypeError("Expected a dictionary for frequency data, got {}".format(type(frequency_dict)))

    return np.array(list(frequency_dict.values()))

def create_distance_matrix(frequency_data):
    num_sequences = len(frequency_data)
    distance_matrix = np.zeros((num_sequences, num_sequences))

    for i in range(num_sequences):
        for j in range(i + 1, num_sequences):  # Matrix is symmetric
            distance = calculate_euclidean_distance(frequency_data[i], frequency_data[j])
            distance_matrix[i][j] = distance
            distance_matrix[j][i] = distance

    return distance_matrix


def format_distance_matrix_phylip(distance_matrix, sequence_names):

    num_sequences = len(sequence_names)
    phylip_formatted = f"{num_sequences}\n"

    for i, name in enumerate(sequence_names):
        distances = " ".join(f"{distance_matrix[i][j]:.3f}" for j in range(num_sequences))
        phylip_formatted += f"{name} {distances}\n"

    return phylip_formatted


def save_to_file(output, filename):
    with open(filename, 'w') as f:
        f.write(output)


def main(input_filenames, output_filename):
    master_amino_acid_freq = initialize_amino_acid_frequencies()
    master_diamino_acid_freq = initialize_diamino_acid_frequencies()

    all_amino_acid_freq = {}
    all_diamino_acid_freq = {}
    sequence_names = []

    for filename in input_filenames:
        sequences = parse_fasta(filename)

        for header, sequence in sequences.items():
            
            if header not in all_amino_acid_freq:
                all_amino_acid_freq[header] = master_amino_acid_freq.copy()
                all_diamino_acid_freq[header] = master_diamino_acid_freq.copy()

            pairs = find_start_stop_pairs(sequence)
            filtered_pairs = filter_pairs_by_length(sequence, pairs)

        for start, end in filtered_pairs:
            dna_segment = sequence[start:end]

            
            protein_sequence = translate_dna_to_protein(dna_segment)
            update_amino_acid_frequencies(protein_sequence, all_amino_acid_freq[header])
            update_dicodon_frequencies(dna_segment, all_diamino_acid_freq[header])

            if header not in sequence_names:
                sequence_names.append(header)


    frequency_data = all_amino_acid_freq  

    frequency_vectors = [frequency_dict_to_vector(freq_dict) for freq_dict in frequency_data.values()]
    distance_matrix = create_distance_matrix(frequency_vectors)
    phylip_output = format_distance_matrix_phylip(distance_matrix, sequence_names)
    save_to_file(phylip_output, output_filename)

input_filenames = ["mamalian1.fasta", "bacterial1.fasta", "mamalian2.fasta", "bacterial2.fasta", "mamalian3.fasta", "bacterial3.fasta", "mamalian4.fasta", "bacterial4.fasta"]
output_filename = "output123.txt"

main(input_filenames, output_filename)

