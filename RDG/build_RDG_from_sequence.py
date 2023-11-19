from RDG import RDG
from progress.bar import Bar
from RDG.plot_RDG import plot

# from profilestats import profile


# table = """TTT F      CTT L      ATT I      GTT V
# TTC F      CTC L      ATC I      GTC V
# TTA L      CTA L      ATA I      GTA V
# TTG L      CTG L      ATG M      GTG V
# TCT S      CCT P      ACT T      GCT A
# TCC S      CCC P      ACC T      GCC A
# TCA S      CCA P      ACA T      GCA A
# TCG S      CCG P      ACG T      GCG A
# TAT Y      CAT H      AAT N      GAT D
# TAC Y      CAC H      AAC N      GAC D
# TAA Stop   CAA Q      AAA K      GAA E
# TAG Stop   CAG Q      AAG K      GAG E
# TGT C      CGT R      AGT S      GGT G
# TGC C      CGC R      AGC S      GGC G
# TGA Stop   CGA R      AGA R      GGA G
# TGG W      CGG R      AGG R      GGG G"""
# table = dict(zip(table.split()[::2],table.split()[1::2]))


def parse_sequence_into_translated_regions(
    sequence, starts=["ATG", "CTG", "GTG"], min_length=10, reinitiation=False
):
    """
    take a nucleotide sequence and return a list of ORFs
    """
    frames = {1: [], 2: [], 3: []}
    for i in range(len(sequence) + 1):
        try:
            codon = sequence[i : i + 3]
            if len(codon) == 3:
                frames[1].append(codon)
        except KeyError:
            continue
        try:
            codon = sequence[i + 1 : i + 4]
            if len(codon) == 3:
                frames[1].append(codon)
        except KeyError:
            continue
        try:
            codon = sequence[i + 2 : i + 5]
            if len(codon) == 3:
                frames[1].append(codon)
        except KeyError:
            continue

    stops = ["TAA", "TAG", "TGA"]

    orf_sequences = []
    for i in range(len(sequence)):
        if sequence[i : i + 3] in starts:
            for j in range(i, len(sequence), 3):
                if sequence[j : j + 3] in stops:
                    orf = [sequence[k : k + 3] for k in range(i, j + 3, 3)]
                    if len("".join(orf)) > min_length:
                        orf_sequences.append("".join(orf))
                    break
    orfs = []
    counter = 1
    for orf in orf_sequences:
        start_codon_position = sequence.find(orf)
        stop_codon_position = sequence.find(orf) + len(orf)
        orfs.append((start_codon_position, stop_codon_position))
        counter += 1
    return orfs


def build_graphs_from_fasta(
    file_path,
    min_lenth=100,
    start_codons=["ATG", "CTG", "GTG"],
    reinitiation=False,
    readthrough=False,
    upstream_limit=1,
):
    file = open(file_path, "r").readlines()
    sequences = {}
    for line in file:
        if line[0] == ">":
            name = line.split(" ")[0][1:]
            sequences[name] = ""
        else:
            sequences[name] += line.strip("\n")

    graphs = []
    for sequence in sequences:
        orfs = parse_sequence_into_translated_regions(
            sequences[sequence], starts=start_codons, min_length=min_lenth
        )
        dg = RDG(name=sequence, locus_stop=len(sequences[sequence]))
        with Bar("building...") as bar:
            for orf in sorted(orfs):
                dg.add_open_reading_frame(
                    orf[0],
                    orf[1],
                    reinitiation=reinitiation,
                    upstream_limit=upstream_limit,
                )
                bar.next()
        graphs.append(dg)
    return graphs


if __name__ == "__main__":
    graphs = build_graphs_from_fasta(
        "/home/jack/projects/decision_graphs/data/PHPT1_transcript_sequence.fa",
        start_codons=["ATG"],
        min_lenth=10,
        reinitiation=False,
        readthrough=False,
        upstream_limit=100,
    )
    # graphs = build_graphs_from_fasta('/home/jack/projects/decision_graphs/data/test_fasta_multi_sequences.fa')
    for dg in graphs:
        orfs = dg.get_orfs()
        longest_orf = max(orfs, key=lambda x: x[1] - x[0])
        print(longest_orf, longest_orf[0] % 3)
        print(dg.describe())
        with open("test.nwk", "w") as f:
            f.write(dg.newick())
        plot(dg)

