from pathlib import Path
from typing import List, Tuple, Set
import ahocorasick
from collections import defaultdict
from RDG import RDG
from rich.progress import Progress
import pandas as pd
import gff2bed
from tempfile import TemporaryDirectory
import tempfile
import gzip
import os


def is_gzipped(file_path: str) -> bool:
    """
    Checks whether the file is gzipped or not

    Inputs:
        file_path: Path to the file to be checked

    Outputs:
        True if gzipped, otherwise False
    """
    try:
        with open(file_path, 'rb') as f:
            # Read the first two bytes of the file
            header = f.read(2)

        # Check if the file starts with the gzip magic number (0x1f 0x8b)
        return header == b'\x1f\x8b'

    except IOError:
        # File not found or unable to open
        return False


def extract_transcript_id(attr_str):
    for attr in attr_str.split(";"):
        # Ensembl GFF3 support
        if attr.startswith("Parent=transcript:") \
                or attr.startswith("ID=transcript:"):
            return attr.split(":")[1]
        # Gencode GFF3 support
        elif attr.startswith("transcript_id="):
            return attr.split("=")[1]
        # Ensembl GTF support
        elif attr.startswith(" transcript_id "):
            return attr.split(" ")[2].replace('"', "")
    return "NA"


def convert_to_transcript_coordinates(gtf_df):
    # Group the dataframe by transcript_id
    grouped = gtf_df.groupby('transcript_id')
    converted_dfs = []

    for transcript_id, transcript_df in grouped:
        # Sort the transcript dataframe by start position
        transcript_df = transcript_df.sort_values('start')
        strand = transcript_df['strand'].iloc[0]

        transcript_coord = 0
        
        new_df = transcript_df.copy()
        for idx, row in transcript_df.iterrows():
            feature_length = row['end'] - row['start'] + 1
            if strand == '+':
                new_df.at[idx, 'start'] = transcript_coord
                new_df.at[idx, 'end'] = transcript_coord + feature_length - 1
            else:  # strand == '-'
                new_df.at[idx, 'start'] = transcript_coord
                new_df.at[idx, 'end'] = transcript_coord + feature_length - 1

            transcript_coord += feature_length

        converted_dfs.append(new_df)

    result_df = pd.concat(converted_dfs)

    return result_df


def parse_gff(gff_path: str, num_transcripts: int) -> pd.DataFrame:
    """
    Read in the gff file at the provided path and return a dataframe

    Inputs:
        gff_path: Path to the gff file

    Outputs:
        gff_df: Dataframe containing the gff information
    """
    if is_gzipped(gff_path):
        with gzip.open(gff_path, 'rt') as f:
            bed_df = gff2bed.convert(
                gff2bed.parse(gff_path)
            )  
    else:
        bed_df = gff2bed.convert(
            gff2bed.parse(gff_path)
        )

    return bed_df


def build_codon_automaton(codons):
    """
    Build an Aho-Corasick automaton for codon pattern matching.

    Args:
        codons (list): List of codon sequences to match

    Returns:
        ahocorasick.Automaton: Compiled automaton for pattern matching
    """
    automaton = ahocorasick.Automaton()
    
    for codon in codons:
        # Store the codon as the value for each pattern
        automaton.add_word(codon, codon)
    
    automaton.make_automaton()
    return automaton


def find_all_positions(sequence, automaton):
    """
    Find positions of all occurrences of patterns from an Aho-Corasick automaton.

    Args:
        sequence (str): Input sequence to search for patterns
        automaton (ahocorasick.Automaton): Aho-Corasick automaton with patterns

    Returns:
        tuple: (frames, codons) where:
            - frames (dict): Maps frame indices (0,1,2) to lists of positions
            - codons (dict): Maps positions to identified codons
    """
    frames = {0: [], 1: [], 2: []}
    codons = {}
    
    # Use automaton.iter() to find all matches
    for end_pos, pattern in automaton.iter(sequence):
        # Calculate the start position (Aho-Corasick gives end position)
        start_pos = end_pos - len(pattern) + 1
        
        # Only accept 3-nucleotide patterns (proper codons)
        if len(pattern) == 3:
            # Determine which frame this codon belongs to
            frame = start_pos % 3
            frames[frame].append(start_pos)
            codons[start_pos] = pattern
            
    return frames, codons


def find_orfs(sequence, startautomaton, stopautomaton, minlength=0, maxlength=1000000):
    """
    Predict Open Reading Frames (ORFs) in a nucleotide sequence.

    Args:
        sequence (str): Nucleotide sequence to analyze
        startautomaton (ahocorasick.Automaton): Automaton for start codons
        stopautomaton (ahocorasick.Automaton): Automaton for stop codons
        minlength (int, optional): Minimum ORF length. Defaults to 0
        maxlength (int, optional): Maximum ORF length. Defaults to 1000000

    Returns:
        list: List of tuples (start, stop) representing ORF positions
    """
    
    orf_list = []
    startpositions, start_codons = find_all_positions(sequence, startautomaton)
    stoppositions, stop_codons = find_all_positions(sequence, stopautomaton)
    
    print("Start positions by frame:", {f: len(p) for f, p in startpositions.items()})
    print("Stop positions by frame:", {f: len(p) for f, p in stoppositions.items()})
    
    for frame, start_pos_list in startpositions.items():
        for start_pos in start_pos_list:
            # Find valid stop codons downstream of start position in the same frame
            valid_stops = [pos for pos in stoppositions[frame] if pos > start_pos]
            
            if valid_stops:
                stop_pos = min(valid_stops)
                # Calculate real positions - start is the first nucleotide, stop is the last
                orf_length = stop_pos - start_pos + 3
                if minlength <= orf_length <= maxlength:
                    orf_list.append((start_pos, stop_pos + 2))  # End position is inclusive of stop codon
            else:
                # If no stop codon, go to the end of the sequence
                remaining_length = ((len(sequence) - start_pos) // 3) * 3
                if minlength <= remaining_length <= maxlength:
                    orf_list.append((start_pos, start_pos + remaining_length - 1))

    return orf_list


def extract_translons(
    sequence: str,
    starts: Set[str] = {"ATG", "CTG", "GTG"},
    min_length: int = 10,
) -> List[Tuple[int, int]]:
    """
    Extract open reading frames (translons) from a nucleotide sequence using Aho-Corasick algorithm.

    Parameters:
    - sequence (str): The input nucleotide sequence.
    - starts (Set[str]): Set of start codons to initiate translon detection.
            Default: {"ATG", "CTG", "GTG"}.
    - min_length (int): Minimum length of translons to be included in
            the result. Default: 10.

    Returns:
    List[Tuple[int, int]]: A list of tuples representing the start and
                      stop positions of detected translons.

    Example:
    ```python
    sequence = "ATGCTAGCATGAATAG"
    translons = extract_translons(sequence)
    print(translons)
    # Output: [(0, 15)]
    ```
    """
    if not sequence:
        return []

    stops = {"TAA", "TAG", "TGA"}
    startautomaton = build_codon_automaton(list(starts))
    stopautomaton = build_codon_automaton(list(stops))

    translons = find_orfs(sequence, startautomaton, stopautomaton, minlength=min_length)

    return sorted([t for t in translons if t[1] - t[0] >= min_length])


def build_graphs_from_fasta(
    file_path: str,
    min_length: int = 100,
    start_codons: Set[str] = {"ATG", "CTG", "GTG"},
    reinitiation: bool = False,
    upstream_limit: int = 1,
    num_starts: int = 3,
) -> List[RDG]:
    """
    Build graphs from a FASTA file.

    Parameters:
    - file_path (str): The path to the FASTA file.
    - min_length (int): Minimum length of translons to be included in the
                result. Default: 100.
    - start_codons (Set[str]): Set of start codons to initiate translon
                detection. Default: {"ATG", "CTG", "GTG"}.
    - reinitiation (bool): Flag indicating whether reinitiation is allowed.
                Default: False.
    - readthrough (bool): Flag indicating whether readthrough is allowed.
                Default: False.
    - upstream_limit (int): Limit for upstream distance in the graph.
                Default: 1.
    - num_starts (int): Number of start codons to search for in the input

    Returns:
    List[RDG]: A list of RDG (RDG) objects representing the graphs.

    Example:
    ```python
    file_path = "example.fasta"
    graphs = build_graphs_from_fasta(file_path)
    for graph in graphs:
        print(f"Graph name: {graph.name}, Locus stop: {graph.locus_stop}")
        for translon in graph.open_reading_frames:
            print(
            f"translon start: {translon.start}, translon stop: {translon.stop}"
                )
    ```
    """
    file_path = Path(file_path)
    sequence_dict = defaultdict(str)

    with file_path.open("r") as file:
        current_name = None
        for line in file:
            if line.startswith(">"):
                current_name = line.split(" ")[0][1:]
            else:
                sequence_dict[current_name] += line.strip()

    result_graphs = []
    with Progress() as progress:
        task = progress.add_task(
            "[cyan]Building graphs...",
            total=min(num_starts, len(sequence_dict)),
        )
        for sequence_name, sequence in sequence_dict.items():
            translons = extract_translons(
                sequence, starts=start_codons, min_length=min_length
            )
            dg = RDG(name=sequence_name, locus_stop=len(sequence))

            for start, stop in sorted(translons)[:num_starts]:
                dg.add_open_reading_frame(
                    start,
                    stop,
                    reinitiation=reinitiation,
                    upstream_limit=upstream_limit,
                )
            progress.update(task, advance=1)
            result_graphs.append(dg)

    return result_graphs


def build_graphs_from_gtf(file_path: str,
    min_length: int = 100,
    reinitiation: bool = False,
    upstream_limit: int = 1,
    num_starts: int = 3,
) -> List[RDG]:
    """
    Build graphs from a GTF file.

    Parameters:
    - file_path (str): The path to the GTF file.

    Returns:
    List[RDG]: A list of RDG (RDG) objects representing the graphs.

    """
    df = parse_gff(file_path)
    transcript_ids = df.iloc[:, 8].apply(
        lambda x: x.split('transcript_id "')[1].split('"')[0])
    start_coords = df.iloc[:, 3]
    stop_coords = df.iloc[:, 4]

    transcripts_df = pd.DataFrame({
        'Start': start_coords,
        'Stop': stop_coords
    })

    result_graphs = []
    for transcript in transcripts_df['Transcript_ID'].unique():
        transcript_df = transcripts_df
        transcript_df = transcript_df.sort_values(by='Start')
        transcript_df = transcript_df.assign(
                            Length=transcript_df['Stop'] - transcript_df['Start']
                            )   
        transcript_df = transcript_df[transcript_df['Length'] >= min_length]

        graph = RDG(transcript)
        for index, row in transcript_df.iterrows()[:num_starts]:
            graph.add_open_reading_frame(
                row['Start'],
                row['Stop'],
                reinitiation=reinitiation,
                upstream_limit=upstream_limit,
            )
        result_graphs.append(graph)

    return result_graphs


def build_graphs_from_bed(file_path: str,
                          min_length: int = 100,
                          reinitiation: bool = False,
                          upstream_limit: int = 1,
                          num_starts: int = 3,
                        ) -> List[RDG]:
    """
    Build graphs from a BED file.

    Parameters:
    - file_path (str): The path to the BED file.

    Returns:
    List[RDG]: A list of RDG (RDG) objects representing the graphs.

    """
    df = pd.read_csv(file_path, sep='\t', header=None)
    transcript_ids = df.iloc[:, 21]
    start_coords = df.iloc[:, 10].apply(
        lambda x: min([int(x) for x in x.split(',')]))
    stop_coords = df.iloc[:, 11].apply(
        lambda x: max([int(x) for x in x.split(',')]))

    transcripts_df = pd.DataFrame({
        'Start': start_coords,
        'Stop': stop_coords,
        'Transcript_ID': transcript_ids
    })

    result_graphs = []
    for transcript in transcripts_df['Transcript_ID'].unique():
        transcript_df = transcripts_df[transcripts_df['Transcript_ID'] == transcript]
        transcript_df = transcript_df.sort_values(by='Start')
        transcript_df = transcript_df.assign(
                            Length=transcript_df['Stop'] - transcript_df['Start']
                            )
        transcript_df = transcript_df[transcript_df['Length'] >= min_length]

        graph = RDG(transcript)
        for index, row in transcript_df.iterrows()[:num_starts]:
            graph.add_open_reading_frame(
                row['Start'],
                row['Stop'],
                reinitiation=reinitiation,
                upstream_limit=upstream_limit,
            )
        result_graphs.append(graph)

    return result_graphs
