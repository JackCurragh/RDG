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


def find_orfs(sequence, start_codons=None, min_length=30):
    """
    Find all Open Reading Frames (ORFs) in a given transcript sequence.
    
    Parameters:
    -----------
    sequence : str
        The DNA/RNA transcript sequence to search for ORFs.
    start_codons : list
        List of start codons to consider (default: ['ATG', 'GTG', 'TTG']).
    min_length : int
        Minimum length of ORF in nucleotides (default: 30).
        
    Returns:
    --------
    list
        List of tuples containing (start_position, end_position).
    """
    # Default start codons if none provided
    if start_codons is None:
        start_codons = ['ATG', 'GTG', 'TTG']
    
    # Stop codons in DNA
    stop_codons = ['TAA', 'TAG', 'TGA']
    
    # Convert sequence to uppercase
    sequence = sequence.upper()
    
    # Store ORFs as (start, end)
    orfs = []
    
    # Scan through the sequence to find start codons
    seq_len = len(sequence)
    i = 0
    
    while i <= seq_len - 3:
        codon = sequence[i:i+3]
        print(codon, start_codons)
        # Check if current position is a start codon
        if codon in start_codons:
            # Start position of potential ORF
            orf_start = i
            print("found start codon")

            
            # Look for the first in-frame stop codon
            j = i
            found_stop = False
            
            while j <= seq_len - 3:
                current_codon = sequence[j:j+3]
                
                if current_codon in stop_codons:
                    # ORF found with start and stop codons
                    orf_end = j + 2  # End position includes the stop codon
                    
                    # Only include ORFs that meet minimum length requirement
                    if orf_end - orf_start + 1 >= min_length:
                        orfs.append((orf_start, orf_end))
                    
                    found_stop = True
                    break
                
                # Move to the next codon (3 nucleotides)
                j += 3
            
            # If no stop codon found, check if it's still worth including
            if not found_stop and seq_len - orf_start >= min_length:
                # ORF runs until the end of sequence
                orf_end = seq_len - 1
                orfs.append((orf_start, orf_end))
        
        # Move to the next position
        i += 1
    
    return orfs


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

    translons = find_orfs(sequence, starts, min_length)

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
