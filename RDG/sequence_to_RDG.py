from pathlib import Path
from typing import List, Tuple, Set
from collections import defaultdict
from RDG import RDG
from rich.progress import Progress


def extract_translons(
    sequence: str,
    starts: Set[str] = {"ATG", "CTG", "GTG"},
    min_length: int = 10,
) -> List[Tuple[int, int]]:
    """
    Extract open reading frames (translons) from a nucleotide sequence.

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
    """
    if not sequence:
        return []

    stops = {"TAA", "TAG", "TGA"}
    translons = []

    for i, _ in enumerate(sequence):
        if sequence[i: i + 3] in starts:
            codon_list = [
                sequence[k: k + 3] for k in range(i, len(sequence), 3)
                ]
            for j, stop_codon in enumerate(codon_list):
                if stop_codon in stops:
                    translon = "".join(codon_list[: j + 1])
                    if len(translon) > min_length:
                        start_codon_position = i
                        stop_codon_position = i + len(translon)
                        translons.append(
                            (start_codon_position, stop_codon_position)
                            )
                    break

    return translons


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
            total=num_starts
        )
        for sequence_name, sequence in sequence_dict.items():
            translons = extract_translons(
                sequence, starts=start_codons, min_length=min_length
            )
            dg = RDG(name=sequence_name, locus_stop=len(sequence))

            for translon_start, translon_stop in sorted(translons)[:num_starts]:
                dg.add_open_reading_frame(
                    translon_start,
                    translon_stop,
                    reinitiation=reinitiation,
                    upstream_limit=upstream_limit,
                )
            progress.update(task, advance=1)
            result_graphs.append(dg)

    return result_graphs
