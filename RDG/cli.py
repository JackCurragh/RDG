"""
Console script for RDG.

Allows interactive command line interface to RDG plot generation 
"""
import click
import pandas as pd

from .sequence_to_RDG import build_graphs_from_fasta, build_graphs_from_bed, build_graphs_from_gtf
from .plot import plot
from .RDG_to_file import save, load, newick_to_file


@click.group()
def rdg_cli():
    pass


@rdg_cli.command()
@click.argument('infile')
@click.option(
    '--start_codons',
    '-s',
    default="ATG,CTG",
    help='Start codons to search for in the input sequence. Default: ATG,CTG'
    )
@click.option(
    '--num_starts',
    '-n',
    default=3,
    help='''
    Number of start codons to search for in the input sequence. Default: 3
    '''
    )
@click.option(
    '--min_length',
    '-m',
    default=10,
    help='''
    Minimum length of the Translon to search for in the input sequence.
    Default: 10
    '''
    )
@click.option(
    '--reinitiation',
    '-r',
    default=True,
    help='''
    Whether to search for reinitiation events in the input sequence.
    Default: False
    ''')
@click.option(
    '--num_visualised',
    '-v',
    default=1,
    help='''
    Number of graphs to visualise. Default: 1
    ''')
def visualise(
        infile,
        start_codons,
        num_starts,
        min_length,
        reinitiation,
        num_visualised
          ):
    graphs = build_graphs_from_fasta(
        infile,
        start_codons=start_codons,
        num_starts=num_starts,
        min_length=min_length,
        reinitiation=reinitiation,
        )
    for graph in graphs[:num_visualised]:
        plot(graph)

@rdg_cli.command()
@click.argument('infile')
@click.option(
    '--input_format',
    '-f',
    type=click.Choice(['fasta', 'gtf', 'bed'], case_sensitive=False),
    required=True,
    help='Input file format (fasta, gtf, bed).'
)
@click.option(
    '--start_codons',
    '-s',
    default="ATG,CTG",
    help='Start codons to search for in the input sequence. Default: ATG,CTG'
    )
@click.option(
    '--num_starts',
    '-n',
    default=3,
    help='''
    Number of start codons to search for in the input sequence. Default: 3
    '''
    )
@click.option(
    '--min_length',
    '-m',
    default=10,
    help='''
    Minimum length of the Translon to search for in the input sequence.
    Default: 10
    '''
    )
@click.option(
    '--reinitiation',
    '-r',
    default=True,
    help='''
    Whether to search for reinitiation events in the input sequence.
    Default: False
    ''')
def construct(
        infile,
        input_format,
        start_codons,
        num_starts,
        min_length,
        reinitiation,
          ):
    if input_format == 'fasta':
        graphs = build_graphs_from_fasta(
            infile,
            start_codons=start_codons,
            num_starts=num_starts,
            min_length=min_length,
            reinitiation=reinitiation,
            )
        for graph in graphs:
            save(graph, f"{graph.name}.sqlite")

    elif input_format == 'gtf':
        graphs = build_graphs_from_gtf(
            infile,
            start_codons=start_codons,
            num_starts=num_starts,
            min_length=min_length,
            reinitiation=reinitiation,
            )
        for graph in graphs:
            save(graph, f"{graph.name}.sqlite")

    elif input_format == 'bed':
        graphs = build_graphs_from_bed(
            infile,
            num_starts=num_starts,
            min_length=min_length,
            reinitiation=reinitiation,
            )
        for graph in graphs:
            save(graph, f"{graph.name}.sqlite")


@rdg_cli.command()
@click.argument('infile')
@click.argument('riboseq_bedgraph')
@click.option(
    '--min_read_support',
    '-m',
    default=10,
    help='Minimum read support for a node to be kept. Default: 10'
    )
def prune(
        infile,
        riboseq_bedgraph,
        min_read_support,
        ):
    '''
    prune graphs based on riboseq data
    '''
    graphs = load(infile)
    out_graphs = []

    riboseq_df = pd.read_csv(riboseq_bedgraph, sep='\t', header=None)
    riboseq_df.columns = ['transcript', 'start', 'end', 'count']
    for graph in graphs:
        riboseq_counts = {} 
        tx_df = riboseq_df[riboseq_df['transcript'] == graph.name]
        for node in graph.nodes:
            riboseq_counts[node] = tx_df[
                (tx_df['start'] <= node.start) & (tx_df['end'] >= node.end)
                ]['count'].sum()

        for node in graph.nodes:
            if riboseq_counts[node] < min_read_support:
                graph.prune(node)
        out_graphs.append(graph)

    for graph in out_graphs:
        save(graph, f"{graph.name}.sqlite")


if __name__ == '__main__':
    rdg_cli()
