"""
Console script for RDG.

Allows interactive command line interface to RDG plot generation 
"""
import sys

import click
import os
from rich.pretty import pprint

from sequence_to_RDG import build_graphs_from_fasta
from plot import plot

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
    Minimum length of the Translon to search for in the input sequence. Default: 10
    '''
    )
@click.option(
    '--reinitiation',
    '-r',
    default=True,
    help='''
    Whether to search for reinitiation events in the input sequence. Default: False
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


if __name__ == '__main__':
    rdg_cli()
