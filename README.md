# RDG
Implementation of the Ribosome Decision Graph Concept

## Background 
  <kbd>![Depiction of typical annotation structure ](https://github.com/JackCurragh/RDG/blob/main/images/current_representation.png)</kbd>
  
This figure is a schematic showing the anatomy of a coding RNA transcript according to canonical annotation approaches. Each coding transcript has three parts. Two regions, 5' UTR and 3' UTR, are considered to be untranslated and a single coding region is annotated as protein coding (CDS).      

<kbd>![ORF plot showing how multiple coding regions per transcript are annotated](https://github.com/JackCurragh/RDG/blob/main/images/ORF_plot.png)</kbd>

The top figure shows the CDS annotation corresponding to the ORF organisation shown on the bottom. 

Current approaches for representing multiple coding regions on the same mRNA transcript rely on the generation of two transcripts with the same sequence that differ in CDS annotation. Such a representation distorts the reality by representing the translation of related coding regions as indpendent of each other.

## Graph Data Structures for Biological Information 

<kbd>![](https://github.com/JackCurragh/RDG/blob/main/images/Other_biological_info)</kbd>

Genomic variation such as single nucleotide polymorphisms, repeat expansions and indels in samples have been shown to lead to reference bias in genomic analysis (Paten et al.,2017). Sequence graphs have been developed to support the representation of genomic variation at the population level. At the genome level, each path through such a graph is a potential haplotype. 

At the transcriptome level, the variation introduced through alternative splicing can be similarly represented in a splice graph where each path through such a graph is a potential isoform.

There is however no such structure that supports the representation of protein coding complexity within single transcripts.     

## Introducing Ribosome Decision Graphs 

<kbd>![](https://github.com/JackCurragh/RDG/blob/main/images/RDG.png)</kbd>

We propose a "Ribosome Decision Graph (RDG)" structure (right) to represent RNA transcripts that encode multiple proteoforms (left). Here, start codons form branch points ("decisions") in the graph with each path representing possible translation fates of individual ribosomes.

In the interest of clear visualisation, branch points are only introduced at translation starts allowing the representation of the mechanism of leaky scanning. Similarly, branch points may be introduced for any other mechanisms described above.

<kbd>![](https://github.com/JackCurragh/RDG/blob/main/images/Proportion_of_peptides.png)</kbd>

A figure depicting a RDG augmented with probabilites of ribosome 'decisions' (initiation, termination, etc.). allowing the prediction of relative synthesis rates of individual proteoforms.


<kbd>![](https://github.com/JackCurragh/RDG/blob/main/images/Impact_of_point_mutationspng)</kbd>


A figure depicting a RDG augmented with probabilites of ribosome 'decisions' (initiation, termination, etc.) allowing the prediction of relative synthesis rates of individual proteoforms. In this graph the first start is a non-AUG (CUG) start that is relatively inefficient, only half of the ribosomes recognises it as a start


# Using RDGs 

## Dependencies 

- matplotlib          3.5.0
- networkx            2.6.3 
- sqlitedict          1.7.0

## Installation 
~~~
pip install RDG
~~~
