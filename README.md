# Gene functional analysis
Repository of scripts and tools for the functional analysis of genes (GO terms, pathway analysis, enrichments)

We store here scripts in `R`, `Python` (and other languages if necessary) that can be used in the different steps of **gene ontology (GO)**, **enrichment** and **pathway analysis**.
This is a tpye of analysis which is typically conducted after GWAS ro gene expression/RNAseq experiments.

GO, enrichment and pathway analysis are highly complex types of analysis, that require interactions with multiple annotation, genomic and metabolic databases, which are constantly being updated.
Therefore, there are usually multiple ways to do the same thing, and there may be multiple points of entry (e.g. input data) that may require preprocessing and conversion steps.
We will try to document and comment steps, operations and scripts in the best possible way, trying to trace as much as possible scripts that must follow in a logical succession of steps.

Details can be found in our [wiki](https://github.com/filippob/gene_functional_analysis/wiki/)

### GO and pathway functional enrichment analysis

A gene set (e.g. from the results of GWAS for a specific phenotype like a disease) is a collection of genes that are supposed to be functionally related. 

GO terms are organized in a **directed acyclic graph** (DAG) (edges between terms represent parent-child relationship).
We can investigate biological processes (BP), molecular functions (MF), and cellular components (CC) associated with a set of genes.

Over Representation Analysis (ORA) [Boyle et al. 2004](https://academic.oup.com/bioinformatics/article/20/18/3710/202612) is used to determine whether GO terms (or pathways) are **over-represented** (**= enriched**) in a gene list (e.g. from GWAS).
The **p-value** is calculated from the **hypergeometric distribution** (distribution of 'successes' *with replacement*):

$$
\text{p-value} = 1- \sum_{i=0}^{k-1} \frac{{M\choose i} \cdot {N-M\choose n-i}}{{N\choose n}}
$$

where $N$ is the number of background genes, $M < N$ is the n. of genes in $[N]$ that are annotated $n$ is the size of the target list of genes,
and $k < n$ is the n. of genes in our list that are annotated.


### KEGG pathways

Collection of pathways of molecular interaction and reaction networks, for a wide range of biochemical processes/categories:

1. Metabolism
2. Genetic information processing
3. Environmental information processing
4. Cellular processes
5. Organismal systems
6. Human diseases
7. Drug development.
