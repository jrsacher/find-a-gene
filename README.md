# Find-A-Gene Project

Joshua Sacher\
BIOT E-100 Fall 2018

## Requirements

### Python 3

Written using Python 3.7

#### Packages

[Biopython](https://biopython.org/)\
[Matplotlib](https://matplotlib.org/)

### Programs

[Clustal Omega](http://clustal.org/)\
included as `clusatlo`

[PhyML](http://www.atgc-montpellier.fr/phyml/)\
included as `PhyML`

### Note

Program assumes MacOS.

To run on Linux, will likely need to download and recompile Clustal and PhyML.\
To run on Windows ... good luck.

## Documentation

### Options

As the BLAST searches can take a great deal of time to run and are performed on NCBI's server, there are options to use the saved outputs of a previous run. If `tblastn_run` or `novels_run` are set to `True`, the saved files are used; if `False`, they are (re)computed.

The options `tree_all` and `tree_novel` are used to build a phylogenetic tree from all searched sequences or novel sequnces only, respectively. Tree file generation is EXTREMELY computationally expensive:
>A rough estimate [of O] can be obtained with the formula : T \* T \* L \* A \* B \* I\
>T (Taxa) : number of taxa\
>L (Length) : sequences length\
>A (Alphabet) : 1 for nucleotides, 12 for amino-acids\
>B (Bootstrap): number of bootstrap replicates + number of random starting trees\
>I (tree Improvement) : 1 for NNI, 4 for SPR

The program is currently configured with 100 bootstrap replicates, which is reportedly the point of diminishing return. Change `bootstrap=100` to a more manageable number (~10) for faster results.

### How the program works

The user must supply a valid accession number as `accession_in`. This is used in an Entrez query to get information on the starting gene.

The gene of interest is then used as the input for a TBLASTN search. By design, it does not search human, mouse, or rat sequences, which are highly likely to be annotated. A results file of all hits is saved as `tblastn_results.xml`.

The results are then filtered based on identity, similarity, and gaps. The thought is that sequences too similar to the orignal (potentially well-studied) sequence are more likely to be annotated. This theory has not been adequately explored (but, hey, it worked out!). The filtered results are converted to a csv file and saved as `filtered_results.csv`.

A number of the filtered results (default is 10; 25 as currenlty configured) are then checked for novelty. This is accomplished by extracting the gi number, then performing a BLASTP search on the sequence in the non-redundant (nr) database.\
By Pevsner's criteria, a gene is novel if there is not a 100% match from the query organism. If the original gene of interest is not found in the results, it is likely a false positive.\
The gi numbers and novelty (`True`/`False`) of the searched genes are added to the original data and saved as `novel_genes.csv`

Multiple sequence alignment for these genes (all or novels only, as mentioned above) is then performed using ClustalOmega. The input FASTA file and results are saved as `msa_[x]_in.fasta` and `msa_[x]_out.aln`, where `[x]` is "all" or "novel," depending.

The alignment file can be used to generate phylogenetic trees using PhyML. The resulting text file is read back in and converted to XML so style can be applied to the image. The final trees are saved as 300 dpi png files--standard white background and transparent--for display in the report.

### TODOs (wishlist)

* Add color to tree text
* Remove hard-coded comparison sequences in tree and let user specify
* Allow user to specify options, file names, etc. via command line
* Convert to web app?
* Automated report generation
* Add additional analyses of novel genes (separate program?)
