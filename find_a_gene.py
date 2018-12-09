# -*- coding: utf-8 -*-
"""
Created on Thu 22 Nov 2018
Last update Sat 08 Dec 2018

@author: jrsacher

Find-a-gene v 0.3
See chapter 4 in Bioinformatics and Functional Genomics (3rd ed)
http://bioinfbook.org/php/C4E3k
"""
import csv
import os
from Bio import AlignIO, Entrez, Phylo, SeqIO
from Bio.Align.Applications import ClustalOmegaCommandline
from Bio.Alphabet import IUPAC
from Bio.Blast import NCBIWWW, NCBIXML
from Bio.Phylo.Applications import PhymlCommandline
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from distutils.util import strtobool
from matplotlib import pyplot as plt

# Change this if you aren't me!
Entrez.email = "joshuasacher@gmail.com"


def main():
    # Flags to indicate if portions of the program have run successfully.
    # Set to True if run, else False. Saves a lot of computation time.
    tblastn_run = False
    novels_run = False
    # MSA doesn't take long, so always run
    # Tree building options
    tree_all = True     # build tree from all searched sequences
    tree_novel = True   # build tree from novel (and initial) sequences only

    # Starting protein. Manually selected accession number
    # Synapsin-1 isoform A SYN1a
    accession_in = "NP_008881.2"

    with Entrez.efetch(db="nucleotide", rettype="gb",
                       retmode="text", id=accession_in) as handle:
        search_seq = SeqIO.read(handle, "gb")

    # Run TBlastN on starting protein
    if tblastn_run:
        tblastn_results = "tblastn_result.xml"
    else:
        tblastn_results = tblastn(accession_in)

    # Filter results list based on identity, similarity, and gap percentages
    filtered_results = inspect_results(tblastn_results)

    # Check a certain number of sequences for novelty
    # Identifies known genes, likely false positives, and novel genes
    # WARNING: THIS CAN TAKE A LOT OF TIME
    if novels_run:
        with open("novel_genes.csv", "r") as infile:
            reader = csv.DictReader(infile)
            novel_search = [row for row in reader]
        # From the CSV, "True" and "False" are read as strings! Convert to bool
        for gene in novel_search:
            gene["novel"] = strtobool(gene["novel"].lower())
    else:
        novel_search = find_novels(filtered_results, search_seq, seqs_to_check=50)

    # Perform multiple sequence alignment using clustalo and builds tree using PhyML
    # Bootstrap of 100 replicates considered more than adequate,
    # but will take a LONG time to run, so default is 10
    if tree_all:
        msa_file = global_msa(novel_search, search_seq, file_name="msa_all")
        build_tree(msa_file, search_seq, file_name="tree_all", bootstrap=25)

    if tree_novel:
        novels = [seq for seq in novel_search if seq["novel"]]
        msa_file = global_msa(novels, search_seq, file_name="msa_novel")
        build_tree(msa_file, search_seq, file_name="tree_novel", bootstrap=100)


def tblastn(seq, db="est", expect=0.05, hitlist_size=500):
    """
    Perform TBLASTN on a string seq, save the output file, and return filename
        db: database to search (default is est);
        expect: expect value threshold (default is 0.05);
        hitlist_size: max hits to return (default is 500)
    """
    # List of organisms to NOT search!
    organisms = ["human", "mouse", "rat"]
    entrez_filter = ("all [filter] NOT(" +
                     " OR ".join([o + "[ORGN]" for o in organisms]) + ")")
    
    print("Performing TBLASTN search on {seq}...")
    results = NCBIWWW.qblast("tblastn", db, seq,
                             expect=expect,
                             hitlist_size=hitlist_size,
                             entrez_query=entrez_filter)
    print("Search complete\n")

    # Save results in XML format
    outfile = "tblastn_result.xml"
    with open(outfile, "w") as out:
        out.write(results.read())
    results.close()
    print(f"Results saved as '{outfile}' to {os.getcwd()}")
    return outfile


def inspect_results(file):
    """
    Filters results based on identity, similarity, and gaps
    """
    with open(file, "r") as result_file:
        # switch to "NCBIXML.parse if more than one search term used
        results = NCBIXML.read(result_file)

    # Upper threshold values
    percent_ident = 0.80
    percent_pos = 0.90
    percent_gaps = 0.10

    print(f"Filtering results based on:\n"
          f"\tpercent identity < {percent_ident}\n"
          f"\tpercent positives < {percent_pos}\n"
          f"\tpercent gaps < {percent_gaps}")
    
    result_data = []
    kept, discarded = 0, 0
    # For diagram of BLAST result object, see
    # http://biopython.org/DIST/docs/tutorial/Tutorial.html#fig:blastrecord
    for alignment in results.alignments:
        for hsp in alignment.hsps:
            if (hsp.identities / hsp.align_length < percent_ident
                    and hsp.positives / hsp.align_length < percent_pos
                    and hsp.gaps / hsp.align_length < percent_gaps):
                result_data.append({"title": alignment.title,
                                    "length": alignment.length,
                                    "expect": hsp.expect,
                                    "identities": hsp.identities,
                                    "positives": hsp.positives,
                                    "gaps": hsp.gaps,
                                    "align_len": hsp.align_length,
                                    "query": hsp.query,
                                    "match": hsp.match,
                                    "subject": hsp.sbjct})
                kept += 1
            else:
                discarded += 1
    print(f"{discarded} records removed; {kept}/{kept + discarded} remain.")

    # Save dict as CSV file in case it's needed again
    outfile = "filtered_results.csv"
    with open(outfile, "w", newline='') as csvfile:
        fieldnames = list(result_data[0].keys())
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(result_data)
    print(f"Filtered data saved as '{outfile}' to {os.getcwd()}\n")

    return result_data


def find_novels(filtered_results, search_seq, seqs_to_check=10):
    """
    Takes in list of dicts containing filtered results, the original search
    sequence, and the maximum number of sequences to check for novelty.
    Returns a subset of the filtered results with the added fields
    "gi" (str) and "novel" (bool)
    """
    # Extract gi numbers from filtered results
    # -- probably not ideal since gi is no longer preferred
    # All titles start with gi|nnnnnnnn|
    for result in filtered_results:
        t = result["title"]
        stop = t.find("|", 3)  # Find inex of end "|"
        result["gi"] = t[3: stop]

    print(f"Checking {min([seqs_to_check, len(filtered_results)])} sequences for novelty...\n")

    # List to store genes after novelty check
    genes = []

    # Check a number of sequences for novelty
    for seq in filtered_results[:min([seqs_to_check, len(filtered_results)])]:
        # Check for duplication before running
        if seq["gi"] not in [gene["gi"] for gene in genes]:
            seq["novel"] = novel_check(seq, search_seq)
            if seq["novel"]:
                seq["gi"] += " ***"     # Make text stand out for novels!
            genes.append(seq)

    # Save results as CSV file
    outfile = "novel_genes.csv"
    with open(outfile, "w", newline='') as csvfile:
        fieldnames = genes[0].keys()
        writer = csv.DictWriter(csvfile, fieldnames)
        writer.writeheader()
        writer.writerows(genes)
    print(f"Gene novelty data saved as '{outfile}' to {os.getcwd()}\n")

    return genes


def novel_check(seq, search_seq):
    """
    Check for novelty by BLASTP search on non-redundant database.
    Return True if novel, else False.
    Novelty criteria:
        NOT novel if 100% identity match from same species
        Likely novel if 100% identity but different species
        Likely novel if <100% identiy matches only
        BUT if original query not found in matches, likely a false positive
    """
    # BLASTP against nr database to confirm novelty
    # Increased hitlist size from 100 to 500 for better false-positive checking
    # Likely should match hitlist_size used in initial search
    print(f"Performing BLASTP search on protein sequence of gi|{seq['gi']}...")
    blast = NCBIWWW.qblast("blastp", "nr", seq["subject"], hitlist_size=500)
    print("Search complete")

    # Save results in XML format
    file = f"BLASTP/{seq['gi']}.xml"
    with open(file, "w") as outfile:
        outfile.write(blast.read())
    blast.close()
    print(f"Results saved as '{file}' to {os.getcwd()}/BLASTP")

    # Read results back in (.read() is destructive for BLAST objects)
    with open(file, "r") as infile:
        # switch to "NCBIXML.parse if more than one search term used
        blast_results = NCBIXML.read(infile)

    is_novel = False
    for alignment in blast_results.alignments:
        # If original accession isn't a hit, it's likely a false positive
        if search_seq.id in alignment.title:
            is_novel = True
        # Known genes have 100% identity in the same species (# identical AAs = length)
        for hsp in alignment.hsps:
            if (search_seq.annotations['organism'] not in alignment.title
                    and hsp.identities == hsp.align_length):
                print("--- Known gene ---\n")
                return False
    if is_novel:
        print("*** Likely novel gene! ***\n")
        print(f"Original: {search_seq.seq[:50]}...{search_seq.seq[-7:]}")
        print(f"Novel   : {seq['subject'][:50]}...{seq['subject'][-7:]}\n")
    else:
        print("--- Potential false postive ---\n")
    return is_novel


def global_msa(matches, search_seq, file_name="msa"):
    """
    Generates a global multiple sequence alignment from a list of seqs
    """
    # Build list of sequences for input FASTA file
    # Start with original search + others for comparison
    # NOTE: synapsin paralogs are currently hard-coded
    seqs = [search_seq,
            SeqIO.read("synapsinIIa.fasta", "fasta"),
            SeqIO.read("synapsinIIb.fasta", "fasta"),
            SeqIO.read("synapsinIII.fasta", "fasta")]
#            SeqIO.read("GFPclover.fasta", "fasta")]
#            SeqIO.read("FireflyLuciferase.fasta", "fasta")]

    # Build Biopython Seq objects from sequences
    for match in matches:
        # Remove gaps and stop from AA sequence to build Seq object's sequence
        seq = SeqRecord(Seq(match["subject"].replace("-", "").replace("*", ""),
                        IUPAC.protein),
                        id="gi|" + match["gi"],
                        description=match["title"])
        seqs.append(seq)

    # Make FASTA file from sequences
    infile = f"{file_name}_in.fasta"
    SeqIO.write(seqs, infile, "fasta")

    # Set up and run ClustalOmega
    clustal = ClustalOmegaCommandline()
    clustal.program_name = "./clustalo"
    outfile = f"{file_name}_out.aln"
    clustal.outfmt = "clustal"
    clustal.infile = infile
    clustal.outfile = outfile
    clustal.force = True    # Allows overwriting file

    # Run MSA, print success/failure
    print(f"Performing multiple sequence alignment on {len(seqs)} sequences")
    stdout, stderr = clustal()
    print(stdout + stderr)
    print(f"Results saved as '{outfile}' to {os.getcwd()}\n")

    # Return MSA file name
    return outfile


def build_tree(msa_file, original, file_name="tree", bootstrap=10):
    """
    Build a phylogenetic tree based on a multiple sequence alignment.
        http://biopython.org/DIST/docs/tutorial/Tutorial.html#htoc217
        PhyML site: http://www.atgc-montpellier.fr/
        "New Algorithms and Methods to Estimate Maximum-Likelihood Phylogenies:
        Assessing the Performance of PhyML 3.0." Guindon S., Dufayard J.F.,
        Lefort V., Anisimova M., Hordijk W., Gascuel O. Systematic Biology,
        59(3):307-21, 2010.
    Color based on novelty
    """
    # See https://biopython.org/wiki/Phylo for general code and settings
    # Convert to phyllip format
    tree_file = file_name   # for convenience
    AlignIO.convert(msa_file, "clustal", tree_file, "phylip-relaxed",
                    alphabet=IUPAC.protein)

    # NOTE: make sure file name is "phyml"
    PhyML = PhymlCommandline("./PhyML")
    PhyML.input = tree_file
    PhyML.datatype = 'aa'  # Specify that amino acids are being input
    PhyML.model = 'LG'  # Amino acid substitution matrix
    PhyML.alpha = 'e'
    # non-parametric bootstrap relplicates; 100 is point of dimiishing returns
    PhyML.bootstrap = bootstrap

    # Run tree generation, print success/failure
    print("Building distance tree from multiple sequence alignment...\n")
    stdout, stderr = PhyML()
    print(stdout + stderr)
    print(f"Newick tree saved as {tree_file + '_phyml_tree.txt'}")

    # Read in tree file, convert to XML (to be able to add color, etc.)
    tree = Phylo.read(tree_file + "_phyml_tree.txt", "newick")
    tree = tree.as_phyloxml()

    # Stylize the tree
    # Colorblind-safe colors can be checked with ColorOracle: http://colororacle.org/
    for clade in tree.find_clades():
        # Bold lines
        clade.width = 3
        # Red if known gene or false positive
        if str(clade.name).startswith("gi|"):
            clade.color = "#e4002b"
        # Blue for originally searched gene
        elif clade.name == original.id:
            clade.color = "#006db6"
        # Black for comparitor nodes
        elif clade.name is not None and not clade.color:
            clade.color = "#000000"
        # Gray for non-terminal nodes
        elif not clade.name:
            clade.color = "#63666a"
        # Green for novel genes
        if str(clade.name).endswith("***"):
            clade.color = "#00bf71"

    # Configure plot. Image size determined based on number of nodes
    tree_len = len(tree.get_terminals())
    plt.rc("font", size=18)  # Bigger font for easier reading
    fig = plt.figure(figsize=(1.6 * tree_len, tree_len), dpi=300)
    axes = fig.add_subplot(1, 1, 1)
    Phylo.draw(tree, axes=axes, do_show=False)
    # Save white background image
    fig.savefig(f"{tree_file}.png", format='png',
                bbox_inches='tight', dpi=300)
    # Save transparent image
    fig.savefig(f"{tree_file}_transparent.png", format='png',
                bbox_inches='tight', dpi=300, transparent=True)
    print(f"Tree images saved as {tree_file + '.png'} "
          f"and {tree_file + '_transparent.png'} to {os.getcwd()}\n")


if __name__ == "__main__":
    main()
