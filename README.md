# thirdkind
Drawing reconciled phylogenetic trees  allowing 1, 2 or 3 reconciliation levels

Build  svg representations of  phylogenetic reconciled (or not) trees with events (loss, duplication, speciation, transfer).

* Input one newick or phyloxml file -> a svg representation of the tree with node events

* Input one recphyloxml file -> a svg representation of the   gene (or symbiot) 'lower' tree(s) inside the associated species (or host) 'upper' tree

* Input a set of recphyloxml files -> a svg representation of all the   gene (or symbiot) 'lower' tree(s) inside the associated species (or host) 'upper' tree, or only 1 gene tree  with the transfers of all genes (or symbiotes)  allowing to visualise transfers for different scenarios or histories.

* Input two nested recphyloxml files -> several svg representations allowing to display 3 level reconciliations (for example gene/symbiot/host)

# Homepage
homepage : https://github.com/simonpenel/thirdkind/wiki

[![thirdkind at crates.io](https://img.shields.io/crates/v/thirdkind.svg)](https://crates.io/crates/thirdkind)
[![License: GPL v2](https://img.shields.io/badge/License-GPL_v2-blue.svg)](https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html)



# Keywords:
phylogeny, reconciled trees, phylogenetic trees, reconciliation, visualisation, svg, recphyloxml, phyloxml, 3 levels reconciliation, gene transfer, speciation, duplication, loss

# Formats:

phyloXML, recPhyloXML, rooted newick ( NHX balises will not be considered ).

# Graphical interface

A web sever dedicated to the graphical interface of _thirdkind_ is available here: http://thirdkind.univ-lyon1.fr/

# Output examples of 2-levels visualisation:

single gene reconciliation with species tree:

https://raw.githubusercontent.com/simonpenel/thirdkind/74b9c84a5b2ed84ff5230fc3a52af856b9aba53d/output_examples/thirdkind_example1.svg


the same gene reconciliation without the species tree:

https://raw.githubusercontent.com/simonpenel/thirdkind/74b9c84a5b2ed84ff5230fc3a52af856b9aba53d/output_examples/thirdkind_example1_bis.svg


multiple genes reconciliation:

https://raw.githubusercontent.com/simonpenel/thirdkind/70a7a11aa89eda61926c5cabf221f47bb44e3409/output_examples/thirdkind_example4.svg

real length with recPhyloXML file:

https://raw.githubusercontent.com/simonpenel/thirdkind/b084872f02a758702e0f90543c715862729166a5/output_examples/thirdkind_real_length.svg

using a  "free living" species:

https://raw.githubusercontent.com/simonpenel/thirdkind/0f8ff64838fc7676b2e98542e953c4ccc9f45c62/output_examples/thirdkind_fl_example2.svg


multiple "spTree":

https://raw.githubusercontent.com/simonpenel/thirdkind/79344589aa8b91a909386d22c95a01bc8c795396/output_examples/thirdkind_multspec.svg

multiple gene trees reconciliation with redundant transfers. Display only 1 gene tree and the transfers in red with an opacity according to the abundance of the transfer:

https://raw.githubusercontent.com/simonpenel/thirdkind/70a7a11aa89eda61926c5cabf221f47bb44e3409/output_examples/thirdkind_example2.svg


# Output examples of 3-levels visualisation:

## Example 1

reconciled symbiote tree with associated gene trees:
https://raw.githubusercontent.com/simonpenel/thirdkind/abfb9e6a28d03860bea43b52312dc706554fd53d/output_examples/thirdkind_example3_mapped_1.svg

symbiote-host reconciliation plus mapped gene transfers:
https://raw.githubusercontent.com/simonpenel/thirdkind/abfb9e6a28d03860bea43b52312dc706554fd53d/output_examples/thirdkind_example3_mapped_2.svg

host tree with associated gene trees:
https://raw.githubusercontent.com/simonpenel/thirdkind/abfb9e6a28d03860bea43b52312dc706554fd53d/output_examples/thirdkind_example3_mapped_3.svg


## Example 2

reconciled symbiote trees with associated gene trees:
https://raw.githubusercontent.com/simonpenel/thirdkind/abfb9e6a28d03860bea43b52312dc706554fd53d/output_examples/thirdkind_example4_mapped_1.svg


symbiote-host reconciliation plus mapped gene transfers:
https://raw.githubusercontent.com/simonpenel/thirdkind/abfb9e6a28d03860bea43b52312dc706554fd53d/output_examples/thirdkind_example4_mapped_2.svg


host tree with associated gene trees:
https://raw.githubusercontent.com/simonpenel/thirdkind/abfb9e6a28d03860bea43b52312dc706554fd53d/output_examples/thirdkind_example4_mapped_3.svg


# Install:

_thirdkind_ is written in Rust. The code is managed using Cargo and published on crates.io.

Install cargo:

    curl https://sh.rustup.rs -sSf | sh

or for Windows see  https://doc.rust-lang.org/cargo/getting-started/installation.html

**Note:**

Since Rust does not include its own linker yet,building _thirdkind_ needs to have a C compiler like gcc installed to act as the linker.
If it is note the case,  install essential build needed by Rust:

    sudo apt install build-essential

Once Cargo is installed just open a new terminal and type:

    cargo install thirdkind
    thirdkind

You may as well install from the sources. This may be useful if you want to use the examples.
Clone or download  the sources here https://github.com/simonpenel/thirdkind and type:

    cargo build --release
    target/release/thirdkind -h



# Run the binary:
Read a newick, phyloxml or recPhyloXML file and create a svg.

Format is guessed according to filename (default is newick)

Usage:

    thirdkind -f input file [-a][-A stArt][-b][-B][-c config file][-d fontsize][-D fontsize][-e][-E][-F format][-g input file][-G #][-h][-H height][-i][-I][-J][-k symbol size][-K Bezier parameter][-l factor][-L][-m][-M][-N eNd][-o output file][-O][-p][-r ratio][-s][-S][-t threshold][-T #][-u threshold][-U #][-v][-W width]|-x][-X][-z thickness][-Z thickness]

    -a : output the redundant transfers analysis
    -A node name : display transfers starting from this node only
    -b : open svg in browser
    -B : with option -l, display branch length
    -c configfile : use a configuration file
    -d fontsize : set font size for gene trees
    -D fontsize : set font size for species trees
    -e : the node associated to FREE_LIVING are drawned in an external tree (free_living option) and superposed in case of multiple genes
    -E : the node associated to FREE_LIVING are drawned in an external tree (free_living option) and shifted in case of multiple genes
    -F phylo/recphylo : force format phyloXML/recPhyloXML
    -g 2nd level input file (for example a gene-symbiote file with -f defining a symbiote-host file)
    -G <n> : display the gene #n in phyloxml style (no species tree)
    -h : help
    -H height : multiply the tree height by factor 'height'
    -i : display internal gene nodes
    -I : display internal species nodes
    -J : with option -t, display the abundance of redudant transfers
    -k size: size of the circles, crosses, squares, etc.
    -K Bezier parameter: curvature of the transfers and branches leading to free living organisms.
    -l factor : use branch length, multiplied by the given factor
    -L : display as landscape
    -m : the input file (-f) is a list of recphyloxml files
    -M : display duplication node at mid-distance in the branch (in progress)  
    -N node name : display transfers ending to this node only
    -o outputfile/prefix : set the name of the output file/set the prefix of the output files
    -O : switching nodes in order to minimise transfer crossings (under development)
    -p : species 'upper' tree uniformisation
    -r ratio : set the ratio between width of species and gene tree
               Default 1.0, you usualy do not need to change it
    -s : drawing species tree only
    -S : display node support
    -t <t> : redudant transfers are displayed as one, with opacity according to abundance and only if abundance is higher tan t
             Only one gene is displayed
    -T <n> : with option -t, select the gene to display. If set to 0, no gene is displayed
    -u <t> : with -g, same as -t, but apply to the '-f' input file, and -t will apply to the '-g' file
    -U <n> : same as -T with -t, but for -u
    -v : verbose
    -W width : multiply the tree width by factor 'width'
    -x : tidy mode (non-layered tidy tree layout)
    -X : tidy mode, avoiding leave names superposition
    -z thickness: thickness of the gene tree
    -Z thickness: thickness of the species tree


    Note on -b option : you must set a browser as default application for opening svg file

    Note on -g option : this will generate 3-levels reconciliation svg files
    For example you may input a gene-symbiote recphyloxml file  with -g and symbiote-host recphyloxml file with -f
    The -t/-u options are not totally implemented for the mapped_1/2/3 svg output files

    Note on -x/-X options : the non-layered tidy tree layout is described in :
                            'van der Ploeg, A. 2014. Drawing non-layered tidy trees in linear time.
                            Software: Practice and Experience, 44(12): 1467–1484.'

`Input format is guessed according to the file name extension:`

    - .phyloxml    => phyloXML
    - .xml         => recPhyloxml
    - .recphyloxml => recPhyloXML
    - .recPhyloXML => recPhyloXML
    - .recphylo    => recPhyloXML
    - any other    => newick

You will find several input file examples in recphylo_examples and xml_examples directories.

# Examples:
    thirdkind -f recphylo_examples/FAM000297_reconciliated.recphylo  -b
    thirdkind -f recphylo_examples/concat.xml -b -t 0
    thirdkind -f recphylo_examples/hote_parasite_page4_BL.recphylo  -b -l 1
    thirdkind -f recphylo_examples/testfiles -m -b -t 3 -J
    thirdkind -f recphylo_examples/test2/hote_parasite_page2.recphylo  -g recphylo_examples/test2/gene_parasite_page2.recphylo  -b  
    thirdkind -f recphylo_examples/test1_mult_parasite/rechp_dtl.recphyloxml -g recphylo_examples/test1_mult_parasite/recgs_mult_host_dtl.recphyloxml -b
    thirdkind -f newick_examples/virus.nhx -l 4  -b
    thirdkind -f newick_examples/virus.nhx -l 4 -x  -b
    thirdkind -f newick_examples/virus.nhx -l 4 -X  -b
# Configuration file:

You may configure some of the features of the svg with the -c option.

The default values are the values of the "config_default.txt" file.

Modify the default values and save it into  "my_config.txt" then type:

```
thirdkind -f recphylo_examples/FAM000600_reconciliated_big.recphylo -c my_config.txt -b

```

Contents of a configuration file :

    species_color:pink
    species_opacity:0.9
    single_gene_color:blue
    gene_opacity:0.9
    species_police_color:orange
    species_police_size:12
    gene_police_size:10
    bezier:1

# Using the API light_phylogeny:

_thirdkind_ uses the light_phylogeny library:

https://github.com/simonpenel/light_phylogeny/wiki


# recPhyloXML documentation

See http://phylariane.univ-lyon1.fr/recphyloxml/

recPhyloXML paper: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6198865/

# phyloXML documentation

See: http://www.phyloxml.org/

phyloXML paper: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2774328/

# Citation

https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/btac062/6525213


# Contact

https://lbbe.univ-lyon1.fr/fr/annuaire-des-membres/penel-simon

# Licence

CECILL: https://choosealicense.com/licenses/cecill-2.1/


