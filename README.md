# thirdkind
Drawing reconciled phylogenetic trees  allowing 1, 2 or 3 reconciliation levels

[![thirdkind at crates.io](https://img.shields.io/crates/v/thirdkind.svg)](https://crates.io/crates/thirdkind)



Build  svg representations of  phylogenetic reconciled (or not) trees with events (loss, duplication, speciation, transfer).

* Input one newick or phyloxml file -> a svg representation of the tree with node events

* Input one recphyloxml file -> a svg representation of the  "path" gene (or symbiot) tree(s) inside the associated "pipe" species (or host) tree

* Input two recphyloxml files -> several svg representations allowing to display 3 level reconciliations (for example gene/symbiot/host)


Keywords:  phylogeny, reconciled trees, phylogenetic trees

# Formats:

phyloXML, recPhyloXML, rooted newick ( NHX balises will not be considered ).

# Output examples of 2-levels visualisation:

single gene reconciliation in recphyloxml:

https://raw.githubusercontent.com/simonpenel/rectree2svg/5e23e92396d44a68337b33c579f2d9d372d18b4d/FAM000696_reconciliated_recphyloxml.svg

the same gene reconciliation in phyloxml:

https://raw.githubusercontent.com/simonpenel/rectree2svg/9244f3136961f909fd7b33818f0a220e3f32c880/FAM000696_reconciliated_xml.svg

multiple genes reconciliation recphyloxml:

https://raw.githubusercontent.com/simonpenel/rectree2svg/684e2ecccf04408cf15815e3fb71bbeebfa19d12/tree2svg.example.svg

multiple gene trees with redundant transfers. Display only 1 gene tree and the transfers according to the abundance of the transfer:

https://raw.githubusercontent.com/simonpenel/rectree2svg/f9ee031fa23b815ff7cc7298fd0dc4fb45707d53/transfers_abundance.svg

# Output examples of 3-levels visualisation::

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

thirdkind is written in Rust. The code is managed using Cargo and published on crates.io.

Install cargo:

    curl https://sh.rustup.rs -sSf | sh

or for Windows see  https://doc.rust-lang.org/cargo/getting-started/installation.html

Once Cargo is installed just open a new terminal and type:

    cargo install thirdkind

You may as well install from the sources. Clone or download  the sources here https://github.com/simonpenel/thirdkind and type:

    cargo build --release
    target/release/thirdkind -h

Note:

Since Rust does not include its own linker yet, you need to have a C compiler like gcc installed to act as the linker.


# Run the binary:
Read a newick, phyloxml or recPhyloXML file and create a svg.

Format is guessed according to filename (default is newick)

Usage:

    thirdkind -f input file [-b][-c config file][-F format][-g input file][-G #][-h][-i][-I][-l factor][-o output file][-p][-r ratio][-s][-v]

    -b : open svg in browser
    -c configfile: use a configuration file    
    -F phylo/recphylo: force format phyloXML/recPhyloXML    
    -g 2nd level input file
    -G <n> : display the gene #n in phyloxml style (no species tree)
    -h : help    
    -H height : multiply the tree height by factor 'height' (default 1.0)
    -i : display internal gene nodes
    -I : display internal species nodes
    -J : with option -t, display the abundance of redudant transfers
    -L : display as landscape
    -l factor: use branch length, using the given factor
    -o outputfile : set name of output file    
    -O switching nodes in order to minimise transfer crossings (under development)
    -p : build a phylogram   
    -r ratio : set the ratio between width of species and gene tree.
               Default 1.0, you usualy do not need to change it.
    -s : drawing species tree only    
    -S : display node support
    -t <t> : redudant transfers are displayed as one, with opacity according to abundance and only if abundance is higher tan t. Only one gene is displayed.
    -T <n> : with option -t, select the gene to display
    -v : verbose   
    -W width : multiply the tree width by factor 'width' (default 1.0)

    Note on -b option : you must set a browser as default application for opening svg file

    Note on -g option : this will generate 3-levels reconciliation svg. For example you may input a gene-symbiote recphyloxml file  with -g and symbiote-host recphyloxml file with -f


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
    thirdkind -f recphylo_examples/test2/hote_parasite_page2.recphylo  -g recphylo_examples/test2/gene_parasite_page2.recphylo  -b  
    thirdkind -f recphylo_examples/test1_mult_parasite/rechp_dtl.recphyloxml -g recphylo_examples/test1_mult_parasite/recgs_mult_host_dtl.recphyloxml -b

# Configuration file:

You may configure some of the features of the svg with the -c option.

The default values are the values of the "config_default.txt" file.

Modify the default values and save it into  "my_config.txt" then type:

```
thirdkind -f recphylo_examples/FAM000600_reconciliated_big.recphylo -c my_config.txt -b

```

# Using the API:

rectree2svg use the light_phylogeny crate:

https://crates.io/crates/light_phylogeny


# recPhyloXML documentation

See http://phylariane.univ-lyon1.fr/recphyloxml/

recPhyloXML paper: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6198865/

# phyloXML documentation

See: http://www.phyloxml.org/

phyloXML paper: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2774328/

# Licence
CECILL: https://choosealicense.com/licenses/cecill-2.1/
