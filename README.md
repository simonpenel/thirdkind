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

single gene reconciliation with species tree:

https://raw.githubusercontent.com/simonpenel/thirdkind/74b9c84a5b2ed84ff5230fc3a52af856b9aba53d/output_examples/thirdkind_example1.svg


the same gene reconciliation without the species tree:

https://raw.githubusercontent.com/simonpenel/thirdkind/74b9c84a5b2ed84ff5230fc3a52af856b9aba53d/output_examples/thirdkind_example1_bis.svg


multiple genes reconciliation:

https://raw.githubusercontent.com/simonpenel/thirdkind/70a7a11aa89eda61926c5cabf221f47bb44e3409/output_examples/thirdkind_example4.svg


example with "free living" species:

https://raw.githubusercontent.com/simonpenel/thirdkind/9eb47ce644998164ff56cc68eb765c0c8a24d389/output_examples/thirdkind_fl_example.svg

multiple gene trees reconciliation with redundant transfers. Display only 1 gene tree and the transfers in red with an opacity according to the abundance of the transfer:

https://raw.githubusercontent.com/simonpenel/thirdkind/70a7a11aa89eda61926c5cabf221f47bb44e3409/output_examples/thirdkind_example2.svg


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
    thirdkind

You may as well install from the sources. This may be useful if you want to use the examples.
Clone or download  the sources here https://github.com/simonpenel/thirdkind and type:

    cargo build --release
    target/release/thirdkind -h

Note:

Since Rust does not include its own linker yet, you need to have a C compiler like gcc installed to act as the linker.


# Run the binary:
Read a newick, phyloxml or recPhyloXML file and create a svg.

Format is guessed according to filename (default is newick)

Usage:

    thirdkind -f input file [-b][-c config file][-d fontsize][-D fontsize][-e][-F format][-g input file][-G #][-h][-H height][-i][-I][-J][-l factor][-L][-m][-o output file][-O][-p][-P][-r ratio][-s][-S][-t threshold][-T #][-u threshold][-U #][-v][-W width]
    
    -b : open svg in browser
    -c configfile: use a configuration file
    -d fontsize: set font size for gene trees
    -D fontsize: set font size for species trees
    -e : the node associated to FREE_LIVING are drawned in an external tree (free_living option)
    -F phylo/recphylo: force format phyloXML/recPhyloXML
    -g 2nd level input file (for example a gene-symbiote file with -f defining a symbiote-host file)
    -G <n> : display the gene #n in phyloxml style (no species tree)
    -h : help
    -H height : multiply the tree height by factor 'height'
    -i : display internal gene nodes
    -I : display internal species nodes
    -J : with option -t, display the abundance of redudant transfers
    -l factor : use branch length, multiplied by the given factor in case of newick/phyloxml,
                and by an optimised factor in case of recphyloxml
    -L : display as landscape
    -m : the input file (-f) is a list of recphyloxml files
    -o outputfile : set name of output file
    -O switching nodes in order to minimise transfer crossings (under development)
    -p : build a phylogram
    -P : pipe species tree uniformisation
    -r ratio : set the ratio between width of species and gene tree.
               Default 1.0, you usualy do not need to change it.
    -s : drawing species tree only
    -S : display node support
    -t <t> : redudant transfers are displayed as one, with opacity according to abundance and only if abundance is higher tan t.
             Only one gene is displayed.
    -T <n> : with option -t, select the gene to display.
    -u <t> : with -g, same as -t, but apply to the '-f' input file, and -t will apply to the '-g' file.
    -U <n> : same as -T with -t, but for -u
    -v : verbose
    -W width : multiply the tree width by factor 'width'

    Note on -b option : you must set a browser as default application for opening svg file

    Note on -g option : this will generate 3-levels reconciliation svg files.
    For example you may input a gene-symbiote recphyloxml file  with -g and symbiote-host recphyloxml file with -f
    The -t/-g options are not totally implemented for the mapped_1/2/3 svg output files.

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

# Configuration file:

You may configure some of the features of the svg with the -c option.

The default values are the values of the "config_default.txt" file.

Modify the default values and save it into  "my_config.txt" then type:

```
thirdkind -f recphylo_examples/FAM000600_reconciliated_big.recphylo -c my_config.txt -b

```

# Using the API light_phylogeny:

thirdkind uses the light_phylogeny library:

https://github.com/simonpenel/light_phylogeny/wiki


# recPhyloXML documentation

See http://phylariane.univ-lyon1.fr/recphyloxml/

recPhyloXML paper: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6198865/

# phyloXML documentation

See: http://www.phyloxml.org/

phyloXML paper: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2774328/

# Licence
CECILL: https://choosealicense.com/licenses/cecill-2.1/
