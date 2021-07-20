/// name = "rectree2svg"
/// version = "2.2.0"
/// authors = ["Simon Penel <simon.penel@univ-lyon1.fr>"]
/// release = "18/04/2021"
/// license = "CECILL-2.1"
///
/// Usage:
/// Draw phylogenetic trees in a svg file.
/// Draw multiple reconciled gene trees with the associated species tree.
/// Draw simple gene or species tree too.
/// Read a newick, phyloxml or recPhyloXML file.
/// Format is guessed according to filename (default is newick).

use std::fs;
use std::env;
use std::process;
use getopt::Opt;
use webbrowser::{Browser};
use light_phylogeny::*;
use log::{info};

// Message d'erreur
// ----------------
fn display_help(programe_name:String) {
    const VERSION: Option<&'static str> = option_env!("CARGO_PKG_VERSION");
    const NAME: Option<&'static str> = option_env!("CARGO_PKG_NAME");
    const DESCRIPTION: Option<&'static str> = option_env!("CARGO_PKG_DESCRIPTION");
    println!("{} v{}", NAME.unwrap_or("unknown"),VERSION.unwrap_or("unknown"));
    println!("{}", DESCRIPTION.unwrap_or("unknown"));
    println!("Usage:");
    println!("{} -f input file [-b][-c config file][-e][-F format][-g input file][-G #][-h]\
    [-H height][-i][-I][-J][-l factor][-L][-m][-o output file][-O][-p][-r ratio][-s][-S]\
    [-t threshold][-t #][-v][-W width]",programe_name);
    println!("    -b : open svg in browser");
    println!("    -c configfile: use a configuration file");
    println!("    -e : the node associated to FREE_LIVING are drawned in an \
    external tree (free_living option).");
    println!("    -F phylo/recphylo: force format phyloXML/recPhyloXML");
    println!("    -g 2nd level input file");
    println!("    -G <n> : display the gene #n in phyloxml style (no species tree)");
    println!("    -h : help");
    println!("    -H height : multiply the tree height by factor 'height' (default 1.0)");
    println!("    -i : display internal gene nodes");
    println!("    -I : display internal species nodes");
    println!("    -J : with option -t, display the abundance of redudant transfers");
    println!("    -l factor : use branch length, using the given factor (default 1.0)");
    println!("    -L : display as landscape");
    println!("    -m : the input file (-f) is a liss of recphyloxml files");
    println!("    -o outputfile : set name of output file");
    println!("    -O switching nodes in order to minimise transfer crossings (under development) ");
    println!("    -p : build a phylogram");
    println!("    -r ratio : set the ratio between width of species and gene tree.");
    println!("               Default 1.0, you usualy do not need to change it. ");

    println!("    -s : drawing species tree only");
    println!("    -S : display node support");
    println!("    -t <t> : redudant transfers are displayed as one, with opacity according \
    to abundance and only if abundance is higher tan t. Only one gene is displayed.");
    println!("    -T <n> : with option -t, select the gene to display");
    println!("    -v : verbose");
    println!("    -W width : multiply the tree width by factor 'width' (default 1.0)");
    println!("");
    println!("    Note on -b option : you must set a browser as default application for opening \
    svg file");
    println!("");
    println!("    Note on -g option : this will generate 3-levels reconciliation svg.");
    println!("    For example you may input a gene-symbiote recphyloxml file  with -g and symbiote-host recphyloxml file with -f");
    println!("");
    println!("Input format is guessed according to the file name extension:");

    println!(".phyloxml    => phyloXML");
    println!(".xml         => recPhyloxml");
    println!(".recphyloxml => recPhyloXML");
    println!(".recPhyloXML => recPhyloXML");
    println!(".recphylo    => recPhyloXML");
    println!("any other    => newick");
    println!("");
    println!("About recPhyloXML format: http://phylariane.univ-lyon1.fr/recphyloxml/");
    println!("recPhyloXML paper: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6198865/");
    println!("About phyloXML format: http://www.phyloxml.org/");
    println!("phyloXML paper: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2774328/");
    println!("");
    println!("Examples:");
    println!("{} -f recphylo_examples/FAM000297_reconciliated.recphylo  -b", programe_name);
    println!("{} -f recphylo_examples/concat.xml -b -t 0 ", programe_name);
    println!("{} -f recphylo_examples/test2/hote_parasite_page2.recphylo  \
    -g recphylo_examples/test2/gene_parasite_page2.recphylo  -b  ", programe_name);
    println!("{} -f recphylo_examples/test1_mult_parasite/rechp_dtl.recphyloxml \
     -g recphylo_examples/test1_mult_parasite/recgs_mult_host_dtl.recphyloxml -b", programe_name);


    process::exit(1);
}
/// enum of the possible input file Formats
#[derive(Debug)]
pub enum  Format {
    Newick,
    Phyloxml,
    Recphyloxml,
}

fn main()  {
    // Initialise les options
    let mut options: Options =  Options::new();
    // Initialise la config
    let mut config: Config = Config::new();
    // Charge la config par deuakt si elle existe
    let fconf = "config_default.txt";
     if fs::metadata(fconf).is_ok() {
         set_config(fconf.to_string(), &mut config);

     }
    // Gestion des arguments et des options
    // ------------------------------------
    let args: Vec<String> = std::env::args().collect();
    if args.len() == 1 {
         display_help(args[0].to_string());
    }
    let mut opts = getopt::Parser::new(&args, "c:bef:F:g:G:hH:iIJl:Lmo:Opr:sSt:T:vW:");
    let mut infile_sh = String::new(); // symbiote host file
    let mut infile_gs = String::new(); // gene symbiote file
    let mut outfile = String::from("thirdkind.svg");
    let mut nb_args = 0;
    let mut level3 = false; // Affichage à 3 niveaux
    let mut multiple_files = false;
    let mut _format = Format::Newick;
    loop {
        match opts.next().transpose() {
            Err(err) => {
                eprintln!("Error : {}",err);
                std::process::exit(1);
            },
            Ok(res) => match res {
                None => break,
                Some(opt) => match opt {
                    Opt('F', Some(string)) => {
                        _format = match string.as_str() {
                            "recphylo" => Format::Recphyloxml,
                            "phyloxml" => Format::Phyloxml,
                            _ => {
                                eprintln!("Error! Please give a correct format (recphylo/phyloxml)");
                                process::exit(1);
                            },
                        };
                    },
                    Opt('g', Some(string)) => {
                        infile_gs = string.clone();
                        level3 = true;
                    },
                    Opt('G', Some(string)) => {
                        options.disp_gene = match string.parse::<usize>(){
                            Ok(valeur) => valeur,
                            Err(_err) => {
                                eprintln!("Error! Please give a integer value with -G option");
                                process::exit(1);
                            },
                        };
                    },
                    Opt('H', Some(string)) => {
                        options.height = match string.parse::<f32>(){
                            Ok(valeur) => valeur,
                            Err(_err) => {
                                eprintln!("Error! Please give a numeric value with -H option");
                                process::exit(1);
                            },
                        };
                    },
                    Opt('e', None) => options.free_living = true,
                    Opt('i', None) => options.gene_internal = true,
                    Opt('I', None) => options.species_internal = true,
                    Opt('J', None) => options.thickness_disp_score = true,
                    Opt('m', None) => multiple_files = true,
                    Opt('b', None) => options.open_browser = true,
                    Opt('r', Some(string)) => {
                        options.ratio = match string.parse::<f32>(){
                            Ok(valeur) => valeur,
                            Err(_err) => {
                                eprintln!("Error! Please give a numeric value with -r option");
                                process::exit(1);
                            },
                        };
                    },
                    Opt('p', None) => options.clado_flag = false,
                    Opt('s', None) => options.species_only_flag = true,
                    Opt('S', None) => options.support = true,
                    Opt('t', Some(string)) => {
                        options.thickness_thresh = match string.parse::<usize>(){
                            Ok(valeur) => valeur,
                            Err(_err) => {
                                eprintln!("Error! Please give a integer value with -t option");
                                process::exit(1);
                            },
                        };
                        options.thickness_flag = true;
                    },
                    Opt('T', Some(string)) => {
                        options.thickness_gene = match string.parse::<usize>(){
                            Ok(valeur) => valeur,
                            Err(_err) => {
                                eprintln!("Error! Please give a integer value with -T option");
                                process::exit(1);
                            },
                        };
                    },
                    Opt('l', Some(string)) => {
                        options.real_length_flag = true;
                        options.scale = match string.parse::<f32>(){
                            Ok(valeur) => valeur,
                            Err(_err) => {
                                eprintln!("Error! Please give a numeric value with -l option");
                                process::exit(1);
                            },
                        };
                    },
                    Opt('L', None) => options.rotate = false,
                    Opt('v', None) => {
                        options.verbose = true;
                        env::set_var("RUST_LOG", "info");
                        env_logger::init();
                        info!("Verbosity set to Info");
                        },
                    Opt('c', Some(string)) => {
                        set_config(string, &mut config);
                    },
                    Opt('f', Some(string)) => {
                        infile_sh = string.clone();
                        nb_args += 1;
                    },
                    Opt('o', Some(string)) => outfile = string.clone(),
                    Opt('O', None) => options.optimisation = true,
                    Opt('h', None) => display_help(args[0].to_string()),
                    Opt('W', Some(string)) => {
                        options.width = match string.parse::<f32>(){
                            Ok(valeur) => valeur,
                            Err(_err) => {
                                eprintln!("Error! Please give a numeric value with -W option");
                                process::exit(1);
                            },
                        };
                    },
                    _ => unreachable!(),
                }
            }
        }
    }
    if nb_args != 1 {
         display_help(args[0].to_string());
    }

    if level3 {
        // Traitement de 2 fichier fichiers recPhyloXML
        println!("Two reconciled files => displaying 3-levels reconciliations. ");
        let  outfile_gene_para = String::from("thirdkind_gene_para.svg");
        let  outfile_para_host = String::from("thirdkind_para_host.svg");
        let  outfile_mapped_1 = String::from("thirdkind_mapped_1.svg");
        let  outfile_mapped_2 = String::from("thirdkind_mapped_2.svg");
        let  outfile_mapped_3 = String::from("thirdkind_mapped_3.svg");

        let transfers = vec![]; // Initialise transfers
// ---------------------------------------------------------
// Create a structure Arena for the global parasite pipe
// tree and a vector of structures Arena for gene path trees
// ---------------------------------------------------------
let mut global_pipe_parasite: ArenaTree<String> = ArenaTree::default();
let mut global_roots: std::vec::Vec<usize> = Vec::new();
let mut path_genes: std::vec::Vec<ArenaTree<String>> = Vec::new();
// ---------------------------------------------------------
// Fill global parasite pipe tree and is roots and path
// genes trees
// ---------------------------------------------------------
read_recphyloxml_multi(infile_gs,&mut global_pipe_parasite,&mut path_genes,&mut global_roots);
let  nb_gntree =  path_genes.len().clone();
println!("Number of gene trees : {}",nb_gntree);
info!("List of gene trees : {:?}",path_genes);
let  nb_parasite_pipe =  global_roots.len().clone();
println!("Number of symbiote trees : {}",nb_parasite_pipe);
println!("List of species trees roots : {:?}",global_roots);
info!("Global symbiote pipe tree : {:?}",global_pipe_parasite);
println!();
println!("Building svg 1: reconciled pipe symbiote tree(s) with gene tree(s) [{}]",
    outfile_mapped_1);
// ---------------------------------------------------------
// Generate svg of the lobal parasite pipe tree and  path
// genes trees
// ---------------------------------------------------------
recphyloxml_processing(&mut global_pipe_parasite,&mut  path_genes, &mut options, &config,true,
        &transfers,outfile_gene_para);
// ---------------------------------------------------------
// Create a structure Arena for the host pipe tree and a
// vector of structures Arena for parasite path trees
// ---------------------------------------------------------
let mut tree_host_pipe: ArenaTree<String> = ArenaTree::default();
let mut path_para_trees:std::vec::Vec<ArenaTree<String>> = Vec::new();
// ---------------------------------------------------------
// Fill  host pipe tree and is roots and path parasite trees
// ---------------------------------------------------------
let mut global_roots: std::vec::Vec<usize> = Vec::new();
read_recphyloxml_multi(infile_sh,&mut tree_host_pipe,&mut path_para_trees, &mut global_roots);
let  nb_parasite_path =  path_para_trees.len().clone();
println!("Number of pipe symbiote trees in gene-symbiote file : {}",nb_parasite_pipe);
println!("Number of path symbiote trees in symbiote-host file : {}",nb_parasite_path);
if nb_parasite_path != nb_parasite_pipe {
    println!();
    println!("==============================================");
    println!("Error! Different number of parasite trees in the 2 files!");
    println!("       Resulting svg will be incomplete.");
    println!("==============================================");
    println!();
}
// ---------------------------------------------------------
// Generate svg of the host pipe tree and path symbiote trees
// ---------------------------------------------------------
recphyloxml_processing(&mut tree_host_pipe,&mut  path_para_trees, &mut options, &config,
    true, &transfers,outfile_para_host);
// ---------------------------------------------------------
// Generation of first 3 levels svg
// ---------------------------------------------------------
info!("Symbiote trees as a 'path tree' : {:?}",path_para_trees);
info!("Symbiote tree as a 'pipe tree' : {:?}",global_pipe_parasite);
println!("==============================================");
println!("Map symbiote as 'path' to symbiote as 'pipe'");
println!("==============================================");
let mut i = 0;
while i < nb_parasite_pipe {
        map_parasite_g2s(&mut global_pipe_parasite, &mut path_para_trees[i]);
i = i + 1;
}
info!("Global symbiote tree wih events : {:?}",global_pipe_parasite);
reset_pos(&mut global_pipe_parasite);
let mut i = 0;
while i < nb_gntree {
    reset_pos(&mut path_genes[i]);
    i = i + 1;
}
println!("==============================================");
println!("Map symbiote as 'pipe' to symbiote as 'path'");
println!("==============================================");
let mut i = 0;
while i < nb_parasite_pipe {
    map_parasite_s2g(&mut global_pipe_parasite, &mut path_para_trees[i],&mut path_genes);
    i = i +  1;
}
info!("Global pipe symbiote after mapping s2g : {:?}",global_pipe_parasite);
println!("==============================================");
println!("Map symbiote as 'path' to symbiote as 'pipe' again");
println!("==============================================");
let mut i = 0;
while i < nb_parasite_pipe {
    map_parasite_g2s(&mut global_pipe_parasite, &mut path_para_trees[i]);
    i = i + 1;
}
reset_pos(&mut global_pipe_parasite);
let mut i = 0;
while i < nb_gntree {
    reset_pos(&mut path_genes[i]);
    i = i + 1;
}
// attention on ne remape pas
recphyloxml_processing(&mut global_pipe_parasite,&mut  path_genes, &mut options, &config,false,
        &transfers,outfile_mapped_1);
let path = env::current_dir().expect("Unable to get current dir");
let url_file = format!("file:///{}/{}", path.display(),"thirdkind_mapped_1.svg".to_string());
if options.open_browser {
    if webbrowser::open_browser(Browser::Default,&url_file).is_ok() {
            info!("Browser OK");
        }
    }

// ---------------------------------------------------------
// Generation of second 3 levels svg
// ---------------------------------------------------------
println!();
println!("Building svg 2:  symbiote tree(s) within host pipe tree and mapped gene transfers [{}]",
    outfile_mapped_2);
let mut i = 0;
let gene_transfers = get_gtransfer(&mut path_genes[i]);
info!("Transfers = {:?}",gene_transfers);
let mut mapped_gene_transfers = map_transfer_mul(gene_transfers, &mut path_para_trees);
info!("Mapped transfers = {:?}",mapped_gene_transfers);
i = i + 1;
while i < nb_gntree {
    let gene_transfers = get_gtransfer(&mut path_genes[i]);
    info!("Transfers = {:?}",gene_transfers);
    let mapped = map_transfer(gene_transfers, &mut path_para_trees[0]);
    info!("Mapped transfers = {:?}",mapped);
    for val in mapped {
        mapped_gene_transfers.push(val);
    }
    i = i + 1;
}
reset_pos(&mut tree_host_pipe);
let mut i = 0;
while i < nb_parasite_pipe {
    reset_pos(&mut path_para_trees[i]);
    i = i + 1;
}
// attention on ne remape pas
recphyloxml_processing(&mut tree_host_pipe, &mut path_para_trees, &mut options, &config,
    false, &mapped_gene_transfers,outfile_mapped_2);
let path = env::current_dir().expect("Unable to get current dir");
let url_file = format!("file:///{}/{}", path.display(),"thirdkind_mapped_2.svg".to_string());
if options.open_browser {
    if webbrowser::open_browser(Browser::Default,&url_file).is_ok() {
        info!("Browser OK");
    }
}

reset_pos(&mut global_pipe_parasite);
phyloxml_processing(&mut global_pipe_parasite, &mut options, &config,"thirdkind_para_simple.svg".to_string());
reset_pos(&mut tree_host_pipe);
phyloxml_processing(&mut tree_host_pipe, &mut options, &config,"thirdkind_host_simple.svg".to_string());
let mut i = 0;
while i < nb_parasite_pipe {
    reset_pos(&mut path_para_trees[i]);
    phyloxml_processing(&mut path_para_trees[i], &mut options, &config,("thirdkind_gene_simple_".to_owned()+&i.to_string()+".svg").to_string());
    i = i + 1;
}

// ---------------------------------------------------------
// Generation of third 3 levels svg UNDER DEVELOPMENT
// --------------------------------------------------------
println!();
println!("Building svg 3: pipe host tree with gene tree(s) inside [{}]",
    outfile_mapped_3);
map_gene_host(&mut path_genes, &mut path_para_trees, &mut tree_host_pipe);
reset_pos(&mut tree_host_pipe);
let mut i = 0;
while i < nb_gntree {
    reset_pos(&mut path_genes[i]);
    i = i + 1;
}
recphyloxml_processing(&mut tree_host_pipe, &mut path_genes, &mut options, &config,
    true, &vec![],outfile_mapped_3);

let url_file = format!("file:///{}/{}", path.display(),"thirdkind_mapped_3.svg".to_string());
    if options.open_browser {
        if webbrowser::open_browser(Browser::Default,&url_file).is_ok() {
            info!("Browser OK");
        }
    }
println!("Output files:");
println!(" - thirdkind_host_simple.svg ...... 1 level:  host tree");
let mut i = 0;
while i < nb_parasite_pipe {
    println!(" - thirdkind_gene_simple_{}.svg .... 2 levels: gene tree(s)",&i);
    i = i + 1;
}
println!(" - thirdkind_symbiote_simple.svg .. 2 levels: symbiote tree(s)");
println!(" - thirdkind_gene_para.svg ........ 2 levels: pipe symbiote tree(s) with gene tree(s) inside");
println!(" - thirdkind_symbiote_host.svg .... 2 levels: pipe host tree with symbiote tree(s) inside");
println!(" - thirdkind_mapped_1.svg ........  3 levels: reconciled pipe symbiote tree(s) with gene tree(s)");
println!(" - thirdkind_mapped_2.svg ........  3 levels: symbiote-host reconciliation plus gene transfers");
println!(" - thirdkind_mapped_3.svg ........  3 levels: pipe host tree with gene tree(s) inside");

if nb_parasite_path != nb_parasite_pipe {
    println!();
    println!("==============================================");
    println!("Error! Different number of symbiote trees in the 2 files!");
    println!("       Resulting svg will be incomplete.");
    println!("==============================================");
    println!();
}

    }
    else if multiple_files {
        // get the url
        let path = env::current_dir().expect("Unable to get current dir");
        let url_file = format!("file:///{}/{}", path.display(),outfile.clone());
        println!("Multiple files processing:");
        let multifilename = &infile_sh.clone();
        let contents = fs::read_to_string(multifilename).expect("Unable to read inout file.");
        let files = contents.lines();
        // let mut sp_tree: ArenaTree<String> = ArenaTree::default();
        // Creation du vecteur de structure ArenaTree pour les genes
        // ---------------------------------------------------------
        let mut gene_trees:std::vec::Vec<ArenaTree<String>> = Vec::new();
        // Empty additional transfers
        let mut transfers = vec![];
        let mut sp_trees:std::vec::Vec<ArenaTree<String>> = Vec::new();

        for filename in files {
            println!("Processing file {}",filename);
            // On cree une structure Arena pour l'arbre d'espece
            // et un vecteur de  structures Arena pour le(s) arbres de gènes
            // -------------------------------------------------------------
            // Creation de la structure ArenaTree pour l'arbre d'espece
            // --------------------------------------------------------
            let mut _sp_tree: ArenaTree<String> = ArenaTree::default();
            // Creation du vecteur de structure ArenaTree pour les genes
            // ---------------------------------------------------------
            let mut _gene_trees:std::vec::Vec<ArenaTree<String>> = Vec::new();
            // Empty additional transfers
            // let mut transfers = vec![];
            let mut _global_roots: std::vec::Vec<usize> = Vec::new();
            read_recphyloxml_multi(filename.to_string(), &mut _sp_tree, &mut _gene_trees,
                &mut _global_roots);
            let  nb_gntree =  _gene_trees.len().clone();
            println!("Number of gene trees : {}",nb_gntree);
            info!("List of gene trees : {:?}",_gene_trees);
            gene_trees.append(&mut _gene_trees);
            sp_trees.push(_sp_tree);
        }
        let  nb_gntree =  gene_trees.len().clone();
        println!("Total number of gene trees : {}",nb_gntree);
        info!("List of all gene trees : {:?}",gene_trees);
        if options.thickness_flag {
            if options.thickness_gene > nb_gntree {
                println!("There are only {} genes in the file, unable to display gene #{}",
                nb_gntree,options.thickness_gene);
                process::exit(1);
            }
            //  Recupere les transferts
            transfers = get_gtransfer(&mut gene_trees[0]);
            let mut i = 1;
            while i < nb_gntree {
                let gene_transfer = get_gtransfer(&mut gene_trees[i]);
                for val in gene_transfer {
                    transfers.push(val);
                }
                i = i + 1;
            }
            info!("Transfers = {:?}",transfers);
            let mut selected_gene_trees:std::vec::Vec<ArenaTree<String>> = Vec::new();
            selected_gene_trees.push(gene_trees.remove(options.thickness_gene));
            recphyloxml_processing(&mut sp_trees[0], &mut selected_gene_trees, &mut options,
                &config, true, &transfers, outfile);
        }
        else {
            if options.disp_gene  > 0 {
                // On traite l'arbre de gene comme un arbre au format phylxoml
                if options.disp_gene > nb_gntree {
                    println!("There are only {} genes in the file, unable to display gene #{}",
                        nb_gntree,options.disp_gene);
                        process::exit(1);
                    }
                let  mut tree = &mut gene_trees[options.disp_gene-1];
                phyloxml_processing(&mut tree, &options, &config, outfile);
            }
            else {
                recphyloxml_processing(&mut sp_trees[0],&mut  gene_trees, &mut options,
                    &config,true, &transfers, outfile);
            }
        }
        if options.open_browser {
            if webbrowser::open_browser(Browser::Default,&url_file).is_ok() {
                info!("Browser OK");
            }
        }
    }
    // Traitement d'un fichier unique qui peu etre newick, phyloXML ou recPhyloXML
    else {
        // Determination du format
        // ------------------------
        let filename = &infile_sh.clone();
        info!("Input filename is {}",filename);
        let dot = filename.rfind('.');
        let format = match dot {
            None => Format::Newick,
            Some(dot) => {
                let suffix = &filename[dot..];
                info!("File suffix is {:?}",suffix);
                match suffix {
                    ".xml" => Format::Recphyloxml,
                    ".phyloxml" => Format::Phyloxml,
                    ".recphyloxml" => Format::Recphyloxml,
                    ".recPhyloXML" => Format::Recphyloxml,
                    ".recphylo" => Format::Recphyloxml,
                    _ => Format::Newick,
                }
            },
        };
        let format = match _format {
            Format::Newick => {
                println!("Assume that format is {:?}",format);
                format
            },
            _ => {
                println!("User defined format {:?}",_format);
                _format
            },
        };
        // get the url
        let path = env::current_dir().expect("Unable to get current dir");
        let url_file = format!("file:///{}/{}", path.display(),outfile.clone());

        // Creation d'une structure ArenaTree (pour phyloxml et newick)
        // -----------------------------------------------------------
        let mut tree: ArenaTree<String> = ArenaTree::default();
        // Charge l'arbre selon le format de fichier
        //  ----------------------------------------
        match format {
            // Phymxoml
            Format::Phyloxml => {
                read_phyloxml(filename.to_string(), &mut tree);
                phyloxml_processing(&mut tree, &options, &config, outfile);
            },
            // Newick
            Format::Newick => {
                read_newick(filename.to_string(), &mut tree);
                phyloxml_processing(&mut tree, &options, &config, outfile);
            },
            // Recphyloxml
            Format::Recphyloxml => {
                // On cree une structure Arena pour l'arbre d'espece
                // et un vecteur de  structures Arena pour le(s) arbres de gènes
                // -------------------------------------------------------------
                // Creation de la structure ArenaTree pour l'arbre d'espece
                // --------------------------------------------------------
                let mut sp_tree: ArenaTree<String> = ArenaTree::default();
                // Creation du vecteur de structure ArenaTree pour les genes
                // ---------------------------------------------------------
                let mut gene_trees:std::vec::Vec<ArenaTree<String>> = Vec::new();
                // Empty additional transfers
                let mut transfers = vec![];
                let mut global_roots: std::vec::Vec<usize> = Vec::new();
                read_recphyloxml_multi(filename.to_string(), &mut sp_tree, &mut gene_trees,
                    &mut global_roots);
                let  nb_gntree =  gene_trees.len().clone();
                println!("Number of gene trees : {}",nb_gntree);
                info!("List of gene trees : {:?}",gene_trees);
                if options.thickness_flag {
                    if options.thickness_gene > nb_gntree {
                        println!("There are only {} genes in the file, unable to display gene #{}",
                        nb_gntree,options.thickness_gene);
                        process::exit(1);
                    }
                    //  Recupere les transferts
                    transfers = get_gtransfer(&mut gene_trees[0]);
                    let mut i = 1;
                    while i < nb_gntree {
                        let gene_transfer = get_gtransfer(&mut gene_trees[i]);
                        for val in gene_transfer {
                            transfers.push(val);
                        }
                        i = i + 1;
                    }
                    info!("Transfers = {:?}",transfers);
                    let mut selected_gene_trees:std::vec::Vec<ArenaTree<String>> = Vec::new();
                    selected_gene_trees.push(gene_trees.remove(options.thickness_gene));
                    recphyloxml_processing(&mut sp_tree, &mut selected_gene_trees, &mut options,
                        &config, true, &transfers, outfile);
                    }
                    else {
                        if options.disp_gene  > 0 {
                            // On traite l'arbre de gene comme un arbre au format phylxoml
                            if options.disp_gene > nb_gntree {
                                println!("There are only {} genes in the file, unable to display gene #{}",
                                nb_gntree,options.disp_gene);
                                process::exit(1);
                            }
                            let  mut tree = &mut gene_trees[options.disp_gene-1];
                            phyloxml_processing(&mut tree, &options, &config, outfile);
                        }
                    else {
                        recphyloxml_processing(&mut sp_tree,&mut  gene_trees, &mut options,
                             &config,true, &transfers, outfile);
                        }
                    }
                },
            }
        if options.open_browser {
            if webbrowser::open_browser(Browser::Default,&url_file).is_ok() {
                info!("Browser OK");
            }
        }
    }
}

fn set_config(configfile: String, config: &mut Config) {
    let contents = fs::read_to_string(configfile)
                .expect("Something went wrong reading the config file");
    let conf = contents.split('\n');
    for line in conf {
        let test: Vec<&str> = line.split(':').collect();
        if test.len() == 2 {
            match test[0] {
                "species_color" => {
                    info!("[set_config] species_color was {}",config.species_color);
                    config.species_color=test[1].to_string();
                    info!("[set_config] species_color is now {}",config.species_color);
                },
                "species_opacity" => {
                    info!("[set_config] species_opacity was {}",config.species_opacity);
                    config.species_opacity=test[1].to_string();
                    info!("[set_config] species_opacity is now {}",config.species_opacity);
                },
                "single_gene_color" => {
                    info!("[set_config] single_gene_color was {}",config.single_gene_color);
                    config.single_gene_color=test[1].to_string();
                    info!("[set_config] single_gene_color is now {}",config.single_gene_color);
                },
                "gene_opacity" => {
                    info!("[set_config] gene_opacity was {}",config.gene_opacity);
                    config.gene_opacity=test[1].to_string();
                    info!("[set_config] gene_opacity is now {}",config.gene_opacity);
                },
                "species_police_color" => {
                    info!("[set_config] species_police_color was {}",config.species_police_color);
                    config.species_police_color=test[1].to_string();
                    info!("[set_config] species_police_color is now {}",config.species_police_color);
                },
                "species_police_size" => {
                    info!("[set_config] species_police_size was {}",config.species_police_size);
                    config.species_police_size=test[1].to_string();
                    info!("[set_config] species_police_size is now {}",config.species_police_size);
                },
                "gene_police_size" => {
                    info!("[set_config] gene_police_size was {}",config.gene_police_size);
                    config.gene_police_size=test[1].to_string();
                    info!("[set_config] gene_police_size is now {}",config.gene_police_size);
                },
                "bezier" => {
                    info!("[set_config] bezier was {}",config.bezier);
                    config.bezier=test[1].to_string();
                    info!("[set_config] bezier is now {}",config.bezier);
                },
                _ => {},
            }
        }

    }
}
