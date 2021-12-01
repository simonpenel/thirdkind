/// name = "thirdkind"
/// version = "2.4.4"
/// authors = ["Simon Penel <simon.penel@univ-lyon1.fr>"]
/// release = "18/11/2021"
/// license = "CECILL-2.1"
///
/// Usage:
/// Build svg representations of phylogenetic reconciled (or not) trees with events (loss, duplication, speciation, transfer).
/// Input one newick or phyloxml file -> a svg representation of the tree with node events
/// Input one recphyloxml file -> a svg representation of the "lower" gene (or symbiot) tree(s) inside the associated "upper" species (or host) tree
/// Input a file describing multiples recphyloxml files -> a svg representation of the "lower" gene (or symbiot) tree(s) inside the associated "upper" species (or host) tree
/// Input two nested recphyloxml files -> several svg representations allowing to display 3 level reconciliations (for example gene/symbiot/host)


use std::fs;
use std::env;
use std::process;
use getopt::Opt;
use webbrowser::{Browser};
use light_phylogeny::*;
use log::{info};

// Message d'aide court
// --------------------
fn display_usage(programe_name:String) {
    const VERSION: Option<&'static str> = option_env!("CARGO_PKG_VERSION");
    const NAME: Option<&'static str> = option_env!("CARGO_PKG_NAME");
    const DESCRIPTION: Option<&'static str> = option_env!("CARGO_PKG_DESCRIPTION");
    println!("{} v{}", NAME.unwrap_or("unknown"),VERSION.unwrap_or("unknown"));
    println!("{}", DESCRIPTION.unwrap_or("unknown"));
    println!();
    println!("Usage:");
    println!("{} -f input file [-b][-c config file][-d fontsize][-D fontsize][-e][-E][-F format][-g input file][-G #][-h]\
    [-H height][-i][-I][-J][-l factor][-L][-m][-o output file][-O][-p][-P][-r ratio][-s][-S]\
    [-t threshold][-T #][-u threshold][-U #][-v][-W width]",programe_name);
        println!();
    println!("Get help:");
    println!("{} -h ",programe_name);
    process::exit(1);
}
// Message d'aide etendu
// ---------------------
fn display_help(programe_name:String) {
    const VERSION: Option<&'static str> = option_env!("CARGO_PKG_VERSION");
    const NAME: Option<&'static str> = option_env!("CARGO_PKG_NAME");
    const DESCRIPTION: Option<&'static str> = option_env!("CARGO_PKG_DESCRIPTION");
    println!("{} v{}", NAME.unwrap_or("unknown"),VERSION.unwrap_or("unknown"));
    println!("{}", DESCRIPTION.unwrap_or("unknown"));
    println!("Usage:");
    println!("{} -f input file [-b][-c config file][-d fontsize][-D fontsize][-e][-E][-F format][-g input file][-G #][-h]\
    [-H height][-i][-I][-J][-l factor][-L][-m][-o output file][-O][-p][-P][-r ratio][-s][-S]\
    [-t threshold][-T #][-u threshold][-U #][-v][-W width]",programe_name);
    println!("    -b : open svg in browser");
    println!("    -c configfile: use a configuration file");
    println!("    -d fontsize: set font size for gene trees");
    println!("    -D fontsize: set font size for species trees");
    println!("    -e : the node associated to FREE_LIVING are drawned in an \
    external tree (free_living option) and superposed in case of multiple genes");
    println!("    -E : the node associated to FREE_LIVING are drawned in an \
    external tree (free_living option) and shifted in case of multiple genes");
    println!("    -F phylo/recphylo: force format phyloXML/recPhyloXML");
    println!("    -g 2nd level input file (for example a gene-symbiote file with -f defining a symbiote-host file)");
    println!("    -G <n> : display the gene #n in phyloxml style (no species tree)");
    println!("    -h : help");
    println!("    -H height : multiply the tree height by factor 'height'");
    println!("    -i : display internal gene nodes");
    println!("    -I : display internal species nodes");
    println!("    -J : with option -t, display the abundance of redudant transfers");
    println!("    -l factor : use branch length, multiplied by the given factor in case of newick/phyloxml,");
    println!("                and by an optimised factor in case of recphyloxml");
    println!("    -L : display as landscape");
    println!("    -m : the input file (-f) is a list of recphyloxml files");
    println!("    -o outputfile : set name of output file");
    println!("    -O switching nodes in order to minimise transfer crossings (under development) ");
    println!("    -p : build a phylogram");
    println!("    -P : species 'upper' tree uniformisation");
    println!("    -r ratio : set the ratio between width of species and gene tree");
    println!("               Default 1.0, you usualy do not need to change it");

    println!("    -s : drawing species tree only");
    println!("    -S : display node support");
    println!("    -t <t> : redudant transfers are displayed as one, with opacity according \
    to abundance and only if abundance is higher tan t\n             Only one gene is displayed");
    println!("    -T <n> : with option -t, select the gene to display");
    println!("    -u <t> : with -g, same as -t, but apply to the '-f' input file, and -t will apply to the '-g' file");
    println!("    -U <n> : same as -T with -t, but for -u");
    println!("    -v : verbose");
    println!("    -W width : multiply the tree width by factor 'width'");
    println!("");
    println!("    Note on -b option : you must set a browser as default application for opening \
    svg file");
    println!("");
    println!("    Note on -g option : this will generate 3-levels reconciliation svg files.");
    println!("    For example you may input a gene-symbiote recphyloxml file  with -g and symbiote-host recphyloxml file with -f");
    println!("    The -t/-u options are not totally implemented for the 3-levels reconciliation svg output files.");
    println!("");
    println!("Input format is guessed according to the file name extension:");

    println!(".phyloxml    => phyloXML");
    println!(".xml         => recPhyloxml");
    println!(".recphyloxml => recPhyloXML");
    println!(".recPhyloXML => recPhyloXML");
    println!(".recphylo    => recPhyloXML");
    println!("any other    => newick");
    println!("");
    println!("Home page: https://github.com/simonpenel/thirdkind/wiki");
    println!("");
    println!("About recPhyloXML format: http://phylariane.univ-lyon1.fr/recphyloxml/");
    println!("recPhyloXML paper: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6198865/");
    println!("About phyloXML format: http://www.phyloxml.org/");
    println!("phyloXML paper: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2774328/");
    println!("");
    println!("Examples with recPhyloXML files (available at https://github.com/simonpenel/thirdkind):");
    println!("{} -f recphylo_examples/FAM000297_reconciliated.recphylo  -b", programe_name);
    println!("{} -f recphylo_examples/concat.xml -b -t 0 ", programe_name);
    println!("{} -f recphylo_examples/hote_parasite_page4_BL.recphylo  -b -l 1", programe_name);
    println!("{} -f recphylo_examples/testfiles -m -b -t 3 -J", programe_name);
    println!("{} -f paramecium_data/liste.txt -m -b -t 25 -J", programe_name);
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
         display_usage(args[0].to_string());
    }
    let mut opts = getopt::Parser::new(&args, "c:bd:D:eEf:F:g:G:hH:iIJl:Lmo:OpPr:sSt:T:u:U:vW:");
    let mut infile_sh = String::new(); // symbiote host file
    let mut infile_gs = String::new(); // gene symbiote file
    let mut outfile = String::from("thirdkind.svg");
    let mut nb_args = 0;
    let mut level3 = false; // Affichage à 3 niveaux
    let mut multiple_files = false;
    let mut _format = Format::Newick;
    // 1st level file option
    let mut thickness_thresh_1st = 0;
    let mut thickness_gene_1st = 1;
    let mut thickness_flag_1st = false;
    // 2nd level file option
    let mut thickness_thresh_2nd = 0;
    let mut thickness_gene_2nd = 1;
    let mut thickness_flag_2nd = false;
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
                    Opt('E', None) => {
                        options.free_living = true;
                        options.free_living_shift = true;
                        },
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
                    Opt('P', None) => options.uniform = true,
                    Opt('s', None) => options.species_only_flag = true,
                    Opt('S', None) => options.support = true,
                    Opt('t', Some(string)) => {
                        thickness_thresh_1st = match string.parse::<usize>(){
                            Ok(valeur) => valeur,
                            Err(_err) => {
                                eprintln!("Error! Please give a integer value with -t option");
                                process::exit(1);
                            },
                        };
                        thickness_flag_1st = true;
                    },
                    Opt('T', Some(string)) => {
                        thickness_gene_1st = match string.parse::<usize>(){
                            Ok(valeur) => valeur,
                            Err(_err) => {
                                eprintln!("Error! Please give a integer value with -T option");
                                process::exit(1);
                            },
                        };
                    },
                    Opt('u', Some(string)) => {
                        thickness_thresh_2nd = match string.parse::<usize>(){
                            Ok(valeur) => valeur,
                            Err(_err) => {
                                eprintln!("Error! Please give a integer value with -u option");
                                process::exit(1);
                            },
                        };
                        thickness_flag_2nd = true;
                    },
                    Opt('U', Some(string)) => {
                        thickness_gene_2nd = match string.parse::<usize>(){
                            Ok(valeur) => valeur,
                            Err(_err) => {
                                eprintln!("Error! Please give a integer value with -U option");
                                process::exit(1);
                            },
                        };
                    },
                    Opt('l', Some(string)) => {
                        options.real_length_flag = true;
                        options.uniform = false; // In case we deal with a recphyloxml
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
                    Opt('d', Some(string)) => {
                        config.gene_police_size = match string.parse::<usize>(){
                            Ok(valeur) => valeur.to_string(),
                            Err(_err) => {
                                eprintln!("Error! Please give a integer value with -d option");
                                process::exit(1);
                            },
                        };
                    },
                    Opt('D', Some(string)) => {
                        config.species_police_size = match string.parse::<usize>(){
                            Ok(valeur) => valeur.to_string(),
                            Err(_err) => {
                                eprintln!("Error! Please give a integer value with -D option");
                                process::exit(1);
                            },
                        };
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
         display_usage(args[0].to_string());
    }

    // Setting options on thickness
    options.thickness_flag = thickness_flag_1st;
    options.thickness_gene = thickness_gene_1st;
    options.thickness_thresh = thickness_thresh_1st;

// ================================================================================================
//          RECONCILIATION A 3 NIVEAUX
// ================================================================================================
    if level3 {
        // Traitement de 2 fichier fichiers recPhyloXML
        println!("Two reconciled files => displaying 3-levels reconciliations. ");
        let  outfile_gene_para = String::from("thirdkind_gene_symbiote.svg");
        let  outfile_para_host = String::from("thirdkind_symbiote_host.svg");
        let  outfile_mapped_1 = String::from("thirdkind_mapped_1.svg");
        let  outfile_mapped_2 = String::from("thirdkind_mapped_2.svg");
        let  outfile_mapped_3 = String::from("thirdkind_mapped_3.svg");
        let transfers = vec![]; // Initialise transfers
        let mut transfers_gene = vec![]; // Transferts de genes
        let mut transfers_para = vec![]; // Transferts de parasites(ou symbiotes)

        // Gestion de l'optiopn free_living dans le cas d'une
        let free_living_3l = options.free_living;
        options.free_living = false;
        //
        // GENE - PARASITE STUFF
        //
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
        println!("\nBuilding 'lower' gene vs 'upper' symbiote reconciliation svg file [{}]",outfile_gene_para.clone());
        read_recphyloxml_multi(infile_gs,&mut global_pipe_parasite,&mut path_genes,&mut global_roots);
        let  nb_gntree =  path_genes.len().clone();
        println!("Number of 'lower' gene trees : {}",nb_gntree);
        info!("List of gene trees : {:?}",path_genes);
        let  nb_parasite_pipe =  global_roots.len().clone();
        println!("Number of 'upper' symbiote trees : {}",nb_parasite_pipe);
        println!("List of 'upper' symbiote tree roots : {:?}",global_roots);
        info!("Global symbiote pipe tree : {:?}",global_pipe_parasite);

        // ---------------------------------------------------------
        // Generate svg of the global parasite pipe tree and  path
        // genes trees (outfile_gene_para)
        // ---------------------------------------------------------

        // If  the option -t is on :
        if thickness_flag_1st {
            // check that gene nb is correct
            if thickness_gene_1st > nb_gntree {
                println!("There are only {} genes in the file, unable to display gene #{}",
                nb_gntree,thickness_gene_1st);
                process::exit(1);
            }

            //  Get the transfers in the genes
            transfers_gene = get_gtransfer(&mut path_genes[0]);
            let mut i = 1;
            while i < nb_gntree {
                let gene_transfer = get_gtransfer(&mut path_genes[i]);
                for val in gene_transfer {
                    transfers_gene.push(val);
                }
                i = i + 1;
            }
            println!("Transfers (genes) = {:?}",transfers_gene);
            // Define the unique gene wich is selected
            let mut selected_gene_trees:std::vec::Vec<ArenaTree<String>> = Vec::new();
            selected_gene_trees.push(path_genes[thickness_gene_1st-1].copie());

            // Define a temporary copy of the species tree
            let mut _global_pipe_parasite = global_pipe_parasite.copie();
            // Set options
            options.thickness_flag = true;
            options.thickness_gene = thickness_gene_1st;
            options.thickness_thresh = thickness_thresh_1st;
            //  Create the svg from temporary variables
            recphyloxml_processing(&mut _global_pipe_parasite, &mut selected_gene_trees, &mut options,
                &config, true, &transfers_gene, outfile_gene_para.clone());
            // We need to run again with the current variables, because we need the upddated
            // global_pipe_parasite path_genes
            recphyloxml_processing(&mut global_pipe_parasite,&mut  path_genes, &mut options, &config,true,
                &transfers,"tmpfile.svg".to_string());
        }
        // No -t option
        else {
            recphyloxml_processing(&mut global_pipe_parasite,&mut  path_genes, &mut options, &config,true,
                &transfers,outfile_gene_para);
        }

        //
        // PARASITE-HOST  STUFF
        //
        // ---------------------------------------------------------
        // Create a structure Arena for the host pipe tree and a
        // vector of structures Arena for parasite path trees
        // ---------------------------------------------------------
        let mut tree_host_pipe: ArenaTree<String> = ArenaTree::default();
        let mut path_para_trees:std::vec::Vec<ArenaTree<String>> = Vec::new();
        println!("\nBuilding 'lower' symbiote vs 'upper' host reconciliation svg file [{}]",outfile_para_host.clone());

        // ---------------------------------------------------------
        // Fill  host pipe tree and is roots and path parasite trees
        // ---------------------------------------------------------
        let mut global_roots: std::vec::Vec<usize> = Vec::new();
        read_recphyloxml_multi(infile_sh,&mut tree_host_pipe,&mut path_para_trees, &mut global_roots);
        let  nb_parasite_path =  path_para_trees.len().clone();
        let  nb_hosts_pipe = global_roots.len();
        println!("Number of 'upper' symbiote trees in gene-symbiote file : {}",nb_parasite_pipe);
        println!("Number of 'lower' symbiote trees in symbiote-host file : {}",nb_parasite_path);
        println!("Number of 'upper' host trees in symbiote-host file : {}",nb_hosts_pipe);
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

        // Reset the option
        options.thickness_flag = false;
        options.free_living = free_living_3l;

        // If  the option -u is on :
        if thickness_flag_2nd {
            options.thickness_flag = true;
            options.thickness_gene = thickness_gene_2nd;
            options.thickness_thresh = thickness_thresh_2nd;

            // check that the number pf the parasite is correct
            if options.thickness_gene > nb_parasite_path {
                println!("There are only {} parasites in the file, unable to display gene #{}",
                nb_parasite_path,options.thickness_gene);
                process::exit(1);
            }

            // Get teh transfers in the parasites
            transfers_para = get_gtransfer(&mut path_para_trees[0]);
            let mut i = 1;
            while i < nb_parasite_path {
                let gene_transfer = get_gtransfer(&mut path_para_trees[i]);
                for val in gene_transfer {
                    transfers_para.push(val);
                }
            i = i + 1;
            }
            println!("Transfers (parasites) = {:?}",transfers_para);
            // Define a temporary copy of the paraistes
            let mut _path_para_trees: std::vec::Vec<ArenaTree<String>> = Vec::new();
            for i in 0 .. path_para_trees.len() {
                _path_para_trees.push(path_para_trees[i].copie());
            }
            // Define the unique parasite  wich is selected
            let mut selected_para_trees:std::vec::Vec<ArenaTree<String>> = Vec::new();
            selected_para_trees.push(_path_para_trees.remove(options.thickness_gene-1));
            // Define a tmprary copy of the host
            let mut _tree_host_pipe = tree_host_pipe.copie();

            //  Create the svg from the temprary variables
            recphyloxml_processing(&mut _tree_host_pipe, &mut selected_para_trees, &mut options,
            &config, true, &transfers_para, outfile_para_host.clone());
            // We need to run again with the current variables, because we need the upddated
            // _tree_host_pipe path_para_trees
            recphyloxml_processing(&mut tree_host_pipe,&mut  path_para_trees, &mut options,
                 &config,true, &transfers,"tmpfile2.svg".to_string());

        }
        // No -u option
        else {
            recphyloxml_processing(&mut tree_host_pipe,&mut  path_para_trees, &mut options,
                &config, true, &transfers,outfile_para_host);
        }

        //
        // GENE-PARASITE-HOST  STUFF : MAP1
        //
        // ---------------------------------------------------------
        // Generation of first 3 levels svg
        // ---------------------------------------------------------
        println!("\nBuilding 'mapped 1': reconciled 'upper' symbiote tree(s) with 'lower' gene tree(s) [{}]",
            outfile_mapped_1);
        info!("Symbiote trees as a 'lower tree' : {:?}",path_para_trees);
        info!("Symbiote tree as a 'upper tree' : {:?}",global_pipe_parasite);
        println!("Map symbiote as 'lower' to symbiote as 'upper'");
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

        //  We need to handel the -t option
        //  C'est complexe il faut eliminer les noeuds associer aux genes qui ne sont pas selectionnes
        //  A developper
        // let mut _path_genes: std::vec::Vec<ArenaTree<String>> = Vec::new();
        // let mut selected_gene_trees:std::vec::Vec<ArenaTree<String>> = Vec::new();
        // let mut _global_pipe_parasite: ArenaTree<String> = ArenaTree::default();
        // if thickness_flag_1st {
        //     for i in 0 .. path_genes.len() {
        //         _path_genes.push(path_genes[i].copie());
        //     }
        // // Define the unique gene wich is selected
        // selected_gene_trees.push(_path_genes.remove(thickness_gene_1st-1));
        // // Define a temporary copy of the species tree
        // _global_pipe_parasite = global_pipe_parasite.copie();
        // }
        println!("Map symbiote as 'upper' to symbiote as 'lower'");
        // if thickness_flag_1st {
        //     let mut i = 0;
        //     while i < nb_parasite_pipe {
        //         map_parasite_s2g(&mut _global_pipe_parasite, &mut path_para_trees[i], &mut selected_gene_trees);
        //         i = i +  1;
        //     }
        // }
        // else {
        let mut i = 0;
        while i < nb_parasite_pipe {
            map_parasite_s2g(&mut global_pipe_parasite, &mut path_para_trees[i],&mut path_genes);
            i = i +  1;
        }
        // }
        info!("Global upper symbiote tree after mapping s2g : {:?}",global_pipe_parasite);
        println!("Map symbiote as 'lower' to symbiote as 'upper' again");
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
        // Reset the option
        options.thickness_flag = false;
        options.free_living = false;
        //  option -t is on
        if thickness_flag_1st {
            options.thickness_flag = true;
            options.thickness_gene = thickness_gene_1st;
            options.thickness_thresh = thickness_thresh_1st;
        }
        // attention on ne remape pas
        recphyloxml_processing(&mut global_pipe_parasite,&mut  path_genes, &mut options, &config,false,
            &transfers_gene,outfile_mapped_1);
        let path = env::current_dir().expect("Unable to get current dir");
        let url_file = format!("file:///{}/{}", path.display(),"thirdkind_mapped_1.svg".to_string());
        if options.open_browser {
            if webbrowser::open_browser(Browser::Default,&url_file).is_ok() {
                info!("Browser OK");
            }
        }
        //
        // GENE-PARASITE-HOST  STUFF : MAP2
        //
        // ---------------------------------------------------------
        // Generation of second 3 levels svg
        // ---------------------------------------------------------
        println!();
        println!("Building 'mapped 2':  'lower' symbiote tree(s) within 'upper' host tree and mapped gene transfers [{}]",
            outfile_mapped_2);
        let mut i = 0;
        // We get the gene transfer here again, but they will be mapped
        let gene_transfers = get_gtransfer(&mut path_genes[i]);
        info!("Transfers = {:?}",gene_transfers);
        let mut mapped_gene_transfers = map_transfer_mul(gene_transfers, &mut path_para_trees);
        info!("Mapped transfers = {:?}",mapped_gene_transfers);
        i = i + 1;
        while i < nb_gntree {
            let gene_transfers = get_gtransfer(&mut path_genes[i]);
            info!("Transfers = {:?}",gene_transfers);
            let mapped = map_transfer_mul(gene_transfers, &mut path_para_trees);
            info!("Mapped transfers = {:?}",mapped);
            for val in mapped {
                mapped_gene_transfers.push(val);
            }
            i = i + 1;
        }
        info!("Mapped transfers = {:?}",mapped_gene_transfers);
        // Reseting the pipe and the paths
        reset_pos(&mut tree_host_pipe);
        let mut i = 0;
        while i < nb_parasite_pipe {
            reset_pos(&mut path_para_trees[i]);
            i = i + 1;
        }
        // Setting option to false because we dont want the parasite transerst to be hidden
        options.thickness_flag = false;

        if thickness_flag_2nd {
            // On ajoute les transferts de paraistes au transferts de gene
            options.thickness_flag = true;
            options.thickness_gene = thickness_gene_2nd;
            options.thickness_thresh = thickness_thresh_2nd;
            //  Create the svg from temporary variables
            for val in transfers_para {
                mapped_gene_transfers.push(val);
            }
        }
        // Reset the option
        options.free_living = free_living_3l;
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
        println!("\nBuilding 'phyloxml style' svg files...");
        //  Simple tree of the parasite
        reset_pos(&mut global_pipe_parasite);
        phyloxml_processing(&mut global_pipe_parasite, &mut options, &config,"thirdkind_symbiote_simple.svg".to_string());
        reset_pos(&mut tree_host_pipe);
        //  Simple tree of the host
        phyloxml_processing(&mut tree_host_pipe, &mut options, &config,"thirdkind_host_simple.svg".to_string());
        //  Simple trees of the genes
        let mut i = 0;
        while i < nb_parasite_pipe {
            reset_pos(&mut path_para_trees[i]);
            phyloxml_processing(&mut path_para_trees[i], &mut options, &config,("thirdkind_gene_simple_".to_owned()+&i.to_string()+".svg").to_string());
            i = i + 1;
        }
        //
        // GENE-PARASITE-HOST  STUFF : MAP3
        //
        // ---------------------------------------------------------
        // Generation of third 3 levels svg UNDER DEVELOPMENT
        // --------------------------------------------------------
        println!();
        println!("Building 'mapped 3': 'upper' host tree with gene tree(s) inside [{}]",outfile_mapped_3);
        map_gene_host(&mut path_genes, &mut path_para_trees, &mut tree_host_pipe);
        reset_pos(&mut tree_host_pipe);
        let mut i = 0;
        while i < nb_gntree {
            reset_pos(&mut path_genes[i]);
            i = i + 1;
        }
        // Reset the option
        options.free_living = false;
        recphyloxml_processing(&mut tree_host_pipe, &mut path_genes, &mut options, &config,
            true, &vec![],outfile_mapped_3);
        let url_file = format!("file:///{}/{}", path.display(),"thirdkind_mapped_3.svg".to_string());
        if options.open_browser {
            if webbrowser::open_browser(Browser::Default,&url_file).is_ok() {
                info!("Browser OK");
            }
        }
        println!("\nOutput summary:");
        println!(" - thirdkind_host_simple.svg ...... 1 level:  host tree");
        let mut i = 0;
        while i < nb_parasite_pipe {
            println!(" - thirdkind_gene_simple_{}.svg .... 2 levels: gene tree(s)",&i);
            i = i + 1;
        }
        println!(" - thirdkind_symbiote_simple.svg .. 2 levels: symbiote tree(s)");
        println!(" - thirdkind_gene_symbiote.svg .... 2 levels: 'upper' symbiote tree(s) with 'lower' gene tree(s) inside");
        println!(" - thirdkind_symbiote_host.svg .... 2 levels: 'upper' host tree with 'lower' symbiote tree(s) inside");
        println!(" - thirdkind_mapped_1.svg ........  3 levels: reconciled 'upper' symbiote tree(s) with 'lower' gene tree(s) inside");
        println!(" - thirdkind_mapped_2.svg ........  3 levels: 'upper' host tree with 'lower' symbiote tree(s) inside plus gene transfers");
        println!(" - thirdkind_mapped_3.svg ........  3 levels: 'upper' host tree with gene tree(s) inside");
        if nb_parasite_path != nb_parasite_pipe {
            println!();
            println!("==============================================");
            println!("Error! Different number of symbiote trees in the 2 files!");
            println!("       Resulting svg will be incomplete.");
            println!("==============================================");
            println!();
        }
    }

    // =================================
    //  RECONCILIATION A 2 DEUX NIVEAUX
    // =================================

    else if multiple_files {
        if options.real_length_flag {
            println!("Note: when using real length option with recPhyloXML, the scaling of branches is automatic.");
        };
        // get the url
        let path = env::current_dir().expect("Unable to get current dir");
        let url_file = format!("file:///{}/{}", path.display(),outfile.clone());
        println!("Multiple files processing:");
        let multifilename = &infile_sh.clone();
        let contents = fs::read_to_string(multifilename);
        let contents = match contents {
            Ok(contents) => contents,
            Err(err) => {
                eprintln!("Something went wrong when reading the  input file {}.\n{}",
                    multifilename,err);
                eprintln!("Please check file name and path.");
                process::exit(1);
            }
        };
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
            selected_gene_trees.push(gene_trees.remove(options.thickness_gene-1));
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
                if options.real_length_flag {
                    println!("Note: when using real length option with recPhyloXML, the scaling of branches is automatic.");
                }
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
                    selected_gene_trees.push(gene_trees.remove(options.thickness_gene-1));
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
