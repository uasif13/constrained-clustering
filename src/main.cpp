#include <iostream>
#include <mpi_telemetry.h>

#include "argparse.h"
#include "constrained.h"
#include "library.h"
#include "mincut_only.h"
/* #include "cm.h" */


int main(int argc, char* argv[]) {
    argparse::ArgumentParser main_program("constrained-clustering");
    /* BEGIN comment out cm */
    /* argparse::ArgumentParser cm("CM"); */
    /* cm.add_description("CM"); */
    /* END comment out cm */

    argparse::ArgumentParser mincut_only("MincutOnly");
    mincut_only.add_description("WCC");

    /* BEGIN comment out cm */
    /* cm.add_argument("--edgelist") */
    /*     .required() */
    /*     .help("Network edge-list file"); */
    /* cm.add_argument("--algorithm") */
    /*     .help("Clustering algorithm to be used (leiden-cpm, leiden-mod, louvain)") */
    /*     .action([](const std::string& value) { */
    /*         static const std::vector<std::string> choices = {"leiden-cpm", "leiden-mod", "louvain"}; */
    /*         if (std::find(choices.begin(), choices.end(), value) != choices.end()) { */
    /*             return value; */
    /*         } */
    /*         throw std::invalid_argument("--algorithm can only take in leiden-cpm, leiden-mod, or louvain."); */
    /*     }); */
    /* cm.add_argument("--resolution") */
    /*     .default_value(double(0.01)) */
    /*     .help("Resolution value for leiden-cpm. Only used if --algorithm is leiden-cpm") */
    /*     .scan<'f', double>(); */
    /* cm.add_argument("--existing-clustering") */
    /*     .default_value("") */
    /*     .help("Existing clustering file"); */
    /* cm.add_argument("--num-processors") */
    /*     .default_value(int(1)) */
    /*     .help("Number of processors") */
    /*     .scan<'d', int>(); */
    /* cm.add_argument("--output-file") */
    /*     .required() */
    /*     .help("Output clustering file"); */
    /* cm.add_argument("--log-file") */
    /*     .required() */
    /*     .help("Output log file"); */
    /* cm.add_argument("--log-level") */
    /*     .default_value(int(1)) */
    /*     .help("Log level where 0 = silent, 1 = info, 2 = verbose") */
    /*     .scan<'d', int>(); */
    /* END comment out cm */

    mincut_only.add_argument("--edgelist")
        .required()
        .help("Network edge-list file");
    mincut_only.add_argument("--existing-clustering")
        .required()
        .help("Existing clustering file");
    mincut_only.add_argument("--num-processors")
        .default_value(int(1))
        .help("Number of processors")
        .scan<'d', int>();
    mincut_only.add_argument("--output-file")
        .required()
        .help("Output clustering file");
    mincut_only.add_argument("--log-header")
        .required()
        .help("Output log header");
    mincut_only.add_argument("--connectedness-criterion")
        .default_value("1log_10(n)")
        .help("String where CC = 0, and otherwise would be in the form of Clog_x(n) or Cn^x for well-connectedness");
    mincut_only.add_argument("--log-level")
        .default_value(int(1))
        .help("Log level where 0 = silent, 1 = info, 2 = verbose")
        .scan<'d', int>();

/*         The two functional forms would be: */
/* F(n) = C log_x(n), where C and x are parameters specified by the user (our default is C=1 and x=10) */
/* G(n) = C n^x, where C and x are parameters specified by the user (here, presumably 0<x<2). Note that x=1 makes it linear. */

    /* BEGIN comment out cm */
    /* main_program.add_subparser(cm); */
    /* END comment out cm */

    main_program.add_subparser(mincut_only);
    try {
        main_program.parse_args(argc, argv);
    } catch (const std::runtime_error& err) {
        std::cerr << err.what() << std::endl;
        std::cerr << main_program;
        std::exit(1);
    }

    /* BEGIN comment out cm */
    /* if(main_program.is_subcommand_used(cm)) { */
    /*     std::string edgelist = cm.get<std::string>("--edgelist"); */
    /*     std::string algorithm = cm.get<std::string>("--algorithm"); */
    /*     double resolution = cm.get<double>("--resolution"); */
    /*     std::string existing_clustering = cm.get<std::string>("--existing-clustering"); */
    /*     int num_processors = cm.get<int>("--num-processors"); */
    /*     std::string output_file = cm.get<std::string>("--output-file"); */
    /*     std::string log_file = cm.get<std::string>("--log-file"); */
    /*     int log_level = cm.get<int>("--log-level") - 1; // so that enum is cleaner */
    /*     ConstrainedClustering* cm = new CM(edgelist, algorithm, resolution, existing_clustering, num_processors, output_file, log_file, log_level); */
    /*     random_functions::setSeed(0); */
    /*     cm->main(); */
    /*     delete cm; */
    /* } else if(main_program.is_subcommand_used(mincut_only)) { */
    /* END comment out cm */
    MPI_Init(&argc, &argv);
    int my_rank;
    int nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    uint64_t * opCount;
    opCount = new uint64_t[nprocs];
    std::string mpi_log_file;
    for (int i = 0; i < nprocs; i++) {
        opCount[i] = 0;
    }
    if(main_program.is_subcommand_used(mincut_only)) {
        std::string edgelist = mincut_only.get<std::string>("--edgelist");
        std::string existing_clustering = mincut_only.get<std::string>("--existing-clustering");
        int num_processors = mincut_only.get<int>("--num-processors");
        std::string output_file = mincut_only.get<std::string>("--output-file");
        std::string log_header = mincut_only.get<std::string>("--log-file");
        std::string log_file = log_header+"_"+std::to_string(my_rank)+".log";
        mpi_log_file = log_header+"_mpi.log";
        int log_level = mincut_only.get<int>("--log-level") - 1; // so that enum is cleaner
        std::string connectedness_criterion = mincut_only.get<std::string>("--connectedness-criterion");
        ConstrainedClustering* mincut_only = new MincutOnly(edgelist, existing_clustering, num_processors, output_file, log_file, connectedness_criterion, log_level, my_rank, nprocs);
        random_functions::setSeed(0);
        mincut_only->main(my_rank, nprocs, opCount);
        delete mincut_only;
    }
    MPI_Finalize(my_rank,nprocs, opCount, mpi_log_file);
}
