#include <mpi.h>
#include "mh.h"
#include "include/sketch/BigInt/BigInt.h"
#include <iostream>
#include "mpiParallelIO.cpp"
#include "calcThreshold.cpp"

// TODO: eventually this whole file should go, move the code other places

int test_one_read() {
    // constants
    int k = 7; // k = 21 is the default for Mash
    std::string fastq_file = "ecsample1.fastq"; // file name of the fastq file in the current directory
    std::string dir = get_current_dir_name();
    std::string filename = dir+"/" + fastq_file;

    // initMPI();

    BigInt threshold = find_threshold(filename, k, SKETCH_SIZE);
    RangeMinHash<std::string> localSketch(LOCAL_SKETCH_SIZE);
    // RangeMinHash<std::string> globalSketch(SKETCH_SIZE);

    // Read in file
    readFile(filename.c_str(), k, localSketch);
    std::cout << localSketch.size() << std::endl;

    // collect local sketches into the global sketch
    // RangeMinHash<std::string> globalSketch = combineSketches(localSketch);

    // TODO: load in all of the local sketches
    // finalizeMPI();
}

int test_expected_unique() {
    // tests that had been used to verify results of expected_unique() and other calculations
        // constants
    int k = 7; // k = 21 is the default for Mash
    // TODO: define k somewhere better

    // size of file
    std::string dir = get_current_dir_name();
    std::fprintf(stderr, "current dir = %s\n", dir.c_str());
    std::string filename = dir+"/ecsample1.fastq";
    std::fprintf(stderr, "file = %s\n", filename.c_str());
    struct stat sb;
    if(stat(filename.c_str(), &sb) == -1) {
        perror("stat");
        exit(EXIT_FAILURE);
    }
    long long size = sb.st_size; // size is in bytes, which is equal to chars
    long long bases = size / 2;
    std::fprintf(stderr, "Size of %s is %lld bytes, or about %lld bases\n", filename.c_str(), size, bases);

    // number of k-mers
    long long kmers = bases - k + 1;
    // this is assuming that the characters used for the names of the reads cancel
    //      out for the k-mers lost at the ends of the reads.
    // TODO: make this more sophisticated

    // estimated number of unique k-mers
    // n is the number of possible k-mers
    BigInt n = BigInt(pow(4, k));
    fprintf(stderr, "estimated number of unique kmers with k of %d and n of %d = %s\n", 3, 20, expected_unique(3, 20).ToString().c_str());
    fprintf(stderr, "estimated number of unique kmers with k of %d and n of %d = %s\n", 3, 10, expected_unique(3, 10).ToString().c_str());
    fprintf(stderr, "estimated number of unique kmers with k of %d and n of %d = %s\n", 3, 10, expected_unique(3, 5).ToString().c_str());
    fprintf(stderr, "estimated number of unique kmers with k of %d and n of %d = %s\n", 20, 100, expected_unique(20, 100).ToString().c_str());
    fprintf(stderr, "double estimated number of unique kmers with k of %d and n of %d = %f\n", 3, 20, expected_unique_dbl(3, 20));
    fprintf(stderr, "double estimated number of unique kmers with k of %d and n of %d = %f\n", 3, 10, expected_unique_dbl(3, 10));
    fprintf(stderr, "double estimated number of unique kmers with k of %d and n of %d = %f\n", 3, 10, expected_unique_dbl(3, 5));
    fprintf(stderr, "double estimated number of unique kmers with k of %d and n of %d = %f\n", 20, 100, expected_unique_dbl(20, 100));

    // BigInt unique_kmers = expected_unique(kmers, n);
    double n_double = pow(4, k);
    BigInt unique_kmers = expected_unique(kmers, n);
    std::cout << "n = " << n_double << " kmers = " << kmers << " unique_kmers = " << unique_kmers << std::endl;
    std::cout << "n = " << n_double << " kmers = " << kmers << " unique_kmers = " << expected_unique_dbl(kmers, n_double) << std::endl;
    fprintf(stderr, "double estimated number of unique kmers with k of %f and n of %f = %f\n", kmers, n_double, expected_unique_dbl(kmers, n_double));

    // fprintf(stderr, "Estimated number of unique kmers with k of %d and n of %d = %d\n", kmers, n_double, unique_kmers);

    // fprintf(stderr, "Estimated number of unique kmers with k of %d is %s.\n", k, unique_kmers.ToString().c_str());

    // TODO: consider header line (at least for testing purposes)

    // std::fprintf(stderr, "The header line may be more than 1\% of the file size, so the estimated k-mers may be overestimated.");

    RangeMinHash<std::string> localSketch(150);
    // NOTE: this assumes that max = 2^64 and min = 0
    if (localSketch.sketch_size() > unique_kmers) {
        std::cout << "NOTE: the size of the sketch is greater than the number of expected unique kmers." << std::endl;
    }
    BigInt threshold = localSketch.sketch_size() / unique_kmers * pow(2, 64) + 1 ;// (std::numeric_limits<uint64_t>::max) - (std::numeric_limits<uint64_t>::min) + 1;

    // Read in file
    readFile(filename.c_str(), k, localSketch);
    std::cout << localSketch.size() << std::endl;
}
