#include <sys/stat.h>
#include <math.h>
#include "mh.h"
#include "include/sketch/BigInt/BigInt.h"
#include <iostream>

const int BASES = 4;
const int SKETCH_SIZE = 150;
const int LOCAL_SKETCH_SIZE = 150;

double expected_unique_dbl(double k, double n) {
    // kmers is the number of picks, which is the number of kmers in the files
    // n is the number of unique kmers that could be chosen

    // return n - pow(n-1, k)*pow(n, 1-k);
    return n - n*pow(n-1, k)/pow(n, k);
    // return n * (1 - pow((n-1)/n, kmers));

    double result = 1;
    for(double i = 0; i < k-1; i += 1) {
        result = 1 + (1 - 1/n)*result;
    }
    return result;
}

// based on answer here: https://math.stackexchange.com/questions/72223/finding-expected-number-of-distinct-values-selected-from-a-set-of-integers
BigInt expected_unique(const BigInt& kmers, const BigInt& n) {
    // kmers is the number of picks, which is the number of kmers in the files
    // n is the number of unique kmers that could be chosen
    return n - n*power(n-1, kmers)/power(n, kmers); // tested
}

BigInt find_threshold(std::string file, int k, int sketch_size) {
    // size of file
    struct stat sb;
    if(stat(file.c_str(), &sb) == -1) {
        perror("stat");
        exit(EXIT_FAILURE);
    }
    long long size = sb.st_size; // size is in bytes, which is equal to chars
    long long bases = size / 2;
    std::fprintf(stderr, "Size of %s is %lld bytes, or about %lld bases\n", file.c_str(), size, bases);

    // number of k-mers in the file. this is assuming that the characters used
    //      for the names of the reads cancels out for the k-mers lost at the
    //      ends of the reads. Future work could improve this calculation.
    long long kmers = bases - k + 1;

    // n is the number of possible k-mers with k bases
    BigInt n = BigInt(pow(BASES, k));

    // number of unique k-mers expected in the file
    BigInt unique_kmers = expected_unique(kmers, n);

    // TODO: do something about the fact that expected_unique takes forever!

    // NOTE: this assumes that max = 2^64 and min = 0, 
    //      which is the case for WangHash used in sketch.
    if (sketch_size > unique_kmers) {
        std::cout << "NOTICE: the size of the sketch is greater than the number of expected unique kmers." << std::endl;
    }
    BigInt threshold = sketch_size / unique_kmers * pow(2, 64) + 1 ;// (std::numeric_limits<uint64_t>::max) - (std::numeric_limits<uint64_t>::min) + 1;

    return threshold;
}
