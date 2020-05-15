#include <mpi.h>
#include "mh.h"
#include "include/sketch/BigInt/BigInt.h"
#include <iostream>
#include "mpiParallelIO.cpp"
#include <sys/stat.h>
#include <math.h>

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

    // BigInt result_a = n - n*power(n-1, kmers)/power(n, kmers);
    return n - n*power(n-1, kmers)/power(n, kmers); // tested

    // std::cout << "kmers = " << kmers << " and n = " << n << std::endl;
    // new strategy
    // BigInt result_a = n * (1 - power((n-1)/n, kmers));
    // std::cout << "Result a of expected_unique = " << result_a << std::endl;
    // return result_a;

    // Strategy produces result
    // while this looks longer (because of the for loop) than the other two, it is actually better because 
    // the one-line versions have two power() calls. 
    // However, this will take longer as we increase the file size and number of kmers
    // TODO: find a way to parallelize???
    // TODO: time this on Cori.
    BigInt result = 1;
    for(BigInt i = 0; i < kmers-1; i += 1) {
        result = 1 + (1 - 1/n)*result;
    }
    std::cout << "Result b of expected_unique = " << result << std::endl;
    return result;

    BigInt unique_kmers = 0;
    BigInt a = 1;
    BigInt b = kmers;
    BigInt c = 1;

    std::cout << "kmers = " << kmers << " and n = " << n << std::endl;

    for(BigInt i = 0; i < kmers-1; i += 1) {
        // a = (-1)^i
        // BigInt a = (i % 2 == 0) ? 1.0 : -1.0; // more efficient than using pow()
        // b = kmers choose (i + 1)
        // BigInt b = nChoosek(kmers, i+1);
        // c = 1 / n^i
        // c *= n;
        // BigInt c = power(n, i);

        unique_kmers += a*b / c;

        // next n choose k
        a = a == -1 ? 1 : -1;
        b *= (kmers-i-1) / (i+2);
        c *= n;

        // std::cout << "a = " << a << " b = " << b << " c = " << c << std::endl;
        // fprintf(stderr, "n = %lld and i = %lld ==> a = %f, b = %lf, c = %lf\n", n, i, a, b, c);
    }
    return unique_kmers;
}

int test_one_read() {
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
