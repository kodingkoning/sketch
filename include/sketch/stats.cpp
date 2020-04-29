#include <math.h>
#include <boost/math/special_functions/binomial.hpp>

// code from https://stackoverflow.com/questions/9330915/number-of-combinations-n-choose-r-in-c
long long unsigned nChoosek( long long unsigned n, long long unsigned k )
{
    if (k > n) return 0;
    if (k * 2 > n) k = n-k;
    if (k == 0) return 1;

    long long unsigned result = n;
    for( int i = 2; i <= k; ++i ) {
        result *= (n-i+1);
        result /= i;
    }
    return result;
}

// code adapted from https://stackoverflow.com/questions/9330915/number-of-combinations-n-choose-r-in-c
// with assistance from Joyce Chew
// currently requires debugging for accuracy
long double nChoosek_logscale( long long unsigned n, long long unsigned k )
{
    long double log_n = log(n);
    // if (k > n) return 0; // by definition cannot be true in this application, so not adapted for log scale
    if (k * 2 > n) k = n-k;
    if (k == 0) return 0;

    long double result = log_n;
    for( int i = 2; i <= k; ++i ) {
        result += log(n-i+1);
        result -= log(i);
    }
    return result;
}

// based on answer here: https://math.stackexchange.com/questions/72223/finding-expected-number-of-distinct-values-selected-from-a-set-of-integers
long double expected_unique(long long unsigned kmers, long long unsigned n) {
    long double unique_kmers = 0;
    for(long long i = 0; i < kmers-1; ++i) {
        // a = (-1)^i
        double a = (i % 2 == 0) ? 1.0 : -1.0; // more efficient than using pow()
        // b = kmers choose (i + 1)
        // long long unsigned b = nChoosek(kmers, i+1);
        long double b = boost::math::binomial_coefficient<long double>(kmers, i+1);

        // long double log_b = nChoosek_logscale(kmers, i+1);
        // assert(log_b == log(b));
        // long double log_b = log(boost::math::binomial_coefficient<long double>(kmers, i+1));
        // long double log_c = -1.0*n*i;
        // c = n^i
        long double c = pow(n, i);
        // unique_kmers += a*exp(log_b+log_c);

        unique_kmers += a*b / c;
        fprintf(stderr, "n = %lld and i = %lld ==> ", n, i);
        fprintf(stderr, "\ta = %f, ", a);
        fprintf(stderr, "\tb = %lf, ", b);
        fprintf(stderr, "\tc = %f\n", c);
        fprintf(stderr, "n = %lld and i = %lld ==> \ta = %f, \tb = %lf, \tc = %f\n", n, i, a, b, c);
    }
    return unique_kmers;
}

int test_one_read() {
    // constants
    int k = 7; // k = 21 is the default for Mash

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
    long long unsigned n = pow(4, k);
    long double unique_kmers = expected_unique(kmers, n);

    fprintf(stderr, "Estimated number of unique kmers with k of %d is %llf.\n", k, unique_kmers);
    test_nChoosek_logscale();
    fprintf(stderr, "estimated number of unique kmers with k of %d and n of %d = %llf\n", 3, 20, expected_unique(3, 20));
    fprintf(stderr, "estimated number of unique kmers with k of %d and n of %d = %llf\n", 3, 10, expected_unique(3, 10));
    fprintf(stderr, "estimated number of unique kmers with k of %d and n of %d = %llf\n", 3, 5, expected_unique(3, 5));

    // TODO: consider header line (at least for testing purposes)

    // std::fprintf(stderr, "The header line may be more than 1\% of the file size, so the estimated k-mers may be overestimated.");
}
