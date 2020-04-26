#include "bcl_mh.h"
#include <random>
#include <mpi.h>
#include "bcl/bcl.hpp"
#include <sys/stat.h>
#include <math.h>
#include <boost/math/special_functions/binomial.hpp>

template<typename T>
int pc(const T &x, const char *s="unspecified") {
    auto it = std::begin(x);
    std::string t = std::to_string(*it);
    while(++it != x.end()) {
        t += ',';
        t += std::to_string(*it);
    }
    std::fprintf(stderr, "Container %s contains: %s\n", s, t.data());
    return 0;
}

template<typename T>
bool forward_sorted(T i1, T i2) {
    auto tmp = *i1++;
    int ret = -1;
    start:
    if(tmp == *i1) throw std::runtime_error("e1");
    if(tmp < *i1) {if(ret >= 0 && ret != 1) throw std::runtime_error("e2"); ret = 1;}
    else  {if(ret >= 0 && ret != 0) throw std::runtime_error("e3"); ret = 0;}
    if(++i1 != i2) goto start;
    if(ret < 0) throw std::runtime_error("e4");
    return ret;
}
using namespace sketch;
using namespace common;
// using namespace mh;
using namespace bcl_mh;


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

int test_nChoosek_logscale() {
    // test cases

    long double result_log = nChoosek_logscale(10, 5);
    double result_ans = exp(result_log);
    long long unsigned result = nChoosek(10, 5);
    double answer = 252;
    fprintf(stderr, "Test 1: 10 choose 5, true value = %f\n", answer);
    fprintf(stderr, "\tnChoosek = %llu\n", result);
    fprintf(stderr, "\tnChoosek_logscale = %f\n\n", result_ans);
    // fprintf(stderr, "test 1: nChoosek_logscale = %f, nChoosek = %llu, answer = %f\n\n", result_ans, result, answer);

    return 0;
}

// based on answer here: https://math.stackexchange.com/questions/72223/finding-expected-number-of-distinct-values-selected-from-a-set-of-integers
long double expected_unique(long long unsigned kmers, long long unsigned n) {
    long double unique_kmers = 0;
    for(long long i = 0; i < kmers-1; ++i) {
        // a = (-1)^i
        double a = (i % 2 == 0) ? 1.0 : -1.0; // more efficient than using pow()
        // b = kmers choose (i + 1)
        // long long unsigned b = nChoosek(kmers, i+1);
        // long double log_b = nChoosek_logscale(kmers, i+1);
        unsigned log_b = log(boost::math::binomial_coefficient<long double>(kmers, i+1));
        long double log_c = -1.0*n*i;
        // c = 1 / n^i
        // double c = 1.0 / pow(n, i);
        unique_kmers += a*exp(log_b+log_c);

        // unique_kmers += a*b*c;
        fprintf(stderr, "n = %lld and i = %lld ==> a = %f, b = %lf, c = %lf\n", n, i, a, log_b, log_c);
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
    fprintf(stderr, "estimated number of unique kmers with k of %d and n of %d = %llf\n", 3, 10, expected_unique(3, 5));

    // TODO: consider header line (at least for testing purposes)

    // std::fprintf(stderr, "The header line may be more than 1\% of the file size, so the estimated k-mers may be overestimated.");
}

int main(int argc, char *argv[]) {

    test_one_read();
    return 0;


    BCL::init(); // includes MPI_Init()
    // These tests are setup well for testing accuracy, but not for time
    double start_time, end_time, setup_time, fill_time, int_time, finalize_time;
    start_time = MPI_Wtime();
   
    // Setup the sketching 
    // rm1 and rm1 of RangeMinHash
    // crmh and crmh2 of Counting RangeMinHash (from mhtest, but add back in later)
    //
    // TODO: remove all print statements from timing sections
    KWiseHasherSet<4> zomg(100);
    std::fprintf(stderr, "hv for 133: %zu\n", size_t(zomg(133, 1)));
    size_t nelem = argc == 1 ? 1000000: size_t(std::strtoull(argv[1], nullptr, 10));
    double olap_frac = argc < 3 ? 0.1: std::atof(argv[2]);
    size_t ss = argc < 4 ? 11: size_t(std::strtoull(argv[3], nullptr, 10));
    RangeMinHash<uint64_t> rm1(1 << ss), rm2(1 << ss);
    //CountingRangeMinHash<uint64_t> crmh(1 << ss), crmh2(1 << ss);
    //KthMinHash<uint64_t> kmh(30, 100);
    std::mt19937_64 mt(1337);
    size_t olap_n = (olap_frac * nelem);
    double true_ji = double(olap_n ) / (nelem * 2 - olap_n);
    olap_frac = static_cast<double>(olap_n) / nelem;

    // Insert items to both minHashes
    std::fprintf(stderr, "Inserting to both\n");
    setup_time = MPI_Wtime();
    //crmh.addh(olap_n); crmh2.addh(olap_n >> 1);
    std::set<uint64_t> z;
    
    // values inserted into both rm1 and rm2
    // TODO: get values that are kmers
    for(size_t i = 0; i < olap_n; ++i) {
        auto v = mt();
        v = WangHash()(v);
        //kmh.addh(v);
        z.insert(v);
        rm1.add(v); rm2.add(v);
    }
    std::fprintf(stderr, "olap_n: %zu. nelem: %zu\n", olap_n, nelem);

    // values inserted into rm1 only
    for(size_t i = nelem - olap_n; i--;) {
        auto v = mt();
        z.insert(v);
        rm1.addh(v);
    }

    // values inserted into rm2 only
    for(size_t i = nelem - olap_n; i--;) {
        auto v = mt();
        rm2.addh(v);
    }
    fill_time = MPI_Wtime();

    // Find the intersections
    // TODO: figure out the difference between intersection size and jaccard index (time differently?)
    
    // intersection between rm1 and rm2
    size_t is = intersection_size(rm1, rm2, typename RangeMinHash<uint64_t>::key_compare());
    double ji = rm1.jaccard_index(rm2);
    std::fprintf(stderr, "sketch is: %zu. sketch ji: %lf. True: %lf\n", is, ji, true_ji);
    assert(std::abs(ji - true_ji) / true_ji < 0.1);

    // intersection between rm1 and rm1
    is = intersection_size(rm1, rm1, typename RangeMinHash<uint64_t>::key_compare());
    ji = rm1.jaccard_index(rm1);
    assert(std::abs(1.0 - ji) < 0.01); // TODO: test sensitivity
    int_time = MPI_Wtime();
    std::fprintf(stderr, "ji for a sketch and itself: %lf\n", ji);

    /*
    // TODO: re-add when ready to test the counting hash
    mt.seed(1337);
    for(size_t i = 0; i < nelem; ++i) {
        auto v = mt();
        crmh.addh(v);
        for(size_t i = 0; i < 9; ++i)
            crmh2.addh(v);
        crmh2.addh(v * v);
    }
    std::fprintf(stderr, "is/jaccard: %zu/%lf. hist intersect: %lf.\n", size_t(crmh.intersection_size(crmh2)), crmh.jaccard_index(crmh2), crmh.histogram_intersection(crmh2));
    auto f1 = crmh.cfinalize();
    auto f2 = crmh2.cfinalize();
    std::fprintf(stderr, "finalized\n");
    std::fprintf(stderr, "crmh sum/sumsq: %zu/%zu\n", size_t(crmh.sum()), size_t(crmh.sum_sq()));
    std::fprintf(stderr, "rmh1 b: %zu. rmh1 rb: %zu. max: %zu. min: %zu\n", size_t(*rm1.begin()), size_t(*rm2.rbegin()), size_t(rm1.max_element()), size_t(rm1.min_element()));
    assert(f1.histogram_intersection(f2) == f2.histogram_intersection(f1));
    assert(f1.histogram_intersection(f2) == crmh.histogram_intersection(crmh2));
    assert(crmh.histogram_intersection(crmh2) ==  f1.tf_idf(f2) || !std::fprintf(stderr, "v1: %f. v2: %f\n", crmh.histogram_intersection(crmh2), f1.tf_idf(f2)));
    std::fprintf(stderr, "tf-idf with equal weights: %lf\n", f1.tf_idf(f2));
    std::fprintf(stderr, "f1 est cardinality: %lf\n", f1.cardinality_estimate());
    std::fprintf(stderr, "f1 est cardinality: %lf\n", f1.cardinality_estimate(ARITHMETIC_MEAN));
    std::fprintf(stderr, "f2 est cardinality: %lf\n", f1.cardinality_estimate());
    std::fprintf(stderr, "f2 est cardinality: %lf\n", f1.cardinality_estimate(ARITHMETIC_MEAN));

    crmh_time = MPI_Wtime();
    */

    // Finalizing rm1 and rm2
    auto m1 = rm1.cfinalize(), m2 = rm2.cfinalize();
    std::fprintf(stderr, "m1 is %s-sorted\n", forward_sorted(m1.begin(), m1.end()) ? "forward": "reverse");
    //auto kmf = kmh.finalize();
    std::fprintf(stderr, "jaccard between finalized MH sketches: %lf, card %lf\n", m1.jaccard_index(m2), m1.cardinality_estimate());
    {
        // TODO: figure out what this part does
        // Part 2: verify merging of final sketches
        size_t n = 10000;
        RangeMinHash<uint64_t> rmh1(20), rmh2(20);
        common::DefaultRNGType gen1(13), gen2(1337);
        for(size_t i = n; i--;) {
            rmh1.addh(gen1());
            rmh2.addh(gen2());
        }
        auto fmh1 = rmh1.cfinalize(), fmh2 = rmh2.cfinalize();
        auto u = fmh1 + fmh2;
        RangeMinHash<uint64_t> rmh3(20);
        gen1.seed(13); gen2.seed(1337);
        for(size_t i = n; i--;)
            rmh3.addh(gen1()), rmh3.addh(gen2());
        auto fmh3 = rmh3.cfinalize();
        assert(fmh3.first == u.first || (pc(fmh3, "fmh3") || pc(u, "u")));
    }
    finalize_time = MPI_Wtime();
    end_time = MPI_Wtime();

    // start_time, end_time, setup_time, fill_time, int_time, hist_time
    std::fprintf(stderr, "Setup time = %f\n", setup_time - start_time);
    std::fprintf(stderr, "Fill time time = %f\n", fill_time - setup_time);
    std::fprintf(stderr, "Intersect time = %f\n", int_time - fill_time);
    //std::fprintf(stderr, "CRMH time = %f\n", crmh_time - int_time);
    std::fprintf(stderr, "Finalizing time = %f\n", finalize_time - int_time);
    std::fprintf(stderr, "Total test time = %f\n", end_time - start_time);

    BCL::finalize();
}
