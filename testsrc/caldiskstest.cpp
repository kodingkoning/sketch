/* caldiskstest.cpp sketches two files using sketch and MPI
* Elizabeth Koning, Spring 2020
* for Senior Project at Calvin University.
*/
#include "mh.h"
#include <random>
#include <mpi.h>
#include "bcl/bcl.hpp"
#include "mpiParallelIO.cpp"

using namespace sketch;
using namespace common;
using namespace mh;

int main(int argc, char *argv[]) {
    int id;
    if (argc < 5) {
        fprintf(stderr, "\n*** Usage: caldiskstest <k> <sketchSize> <inputFile1> <inputFile2>\n\n");
        exit(1);
    }
    unsigned k = atoi(argv[1]); // k = 21 is the default for Mash. It should not go above 32 because it must be represented by an 64 bit unsigned int.
    unsigned sketchSize = atoi(argv[2]);
    std::string filename1 = argv[3];
    std::string filename2 = argv[4];
    if (k == 0 || sketchSize == 0) {
        fprintf(stderr, "\n*** k and sketchSize must be integers.\n\n");
        exit(1);
    }

    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    RangeMinHash<uint64_t> sketch(sketchSize);
    sketchFromFile(filename1, sketch, k);

    RangeMinHash<uint64_t> sketch2(sketchSize);
    sketchFromFile(filename2, sketch2, k);
    
    if (id == 0) {
        double startTime = MPI_Wtime();
    	auto s1 = sketch.cfinalize();
    	auto s2 = sketch2.cfinalize();
    	double similarity = s1.jaccard_index(s2);
        double compareTime = MPI_Wtime() - startTime;
    	std::cout << "* similarity between sketches = " << similarity << std::endl;
        std::cout << "* compare time = " << compareTime << std::endl;

        // vector<uint64_t> a = sketch2.mh2vec();
        // vector<uint64_t> b = sketch.mh2vec();

        // std::cout << "C. elegans sketch = " << std::endl;
        // for(auto ir = a.cbegin(); ir != a.cend(); ++ir) {
        //     std::cout << *ir << "\n";
        // }
        // std::cout << std::endl;
        // std::cout << "\nE. coli sketch = " << std::endl;
        // for(auto ir = b.cbegin(); ir != b.cend(); ++ir) {
        //     std::cout << *ir << "\n";
        // }
        // std::cout << std::endl;
    }
    
    MPI_Finalize();
    return 0;
}
