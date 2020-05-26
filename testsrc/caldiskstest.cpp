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
    if (argc < 3) {
        fprintf(stderr, "\n*** Usage: caldiskstest <inputFile1> <inputFile2>\n\n");
        exit(1);
    }
    std::string fastq_file1 = argv[1];
    std::string fastq_file2 = argv[2];
    std::string dir = get_current_dir_name();
    std::string filename1 = dir+"/" + fastq_file1;
    std::string filename2 = dir + "/" + fastq_file2;

    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    RangeMinHash<uint64_t> sketch(LOCAL_SKETCH_SIZE);
    sketchFromFile(filename1, sketch);

    RangeMinHash<uint64_t> sketch2(LOCAL_SKETCH_SIZE);
    sketchFromFile(filename2, sketch2);
    
    vector<uint64_t> a = sketch2.mh2vec();
    vector<uint64_t> b = sketch.mh2vec();

    if (id == 0) {
        double startTime = MPI_Wtime();
    	auto s1 = sketch.cfinalize();
    	auto s2 = sketch2.cfinalize();
    	double similarity = s1.jaccard_index(s2);
        double compareTime = MPI_Wtime() - startTime;
    	std::cout << "* similarity between sketches = " << similarity << std::endl;
        std::cout << "* compare time = " << compareTime << std::endl;
    }
    
    MPI_Finalize();
    return 0;
}
