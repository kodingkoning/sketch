/* mpiParallelIO.cpp handles the parallel I/O for Cal-DisKS
* Elizabeth Koning, Spring 2020
* for Senior Project at Calvin University.
*/

#include <stdio.h>	 /* I/O stuff */
#include <stdlib.h>	 /* calloc, etc. */
#include <mpi.h>	 /* MPI calls */
#include <string.h>	 /* strlen() */
#include <stdbool.h> /* bool */
#include <sys/stat.h>
#include <iostream>
#include "mh.h"
#include "calcThreshold.cpp"

using namespace sketch;

int readFile(const char *fileName, int k, RangeMinHash<uint64_t>& localSketch, int nProcs, int id);
void readArray(const char * fileName, char ** a, int * n);
void parallelReadArray(const char * fileName, char ** a, int * n, int id, int nProcs);
void scatterArray(char ** a, char ** allA, int * total, int * n, int id, int nProcs);
void sketchKmers(char* a, int numValues, int k, RangeMinHash<uint64_t> & kmerSketch);
void combineSketches(RangeMinHash<uint64_t> & localSketch, RangeMinHash<uint64_t> & globalSketch, int nProcs, int id);

void sketchFromFile(std::string filename, RangeMinHash<uint64_t>& globalSketch) {
    int k = 21; // k = 21 is the default for Mash. It should not go above 32 because it must be represented by an 64 bit unsigned int.
	int nProcs, id;
    double startTime, totalTime, threshTime, sketchTime, gatherTime;

	MPI_Comm_rank(MPI_COMM_WORLD, &id);
	MPI_Comm_size(MPI_COMM_WORLD, &nProcs);

    startTime = MPI_Wtime();

	// TODO : would be better to only calculate this on one process
    // BigInt threshold = find_threshold(filename, k, SKETCH_SIZE);
    threshTime = MPI_Wtime();

    RangeMinHash<uint64_t> localSketch(LOCAL_SKETCH_SIZE);

	std::cout << "Created sketch objects..." << std::endl;

    readFile(filename.c_str(), k, localSketch, nProcs, id);
	std::cout << "Read FASTQ file and made local sketches..." << std::endl;

    sketchTime = MPI_Wtime();

    combineSketches(localSketch, globalSketch, nProcs, id);
	std::cout << "Combined the sketches..." << std::endl;

    gatherTime = MPI_Wtime();

    totalTime = MPI_Wtime() - startTime;
    std::cout << "\n* Threshold calculation time = \t" << (threshTime - startTime) << std::endl;
    std::cout << "\n* Local sketching time = \t" << (sketchTime - threshTime) << std::endl;
    std::cout << "\n* Sketch gather time = \t\t" << (gatherTime - sketchTime) << std::endl;
    std::cout << "\n* Total Cal_DisKS time = \t" << (totalTime - startTime) << std::endl;
}

int readFile(const char *fileName, int k, RangeMinHash<uint64_t>& localSketch, int nProcs, int id)
{
	int allCount, localCount;
	char *a;
	double startTime, sketchTime, ioTime, totalTime;

	startTime = MPI_Wtime();

	parallelReadArray(fileName, &a, &localCount, id, nProcs);
	ioTime = MPI_Wtime() - startTime;

	sketchKmers(a, localCount, k, localSketch);

	sketchTime = MPI_Wtime() - startTime - ioTime;

	totalTime = MPI_Wtime() - startTime;

	if (id == 0)
	{
		printf("For %d processes:\nioTime\t\tsketchTime\t\ttotalTime\n%f\t%f\t%f\t%f\n\n", nProcs, ioTime, sketchTime, totalTime);
	}

	free(a);
	return 0;
}

/* parallelReadArray fills an array with values from a file.
 * Receive:	fileName, a char*
 * 		a, the address of a pointer to an array,
 * 		n, the address of an int,
 * 		id, an int id of the current process,
 * 		nProcs, an int number of MPI processes.
 * PRE: fileName contains k-mers, and may contain other characters.
 * POST: a points to a dynamically allocated array 
 * 	containing file size / nProcs values from fileName.
 */
void parallelReadArray(const char *fileName, char **a, int *n, int id, int nProcs)
{
	int count, howMany, offset, chunkSize, remainder, headerLen, numLen, chunkChars;
	const int DEFAULT_BUF_LEN = 10;
	int error;
	MPI_File file;
	MPI_Status status;
	MPI_Offset fileSize;
	char *buffer;

	// open MPI file for parallel I/O
	error = MPI_File_open(MPI_COMM_WORLD, fileName, MPI_MODE_RDONLY, MPI_INFO_NULL, &file);
	if (error != MPI_SUCCESS)
	{
		fprintf(stderr, "\n*** Unable to open input file '%s'\n\n", fileName);
	}

	// get the size of the file
	error = MPI_File_get_size(file, &fileSize);
	if (error != MPI_SUCCESS)
	{
		fprintf(stderr, "\n*** Unable to get size of file '%s'\n\n", fileName);
	}

	// find size of each process's chunk
	chunkSize = fileSize / nProcs;
	offset = id * chunkSize; // TODO: add room on either end to handle not losing k-mers between reads
	remainder = chunkSize % nProcs;
	if (remainder && id == nProcs - 1)
	{
		chunkSize += remainder;
	}

	buffer = (char *)calloc(chunkSize + 1, sizeof(char));
	if (buffer == NULL)
	{
		fprintf(stderr, "\n** Unable to allocate %d-length array", chunkSize);
	}
	MPI_File_read_at(file, offset, buffer, chunkSize, MPI_CHAR, &status);
	// buffer will contain all of the chars
	// TODO: add buffer room at the ends so that we don't lose some of the k-mers

	MPI_File_close(&file);

	*n = chunkSize;
	*a = buffer;
}

/* complemenntbase() and reversecomplement() come from BELLA code
 */
char
complementbase(char n) {
	switch(n)
	{
	case 'A':
		return 'T';
	case 'T':
		return 'A';
	case 'G':
		return 'C';
	case 'C':
		return 'G';
	}
	assert(false);
	return ' ';
}

std::string
reversecomplement(const std::string& seq) {

	std::string cpyseq = seq;
	std::reverse(cpyseq.begin(), cpyseq.end());

	std::transform(
		std::begin(cpyseq),
		std::end  (cpyseq),
		std::begin(cpyseq),
	complementbase);

	return cpyseq;
}

// modified from 32 bit version at from https://github.com/Ensembl/treebest/blob/master/common/hash_com.h
inline uint64_t kmer_int(const char *s) {
	uint64_t h = 0;
	for ( ; *s; s++)
		h = (h << 5) - h + *s;
	return h;
} 


/* sketchKmers adds the kmers in the data read to a Minhash sketch
 * Receive: a, a pointer to the head of an array;
 * 			numValues, the number of chars in the array;
 * 			k, the number of bases in a k-mer;
 * 			kmerSketch, the empty sketch to fill with k-mers;
 * Postcondition: kmerSketch is filled with k-mers from a.
 */
void sketchKmers(char* a, int numValues, int k, RangeMinHash<uint64_t> & kmerSketch) {
	std::string kmer = "";
	for(int i = 0; i < numValues; i++) {
		if(a[i] == 'A' || a[i] == 'T' || a[i] == 'C' || a[i]== 'G') {
			if(kmer.length() < k) {
				kmer.push_back(a[i]);
			} else {
				// TODO: check against threshold for the hash values (will need to send the hash value to the sketch for confirmation)
				std::string twin = reversecomplement(kmer);
				if (twin < kmer) {
					kmerSketch.add(kmer_int(twin.c_str()));
				} else {
					kmerSketch.add(kmer_int(kmer.c_str()));
				}
				kmer = kmer.substr(1, k-1) + a[i];
			}
		} else {
			kmer = "";
		}
	}
}

/* combineSketches adds the kmers in the data read to a Minhash sketch
 * Receive: localSketch, the local sketch for easy MPI process;
 * 			globalSketch, the global sketch for process 0 to gather the hash values;
 * 			nProcs, the number of MPI processes;
 * 			id, the id of the current MPI process;
 * Postcondition: globalSketch for process 0 has the minimum values from the local sketches.
 */
void combineSketches(RangeMinHash<uint64_t> & localSketch, RangeMinHash<uint64_t> & globalSketch, int nProcs, int id) {
	unsigned num_vals = localSketch.size();
	uint64_t * local_data = localSketch.mh2vec().data();
	uint64_t * global_data = NULL;
	if(id == 0) {
		global_data = new uint64_t[num_vals*nProcs];
	}

	MPI_Gather(local_data, num_vals,  MPI_UNSIGNED_LONG_LONG, global_data, num_vals, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);

	if(id == 0) {
		for (unsigned i = 0; i < num_vals*nProcs; i++) {
			globalSketch.add(global_data[i]);
		}
	}

	delete [] global_data;
}
