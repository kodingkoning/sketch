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

bool debug = true;

void readArray(const char * fileName, char ** a, int * n);
int parallelReadArray(const char * fileName, char ** a, int * n, int id, int nProcs, unsigned k);
void scatterArray(char ** a, char ** allA, int * total, int * n, int id, int nProcs);
void sketchKmers(char* a, int numValues, unsigned k, RangeMinHash<uint64_t> & kmerSketch, int id);
void combineSketches(RangeMinHash<uint64_t> & localSketch, RangeMinHash<uint64_t> & globalSketch, int nProcs, int id);
void sketchReduction(RangeMinHash<uint64_t> & localSketch, RangeMinHash<uint64_t> & globalSketch, int id, int nProcs);

void sketchFromFile(std::string filename, RangeMinHash<uint64_t>& globalSketch, unsigned k) {
	int nProcs, id, error, minChunks, chunksPerProc;
	double startTime, totalTime, threshTime, ioTime, sketchTime, gatherTime, tempTime;
	int localCount;
	int readChunks = 1;
	char *a;
	MPI_File file;
	MPI_Offset fileSize;

	MPI_Comm_rank(MPI_COMM_WORLD, &id);
	MPI_Comm_size(MPI_COMM_WORLD, &nProcs);

    startTime = MPI_Wtime();

	// TODO : would be better to only calculate this on one process
    // BigInt threshold = find_threshold(filename, k, SKETCH_SIZE);
    threshTime = MPI_Wtime();

    RangeMinHash<uint64_t> localSketch(globalSketch.sketch_size());
    
    // open MPI file for parallel I/O
    error = MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &file);
    if (error != MPI_SUCCESS) {
	    fprintf(stderr, "\n*** Unable to open input file '%s'\n\n", filename);
    }

    // get the size of the file
    error = MPI_File_get_size(file, &fileSize);
    if (error != MPI_SUCCESS) {
	    fprintf(stderr, "\n*** Unable to get size of file '%s'\n\n", filename);
    }

    minChunks = fileSize / INT_MAX + (fileSize % INT_MAX != 0);
    if(debug) std::cout << "Minimum chunks = " << minChunks << std::endl;
    if(minChunks <= nProcs) {
	    chunksPerProc = 1;
	    parallelReadArray(filename.c_str(), &a, &localCount, id, nProcs, k);
	    ioTime = MPI_Wtime();
	    sketchKmers(a, localCount, k, localSketch, id);
	    free(a);
	    sketchTime = MPI_Wtime();
    } else {
	    ioTime = 0;
	    sketchTime = 0;
	    chunksPerProc = minChunks / nProcs + (minChunks % INT_MAX != 0);
	    for(int chunk = id*chunksPerProc; chunk < (id+1)*chunksPerProc; ++chunk) {
		    tempTime = MPI_Wtime();
		    parallelReadArray(filename.c_str(), &a, &localCount, chunk, nProcs*chunksPerProc, k);
		    ioTime += MPI_Wtime() - tempTime;
		    sketchKmers(a, localCount, k, localSketch, id);
		    free(a);
		    sketchTime += MPI_Wtime() - tempTime;
	    }
    }
	/*
	int readStatus = parallelReadArray(filename.c_str(), &a, &localCount, id, nProcs, k);
	if(debug) std::cout << "parallelReadArray() done" << std::endl;
	if(readStatus) {
		ioTime = 0;
		sketchTime = 0;
		if(debug) std::cout << "Process " << id << ": splitting into chunks." << std::endl;
		int newChunks = 2;
		for(int chunkIndex = 0; chunkIndex < newChunks; ++chunkIndex) {
			tempTime = MPI_Wtime();
			readStatus = parallelReadArray(filename.c_str(), &a, &localCount, id+newChunks, newChunks*nProcs, k);
			if(readStatus) {
				fprintf(stderr, "\n*** Unable to allocate memory to read the array.\n\n");
			return;
			}
			ioTime += MPI_Wtime() - tempTime;
			sketchKmers(a, localCount, k, localSketch, id);
			free(a);
			sketchTime += MPI_Wtime() - tempTime;
		}
	} else {
		ioTime = MPI_Wtime();
		sketchKmers(a, localCount, k, localSketch, id);
		free(a);
		sketchTime = MPI_Wtime();
	}
	*/ /*
	//TODO: forcibly test if this actually works, at least if done once
	ioTime = 0;
	sketchTime = 0;
	while(readStatus) {
		readChunks *= 2;
		if(debug) std::cout << "splitting read chunks to " << readChunks << std::endl;
		for(int chunkIndex = 0; chunkIndex < readChunks; chunkIndex++) {
			tempTime = MPI_Wtime();
			readStatus = parallelReadArray(filename.c_str(), &a, &localCount, id*readChunks + chunkIndex, nProcs*readChunks, k);
			ioTime += MPI_Wtime() - tempTime;
			if(readStatus) {
				chunkIndex = readChunks;
			}
			else {
				tempTime = MPI_Wtime();
				sketchKmers(a, localCount, k, localSketch, id);
				free(a);
				if(debug) std::cout << "Process " << id << " sketched chunk " << id*readChunks+chunkIndex << " of " << nProcs*readChunks << std::endl;
				sketchTime += MPI_Wtime() - tempTime;
			}
		}
		ioTime = threshTime + ioTime;
		sketchTime = ioTime + sketchTime;
	}
	if(readChunks == 1) {
		ioTime = MPI_Wtime();
		sketchKmers(a, localCount, k, localSketch, id);
		free(a);
    	sketchTime = MPI_Wtime();
	}
	*/

	/*
	ioTime = MPI_Wtime();
	sketchKmers(a, localCount, k, localSketch, id);
	free(a);
	sketchTime = MPI_Wtime();
	*/

	if(debug) std::cout << "Process " << id << ": Local sketching complete." << std::endl;

    combineSketches(localSketch, globalSketch, nProcs, id); // TODO: pick which approach is what we want
	//sketchReduction(localSketch, globalSketch, id, nProcs);
	if(debug) std::cout << "Process " << id << ": Sketches combined." << std::endl;

    gatherTime = MPI_Wtime();

    totalTime = MPI_Wtime() - startTime;

    if (id == 0) {
		std::cout << "For file " << filename << " with " << nProcs << " processes: " << std::endl;
    	std::cout << " * Used " << chunksPerProc << "  chunks per process" << std::endl;
    	std::cout << " * Threshold calculation time = " << (threshTime - startTime) << std::endl;
    	std::cout << " * Parallel read from file time = " << (ioTime - threshTime) << std::endl;
    	std::cout << " * Local sketching time = \t" << (sketchTime - ioTime) << std::endl;
    	std::cout << " * Sketch combine time = \t" << (gatherTime - sketchTime) << std::endl;
    	std::cout << " * Total Cal_DisKS time = \t" << (totalTime) << std::endl;
    }
}

/* parallelReadArray fills an array with values from a file.
 * Receive:	fileName, a char*
 * 		a, the address of a pointer to an array,
 * 		n, the address of an int,
 * 		id, an int id of the current process,
 * 		nProcs, an int number of MPI processes
 * 		k, an int for the length of the k-mers.
 * PRE: fileName contains k-mers, and may contain other characters.
 * POST: a points to a dynamically allocated array 
 * 	containing file size / nProcs values from fileName.
 */
int parallelReadArray(const char *fileName, char **a, int *n, int id, int nProcs, unsigned k)
{
	int error;
	MPI_File file;
	MPI_Status status;
	int chunkSize;
	MPI_Offset fileSize, offset, remainder;
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
	offset = id * chunkSize;
	remainder = fileSize % nProcs;
	if (remainder && id == nProcs - 1)
	{
		chunkSize += remainder;
	}
	if (id < nProcs -1) { // adds room on end to account for k-mers in part of each I/O section
		chunkSize += k + 1;
	}

	if(chunkSize < 0) {
		fprintf(stderr, "Process %3d: chunkSize < 0\n", id);
		return 1;
	}

	buffer = (char *)calloc(chunkSize + 1, sizeof(char));
	if (buffer == NULL)
	{
		fprintf(stderr, "\n*** Unable to allocate memory to read \n\n");
		return 1;
	}
	if(debug) std::cout << "Process " << id << ": file = " << file << ", offset = " << offset << ", chunkSize = " << chunkSize << std::endl;
	error = MPI_File_read_at(file, offset, buffer, chunkSize, MPI_CHAR, &status);
	if (error != MPI_SUCCESS) {
		fprintf(stderr, "\n*** Unable to read from input file\n\n");
		char error_string[BUFSIZ];
                int length_of_error_string;
                MPI_Error_string(error, error_string, &length_of_error_string);
                fprintf(stderr, "%3d: %s\n", id, error_string);
	} else {
	}

	MPI_File_close(&file);

	*n = chunkSize;
	*a = buffer;
	return 0;
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
void sketchKmers(char* a, int numValues, unsigned k, RangeMinHash<uint64_t> & kmerSketch, int id) {
	std::string kmer = "";
	for(int i = 0; i < numValues; ++i) {
		if(a[i] == 'A' || a[i] == 'T' || a[i] == 'C' || a[i]== 'G') {
			if(kmer.length() < k) {
				kmer.push_back(a[i]);
			}
			if(kmer.length() == k) {
				// TODO: check against threshold for the hash values (will need to send the hash value to the sketch for confirmation)
				std::string twin = reversecomplement(kmer);
				if (twin < kmer) {
					kmerSketch.addh(kmer_int(twin.c_str()));
				} else {
					kmerSketch.addh(kmer_int(kmer.c_str()));
				}
				kmer = kmer.substr(1, k-1) + a[i]; // start at 1 and get k-1 chars
			}
		} else {
			kmer = "";
		}
	}
	if(debug) {
		if (kmerSketch.size() == 0) {
			int counter = 0;
			std::cout << "Process " << id << " has 0 kmers" << std::endl;
		for(int i = 0; i < numValues; ++i) if(!(a[i] == 'A' || a[i] == 'T' || a[i] == 'C' || a[i]== 'G')) ++counter;
			std::cout << "Process " << id << " has " << counter << " non-bases " <<  std::endl;
		}
		std::cout << "Process " << id << ": first char = " << (a[0] == 0) << std::endl;
	}
	if(debug) std::cout << "Process " << id << ": size = " << kmerSketch.size() << std::endl;
}

/* combineSketches adds the kmers in the data read to a Minhash sketch
 * Receive: localSketch, the local sketch for easy MPI process;
 * 			globalSketch, the global sketch for process 0 to gather the hash values;
 * 			nProcs, the number of MPI processes;
 * 			id, the id of the current MPI process;
 * Postcondition: globalSketch for process 0 has the minimum values from the local sketches.
 */
void combineSketches(RangeMinHash<uint64_t> & localSketch, RangeMinHash<uint64_t> & globalSketch, int nProcs, int id) {

	if(nProcs == 1) {
		// TODO: decide if this simplification should stay
		globalSketch += localSketch;
		return;
	}

	unsigned num_vals = localSketch.sketch_size();
	vector<uint64_t> local_data = localSketch.mh2vec();
	uint64_t * global_data = NULL;
	int error_code;

	if(id == 0) {
		global_data = (uint64_t *)calloc(num_vals*nProcs, sizeof(uint64_t));
		if(global_data == NULL) {
			std::cout << "Unable to allocate array for global data, so cannot gather" << std::endl;
			return;
		}
	}

	if(debug) {
		std::cout << "Process " << id << ": " << localSketch.size() << std::endl;
		if(localSketch.min_element() == 0) {
			std::cout << "Process " << id << ": minimum element is 0"<< std::endl;
		} else {
			std::cout << "Process " << id << ": minimum element is not 0!"<< std::endl;
		}
	}

	if(debug) std::cout << "Process " << id << ": size = " << localSketch.size() << " and capacity = " << localSketch.sketch_size() << std::endl;

	error_code = MPI_Gather(local_data.data(), num_vals,  MPI_UNSIGNED_LONG_LONG, global_data, num_vals, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);

	if(error_code != MPI_SUCCESS) {
		char error_string[BUFSIZ];
		int length_of_error_string;
		MPI_Error_string(error_code, error_string, &length_of_error_string);
		fprintf(stderr, "%3d: %s\n", id, error_string);
	}

	if(id == 0) {
		for (unsigned i = 0; i < num_vals*nProcs; i++) {
			globalSketch.add(global_data[i]);
		}
	}

	free(global_data);
}

void sketchReduction(RangeMinHash<uint64_t> & localSketch, RangeMinHash<uint64_t> & globalSketch, int id, int nProcs) {
	int n;
	int n_vals = localSketch.sketch_size();
	vector<uint64_t> local_data = localSketch.mh2vec();
	uint64_t * buffer = NULL; // (uint64_t *)calloc(n_vals, sizeof(uint64_t));
	if (id%2 == 0) {
		buffer = (uint64_t *)calloc(n_vals, sizeof(uint64_t));
		if(buffer == NULL) { std::cout << "Process " << id << " unable to receive sketches." << std::endl; return; }
	}
	for(n = 1; n < nProcs; n *= 2) {
		// TODO: id receives, id+n sends data
		if( id%(n*2) == 0 && id+n < nProcs) {
			// TODO: receive data from id + n
			if(debug) std::cout << "receiving to " << id << " from " << id+n << std::endl;
			MPI_Recv(buffer, n_vals, MPI_UNSIGNED_LONG_LONG, id+n, 0, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			for(int i = 0; i < n_vals; ++i) {
				localSketch.add(buffer[i]);
			}
		} else if((id-n)%(n*2) == 0) {
			// TODO: send data to id - n
			if(debug) std::cout << "sending from " << id << " to " << id-n << std::endl;
			MPI_Send(local_data.data(), local_data.size(), MPI_UNSIGNED_LONG_LONG, id-n, 0, MPI_COMM_WORLD);
		}
	}
	free(buffer);

	if(id == 0) {
		globalSketch += localSketch;
	}
}
