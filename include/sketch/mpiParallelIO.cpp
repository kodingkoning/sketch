/* arraySum.c uses an array to sum the values in an input file,
*  whose name is specified on the command-line.
* Joel Adams, Fall 2005
* for CS 374 (HPC) at Calvin College.
*
* MPI parallelism added by Elizabeth Koning, Fall 2019
* for CS 374 (HPC) at Calvin University.
*
* Modified to read chars for k-mers by Elizabeth Koning, Spring 2020
* for Senior Project at Calvin University.
*/

// TODO: test all of this

#include <stdio.h>	 /* I/O stuff */
#include <stdlib.h>	 /* calloc, etc. */
#include <mpi.h>	 /* MPI calls */
#include <string.h>	 /* strlen() */
#include <stdbool.h> /* bool */
#include "mh.h"
#include <sys/stat.h>

using namespace sketch;

void readArray(const char * fileName, char ** a, int * n);
void parallelReadArray(const char * fileName, char ** a, int * n, int id, int nProcs);
void scatterArray(char ** a, char ** allA, int * total, int * n, int id, int nProcs);
void sketchKmers(char* a, int numValues, int k, RangeMinHash<std::string> & kmerSketch);

int readFile(const char *fileName, int k, RangeMinHash<std::string>& localSketch)
{
	int nProcs, id, allCount, localCount;
	double sum;
	char *a;
	char *allA;
	double startTime, sumTime, ioTime, scatterTime, totalTime;
	bool parallelIO = true;

	MPI_Init(NULL, NULL);

	MPI_Comm_rank(MPI_COMM_WORLD, &id);
	MPI_Comm_size(MPI_COMM_WORLD, &nProcs);

	startTime = MPI_Wtime();

	if (parallelIO)
	{
		parallelReadArray(fileName, &a, &localCount, id, nProcs);
		ioTime = MPI_Wtime() - startTime;
		scatterTime = 0;
	}
	else
	{
		if (id == 0)
		{
			readArray(fileName, &allA, &allCount);
		}
		ioTime = MPI_Wtime() - startTime;
		scatterArray(&a, &allA, &allCount, &localCount, id, nProcs);
		scatterTime = MPI_Wtime() - ioTime - startTime;
	}

	// addToSketch(kmers);
	sketchKmers(a, localCount, k, localSketch);
	// instead of finding the sum of numbers, we will be adding the values to a MinHash sketch
	//   sum = parallelSumArray(a, localCount);

	sumTime = MPI_Wtime() - startTime - ioTime - scatterTime;

	totalTime = MPI_Wtime() - startTime;

	if (id == 0)
	{
		printf("The sum of the values in the input file '%s' is %g\n",
			   fileName, sum);

		printf("For %d processes:\nioTime\t\tscatterTime\tsumTime\t\ttotalTime\n%f\t%f\t%f\t%f\n\n", nProcs, ioTime, scatterTime, sumTime, totalTime);
	}

	MPI_Finalize();

	if (id == 0 && !parallelIO)
		free(allA);
	free(a);
	return 0;
}

/* parallelReadArray fills an array with values from a file.
 * Receive:	fileName, a char*
 * 		a, the address of a pointer to an array,
 * 		n, the address of an int,
 * 		id, an int id of the current process,
 * 		nProcs, an int number of MPI processes.
 * PRE: fileName contains N, followed by N double values.
 * POST: a points to a dynamically allocated array 
 * 	containing N / nProcs values from fileName and n * nProcs == N.
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

/* readArray fills an array with values from a file.
* Receive: fileName, a char*,
*          a, the address of a pointer to an array,
*          n, the address of an int.
* PRE: fileName contains N, followed by N double values.
* POST: a points to a dynamically allocated array
*        containing the N values from fileName
*        and n == N.
*/
void readArray(const char *fileName, char **a, int *n)
{
	// TODO: handle the case where we can't store all of it in memory and need to divide it up as we read
	int count, howMany;
	char *tempA;
	FILE *fin;
	struct stat sb;

	fin = fopen(fileName, "r");
	if (fin == NULL)
	{
		fprintf(stderr, "\n*** Unable to open input file '%s'\n\n",
				fileName);
		exit(1);
	}

    if(stat(fileName, &sb) == -1) {
        perror("stat");
        exit(EXIT_FAILURE);
    }
    long long size = sb.st_size; // size is in bytes, which is equal to chars
    long long bases = size / 2;
    std::fprintf(stderr, "Size of %s is %lld bytes, or about %lld bases\n", fileName, size, bases);
	howMany = size;

	tempA = (char*)calloc(howMany, sizeof(char));
	if (tempA == NULL)
	{
		fprintf(stderr, "\n*** Unable to allocate %d-length array",
				howMany);
		exit(1);
	}

	// TODO: improve this section, initially designed for doubles and can be simpler with chars
	for (count = 0; count < howMany; count++)
		fscanf(fin, "%s", &tempA[count]);

	fclose(fin);

	*n = howMany;
	*a = tempA;
}

/* scatterArray scatters the results to each process.
 * Receive:	a, the address of a pointer to an array,
 * 		allA, the address of a pointer to an array,
 * 		total, an address to an int number of elements in allA,
 * 		n, the address of an int,
 * 		id, an int id of the current process,
 * 		nProcs, an int number of MPI processes.
 * PRE: allA is filled with total number of values.
 * POST: the values in allA are scattered to a, which will contain n items.
 */
void scatterArray(char **a, char **allA, int *total, int *n, int id, int nProcs)
{
	int chunkSize, remainder, i;
	int *chunkSizes;
	int *displacements;
	char *tempA;

	MPI_Bcast(total, 1, MPI_INT, 0, MPI_COMM_WORLD);

	chunkSize = *total / nProcs;
	chunkSizes = (int*)calloc(nProcs, sizeof(int));
	displacements = (int*)calloc(nProcs, sizeof(int));
	for (i = 0; i < nProcs; ++i)
	{
		chunkSizes[i] = chunkSize;
		displacements[i] = chunkSize * i;
	}
	remainder = *total % nProcs;
	if (remainder)
	{
		chunkSizes[nProcs - 1] += remainder;
	}
	tempA = (char*)calloc(chunkSizes[id], sizeof(char));

	MPI_Scatterv(*allA, chunkSizes, displacements, MPI_CHAR, tempA, chunkSizes[id], MPI_CHAR, 0, MPI_COMM_WORLD);

	*n = chunkSizes[id];
	*a = tempA;
	free(chunkSizes);
	free(displacements);
}

/* sketchKmers adds the kmers in the data read to a Minhash sketch
 * Receive: a, a pointer to the head of an array;
 * 			numValues, the number of chars in the array.
 * Return: the MinHash sketch with the k-mers
 */
void sketchKmers(char* a, int numValues, int k, RangeMinHash<std::string> & kmerSketch) {
	// RangeMinHash<std::string> kmerSketch();
	// TODO: analyze values
	char * kmer = (char*)calloc(k, sizeof(char));
	// sketchedKmers = *kmerSketch;
	// return kmerSketch;
	int kmer_len = 0;
	for(int i = 0; i < numValues; i++) {
		if(a[i] == 'A' || a[i] == 'T' || a[i] == 'C' || a[i]== 'G') {
			if(kmer_len < k) {
				kmer_len++;
				kmer[kmer_len] = a[i];
			} else {
				// TODO: check against threshold for the hash values (will need to send the hash value to the sketch for confirmation)
				kmerSketch.add(kmer);
				for(int j = 0; j < k-1; j++) {
					// TODO: make this more efficient. Maybe keep an index stored of the first...
					kmer[j] = kmer[j+1];
				}
				kmer[k-1] = a[i];
			}
		} else {
			kmer_len = 0;
		}
	}
}

/* sumArray sums the values in an array of doubles.
* Receive: a, a pointer to the head of an array;
*          numValues, the number of values in the array.
* Return: the sum of the values in the array.
*/

double sumArray(double *a, int numValues)
{
	int i;
	double result = 0.0;

	for (i = 0; i < numValues; i++)
	{
		result += *a;
		a++;
	}

	return result;
}

/* parallelSumArray sums the values in an array of doubles.
* Receive: a, a pointer to the head of an array;
*          numValues, the number of values in the array.
* Return: the sum of the values in the array.
*/

double parallelSumArray(double *a, int numValues)
{
	int i;
	double result = 0.0;
	double resultSum = 0.0;

	for (i = 0; i < numValues; i++)
	{
		result += *a;
		a++;
	}

	MPI_Reduce(&result, &resultSum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	return resultSum;
}
