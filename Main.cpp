// PR_projekt1.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <stdio.h>
#include <time.h>
#include <omp.h>
#include <cmath>

#include <iostream>
#include <fstream>
#include <ctime>

// MIN >= 2
const int MIN = 2;
const int MAX = 250000000;

// W podejsciu domenowym watki dziela prace na THREAD_NUM * DOMAIN_MULTIPLIER przedzialow
const int DOMAIN_MULTIPLIER = 2;

#define THREAD_NUM 4


const bool eratostenes = false;
const bool printPrimes = false;

enum version { SEQ_DIV, SEQ_DIV_PRIM, SEQ_ERA, SEQ_ERA_PRIM, PAR_DIV, PAR_DIV_PRIM, PAR_ERA_FUN, PAR_ERA_FUN_PRIM, PAR_ERA_DOM, PAR_ERA_DOM_PRIM, PAR_DIV_2 };

// prototypy funkcji pomocniczych
void showResult(bool* matrix);
// prototypy funkcji do przypisania dla tablicy wskaznikow do funkcji
void sequentialDivision(bool* matrix);
void sequentialDivisionWithPrimes(bool* matrix);
void sequentialEratostenes(bool* matrix);
void sequentialEratostenesWithPrimes(bool* matrix);
void parallelDivision(bool* matrix);
void parallelDivisionWithPrimes(bool* matrix);
void parallelEratostenesFunctional(bool* matrix);
void parallelEratostenesFunctionalWithPrimes(bool* matrix);
void parallelEratostenesDomain(bool* matrix);
void parallelEratostenesDomainWithPrimes(bool* matrix);
void parallelDivisionLivePrimes(bool* matrix);

void (*funcPtr[])(bool*) = { sequentialDivision, sequentialDivisionWithPrimes, sequentialEratostenes, sequentialEratostenesWithPrimes, parallelDivision, parallelDivisionWithPrimes, parallelEratostenesFunctional,
							  parallelEratostenesFunctionalWithPrimes, parallelEratostenesDomain, parallelEratostenesDomainWithPrimes, parallelDivisionLivePrimes };


//glowne funkcje
void cleanMatrix(bool* matrix, bool val, int maxRange) {
	for (int i = 0; i <= maxRange; i++) {
		matrix[i] = val;
	}
}

bool isPrime(int n) {
	if (n == 2)
	{
		return true;
	}
	else if (n % 2 == 0) {
		return false;
	}
	int squareRoot = sqrt(n);
	for (int i = 3; i <= squareRoot; i += 2) {
		if (n % i == 0) {
			return false;
		}
	}
	return true;
}

bool isPrimeExtended(int n, int* primes, int prime_counter) {
	if (n == 2)
	{
		return true;
	}
	else if (n % 2 == 0) {
		return false;
	}
	int squareRoot = sqrt(n);

	for (int i = 0; i < prime_counter; i++) {
		if (primes[i] > squareRoot)
		{
			break;
		}
		if (n % primes[i] == 0) {
			return false;
		}
	}

	if (prime_counter > 1)
	{
		for (int i = primes[prime_counter - 1] + 2; i <= squareRoot; i += 2) {
			if (n % i == 0) {
				return false;
			}
		}
	}
	else
	{
		for (int i = 3; i <= squareRoot; i += 2) {
			if (n % i == 0) {
				return false;
			}
		}
	}
	return true;
}



// Wersje sekwencyjne
void sequentialDivision(bool* matrix) {
	for (int i = MIN; i <= MAX; i++) {
		if (isPrime(i)) {
			matrix[i - MIN] = true;
		}
		else
		{
			matrix[i - MIN] = false;
		}
	}
}

void sequentialDivisionWithPrimes(bool* matrix) {

	int squareRoot = sqrt(MAX);
	int* primes = new int[squareRoot];
	int prime_counter = 0;

	for (int i = 2; i <= squareRoot; i++)
	{
		if (isPrimeExtended(i, primes, prime_counter)) {
				primes[prime_counter] = i;
				prime_counter++;

				if (i >= MIN)
				{
					matrix[i - MIN] = true;
				}
			
		}
		else if (i >= MIN) {
			matrix[i - MIN] = false;
		}
	}
	for (int i = (MIN > squareRoot + 1) ? MIN : squareRoot + 1; i <= MAX; i++) {
		if (isPrimeExtended(i, primes, prime_counter)) {
			matrix[i - MIN] = true;
			if (prime_counter < squareRoot) {
				primes[prime_counter] = i;
				prime_counter++;
			}
		}
		else
		{
			matrix[i - MIN] = false;
		}


	}

	delete[] primes;

}


void sequentialEratostenes(bool* matrix) {
	cleanMatrix(matrix, true, MAX - 2);
	int squareRoot = sqrt(MAX);
	for (int counter = 2; counter <= squareRoot; counter++)
	{
		if (matrix[counter - 2] == false) {
			continue;
		}
		for (int i = counter * 2; i <= MAX; i += counter) {
			matrix[i - 2] = false;
		}
	}
}

void sequentialEratostenesWithPrimes(bool* matrix) {
	cleanMatrix(matrix, true, MAX - 2);
	int squareRoot = sqrt(MAX);
	int* primes = new int[squareRoot];
	int primeCounter = 0;
	int secondSquareRoot = sqrt(squareRoot);
	for (int counter = 2; counter <= secondSquareRoot; counter++)
	{
		if (matrix[counter - 2] == false) {
			continue;
		}
		for (int i = counter * 2; i <= squareRoot; i += counter) {
			matrix[i - 2] = false;
		}
	}
	for (int i = 0; i <= squareRoot; i++) {
		if (matrix[i] == true) {
			primes[primeCounter] = i + 2;
			primeCounter++;
		}
	}
	for (int i = 0; i < primeCounter; i++) {
		for (int j = primes[i] * 2; j <= MAX; j += primes[i]) {
			matrix[j - 2] = false;
		}
	}
	delete[] primes;
}

// Wersje równoległe
void parallelDivision(bool* matrix) {
#pragma omp parallel
	{
#pragma omp for
		for (int i = MIN; i <= MAX; i++) {
			if (isPrime(i)) {
				matrix[i - MIN] = true;
			}
			else
			{
				matrix[i-MIN] = true;
			}
		}
	}
}


void parallelDivisionLivePrimes(bool* matrix) {

	cleanMatrix(matrix, false, MAX - MIN);
	int squareRoot = sqrt(MAX);

	int* primes = new int[squareRoot];
	int prime_counter = 0;

	int last_j = 0;

#pragma omp parallel
	{
		int threadId = omp_get_thread_num();
#pragma omp for schedule(static)
		for (int i = MIN; i <= MAX; i++) {

			if (isPrimeExtended(i, primes, prime_counter)) {
				matrix[i - MIN] = true;

				if (threadId == 0)
				{

					if (prime_counter < squareRoot)
					{

						for (int j = last_j + 1; j < i; j++)
						{
							if (matrix[j - MIN] == true) {

								primes[prime_counter] = j;
								prime_counter++;

								if (prime_counter >= squareRoot)
								{
									break;
								}
							}
						}
						last_j = i;
					}

				}

			}


		}
	}

	delete[] primes;
}

void parallelDivisionWithPrimes(bool* matrix) {

	//cleanMatrix(matrix, false, MAX - MIN);
	int squareRoot = sqrt(MAX);

	int* primes = new int[int(squareRoot / 2) + 10];
	int prime_counter = 0;

	bool* matrixUnderMin = new bool[MIN];
	//cleanMatrix(matrixUnderMin, false, MIN - 1);


#pragma omp parallel
	{

#pragma omp for 
		for (int i = 2; i <= squareRoot; i++) {

			if (isPrime(i)) {

				if (i < MIN) {
					matrixUnderMin[i] = true;
				}
				else {
					matrix[i - MIN] = true;
				}
			}
			else
			{
				if (i < MIN) {
					matrixUnderMin[i] = false;
				}
				else {
					matrix[i - MIN] = false;
				}

			}
		}

#pragma omp single
		{
			for (int i = 2; i <= squareRoot; i++)
			{
				if (i < MIN)
				{
					if (matrixUnderMin[i] == true)
					{
						primes[prime_counter] = i;
						prime_counter += 1;
					}
				}
				else
				{
					if (matrix[i - MIN] == true)
					{
						primes[prime_counter] = i;
						prime_counter += 1;
					}
				}
			}
		}

#pragma omp for
		for (int i = (squareRoot + 1 > MIN) ? squareRoot + 1 : MIN; i <= MAX; i++) {
			if (isPrimeExtended(i, primes, prime_counter)) {
				matrix[i - MIN] = true;
			}
			else {
				matrix[i - MIN] = false;
			}
		}

	}

	delete[] primes;
	delete[] matrixUnderMin;
}

/*void parallelEratostenesFunctional(bool* matrix) {
	cleanMatrix(matrix, true);
	if (MAX > 0) {
		matrix[1] = false;
	}
	matrix[0] = false;
	int counter = 2;
	int squareRoot = sqrt(MAX);
#pragma omp parallel
	{
		int threadId = omp_get_thread_num();
		int numThreads = omp_get_num_threads();
		for (int i = counter + threadId; i <= squareRoot; i+=numThreads)
		{
			printf("Thread: %d, iter: %d\n", threadId, i);
			if (matrix[i] == false) {
				continue;
			}
			for (int j = i * 2; j <= MAX; j += i) {
				matrix[j] = false;
			}
		}
	}
}*/

void parallelEratostenesFunctional(bool* matrix) {
	cleanMatrix(matrix, true, MAX - 2);
	int counter = 2;
	int squareRoot = sqrt(MAX);
#pragma omp parallel
	{
		//int threadId = omp_get_thread_num();
#pragma omp for
		for (int i = counter; i <= squareRoot; i++)
		{
			//printf("Thread: %d, iter: %d\n", threadId, i);
			if (matrix[i - 2] == true) {
				for (int j = i * 2; j <= MAX; j += i) {
					matrix[j - 2] = false;
				}
			}
		}
	}
}

void parallelEratostenesFunctionalWithPrimes(bool* matrix) {
	cleanMatrix(matrix, true, MAX - 2);
	int counter = 2;
	int squareRoot = sqrt(MAX);
	int* primes = new int[squareRoot];
	int primeCounter = 0;
	int secondSquareRoot = sqrt(squareRoot);
#pragma omp parallel
	{
		//int threadId = omp_get_thread_num();
#pragma omp for
		for (int i = counter; i <= secondSquareRoot; i++)
		{
			//printf("Thread: %d, iter: %d\n", threadId, i);
			if (matrix[i - 2] == true) {
				for (int j = i * 2; j <= squareRoot; j += i) {
					matrix[j - 2] = false;
				}
			}
		}
	}
	for (int i = 0; i <= squareRoot; i++) {
		if (matrix[i] == true) {
			primes[primeCounter] = i + 2;
			primeCounter++;
		}
	}
#pragma omp for
	for (int i = 0; i < primeCounter; i++) {
		for (int j = primes[i] * 2; j <= MAX; j += primes[i]) {
			matrix[j - 2] = false;
		}
	}
	delete[] primes;
}
// STARA WERSJA, DZIALA, ALE NIE KORZYSTA Z PRAGMA OMP FOR
/*void parallelEratostenesDomain(bool* matrix) {
	cleanMatrix(matrix, true);
	int squareRoot = sqrt(MAX);
	int numThreads = omp_get_num_threads();
#pragma omp parallel
	{
		int threadId = omp_get_thread_num();
		int lowerBound = threadId * (squareRoot / numThreads);
		if (lowerBound < 2) lowerBound = 2;
		int upperBound = (threadId + 1) * (squareRoot / numThreads);
		for (int i = lowerBound; i <= upperBound; i++)
		{
			if (matrix[i - 2] == false) {
				continue;
			}
			for (int j = i * 2; j <= MAX; j += i) {
				matrix[j - 2] = false;
			}
		}
	}
}*/

void parallelEratostenesDomain(bool* matrix) {
	cleanMatrix(matrix, true, MAX - 2);
	int squareRoot = sqrt(MAX);
	int numThreads = THREAD_NUM;
	printf("Num threads: %d\n", numThreads);
	int range = MAX / (numThreads * DOMAIN_MULTIPLIER);
	printf("Range: %d\n", range);
#pragma omp parallel
	{
#pragma omp for
		for (int i = 2; i <= MAX; i += range)
		{
			//int threadId = omp_get_thread_num();
			//printf("Thread id: %d, ", threadId);
			int upperBound = (i + range > MAX) ? MAX : i + range;
			for (int j = 2; j <= squareRoot; j++)
			{
				if (matrix[j - 2] == true) {
					for (int k = j * 2; k <= upperBound; k += j) {
						matrix[k - 2] = false;
					}
				}
			}
		}
	}
}

void parallelEratostenesDomainWithPrimes(bool* matrix) {
	cleanMatrix(matrix, true, MAX - 2);
	int squareRoot = sqrt(MAX);
	int secondSquareRoot = sqrt(squareRoot);
	int numThreads = THREAD_NUM;
	int* primes = new int[squareRoot];
	int primeCounter = 0;
	int range = squareRoot / (numThreads * DOMAIN_MULTIPLIER);
#pragma omp parallel
	{
		// najpierw watki dziela sie praca wykonywana na przedziale 2 .. sqrt(MAX) 
#pragma omp for
		for (int i = 2; i <= squareRoot; i += range)
		{
			//int threadId = omp_get_thread_num();
			//printf("Thread id: %d, ", threadId);
			int upperBound = (i + range > squareRoot) ? squareRoot : i + range;
			for (int j = 2; j <= secondSquareRoot; j++)
			{
				if (matrix[j - 2] == true) {
					for (int k = j * 2; k <= upperBound; k += j) {
						matrix[k - 2] = false;
					}
				}
			}
		}
#pragma omp single
		{
			for (int i = 2; i < squareRoot; i++) {
				if (matrix[i - 2] == true) {
					primes[primeCounter] = i;
					primeCounter++;
				}
			}
			// watki dziela sie pozostala praca na przedziale sqrt(MAX) .. MAX
			range = (MAX - squareRoot) / (numThreads * DOMAIN_MULTIPLIER);
		}
#pragma omp for
		for (int i = squareRoot + 1; i <= MAX; i += range)
		{
			int upperBound = (i + range > MAX) ? MAX : i + range;
			for (int j = 0; j < primeCounter; j++)
			{
				for (int k = primes[j] * 2; k <= upperBound; k += primes[j]) {
					matrix[k - 2] = false;
				}
			}
		}
	}
	delete[] primes;
}


int main(int argc, char* argv[])
{

	bool* matrix;
	if (eratostenes) { // usuwamy 2, gdyz nie przetwarzamy 0 i 1
		matrix = new bool[MAX + 1 - 2];
	}
	else {
		matrix = new bool[MAX - MIN + 1];
	}


	double start, stop;
	int i;

	omp_set_num_threads(THREAD_NUM);

	start = omp_get_wtime();

	//--- POCZĄTEK PRZETWARZANIA ---

	//{ SEQ_DIV, SEQ_DIV_PRIM, SEQ_ERA, SEQ_ERA_PRIM, PAR_DIV, 
	// PAR_DIV_PRIM, PAR_ERA_FUN, PAR_ERA_FUN_PRIM, PAR_ERA_DOM, 
	// PAR_ERA_DOM_PRIM, PAR_DIV_2 }
	funcPtr[SEQ_DIV_PRIM](matrix);

	//--- KONIEC PRZETWARZANIA ---

	stop = omp_get_wtime();

	printf("Czas: %f\n", (double)(stop - start));
	if (printPrimes)
		showResult(matrix);

	delete[] matrix;
	return 0;
}


// funkcje pomocnicze
void showResult(bool* matrix)
{

	std::ofstream f_results;
	f_results.open("wyniki.txt", std::ios::app);


	int counter = 0;
	int upperLimit = MAX - MIN;
	if (eratostenes) {
		upperLimit = MAX - 2;
	}
	int i = 0;
	if (eratostenes) {
		i = MIN - 2;
	}
	for (i; i <= upperLimit; i++) {
		if (matrix[i]) {
			counter++;
			if (eratostenes) {
				printf("%d, ", i + 2);
			}
			else {
				printf("%d, ", i + MIN);
			}
			if (counter % 10 == 0) {
				printf("\n");
			}
		}
	}

	printf("\n\nW zakresie od %d do %d znaleziono %d liczb pierwszych\n", MIN, MAX, counter);

	time_t result = time(NULL);

	char str[26];
	ctime_s(str, sizeof str, &result);
	f_results << str;
	f_results << "Algorytm: ";
	f_results << "?";
	f_results << "\n";
	f_results << "Liczba watkow: ";
	f_results << THREAD_NUM;
	f_results << "\n";
	f_results << "W zakresie od ";
	f_results << MIN;
	f_results << " do ";
	f_results << MAX;
	f_results << " znaleziono ";
	f_results << counter;
	f_results << " liczb pierwszych\n\n";


	f_results.close();

}
