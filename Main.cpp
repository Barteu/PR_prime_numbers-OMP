#include <stdio.h>
#include <time.h>
#include <omp.h>
#include <cmath>
#include <iostream>
#include <fstream>
#include <ctime>

//############### PARAMETRY DO USTAWIENIA #################

// MIN >= 2
const int MIN = 2;
const int MAX = 350000000; // 175000000 350000000

#define THREAD_NUM 4
const bool printPrimes = false;   // Czy wypisac wyznaczone liczby pierwsze wraz z ich ilością

// W podejsciu domenowym watki dziela prace na THREAD_NUM * DOMAIN_MULTIPLIER przedzialow
const int DOMAIN_MULTIPLIER = 2;

//##########################################################

enum version { SEQ_DIV, SEQ_DIV_PRIM, SEQ_ERA, SEQ_ERA_OPT, PAR_DIV, PAR_DIV_PRIM, PAR_ERA_FUN, PAR_ERA_DOM, PAR_ERA_FUN_LOC };

// prototypy funkcji pomocniczych
void showResult(bool* matrix);
void cleanMatrix(bool* matrix, bool val, int maxRange);
bool isPrime(int n);
bool isPrimeExtended(int n, int* primes, int prime_counter);
// prototypy funkcji do przypisania dla tablicy wskaznikow do funkcji
void sequentialDivision(bool* matrix);
void sequentialDivisionWithPrimes(bool* matrix);
void sequentialEratostenes(bool* matrix);
void sequentialEratostenesOpt(bool* matrix);
void parallelDivision(bool* matrix);
void parallelDivisionWithPrimes(bool* matrix);
void parallelEratostenesFunctional(bool* matrix);
void parallelEratostenesFunctionalLocal(bool* matrix);
void parallelEratostenesDomain(bool* matrix);

void (*funcPtr[])(bool*) = { sequentialDivision, sequentialDivisionWithPrimes, sequentialEratostenes, sequentialEratostenesOpt, parallelDivision, parallelDivisionWithPrimes, parallelEratostenesFunctional,
							   parallelEratostenesDomain, parallelEratostenesFunctionalLocal};


int main(int argc, char* argv[])
{

	bool* matrix;
	matrix = new bool[MAX - MIN + 1];
	double start, stop;
	omp_set_num_threads(THREAD_NUM);

	start = omp_get_wtime();

	//-------- POCZĄTEK PRZETWARZANIA --------- 

	//SEQ_DIV, SEQ_DIV_PRIM, SEQ_ERA, SEQ_ERA_OPT, PAR_DIV, PAR_DIV_PRIM, PAR_ERA_FUN, PAR_ERA_DOM, PAR_ERA_FUN_LOC
	funcPtr[SEQ_ERA_OPT](matrix);    // <<------ w tym miejscu prosze wybrac uzywany algorytm
	
	//--------- KONIEC PRZETWARZANIA ----------

	stop = omp_get_wtime();

	printf("Czas: %f\n", (double)(stop - start));
	if(printPrimes)
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
	for (int i = 0 ; i <= upperLimit; i++) {
		if (matrix[i]) {
			counter++;
			printf("%d, ", i + MIN);
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
	printf("SEQ DIV\n");
	cleanMatrix(matrix, false, MAX - MIN);
	for (int i = MIN; i <= MAX; i++) {
		if (isPrime(i)) {
			matrix[i - MIN] = true;
		}
	}
}

void sequentialDivisionWithPrimes(bool* matrix) {
	printf("SEQ DIV PRIM\n");
	cleanMatrix(matrix, false, MAX - MIN);

	int squareRoot = sqrt(MAX);
	int* primes = new int[squareRoot];
	int prime_counter = 0;

	for (int i = 2; i <= squareRoot; i++)
	{
		if (isPrimeExtended(i, primes, prime_counter)) {
			if (prime_counter < squareRoot) {
				primes[prime_counter] = i;
				prime_counter++;

				if (i >= MIN)
				{
					matrix[i - MIN] = true;
				}
			}
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
	}

	delete[] primes;

}

void sequentialEratostenes(bool* matrix) {
	printf("SEQ ERA\n");
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

// FIN ver
void sequentialEratostenesOpt(bool* matrix) {
	printf("SEQ ERA OPT\n");
	cleanMatrix(matrix, true, MAX - MIN);
	int squareRoot = sqrt(MAX);
	int secondSquareRoot = sqrt(squareRoot);
	// W celu optymalizacji algorytmu
	bool* matrixSqrt = new bool[squareRoot];
	for (int i = 0; i < squareRoot; i++) {
		matrixSqrt[i] = true;
	}
	// Dla 2 ... sqrt(MAX)
	for (int i = 2; i < secondSquareRoot; i++) {
		if (matrixSqrt[i - 2] == false) {
			continue;
		}
		for (int j = i * i; j <= squareRoot; j += i) {
			matrixSqrt[j - 2] = false;
			if (j >= MIN)
			{
				matrix[j - MIN] = false;
			}
		}
	}
	// Dla MIN ... MAX
	for (int counter = 2; counter < squareRoot; counter++) {
		//printf("%d\n", counter);
		if (matrixSqrt[counter - 2] == false) {
			continue;
		}
		int i = counter * counter;
		if (i < MIN && (counter * ceil(float(MIN) / float(counter)))>counter * counter) {
			i = counter * ceil(float(MIN) / float(counter));
		}
		for (i; i <= MAX; i += counter) {
			matrix[i - MIN] = false;
		}
	}
	delete[] matrixSqrt;
}

// FIN ver
void sequentialEratostenesWithPrimes(bool* matrix) {
	printf("SEQ ERA PRIM\n");
	cleanMatrix(matrix, true, MAX - MIN);
	int squareRoot = sqrt(MAX);
	bool* matrixSqrt = new bool[squareRoot];
	for (int i = 0; i < squareRoot; i++) {
		matrixSqrt[i] = true;
	}
	int* primes = new int[squareRoot];
	int primeCounter = 0;
	int secondSquareRoot = sqrt(squareRoot);
	for (int counter = 2; counter <= secondSquareRoot; counter++)
	{
		if (matrixSqrt[counter - 2] == false) {
			continue;
		}
		for (int i = counter * counter; i <= squareRoot; i += counter) {
			matrixSqrt[i - 2] = false;
			if (i >= MIN) {
				matrix[i - MIN] = false;
			}
		}
	}
	for (int i = 2; i <= squareRoot; i++) {
		if(matrixSqrt[i - 2] == true){
			primes[primeCounter] = i;
			primeCounter++;
		}
	}
	for (int i = 0; i < primeCounter; i++) {
		int j = primes[i] * primes[i];
		if (j < MIN && (primes[i] * ceil(float(MIN) / primes[i]))>primes[i]* primes[i]) {
			j = primes[i] * ceil(float(MIN) / primes[i]);
		}
		for (j; j <= MAX; j += primes[i]) {
			matrix[j - MIN] = false;
		}
	}
	delete[] matrixSqrt;
	delete[] primes;
}

// Wersje równoległe
void parallelDivision(bool* matrix) {
	printf("PAR DIV\n");
	cleanMatrix(matrix, false, MAX - MIN);
	int squareRoot = sqrt(MAX);
	int chunk_size = ceil(MAX / 100.0);
#pragma omp parallel
	{
#pragma omp for schedule(dynamic,chunk_size)
		for (int i = MIN; i <= MAX; i++) {
			if (isPrime(i)) {
				matrix[i - MIN] = true;
			}
		}
	}
}

void parallelDivisionWithPrimes(bool* matrix) {
	printf("PAR DIV PRIM\n");
	cleanMatrix(matrix, false, MAX - MIN);
	int squareRoot = sqrt(MAX);

	int* primes = new int[int(squareRoot / 2) + 10];
	int prime_counter = 0;

	bool* matrixUnderMin = new bool[MIN];
	cleanMatrix(matrixUnderMin, false, MIN - 1);

	int chunk_size = ceil(squareRoot / 1000.0);
	int chunk2_size = ceil((MAX - MIN) / 1000.0);
	
#pragma omp parallel
	{
		int threadId = omp_get_thread_num();

#pragma omp for schedule(dynamic,chunk_size)
		for (int i = 2; i <= squareRoot; i++) {

			if (isPrime(i)) {
				if (i < MIN) {
					matrixUnderMin[i] = true;
				}
				else {
					matrix[i - MIN] = true;
				}
			}
		}

#pragma omp master
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

#pragma omp for schedule(dynamic,chunk2_size)
		for (int i = (squareRoot + 1 > MIN) ? squareRoot + 1 : MIN; i <= MAX; i++) {
			if (isPrimeExtended(i, primes, prime_counter)) {
				matrix[i - MIN] = true;
			}
		}


	}

	delete[] primes;
	delete[] matrixUnderMin;
}


// FIN ver
void parallelEratostenesFunctional(bool* matrix) {
	printf("PAR ERA FUN\n");
	cleanMatrix(matrix, true, MAX - MIN);
	int squareRoot = sqrt(MAX);
	int secondSquareRoot = sqrt(squareRoot);
	// W celu optymalizacji algorytmu
	bool* matrixSqrt = new bool[squareRoot];
	for (int i = 0; i < squareRoot; i++) {
		matrixSqrt[i] = true;
	}
	int chunk_size = ceil(squareRoot / 1000.0);
	int chunk2_size = ceil((MAX - MIN) / 1000.0);

	// Dla 2 ... sqrt(MAX)
#pragma omp parallel
	{
#pragma omp for schedule(dynamic,chunk_size)
		for (int i = 2; i < secondSquareRoot; i++) {
			if (matrixSqrt[i - 2] == false) {
				continue;
			}
			for (int j = i * i; j <= squareRoot; j += i) {
				matrixSqrt[j - 2] = false;
				if (j >= MIN)
				{
					matrix[j - MIN] = false;
				}
			}
		}
		// Dla MIN ... MAX
#pragma omp for schedule(dynamic,chunk2_size)
		for (int counter = 2; counter < squareRoot; counter++) {
			if (matrixSqrt[counter - 2] == false) {
				continue;
			}
			int i = counter * counter;
			if (i < MIN && (counter * ceil(float(MIN) / float(counter)))>counter* counter) {
				i = counter * ceil(float(MIN) / float(counter));
			}
			for (i; i <= MAX; i += counter) {
				matrix[i - MIN] = false;
			}
		}
	}
	delete[] matrixSqrt;
}

// wersja hybrydowa, podzial funkcjonalny z lokalnoscia
void parallelEratostenesFunctionalLocal(bool* matrix) {
	printf("PAR ERA FUN LOC\n");
	cleanMatrix(matrix, true, MAX - MIN);
	int squareRoot = sqrt(MAX);
	int secondSquareRoot = sqrt(squareRoot);
	// W celu optymalizacji algorytmu
	bool* matrixSqrt = new bool[squareRoot];
	for (int i = 0; i < squareRoot; i++) {
		matrixSqrt[i] = true;
	}
	int chunk_size = ceil(squareRoot / 1000.0);
	int range = 1024;//squareRoot / (numThreads * DOMAIN_MULTIPLIER);
	int range2 = 32768;//squareRoot / (numThreads * DOMAIN_MULTIPLIER);
	// Dla 2 ... sqrt(MAX)
#pragma omp parallel 
	{
		for (int x = 2; x < squareRoot; x += range) {
			int upperBound = (x + range > squareRoot) ? squareRoot : x + range;
#pragma omp for schedule(dynamic,chunk_size)
			for (int i = 2; i < secondSquareRoot; i++) {
				if (matrixSqrt[i - 2] == false) {
					continue;
				}
				int j = i * i;
				if (j < x && (i * ceil(float(x) / float(i)))>i* i) {
					j = i * ceil(float(x) / float(i));
				}
				for (j; j <= upperBound; j += i) {
					matrixSqrt[j - 2] = false;
					if (j >= MIN)
					{
						matrix[j - MIN] = false;
					}
				}
			}
		}

		// Dla MIN ... MAX
for(int x = MIN; x < MAX; x+= range2){
	int upperBound = (x + range> MAX) ? MAX : x + range2;

#pragma omp for schedule(dynamic,chunk_size)
		for (int counter = 2; counter < squareRoot; counter++) {
			if (matrixSqrt[counter - 2] == false) {
				continue;
			}
			int i = counter * counter;
			if (i < x && (counter * ceil(float(x) / float(counter)))>counter* counter) {
				i = counter * ceil(float(x) / float(counter));
			}
			for (i; i <= upperBound; i += counter) {
				matrix[i - MIN] = false;
			}
		}
	}
}
	delete[] matrixSqrt;
}

void parallelEratostenesDomain(bool* matrix) {
	printf("PAR ERA DOM\n");
	cleanMatrix(matrix, true, MAX - MIN);
	int squareRoot = sqrt(MAX);
	int numThreads = THREAD_NUM;
	printf("Num threads: %d\n", numThreads);

	// W celu optymalizacji algorytmu
	bool* matrixSqrt = new bool[squareRoot];
	
	for (int i = 0; i < squareRoot; i++) {
		matrixSqrt[i] = true;
	}
	int secondSquareRoot = sqrt(squareRoot);
	// Dla 2 ... sqrt(MAX)
	int range = 1024;//squareRoot / (numThreads * DOMAIN_MULTIPLIER);
	printf("Range sqrt: %d\n", range);
	//  range2
	//range2 = 32768;// (MAX - MIN) / (numThreads * DOMAIN_MULTIPLIER); // range = 4000 - 64000 gwarantuje dobry MEMORY BOUND 
	int range2 = 32768;
	// W sumie dla wiekszych instacji wystarczy dosc mocno zwiekszac DOMAIN_MULTIPLIER, tak zeby (MAX-MIN)/DOMAIN_MULTIPLIER = ~64000
	printf("Range2 (MAX-MIN): %d\n", range2);
	int chunk2_size = ceil(((MAX - MIN)/range2) / 1000.0);

#pragma omp parallel
	{
#pragma omp for schedule(dynamic,1)
		for (int i = 2; i <= squareRoot; i += range)
		{
			int upperBound = (i + range > squareRoot) ? squareRoot : i + range;
			for (int j = 2; j <= secondSquareRoot; j++)
			{
				if (matrixSqrt[j - 2] == true) {
					int k = j * j;
					if (k < i && (j * ceil(float(i) / j))>j * j) {
						k = j * ceil(float(i) / j);
					}
					for (k; k <= upperBound; k += j) {
						matrixSqrt[k - 2] = false;
						if (k >= MIN) {
							matrix[k - MIN] = false;
						}
					}
				}
			}
		}
		// Dla MIN ... MAX
#pragma omp for schedule(dynamic,chunk2_size)
		for (int i = MIN; i <= MAX; i += range2)
		{
			int upperBound = (i + range2 > MAX) ? MAX : i + range2;
			for (int j = 2; j <= squareRoot; j++)
			{
				if (matrixSqrt[j - 2] == true) {
					int k = j * j;
					if (k < i && (j * ceil(float(i) / j))>j * j) {
						k = j * ceil(float(i) / j);
					}
					for (k; k <= upperBound; k += j) {
						matrix[k - MIN] = false;
					}
				}
			}
		}
	}
	delete[] matrixSqrt;
}


