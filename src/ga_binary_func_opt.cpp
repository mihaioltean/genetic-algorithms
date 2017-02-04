// Genetic Algorithm with binary encoding for function optimization problems
// (finding the minum or maximum of a given function)
// Steady state - one population with replacement

// author: Mihai Oltean
// web: www.tcreate.org
// github.com/mihaioltean
// email: mihai.oltean@gmail.com


// MIT License

// compiled with Visual Studio 2013 Express Edition
// last modified on: 2017.02.04

//--------------------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//--------------------------------------------------------------------
struct b_chromosome{
	// also called individual, or potential solution
	char *x;      // an array of genes; one for each dimension
	double fitness; // the quality of an individual
};
//--------------------------------------------------------------------
void generate_random(b_chromosome c, int num_dims, int num_bits_per_dimension)
{
	// num_dims is the number of dimensions of the function to be optimized
	// min_x, max_x is the range for variables
	int _length = num_dims * num_bits_per_dimension;
	for (int i = 0; i < _length; i++)
		c.x[i] = rand() % 2;
}
//--------------------------------------------------------------------
void copy_individual(b_chromosome *dest, b_chromosome source, int num_dims, int num_bits_per_dimension)
{
	int _length = num_dims * num_bits_per_dimension;
	for (int i = 0; i < _length; i++)
		dest->x[i] = source.x[i];
	dest->fitness = source.fitness;
}
//--------------------------------------------------------------------
double binary_to_real(char* b_string, int num_bits_per_dimension, double min_x, double max_x)
{
	// transform a binary string of num_bits_per_dimension size into a real number in [min_x ... max_x] interval
	double x_real = 0;
	for (int j = 0; j < num_bits_per_dimension; j++)
		x_real = x_real * 2 + (int)b_string[j]; // now I have them in interval [0 ... 2 ^ num_bits_per_dimension - 1]
	x_real /= ((1 << num_bits_per_dimension) - 1); // now I have them in [0 ... 1] interval
	x_real *= (max_x - min_x); // now I have them in [0 ... max_x - min_x] interval
	x_real += min_x; // now I have them in [min_x ... max_x] interval

	return x_real;
}
//--------------------------------------------------------------------
void compute_fitness(b_chromosome *c, int num_dims, int num_bits_per_dimension, double min_x, double max_x)
{
	// first of all transform from binary encoding to real encoding
	// each consecutive num_bits_per_dimension are translated in base 10 and then in [min_x, max_x] interval
	double *x_real = new double[num_dims];
	for (int i = 0; i < num_dims; i++) 
		x_real[i] = binary_to_real(c->x + i * num_bits_per_dimension, num_bits_per_dimension, min_x, max_x);

	// rastrigin
	c->fitness = 10 * num_dims;
	for (int i = 0; i < num_dims; i++)
		c->fitness += x_real[i] * x_real[i] - 10 * cos(2 * 3.1415 * x_real[i]);

	delete[] x_real;
}
//--------------------------------------------------------------------
void mutation(b_chromosome c, int num_dims, int num_bits_per_dimension, double pm)
// mutate the individual
{
	// mutate each symbol with the same pm probability

	int _length = num_dims * num_bits_per_dimension;
	for (int i = 0; i < _length; i++){
		double p = rand() / (double)RAND_MAX;
		if (p < pm)
			c.x[i] = 1 - c.x[i];
	}
}
//--------------------------------------------------------------------
void one_cut_point_crossover(b_chromosome parent1, b_chromosome parent2, b_chromosome offspring1, b_chromosome offspring2, int num_dims, int num_bits_per_dimension)
{
	int pct;
	pct = 1 + rand() % (num_dims * num_bits_per_dimension - 1);
	// the cutting point should be after first gene 
	//and before the last one
	for (int i = 0; i < pct; i++) {
		offspring1.x[i] = parent1.x[i];
		offspring2.x[i] = parent2.x[i];
	}
	for (int i = pct; i < num_dims * num_bits_per_dimension; i++) {
		offspring1.x[i] = parent2.x[i];
		offspring2.x[i] = parent1.x[i];
	}
}
//--------------------------------------------------------------------
void uniform_crossover(b_chromosome parent1, b_chromosome parent2, b_chromosome offspring1, b_chromosome offspring2, int num_dims, int num_bits_per_dimension)
{
	// uniform crossover can also be used
	// for each gene we decide randomly where it goes
	// (to the first or second offspring)
	int _length = num_dims * num_bits_per_dimension;
	for (int i = 0; i < _length; i++){
		if (rand() % 2) {// flip
			offspring1.x[i] = parent2.x[i];
			offspring2.x[i] = parent1.x[i];
		}
		else {
			offspring1.x[i] = parent1.x[i];
			offspring2.x[i] = parent2.x[i];
		}
	}
}
//--------------------------------------------------------------------
int sort_function(const void *a, const void *b)
{
	if (((b_chromosome *)a)->fitness >((b_chromosome *)b)->fitness)
		return 1;
	else
		if (((b_chromosome *)a)->fitness < ((b_chromosome *)b)->fitness)
			return -1;
		else
			return 0;
}
//--------------------------------------------------------------------
void print_chromosome(b_chromosome *c, int num_dims, int num_bits_per_dimension, double min_x, double max_x)
{
	printf("x = (");

	for (int i = 0; i < num_dims; i++) {
		double x_real = binary_to_real(c->x + i * num_bits_per_dimension, num_bits_per_dimension, min_x, max_x);
		printf("%lf ", x_real);
	}
	printf(") ");
	printf("fitness = %lf\n", c->fitness);
}
//--------------------------------------------------------------------
int tournament_selection(int tournament_size, b_chromosome *pop, int pop_size)
{
	// randomly pick several individuals
	// and return the best of them
	int selected_index;
	selected_index = rand() % pop_size;
	for (int i = 1; i < tournament_size; i++) {
		int r = rand() % pop_size;
		selected_index = pop[r].fitness < pop[selected_index].fitness ? r : selected_index;
	}
	return selected_index;
}
//--------------------------------------------------------------------
void start_steady_state_ga(int pop_size, int num_gens, int num_dims, int num_bits_per_dimension, double pcross, double pm, double min_x, double max_x)
// Steady-State genetic algorithm
// each step:
// pick 2 parents, mate them, 
// mutate the offspring
// and replace the worst in the population
// (only if offspring are better)
{
	// allocate memory
	b_chromosome *population;
	population = (b_chromosome*)malloc(pop_size * sizeof(b_chromosome));
	for (int i = 0; i < pop_size; i++)
		population[i].x = (char*)malloc(num_dims * num_bits_per_dimension);

	b_chromosome offspring1, offspring2;
	offspring1.x = (char*)malloc(num_dims * num_bits_per_dimension);
	offspring2.x = (char*)malloc(num_dims * num_bits_per_dimension);

	// initialize
	for (int i = 0; i < pop_size; i++) {
		generate_random(population[i], num_dims, num_bits_per_dimension);
		compute_fitness(&population[i], num_dims, num_bits_per_dimension, min_x, max_x);
	}

	qsort((void *)population, pop_size, sizeof(population[0]), sort_function);

	printf("generation 0\n");
	// print the best from generation 0
	print_chromosome(&population[0], num_dims, num_bits_per_dimension, min_x, max_x);

	for (int g = 1; g < num_gens; g++) {
		for (int k = 0; k < pop_size; k += 2) {
			// choose the parents using binary tournament
			int r1 = tournament_selection(2, population, pop_size);
			int r2 = tournament_selection(2, population, pop_size);
			// crossover
			double p = rand() / double(RAND_MAX);
			if (p < pcross)
				one_cut_point_crossover(population[r1], population[r2], offspring1, offspring2, num_dims, num_bits_per_dimension);
			else {
				copy_individual(&offspring1, population[r1], num_dims, num_bits_per_dimension);
				copy_individual(&offspring2, population[r2], num_dims, num_bits_per_dimension);
			}
			// mutate the result and compute its fitness
			mutation(offspring1, num_dims, num_bits_per_dimension, pm);
			compute_fitness(&offspring1, num_dims, num_bits_per_dimension, min_x, max_x);
			mutation(offspring2, num_dims, num_bits_per_dimension, pm);
			compute_fitness(&offspring2, num_dims, num_bits_per_dimension, min_x, max_x);

			// are offspring better than the worst ?
			if (offspring1.fitness < population[pop_size - 1].fitness) {

				copy_individual(&population[pop_size - 1], offspring1, num_dims, num_bits_per_dimension);
				qsort((void *)population, pop_size, sizeof(population[0]), sort_function);
			}
			if (offspring2.fitness < population[pop_size - 1].fitness) {
				copy_individual(&population[pop_size - 1], offspring2, num_dims, num_bits_per_dimension);
				qsort((void *)population, pop_size, sizeof(population[0]), sort_function);
			}
		}
		printf("generation %d\n", g);
		print_chromosome(&population[0], num_dims, num_bits_per_dimension, min_x, max_x);
	}
	// free memory
	free(offspring1.x);
	free(offspring2.x);

	for (int i = 0; i < pop_size; i++)
		free(population[i].x);
	free(population);
}
//--------------------------------------------------------------------
int main(void)
{
	int pop_size = 100;    // the number of individuals in population 
	// must be an even number
	int num_gens = 100;   // the number of generations
	double pm = 0.03;      // mutation probability
	double pcross = 0.9;  // crossover probability


	int num_dims = 5;  // number of dimensions of the function to be optimized
	double min_x = -10;  // definition domain for functions
	double max_x = 10;

	int num_bits_per_dimension = 30; // 10 bits precision

	srand(0);

	start_steady_state_ga(pop_size, num_gens, num_dims, num_bits_per_dimension, pcross, pm, min_x, max_x);

	printf("Press enter ...");
	getchar();

	return 0;
}
//--------------------------------------------------------------------