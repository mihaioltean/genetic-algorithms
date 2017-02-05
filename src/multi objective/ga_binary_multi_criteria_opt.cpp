// Genetic Algorithm with binary encoding for multicriteria optimization
// Steady state - one population with replacement (the new individuals will replace dominated ones)

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

#include "lista_voidp.h"
//--------------------------------------------------------------------
struct b_chromosome{
	               // also called individual, or potential solution
	char *x;       // an array of genes; one for each dimension
	double f[2];   // the quality of an individual for each objective function (we work with 2 objectives)
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
	dest->f[0] = source.f[0];
	dest->f[1] = source.f[1];
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
//---------------------------------------------------------------------------
// function to be optimized
//---------------------------------------------------------------------------
//ZDT1
//---------------------------------------------------------------------------
double f1(double *p, int num_dimensions)
{
	return p[0];
}
//---------------------------------------------------------------------------
double f2(double* p, int num_dimensions)
{
	// test function T1

	double sum = 0;
	for (int i = 1; i < num_dimensions; i++)
		sum += p[i];
	double g = 1 + (9 * sum) / (double)(num_dimensions - 1);

	double h = 1 - sqrt(f1(p, num_dimensions) / g);

	return g * h;
}
//---------------------------------------------------------------------------
void compute_fitness(b_chromosome *c, int num_dims, int num_bits_per_dimension, double min_x, double max_x)
{
	// first of all transform from binary encoding to real encoding
	// each consecutive num_bits_per_dimension are translated in base 10 and then in [min_x, max_x] interval
	double *x_real = new double[num_dims];
	for (int i = 0; i < num_dims; i++) 
		x_real[i] = binary_to_real(c->x + i * num_bits_per_dimension, num_bits_per_dimension, min_x, max_x);

	c->f[0] = f1(x_real, num_dims);
	c->f[1] = f2(x_real, num_dims);

	delete[] x_real;
}
//--------------------------------------------------------------------
void mutation(b_chromosome c, int num_dims, int num_bits_per_dimension, double pm)
// mutate the individual
{
	// mutate each symbol with the same pm probability
	int _length = num_dims * num_bits_per_dimension;
	for (int i = 0; i < _length; i++) {
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
	for (int i = 0; i < _length; i++) {
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
void print_chromosome(b_chromosome *c, int num_dims, int num_bits_per_dimension, double min_x, double max_x)
{
	printf("x = (");

	for (int i = 0; i < num_dims; i++) {
		double x_real = binary_to_real(c->x + i * num_bits_per_dimension, num_bits_per_dimension, min_x, max_x);
		printf("%lf ", x_real);
	}
	printf(") ");
	printf("fitness = %lf %lf\n", c->f[0], c->f[0]);
}
//--------------------------------------------------------------------
bool dominates(double* p1, double* p2) // returns true if p1 dominates p2
{
	for (int i = 0; i < 2; i++)
		if (p2[i] < p1[i])
			return false;
	for (int i = 0; i < 2; i++)
		if (p1[i] < p2[i])
			return true;
	return false;
}
//---------------------------------------------------------------------------
void sort_list(TLista &nondominated)
{
	bool sorted = false;
	while (!sorted) {
		sorted = true;
		for (node_double_linked * node_p = nondominated.head; node_p->next; node_p = node_p->next) {
			double *p = (double*)nondominated.GetCurrentInfo(node_p);
			double *p_next = (double*)nondominated.GetCurrentInfo(node_p->next);
			if (p[0] > p_next[0]) {
				void *tmp_inf = node_p->inf;
				node_p->inf = node_p->next->inf;
				node_p->next->inf = tmp_inf;
				sorted = false;
			}
		}
	}
}
//---------------------------------------------------------------------------
double compute_hypervolume(b_chromosome* population, int pop_size, TLista& nondominated)
{
	// create a list with nondominated
	nondominated.Clear();
	nondominated.Add(population[0].f);

	for (int i = 1; i < pop_size; i++) {
		node_double_linked * node_p = nondominated.head;
		bool dominated = false;
		while (node_p) {
			double *p = (double*)nondominated.GetCurrentInfo(node_p);
			if (dominates(p, population[i].f)) {
				dominated = true;
				break;
			}
			else
				if (dominates(population[i].f, p))
					node_p = nondominated.DeleteCurrent(node_p);
				else// move to the next one
					node_p = node_p->next;
		}
		if (!dominated)
			nondominated.Add(population[i].f);
	}

	// compute the distance to front
	double reference[2] = { 11, 11 };

	sort_list(nondominated);

	double hyper_volume = 0;
	for (node_double_linked * node_p = nondominated.head; node_p->next; node_p = node_p->next) {
		double *p = (double*)nondominated.GetCurrentInfo(node_p);
		double *p_next = (double*)nondominated.GetCurrentInfo(node_p->next);
		hyper_volume += (p_next[0] - p[0]) * (reference[1] - p[1]);
	}
	double *p = (double*)nondominated.GetTailInfo();
	hyper_volume += (reference[0] - p[0]) * (reference[1] - p[1]);

	return hyper_volume;
}
//---------------------------------------------------------------------------
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

	FILE *f = fopen("pareto.txt", "w");

	TLista nondominated; 

	double hv = compute_hypervolume(population, pop_size, nondominated);

	printf("generation %d hypervolume = %lf num nondominated = %d\n", 0, hv, nondominated.count);
	fprintf(f, "%d ", nondominated.count);
	for (node_double_linked * node_p = nondominated.head; node_p->next; node_p = node_p->next) {
		double *p = (double*)nondominated.GetCurrentInfo(node_p);
		fprintf(f, "%lf %lf ", p[0], p[1]);
	}
	fprintf(f, "\n");

	for (int g = 1; g < num_gens; g++) {
		for (int k = 0; k < pop_size; k += 2) {
			// choose the parents using binary tournament
			int r1 = rand() % pop_size;
			int r2 = rand() % pop_size;

			// crossover
			double p = rand() / double(RAND_MAX);
			if (p < pcross)
				uniform_crossover(population[r1], population[r2], offspring1, offspring2, num_dims, num_bits_per_dimension);
			else {
				copy_individual(&offspring1, population[r1], num_dims, num_bits_per_dimension);
				copy_individual(&offspring2, population[r2], num_dims, num_bits_per_dimension);
			}
			// mutate the result and compute its fitness
			mutation(offspring1, num_dims, num_bits_per_dimension, pm);
			compute_fitness(&offspring1, num_dims, num_bits_per_dimension, min_x, max_x);
			mutation(offspring2, num_dims, num_bits_per_dimension, pm);
			compute_fitness(&offspring2, num_dims, num_bits_per_dimension, min_x, max_x);

			// are offspring dominating something?

			for (int i = 0; i < pop_size; i++)
				if (dominates(offspring1.f, population[i].f)) {
					copy_individual(&population[i], offspring1, num_dims, num_bits_per_dimension);
					break;
				}

			for (int i = 0; i < pop_size; i++)
				if (dominates(offspring2.f, population[i].f)) {
					copy_individual(&population[i], offspring2, num_dims, num_bits_per_dimension);
					break;
				}
		}
		
		double hv = compute_hypervolume(population, pop_size, nondominated);

		printf("generation %d hypervolume = %lf num nondominated = %d\n", g, hv, nondominated.count);

		fprintf(f, "%d ", nondominated.count);
		for (node_double_linked * node_p = nondominated.head; node_p->next; node_p = node_p->next) {
			double *p = (double*)nondominated.GetCurrentInfo(node_p);
			fprintf(f, "%lf %lf ", p[0], p[1]);
		}
		fprintf(f, "\n");

		printf("");
	}
	// free memory
	free(offspring1.x);
	free(offspring2.x);

	for (int i = 0; i < pop_size; i++)
		free(population[i].x);
	free(population);

	fclose(f);
}
//--------------------------------------------------------------------
int main(void)
{
	int pop_size = 100;    // the number of individuals in population 
	// must be an even number
	int num_gens = 100;   // the number of generations
	double pm = 0.001;      // mutation probability
	double pcross = 0.9;  // crossover probability


	int num_dims = 30;  // number of dimensions of the function to be optimized
	double min_x = 0;  // definition domain for functions
	double max_x = 1;

	int num_bits_per_dimension = 30; // 10 bits precision

	srand(0);

	start_steady_state_ga(pop_size, num_gens, num_dims, num_bits_per_dimension, pcross, pm, min_x, max_x);

	printf("Press enter ...");
	getchar();

	return 0;
}
//--------------------------------------------------------------------