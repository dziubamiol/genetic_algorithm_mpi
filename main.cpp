#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <chrono>
#include <algorithm>
#include "mpi.h"
#define POPULATION_SIZE 2000

#define SIZE 100
#define NTHREADS 4
#define OFFSET (SIZE / NTHREADS)
/*
struct Timer {
public:
    std::chrono::high_resolution_clock::time_point start_time;

    void start_timer() {
        start_time = std::chrono::steady_clock::now();
    }

    double get_time_elapsed() {
        std::chrono::high_resolution_clock::time_point now = std::chrono::steady_clock::now();;
        std::chrono::duration<double> time_span =  now - start_time;
        return time_span.count();
    }

};
*/
/* --- Calculation functions --- */

double target_function(double x) {
    return std::sin(x);
}

double get_random(double start, double end, std::mt19937 &rand) {
    std::uniform_real_distribution<> dist(start, end);

    return dist(rand);
}

unsigned int get_random_int(double end, std::mt19937 &rand) {
    std::uniform_int_distribution<> dist(0, end);

    return dist(rand);
}

struct XSort {
    bool operator() (const double &l, const double &r) const {
        return target_function(r) < target_function(l);
    }
};

/* --- Genetic algorithms functions --- */

void mutate(std::vector<double> &population, std::mt19937 &rand) {
    unsigned int population_size = population.size();

    for (unsigned int i = 0; i < (population_size * 2 / 10); i++) {
        unsigned int person = get_random_int(population_size, rand);

        population[person] = population[person] + get_random(-0.125, 0.125, rand);
    }
}

void spawn(std::vector<double> &population, std::mt19937 &rand) {
    unsigned int population_size = population.size();

    for (unsigned int i = 0; i < population_size; i++) {
        double father = population[i];
        double mother = population[get_random_int(population_size, rand)];

        population.push_back((father + mother) / 2);
    }

    while (population.size() < 500) {
        double father = population[get_random_int(population_size, rand)];
        double mother = population[get_random_int(population_size, rand)];

        population.push_back((father + mother) / 2);
    }
}

void kill(std::vector<double> &population) {
    population.erase(population.begin() + population.size() / 2, population.end());
}

double get_effectiveness(std::vector<double> &population) {
    double sum = 0;

    for (double i : population) {
        sum = target_function(i);
    }

    return sum / population.size();
}

/* --- Region calculation functions --- */

double calculate_region(double start, double end, std::mt19937 &rand) {
    std::vector<double> population(POPULATION_SIZE);

    for (int i = 0; i < POPULATION_SIZE; i++) {
        population.push_back(get_random(start, end, rand));
    }

    double prev_effectiveness = get_effectiveness(population);

    while (true) {
        double effectiveness = get_effectiveness(population);

        if (std::abs(effectiveness - prev_effectiveness) < 0.0025) {
            // return here best of population
            std::sort(population.begin(), population.end(), XSort());

            return population[0];
        }
        prev_effectiveness = effectiveness;

        // sort here
        std::sort(population.begin(), population.end(), XSort());
        kill(population);

        mutate(population, rand);
        spawn(population, rand);
    }
 }

int main(int argc, char *argv[]) {
    int rc, i, numprocs, myid;	
    std::mt19937 rand { std::random_device{}() };
    std::vector<double> picks;
     
    if (rc = MPI_Init(&argc, &argv)) {
	std::cout << "Startup failure, aborting!\n";
        MPI_Abort(MPI_COMM_WORLD, rc);	
    } 
    
    MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD,&myid); 
    
    if (myid == 0) {
        std::cout << "Main process has been started\n";
	std::cout << "Number of processes " << numprocs << "\n";
	std::cout << "OFFSET " << OFFSET << "\n";
	picks.resize(1000);
    }

    std::vector<double> partialPicks;
   
    double *sendbuf;
    double rbuf[OFFSET];

    sendbuf = (double *)malloc(numprocs * OFFSET * sizeof(double));

    MPI_Scatter(sendbuf, OFFSET, MPI_DOUBLE, rbuf, OFFSET, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    {
    	for (int i = (myid * OFFSET); i < (myid + 1) * OFFSET; i++) {
        	double pick = calculate_region(i * 2 * 3.14, i * 2 * 3.14 + 2 * 3.14, rand);
                int index = i - (myid * OFFSET);
		rbuf[index] = pick;
    	}
    }

    MPI_Gather(rbuf, OFFSET, MPI_DOUBLE, sendbuf, OFFSET, MPI_DOUBLE, 0, MPI_COMM_WORLD);	

    if (myid == 0) {
	std::cout << "Receive data\n";    
	for (int i = 0; i < numprocs * OFFSET; i++) {
	        std::cout << "Pick #" << i << " x: " << sendbuf[i] << "\n";
	}
    }

    
    MPI_Finalize();
    return 0;
}
