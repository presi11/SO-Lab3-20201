#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <assert.h>
#include <pthread.h>
#include <semaphore.h>
#include <sys/time.h>

typedef struct Args{
	int curr_index;
	int end_index;
	int n_threads;
	int max_iters;
	int *n_threads_finished;
	int *n_iter;
	double scalar;
	double *X;
	double *Y;
	double *Y_avgs;
}Saxpy_args;

sem_t sem;
sem_t sem1;
pthread_mutex_t lock;
int p = 10000000;
Saxpy_args **create_args(int n_elements, int n_threads, int max_iters, int *n_threads_finished, int *n_iter,
double scalar,  double *X, double *Y, double *Y_avgs);
double calculate_saxpy_n_time(int curr_index, int end_index, int n_time, double scalar, double *X, double *Y);
void update_Y_n_time(int curr_index, int end_index, int n_time, double scalar, double *X, double *Y);
void *calculate_saxpy(void *args);
pthread_t *create_thread_arr(int n_threads, Saxpy_args **args);
void wait_for_threads(pthread_t *threads, int n_threads);

int main(int argc, char* argv[]){
	unsigned int value = 0;
	sem_init(&sem, 0, value);
	sem_init(&sem1, 0, 0);
	// Variables to obtain command line parameters
	unsigned int seed = 1;
	int n_threads = 16;
  	int max_iters = 1000;
  	// Variables to perform SAXPY operation
	double* X;
	double a;
	double* Y;
	double* Y_avgs;
	int i, it;
	// Variables to get execution time
	struct timeval t_start, t_end;
	double exec_time;

	// Getting input values
	int opt;
	while((opt = getopt(argc, argv, ":p:s:n:i:")) != -1){  
		switch(opt){  
			case 'p':  
			printf("vector size: %s\n", optarg);
			p = strtol(optarg, NULL, 10);
			assert(p > 0 && p <= 2147483647);
			break;  
			case 's':  
			printf("seed: %s\n", optarg);
			seed = strtol(optarg, NULL, 10);
			break;
			case 'n':  
			printf("threads number: %s\n", optarg);
			n_threads = strtol(optarg, NULL, 10);
			break;  
			case 'i':  
			printf("max. iterations: %s\n", optarg);
			max_iters = strtol(optarg, NULL, 10);
			break;  
			case ':':  
			printf("option -%c needs a value\n", optopt);  
			break;  
			case '?':  
			fprintf(stderr, "Usage: %s [-p <vector size>] [-s <seed>] [-n <threads number>]\n", argv[0]);
			exit(EXIT_FAILURE);
		}  
	}  
	srand(seed);

	printf("p = %d, seed = %d, n_threads = %d, max_iters = %d\n", \
	 p, seed, n_threads, max_iters);	

	// initializing data
	X = (double*) malloc(sizeof(double) * p);
	Y = (double*) malloc(sizeof(double) * p);
	Y_avgs = (double*) malloc(sizeof(double) * max_iters);

	for(i = 0; i < p; i++){
		X[i] = (double)rand() / RAND_MAX;
		Y[i] = (double)rand() / RAND_MAX;
	}
	
	for(i = 0; i < max_iters; i++){
		Y_avgs[i] = 0.0;
	}
	a = (double)rand() / RAND_MAX;
#ifdef DEBUG
	printf("vector X= [ ");
	for(i = 0; i < p-1; i++){
		printf("%f, ",X[i]);
	}
	printf("%f ]\n",X[p-1]);

	printf("vector Y= [ ");
	for(i = 0; i < p-1; i++){
		printf("%f, ", Y[i]);
	}
	printf("%f ]\n", Y[p-1]);

	printf("a= %f \n", a);
#endif
	/*
	 *	Function to parallelize 
	 */
	
	gettimeofday(&t_start, NULL);
	//SAXPY iterative SAXPY mfunction
	it = 0;
	int *n_threads_finished = (int*) calloc(1, sizeof(int));
	Saxpy_args **args = create_args(p, n_threads, max_iters, n_threads_finished, &it, a,  X, Y, Y_avgs);
	pthread_t *threads;
	threads = create_thread_arr(n_threads, args); // n_threads are created
	wait_for_threads(threads, n_threads);
	gettimeofday(&t_end, NULL);
#ifdef DEBUG
	printf("RES: final vector Y= [ ");
	for(i = 0; i < p-1; i++){
		printf("%f, ", Y[i]);
	}
	printf("%f ]\n", Y[p-1]);
#endif
	
	// Computing execution time
	exec_time = (t_end.tv_sec - t_start.tv_sec) * 1000.0;  // sec to ms
	exec_time += (t_end.tv_usec - t_start.tv_usec) / 1000.0; // us to ms
	printf("Execution time: %f ms \n", exec_time);
	printf("Last 3 values of Y: %f, %f, %f \n", Y[p-3], Y[p-2], Y[p-1]);
	printf("Last 3 values of Y_avgs: %f, %f, %f \n", Y_avgs[max_iters-3], Y_avgs[max_iters-2], Y_avgs[max_iters-1]);
	free(threads);
	return 0;
}

Saxpy_args **create_args(int n_elements, int n_threads, int max_iters, int *n_threads_finished, int *n_iter,
double scalar,  double *X, double *Y, double *Y_avgs){
	Saxpy_args **args = (Saxpy_args**)malloc(n_threads*sizeof(Saxpy_args*));
	int arr_chunk = n_elements/n_threads;
	int n_leftovers = n_elements%n_threads;
	int curr_index = 0;
	int end_index = 0;
	int n_threads_copy = n_threads;
	int n; 
	while(n_threads--){
		args[n_threads] = (Saxpy_args*)malloc(sizeof(Saxpy_args));
		if(n_leftovers){
			n = arr_chunk + 1;
			n_leftovers--;
		}else
			n = arr_chunk;
		args[n_threads] -> curr_index = curr_index;
		end_index = curr_index + n - 1;
		args[n_threads] -> end_index = end_index;
		args[n_threads] -> n_threads = n_threads_copy;
		args[n_threads] -> max_iters = max_iters;
		args[n_threads] -> n_threads_finished = n_threads_finished;
		args[n_threads] -> n_iter = n_iter;
		args[n_threads] -> scalar = scalar;
		args[n_threads] -> X = X;
		args[n_threads] -> Y = Y;
		args[n_threads] -> Y_avgs = Y_avgs;
		curr_index = end_index + 1;
	}
	return args;
}

double calculate_saxpy_n_time(int curr_index, int end_index, int n_time, double scalar, double *X, double *Y){
	double result = 0.0;
	while(curr_index <= end_index){
		result += scalar*n_time*X[curr_index] + Y[curr_index];
		curr_index ++;
	}
	return result;
}

void update_Y_n_time(int curr_index, int end_index, int n_time, double scalar, double *X, double *Y){
	while(curr_index <= end_index){
		Y[curr_index] = scalar*n_time*X[curr_index] + Y[curr_index];
		curr_index++;
	}
}
void *calculate_saxpy(void *args){
	Saxpy_args *saxpy_args = (Saxpy_args*)args;
	int curr_index = saxpy_args -> curr_index;
	int end_index = saxpy_args -> end_index;
	int n_threads = saxpy_args-> n_threads;
	int max_iters = saxpy_args -> max_iters;
	int *n_threads_finished =  saxpy_args -> n_threads_finished;
	int *n_iter = saxpy_args -> n_iter;
	double scalar = saxpy_args -> scalar;
	double *X = saxpy_args -> X;
	double *Y = saxpy_args -> Y;
	double *Y_avgs = saxpy_args -> Y_avgs;
	double result;
	while(*n_iter < max_iters){
		result = calculate_saxpy_n_time(curr_index, end_index, *n_iter+1, scalar, X, Y)/p;
		pthread_mutex_lock(&lock);
		Y_avgs[*n_iter] += result;
		(*n_threads_finished)++;
		if(*n_threads_finished == n_threads){
			(*n_iter)++;
			sem_init(&sem1, 0, 0);
			sem_post(&sem);
		}
		pthread_mutex_unlock(&lock);
		sem_wait(&sem);
		sem_post(&sem);
		pthread_mutex_lock(&lock);
		(*n_threads_finished)--;
		if(*n_threads_finished == 0){
			sem_init(&sem, 0, 0);
			sem_post(&sem1);
		}
		pthread_mutex_unlock(&lock);
		sem_wait(&sem1);
		sem_post(&sem1);		
	}
	update_Y_n_time(curr_index, end_index, max_iters, scalar, X, Y);
	pthread_exit(NULL);
}

pthread_t *create_thread_arr(int n_threads, Saxpy_args **args){
	pthread_t *threads = (pthread_t*)malloc(n_threads*sizeof(pthread_t));
	while(n_threads--)
		pthread_create(&threads[n_threads], NULL, calculate_saxpy, args[n_threads]);
	return threads;
}

void wait_for_threads(pthread_t *threads, int n_threads){
	while(n_threads--){
		pthread_join(threads[n_threads], NULL);
	}
}