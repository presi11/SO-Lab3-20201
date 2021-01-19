#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <assert.h>
#include <pthread.h>
#include <semaphore.h>
#include <sys/time.h>

typedef struct Args{
	// Argumentos que se les pasa a cada hilo, por que las funciones para crear un hilo solo acepta un argumento y ademas es apuntador
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

sem_t sem1;
sem_t sem2;
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
	sem_init(&sem1, 0, value);
	sem_init(&sem2, 0, 1);
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
	it = 0; //Numero de iteracciones que llevamos hasta el momento, inicializadas en 0
	int *n_threads_finished = (int*) calloc(1, sizeof(int)); //Apuntador al entero de numero de hijos que termino de hacer el trabajo
	Saxpy_args **args = create_args(p, n_threads, max_iters, n_threads_finished, &it, a,  X, Y, Y_avgs); //Llamadoa la creacion de argumentos
	pthread_t *threads;// Declaracion de arreglo de hilos
	threads = create_thread_arr(n_threads, args); // n_threads crecion del nuemro de hilos ingresado
	wait_for_threads(threads, n_threads); // Funcion para esperar que todos los hilos acaben
	gettimeofday(&t_end, NULL);// Toma de tiempo final
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

//Creacion de un arreglo de argumentos, apuntadores a la estructura de los argumentos 
Saxpy_args **create_args(int n_elements, int n_threads, int max_iters, int *n_threads_finished, int *n_iter,
double scalar,  double *X, double *Y, double *Y_avgs){
	Saxpy_args **args = (Saxpy_args**)malloc(n_threads*sizeof(Saxpy_args*));//Grabamos en memoria los apuntadores
	int arr_chunk = n_elements/n_threads;//Dividimos el arreglo en el numero de hilos para asignar la parte del arreglo que trabajara
	int n_leftovers = n_elements%n_threads;//Elemntos que sobran de uan divison incompleta, se reparten equitatuivamente entre los hilos existentes
	int curr_index = 0; //Apuntamos a la posicion inicial del arreglo
	int end_index = 0; // Declaramos la posicion final, la cual se actulizara posteriormente
	int n_threads_copy = n_threads;//Copia de numero de hilos, que se le da cada hilo, para saber cuantos hilos existen en total
	int n; //Cuantos elementos existene ne cada reparticion de elemntos que trabajara un hilo
	//REpartimos las partes del arreglo para cada hilo
	while(n_threads--){
		args[n_threads] = (Saxpy_args*)malloc(sizeof(Saxpy_args));//Declaramos el apuntador al argumento
		//Verificamos si hay sobras de la divion para repartirlas
		if(n_leftovers){
			n = arr_chunk + 1;
			n_leftovers--;
		}else //No se ele asigna trabajo extra a los hilos
			n = arr_chunk;
		args[n_threads] -> curr_index = curr_index; //-> Acceder al elemnto dentro del apuntador
		end_index = curr_index + n - 1; //Se actuliza el final del arreglo
		args[n_threads] -> end_index = end_index; //Se asigna el final del arreglo para dicho hilo
		args[n_threads] -> n_threads = n_threads_copy;// Le asignamos el numero de hilos que existen
		args[n_threads] -> max_iters = max_iters;//Maximo de interaciones totales que hara el hilo
		args[n_threads] -> n_threads_finished = n_threads_finished;//Saber cuantos hilos a termiando, se comparte la direcciond memoria, sin lo cambio en uno todos lo saben
		args[n_threads] -> n_iter = n_iter;//Numero de iteraciones hasta el momento, todos ceros en su creacion
		args[n_threads] -> scalar = scalar;//Scalar sera el valor de A*x[i] de la funcion
		args[n_threads] -> X = X;//El apuntador de X de la funcion
		args[n_threads] -> Y = Y;//El apuntador de Y de la funcion
		args[n_threads] -> Y_avgs = Y_avgs;//Promedios de valores de Y que tendra cada hilo
		curr_index = end_index + 1;//Actulizamos para que reparta el proximo grupo del arreglo para el siguiente hilo
	}
	return args;//Retornamos los argumentos
}

//Calculo del saxpy, Void sirve para cualqueir tipo de apuntador
void *calculate_saxpy(void *args){
	Saxpy_args *saxpy_args = (Saxpy_args*)args;//casting, queremos saber que si esta a la estructura creda anteriormente
	int curr_index = saxpy_args -> curr_index;//Sacamos los argumentos para manejarlos de forma local
	int end_index = saxpy_args -> end_index;//Sacamos los argumentos para manejarlos de forma local
	int n_threads = saxpy_args-> n_threads;//Sacamos los argumentos para manejarlos de forma local
	int max_iters = saxpy_args -> max_iters;//Sacamos los argumentos para manejarlos de forma local
	int *n_threads_finished =  saxpy_args -> n_threads_finished;//Sacamos los argumentos para manejarlos de forma local
	int *n_iter = saxpy_args -> n_iter;//Sacamos los argumentos para manejarlos de forma local
	double scalar = saxpy_args -> scalar;//Sacamos los argumentos para manejarlos de forma local
	double *X = saxpy_args -> X;//Sacamos los argumentos para manejarlos de forma local
	double *Y = saxpy_args -> Y;//Sacamos los argumentos para manejarlos de forma local
	double *Y_avgs = saxpy_args -> Y_avgs;//Sacamos los argumentos para manejarlos de forma local
	double result;//Declaramos el resultado
	//Iteraremos hasta que le maximo de iteraciones de dicho hilo se cumpla
	while(*n_iter < max_iters){
		result = calculate_saxpy_n_time(curr_index, end_index, *n_iter+1, scalar, X, Y)/p; //CAlculamos ell Y final de cada hilo
		pthread_mutex_lock(&lock);//Cerrojo, bloqueo las cosas criticas, para que no sean modificados por los demas hilos
		Y_avgs[*n_iter] += result;//El trabajo que realizo cada hilo
		(*n_threads_finished)++;//Este hilo ya termino su trabajo
		//EL ultimo hilo en terminar reinica el segundo semaforo y aumenta el numero de iteracciones hechas
		//Segundo semaforo lo utilizmaos para reinicia rlas condiciones inciales y poder realizar nuevamente las iteracciones
		//El primer semaforo lo usamos para saber que todos los hilos terminaron
		if(*n_threads_finished == n_threads){
			(*n_iter)++;//aumentamos el numero de iteraciones que haramos el calculo de SAXPY
			sem_init(&sem2, 0, 0);//Se inicliaza asi para que todos los hilos tengan qu esperar
			sem_post(&sem1);//Activa un hilo que este en reposo
		}
		pthread_mutex_unlock(&lock);//Desbloquemos los valores para ser utilziados
		sem_wait(&sem1);//Ponemos a esperar a los hilos
		sem_post(&sem1);//Activamos un hilo y se hace un efecto cascada  a traves de los hilos activandoce
		pthread_mutex_lock(&lock);//reiniciar
		(*n_threads_finished)--;//Reiniciar uno a uno los hilos para decirles que no han terminado y ponerlos a esperar hasta que ninguna haya acabado
		if(*n_threads_finished == 0){
			sem_init(&sem1, 0, 0);//reinicmos el primer semaforo  
			sem_post(&sem2);//Activo un hilo en reposo
		}
		pthread_mutex_unlock(&lock);// Desbloqueamos Nuevamente los valores para ser utilizados 
		sem_wait(&sem2);//Ponemos a esperar al hilo
		sem_post(&sem2);//Activamos un hilo y se hace un efecto cascada  a traves de los hilos activandoce		
	}
	//Terminacion de la iteracion de un hilo

	update_Y_n_time(curr_index, end_index, max_iters, scalar, X, Y);//ACtulizamos el vector Y con sus valores resultantes depues de las n iteracione de un hilo
	pthread_exit(NULL);//Acabo 
}

//Creamos el arreglo de apuntadores
pthread_t *create_thread_arr(int n_threads, Saxpy_args **args){
	pthread_t *threads = (pthread_t*)malloc(n_threads*sizeof(pthread_t));//Arreglo de hilos instanciado
	while(n_threads--)
		pthread_create(&threads[n_threads], NULL, calculate_saxpy, args[n_threads]);// Funcion para crear un hilo
	return threads;
}

//Ponemos a esperar los hilos para que Actualicen el valor de Y
void wait_for_threads(pthread_t *threads, int n_threads){
	while(n_threads--){
		pthread_join(threads[n_threads], NULL);
	}
}

//Modificiaciond del arreglo final tras haber relizado todas las iteracione por hilo
void update_Y_n_time(int curr_index, int end_index, int n_time, double scalar, double *X, double *Y){
	while(curr_index <= end_index){
		Y[curr_index] = scalar*n_time*X[curr_index] + Y[curr_index];
		curr_index++;
	}
}

//Calculamos el valor de Y final actual de dicho hilo
double calculate_saxpy_n_time(int curr_index, int end_index, int n_time, double scalar, double *X, double *Y){
	double result = 0.0; // inicliazamos el resultado en 0.0
	//Recorremso el grupo del arreglo asigando para un hilo en particular
	while(curr_index <= end_index){
		//Asigamos el valor que tendra final despues de recorrer todo el grupo del hilo
		result += scalar*n_time*X[curr_index] + Y[curr_index];//n_time sera el numero de iteraciones en ese hilo
		curr_index ++; //Actulizamos el index de l while
	}
	return result;//retorno
}