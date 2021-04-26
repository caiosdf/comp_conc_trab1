/**
 * 
 * Cálculo da propagação de calor em uma placa (reprsentada por um matriz)
 * O calor se propaga a partir da fonte que localizada na extremidade
 * 
 * Uso do programa:
 *      ./trabalho1 <dimensao da matriz> <numero de threads>
 * 
 */

#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include "timer.h"

#define MAX_ITERATIONS 3000

float *M1, *M2,*M1Seq, *M2Seq, *lastItr;
pthread_mutex_t barrier;
pthread_cond_t go;
int narrived = 0, nthreads, sequencial_executed = 0; 
float epsilon = 0.001f;

typedef struct
{
    int id;
    int dim;
} argType;
//função utilizada para comparar floats com a precisão de 10
int compare_float(float x, float y, float epsilon){
    if(x > y){
        if((x - y) < epsilon){
            return 1;
        }
        else return 0;
    }
    else{
        if((y - x) < epsilon){
            return 1;
        }
        else return 0;
    }
}
void Barrier(){
    pthread_mutex_lock(&barrier);
    narrived++;
    if(narrived == nthreads){
        narrived = 0;
        pthread_cond_broadcast(&go);
    }else {
        pthread_cond_wait(&go, &barrier);
    }
    pthread_mutex_unlock(&barrier);
}
void *task(void *arg)
{
    argType *args = (argType *)arg;
    float *old, *new;
    int itr;
    for (itr = 1; itr <= MAX_ITERATIONS; itr++)
    {
        if (itr % 2 == 0)
        {
            old = M1;
            new = M2;
        }
        else
        {
            old = M2;
            new = M1;
        }
        for (int row = args->id + 1; row < args->dim - 1; row += nthreads)
        {
            for (int col = 1; col < args->dim - 1; col++)
            {
                new[row * (args->dim) + col] = 0.25 * old[((row - 1) * (args->dim) + col)];
                new[row * (args->dim) + col] += 0.25 * old[((row + 1) * (args->dim) + col)];
                new[row * (args->dim) + col] += 0.25 * old[(row * (args->dim) + (col - 1))];
                new[row * (args->dim) + col] += 0.25 * old[(row * (args->dim) + (col + 1))];
            }
        }
        Barrier();
    }
    pthread_exit(NULL);
}

void initiator(int dim)
{
    int i, j;

    /* Preenche toda a matriz com zeros */
    for (i = 0; i < dim; i++)
    {
        for (j = 0; j < dim; j++)
        {
            M1[i * dim + j] = 0.0;
            M2[i * dim + j] = 0.0;
            M1Seq[i * dim + j] = 0.0;
            M2Seq[i * dim + j] = 0.0;
            lastItr[i * dim + j] = 0.0;
        }
    }

    /* Preenche a equerda da matriz com 0 e a borda direita com um degradê de 0 a 100 */
    for (i = 0; i < dim; i++)
    {
        M1[i * dim + 0] = 0.0;
        M1[i * dim + (dim - 1)] = (100.0 / (dim - 1)) * i;
        M2[i * dim + 0] = 0.0;
        M2[i * dim + (dim - 1)] = (100.0 / (dim - 1)) * i;
        M1Seq[i * dim + 0] = 0.0;
        M1Seq[i * dim + (dim - 1)] = (100.0 / (dim - 1)) * i;
        M2Seq[i * dim + 0] = 0.0;
        M2Seq[i * dim + (dim - 1)] = (100.0 / (dim - 1)) * i;
        lastItr[i * dim + 0] = 0.0;
        lastItr[i * dim + (dim - 1)] = (100.0 / (dim - 1)) * i;
    }

    /* Preenche a borda superios com 0 e a borda inferior com um degradê de 0 a 100 */
    for (j = 0; j < dim; j++)
    {
        M1[0 * dim + j] = 0.0;
        M1[(dim - 1) * dim + j] = (100.0 / (dim - 1)) * j;
        M2[0 * dim + j] = 0.0;
        M2[(dim - 1) * dim + j] = (100.0 / (dim - 1)) * j;
        M1Seq[0 * dim + j] = 0.0;
        M1Seq[(dim - 1) * dim + j] = (100.0 / (dim - 1)) * j;
        M2Seq[0 * dim + j] = 0.0;
        M2Seq[(dim - 1) * dim + j] = (100.0 / (dim - 1)) * j;
        lastItr[0 * dim + j] = 0.0;
        lastItr[(dim - 1) * dim + j] = (100.0 / (dim - 1)) * j;
    }
}

void print_matrix(float *M, int dim)
{
    for (int row = 0; row < dim; row++)
    {
        for (int col = 0; col < dim; col++)
        {
            printf("%.2f\t", M[row * dim + col]);
        }
        printf("\n");
    }
    printf("\n");
}

void destructor(){
    free(M1);
    free(M2);
    free(M1Seq);
    free(M2Seq);
    free(lastItr);
}

int task_sequential(int dim){
    float *oldSeq, *newSeq;
    int i;
    if (M1Seq == NULL || M2Seq == NULL)
    {
        printf("ERRO--malloc\n");
        return 2;
    }
    for (i = 1; i <= MAX_ITERATIONS; i++)
    {
        if (i % 2 == 0)
        {
            oldSeq = M1Seq;
            newSeq = M2Seq;
        }
        else
        {
            oldSeq = M2Seq;
            newSeq = M1Seq;
        }

        if(i == MAX_ITERATIONS - 1){
            for (int row = 1; row < dim - 1; row ++)
            {
                for (int col = 1; col < dim - 1; col++)
                {
                    lastItr[row * dim + col] = oldSeq[row * dim + col];
                    newSeq[row * dim + col] = 0.25 * oldSeq[((row - 1) * dim + col)];
                    newSeq[row * dim + col] += 0.25 * oldSeq[((row + 1) * dim + col)];
                    newSeq[row * dim + col] += 0.25 * oldSeq[(row * dim + (col - 1))];
                    newSeq[row * dim + col] += 0.25 * oldSeq[(row * dim + (col + 1))];
                }
            }
            continue;
        }

        for (int row = 1; row < dim - 1; row ++)
        {
            for (int col = 1; col < dim - 1; col++)
            {
                newSeq[row * dim + col] = 0.25 * oldSeq[((row - 1) * dim + col)];
                newSeq[row * dim + col] += 0.25 * oldSeq[((row + 1) * dim + col)];
                newSeq[row * dim + col] += 0.25 * oldSeq[(row * dim + (col - 1))];
                newSeq[row * dim + col] += 0.25 * oldSeq[(row * dim + (col + 1))];
            }
        }
        
    }
    sequencial_executed = 1;
    return 0;
}

//Verifica se os valores da execução sequencial e paralela são iguais.
int consistency_test(float *A, float *B, int dim){
    for(int i = 0; i < dim; i++){
        for(int j = 0; j < dim; j++){
            if(A[i*dim + j] != B[i*dim + j]){
                printf("ERRO: inconsistência nos resultados\n");
                destructor();
                return 4;
            }
        }
    }
    return 0;
}

int coherence_test(float *A, int dim){
    if(!sequencial_executed) task_sequential(dim);
    for(int i = 1; i < dim - 1; i++){
        for(int j = 1; j < dim - 1; j++){
            float expected_value =  (0.25 * lastItr[((i - 1) * dim + j)]) + (0.25 * lastItr[((i + 1) * dim + j)]) + (0.25 * lastItr[(i * dim + (j - 1))]) + (0.25 * lastItr[(i * dim + (j + 1))]);
            if(!compare_float(A[i*dim+j], expected_value, epsilon)){
                printf("%d\n", i*dim + j);
                printf("%.10f\n", A[i*dim + j]);
                printf("%.10f\n", expected_value);
                printf("ERRO: incoerência nos resultados\n");
                destructor();
                return 5;
            }
            expected_value = 0;
        }
    }
    return 0;
}

int main(int argc, char **argv)
{
    int dim;
    pthread_t *tid;
    argType *args;
	double start, end, delta;

    if (argc < 3)
    {
        printf("Digite: %s <dimensao da matriz> <numero de threads>\n", argv[1]);
    }

    dim = atoi(argv[1]);
    nthreads = atoi(argv[2]);

    if (nthreads > dim - 2)
    {
        nthreads = dim - 2;
    }

    if (dim == 0 || nthreads == 0)
    {
        printf("Digite: %s <dimensao da matriz> <numero de threads>\n", argv[1]);
        return 1;
    }

    M1 = (float *)malloc(sizeof(float) * dim * dim);
    M2 = (float *)malloc(sizeof(float) * dim * dim);
    M1Seq = (float *)malloc(sizeof(float) * dim * dim);
    M2Seq = (float *)malloc(sizeof(float) * dim * dim);
    lastItr = (float *)malloc(sizeof(float) * dim * dim);
    
    if (M1 == NULL || M2 == NULL || M1Seq == NULL || M2Seq == NULL || lastItr == NULL)
    {
        printf("ERRO--malloc\n");
        return 2;
    }

    initiator(dim);

    tid = (pthread_t *)malloc(sizeof(pthread_t) * nthreads);
    if (tid == NULL)
    {
        puts("ERRO--malloc");
        return 2;
    }
    args = (argType *)malloc(sizeof(argType) * nthreads);
    if (args == NULL)
    {
        puts("ERRO--malloc");
        return 2;
    }

    GET_TIME(start);
    for (int i = 0; i < nthreads; i++)
    {
        (args + i)->id = i;
        (args + i)->dim = dim;
        if (pthread_create(tid + i, NULL, task, (void *)(args + i)))
        {
            puts("ERRO--pthread_create");
            return 3;
        }
    }

    for (int i = 0; i < nthreads; i++)
    {
        pthread_join(*(tid + i), NULL);
    }
	GET_TIME(end)
	delta = end - start;
    printf("Tempo de execução concorrente: %lf segundos.\n", delta);
    
    GET_TIME(start);
    task_sequential(dim);
	GET_TIME(end)
	delta = end - start;
    printf("Tempo de execução sequencial: %lf segundos.\n", delta);
    
    free(tid);
    free(args);
    int test1 = consistency_test(M1,M1Seq, dim);
    int test2 = coherence_test(M1, dim);
    if (test1 != 0) return test1;
    if (test2 != 0) return test2;
    //print_matrix(M1, dim); //printa a aproximação calculada.
    //print_matrix(lastItr, dim);
    destructor();

    return 0;
}
