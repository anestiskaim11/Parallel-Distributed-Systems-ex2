/*
 *************************************************************************
 authors: Anestis Kaimakamidis 9627 Konstantinos Kalamaras 9716
 *************************************************************************
 MPI DISTRTIBUTE BY MEDIAN
 */
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <mpi.h>
#include <unistd.h>
#include "float.h"

//time values
struct timeval startwtime, endwtime;
static clock_t st_time;
static clock_t en_time;
double p_time;


//takes a 2d table and returns an 1d array
int* D2_to_D(int ** table, int rows, int cols){
    int *array = (int*)malloc(rows*cols*sizeof(int));
    for(int i = 0; i < rows; i++){
        for(int j = 0; j < cols; j++) array[i*cols + j] = table[i][j];
    }
    return array;
}

//takes an 1d array and returns a 2d table
int** D_to_2D(int *array, int rows, int cols){
    int **table = (int**)malloc(rows*sizeof(int*));
    for(int i = 0; i < rows; i++) table[i] = (int*)malloc(cols*sizeof(int));
    for(int i = 0; i < rows; i++){
        for(int j = 0; j < cols; j++) table[i][j] = array[i*cols + j];
    }
    return table;
}

/*
 function readArray(): reads a binary file of points and returns the data as an array
 */
double * readArray(){
    FILE *data;
    data=fopen("mnist.bin", "rb");
   
    long int *size = (long int *)malloc(2*sizeof(long int));
    fread(size,sizeof(long int),2,data);
    int rows = (int) size[0];
    int cols = (int) size[1];
    
    double *myArray=(double*)malloc(rows*cols*sizeof(double*));
    fread(myArray,sizeof(double),rows*cols,data);
    return myArray;
}
/*
 function find_final(int** final, int* results, int p, int n, int flag, int start): gets the results array and writes the exchange board (final)
 flag is 0 or 1 and determines if we are dealing with small or big distances
 for more information check the report
 */
void find_final(int** final, int* results, int p, int n, int flag, int start){
    int index = 0;
    int* given = (int*)malloc(p*sizeof(int));
    for(int i = 0; i < p; i++) given[i] = 0;
    if(flag){
        for(int i = 0; i < p; i++){
            for(int j = p - 1; j >= i + 1; j--){
                if(results[i] == n/p) continue;
                if(results[j] <= (n/p) - results[i] && results[j] > 0){
                    while(final[j][index + 2] != -1) index += 3;
                    final[j][index] = given[j];
                    final[j][1 + index] = results[j] + given[j] - 1;
                    final[j][2 + index] = i + start;
                    results[i] += results[j];
                    given[j] += results[j];
                    results[j] = 0;
                }
                else if(results[j] > (n/p) - results[i]){
                    while(final[j][2 + index] != -1) index += 3;
                    final[j][index] = given[j];
                    final[j][index + 1] = (n/p) - results[i] + given[j] - 1;
                    final[j][index + 2] = i + start;
                    given[j] += (n/p) - results[i];
                    results[j] -= (n/p) - results[i];
                    results[i] = n/p;
                }
                index = 0;
            }
        }
    }
    else{
        for(int i = p - 1; i > 0; i--){
            for(int j = 0; j < i; j++){
                if(results[i] == n/p) continue;
                if(results[j] <= (n/p) - results[i] && results[j] > 0){
                    while(final[j][index + 2] != -1) index += 3;
                    final[j][index] = (n/p) - results[j] - given[j];
                    final[j][1 + index] = (n/p) - given[j] - 1;
                    final[j][2 + index] = i + start;
                    results[i] += results[j];
                    given[j] += results[j];
                    results[j] = 0;
                }
                else if(results[j] > (n/p) - results[i]){
                    while(final[j][index + 2] != -1) index += 3;
                    final[j][index] = results[i] - given[j];
                    final[j][index + 1] = (n/p) - given[j] - 1;
                    final[j][index + 2] = i + start;
                    given[j] += (n/p) - results[i];
                    results[j] -= (n/p) - results[i];
                    results[i] = n/p;
                }
                index = 0;
            }
        }
    }
    free(given);
}


/*
 function swap(double *a, int s1, int s2, int d): swaps s1 and s2 d-dimentional points
 */
void swap(double *a, int s1, int s2, int d){
    double tmp;
    for(int i = 0; i < d; i++){
        tmp = a[s1*d + i];
        a[s1*d + i] = a[s2*d + i];
        a[s2*d + i] = tmp;
    }
}

/*
 function distribute(double *a, double *distances, int l, int r, double pivot, int d): distributes points with distance from pivot lower and higher than the median
 note that the median is referred as pivot in this function
 */
int distribute(double *a, double *distances, int l, int r, double pivot, int d){
    double tmp;
    while(l < r){
        while(distances[l] <= pivot) l++;
        while(distances[r] > pivot) r--;
        if(l < r){
            tmp = distances[l];
            distances[l] = distances[r];
            distances[r] = tmp;
            swap(a, l, r, d);
        }
        else if(l > r) return r;
    }
}


/*
 functions partition(double *arr, int l, int r) and quick_select(double *arr, int l, int r, int k): typical quickselect routine
 */
int partition(double *arr, int l, int r){
    double x = arr[r], tmp = 0;
    int i = l;
    for (int j = l; j <= r - 1; j++) {
        if (arr[j] <= x) {
            tmp = arr[i];
            arr[i] = arr[j];
            arr[j] = tmp;
            i++;
        }
    }
    tmp = arr[i];
    arr[i] = arr[r];
    arr[r] = tmp;
    return i;
}
double quick_select(double *arr, int l, int r, int k){
    if (k > 0 && k <= r - l + 1) {
 
        int index = partition(arr, l, r);
 
        if (index - l == k - 1)
            return arr[index];
 
        if (index - l > k - 1)
            return quick_select(arr, l, index - 1, k);
 
        return quick_select(arr, index + 1, r,
                            k - index + l - 1);
    }

}


/*
 function find_dist(double *a, double *distances, int size, double *c, int d): finds the distance between d-dimentional point c (pivot) and all the points in a array
 and stores the results in distances array
 */
void find_dist(double *a, double *distances, int size, double *c, int d){
    for(int i = 0; i < size; i++){
        for(int j = 0; j < d; j++){
            if(j == 0) distances[i] = 0;
            distances[i] += (a[i*d + j] - c[j]) * (a[i*d + j] - c[j]);
        }
    }
}

/*
 function distributeByMedian(int n, int p, int d, int SelfTID, double* a, double* pivot, int start):
 *n: total number of points
 *p: total number of processes
 *d: number of dimentions
 *SelfTID: rank of each process
 *a: data of each process
 *pivot: pivot point
 *start: start of a group of processes (needed for recursion)
 */
void distributeByMedian(int n, int p, int d, int SelfTID, double* a, double* pivot, int start){
    int t, err;
    MPI_Status mpistat;
    double *distances = (double*)malloc(n/p*sizeof(double));
    double median;
    if(SelfTID == start){ //master process
        find_dist(a, distances, n/p, pivot, d); //find distances
        //receive distances of other processes
        double **distances1 = (double**)malloc((p-1)*sizeof(double*));
        for(int i = 0; i < p - 1; i++) distances1[i] = (double*)malloc(n/p*sizeof(double));
        for(t= 1; t < p; t++){
            MPI_Recv(distances1[t-1],n/p,MPI_DOUBLE,t+start,MPI_ANY_TAG,MPI_COMM_WORLD,&mpistat);
        }
        double *array = (double*)malloc(n*sizeof(double));
        for(int i = 0; i < n; i++){
            if(i < n/p) array[i] = distances[i];
            else array[i] = distances1[i*p/n - 1][i % (n/p)];
        }
        free(distances1);
        //quickselect for median
        median = (quick_select(array, 0, n - 1, n/2) + quick_select(array, 0, n - 1, n/2 + 1))/2;
        //send median
        for(t = start + 1; t < p + start; t++){
            err = MPI_Send(&median,1,MPI_DOUBLE,t,100+t,MPI_COMM_WORLD);
            if( err ) {
                printf("Error=%i in MPI_Send to %i\n",err,t);
            }
        }
        free(array);
        //receive all the results (results: number of "small" distances results_max: number of "big distances")
        int *results = (int*)malloc(p*sizeof(int*));
        int *results_max = (int*)malloc(p*sizeof(int*));
        results[0] = distribute(a, distances, 0, n/p - 1, median, d);
        results[0]++;
        results_max[0] = (n/p) - results[0];
        for(t = 1; t < p; t++){
            MPI_Recv(&results[t],1,MPI_INT,t+start,MPI_ANY_TAG,MPI_COMM_WORLD,&mpistat);
            results[t]++;
            results_max[t] = (n/p) - results[t];
        }
        //find and send the final board
        int **final = (int**)malloc(p*sizeof(int*));
        for(int i = 0; i < p; i++) final[i] = (int*)malloc(3*p/2*sizeof(int));
        for(int i = 0; i < p; i++){
            for(int j = 2;  j < 3*p/2; j += 3) final[i][j] = -1;
        }
        find_final(final, results, p, n, 1, start);
        find_final(final, results_max, p, n, 0, start);
        int* final_1d = D2_to_D(final, p, 3*p/2);
        for(t = start + 1; t < p + start; t++){
            MPI_Send(final_1d,3*p*p/2,MPI_INT,t,100+t,MPI_COMM_WORLD);
        }
        free(final_1d);
        free(results);
        free(results_max);
        //make exchanges
        double *received = (double*)malloc(d*n/p*sizeof(double));
        int index = 0;
        for(int i = 0; i < p; i++){
            for(int j = 2; j < 3*p/2; j += 3){
                int count = final[i][j-1] - final[i][j-2] + 1;
                if(i + start == SelfTID && final[i][j] != -1){
                    MPI_Send(&(a[d*final[i][j-2]]),d*count,MPI_DOUBLE,final[i][j],200+SelfTID,MPI_COMM_WORLD);
                }
                if(final[i][j] == SelfTID){
                    MPI_Recv(&(received[index]),d*count,MPI_DOUBLE,i+start,MPI_ANY_TAG,MPI_COMM_WORLD,&mpistat);
                    index += (final[i][j-1] - final[i][j-2] + 1)*d;
                }
            }
        }
        //store received values
        for(int i = n*d/p - index; i < n*d/p; i++) a[i] = received[i - n*d/p + index];
        free(received);
        free(final);
        free(distances);
    }
    else{ //other processes
        find_dist(a, distances, n/p, pivot, d); //find distances
        MPI_Send(distances,n/p,MPI_DOUBLE,start,100+SelfTID,MPI_COMM_WORLD); //send distances
        MPI_Recv(&median,1,MPI_DOUBLE,start,MPI_ANY_TAG,MPI_COMM_WORLD,&mpistat); //get median
        int results = distribute(a, distances, 0, n/p - 1, median, d);//calculate results
        MPI_Send(&results,1,MPI_INT,start,100+SelfTID,MPI_COMM_WORLD);//send results
        int *array = (int *)malloc(3*p*p/2*sizeof(int));
        MPI_Recv(array,3*p*p/2,MPI_INT,start,MPI_ANY_TAG,MPI_COMM_WORLD,&mpistat);//receive exchange table
        int** table = D_to_2D(array, p, 3*p/2);
        free(array);
        double *received = (double*)malloc(d*n/p*sizeof(double));
        //make exchanges
        int index = 0;
        for(int i = 0; i < p; i++){
            for(int j = 2; j < 3*p/2; j += 3){
                int count = table[i][j-1] - table[i][j-2] + 1;
                if(i + start == SelfTID && table[i][j] != -1) MPI_Send(&(a[d*table[i][j-2]]),d*count,MPI_DOUBLE,table[i][j],200+SelfTID,MPI_COMM_WORLD);
                if(table[i][j] == SelfTID){
                    MPI_Recv(&(received[index]),d*count,MPI_DOUBLE,i+start,MPI_ANY_TAG,MPI_COMM_WORLD,&mpistat);
                    index += (table[i][j-1] - table[i][j-2] + 1)*d;
                }
            }
        }
        //store results
        if(SelfTID - start < p/2){
            for(int i = n*d/p - index; i < n*d/p; i++) a[i] = received[i - n*d/p + index];
        }
        else if(SelfTID - start >= p/2){
            for(int i = 0; i < index; i++) a[i] = received[i];
        }
        free(received);
        free(table);
        free(distances);
    }
}


/*
 function control_recursion(int start, int n, int p, int d, int SelfTID, double* a, double *pivot):
 *start: start of a group of processes
 *n: total points in a group
 *p: total processes in a group
 *d: number of dimentions
 *SelfTID: rank of the process
 *a: data of each process
 *pivot: pivot point
 */
void control_recursion(int start, int n, int p, int d, int SelfTID, double* a, double *pivot){
    if(p == 1) return;
    distributeByMedian(n, p, d, SelfTID, a, pivot, start);
    //split into two groups
    if(SelfTID - start >= p/2) start += p/2;
    control_recursion(start, n/2, p/2, d, SelfTID, a, pivot);
}


int main(int argc, char** argv){
    int n = 65536, d = 784, p;
    int SelfTID, NumTasks, err, t;
    err = MPI_Init( &argc, &argv );
    MPI_Status mpistat;
    if( err ) {
        printf("Error=%i in MPI_Init\n",err);
    }
    MPI_Comm_size( MPI_COMM_WORLD, &NumTasks );
    MPI_Comm_rank( MPI_COMM_WORLD, &SelfTID );
    p = NumTasks;
    double *a = (double*)malloc(n/p*d*sizeof(double));
    double *pivot = (double*)malloc(d*sizeof(double));
    if(SelfTID == 0){
        gettimeofday (&startwtime, NULL);
        double *b = readArray(); //read data
        for(int i = 0; i < n*d/p; i++) a[i] = b[i];
        //send data
        for(t=1; t<p; t++){
            err = MPI_Send(&(b[t*d*n/p]),n*d/p,MPI_DOUBLE,t,100+t,MPI_COMM_WORLD);
            if( err ) {
                printf("Error=%i in MPI_Send to %i\n",err,t);
            }
        }
        //find and send random pivot
        int random = rand() % n/p;
        for(int i = random*d; i < random*d+d; i++) pivot[i - random*d] = a[i];
        for(t=1; t<p; t++){
            err = MPI_Send(pivot,d,MPI_DOUBLE,t,100+t,MPI_COMM_WORLD);
            if( err ) {
                printf("Error=%i in MPI_Send to %i\n",err,t);
            }
        }
        free(b);
    }
    else{
        MPI_Recv(a,n*d/p,MPI_DOUBLE,0,MPI_ANY_TAG,MPI_COMM_WORLD,&mpistat);
        MPI_Recv(pivot,d,MPI_DOUBLE,0,MPI_ANY_TAG,MPI_COMM_WORLD,&mpistat);
    }
    //controll recursion
    control_recursion(0, n, p, d, SelfTID, a, pivot);
    //find new distances
    double* distances = (double*)malloc(n/p*sizeof(double));
    find_dist(a, distances, n/p, pivot, d);
    //check the result
    double max = DBL_MIN, min = DBL_MAX;
    for(int i = 0; i < n/p; i++){
        if(distances[i] > max) max = distances[i];
        if(distances[i] < min) min = distances[i];
    }
    if(!SelfTID){
        double** minmax = (double**)malloc(p*sizeof(double*));
        for(int i = 0; i < p; i++) minmax[i] = (double*)malloc(2*sizeof(double));
        minmax[0][0] = min; minmax[0][1] = max;
        for(t = 1; t < p; t++){
            MPI_Recv(&(minmax[t][0]),1,MPI_DOUBLE,t,MPI_ANY_TAG,MPI_COMM_WORLD,&mpistat);
            MPI_Recv(&(minmax[t][1]),1,MPI_DOUBLE,t,MPI_ANY_TAG,MPI_COMM_WORLD,&mpistat);
        }
        int ok = 1;
        for(int i = 1; i < p; i++){
            if(minmax[i][0] < minmax[i-1][1]) ok = 0;
        }
        if(ok) printf("CHECK OK\n");
        else printf("WRONG RESULT\n");
        //print runtime
        gettimeofday (&endwtime, NULL);
        p_time = (double)((endwtime.tv_usec - startwtime.tv_usec)/1.0e6
                              + endwtime.tv_sec - startwtime.tv_sec);
        printf("TIME:%f\n", p_time);
    }
    else{
        MPI_Send(&min,1,MPI_DOUBLE,0,100+SelfTID,MPI_COMM_WORLD);
        MPI_Send(&max,1,MPI_DOUBLE,0,100+SelfTID,MPI_COMM_WORLD);
    }
    MPI_Finalize();
    free(distances);
    free(a);
    free(pivot);
}


