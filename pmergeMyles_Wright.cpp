/*  Myles Wright
    Butler University 2021
    CS452 PMERGE
*/
#include <iostream>
#include <stdio.h>
#include <string.h>
#include <ctime>
#include <math.h>
#include "mpi.h" // message passing interface
using namespace std;

// New compile and run commands for MPI!
// mpicxx -o blah file.cpp
// mpirun -q -np 32 blah

int bSearch(int * arr, int x, int p, int r){ //bSearch function
    if(p <= r){                                //x is what we are searching for //RANK is the bSearch array  //p is the partition                                    
        int m = (p+r)/2;    //index of middle element   //r is the index of the last element of the array
        if(arr[m] == x){
            return m;
        }
        if(arr[m] > x){
            return bSearch(arr, x, p, m-1);
        }
        if(arr[m] < x){
            return bSearch(arr, x, m+1, r);
        }
        if(x < arr[m] && x > arr[m-1]){
            //x is inbetween these two numbers but not in arr
            //the point is that the rank is still noted, it doesn't have to be inside of arr, we just need to know where it fits in.
            return m-1;
        }
    }
    return -1; //there was a major issue                      

}


void rank(int * a, int * RANK, int x, int i){ //rank function, uses binary search. i is the index value of x
    int rank = 0; //rank value
    int r = sizeof(a) - 1;
    rank = bSearch(a, x, 0, r); //using the array, what we are looking for, 0 as p, and the last index as r
    RANK[i] = rank;
}

void smerge(int * a, int first, int lasta, int lastb, int * output = NULL){ //smerge function
    
    int arrA = lasta - first + 1;
    int arrB = lastb - lasta; //size of b

    //temp arrays
    int * arrL = new int[arrA];
    int * arrR = new int[arrB];

    //copy
    for (int i = 0; i < arrA; i++){
        arrL[i] = a[first + i];
    }
    for (int j = 0; j < arrB; j++){
        arrR[j] = a[lasta + 1 + j];
    }

    int i = 0; //index of first array
    int j = 0; //index of second array
    int k = first; //index of merged arrays

    while(i < arrA && j < arrB){ //copy to arr
        if(arrL[i] <= arrR[j]){
            a[k] = arrL[i];
            i++;    
        }
        else{
            a[k] = arrR[j];
            j++;
        }
        k++;
    }

    while(i < arrA){ //remaining elements of arrL[]
        a[k] = arrL[i];
        i++;
        k++;
    }
    
    while(j < arrB){
        a[k] = arrR[j];
        j++;
        k++;
    }
}


void mergeSort(int * A){ //mergeSort function, basically the same as in smerge
   // mergeSort(A);
   // mergeSort(A);
   // smerge();
}

void printArray(int * a, int size){
    for(int i = 0; i < size; i++){
        cout << a[i] << " ";
    }
    cout << "\n";
}


int main (int argc, char * argv[]) {

	int my_rank;			// my CPU number for this process
	int p;					// number of CPUs that we have
	int source;				// bSearch of the sender
	int dest;				// bSearch of destination
	int tag = 0;			// message number
	char message[100];		// message itself
	MPI_Status status;		// return status for receive
	
	// Start MPI
	MPI_Init(&argc, &argv);
	
	// Find out my bSearch!
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	
	// Find out the number of processes!
	MPI_Comm_size(MPI_COMM_WORLD, &p);
	

    int * input = NULL;
    int * a = NULL;
    int * b = NULL;                //MYLES: Create global arrays
    int * RANKA = NULL;
    int * RANKB = NULL;
    int * WIN = NULL;
    int size = 0;

	// THE REAL PROGRAM IS HERE
    if(my_rank == 0){ //MYLES: create and broadcast the array
        cout << "What in gods name is the size of this thing?";

        cin >> size;
        input = new int[size];

        int size_div = size/2;
        int size_log = (size_div)/log2(size_div);

        a = &input[0];
        b = &input[size_div]; //MYLES: cut this hoe in half
        
        RANKA = new int[size_log]; //correct?
        RANKB = new int[size_log]; //please be correct
        

        //please sort a and b before broadcast
        

        /*MPI_Bcast(a, size/2, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(b, size/2, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(RANKA, (size/2)/log2(size/2), MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(RANKB, (size/2)/log2(size/2), MPI_INT, 0, MPI_COMM_WORLD);
        */

        
    }

    int * recv = new int[1]; //receiving buffer
    int s = 0;
    MPI_Scatter(a, size/p, MPI_INT, recv, 1, MPI_INT, 0, MPI_COMM_WORLD); //scatter the arrays
    MPI_Scatter(b, size/p, MPI_INT, recv, 1, MPI_INT, 0, MPI_COMM_WORLD);
    cout << "Process " << my_rank << "has data" << recv[0] << "\n"; //print test

	// Shut down MPI
	MPI_Finalize();

	return 0;
}