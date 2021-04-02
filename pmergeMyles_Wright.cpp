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

int srank(int * arr, int x, int p, int r){ //bSearch function
    if(p <= r){                                //x is what we are searching for //RANK is the bSearch array  //p is the partition                                    
        int m = (p+r)/2;    //index of middle element   //r is the index of the last element of the array
        if(arr[m] == x){
            return m;
        }
        else if(arr[m] > x){
            return srank(arr, x, p, m-1);
        }
        else if(arr[m] < x){
            return srank(arr, x, m+1, r);
        }
        
    }
    else{
        return p;
    }
    return -1; //there was a major issue                      

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


void mergesort(int * a, int first, int last){ //mergesort. yee yee
    if(first >= last){
        return;
    }
    int p = first + (last-first)/2;
    mergesort(a, first, p);
    mergesort(a, p+1, last);
    smerge(a, first, p, last);
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
	

    int * input; int * a; int * b; 
    int RANKA;
    int RANKB; 
    int * WIN; int * recvA; int * recvB;              //MYLES: Create global arrays
    
    int size;
    

	// THE REAL PROGRAM IS HERE
    if(my_rank == 0){ //MYLES: create and broadcast the array
        cout << "What in gods name is the size of this thing?\n";
        cin >> size;
        input = new int[size];
        for(int i = 0; i < size; i++){
            input[i] = rand() %(500-0+1)+0; // :)
        }
        int size_div = size/2;
        a = &input[0];
        b = &input[size_div]; //MYLES: cut this hoe in half
    }
    
    MPI_Bcast(&size, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    
    int recv_size = (size/2)/p;
    recvA = new int[recv_size];
    recvB = new int[recv_size];
    MPI_Scatter(a, recv_size, MPI_INT, recvA, recv_size, MPI_INT, 0, MPI_COMM_WORLD); //scatter the arrays
    MPI_Scatter(b, recv_size, MPI_INT, recvB, recv_size, MPI_INT, 0, MPI_COMM_WORLD);
    
    mergesort(recvA, 0, (recv_size-1));
    mergesort(recvB, 0, (recv_size-1));

    cout << "Process " << my_rank << " Array A: ";
    printArray(recvA, recv_size); //print test
    
    cout << "Process " << my_rank << " Array B: ";
    printArray(recvB, recv_size); //print test

    RANKA = 0;
    RANKB = 0;
    RANKA = srank(b, recvA[0], 0, (recv_size-1));
    RANKB = srank(a, recvB[0], 0, (recv_size-1));

    cout << "Process " << my_rank << " RANKA: " << RANKA;
    cout << "Process " << my_rank << " RANKB: " << RANKB;

	// Shut down MPI
	MPI_Finalize();

	return 0;
}