/*  Myles Wright
    Butler University 2021
    CS452 PMERGE
*/
#include <iostream>
#include <stdio.h>
#include <string.h>
#include <ctime>
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
            return bSearch(arr, p, m-1, num);
        }
        if(arr[m] < x){
            return bSearch(arr, m+1, r, num);
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


void mergeSort(int * A, ){ //mergeSort function, basically the same as in smerge
    mergeSort(A, );
    mergeSort(A, );
    smerge();


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
	
	// THE REAL PROGRAM IS HERE
    if(my_rank == 0){ //create and broadcast the array, array creation taken from mergesortMylesWright.cpp
        cout << "How large is the first array hoe?\n"; 
        int size = 0; //size of array
        cin >> size; //get input
        int * a = new (nothrow) int [size]; //assign size in array a decclaration
                                             //nothrow will not throw if declaration fails :)
        int * RANKA = new int[size];
        cout << "How large is the second array?";
        int size_2 = 0;
        cin >> size_2;
        int * b = new int [size_2];
        int * RANKB = new int[size_2];

        if(a==nullptr){ //if there is a memory allocation failure, deal with it
            cout << "OH NO ABORT ABORT MISSON WE ARE GOING DOWN HOLY SHIT I HAVE TWO KIDS IT CANT END LIKE-";
            exit(EXIT_FAILURE);
        }

        srand(71911); //random number gen
        int x = 1000000; //prompt said to use this number

        for(int i = 0; i < size; i++){
            a[i] = rand() % 8; // :)
        }
        srand(time());
        for(int i = 0; i < size; i++){
            b[i] = (rand() % 500)+1; // :) between 0 and 500
        }
        
        MPI_Bcast(b, size_2, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(a, size, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(SRANKA, size, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(SRANKB, size_2, MPI_INT, 0, MPI_COMM_WORLD);
        int i = 1;
        MPI_Send( i,  1, MPI_INT, 1, 0, MPI_COMM_WORLD);
        

    }
    else{
        MPI_Recv(i, 1, MPI_INT, (my_rank-1), 0, MPI_COMM_WORLD);

        while(i<p){
            if(i==my_rank){
                //recieve all of the goods ;)
                MPI_Recv(a, size, 0, MPI_INT, 0, MPI_COMM_WORLD); //the size variable probably wont work
                MPI_Recv(b, size_2, 0, MPI_INT, 0, MPI_COMM_WORLD);
                MPI_Recv(SRANKA, size, 0, MPI_INT, 0, MPI_COMM_WORLD);
                MPI_Recv(SRANKB, size_2, 0, MPI_INT, 0, MPI_COMM_WORLD);

                //send to rank function
                //give element at logi

                //rank of a in b
                //define the element that were looking for
                int log = log2(i);
                rank(b, SRANKA, a[log], log); //pass in the array we are looking through, the rankarrayA, the actual value we are looking for(every lognth element cause why would we want to do every fucking one), and the index value of log(i)
                rank(a, SRANKB, b[log], log); //rank of b in a

                MPI_Send{}
                

            }
        }
    }
    /*
    FIRST: Gotta sort both arrays.
    SECOND: MERGE THE FUCKERS IN PARALLEL
    THIRD: ???
    FOURTH: PROFIT

    */
	// Shut down MPI
	MPI_Finalize();

	return 0;
}