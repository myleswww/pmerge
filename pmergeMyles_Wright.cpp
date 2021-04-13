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

int srank(int * arr, int x, int valtofind){ //MYLES: a fun binary search style algorithm that returns the value of where the value should be
    if (x==1){
        if(valtofind<arr[0]){
            return 0;
        }
        else{
            return 1;
        }
    }
    else{
        if(valtofind < arr[x/2]){
            return srank(arr, x/2, valtofind);
        }
        else{
            return x/2 + srank(&arr[x/2], x/2, valtofind);
        }
    }

}

void insertionSort(int* array, int size) //CATHAL: hehe
{
    int k;
    int j;
    for (int i = 1; i < size; i++)
    {
        k = array[i];
        j = i;
        while (j > 0 && array[j - 1] > k)
        {
            array[j] = array[j - 1];
            j--;
        }
        array[j] = k;
    }
}

void smerge(int* A, int a_start, int a_end, int* B, int b_start, int b_end, int* C, int c_start, int c_end)
{
	for(int i = c_start; i < c_end; i++)
	{
		//If there are elements left in A and either there are no elements left in B, 
		//or the next A element is less than the next B element, place the next A element
		//In the next C spot
		if(a_start < a_end && (b_start == b_end || A[a_start] < B[b_start])){
			C[i] = A[a_start++];
		}
		else{		//Else place the next B element in the next C spot
			C[i] = B[b_start++];
		}
	}
}

void printArray(int * a, int size){
    for(int i = 0; i < size; i++){
        cout << a[i] << " ";
    }
    cout << "\n";
}

void clear(int*a, int n){
    for(int i = 0; i<n; i++){
        a[i] = 0;
    }
}


void pmerge(int* A, int* B, int* C, int n, int my_rank, int p)		//Merge Array A and B, both size n/2, into array C, size n
{
	int logn = log2(n/2);
	int x = ceil((n/2)/logn);		//Number of sampled elements
	
	//Arrays to hold positions and ranks used to determine end points of sub arrays to merge using sequential merge
	//lower case is used for local arrays to be filled out by individual processes and then all reduced into capital arrays
	//i.e Allreduce(ar, AR, n, ...)
	int* ar = new int[2*x + 2];
	int* br = new int[2*x + 2];
	int* Ar = new int[2*x + 2];
	int* Br = new int[2*x + 2];
	clear(ar,2*x+2);
	clear(br,2*x+2);
	clear(Ar,2*x+2);
	clear(Br,2*x+2);
	
	for(int i = my_rank; i < x; i+=p){		//Fill in positions of sampled elements(in parallel)
		ar[i] = 1 + i*logn;
		br[i] = 1 + i*logn;
	}		
	for(int i = my_rank; i < x; i+=p){		//Fill in ranks of sampled elements(in parallel)
		br[i+x] = srank(B, n/2, A[0+i*logn]);
		ar[i+x] = srank(A, n/2, B[0+i*logn]);
	}
	MPI_Allreduce(ar, Ar, 2*x+2, MPI_INT, MPI_SUM, MPI_COMM_WORLD);		//combine local arrays and broadcast
	MPI_Allreduce(br, Br, 2*x+2, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	
	Ar[2*x] = 0;				//Fill in First Position(i = 0)
	Ar[2*x + 1] = n/2;			//Fill in Last Position(i = n)
	Br[2*x] = 0;				//Fill in First Position(i = 0)
	Br[2*x + 1] = n/2;			//Fill in Last Position(i = n)
	insertionSort(Ar, 2*x+2);			//Sort these to make sure they are in order
	insertionSort(Br, 2*x+2);			//Sorted using insertion sort
	
	int* localC = new int[n];		//Local C array for the sequnetial merges performed by individual process by striping.
	clear(localC, n);
	
	//Striping the sequential merges
	for(int i = my_rank; i < 2*x+1; i+=p){
		smerge(A, Ar[i], Ar[i+1], B, Br[i], Br[i+1], localC, Ar[i] + Br[i], Ar[i+1]+Br[i+1] );
	}
	
	//Combine each process merges into the resulting C array
	MPI_Allreduce(localC, C, n, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	
}

void merge_sort(int*a, int*b, int n, int my_rank, int p){
    //MYLES: We have made a base for arrays less than 4 just use insertion because its better, could use another sort but I saw this one in a sort visualizer video on youtube and thought it looked cool so why not?
    if(n==4){
        for(int i = 0; i < n; i++){
            b[i] = a[i];
        }
        insertionSort(b, n);
    }
    else{
        int*c = new int[n]; //basically we are creating a temp array to store the two halfs then run mergesort recursively
        clear(c,n);
        merge_sort(&a[0], &c[0], n/2, my_rank, p);
        merge_sort(&a[n/2], &c[n/2], n - n/2, my_rank, p);
        pmerge(&c[0], &c[n/2], &b[0], n, my_rank, p);
    }
}

int main(int argc, char * argv[]) {

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
    bool results = true;

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
    MPI_Bcast(a, (size/2), MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(b, size/2, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    
    merge_sort(a, b, size, my_rank, p);
    if(results && my_rank==0){
        cout<<"\nA:" << endl;
        printArray(a, size);
        cout << "\nB:"<< endl;
        printArray(b, size);
    }

	// Shut down MPI
	MPI_Finalize();

	return 0;
}