/*  Myles Wright, Cathal Osullivan
    Butler University 2021
    CS452 PMERGE
*/


//GUIDES THAT WE USED:
//https://github.com/RachelBurke/CS452/blob/master/Project3/cloneWars.cpp
//https://cse.buffalo.edu/faculty/miller/Courses/CSE702/Swati.Nair-Fall-2018.pdf
//https://www.youtube.com/watch?v=_XOZ2IiP2nw
//http://selkie.macalester.edu/csinparallel/modules/ParallelSorting/build/html/MergeSort/MergeSort.html
//https://stanford.edu/~rezab/classes/cme323/S16/notes/Lecture04/cme323_lec4.pdf

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
int srank(int * arr, int x, int valtofind);
void insertionSort(int* array, int size);
void smerge(int* A, int a_start, int a_end, int* B, int b_start, int b_end, int* C, int c_start, int c_end);
void printArray(int * a, int size);
void pmerge(int* A, int* B, int* C, int n, int my_rank, int p);
void mergesort(int*a, int*b, int n, int my_rank, int p);

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
        cout << "array A:" << endl;
        printArray(a, size_div);
        cout << "array B:" << endl;
        printArray(b, size_div);
        
    }
    
    MPI_Bcast(&size, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(a, (size/2), MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(b, size/2, MPI_INT, 0, MPI_COMM_WORLD);
    
    cout << "process " << my_rank << "a:" << endl;
    cout << "boots withthe fur" << endl;
    merge_sort(a, b, size, my_rank, p);
    cout << "ASS" << endl;
    if(my_rank==0){
        cout<<"\nA:" << endl;
        printArray(a, size);
        cout << "\nB:"<< endl;
        printArray(b, size);
    }

	// Shut down MPI
	MPI_Finalize();

	return 0;
}



void insertionSort(int* array, int size){//CATHAL: hehe
    cout << "inside insertion" << endl;
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

int srank(int * arr, int x, int valtofind){ //MYLES: a fun binary search style algorithm that returns the value of where the value should be
    cout << "inside of srank" << endl;
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

void smerge(int* a, int a_start, int a_end, int* b, int b_start, int b_end, int* c, int c_start, int c_end)
{
    cout << "Inside of smerge" << endl;
	for(int i = c_start; i < c_end; i++)
	{
		if(a_start < a_end && (b_start == b_end || a[a_start] < b[b_start])){
			c[i] = a[a_start++];
		}
		else{		
			c[i] = b[b_start++];
		}
	}
}

void printArray(int * a, int size){
    for(int i = 0; i < size; i++){
        cout << a[i] << " ";
    }
    cout << "\n";
}

void clear(int*arr, int x){ //fills with
    for(int i = 0; i<x; i++){
        arr[i] = 0;
    }
}


void pmerge(int* a, int* b, int* c, int x, int my_rank, int p)		//Merge Array A and B, both size n/2, into array C, size n
{
    cout << "inside of pmerge" << endl;
	int logn = log2(x/2);
	int y = ceil((x/2)/logn);		//Number of sampled elements
	
	
    int rsize = 2*y+2;
	int* rankA = new int[rsize];
	int* rankB = new int[rsize];
	int* Ar = new int[rsize];
	int* Br = new int[rsize];
	clear(rankA,rsize);
	clear(rankB,rsize);
	clear(Ar,rsize);
	clear(Br,rsize);
	
    

	for(int i = my_rank; i < x; i+=p){		
		rankA[i] = 1 + i*logn;
		rankB[i] = 1 + i*logn;
	}
    cout << "test 2" << endl;		
	for(int i = my_rank; i < y; i+=p){
        cout << "going to srank ONE" << endl;		
		rankA[i+x] = srank(a, x/2, A[0+i*logn]);
        cout << "going to srank TWO" << endl;
		rankB[i+x] = srank(b, x/2, B[0+i*logn]);
        cout << "done with srank" << endl;
	}

    cout << "rankA: " << endl;
    printArray(rankA, rsize);
    cout << "rankB: " << endl;
    printArray(rankB, rsize);
    cout << "Ar: " << endl;
    printArray(Ar, rsize);
    cout << "Br: " << endl;
    printArray(Br, rsize);

    cout << "Test 3" << endl;
	MPI_Allreduce(rankA, Ar, rsize, MPI_INT, MPI_SUM, MPI_COMM_WORLD);		
    cout << "Allreduce 1 done" << endl;
	MPI_Allreduce(rankB, Br, rsize, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	cout << "Allreduce done" << endl;

	Ar[2*y] = 0;				
	Ar[2*y + 1] = x/2;			
	Br[2*y] = 0;				
	Br[2*y + 1] = x/2;


    cout << "rank arrays before sort: " << endl;
    printArray(Ar, 2*y);
    printArray(Br, 2*y);
	insertionSort(Ar, rsize);			
	insertionSort(Br, rsize);			
	cout << "rank arrays after sort: " << endl;
    printArray(Ar, 2*y);
    printArray(Br, 2*y);
	int* localC = new int[x];		
	clear(localC, y);
	

	for(int i = my_rank; i < 2*x+1; i+=p){
		smerge(a, Ar[i], Ar[i+1], b, Br[i], Br[i+1], localC, Ar[i] + Br[i], Ar[i+1]+Br[i+1]);
	}
	
	//Combine each process merges into the resulting C array
	MPI_Allreduce(localC, c, x, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    cout << "array C:" << endl;
    printArray(c, x);
	
}

void mergesort(int*a, int*b, int n, int my_rank, int p){
    //MYLES: We have made a base for arrays less than 4 just use insertion because its better, could use another sort but I saw this one in a sort visualizer video on youtube and thought it looked cool so why not?
    cout << "Inside of merge_sort" << endl;
    if(n==4){
        for(int i = 0; i < n; i++){
            b[i] = a[i];
        }
        insertionSort(b, n);
    }
    else{
        int*c = new int[n]; //basically we are creating a temp array to store the two halfs then run mergesort recursively
        clear(c,n);
        mergesort(&a[0], &c[0], n/2, my_rank, p);
        mergesort(&a[n/2], &c[n/2], n - n/2, my_rank, p);
        pmerge(&c[0], &c[n/2], &b[0], n, my_rank, p);
    }
}