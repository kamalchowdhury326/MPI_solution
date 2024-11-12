/* File:     hw4_kamal.c
 *Name: Md Kamal Hossain Chowdhury
 * 
 * Homework #:4
 * 
 *
 * Compile:  mpicc -g -Wall -o hw4_Isendrecv hw4_Isendrecv.c
 * Run:      mpiexec -n <comm_sz> ./hw4_Isendrecv 4 4 /scratch/ualmkc001/
 *
 * 
 */
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <sys/time.h>
#include <string.h>

#define DIES   0
#define ALIVE  1

void Check_for_error(int local_ok, char fname[], char message[], 
      MPI_Comm comm);

void Print_vector(int local_b[], int *counts, int *displs, int n, char title[], 
      int my_rank, MPI_Comm comm);
void Parallel_vector_sum(int local_x[], int local_y[], 
      int local_z[], int local_n);
void Print_local_vector( int local_b[],int *counts,char title[],int my_rank,      MPI_Comm  comm);
void compute_local(int local_x[], int n, int NTIMES,int counts[], int displs[], int my_rank,int comm_sz, MPI_Comm  comm);



char* filePath;
double gettime(void) {
  struct timeval tval;

  gettimeofday(&tval, NULL);

  return( (double)tval.tv_sec + (double)tval.tv_usec/1000000.0 );
}
void printarray(int **a, int row,int col, int k) {
  int i, j;
  printf("Life after %d iterations:\n", k) ;
  for (i = 1; i < row+1; i++) {
    for (j = 1; j< col+1; j++)
      printf("%d ", a[i][j]);
    printf("\n");
  }
  printf("\n");
}
void printarray_1d(int *a, int N, int k) {
  int i, j;
  printf("Life after %d iterations:\n", k) ;
  for (i = 0; i < N; i++) {
    for (j = 0; j< N; j++)
      printf("%d ", a[i*N+j]);
    printf("\n");
  }
  printf("\n");
}

void allocate_vector_for_x(
    int **local_x_pp /* out */,
    int local_n /* in  */,
    int n /* in*/,
    MPI_Comm comm /* in  */)
{
  int local_ok = 1;
  char *fname = "Allocate_vectors";

  *local_x_pp = malloc(local_n * sizeof(int));

  if (*local_x_pp == NULL)
    local_ok = 0;
  Check_for_error(local_ok, fname, "Can't allocate local vector(s)",
                  comm);
}
int **allocarray(int P, int Q) {
  int i, *p, **a;

  p = (int *)malloc(P*Q*sizeof(int));
  a = (int **)malloc(P*sizeof(int*));
  for (i = 0; i < P; i++)
    a[i] = &p[i*Q]; 

  return a;
}
void Read_vector(
      int    local_a[]   /* out */, 
      int       counts[]    /* in  */, 
      int       displs[]    /* in  */, 
      int       n           /* in  */,
      char      vec_name[]  /* in  */,
      int       my_rank     /* in  */, 
      MPI_Comm  comm        /* in  */) {

   int* a = NULL;
   int i, local_n;
   int local_ok = 1;
   char* fname = "Read_vector";

   local_n = counts[my_rank];
   if (my_rank == 0) {
      a = malloc(n*n*sizeof(int));
      
      if (a == NULL) local_ok = 0;
      Check_for_error(local_ok, fname, "Can't allocate temporary vector", 
            comm);
     
        /* Initialize the life matrix */
        for (i = 0; i < n*n; i++) {
        srand(54321|i);
        if (drand48() < 0.5) 
            a[i] = ALIVE ;
        else
            a[i] = DIES ;
        }
      
    #ifdef DEBUG2
        printf("rank=%d\n",my_rank);
        /* Display the life matrix */
        printarray_1d(a, n, 0);
    #endif
      
      MPI_Scatterv(a, counts, displs, MPI_INT, local_a, local_n, MPI_INT, 0,comm);
     
      free(a);
   } else {
        Check_for_error(local_ok, fname, "Can't allocate temporary vector", 
              comm);
        MPI_Scatterv(a, counts, displs, MPI_INT, local_a, local_n, MPI_INT, 0,comm);
   }
  
}  /* Read_vector */  


void freearray(int **a) {
  free(&a[0][0]);
  free(a);
}

void file_write(char* path, int **A, int m, int n,int iterations, int P){
   
    const char* name = path;
    const char* extension = ".txt";
    char size[5];
    snprintf(size, 5, "%d", n);
    char iter[5];
    snprintf(iter, 5, "%d", iterations);
    char p_char[5];
    snprintf(p_char, 5, "%d", P);
    

    char* name_with_extension;
    name_with_extension = malloc(strlen(name)+1+25); /* make space for the new string (should check the return value ...) */
    strcpy(name_with_extension, name); /* copy name into the new var */
    strcat(name_with_extension,"output");
    strcat(name_with_extension,".");
    strcat(name_with_extension,size);
    strcat(name_with_extension,".");
    strcat(name_with_extension,iter);
  
    strcat(name_with_extension,".");
    strcat(name_with_extension,p_char);

    strcat(name_with_extension, extension); /* add the extension */
    
    
    // Open file in write mode
    FILE *file = fopen(name_with_extension, "w");

    // Check if file creation is successful
    if (file == NULL) {
        printf("Failed to create the file.\n");
        return ;
    }
   
    for(int i=0;i<n;i++){
       for(int j=0;j<n;j++){
          // Write something to the file
         fprintf(file, " %d ",A[i][j]);
       }
       fprintf(file, "\n");
    }
    

    // Close the file
    fclose(file);
    printf("File created successfully at %s\n", name_with_extension);
   // free(name_with_extension);
}

void print_life(int **life,int nRowsGhost,int nColsGhost,int my_rank,int local_n){
	int i,j;
	 printf("\n----life print----\n");
   
   printf("my_rank=%d local_n=%d \n",my_rank,local_n);
      for (i = 0; i < nRowsGhost; i++)
	  {
		  for(j=0;j<nColsGhost;j++){
			  
			printf(" %3d ",life[i][j]);  
		  }
         
      printf("\n");
	  }
}

int compute(int **life, int **temp, int nRows,int nCols) {
  int i, j, value, flag=0;

  for (i = 1; i < nRows+1; i++) {
    for (j = 1; j < nCols+1; j++) {
      /* find out the value of the current cell */
      value = life[i-1][j-1] + life[i-1][j] + life[i-1][j+1]
	+ life[i][j-1] + life[i][j+1]
	+ life[i+1][j-1] + life[i+1][j] + life[i+1][j+1] ;
      
      /* check if the cell dies or life is born */
      if (life[i][j]) { // cell was alive in the earlier iteration
	if (value < 2 || value > 3) {
	  temp[i][j] = DIES ;
	  flag++; // value changed 
	}
	else // value must be 2 or 3, so no need to check explicitly
	  temp[i][j] = ALIVE ; // no change
      } 
      else { // cell was dead in the earlier iteration
	if (value == 3) {
	  temp[i][j] = ALIVE;
	  flag++; // value changed 
	}
	else
	  temp[i][j] = DIES; // no change
      }
    }
  }

  return flag;
}

/*-------------------------------------------------------------------*/
int main(int argc, char **argv) {
   int n, local_n, i, remain;
   int comm_sz, my_rank;
   int *local_x;
   
   
   
   MPI_Comm comm;
   int *counts, *displs;

   MPI_Init(NULL, NULL);
   comm = MPI_COMM_WORLD;
   MPI_Comm_size(comm, &comm_sz);
   MPI_Comm_rank(comm, &my_rank);
   
   if (argc != 4) {
      printf("Usage: mpiexec -n <number_of_process> %s <N> <Generations> <filePath>\n", argv[0]);
      exit(-1);
    }
   n = atoi(argv[1]);
   int NTIMES = atoi(argv[2]);
   filePath=argv[3];
   //n=4;
  
  
   /* compute counts and displacements */
   counts = (int *)malloc(comm_sz*sizeof(int));
   displs = (int *)malloc(comm_sz*sizeof(int));
   remain = n % comm_sz;
   for (i = 0; i < comm_sz; i++){
       counts[i] = n/comm_sz + ((i < remain)? 1 : 0);
       counts[i]=counts[i]*n;
       #ifdef DEBUG3
       printf("counts[%d]=%d\n",i,counts[i]);
       #endif
   }
   displs[0] = 0;
   for (i = 1; i < comm_sz; i++){
       displs[i] = displs[i-1] + counts[i-1];
       #ifdef DEBUG3
       printf("displs[%d]=%d\n",i,displs[i]);
       #endif
   }
   local_n = counts[my_rank];

#  ifdef DEBUG
   printf("Proc %d > n = %d, local_n = %d\n", my_rank, n, local_n);
#  endif
    
   // allocate1DLife(&local_x, local_n,n, comm);
   allocate_vector_for_x(&local_x, local_n, n, comm);
   
    Read_vector(local_x, counts, displs, n, "x", my_rank, comm);
   
    compute_local(local_x, n,NTIMES,counts,displs,my_rank,comm_sz,comm);
   
    free(local_x);
  
    free(counts);
    free(displs);
  
  
   MPI_Finalize();

   return 0;
}  /* main */

/*-------------------------------------------------------------------
 * Function:  Check_for_error
 * Purpose:   Check whether any process has found an error.  If so,
 *            print message and terminate all processes.  Otherwise,
 *            continue execution.
 * In args:   local_ok:  1 if calling process has found an error, 0
 *               otherwise
 *            fname:     name of function calling Check_for_error
 *            message:   message to print if there's an error
 *            comm:      communicator containing processes calling
 *                       Check_for_error:  should be MPI_COMM_WORLD.
 *
 * Note:
 *    The communicator containing the processes calling Check_for_error
 *    should be MPI_COMM_WORLD.
 */
void Check_for_error(
      int       local_ok   /* in */, 
      char      fname[]    /* in */,
      char      message[]  /* in */, 
      MPI_Comm  comm       /* in */) {
   int ok;

   MPI_Allreduce(&local_ok, &ok, 1, MPI_INT, MPI_MIN, comm);
   if (ok == 0) {
      int my_rank;
      MPI_Comm_rank(comm, &my_rank);
      if (my_rank == 0) {
         fprintf(stderr, "Proc %d > In %s, %s\n", my_rank, fname, 
               message);
         fflush(stderr);
      }
      MPI_Finalize();
      exit(-1);
   }
}  /* Check_for_error */






void compute_local(
      int    local_x[]  /* in */,
      int n,	          /* in */
      int NTIMES, 
	    int       counts[]    /* in  */, 
      int       displs[]    /* in  */,
      int       my_rank    /* in */,
      int       comm_sz    /* in */,
      MPI_Comm  comm       /* in */) {

   
   int i, j,local_n;
   int **life=NULL, **temp=NULL, **ptr ;
   int **final_board2D=NULL;
   local_n = counts[my_rank];
   int nCols=n;
   int nRows=local_n/n;
   int upper_rank,down_rank;
   double t1, t2;
   int flag=1,k;
   
  int *final_board=NULL;
   if(my_rank==0){
     
     final_board=malloc(n*n*sizeof(int));
   } 
   
   MPI_Status status;
   MPI_Request request = MPI_REQUEST_NULL;


   int nRowsGhost=nRows+2;
   int nColsGhost=nCols+2;
   life = allocarray(nRowsGhost,nColsGhost);
   temp = allocarray(nRowsGhost,nColsGhost);
  
   
   int row=0;
   int col=0;
   
  for (i = 0; i < nRowsGhost; i++) {
	  for(j=0;j<nColsGhost;j++){
     if(i==0 || j==0 ||i==nRowsGhost-1 ||j==nColsGhost-1)
	 {     life[i][j] = DIES ;
         temp[i][j]  = DIES ;
	 }
  }
  }
  
  
   for (i = 0; i < local_n; i++){
	   row=i/n;
	   col=i%n;
	   // give space for ghost cell
	   row=row+1;
	   col=col+1;
    #ifdef DEBUG3
	   printf("\n----local_x[%d]=%d-----",i,local_x[i]);
    #endif
	   life[row][col]=local_x[i];
    #ifdef DEBUG3
	   printf(" life[%d][%d]= %d \n",row,col,life[row][col]);
    #endif
   }
   #ifdef DEBUG2
    printf("rank=%d\n",my_rank);
    /* Display the life matrix */
    printarray(life, nRows,nCols, 0);
   #endif
   
   upper_rank = my_rank + 1;
   if (upper_rank >= comm_sz) upper_rank = MPI_PROC_NULL;
   down_rank = my_rank - 1;
   if (down_rank < 0) down_rank = MPI_PROC_NULL;
   
   if(my_rank==0){
      t1 = gettime();
      }
   /* Play the game of life for given number of iterations */
    for (k = 0; k < NTIMES; k++) {
      flag = 0;
     
      if ((my_rank % 2) == 0) {
        /* exchange up */
        MPI_Isend(&(life[nRows][0]),nColsGhost,MPI_INT,upper_rank,0,comm,&request);
        MPI_Irecv(&(life[nRows+1][0]), nColsGhost, MPI_INT, upper_rank, 0,comm, &request );
        
        MPI_Wait(&request, &status);
        // MPI_Sendrecv( &(life[nRows][0]), nColsGhost, MPI_INT, upper_rank, 0, 
        //         &(life[nRows+1][0]), nColsGhost, MPI_INT, upper_rank, 0, 
        //         comm, &status );
          
          }
          else {
        /* exchange down */
        MPI_Isend(&(life[1][0]),nColsGhost,MPI_INT,down_rank,0,comm,&request);
        MPI_Irecv(&(life[0][0]), nColsGhost, MPI_INT, down_rank, 0,comm, &request);
        MPI_Wait(&request, &status); //blocks and waits for destination process to receive data
        
        // MPI_Sendrecv( &(life[1][0]), nColsGhost, MPI_INT, down_rank, 0,
        //         &(life[0][0]), nColsGhost, MPI_INT, down_rank, 0, 
        //         comm, &status );
          }

      /* Do the second set of exchanges */
      if ((my_rank % 2) == 1) {
        /* exchange up */
        MPI_Isend(&(life[nRows][0]),nColsGhost,MPI_INT,upper_rank,0,comm,&request);
	      MPI_Irecv(&(life[nRows+1][0]), nColsGhost, MPI_INT, upper_rank, 0,comm, &request );
	      MPI_Wait(&request, &status); //blocks and waits for destination process to receive data
	
        // MPI_Sendrecv( &(life[nRows][0]), nColsGhost, MPI_INT, upper_rank, 1, 
        //         &(life[nRows+1][0]), nColsGhost, MPI_INT, upper_rank, 1, 
        //         comm, &status );
          }
      else {
        /* exchange down */
        MPI_Isend(&(life[1][0]),nColsGhost,MPI_INT,down_rank,0,comm,&request);
        MPI_Irecv(&(life[0][0]), nColsGhost, MPI_INT, down_rank, 0,comm, &request);
        MPI_Wait(&request, &status); //blocks and waits for destination process to receive data
        
        // MPI_Sendrecv( &(life[1][0]), nColsGhost, MPI_INT, down_rank, 1,
        //         &(life[0][0]), nColsGhost, MPI_INT, down_rank, 1, 
        //         comm, &status );
        //   
        }
   
  
   flag=compute(life,temp,nRows,nCols);
  
  
  // Each MPI process sends its rank to reduction, root MPI process collects the result
    int reduction_flag = 0;
    MPI_Allreduce(&flag, &reduction_flag, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
 
    if(my_rank == 0)
    {   
      if(!reduction_flag)
        printf("The sum of all flag is %d after k=%d.\n", reduction_flag,k);
    }
    if(!reduction_flag){
      break;
    }

    MPI_Barrier(comm);
    
    
    /* copy the new values to the old array */
    ptr = life;
    life = temp;
    temp = ptr;
   
  
#ifdef DEBUG2
    /* Print no. of cells alive after the current iteration */
    printf("No. of cells whose value changed in iteration %d = %d\n",k+1,flag) ;
    printf("rank=%d\n",my_rank);
    /* Display the life matrix */
    printarray(life, row,col, k+1);
#endif
  }

  for(int i=1;i<nRows+1;i++){
    for(int j=1;j<nCols+1;j++){
      local_x[(i-1)*n+(j-1)]=life[i][j];
      // printf("local_x[%d]=%d\n",(i-1)*n+(j-1),local_x[(i-1)*n+(j-1)]);

    }

  }
  freearray(life);
  freearray(temp);
  // MPI_Barrier(comm);
  
  MPI_Gatherv(local_x,local_n,MPI_INT,final_board,counts,displs,MPI_INT,0,comm);
  
  
  if(my_rank==0){
    #ifdef DEBUG2
    printf("Final Board:\n");
    printarray_1d(final_board, n, k+1);
    #endif
    final_board2D = allocarray(n,n);
    for(int i=0;i<n;i++){
      for(int j=0;j<n;j++){
        final_board2D[i][j]=final_board[i*n+j];
      }
    }
  }
  
  if(my_rank==0){
    t2 = gettime();
    printf("Time taken %f seconds for size=%dx%d  generations=%d after %d iterations\n", t2 - t1, n,n,NTIMES,k);
    //printf("%s\n",filePath);
    file_write(filePath,final_board2D,n,n,NTIMES,comm_sz);
    free(final_board);
    freearray(final_board2D);
  }
  
  
  // freearray(life);
  // freearray(temp);
  
}  /* Print_life */

