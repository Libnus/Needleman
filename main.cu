#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// #include <cuda.h>

#define MATCH 1
#define MISMATCH -1
#define GAP -2

#define NUM_THREADS 256

__global__ void NW_kernel(int* d_data, int* d_dataComputed, const char* sequence_a, const char* sequence_b){
    // get index
    for( int index = blockIdx.x * blockDim.x + threadIdx.x; index < rows*cols; index += blockDim.x*gridDim.x ){
        int y = floor(index/(double)rows);
        int x = index % cols;
        if((y == 0 || x == 0) || (y > rows || x > cols)) continue;


        int top_index = index - cols;
        int left_index = index-1;
        int diagonal_index = top_index-1;

        while( !d_dataComputed[top_index] || !d_dataComputed[left_index] || !d_dataComputed[diagonal_index] ) { continue; } // wait for values to be filled in d_data

        int top_score = d_data[top_index] * GAP;
        int left_score = d_data[left_index] * GAP;
        int diagonal_score = sequence_a[x] == sequence_b[x] ? d_data[diagonal_index] * MATCH : d_data[diagonal_index] * MISMATCH;

        d_dataComputed[index] = 1;
    }

}


void NW_kernelLaunch( int* d_data, int* d_dataComputed, int rows, int cols, int grid_length, char sequence_a, char sequence_b ){
    // controller for thread launching
    int diagonals = grid_length - 1;

    // #ifdef DIAGONAL_SERIAL
    // // for each diagonal launch a kernel
    // for(int i = 0; i < diagonals; i++){
    //     NW_kernel<<<1, NUM_THREADS>>>(d_data, d_dataComputed, seqeunce_a, sequence_b, rows, cols);

    //     cudaDeviceSynchronize(); // sync before moving onto next diagonal
    // }
    // #endif

    #ifdef DIAGONAL_PARAALLEL
    NW_kernel<<<diagonals, NUM_THREADS>>>(d_data, sequence_a, sequence_b, rows, cols);
    #endif
}

static inline int min(int a, int b) { return a < b ? a : b; }
static inline int max(int a, int b) { return a > b ? a : b; }

static inline void NW_init(int* d_data, int* d_dataComputed, int rows, int cols){
    int diff_min = min(rows, cols);
    int diff_max = max(rows, cols) - diff_min;

    d_data[0] = 0;
    d_dataComputed[0] = 1;


    // o(n) despite two loops :)
    int i;
    for(i = 1; i < diff_min; i++){
        // d_data[0][i] = d_data[0][i-1]-1;
        // d_data[i][0] = d_data[i-1][0]-1;
        int row_index = i*cols;

        printf("row_index %d, cols %d\n", row_index, cols);
        d_data[row_index] = d_data[row_index-cols] - 1; // rows
        d_data[i] = d_data[i-1] - 1; //cols

        d_dataComputed[row_index] = 1;
        d_dataComputed[i] = 0;
    }

    if (rows == cols) return;

    int rows_bigger = 0;
    if(rows > cols) rows_bigger = 1;

    for(; i < diff_max; i++){
        if(rows_bigger) {
            int row_index = i*cols;

            d_data[i*cols] = d_data[row_index-cols] - 1; // rows
            d_dataComputed[row_index] = 1;
        }
        else{
            d_data[i] = d_data[i-1] - 1; // cols
            d_dataComputed[i] = 1;
        }
    }
}

// needleman-wunsch algorithm
int main(int argc, char* argv[]){
    if(argc > 1){
        printf("Incorrect number of inputs\n");
        return EXIT_FAILURE;
    }

    // test sequences
    char* sequence_a = "GCATGCG";
    char* sequence_b = "GATTACA";

    // construct a 1d array to represent the 2d grid for needleman
    int grid_length = (strlen(sequence_a)+1) + (strlen(sequence_b)+1);
    int cols = strlen(sequence_a)+1;
    int rows = strlen(sequence_b)+1;

    int* d_data = calloc(grid_length, sizeof(int));
    int* d_dataComputed = calloc(grid_length, sizeof(int));

    NW_init(d_data, d_dataComputed, rows, cols);

    for(size_t i = 0; i < rows; i++){
        for(size_t j = 0; j < cols; j++){
            int index = i*cols + j;
            printf("%d ", d_data[index]);
        }
        printf("\n");
    }

    NW_kernelLaunch( d_data, d_dataComputed, rows, cols, grid_length, sequence_a, sequence_b );


    return EXIT_SUCCESS;
}
