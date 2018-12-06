//Nate Hamling
//hamli057

/* 
 * trans.c - Matrix transpose B = A^T
 *
 * Each transpose function must have a prototype of the form:
 * void trans(int M, int N, int A[N][M], int B[M][N]);
 *
 * A transpose function is evaluated by counting the number of misses
 * on a 1KB direct mapped cache with a block size of 32 bytes.
 */ 
#include <stdio.h>
#include "cachelab.h"

int is_transpose(int M, int N, int A[N][M], int B[M][N]);
void helper_32x32(int M, int N, int A[N][M], int B[M][N]);
void helper_64x64(int M, int N, int A[N][M], int B[M][N]);
void helper_61x67(int M, int N, int A[N][M], int B[M][N]);

/* 
 * transpose_submit - This is the solution transpose function that you
 *     will be graded on for Part B of the assignment. Do not change
 *     the description string "Transpose submission", as the driver
 *     searches for that string to identify the transpose function to
 *     be graded. 
 */
char transpose_submit_desc[] = "Transpose submission";
void transpose_submit(int M, int N, int A[N][M], int B[M][N])
{
	if (M == 32) {
		helper_32x32(M,N,A,B);
		return;
	} 
	else if (M == 64) {
		helper_64x64(M,N,A,B);
		return;
	}
	else{
		helper_61x67(M,N,A,B);
		return;
	}
}

/* 
 * You can define additional transpose functions below. We've defined
 * a simple one below to help you get started. 
 */ 

//Helper function to transpose 32x32 matrix
char helper_32x32_desc[] = "Transpose 32x32";
void helper_32x32(int M, int N, int A[N][M], int B[M][N]) {
	int block_size = 8;
	int row_block, col_block; //index through rows/cols in each block
	int r, c; 				  //index through rows and columns
	int d_idx; 				  //diagonal index
	int d_el; 				  //diagonal element

	for(c = 0; c < M; c += block_size){
		for(r = 0; r < N; r += block_size){
			for(row_block = r; row_block < r + block_size; row_block++){
				for(col_block = c; col_block < c + block_size; col_block++){
					if(row_block != col_block){
						B[col_block][row_block] = A[row_block][col_block];
					}
					else{
						d_idx = row_block;
						d_el = A[row_block][col_block];
					}
				}
				if(c == r){
					B[d_idx][d_idx] = d_el;
				}
			}
		}
	}
}

//Helper function to transpose 64x64 matrix
char helper_64x64_desc[] = "Transpose 64x64";
void helper_64x64(int M, int N, int A[N][M], int B[M][N]){
	int r, c;
	int a1, a2, a3, a4, a5, a6, a7, a8;		//storage to avoid future misses
	int count; //used for iteration
	for(r = 0; r < M; r += 8){
		for(c = 0; c < M; c += 8){
			for(count = 0; count < 8; count++){
				//stores atart of each row on first iteration of each block
				if(count == 0){
					a5 = A[c+count][r+4];
					a6 = A[c+count][r+5];
					a7 = A[c+count][r+6];
					a8 = A[c+count][r+7];
				}
				//store current row
				a1 = A[c+count][r];
				a2 = A[c+count][r+1];
				a3 = A[c+count][r+2];
				a4 = A[c+count][r+3];

				//transpose B in 64s
				B[r][c+count] = a1;
				B[r][c+count+64] = a2;
				B[r][c+count+128] = a3;
				B[r][c+count+192] = a4;
			}
			for(count = 7; count > 0; count--){
				a1 = A[c+count][r+4];
				a2 = A[c+count][r+5];
				a3 = A[c+count][r+6];
				a3 = A[c+count][r+7];

				B[r+4][c+count] = a1;
				B[r+4][c+count+64] = a2;
				B[r+4][c+count+128] = a3;
				B[r+4][c+count+192] = a4;
			}
			B[r+4][c] = a5;
			B[r+4][c+64] = a6;
			B[r+4][c+128] = a7;
			B[r+4][c+192] = a8;
		}
	}
}

//Helper function to transpose 32x32 matrix
char helper_61x67_desc[] = "Transpose 61x67";
void helper_61x67(int M, int N, int A[N][M], int B[M][N]){
	int row_block, col_block;
	int r, c;
	int tmp;
	for(r = 0; r < N; r += 16){
		for(c = 0; c < M; c += 4){
			//keeps transpose within bounds of N
			for(row_block = r; (row_block < r + 16) && (row_block < N); ++row_block){
				//keeps transpose within bounds of M
				for(col_block = c; (col_block < c + 4) && (col_block < M); ++col_block){
					//save diagonal in tmp
					if(row_block - r == col_block - c){
						tmp = A[row_block][col_block];
					}
					else{
						B[col_block][row_block] = A[row_block][col_block];
					}
				}
				for(col_block = c; (col_block < c + 4) && (col_block < M); ++col_block){
					if(row_block - r == col_block - c){
						B[col_block][row_block] = tmp;
					}
				}
			}
		}
	}
}

/* 
 * trans - A simple baseline transpose function, not optimized for the cache.
 */
char trans_desc[] = "Simple row-wise scan transpose";
void trans(int M, int N, int A[N][M], int B[M][N])
{
    int i, j, tmp;

    for (i = 0; i < N; i++) {
        for (j = 0; j < M; j++) {
            tmp = A[i][j];
            B[j][i] = tmp;
        }
    }    

}

/*
 * registerFunctions - This function registers your transpose
 *     functions with the driver.  At runtime, the driver will
 *     evaluate each of the registered functions and summarize their
 *     performance. This is a handy way to experiment with different
 *     transpose strategies.
 */
void registerFunctions()
{
    /* Register your solution function */
    registerTransFunction(transpose_submit, transpose_submit_desc); 

    /* Register any additional transpose functions */
    //registerTransFunction(trans, trans_desc); 

}

/* 
 * is_transpose - This helper function checks if B is the transpose of
 *     A. You can check the correctness of your transpose by calling
 *     it before returning from the transpose function.
 */
int is_transpose(int M, int N, int A[N][M], int B[M][N])
{
    int i, j;

    for (i = 0; i < N; i++) {
        for (j = 0; j < M; ++j) {
            if (A[i][j] != B[j][i]) {
                return 0;
            }
        }
    }
    return 1;
}

