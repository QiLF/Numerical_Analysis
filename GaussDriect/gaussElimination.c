/**
 * HW6
 * qlf
 * result:
x=(
-0.289233816015754,
0.345435715779115,
-0.712811731086879,
-0.220608510570529,
-0.430400432704022,
0.154308739838311,
-0.057822873289041,
0.201053894823681,
0.290228661879745,
)
 */

#include <stdio.h>
#include <math.h>

#define MATRIX_SIZE 9

/**
 * @brief  Gauss Elimination Algorithm with Scaled Row Pivoting for Ax = b
 * @param A input: n*n Matrix
 * @param b input: b = (b1, b2,...bn)    output: x = (x1, x2,...xn)
 * @param n input: size of A
 */
void gaussElimination(double A[][MATRIX_SIZE], double *b, int n) {
    for (int i = 0; i < n; i++) {
        //1.找列主元
        int k = i;
        for (int j = i + 1; j < n; j++) {
            if (fabs(A[j][i]) > fabs(A[k][i])) {
                k = j;
            }
        }
        if(k!=i){
            //2.交换扩展矩阵(A,b)两行
            for (int j = i; j < n; j++) {
                //swap(A[i][j], A[k][j])
                double tmp = A[i][j];
                A[i][j] = A[k][j];
                A[k][j] = tmp;
            }
            //swap(b[k], b[i])
            double tmp = b[i];
            b[i] = b[k];
            b[k] = tmp;
        }
        //3.消元
        for (int j = i + 1; j < n; j++) {
            double tmp = - A[j][i] / A[i][i];

            for (k = i + 1; k < n; k++) {
                A[j][k] = A[j][k] + tmp * A[i][k];
            }
            b[j] = b[j] + tmp * b[i];
        }
    }
    //4.回代求解
    for (int i = n - 1; i >= 0; i--) {
        for (int j = i + 1; j < n; j++) {
            b[i] = b[i] - A[i][j] * b[j];
        }
        b[i] = b[i] / A[i][i];
    }
}


int main() {

    //
    double A[9][9] = {
            {31,  -13, 0,   0,   0,   -10, 0,   0,  0},
            {-13, 35,  -9,  0,   -11, 0,   0,   0,  0},
            {0,   -9,  31,  -10, 0,   0,   0,   0,  0},
            {0,   0,   -10, 79,  -30, 0,   0,   0,  -9},
            {0,   0,   0,   -30, 57,  -7,  0,   -5, 0},
            {0,   0,   0,   0,   -7,  47,  -30, 0,  0},
            {0,   0,   0,   0,   0,   -30, 41,  0,  0},
            {0,   0,   0,   0,   -5,  0,   0,   27, -2},
            {0,   0,   0,   -9,  0,   0,   0,   -2, 29}
    };
    double b[9] = {-15, 27, -23, 0, -20, 12, -7, 7, 10};
    int n = 9;
    gaussElimination(A, b, n);
    //
    for (int i = 0; i < n; i++) {
        if (i == 0) {
            printf("x=(\n%.15f,\n", b[0]);
        } else {
            printf("%.15f,\n", b[i]);
        }
    }
    printf("\b)\n");
    return 0;
}