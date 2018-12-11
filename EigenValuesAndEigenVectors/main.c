/*
* HW8
* Description: Get eigenvalues and eigenvectors of matrix by Jacobi method
* Author: qlf
* Date: 2018-12-11
* Result:
r1=-4.061950892717061
v1=(0.768639921841930,-0.221691402575936,-0.369283259380001,-0.395445295112373,-0.259419516401363)
r2=-12.760554989645305
v2=(0.089240149713038,0.591435168953867,-0.440770459627694,0.525458539215943,-0.414554413905271)
r3=11.023392437654190
v3=(0.411866187348695,-0.342988626053424,0.542906882860872,0.615543660391281,-0.197694625145014)
r4=4.574951567326774
v4=(-0.103969141997373,-0.546298068116902,-0.607111450341849,0.422399393283991,0.379134906993808)
r5=16.224161877381380
v5=(0.469878204922516,0.430086719299899,0.077582727893399,0.101000479424375,0.760276074797378)
*/
#include <stdio.h>
#include <math.h>

#define MATRIX_SIZE 5
#define ERROR_THRESHOLD 1E-7
#define MY_SUCCESS 1
#define MY_FAIL -1

// get the sum of absolute value of the matrix (A - diag(A))
double off(double A[][MATRIX_SIZE]) {
    double sum = 0;
    for (int i = 0; i < MATRIX_SIZE; i++) {
        for (int j = 0; j < MATRIX_SIZE; j++) {
            if (i == j) {
                continue;
            }
            sum += fabs(A[i][j]);
        }
    }
    return sum;
}

/**
 * @brief get the index (p,q) of the matrix elements whose norm is the maximum s.t. p!=q
 * @param A input matrix
 * @param row pointer of row index
 * @param col pointer of col index
 */
double getMaxElemIndex(double A[][MATRIX_SIZE], int *row, int *col) {
    //in case of A[p][q]==0 for all p!=q
    double temp = -1;
    int p = -1;
    int q = -1;
    for (int i = 0; i < MATRIX_SIZE; i++) {
        for (int j = 0; j < MATRIX_SIZE; j++) {
            if (i == j) {
                continue;
            }
            if (temp < fabs(A[i][j])) {
                p = i;
                q = j;
                temp = fabs(A[i][j]);
            }
        }
    }
    *row = p;
    *col = q;
}

/**
 * @brief initialize Givens matrix Q
 * @param Q output, the Givens matrix Q
 * @param p row index
 * @param q col index
 * @param cos_theta the value of cos(theta)
 * @param sin_theta the value of sin(theta)
 */
void initQMatrix(double Q[][MATRIX_SIZE], int p, int q, double cos_theta, double sin_theta) {
    for (int i = 0; i < MATRIX_SIZE; i++) {
        for (int j = 0; j < MATRIX_SIZE; j++) {
            Q[i][j] = 0;
        }
        Q[i][i] = 1;
    }
    Q[p][p] = cos_theta;
    Q[q][q] = cos_theta;
    Q[p][q] = sin_theta;
    Q[q][p] = -sin_theta;
}

/**
 * @brief Calculate matrix multiplication
 * @param A The input matrix
 * @param B The input matrix
 * @param res output res=AB
 */
void matrixMultiply(double A[][MATRIX_SIZE], double B[][MATRIX_SIZE], double res[][MATRIX_SIZE]) {
    double C[MATRIX_SIZE][MATRIX_SIZE] = {0};
    for (int i = 0; i < MATRIX_SIZE; i++) {
        for (int j = 0; j < MATRIX_SIZE; j++) {
            double temp = 0;
            for (int k = 0; k < MATRIX_SIZE; k++) {
                temp += A[i][k] * B[k][j];
            }
            C[i][j] = temp;
        }
    }//end for
    for (int i = 0; i < MATRIX_SIZE; i++) {
        for (int j = 0; j < MATRIX_SIZE; j++) {
            res[i][j] = C[i][j];
        }
    }
}

/**
 * @brief Transpose @param Q, output the result to @param res
 * @param Q The input matrix
 * @param res output: res=Q'
 */
void matrixTranspose(double Q[][MATRIX_SIZE], double res[][MATRIX_SIZE]) {
    for (int i = 0; i < MATRIX_SIZE; i++) {
        for (int j = 0; j < MATRIX_SIZE; j++) {
            res[i][j] = Q[j][i];
        }
    }
}

/**
 * @brief Get eigenvalues and eigenvectors of matrix by Jacobi method
 * @param A :The input matrix
 * @param eigenvalues :output, eigenvalues of A
 * @param  V output, eigenvectors corresponding to the eigenvalues
 * @retval MY_SUCCESS:The function runs successfully; MY_FAIL: An error occurred and exit
 */
int JacobiMethod(double A[][MATRIX_SIZE], double eigenvalues[MATRIX_SIZE], double V[][MATRIX_SIZE]) {
    //检验A是否为对称阵
    for (int i = 0; i < MATRIX_SIZE; i++) {
        for (int j = i + 1; j < MATRIX_SIZE; j++) {
            if (A[i][j] != A[j][i]) {
                return MY_FAIL;
            }
        }
    }
    //初始化V为单位阵I
    for (int i = 0; i < MATRIX_SIZE; i++) {
        for (int j = 0; j < MATRIX_SIZE; j++) {
            V[i][j] = 0;
        }
        V[i][i] = 1;
    }
    //Jacobi迭代
    double c = 0;
    double d = 0;
    double t = 0;
    while (1) {
        int p = 0;
        int q = 0;
        getMaxElemIndex(A, &p, &q);// p!=q
        double s = 0;
        if (fabs(A[p][q]) > 1E-14) {//A[p][q]!=0 保证A为对角阵也能输出正确结果
            s = (A[q][q] - A[p][p]) / (2 * A[p][q]);
            if (s >= 0) {
                t = 1 / (s + sqrt(1 + s * s));
            } else {
                t = 1 / (s - sqrt(1 + s * s));
            }//end if
            c = 1 / sqrt(1 + t * t);
            d = t * c;
        } else {//A[p][q]==0
            c = 1;
            d = 0;
        }//end if
        //实际上矩阵乘法只涉及A[i][p],A[i][q]可以优化,这里简单起见直接暴力
        double Q[MATRIX_SIZE][MATRIX_SIZE] = {0};
        initQMatrix(Q, p, q, c, d);
        double Q_transpose[MATRIX_SIZE][MATRIX_SIZE] = {0};
        matrixTranspose(Q, Q_transpose);
        matrixMultiply(Q_transpose, A, A);// A<-Q'A
        matrixMultiply(A, Q, A);// A<-AQ
        matrixMultiply(V, Q, V);// V<-VQ
        if (off(A) < ERROR_THRESHOLD) {//检验迭代是否结束
            break;
        }//end if
    }//end for
    for (int i = 0; i < MATRIX_SIZE; i++) {//A已经近似为对角阵，输出对角线元素为特征值
        eigenvalues[i] = A[i][i];
    }
    return MY_SUCCESS;
}

int main() {
    double A[][MATRIX_SIZE] = {
            {3, 2,  5,  4,  6},
            {2, 1,  3,  -7, 8},
            {5, 3,  2,  5,  -4},
            {4, -7, 5,  1,  3},
            {6, 8,  -4, 3,  8}
    };
    double r[MATRIX_SIZE] = {0};
    double v[MATRIX_SIZE][MATRIX_SIZE] = {0};
    JacobiMethod(A, r, v);
    for (int i = 0; i < MATRIX_SIZE; i++) {
        printf("r%d=%.15f\nv%d=(", i + 1, r[i], i + 1);
        for (int j = 0; j < MATRIX_SIZE; j++) {
            printf("%.15f,", v[j][i]);
        }
        printf("\b)\n");
    }
    return 0;
}