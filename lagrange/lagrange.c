/* HW2
 * 
 * 2018-9-25
 * RESULT:
   interpolation points: xi=-5+10i/n  where i=1,2...,n
   n:  5, maxError:4.326923076923077e-001
   n: 10, maxError:1.915643050219250e+000
   n: 20, maxError:5.976568477453244e+001
   n: 40, maxError:1.039408117697524e+005

   interpolation points: xi=-5cos[(2i+1)PI/(2n+2)] where i=1,2,...,n
   n:  5,maxError:5.559113388123953e-001
   n: 10,maxError:1.091467246497665e-001
   n: 20,maxError:1.532508854382753e-002
   n: 40,maxError:2.889123107672198e-004
 */
#include <stdio.h>
#include <malloc.h>
#include <math.h>

#define MY_FAIL 0
#define EQUAL ==
#define MODE_ONE 1
#define MODE_TWO 2
int nList[] = {5, 10, 20, 40};  // a list of Num of interpolation points

/* Create a list of x by Num of interpolation points and Mode
 * n   : size of the list
 * mode: MODE_ONE is xi = -5+10i/n  where i = 1,2...,n
 *       MODE_TWO is xi = -5cos[(2i+1)PI/(2n+2)] where i = 1,2,...,n
 */
double *create_xList(int n, int mode) {
    if (n <= 0)     // in case of n < 0
    {
        printf("Error: n should be a positive num");
        return MY_FAIL;
    }
    double *xList = (double *) malloc((n + 1) * sizeof(double));
    if (mode EQUAL MODE_ONE)     // MODE1 -5+10i/n
    {
        for (int i = 0; i <= n; i++) {
            xList[i] = -5.0 + 10.0 / n * i;
        }
    }
    else if (mode EQUAL MODE_TWO)    // MODE2 -5cos[(2i+1)PI/(2n+2)]
    {
        for (int i = 0; i <= n; i++) {
            xList[i] = -5 * cos((2 * i + 1) * M_PI / (2 * n + 2));
        }
    }
    return xList;
}

/* Calculate the value of function f(x)=1/(1+x^2) by a list of x-coordinate
 * xList: the x-coordinate of interpolation points that we care about
 * n: the num of interpolation points
 */
double *create_fxList(double *xList, int n) {
    double *fxList = (double *) malloc((n + 1) * sizeof(double));
    for (int i = 0; i <= n; i++) {
        fxList[i] = 1 / (1 + xList[i] * xList[i]);
    }
    return fxList;
}

/* Create a list of y so that we can estimate the max error of Ln
 * yi=-5+10*i/500 where i=0,1,2,...,500
 * the size of this list is 501
 */
double *create_yList() {
    double *yList = (double *) malloc((500 + 1) * sizeof(double));
    for (int i = 0; i <= 500; i++) {
        yList[i] = 10.0 * i / 500.0 - 5;
    }
    return yList;
}

/* Calculate the value of base function
 * xList: the x-coordinate of interpolation points that we care about
 * i: the label of the base func  i belongs to {0,1,...,n}
 * n: num of interpolation points
 * x: the x-coordinate of the point we need to calculate
 */
double get_li(double *xList, int i, int n, double x) {
    double val = 1;
    for (int j = 0; j <= n; j++) {
        if (j EQUAL i) {
            continue;
        }
        val *= (x - xList[j]) / (xList[i] - xList[j]);
    }
    return val;
}

/* Calculate the value of Lagrange Function
 * x: the x-coordinate of the point we need to calculate
 * xList: the x-coordinate of interpolation points
 * n: the num of interpolation points
 */
double get_Lnx(double x, double *xList, int n) {
    double val = 0;
    double *fxList = create_fxList(xList, n);
    for (int i = 0; i <= n; i++) {
        val += fxList[i] * get_li(xList, i, n, x);
    }
    free(fxList);
    return val;
}

/* Estimate the max error
 * yList: the x-coordinate of test points
 * xList: the x-coordinate of interpolation points
 * n: the num of interpolation points
 */
double estimate_maxError(double *yList, double *xList, int n) {
    double max_error = 0;
    for (int i = 0; i <= 500; i++) {
        double f_yi = 1 / (1 + yList[i] * yList[i]);
        double Ln_yi = get_Lnx(yList[i], xList, n);
        double abs_error = fabs(f_yi - Ln_yi);
        if (max_error < abs_error) {
            max_error = abs_error;
        }
    }
    return max_error;
}

int main() {
    double *yList = create_yList();
    printf("interpolation points: xi=-5+10i/n  where i=1,2...,n\n");
    for (int i = 0; i < 4; i++) {
        double *xList = create_xList(nList[i], MODE_ONE);
        double maxError = estimate_maxError(yList, xList, nList[i]);
        printf("n: %2d, maxError:%.15e\n", nList[i], maxError);
        free(xList);
    }
    printf("\n");
    printf("interpolation points: xi=-5cos[(2i+1)PI/(2n+2)] where i=1,2,...,n\n");
    for (int i = 0; i < 4; i++) {
        double *xList = create_xList(nList[i], MODE_TWO);
        double maxError = estimate_maxError(yList, xList, nList[i]);
        printf("n: %2d,maxError:%.15e\n", nList[i], maxError);
        free(xList);
    }
    free(yList);
    return 0;
}
