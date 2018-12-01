/**
 * HW5
 * nonlinear equation f(x) = 1/3 * x^3 - x
 * author: qlf
 * res:
    Newton method result:
    input:x0=0.100000000000000  res: root=-0.000000000000000, recursive_num=3,convergenceDegree=3.055519
    input:x0=0.200000000000000  res: root=0.000000000000000, recursive_num=4,convergenceDegree=3.025366
    input:x0=0.900000000000000  res: root=-1.732050807568877, recursive_num=7,convergenceDegree=2.011435
    input:x0=9.000000000000000  res: root=1.732050807568877, recursive_num=10,convergenceDegree=1.966333
    Scant method result:
    x0=-0.100000000000000,x1=0.100000000000000 res: root=-0.000000000000000, recursive_num=1,convergenceDegree=1.000000
    x0=-0.200000000000000,x1=0.200000000000000 res: root=0.000000000000000, recursive_num=1,convergenceDegree=1.000000
    x0=-2.000000000000000,x1=0.900000000000000 res: root=1.732050807753635, recursive_num=11,convergenceDegree=1.588367
    x0=0.900000000000000,x1=9.000000000000000 res: root=1.732050807453568, recursive_num=13,convergenceDegree=1.626905
 */
#include <stdio.h>
#include <math.h>

#define ERROR_THRESHOLD 1E-8

typedef struct {
    double root;
    int recursive_num;
    double convergenceDegree;
} Result;

double initValsForNetween[] = {0.1, 0.2, 0.9, 9.0};
double initValsForaSecant[][2] = {{-0.1, 0.1},
                                  {-0.2, 0.2},
                                  {-2.0, 0.9},
                                  {0.9,  9.0}};

/**
 * @brief Newton method for solving nonlinear equations
 * @param x0 the initial value for newton method
 * @retval the root of f(x)
 */
Result newtonMethod(double x0) {
    Result res = {0, 1, 0};
    double x1 = x0;
    double next = x1 - (x1 * x1 * x1 / 3 - x1) / (x1 * x1 - 1);
    while (fabs(next - x1) >= ERROR_THRESHOLD) {
        x0 = x1;
        x1 = next;
        next = x1 - (x1 * x1 * x1 / 3 - x1) / (x1 * x1 - 1);
        res.recursive_num++;
    }
    res.root = next;
    res.convergenceDegree = log(fabs(x1 - res.root)) / log(fabs(x0 - res.root));
    return res;
}

/**
 * @brief Secant method for solving nonlinear equations
 * @param x0 the initial value for secant method
 * @param x1 the initial value for secant method
 * @retval the root of f(x).   type: ValsForSecant
 */
Result secantMethod(double x0, double x1) {
    Result res = {0, 1, 0};
    double fx1 = x1 * x1 * x1 / 3 - x1;
    double fx0 = x0 * x0 * x0 / 3 - x0;
    double next = x1 - fx1 * (x1 - x0) / (fx1 - fx0);
    while (fabs(next - x1) >= ERROR_THRESHOLD && fabs(next * next * next / 3 - next) >= ERROR_THRESHOLD) {
        x0 = x1;
        fx0 = fx1;
        x1 = next;
        fx1 = x1 * x1 * x1 / 3 - x1;
        next = x1 - fx1 * (x1 - x0) / (fx1 - fx0);
        res.recursive_num++;
    }
    res.root = next;
    res.convergenceDegree = log(fabs(x1 - res.root)) / log(fabs(x0 - res.root));
    return res;
}


int main() {
    puts("Newton method result:");
    for (int i = 0; i < 4; i++) {
        Result res = newtonMethod(initValsForNetween[i]);
        printf("input:x0=%.15f  res: root=%.15f, recursive_num=%d,convergenceDegree=%f\n",
               initValsForNetween[i], res.root, res.recursive_num, res.convergenceDegree);
    }
    puts("Scant method result:");
    for (int i = 0; i < 4; i++) {
        Result res = secantMethod(initValsForaSecant[i][0], initValsForaSecant[i][1]);
        printf("x0=%.15f,x1=%.15f res: root=%.15f, recursive_num=%d,convergenceDegree=%f\n",
               initValsForaSecant[i][0], initValsForaSecant[i][1], res.root, res.recursive_num, res.convergenceDegree);
    }
    return 0;
}