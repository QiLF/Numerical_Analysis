/**
 * HW5
 * nonlinear equation f(x) = 1/3 * x^3 - x
 * author: lf
 * res:
 */
#include <stdio.h>
#include <math.h>

#define ERROR_THRESHOLD 1E-8

typedef struct {
    double root;
    int recursive_num;
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
    Result res = {0, 1};
    double next = x0 - (x0 * x0 * x0 / 3 - x0) / (x0 * x0 - 1);
    while (fabs(next - x0) >= ERROR_THRESHOLD) {
        x0 = next;
        next = x0 - (x0 * x0 * x0 / 3 - x0) / (x0 * x0 - 1);
        res.recursive_num++;
    }
    res.root = next;
    return res;
}

/**
 * @brief Secant method for solving nonlinear equations
 * @param x0 the initial value for secant method
 * @param x1 the initial value for secant method
 * @retval the root of f(x).   type: ValsForSecant
 */
Result secantMethod(double x0, double x1) {
    Result res = {0, 1};
    double next = x1 - (x1 * x1 * x1 / 3 - x1)*(x1-x0)/((x1 * x1 * x1 / 3 - x1)-(x0 * x0 * x0 / 3 - x0));
    while (fabs(next - x1) >= ERROR_THRESHOLD) {
        x0 = x1;
        x1 = next;
        next = x1 - (x1 * x1 * x1 / 3 - x1)*(x1-x0)/((x1 * x1 * x1 / 3 - x1)-(x0 * x0 * x0 / 3 - x0));
        res.recursive_num++;
    }
    res.root = next;
    return res;
}


int main() {
    puts("Newton method result:");
    for (int i = 0; i < 4; i++) {
        Result res = newtonMethod(initValsForNetween[i]);
        printf("input:x0=%.15f  res: root=%.15f, recursive_num=%d\n", initValsForNetween[i], res.root,
               res.recursive_num);
    }
    puts("Scant method result:");
    for (int i = 0; i < 4; i++) {
        Result res = secantMethod(initValsForaSecant[i][0], initValsForaSecant[i][1]);
        printf("x0=%.15f,x1=%.15f res: root=%.15f, recursive_num=%d\"\n", initValsForaSecant[i][0],
               initValsForaSecant[i][1], res.root, res.recursive_num);
    }
    return 0;
}