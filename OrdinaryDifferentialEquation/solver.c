/**
 * HW8
 * 2018-12-28
 * author: qlf
 * Result
    Runge-Kutta Result:
    h=0.100000,err=0.000017800000795,ok=4.094116865347997
    h=0.050000,err=0.000001042240726,ok=4.050319982837565
    h=0.025000,err=0.000000062907182,ok=4.025859479718573
    h=0.012500,err=0.000000003861853
    Adams Result:
    h=0.100000,err=0.000081825365899,ok=4.342281405160378
    h=0.050000,err=0.000004033954575,ok=4.370487209208222
    h=0.025000,err=0.000000195021844,ok=4.239952930241303
    h=0.012500,err=0.000000010321201
 */
#include <stdio.h>
#include <math.h>
#include <malloc.h>

#define HIGH 1.5
#define LOW 0
#define ERROR_THRESHOLD 1e-7

//get the value of f(x,y) at (x,y)
double fxy(double x, double y) {
    return (-x * x * y * y);
}

/**
 * @brief solve dy/dx = f(x,y) by 4-order Runge-Kutta method
 * @param y0 initial val of y
 * @param h step length
 * @param y a double array: vals of y on [0,1.5]
 */
void solveByRungeKutta(double y0, double h, double y[]) {
    double xi = LOW;
    y[0] = y0;
    int steps = (int) ((HIGH - LOW) / h + 0.5) + 1;
    for (int i = 1; i < steps; i++) {
        double k1 = fxy(xi, y[i - 1]);
        double k2 = fxy(xi + h / 2, y[i - 1] + h * k1 / 2);
        double k3 = fxy(xi + h / 2, y[i - 1] + h * k2 / 2);
        double k4 = fxy(xi + h, y[i - 1] + h * k3);
        xi += h;
        y[i] = y[i - 1] + (k1 + 2 * k2 + 2 * k3 + k4) * h / 6;
    }
}

/**
 * @brief solve dy/dx = f(x,y) by 4-order Adams method
 * @param y0 initial val of y
 * @param h step length
 * @param y a double array: vals of y on [0,1.5]
 */
void solveByAdams(double y0, double h, double x[], double y[]) {
    //三阶Runge-Kutta公式计算起始值 y1,y2
    int steps = (int) ((HIGH - LOW) / h + 0.5) + 1;
    x[0] = LOW;
    y[0] = y0;
    for (int i = 1; i < steps; i++) {
        x[i] = x[i - 1] + h;
    }
    for (int i = 1; i < 4; i++) {
        double k1 = fxy(x[i - 1], y[i - 1]);
        double k2 = fxy(x[i - 1] + h / 2, y[i - 1] + h * k1 / 2);
        double k3 = fxy(x[i - 1] + h / 2, y[i - 1] + h * k2 / 2);
        double k4 = fxy(x[i - 1] + h, y[i - 1] + h * k3);
        y[i] = y[i - 1] + (k1 + 2 * k2 + 2 * k3 + k4) * h / 6;
    }
    //Admas求解
    for (int i = 4; i < steps; i++) {
        //预估
        y[i] = y[i - 1] + h / 24 * (55 * fxy(x[i - 1], y[i - 1]) - 59 * fxy(x[i - 2], y[i - 2]) +
                                    37 * fxy(x[i - 3], y[i - 3]) - 9 * fxy(x[i - 4], y[i - 4]));
        //校正
        y[i] = y[i - 1] + h / 24 *
                          (9 * fxy(x[i], y[i]) + 19 * fxy(x[i - 1], y[i - 1]) - 5 * fxy(x[i - 2], y[i - 2]) +
                           fxy(x[i - 3], y[i - 3]));
    }
}

// get the real value of y(1.5)
double realRes(double x) {
    return 3 / (1 + pow(x, 3));
}


int main() {
    double error[4] = {0};
    double h[4] = {0};
    //Runge-Kutta
    for (int i = 0; i < 4; i++) {
        h[i] = 0.1 / pow(2, i);
        double y0 = 3;
        int length = (int) ((HIGH - LOW) / h[i] + 0.5) + 1;
        double *y = malloc(sizeof(double) * length);
        solveByRungeKutta(y0, h[i], y);
        error[i] = fabs(realRes(1.5) - y[length - 1]);
        free(y);
    }
    printf("Runge-Kutta Result:\n");
    for (int i = 0; i < 3; i++) {
        printf("h=%f,err=%.15f,ok=%.15f\n", h[i], error[i], log(error[i] / error[i + 1]) / log(2));
    }
    printf("h=%f,err=%.15f\n", h[3], error[3]);
    //Adams
    for (int i = 0; i < 4; i++) {
        h[i] = 0.1 / pow(2, i);
        double y0 = 3;
        int length = (int) ((HIGH - LOW) / h[i] + 0.5) + 1;
        double *x = malloc(sizeof(double) * length);
        double *y = malloc(sizeof(double) * length);
        solveByAdams(y0, h[i], x, y);
        double errori = fabs(realRes(1.5) - y[length - 1]);
        error[i] = errori;
        free(y);
        free(x);
    }
    printf("Adams Result:\n");
    for (int i = 0; i < 3; i++) {
        printf("h=%f,err=%.15f,ok=%.15f\n", h[i], error[i], log(error[i] / error[i + 1]) / log(2));
    }
    printf("h=%f,err=%.15f\n", h[3], error[3]);
    return 0;
}