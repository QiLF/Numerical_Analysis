/*
 * HW3
 * PB16000953 祁琳峰
 * 2018-10-30
 */

#include <stdio.h>
#include <math.h>
#include <malloc.h>

typedef struct {
    double value;
    double error;
    double error_order;
}Res;


/**
 * @param1 power: determine the num of interpolation points;  #points = 2^power + 1
 * @return: integration value
 */
double simpson_integration(int power) {
    double step = (5 - 1) / pow(2, power);
    int size = (int) pow(2, power) + 1;
    double *values = malloc(sizeof(double) * size);
    for (int i = 0; i < size; i++) {
        values[i] = sin(1 + i * step);
    }
    double res = values[0];
    for (int i = 1; i < size - 1; i += 2) {
        res += 4 * values[i] + 2 * values[i + 1];
    }
    res += values[size - 1];
    res = res * step / 3;
    free(values);
    return res;
}

/**
 * @param1 power: determine the num of interpolation points;  #points = 2^power + 1
 * @return: integration value
 */
double trapezoidal_integration(int power) {
    double step = (5 - 1) / pow(2, power);
    int size = (int) pow(2, power) + 1;
    double *values = malloc(sizeof(double) * size);
    for (int i = 0; i < size; i++) {
        values[i] = sin(1 + i * step);
    }
    double res = 0.5 * values[0];
    for (int i = 1; i < size - 1; i++) {
        res += values[i];
    }
    res += 0.5 * values[size - 1];
    res *= step;
    free(values);
    return res;
}

/**
 * @param val: the result of numerical integration we calculate
 * @return: the error between numerical integration and the real result
 */
double get_error(double val) {
    double real_val = -cos(5) + cos(1);
    return fabs(real_val - val);
}

/*
 * @param res*: pointer of Res
 * @param n: length of results - 1
 * @return: void
 */
void calc_error_order(Res *res,int n)
{
    for(int i=0;i<n;i++)
    {
        res[i].error_order=log(res[i].error/res[i+1].error)/log(2);
    }
}

int main() {
    printf("simpson integration result:\n");
    Res results[13];
    for (int i = 1; i <= 13; i++) {
        results[i-1].value = simpson_integration(i);
        results[i-1].error=get_error(results[i-1].value);
    }
    calc_error_order(results,12);
    for(int i=1;i<=12;i++)
    {
        printf("input:%d res:%.15f error:%.15f error_order:%f \n", i, results[i-1].value, results[i-1].error,
                results[i-1].error_order);
    }
    printf("trapezoidal integration result:\n");
    for (int i = 1; i <= 13; i++) {
        results[i-1].value = trapezoidal_integration(i);
        results[i-1].error=get_error(results[i-1].value);
    }
    calc_error_order(results,12);
    for(int i=1;i<=12;i++)
    {
        printf("input:%d res:%.15f error:%.15f error_order:%f \n", i, results[i-1].value, results[i-1].error,
               results[i-1].error_order);
    }
    return 0;
}