/* HW1
 * 
 * 2018-9-11
 * input:0.000000000000000 res:1.644933566848396e+000
 * input:0.500000000000000 res:1.227410944427013e+000
 * input:1.000000000000000 res:9.999995000002929e-001
 * input:1.414213562373095 res:8.749824102346320e-001
 * input:10.000000000000000 res:2.928959163104470e-001
 * input:100.000000000000000 res:5.187278512688917e-002
 * input:300.000000000000000 res:2.094121640573054e-002
*/
#include <stdio.h>
#include <math.h>

#define EQUAL ==

/* estimate the number k so that the error of sum can be less than 1e-6 */
unsigned long estimateK(double x) {
    unsigned long k;
    if (x EQUAL 0) {    //  x=0 need to be considered separately
        k = 2e6;
    } else {
        k = (1 + 1 / x) * 1e6;
    }
    return k;
}

/* calculate the sum of the series with error less than 1e-6 */
double getSum(double x) {
    double sum = 0;
    unsigned long k = estimateK(x);
    for (int i = 1; i <= k; i++) {
        double temp = i * (i + x);
        sum += 1 / temp;
    }
    return sum;
}

int main() {
    // inputs
    double x_array[] = {0.0, 0.5, 1.0, 0, 10.0, 100.0, 300.0};
    x_array[3] = sqrt(2);
    for (int i = 0; i < 7; i++) {
        printf("input:%.15f res:%.15e\n", x_array[i], getSum(x_array[i]));
    }
    return 0;
}
