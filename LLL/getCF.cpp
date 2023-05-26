#include <stdio.h>
#include <iostream>
#include <fstream>
#include <math.h>

using namespace std;

const double E = 2.718281828459045235360287471352662497757247093699959574966;
const double PI = 3.141592653589793238462643383279502884197169399375105820974;
const double GOLDENRATIO = (1.0 + sqrt(5))/2.0;

void find_cf(double x, int y);

int main() {
    
    find_cf(log2(3),35);
    
    return 0;
}

void find_cf(double x, int y) {
    
    double realNumber = x;
    int MAX = y;
    long p[MAX];
    long q[MAX];
    long a[MAX];
    int i = 0;
    
    // The first two convergents are 0/1 and 1/0
    p[0] = 0; q[0] = 1;
    p[1] = 1; q[1] = 0;
    
    // The rest of the convergents (and continued fractions)
    for (i = 2; i < MAX + 2; i++) {
        a[i] = lrint(floor(x));
        p[i] = a[i] * p[i - 1] + p[i - 2];
        q[i] = a[i] * q[i - 1] + q[i - 2];
//        cout << a[i] << ":  " << p[i] << " " << q[i] << "   " << fabs((double)p[i]/q[i] - realNumber) << endl;
        if (fabs(x - a[i]) < 1.19209e-07) {
            return;
        }
        x = 1.0 / (x - a[i]);
        cout << a[i] << ":  " << p[i] << " " << q[i] << endl;
    }
    
}
