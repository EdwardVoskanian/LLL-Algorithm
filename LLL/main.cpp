#include <iostream>
#include <fstream>
#include <math.h>
#include <string>
#include <iomanip>
#include "test.h"

const double ALPHA = 0.75; // Reduction parameter
const double MIN = 0.0000000000001;
const int MAX = 20;

int main(int argc, char* argv[]) {
    
    double delta = atof(argv[1]);
    int nSteps = atoi(argv[2]);
    double stepSize = (delta - MIN)/nSteps;
    int temp = 0;
    int done = 0;
    
    cout << endl;
    cout << "range: (" << MIN << "," << delta << "]" << endl;
    cout << "stepsize: " << stepSize << endl;
    while (delta > MIN) {
        getApproximation(delta,stepSize,ALPHA,MAX,temp,done);
    }
    
    return 0;
}
