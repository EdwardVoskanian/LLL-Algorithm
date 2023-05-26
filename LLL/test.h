using namespace std;

void getNumQuotients(int &n) {
    
    fstream file;
    int input;
    string fileName = "input.txt";
    
    file.open(fileName.c_str());
    n = -1;
    while (file >> input) {
        n = n + 1;
    }
    file.close();
    
    return;
}

void getDelta(int &n) {
    
    fstream file;
    int input;
    string fileName = "input.txt";
    
    file.open(fileName.c_str());
    n = -1;
    while (file >> input) {
        n = n + 1;
    }
    file.close();
    
    return;
}

void getBase(double &base) {

    fstream file;
    string fileName = "input.txt";

    file.open(fileName.c_str());
    file >> base;
    file.close();

    return;
}

void getQuotients(double * quotients,int base) {

    int i = 0;
    fstream file;
    double input;
    string fileName = "input.txt";
    
    file.open(fileName.c_str());
    //cout << endl;
    while (file >> input) {
        if (i > 0) {
            quotients[i - 1] = log(input)/log(base);
            i = i + 1;
        }
        else {
            i = i + 1;
        }
    }
    file.close();
    quotients[0] = sqrt(2);
    quotients[1] = sqrt(3);
    quotients[1] = sqrt(5);
    quotients[1] = sqrt(7);
    
    return;
}

void getConvergent(double alpha, int m, int &a, int &b) {

    int i;
    double p[m + 2];
    double q[m + 2];
    double A[m + 2];
    char answer;
    
    p[0] = 0; q[0] = 1;
    p[1] = 1; q[1] = 0;

    for (i = 2; i < m + 2; i++) {
        A[i] = lrint(floor(alpha));
        p[i] = A[i] * p[i - 1] + p[i - 2];
        q[i] = A[i] * q[i - 1] + q[i - 2];
        if (fabs(alpha - A[i]) < 1.19209e-07) {
            a = p[i];
            b = q[i];
            return;
        }
        alpha = 1.0 / (alpha - A[i]);
    }
    a = p[i - 1];
    b = q[i - 1];
    
    return;
}

void updateConvergents(int n, double * quotients, int ** convergents, int ** ticket, int MAX) {
    
    int i, j = 0;
    int a, b;
    char answer;
    
    for (i = 0; i < n; i++) {
        if (ticket[0][i] == 1 && ticket[1][i] < MAX + 1) {
            ticket[1][i] = ticket[1][i] + 1;
            getConvergent(quotients[i], ticket[1][i], a, b);
            convergents[0][i] = a;
            convergents[1][i] = b;
            ticket[0][i] = 0;
        }
        else {
            if (ticket[1][i] == MAX) exit(0);
        }
    }
//    cout << endl;
//    for (j = 0; j < 2 * n; j++) {
//        if (j < n) {
//            cout << ticket[1][j] << " ";
//        }
//        else {
//            cout << " ";
//        }
//    }
//    cout << endl;
//    for (j = 0; j < n; j++) {
//        cout << convergents[0][j] << "/" << convergents[1][j] << "  ";
//    }
//    cout << endl;
	
    return;
}

void gramSchmidt(int n, double ** Y, double ** Y_star, double ** M, double * gamma_star) {

    int i, j, k = 0;
    double sum = 0.0;

    for (i = 0; i < n + 1; i++) {
        for (j = 0; j < n + 1; j++) {
            Y_star[i][j] = Y[i][j];
        }
        for (j = 0; j < i; j++) {
            sum = 0.0;
            for (k = 0; k < n + 1; k++) {
                sum = sum + Y[i][k] * Y_star[j][k];
            }
            M[i][j] = sum / gamma_star[j];
            for (k = 0; k < n + 1; k++) {
                Y_star[i][k] = Y_star[i][k] - M[i][j] * Y_star[j][k];
            }
        }
        for (j = 0; j < n + 1; j++) {
            gamma_star[i] = gamma_star[i] + Y_star[i][j] * Y_star[i][j];
        }
    }

    return;
}

void reduce(int n, int k, int l, double ** Y, double ** M, double ** C) {

    int i, j = 0;

    if (fabs(M[k][l]) >= 0.5) {
        for (i = 0; i < n + 1; i++) {
            Y[k][i] = Y[k][i] - (round(M[k][l]) * Y[l][i]);
            C[k][i] = C[k][i] - (round(M[k][l]) * C[l][i]);
        }
        for (j = 0; j < l; j++) {
            M[k][j] = M[k][j] - (round(M[k][l]) * M[l][j]);
        }
        M[k][l] = M[k][l] - round(M[k][l]);
    }

    return;
}

void exchange(int n, int k, double ** Y, double ** M, double * gamma_star, double ** C) {

    int i = 0;
    double temp = 0.0;
    double temp_1 = 0.0;
    double temp_2 = 0.0;
    double nu = 0.0;
    double delta = 0.0;
    double xsi = 0.0;

    /* exchange y_k-1 and y_k */
    for (i = 0; i < n + 1; i++) {
        temp_1 = Y[k - 1][i];
        temp_2 = C[k - 1][i];
        Y[k - 1][i] = Y[k][i];
        C[k - 1][i] = C[k][i];
        Y[k][i] = temp_1;
        C[k][i] = temp_2;
    }
    nu = M[k][k - 1];
    delta = gamma_star[k] + (nu * nu) * gamma_star[k - 1];
    M[k][k - 1] = nu * gamma_star[k - 1] / delta;
    gamma_star[k] = gamma_star[k] * gamma_star[k - 1] / delta;
    gamma_star[k - 1] = delta;
    /* exchange mu(k-1,t) and mu(k,t) */
    for (i = 0; i < k - 1; i++) {
        temp = M[k - 1][i];
        M[k - 1][i] = M[k][i];
        M[k][i] = temp;
    }
    for (i = k + 1; i < n + 1; i++) {
        xsi = M[i][k];
        M[i][k] = M[i][k - 1] - (nu * M[i][k]);
        M[i][k - 1] = (M[k][k - 1] * M[i][k]) + xsi;
    }

    return;
}

void getC(double ALPHA, int n, double delta, int ** convergents, double ** C) {

    int i, j = 0;
    
    double ** X = new double * [n + 1];
        for (i = 0; i < n + 1; i++) {
            X[i] = new double[n + 1];
            for (j = 0; j < n + 1; j++) {
                X[i][j] = 0.0;
            }
        }
        for (i = 1; i < n + 1; i++) {
            X[0][i] = 0.0;
        }
        for (i = 1; i < n + 1; i++) {
            X[i][i] = -1.0;
        }
    double ** Y = new double * [n + 1];
        for (i = 0; i < n + 1; i++) {
            Y[i] = new double[n + 1];
            for (j = 0; j < n + 1; j++) {
                Y[i][j] = 0.0;
            }
        }
    double ** Y_star = new double * [n + 1];
        for (i = 0; i < n + 1; i++) {
            Y_star[i] = new double[n + 1];
            for (j = 0; j < n + 1; j++) {
                Y[i][j] = 0.0;
            }
        }
    double ** M = new double * [n + 1];
        for (i = 0; i < n + 1; i++) {
            M[i] = new double[n + 1];
            for (j = 0; j < n + 1; j++) {
                M[i][j] = 0.0;
            }
        }
        for (i = 0; i < n + 1; i++) {
            M[i][i] = 1.0;
        }
    double * gamma_star = new double[n + 1];
        for (i = 0; i < n + 1; i++) {
            gamma_star[i] = 0.0;
        }
    
    // INTIALIZE <X>
    X[0][0] = pow(2.0, -n * (n + 1.0) / 4.0) * pow(delta, n + 1);
    for (i = 0; i < n; i++) {
        X[0][i + 1] = (double)convergents[0][i] / convergents[1][i];
    }
  
    // STEP 1: MAKE A COPY OF <X>
    for (i = 0; i < n + 1; i++) {
        for (j = 0; j < n + 1; j++) {
            Y[i][j] = X[i][j];
        }
    }

    // STEP 2: COMPUTE THE GRAM-SCHMIDT ORTHOGONALIZATION
    gramSchmidt(n, Y, Y_star, M, gamma_star);

    // STEP 3: PERFORM THE BASIS REDUCTION
    int k = 1;
    while (k < n + 1) {
        reduce(n, k, k - 1, Y, M, C);
        if (gamma_star[k] >= (ALPHA - pow(M[k][k - 1], 2)) * gamma_star[k - 1]) {
            for (int l = k - 2; l >= 0; l--) {
                reduce(n, k, l, Y, M, C);
            }
            k = k + 1;
        }
        else {
            exchange(n, k, Y, M, gamma_star, C);
            if (k > 1) {
                k = k - 1;
            }
        }
    }

    return;
}

void check(int n, double ** C, double * quotients, int ** convergents, int ** ticket) {

    double q = abs(C[0][0]);
    double k;
    double error_1; // positive distance between the real number and the convergent
    double error_2; // positive distance between simultaneous Diophanitine approximation and the convergent
    int i, j = 0;

    for (i = 0; i < n; i++) {
        k = abs(C[0][i + 1]);
        error_1 = fabs(quotients[i] - (double)convergents[0][i] / convergents[1][i]);
        error_2 = fabs(k/q - (double)convergents[0][i] / convergents[1][i]);
        if (error_1 / error_2 > .5) {
            ticket[0][i] = ticket[0][i] + 1;
        }
    }
    
    return;
}

void results(int n, double ** C, double * quotients, int ** convergents, int &temp) {

    int i, j = 0;
    double q = abs(C[0][0]);
        if (q == temp) {
            return;
        }
        temp = q;
    double k;
    double error_1; // positive distance between the real number and the convergent
    double error_2; // positive distance between simultaneous Diophanitine approximation and the convergent
    
    cout << endl;
    cout << "A new approximation:" << endl << endl;
    for (i = 0; i < n; i++) {
        cout << "a_" << i + 1 << " / b_" << i + 1 << " = " << convergents[0][i] << " / " << convergents[1][i] << endl;
    }
    cout << endl;
    for (i = 0; i < n; i++) {
        k = abs(C[0][i + 1]);
        error_1 = fabs(k/q - (double)convergents[0][i] / convergents[1][i]);
        error_2 = fabs(k/q - quotients[i]);
        cout << "k_" << i + 1 << " / q = " << k << " / " << q << "    " << error_1 << "  " << error_2 << "  " << error_1 / error_2 << endl;
    }
//    cout << endl;
    
    return;
}

void getApproximation(double &delta,double stepSize,double ALPHA,int MAX,int &temp,int &done) {
    
    int n;
    int i, j = 0;
    int S = 0;
    char answer;
    double base;
        getNumQuotients(n);
        getBase(base);
    double * quotients = new double[n];
    getQuotients(quotients, base);
//            for (i = 0; i < n; i ++) {
//                cout << quotients[i] << endl;
//            }
    int ** convergents = new int * [2];
        for (i = 0; i < 2; i++) {
            convergents[i] = new int[n];
        }
    int ** ticket = new int * [2];
        for (i = 0; i < 2; i++) {
            ticket[i] = new int[n];
        }
        for (i = 0; i < n; i++) {
            ticket[0][i] = 1; // 1 means get the next convergent
            ticket[1][i] = 0;
        }
    double ** C = new double * [n + 1];
        for (i = 0; i < n + 1; i++) {
            C[i] = new double[n + 1];
        }
    
    // PROGRAM STARTS HERE
//    cout << endl;
//    cout << "Starting ticket:" << endl;
//    for (i = 0; i < 2; i++) {
//        if(i == 1) {
//            for (j = 0; j < n; j++) {
//                cout << ticket[i][j] << " ";
//            }
//            cout << endl;
//        }
//        else {
//            for (j = 0; j < n; j++) {
//                cout << ticket[i][j] << " ";
//            }
//            cout << endl;
//        }
//    }
//    cout << "Would you like to continue: [y/n] ";
//    cin >> answer;
//    if (answer != 'y') exit (0);
    do {
        updateConvergents(n,quotients,convergents,ticket,MAX);
        for (j = 0; j < n; j++) {
            if (ticket[1][j] > MAX - 1) {
                cout << "Increase MAX to compute more convergents." << endl << endl;
                exit(0);
            }
        }
//        cout << "Would you like to continue: [y/n] ";
//        cin >> answer;
//        if (answer != 'y') exit (0);
        for (i = 0; i < n + 1; i++) {
            for (j = 0; j < n + 1; j++) {
                C[i][j] = 0.0;
            }
        }
        for (i = 0; i < n + 1; i++) {
            C[i][i] = 1.0;
        }
        getC(ALPHA,n,delta,convergents,C);  // Update <C> with the LLL algorithm
        check(n,C,quotients,convergents,ticket); // Update the ticket if the convergents aren't good enough
        S = 0;
        for (i = 0; i < n; i++) {
            S = S + ticket[0][i];
        }
    } while (S != 0);
    
    results(n, C, quotients, convergents, temp);
    delta = delta - stepSize;
    
    return;
}
