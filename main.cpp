#include <iostream>
#include <cmath>
#include <vector>
#include <random>
#include <iomanip>
using namespace std;
uint32_t seed_val = 354782;
mt19937 rng(seed_val);
double f(double x){
    return (x*x*sin(x));
}
vector<double> getcolumn (vector<vector<double>> A, int k){ //get a column of number k as a vector from a matrix
    vector<double> result;
    for(int i = 0; i < A.size(); i++){
        result.push_back(A[i][k]);
    }
    return result;
}
vector<vector<double>> transpose (vector<vector<double>> A){ //returns the transposed matrix
    vector<vector<double>> result (A[0].size());
    for(int i = 0; i < result.size(); i++){
        result[i] = getcolumn(A,i);
    }
    return result;
}
vector<double> vnumm(double a, vector<double> b){ //multiply a vector by a number
    for (int i = 0; i < b.size(); i++){
        b[i] = b[i]*a;
    }
    return b;
}
vector <vector<double>> mmultiply (vector <vector<double>> A, vector <vector<double>> B){ //matrix multiplication
    if (A[0].size() == B.size()){
        vector <vector<double>> C (A.size());
        for (int i = 0; i < C.size(); i++){
            C[i].resize(B[0].size());
        }
        for (int i = 0; i < A.size(); i++){
            for (int j = 0; j < B[0].size(); j++){
                for (int k = 0; k < A[0].size(); k++){
                    C[i][j] += A[i][k]*B[k][j];
                }
            }
        }
        return C;
    }
}
pair<vector<vector<double>>,vector<vector<double>>> LUP (vector <vector<double>> A){ //returns a matrix that stores L,U and a matrix P (transposition matrix)
    vector <vector<double>> P(A.size());
    for(int i = 0; i < P.size(); i++){
        P[i].resize(A.size(),0);
        P[i][i] = 1;
    }
    for(int i = 0; i<A.size()-1; i++){
        double lead = INT64_MIN;
        double nlead = -1;
        for (int j = i; j < A.size(); j++){
            if (abs(A[j][i]) > lead){
                lead = abs(A[j][i]);
                nlead = j;
            }
        }
        swap(A[i],A[nlead]);
        swap(P[i],P[nlead]);
        for (int j = i+1; j < A.size(); j++){
            A[j][i] = A[j][i]/A[i][i];
            for (int k = i+1; k<A.size(); k++){
                A[j][k] = A[j][k]-A[j][i]*A[i][k];
            }
        }
    }
    return make_pair(A,P);
}
vector <double> LUPsolve(vector <vector<double>> A, vector<vector<double>> b){ //solves the equation system by using the results of LUP function
    pair<vector<vector<double>>,vector<vector<double>>> LpUaP = LUP (A);
    vector<vector<double>> LU = LpUaP.first;
    b = mmultiply(LpUaP.second,b);
    vector<double> y(b.size());
    for(int i = 0; i<b.size(); i++){
        y[i] = b[i][0];
    }
    for (int i = 0; i < A.size(); i++){
        for (int k = 0; k<i;k++){
            y[i]-=LU[i][k]*y[k];
        }
    }
    vector<double> x(b.size());
    for(int i = b.size()-1; i>=0; i--){
        x[i] = y[i];
        for (int k = i+1; k<b.size(); k++){
            x[i] -= LU[i][k]*x[k];
        }
        x[i] = x[i]/LU[i][i];
    }
    return x;
}
vector <vector<double>> data_gen(double a, double b, int n, int k, double alpha){
    std::normal_distribution<double> normal_dist(0, alpha);
    vector<vector<double>> result (n*k,vector<double>(2));
    double step = (b-a)/(n-1), init = a;
    for (int i = 0; i < n; i++){
        for (int j = 0; j < k; j++){
            result[i*k+j][0] = init;
            result[i*k+j][1] = f(init) + normal_dist(rng);
        }
        init += step;
    }
    return result;
}
double poly(vector<double> coeffs, double x){
    double result = 0;
    for (int i = 0; i < coeffs.size(); i++){
        result+=coeffs[i]*pow(x,i);
    }
    return result;
}
vector <double> vsub(vector<double>a, vector<double>b){
    if (a.size()<b.size()){
        swap(a,b);
    }
    for (int i = 0; i < b.size(); i++){
        a[i]-=b[i];
    }
    return a;
}
vector <double> vsumm(vector<double>a, vector<double>b){
    if (a.size()<b.size()){
        swap(a,b);
    }
    for (int i = 0; i < b.size(); i++){
        a[i]+=b[i];
    }
    return a;
}
double finderr (vector <double> coeffs, vector<vector<double>> table){
    double sum = 0;
    for (int i = 0; i < table.size(); i++){
        sum += pow(poly(coeffs,table[i][0]) - table[i][1],2);
    }
    return sum;
}
int main() {
    double a,b;
    int n,k,m;
    cin >> a >> b >> m >> k; //левая граница отрезка, правая граница отрезка, количество точек, количество измерений на точку.
    vector<vector<double>> table = data_gen(a,b,m,k,0.05), npolys, opolys;
    for (int i = 0; i < m*k;i++){
        cout << fixed << setprecision(16) <<table[i][0] << ' ' << table[i][1] << endl;
    }
    m = m*k; //после генерации данных смысл переменной меняю на общее количество измерений.
    for (n = 1; n <= 5; n++) {
        vector<vector<double>> E(m,vector<double>(n+1)), Et, vals(0),q;
        vals.push_back(getcolumn(table, 1));
        vals = transpose(vals);
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n + 1; j++) {
                E[i][j] = pow(table[i][0], j);
            }
        }
        Et = transpose(E);
        vector<double> n_coeffs = LUPsolve(mmultiply(Et, E), mmultiply(Et, vals)), o_coeffs(n + 1), op_coeffs;
        q.push_back({1});
        double sumxi = 0;
        for (int i = 0; i < m; i++) {
            sumxi += table[i][0];
        }
        q.push_back({-sumxi / m, 1});
        for (int i = 1; i <= n - 1; i++) {
            q.push_back(q[i]);
            q[i + 1].insert(q[i + 1].begin(), 0);
            double alpha = 0, beta = 0, sumsq1 = 0, sumsq2 = 0;
            for (int j = 0; j < m; j++) {
                alpha += table[j][0] * pow(poly(q[i], table[j][0]), 2);
                beta += table[j][0] * poly(q[i], table[j][0]) * poly(q[i - 1], table[j][0]);
                sumsq1 += pow(poly(q[i], table[j][0]), 2);
                sumsq2 += pow(poly(q[i - 1], table[j][0]), 2);
            }
            alpha = alpha / sumsq1;
            beta = beta / sumsq2;
            q[i + 1] = vsub(q[i + 1], vnumm(alpha, q[i]));
            q[i + 1] = vsub(q[i + 1], vnumm(beta, q[i - 1]));
        }
        for (int i = 0; i <= n; i++) {
            double sum = 0, sumsq = 0;
            for (int j = 0; j < m; j++) {
                sum += poly(q[i], table[j][0]) * table[j][1];
                sumsq += pow(poly(q[i], table[j][0]), 2);
            }
            o_coeffs[i] = sum / sumsq;
        }
        op_coeffs = vnumm(o_coeffs[n], q[n]);
        for (int i = n - 1; i >= 0; i--) {
            op_coeffs = vsumm(op_coeffs, vnumm(o_coeffs[i], q[i]));
        }
        cout << "Polynomial of order " << n << " calculated using normal equations:" << endl;
        cout << "y = " <<fixed << setprecision(16) << n_coeffs[0];
        for (int i = 1; i < n_coeffs.size();i++){
            cout << " + " << fixed <<setprecision(16) << n_coeffs[i] << "x^" << i;
        }
        cout << endl;
        cout << "Polynomial of order " << n << " calculated using orthogonal polynomials:" << endl;
        cout<< "y = " <<fixed << setprecision(16) << op_coeffs[0];
        for (int i = 1; i < op_coeffs.size();i++){
            cout << " + " <<fixed <<  setprecision(16) <<op_coeffs[i] << "x^" << i;
        }
        npolys.push_back(n_coeffs);
        opolys.push_back(op_coeffs);
        cout << endl;
        E.resize(0);
        Et.resize(0);
        vals.resize(0);
        q.resize(0);
        n_coeffs.resize(0);
        o_coeffs.resize(0);
        op_coeffs.resize(0);
    }
    cout << "Table:" << endl;
    cout << "n for NEM            for OPM" << endl;
    for (int i = 0 ;i < 5; i++){
        cout << i+1 << ' ' << finderr(npolys[i],table) << ' ' << finderr(opolys[i],table) << endl;
    }
    return 0;
}