#include "MatrixLibrary.h"
#include <iostream>
#include <iomanip>
#include <omp.h>
#include <vector>
#include <cmath>
#include <string>
#include <ctime>

using namespace std;

void main() {
	Matrix<float> M(3,3);
	cin >> M;
	cout << M.R_Inverse(0.01) << endl;
	system("pause");
}