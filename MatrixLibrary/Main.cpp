#include "MatrixLibrary.h"
#include <iostream>
#include <iomanip>
#include <omp.h>
#include <vector>
#include <cmath>
#include <string>
#include <ctime>

using namespace std;

float Inverse(float Value) {
	return Value * (-1);
}

void main() {
	/*TESTING MATRIX CLASS (POINTERS)*/
	cout << "Pointers\n" << endl;
	Matrix<float>* Test = new Matrix<float>(5,5);
	cout << "Randomizing\n" << endl;
	Test->P_Randomize()->P_Print();
	Matrix<float>* Pointers1 = new Matrix<float>(Test);
	cout << "Cut\n" << endl;
	Pointers1->P_Cut(0,0,2,2)->P_Print();
	cout << "Scalar Product\n" << endl;
	cout << Pointers1->P_ScalarProduct(Pointers1) << endl;
	cout << "Subtract\n" << endl;
	Pointers1->P_Subtract(Pointers1)->P_Print();
	cout << "Multiply\n" << endl;
	Pointers1->P_Multiply(Pointers1)->P_Print();
	cout << "Add\n" << endl;
	Pointers1->P_Add(Pointers1)->P_Print();
	cout << "Multiply\n" << endl;
	Pointers1->P_Multiply(3)->P_Print();
	cout << "Inverse\n" << endl;
	Pointers1->P_Inverse(0.01)->P_Print();
	cout << "Determinant\n" << endl;
	cout << Pointers1->P_Determinant() << endl;
	cout << "Use Function\n" << endl;
	Pointers1->P_UseFunction(Inverse)->P_Print();
	cout << "Reshape\n" << endl;
	Pointers1->P_Reshape(1, 25)->P_Print();
	cout << "Exponentiate\n" << endl;
	Pointers1->P_Exponentiate(2)->P_Print();
	cout << "Transpose\n" << endl;
	Pointers1->P_Transpose()->P_Print();
	Matrix<float>* Pointers2 = new Matrix<float>(Test);
	cout << "Cut\n" << endl;
	Matrix<float>::P_Cut(Pointers2,0, 0, 2, 2)->P_Print();
	cout << "Scalar Product\n" << endl;
	cout << Matrix<float>::P_ScalarProduct(Pointers2,Pointers2) << endl;
	cout << "Subtract\n" << endl;
	Matrix<float>::P_Subtract(Pointers2,Pointers2)->P_Print();
	cout << "Multiply\n" << endl;
	Matrix<float>::P_Multiply(Pointers2,Pointers2)->P_Print();
	cout << "Add\n" << endl;
	Matrix<float>::P_Add(Pointers2,Pointers2)->P_Print();
	cout << "Multiply\n" << endl;
	Matrix<float>::P_Multiply(Pointers2, 3)->P_Print();
	cout << "Inverse\n" << endl;
	Matrix<float>::P_Inverse(Pointers2, 0.01)->P_Print();
	cout << "Determinant\n" << endl;
	cout << Matrix<float>::P_Determinant(Pointers2) << endl;
	cout << "Use Function\n" << endl;
	Matrix<float>::P_UseFunction(Pointers2, Inverse)->P_Print();
	cout << "Reshape\n" << endl;
	Matrix<float>::P_Reshape(Pointers2, 1, 25)->P_Print();
	cout << "Exponentiate\n" << endl;
	Matrix<float>::P_Exponentiate(Pointers2, 2)->P_Print();
	cout << "Transpose\n" << endl;
	Matrix<float>::P_Transpose(Pointers2)->P_Print();
	cout << "References\n" << endl;
	Matrix<float> Reference(5,5);
	cout << "Randomizing\n" << endl;
	Reference.R_Randomize().R_Print();
	Matrix<float> Reference1(Test);
	cout << "Cut\n" << endl;
	Reference1.R_Cut(0, 0, 2, 2).R_Print();
	cout << "Scalar Product\n" << endl;
	cout << Reference1.R_ScalarProduct(Reference1) << endl;
	cout << "Subtract\n" << endl;
	Reference1.R_Subtract(Reference1).R_Print();
	cout << "Multiply\n" << endl;
	Reference1.R_Multiply(Reference1).R_Print();
	cout << "Add\n" << endl;
	Reference1.R_Add(Reference1).R_Print();
	cout << "Multiply\n" << endl;
	Reference1.R_Multiply(3).R_Print();
	cout << "Inverse\n" << endl;
	Reference1.R_Inverse(0.01).R_Print();
	cout << "Determinant\n" << endl;
	cout << Reference1.R_Determinant() << endl;
	cout << "Use Function\n" << endl;
	Reference1.R_UseFunction(Inverse).R_Print();
	cout << "Reshape\n" << endl;
	Reference1.R_Reshape(1, 25).R_Print();
	cout << "Exponentiate\n" << endl;
	Reference1.R_Exponentiate(2).R_Print();
	cout << "Transpose\n" << endl;
	Reference1.R_Transpose().R_Print();
	Matrix<float> Reference2(Test);
	cout << "Cut\n" << endl;
	Matrix<float>::R_Cut(Reference2, 0, 0, 2, 2).R_Print();
	cout << "Scalar Product\n" << endl;
	cout << Matrix<float>::R_ScalarProduct(Reference2, Reference2) << endl;
	cout << "Subtract\n" << endl;
	Matrix<float>::R_Subtract(Reference2, Reference2).R_Print();
	cout << "Multiply\n" << endl;
	Matrix<float>::R_Multiply(Reference2, Reference2).R_Print();
	cout << "Add\n" << endl;
	Matrix<float>::R_Add(Reference2, Reference2).R_Print();
	cout << "Multiply\n" << endl;
	Matrix<float>::R_Multiply(Reference2, 3).R_Print();
	cout << "Inverse\n" << endl;
	Matrix<float>::R_Inverse(Reference2, 0.01).R_Print();
	cout << "Determinant\n" << endl;
	cout << Matrix<float>::R_Determinant(Reference2) << endl;
	cout << "Use Function\n" << endl;
	Matrix<float>::R_UseFunction(Reference2, Inverse).R_Print();
	cout << "Reshape\n" << endl;
	Matrix<float>::R_Reshape(Reference2, 1, 25).R_Print();
	cout << "Exponentiate\n" << endl;
	Matrix<float>::R_Exponentiate(Reference2, 2).R_Print();
	cout << "Transpose\n" << endl;
	Matrix<float>::R_Transpose(Reference2).R_Print();
	cout << "END" << endl;
	system("pause");
}
