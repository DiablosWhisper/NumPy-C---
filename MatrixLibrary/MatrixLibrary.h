#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <omp.h>
#include <vector>
#include <cmath>
#include <string>
#include <ctime>
#pragma once

using namespace std;

template<typename Type>
class Matrix {

#pragma region
private:
	/*[Matrix] Inner Static Exceptions*/

	static void CheckSummationAndSubtraction(const Matrix* FirstMatrix, const Matrix* SecondMatrix) {

		if (FirstMatrix->Length != SecondMatrix->Length || FirstMatrix->Height != SecondMatrix->Height) {

			throw "The length/height of the first matrix doesn't equal the length/height of the second matrix";
		}
	};
	static void CheckMultiplication(const Matrix* FirstMatrix, const Matrix* SecondMatrix) {

		if (FirstMatrix->Length != SecondMatrix->Height) {

			throw "The length of the first matrix doesn't equal the height of the second matrix";
		}
	};
	static void CheckSquareness(const Matrix* SomeMatrix) {

		if (SomeMatrix->Length != SomeMatrix->Height) {

			throw "The length doesn't equal the height";
		}
	};
	static void CheckParameters(const Matrix* SomeMatrix) {

		if (SomeMatrix->Length == 0 || SomeMatrix->Height == 0) {

			throw "The length/height equals zero";
		}
	};

	/*[Matrix] Inner Exceptions*/

	void CheckSummationAndSubtraction(const Matrix* SomeMatrix) {

		if (Length != SomeMatrix->Length || Height != SomeMatrix->Height) {

			throw "The length/height of the first matrix doesn't equal the length/height of the second matrix";
		}
	};
	void CheckMultiplication(const Matrix* SomeMatrix) {

		if (Length != SomeMatrix->Height) {

			throw "The length of the first matrix doesn't equal the height of the second matrix";
		}
	};
	void CheckSquareness() {

		if (Length != Height) {

			throw "The length doesn't equal the height";
		}
	};
	void CheckParameters() {

		if (Length == 0 || Height == 0) {

			throw "The length/height equals zero";
		}
	};

	/*[Matrix] Properties*/

	Type** InnerMatrix;
	int Length;
	int Height;
#pragma endregion

#pragma region
public:
	/*[Matrix] Outter Static Matrix Methods*/

	static Matrix R_Cut(const Matrix& SomeMatrix, int StartX, int StartY, int EndX, int EndY) {

		CheckParameters(&SomeMatrix);

		Matrix<Type> Temporary(EndY - StartY, EndX - StartX);

		for (int FirstIndex = StartY; FirstIndex < EndY; FirstIndex++) {

			for (int SecondIndex = StartX; SecondIndex < EndX; SecondIndex++) {

				Temporary.InnerMatrix[FirstIndex - StartY][SecondIndex - StartX] = SomeMatrix.InnerMatrix[FirstIndex][SecondIndex];
			}
		}

		return Temporary;
	};
	static Type R_ScalarProduct(const Matrix& FirstMatrix, const Matrix& SecondMatrix) {

		CheckSummationAndSubtraction(&FirstMatrix, &SecondMatrix);

		Type Scalar = 0;

		for (int FirstIndex = 0; FirstIndex < FirstMatrix.Height; FirstIndex++) {

			for (int SecondIndex = 0; SecondIndex < FirstMatrix.Length; SecondIndex++) {

				Scalar += FirstMatrix.InnerMatrix[FirstIndex][SecondIndex] * SecondMatrix.InnerMatrix[FirstIndex][SecondIndex];
			}
		}

		return Scalar;
	};
	static Matrix R_Subtract(const Matrix& FirstMatrix, const Matrix& SecondMatrix) {

		CheckSummationAndSubtraction(&FirstMatrix, &SecondMatrix);

		Matrix<Type> Temporary(FirstMatrix.Height, FirstMatrix.Length);

		omp_set_num_threads(4);

		#pragma omp parallel
		{
			#pragma omp for

			for (int FirstIndex = 0; FirstIndex < FirstMatrix.Height; FirstIndex++) {

				for (int SecondIndex = 0; SecondIndex < FirstMatrix.Length; SecondIndex++) {

					Temporary.InnerMatrix[FirstIndex][SecondIndex] = FirstMatrix.InnerMatrix[FirstIndex][SecondIndex] - SecondMatrix.InnerMatrix[FirstIndex][SecondIndex];
				}
			}
		}

		return Temporary;
	};
	static Matrix R_Multiply(const Matrix& FirstMatrix, const Matrix& SecondMatrix) {

		CheckMultiplication(&FirstMatrix, &SecondMatrix);

		const int ThatLength = SecondMatrix.Length;
		const int ThatHeight = SecondMatrix.Height;

		const int ThisLength = FirstMatrix.Length;
		const int ThisHeight = FirstMatrix.Height;

		Matrix<Type> Temporary(ThisHeight, ThatLength);

		Type* ThatColumn = new Type[ThatHeight];

		omp_set_num_threads(4);

		#pragma omp parallel
		{
			#pragma omp for

			for (int SecondIndex = 0; SecondIndex < ThatLength; SecondIndex++) {

				for (int AuxiliaryIndex = 0; AuxiliaryIndex < ThisLength; AuxiliaryIndex++) {

					ThatColumn[AuxiliaryIndex] = SecondMatrix.InnerMatrix[AuxiliaryIndex][SecondIndex];
				}

				for (int FirstIndex = 0; FirstIndex < ThisHeight; FirstIndex++) {

					Type* ThisRow = FirstMatrix.InnerMatrix[FirstIndex];

					Type Summation = 0;

					for (int AuxiliaryIndex = 0; AuxiliaryIndex < ThisLength; AuxiliaryIndex++) {

						Summation += ThisRow[AuxiliaryIndex] * ThatColumn[AuxiliaryIndex];
					}

					Temporary.InnerMatrix[FirstIndex][SecondIndex] = Summation;
				}
			}
		}

		delete[] ThatColumn;

		return Temporary;
	};
	static Matrix R_Product(const Matrix& FirstMatrix, const Matrix& SecondMatrix) {

		CheckSummationAndSubtraction(&FirstMatrix, &SecondMatrix);

		Matrix<Type> Temporary(FirstMatrix.Height, FirstMatrix.Length);

		omp_set_num_threads(4);

		#pragma omp parallel
		{
			#pragma omp for

			for (int FirstIndex = 0; FirstIndex < FirstMatrix.Height; FirstIndex++) {

				for (int SecondIndex = 0; SecondIndex < FirstMatrix.Length; SecondIndex++) {

					Temporary.InnerMatrix[FirstIndex][SecondIndex] = FirstMatrix.InnerMatrix[FirstIndex][SecondIndex] * SecondMatrix.InnerMatrix[FirstIndex][SecondIndex];
				}
			}
		}

		return Temporary;
	};
	static Matrix R_Add(const Matrix& FirstMatrix, const Matrix& SecondMatrix) {

		CheckSummationAndSubtraction(&FirstMatrix, &SecondMatrix);

		Matrix<Type> Temporary(FirstMatrix.Height, FirstMatrix.Length);

		omp_set_num_threads(4);

		#pragma omp parallel
		{
			#pragma omp for

			for (int FirstIndex = 0; FirstIndex < FirstMatrix.Height; FirstIndex++) {

				for (int SecondIndex = 0; SecondIndex < FirstMatrix.Length; SecondIndex++) {

					Temporary.InnerMatrix[FirstIndex][SecondIndex] = FirstMatrix.InnerMatrix[FirstIndex][SecondIndex] + SecondMatrix.InnerMatrix[FirstIndex][SecondIndex];
				}
			}
		}

		return Temporary;
	};
	static Matrix R_Multiply(const Matrix& SomeMatrix, Type Coefficient) {

		CheckParameters(&SomeMatrix);

		Matrix<Type> Temporary(SomeMatrix.Height, SomeMatrix.Length);

		omp_set_num_threads(4);

		#pragma omp parallel
		{
			#pragma omp for

			for (int FirstIndex = 0; FirstIndex < SomeMatrix.Height; FirstIndex++) {

				for (int SecondIndex = 0; SecondIndex < SomeMatrix.Length; SecondIndex++) {

					Temporary.InnerMatrix[FirstIndex][SecondIndex] = SomeMatrix.InnerMatrix[FirstIndex][SecondIndex] * Coefficient;
				}
			}
		}

		return Temporary;
	};
	static Matrix R_Inverse(const Matrix& SomeMatrix, double Epsilon) {

		CheckSquareness(&SomeMatrix);

		Type SecondRate = 0, FirstRate = 0;

		Matrix<Type> SecondTemporary(SomeMatrix.Height, SomeMatrix.Length);

		SecondTemporary = SecondTemporary.R_UnitFilling().R_Multiply(2);

		Matrix<Type> FirstTemporary(SomeMatrix);

		for (int FirstIndex = 0; FirstIndex < SomeMatrix.Height; FirstIndex++) {

			Type ColumnSummation = 0, RowSummation = 0;

			for (int SecondIndex = 0; SecondIndex < SomeMatrix.Length; SecondIndex++) {

				ColumnSummation += fabs(FirstTemporary.InnerMatrix[SecondIndex][FirstIndex]);

				RowSummation += fabs(FirstTemporary.InnerMatrix[FirstIndex][SecondIndex]);
			}

			SecondRate = fmax(RowSummation, SecondRate);

			FirstRate = fmax(ColumnSummation, FirstRate);
		}

		FirstTemporary = FirstTemporary.R_Transpose().R_Multiply(1 / (FirstRate * SecondRate));

		Matrix<Type> InverseMatrix(FirstTemporary);

		if (SomeMatrix.R_Determinant() != 0) {

			for (; fabs((SomeMatrix.R_Multiply(InverseMatrix))->R_Determinant() - 1) >= Epsilon;) {

				Matrix<Type> PreviousStep(InverseMatrix);

				InverseMatrix = SomeMatrix.R_Multiply(PreviousStep);

				InverseMatrix = InverseMatrix.R_Multiply(-1).R_Add(SecondTemporary);

				InverseMatrix = PreviousStep.R_Multiply(InverseMatrix);
			}

			return InverseMatrix;
		}

		return nullptr;
	};
	static Type R_Determinant(const Matrix& SomeMatrix) {

		CheckSquareness(&SomeMatrix);

		Matrix<Type> Temporary(SomeMatrix);

		double Value = 1;

		omp_set_num_threads(4);

		#pragma omp parallel
		{
			#pragma omp for

			for (int FirstIndex = 0; FirstIndex < SomeMatrix.Height - 1; FirstIndex++) {

				for (int SecondIndex = FirstIndex + 1; SecondIndex < SomeMatrix.Length; SecondIndex++) {

					double Coefficient = -Temporary.InnerMatrix[SecondIndex][FirstIndex] / Temporary.InnerMatrix[FirstIndex][FirstIndex];

					for (int AuxiliaryIndex = FirstIndex; AuxiliaryIndex < SomeMatrix.Height; AuxiliaryIndex++) {

						Temporary.InnerMatrix[SecondIndex][AuxiliaryIndex] += Temporary.InnerMatrix[FirstIndex][AuxiliaryIndex] * Coefficient;
					}
				}
			}
		}

		for (int Index = 0; Index < SomeMatrix.Height; Index++) {

			Value *= Temporary.InnerMatrix[Index][Index];
		}

		return Value;
	};

	static Matrix* P_Cut(Matrix* SomeMatrix,int StartX, int StartY, int EndX, int EndY) {

		CheckParameters(SomeMatrix);

		Matrix<Type>* Temporary = new Matrix<Type>(EndY - StartY, EndX - StartX);

		for (int FirstIndex = StartY; FirstIndex < EndY; FirstIndex++) {

			for (int SecondIndex = StartX; SecondIndex < EndX; SecondIndex++) {

				Temporary->InnerMatrix[FirstIndex - StartY][SecondIndex - StartX] = SomeMatrix->InnerMatrix[FirstIndex][SecondIndex];
			}
		}

		return Temporary;
	};
	static Type P_ScalarProduct(Matrix* FirstMatrix, Matrix* SecondMatrix) {

		CheckSummationAndSubtraction(FirstMatrix, SecondMatrix);

		Type Scalar = 0;

		for (int FirstIndex = 0; FirstIndex < FirstMatrix->Height; FirstIndex++) {

			for (int SecondIndex = 0; SecondIndex < FirstMatrix->Length; SecondIndex++) {

				Scalar += FirstMatrix->InnerMatrix[FirstIndex][SecondIndex] * SecondMatrix->InnerMatrix[FirstIndex][SecondIndex];
			}
		}

		return Scalar;
	};
	static Matrix* P_Subtract(Matrix* FirstMatrix, Matrix* SecondMatrix) {

		CheckSummationAndSubtraction(FirstMatrix, SecondMatrix);

		Matrix<Type>* Temporary = new Matrix<Type>(FirstMatrix->Height, FirstMatrix->Length);

		omp_set_num_threads(4);

		#pragma omp parallel
		{
			#pragma omp for

			for (int FirstIndex = 0; FirstIndex < FirstMatrix->Height; FirstIndex++) {

				for (int SecondIndex = 0; SecondIndex < FirstMatrix->Length; SecondIndex++) {

					Temporary->InnerMatrix[FirstIndex][SecondIndex] = FirstMatrix->InnerMatrix[FirstIndex][SecondIndex] - SecondMatrix->InnerMatrix[FirstIndex][SecondIndex];
				}
			}
		}

		return Temporary;
	};
	static Matrix* P_Multiply(Matrix* FirstMatrix, Matrix* SecondMatrix) {

		CheckMultiplication(FirstMatrix, SecondMatrix);

		const int ThatLength = SecondMatrix->Length;
		const int ThatHeight = SecondMatrix->Height;

		const int ThisLength = FirstMatrix->Length;
		const int ThisHeight = FirstMatrix->Height;

		Matrix<Type>* Temporary = new Matrix<Type>(ThisHeight, ThatLength);

		Type* ThatColumn = new Type[ThatHeight];

		omp_set_num_threads(4);

		#pragma omp parallel
		{
			#pragma omp for

			for (int SecondIndex = 0; SecondIndex < ThatLength; SecondIndex++) {

				for (int AuxiliaryIndex = 0; AuxiliaryIndex < ThisLength; AuxiliaryIndex++) {

					ThatColumn[AuxiliaryIndex] = SecondMatrix->InnerMatrix[AuxiliaryIndex][SecondIndex];
				}

				for (int FirstIndex = 0; FirstIndex < ThisHeight; FirstIndex++) {

					Type* ThisRow = FirstMatrix->InnerMatrix[FirstIndex];

					Type Summation = 0;

					for (int AuxiliaryIndex = 0; AuxiliaryIndex < ThisLength; AuxiliaryIndex++) {

						Summation += ThisRow[AuxiliaryIndex] * ThatColumn[AuxiliaryIndex];
					}

					Temporary->InnerMatrix[FirstIndex][SecondIndex] = Summation;
				}
			}
		}

		delete[] ThatColumn;

		return Temporary;
	};
	static Matrix* P_Product(Matrix* FirstMatrix, Matrix* SecondMatrix) {

		CheckSummationAndSubtraction(FirstMatrix, SecondMatrix);

		Matrix<Type>* Temporary = new Matrix<Type>(FirstMatrix->Height, FirstMatrix->Length);

		omp_set_num_threads(4);

		#pragma omp parallel
		{
			#pragma omp for

			for (int FirstIndex = 0; FirstIndex < FirstMatrix->Height; FirstIndex++) {

				for (int SecondIndex = 0; SecondIndex < FirstMatrix->Length; SecondIndex++) {

					Temporary->InnerMatrix[FirstIndex][SecondIndex] = FirstMatrix->InnerMatrix[FirstIndex][SecondIndex] * SecondMatrix->InnerMatrix[FirstIndex][SecondIndex];
				}
			}
		}

		return Temporary;
	};
	static Matrix* P_Add(Matrix* FirstMatrix, Matrix* SecondMatrix) {

		CheckSummationAndSubtraction(FirstMatrix, SecondMatrix);

		Matrix<Type>* Temporary = new Matrix<Type>(FirstMatrix->Height, FirstMatrix->Length);

		omp_set_num_threads(4);

		#pragma omp parallel
		{
			#pragma omp for

			for (int FirstIndex = 0; FirstIndex < FirstMatrix->Height; FirstIndex++) {

				for (int SecondIndex = 0; SecondIndex < FirstMatrix->Length; SecondIndex++) {

					Temporary->InnerMatrix[FirstIndex][SecondIndex] = FirstMatrix->InnerMatrix[FirstIndex][SecondIndex] + SecondMatrix->InnerMatrix[FirstIndex][SecondIndex];
				}
			}
		}

		return Temporary;
	};
	static Matrix* P_Multiply(Matrix* SomeMatrix, Type Coefficient) {

		CheckParameters(SomeMatrix);

		Matrix<Type>* Temporary = new Matrix<Type>(SomeMatrix->Height, SomeMatrix->Length);

		omp_set_num_threads(4);

		#pragma omp parallel
		{
			#pragma omp for

			for (int FirstIndex = 0; FirstIndex < SomeMatrix->Height; FirstIndex++) {

				for (int SecondIndex = 0; SecondIndex < SomeMatrix->Length; SecondIndex++) {

					Temporary->InnerMatrix[FirstIndex][SecondIndex] = SomeMatrix->InnerMatrix[FirstIndex][SecondIndex] * Coefficient;
				}
			}
		}

		return Temporary;
	};
	static Matrix* P_Inverse(Matrix* SomeMatrix, double Epsilon) {

		CheckSquareness(SomeMatrix);

		Type SecondRate = 0, FirstRate = 0;

		Matrix<Type>* SecondTemporary = new Matrix<Type>(SomeMatrix->Height, SomeMatrix->Length);

		SecondTemporary = SecondTemporary->P_UnitFilling()->P_Multiply(2);

		Matrix<Type>* FirstTemporary = new Matrix<Type>(SomeMatrix);

		for (int FirstIndex = 0; FirstIndex < SomeMatrix->Height; FirstIndex++) {

			Type ColumnSummation = 0, RowSummation = 0;

			for (int SecondIndex = 0; SecondIndex < SomeMatrix->Length; SecondIndex++) {

				ColumnSummation += fabs(FirstTemporary->InnerMatrix[SecondIndex][FirstIndex]);

				RowSummation += fabs(FirstTemporary->InnerMatrix[FirstIndex][SecondIndex]);
			}

			SecondRate = fmax(RowSummation, SecondRate);

			FirstRate = fmax(ColumnSummation, FirstRate);
		}

		FirstTemporary = FirstTemporary->P_Transpose()->P_Multiply(1 / (FirstRate * SecondRate));

		Matrix<Type>* InverseMatrix = new Matrix(FirstTemporary);

		delete FirstTemporary;

		if (SomeMatrix->P_Determinant() != 0) {

			for (; fabs((SomeMatrix->P_Multiply(InverseMatrix))->P_Determinant() - 1) >= Epsilon;) {

				Matrix<Type>* PreviousStep = new Matrix<Type>(InverseMatrix);

				InverseMatrix = SomeMatrix->P_Multiply(PreviousStep);

				InverseMatrix = InverseMatrix->P_Multiply(-1)->P_Add(SecondTemporary);

				InverseMatrix = PreviousStep->P_Multiply(InverseMatrix);

				delete PreviousStep;
			}

			delete SecondTemporary;

			return InverseMatrix;
		}

		return nullptr;
	};
	static Type P_Determinant(Matrix* SomeMatrix) {

		CheckSquareness(SomeMatrix);

		Matrix<Type>* Temporary = new Matrix(SomeMatrix);

		double Value = 1;

		omp_set_num_threads(4);

		#pragma omp parallel
		{
			#pragma omp for

			for (int FirstIndex = 0; FirstIndex < SomeMatrix->Height - 1; FirstIndex++) {

				for (int SecondIndex = FirstIndex + 1; SecondIndex < SomeMatrix->Length; SecondIndex++) {

					double Coefficient = -Temporary->InnerMatrix[SecondIndex][FirstIndex] / Temporary->InnerMatrix[FirstIndex][FirstIndex];

					for (int AuxiliaryIndex = FirstIndex; AuxiliaryIndex < SomeMatrix->Height; AuxiliaryIndex++) {

						Temporary->InnerMatrix[SecondIndex][AuxiliaryIndex] += Temporary->InnerMatrix[FirstIndex][AuxiliaryIndex] * Coefficient;
					}
				}
			}
		}

		for (int Index = 0; Index < SomeMatrix->Height; Index++) {

			Value *= Temporary->InnerMatrix[Index][Index];
		}

		delete Temporary;

		return Value;
	};

	/*[Matrix] Outter Static Methods*/

	static Matrix R_LoadFromFIle(const Matrix& SomeMatrix, string Path, char Delimiter) {

		ifstream Input(Path);

		if (!Input.is_open()) {

			return nullptr;
		}

		int Height = 0, Length = 0;

		string Item, Text;

		getline(Input, Text, '\n');

		Height++;

		stringstream Stream(Text);

		for (; getline(Stream, Item, Delimiter);) {

			Length++;
		}

		for (; getline(Input, Text, '\n');) {

			Height++;
		}

		Input.clear();

		Input.seekg(0, Input.beg);

		SomeMatrix(Height, Length);

		Height = 0, Length = 0;

		if (Input.is_open()) {

			for (; getline(Input, Text, '\n');) {

				stringstream FirstStream(Text);

				for (; getline(FirstStream, Item, Delimiter);) {

					stringstream SecondStream(Item);

					SecondStream >> SomeMatrix->InnerMatrix[Height][Length];

					SecondStream.clear();

					Length++;
				}

				Length = 0;

				Height++;

				FirstStream.clear();
			}

			Input.close();

			Height = 0;
		}

		return SomeMatrix;
	};
	static Matrix R_UseFunction(const Matrix& SomeMatrix, Type(*Function)(Type)) {

		CheckParameters(SomeMatrix);

		Matrix<Type> Temporary(SomeMatrix.Height, SomeMatrix.Length);

		for (int FirstIndex = 0; FirstIndex < SomeMatrix.Height; FirstIndex++) {

			for (int SecondIndex = 0; SecondIndex < SomeMatrix.Length; SecondIndex++) {

				Temporary.InnerMatrix[FirstIndex][SecondIndex] = Function(SomeMatrix.InnerMatrix[FirstIndex][SecondIndex]);
			}
		}

		return Temporary;
	};
	static Matrix R_Reshape(const Matrix& SomeMatrix, int Height, int Length) {

		CheckParameters(&SomeMatrix);

		Type* DeconvolutedMatrix = new Type[SomeMatrix.Length * SomeMatrix.Height];

		Matrix<Type> Temporary(Height, Length);

		int InnerHeight = 0;

		for (int FirstIndex = 0; FirstIndex < SomeMatrix.Height; FirstIndex++) {

			for (int SecondIndex = 0; SecondIndex < SomeMatrix.Length; SecondIndex++) {

				DeconvolutedMatrix[InnerHeight] = SomeMatrix.InnerMatrix[FirstIndex][SecondIndex];

				InnerHeight++;
			}
		}

		InnerHeight = 0;

		int Counter = 0;

		for (int Index = 0; Index < SomeMatrix.Length * SomeMatrix.Height; Index++) {

			if (Counter == Length) {

				Counter = 0;

				InnerHeight++;
			}

			Temporary.InnerMatrix[InnerHeight][Counter] = DeconvolutedMatrix[Index];

			Counter++;
		}

		delete[] DeconvolutedMatrix;

		return Temporary;
	};
	static Matrix R_Exponentiate(const Matrix& SomeMatrix, int Power) {

		CheckParameters(&SomeMatrix);

		Matrix<Type> Temporary(SomeMatrix.Height, SomeMatrix.Length);

		for (int Index = 0; Index < Power; Index++) {

			Temporary = Temporary.R_Add(SomeMatrix.R_Multiply(SomeMatrix));
		}

		return Temporary;
	};
	static Matrix R_Fill(const Matrix& SomeMatrix, Type Value) {

		CheckParameters(&SomeMatrix);

		for (int FirstIndex = 0; FirstIndex < SomeMatrix.Height; FirstIndex++) {

			for (int SecondIndex = 0; SecondIndex < SomeMatrix.Length; SecondIndex++) {

				SomeMatrix.InnerMatrix[FirstIndex][SecondIndex] = Value;
			}
		}

		return SomeMatrix;
	};
	static Matrix R_Randomize(const Matrix& SomeMatrix) {

		CheckParameters(&SomeMatrix);

		srand(1 + rand() % 100);

		for (int FirstIndex = 0; FirstIndex < SomeMatrix.Height; FirstIndex++) {

			for (int SecondIndex = 0; SecondIndex < SomeMatrix.Length; SecondIndex++) {

				SomeMatrix.InnerMatrix[FirstIndex][SecondIndex] = static_cast <Type> (rand()) / static_cast <Type> (RAND_MAX) - 0.5f;
			}
		}

		return SomeMatrix;
	};
	static Matrix R_UnitFilling(const Matrix& SomeMatrix) {

		CheckSquareness(&SomeMatrix);

		for (int Index = 0; Index < SomeMatrix.Height; Index++) {

			SomeMatrix.InnerMatrix[Index][Index] = 1;
		}

		return SomeMatrix;
	};
	static Matrix R_Transpose(const Matrix& SomeMatrix) {

		CheckParameters(&SomeMatrix);

		Matrix<Type> Temporary(SomeMatrix.Length, SomeMatrix.Height);

		omp_set_num_threads(4);

		#pragma omp parallel
		{
			#pragma omp for

			for (int FirstIndex = 0; FirstIndex < SomeMatrix.Length; FirstIndex++) {

				for (int SecondIndex = 0; SecondIndex < SomeMatrix.Height; SecondIndex++) {

					Temporary.InnerMatrix[FirstIndex][SecondIndex] = SomeMatrix.InnerMatrix[SecondIndex][FirstIndex];
				}
			}
		}

		return Temporary;
	};
	static Type R_Average(const Matrix& SomeMatrix) {

		CheckParameters(&SomeMatrix);

		Type Value = 0;

		for (int FirstIndex = 0; FirstIndex < SomeMatrix.Height; FirstIndex++) {

			for (int SecondIndex = 0; SecondIndex < SomeMatrix.Length; SecondIndex++) {

				Value += SomeMatrix.InnerMatrix[FirstIndex][SecondIndex];
			}
		}

		return Value / (SomeMatrix.Height * SomeMatrix.Length);
	};
	static Matrix R_Clear(const Matrix& SomeMatrix) {

		CheckParameters(&SomeMatrix);

		for (int FirstIndex = 0; FirstIndex < SomeMatrix.Height; FirstIndex++) {

			for (int SecondIndex = 0; SecondIndex < SomeMatrix.Length; SecondIndex++) {

				SomeMatrix.InnerMatrix[FirstIndex][SecondIndex] = 0;
			}
		}

		return SomeMatrix;
	};
	static Matrix R_Print(const Matrix& SomeMatrix) {

		CheckParameters(&SomeMatrix);

		for (int FirstIndex = 0; FirstIndex < SomeMatrix.Height; FirstIndex++) {

			for (int SecondIndex = 0; SecondIndex < SomeMatrix.Length; SecondIndex++) {

				cout << " " << SomeMatrix.InnerMatrix[FirstIndex][SecondIndex];
			}

			cout << endl;
		}

		cout << endl;

		return SomeMatrix;
	};
	static Type R_Rate(const Matrix& SomeMatrix) {

		CheckParameters(&SomeMatrix);

		Type Value = 0;

		for (int FirstIndex = 0; FirstIndex < SomeMatrix.Height; FirstIndex++) {

			for (int SecondIndex = 0; SecondIndex < SomeMatrix.Length; SecondIndex++) {

				Value += pow(SomeMatrix.InnerMatrix[FirstIndex][SecondIndex], 2);
			}
		}

		return sqrt(Value);
	};
	static Type R_Max(const Matrix& SomeMatrix) {

		CheckParameters(&SomeMatrix);

		Type Value = 0;

		for (int FirstIndex = 0; FirstIndex < SomeMatrix.Height; FirstIndex++) {

			for (int SecondIndex = 0; SecondIndex < SomeMatrix.Length; SecondIndex++) {

				Value = fmax(Value, SomeMatrix.InnerMatrix[FirstIndex][SecondIndex]);
			}
		}

		return Value;
	};

	static Matrix* P_LoadFromFIle(Matrix* SomeMatrix ,string Path, char Delimiter) {

		ifstream Input(Path);

		if (!Input.is_open()) {

			return nullptr;
		}

		int Height = 0, Length = 0;

		string Item, Text;

		getline(Input, Text, '\n');

		Height++;

		stringstream Stream(Text);

		for (; getline(Stream, Item, Delimiter);) {

			Length++;
		}

		for (; getline(Input, Text, '\n');) {

			Height++;
		}

		Input.clear();

		Input.seekg(0, Input.beg);

		SomeMatrix = new Matrix<Type>(Height, Length);

		Height = 0, Length = 0;

		if (Input.is_open()) {

			for (; getline(Input, Text, '\n');) {

				stringstream FirstStream(Text);

				for (; getline(FirstStream, Item, Delimiter);) {

					stringstream SecondStream(Item);

					SecondStream >> SomeMatrix->InnerMatrix[Height][Length];

					SecondStream.clear();

					Length++;
				}

				Length = 0;

				Height++;

				FirstStream.clear();
			}

			Input.close();

			Height = 0;
		}

		return SomeMatrix;
	};
	static Matrix* P_UseFunction(Matrix* SomeMatrix, Type(*Function)(Type)) {

		CheckParameters(SomeMatrix);

		Matrix<Type>* Temporary = new Matrix<Type>(SomeMatrix->Height, SomeMatrix->Length);

		for (int FirstIndex = 0; FirstIndex < SomeMatrix->Height; FirstIndex++) {

			for (int SecondIndex = 0; SecondIndex < SomeMatrix->Length; SecondIndex++) {

				Temporary->InnerMatrix[FirstIndex][SecondIndex] = Function(SomeMatrix->InnerMatrix[FirstIndex][SecondIndex]);
			}
		}

		return Temporary;
	};
	static Matrix* P_Reshape(Matrix* SomeMatrix, int Height, int Length) {

		CheckParameters(SomeMatrix);

		Type* DeconvolutedMatrix = new Type[SomeMatrix->Length * SomeMatrix->Height];

		Matrix<Type>* Temporary = new Matrix<Type>(Height, Length);

		int InnerHeight = 0;

		for (int FirstIndex = 0; FirstIndex < SomeMatrix->Height; FirstIndex++) {

			for (int SecondIndex = 0; SecondIndex < SomeMatrix->Length; SecondIndex++) {

				DeconvolutedMatrix[InnerHeight] = SomeMatrix->InnerMatrix[FirstIndex][SecondIndex];

				InnerHeight++;
			}
		}

		InnerHeight = 0;

		int Counter = 0;

		for (int Index = 0; Index < SomeMatrix->Length * SomeMatrix->Height; Index++) {

			if (Counter == Length) {

				Counter = 0;

				InnerHeight++;
			}

			Temporary->InnerMatrix[InnerHeight][Counter] = DeconvolutedMatrix[Index];

			Counter++;
		}

		delete[] DeconvolutedMatrix;

		return Temporary;
	};
	static Matrix* P_Exponentiate(Matrix* SomeMatrix, int Power) {

		Matrix<Type>* Temporary = new Matrix<Type>(SomeMatrix->Height, SomeMatrix->Length);

		for (int Index = 0; Index < Power; Index++) {

			Temporary = Temporary->P_Add(SomeMatrix->P_Multiply(SomeMatrix));
		}

		return Temporary;
	};
	static Matrix* P_Fill(Matrix* SomeMatrix, Type Value) {

		CheckParameters(SomeMatrix);

		for (int FirstIndex = 0; FirstIndex < SomeMatrix->Height; FirstIndex++) {

			for (int SecondIndex = 0; SecondIndex < SomeMatrix->Length; SecondIndex++) {

				SomeMatrix->InnerMatrix[FirstIndex][SecondIndex] = Value;
			}
		}

		return SomeMatrix;
	};
	static Matrix* P_Randomize(Matrix* SomeMatrix) {

		CheckParameters(SomeMatrix);

		srand(1 + rand() % 100);

		for (int FirstIndex = 0; FirstIndex < SomeMatrix->Height; FirstIndex++) {

			for (int SecondIndex = 0; SecondIndex < SomeMatrix->Length; SecondIndex++) {

				SomeMatrix->InnerMatrix[FirstIndex][SecondIndex] = static_cast <Type> (rand()) / static_cast <Type> (RAND_MAX) - 0.5f;
			}
		}

		return SomeMatrix;
	};
	static Matrix* P_UnitFilling(Matrix* SomeMatrix) {

		CheckSquareness(SomeMatrix);

		for (int Index = 0; Index < SomeMatrix->Height; Index++) {

			SomeMatrix->InnerMatrix[Index][Index] = 1;
		}

		return SomeMatrix;
	};
	static Matrix* P_Transpose(Matrix* SomeMatrix) {

		CheckParameters(SomeMatrix);

		Matrix<Type>* Temporary = new Matrix<Type>(SomeMatrix->Length, SomeMatrix->Height);

		omp_set_num_threads(4);

		#pragma omp parallel
		{
			#pragma omp for

			for (int FirstIndex = 0; FirstIndex < SomeMatrix->Length; FirstIndex++) {

				for (int SecondIndex = 0; SecondIndex < SomeMatrix->Height; SecondIndex++) {

					Temporary->InnerMatrix[FirstIndex][SecondIndex] = SomeMatrix->InnerMatrix[SecondIndex][FirstIndex];
				}
			}
		}

		return Temporary;
	};
	static Matrix* P_Clear(Matrix* SomeMatrix) {

		CheckParameters(SomeMatrix);

		for (int FirstIndex = 0; FirstIndex < SomeMatrix->Height; FirstIndex++) {

			for (int SecondIndex = 0; SecondIndex < SomeMatrix->Length; SecondIndex++) {

				SomeMatrix->InnerMatrix[FirstIndex][SecondIndex] = 0;
			}
		}

		return SomeMatrix;
	};
	static Matrix* P_Print(Matrix* SomeMatrix) {

		CheckParameters(SomeMatrix);

		for (int FirstIndex = 0; FirstIndex < SomeMatrix->Height; FirstIndex++) {

			for (int SecondIndex = 0; SecondIndex < SomeMatrix->Length; SecondIndex++) {

				cout << " " << SomeMatrix->InnerMatrix[FirstIndex][SecondIndex];
			}

			cout << endl;
		}

		cout << endl;

		return SomeMatrix;
	};
	static Type P_Average(Matrix* SomeMatrix) {

		CheckParameters(SomeMatrix);

		Type Value = 0;

		for (int FirstIndex = 0; FirstIndex < SomeMatrix->Height; FirstIndex++) {

			for (int SecondIndex = 0; SecondIndex < SomeMatrix->Length; SecondIndex++) {

				Value += SomeMatrix->InnerMatrix[FirstIndex][SecondIndex];
			}
		}

		return Value / (SomeMatrix->Height * SomeMatrix->Length);
	};
	static Type P_Rate(Matrix* SomeMatrix) {

		CheckParameters(SomeMatrix);

		Type Value = 0;

		for (int FirstIndex = 0; FirstIndex < SomeMatrix->Height; FirstIndex++) {

			for (int SecondIndex = 0; SecondIndex < SomeMatrix->Length; SecondIndex++) {

				Value += pow(SomeMatrix->InnerMatrix[FirstIndex][SecondIndex], 2);
			}
		}

		return sqrt(Value);
	};
	static Type P_Max(Matrix* SomeMatrix) {

		CheckParameters(SomeMatrix);

		Type Value = 0;

		for (int FirstIndex = 0; FirstIndex < SomeMatrix->Height; FirstIndex++) {

			for (int SecondIndex = 0; SecondIndex < SomeMatrix->Length; SecondIndex++) {

				Value = fmax(Value, SomeMatrix->InnerMatrix[FirstIndex][SecondIndex]);
			}
		}

		return Value;
	};

	/*[Matrix] Outter Matrix Methods*/

	Matrix R_Cut(int StartX, int StartY, int EndX, int EndY) {

		CheckParameters();

		Matrix<Type> Temporary(EndY - StartY, EndX - StartX);

		for (int FirstIndex = StartY; FirstIndex < EndY; FirstIndex++) {

			for (int SecondIndex = StartX; SecondIndex < EndX; SecondIndex++) {

				Temporary.InnerMatrix[FirstIndex - StartY][SecondIndex - StartX] = InnerMatrix[FirstIndex][SecondIndex];
			}
		}

		return Temporary;
	};
	Type R_ScalarProduct(const Matrix& SomeMatrix) {

		CheckSummationAndSubtraction(&SomeMatrix);

		Type Value = 0;

		for (int FirstIndex = 0; FirstIndex < Height; FirstIndex++) {

			for (int SecondIndex = 0; SecondIndex < Length; SecondIndex++) {

				Value += InnerMatrix[FirstIndex][SecondIndex] * SomeMatrix.InnerMatrix[FirstIndex][SecondIndex];
			}
		}

		return Value;
	};
	Matrix R_Subtract(const Matrix& SomeMatrix) {

		CheckSummationAndSubtraction(&SomeMatrix);

		Matrix<Type> Temporary(Height, Length);

		omp_set_num_threads(4);

		#pragma omp parallel
		{
			#pragma omp for

			for (int FirstIndex = 0; FirstIndex < Height; FirstIndex++) {

				for (int SecondIndex = 0; SecondIndex < Length; SecondIndex++) {

					Temporary.InnerMatrix[FirstIndex][SecondIndex] = InnerMatrix[FirstIndex][SecondIndex] - SomeMatrix.InnerMatrix[FirstIndex][SecondIndex];
				}
			}
		}

		return Temporary;
	};
	Matrix R_Multiply(const Matrix& SomeMatrix) {

		CheckMultiplication(&SomeMatrix);

		const int ThatLength = SomeMatrix.Length;
		const int ThatHeight = SomeMatrix.Height;

		const int ThisLength = Length;
		const int ThisHeight = Height;

		Matrix<Type> Temporary(ThisHeight, ThatLength);

		Type* ThatColumn = new Type[ThatHeight];

		omp_set_num_threads(4);

		#pragma omp parallel
		{
			#pragma omp for

			for (int SecondIndex = 0; SecondIndex < ThatLength; SecondIndex++) {

				for (int AuxiliaryIndex = 0; AuxiliaryIndex < ThisLength; AuxiliaryIndex++) {

					ThatColumn[AuxiliaryIndex] = SomeMatrix.InnerMatrix[AuxiliaryIndex][SecondIndex];
				}

				for (int FirstIndex = 0; FirstIndex < ThisHeight; FirstIndex++) {

					Type* ThisRow = InnerMatrix[FirstIndex];

					Type Summation = 0;

					for (int AuxiliaryIndex = 0; AuxiliaryIndex < ThisLength; AuxiliaryIndex++) {

						Summation += ThisRow[AuxiliaryIndex] * ThatColumn[AuxiliaryIndex];
					}

					Temporary.InnerMatrix[FirstIndex][SecondIndex] = Summation;
				}
			}
		}

		delete[] ThatColumn;

		return Temporary;
	};
	Matrix R_Product(const Matrix& SomeMatrix) {

		CheckSummationAndSubtraction(&SomeMatrix);

		Matrix<Type> Temporary(Height, Length);

		omp_set_num_threads(4);

		#pragma omp parallel
		{
			#pragma omp for

			for (int FirstIndex = 0; FirstIndex < Height; FirstIndex++) {

				for (int SecondIndex = 0; SecondIndex < Length; SecondIndex++) {

					Temporary.InnerMatrix[FirstIndex][SecondIndex] = InnerMatrix[FirstIndex][SecondIndex] * SomeMatrix.InnerMatrix[FirstIndex][SecondIndex];
				}
			}
		}

		return Temporary;
	};
	Matrix R_Add(const Matrix& SomeMatrix) {

		CheckSummationAndSubtraction(&SomeMatrix);

		Matrix<Type> Temporary(Height, Length);

		omp_set_num_threads(4);

		#pragma omp parallel
		{
			#pragma omp for

			for (int FirstIndex = 0; FirstIndex < Height; FirstIndex++) {

				for (int SecondIndex = 0; SecondIndex < Length; SecondIndex++) {

					Temporary.InnerMatrix[FirstIndex][SecondIndex] = InnerMatrix[FirstIndex][SecondIndex] + SomeMatrix.InnerMatrix[FirstIndex][SecondIndex];
				}
			}
		}

		return Temporary;
	};
	Matrix R_Multiply(Type Coefficient) {

		CheckParameters();

		Matrix<Type> Temporary(Height, Length);

		omp_set_num_threads(4);

		#pragma omp parallel
		{
			#pragma omp for

			for (int FirstIndex = 0; FirstIndex < Height; FirstIndex++) {

				for (int SecondIndex = 0; SecondIndex < Length; SecondIndex++) {

					Temporary.InnerMatrix[FirstIndex][SecondIndex] = InnerMatrix[FirstIndex][SecondIndex] * Coefficient;
				}
			}
		}

		return Temporary;
	};
	Matrix R_Inverse(double Epsilon) {

		CheckSquareness();

		Type SecondRate = 0, FirstRate = 0;

		Matrix<Type> SecondTemporary(Height, Length);

		SecondTemporary = SecondTemporary.R_UnitFilling().R_Multiply(2);

		Matrix<Type> FirstTemporary(this);

		for (int FirstIndex = 0; FirstIndex < Height; FirstIndex++) {

			Type ColumnSummation = 0, RowSummation = 0;

			for (int SecondIndex = 0; SecondIndex < Length; SecondIndex++) {

				ColumnSummation += fabs(FirstTemporary.InnerMatrix[SecondIndex][FirstIndex]);

				RowSummation += fabs(FirstTemporary.InnerMatrix[FirstIndex][SecondIndex]);
			}

			SecondRate = fmax(RowSummation, SecondRate);

			FirstRate = fmax(ColumnSummation, FirstRate);
		}

		FirstTemporary = FirstTemporary.R_Transpose().R_Multiply(1 / (FirstRate * SecondRate));

		Matrix<Type> InverseMatrix(FirstTemporary);

		if (this->R_Determinant() != 0) {

			for (; fabs((this->R_Multiply(InverseMatrix)).R_Determinant() - 1) >= Epsilon;) {

				Matrix<Type> PreviousStep(InverseMatrix);

				InverseMatrix = this->R_Multiply(PreviousStep);

				InverseMatrix = InverseMatrix.R_Multiply(-1).R_Add(SecondTemporary);

				InverseMatrix = PreviousStep.R_Multiply(InverseMatrix);
			}

			return InverseMatrix;
		}

		return nullptr;
	};
	Type R_Determinant() {

		CheckSquareness();

		Matrix<Type> Temporary(this);

		double Value = 1;

		omp_set_num_threads(4);

		#pragma omp parallel
		{
			#pragma omp for

			for (int FirstIndex = 0; FirstIndex < Height - 1; FirstIndex++) {

				for (int SecondIndex = FirstIndex + 1; SecondIndex < Length; SecondIndex++) {

					double Coefficient = -Temporary.InnerMatrix[SecondIndex][FirstIndex] / Temporary.InnerMatrix[FirstIndex][FirstIndex];

					for (int AuxiliaryIndex = FirstIndex; AuxiliaryIndex < Height; AuxiliaryIndex++) {

						Temporary.InnerMatrix[SecondIndex][AuxiliaryIndex] += Temporary.InnerMatrix[FirstIndex][AuxiliaryIndex] * Coefficient;
					}
				}
			}
		}

		for (int Index = 0; Index < Height; Index++) {

			Value *= Temporary.InnerMatrix[Index][Index];
		}

		return Value;
	};

	Matrix* P_Cut(int StartX, int StartY, int EndX, int EndY) {

		CheckParameters();

		Matrix<Type>* Temporary = new Matrix<Type>(EndY - StartY, EndX - StartX);

		for (int FirstIndex = StartY; FirstIndex < EndY; FirstIndex++) {

			for (int SecondIndex = StartX; SecondIndex < EndX; SecondIndex++) {

				Temporary->InnerMatrix[FirstIndex - StartY][SecondIndex - StartX] = InnerMatrix[FirstIndex][SecondIndex];
			}
		}

		return Temporary;
	};
	Type P_ScalarProduct(Matrix* SomeMatrix) {

		CheckSummationAndSubtraction(SomeMatrix);

		Type Scalar = 0;

		for (int FirstIndex = 0; FirstIndex < Height; FirstIndex++) {

			for (int SecondIndex = 0; SecondIndex < Length; SecondIndex++) {

				Scalar += InnerMatrix[FirstIndex][SecondIndex] * SomeMatrix->InnerMatrix[FirstIndex][SecondIndex];
			}
		}

		return Scalar;
	};
	Matrix* P_Subtract(Matrix* SomeMatrix) {

		CheckSummationAndSubtraction(SomeMatrix);

		Matrix<Type>* Temporary = new Matrix<Type>(Height, Length);

		omp_set_num_threads(4);

		#pragma omp parallel
		{
			#pragma omp for

			for (int FirstIndex = 0; FirstIndex < Height; FirstIndex++) {

				for (int SecondIndex = 0; SecondIndex < Length; SecondIndex++) {

					Temporary->InnerMatrix[FirstIndex][SecondIndex] = InnerMatrix[FirstIndex][SecondIndex] - SomeMatrix->InnerMatrix[FirstIndex][SecondIndex];
				}
			}
		}

		return Temporary;
	};
	Matrix* P_Multiply(Matrix* SomeMatrix) {

		CheckMultiplication(SomeMatrix);

		const int ThatLength = SomeMatrix->Length;
		const int ThatHeight = SomeMatrix->Height;

		const int ThisLength = Length;
		const int ThisHeight = Height;

		Matrix<Type>* Temporary = new Matrix<Type>(ThisHeight, ThatLength);

		Type* ThatColumn = new Type[ThatHeight];

		omp_set_num_threads(4);

		#pragma omp parallel
		{
			#pragma omp for

			for (int SecondIndex = 0; SecondIndex < ThatLength; SecondIndex++) {

				for (int AuxiliaryIndex = 0; AuxiliaryIndex < ThisLength; AuxiliaryIndex++) {

					ThatColumn[AuxiliaryIndex] = SomeMatrix->InnerMatrix[AuxiliaryIndex][SecondIndex];
				}

				for (int FirstIndex = 0; FirstIndex < ThisHeight; FirstIndex++) {

					Type* ThisRow = InnerMatrix[FirstIndex];

					Type Summation = 0;

					for (int AuxiliaryIndex = 0; AuxiliaryIndex < ThisLength; AuxiliaryIndex++) {

						Summation += ThisRow[AuxiliaryIndex] * ThatColumn[AuxiliaryIndex];
					}

					Temporary->InnerMatrix[FirstIndex][SecondIndex] = Summation;
				}
			}
		}

		delete[] ThatColumn;

		return Temporary;
	};
	Matrix* P_Product(Matrix* SomeMatrix) {

		CheckSummationAndSubtraction(SomeMatrix);

		Matrix<Type>* Temporary = new Matrix<Type>(Height, Length);

		omp_set_num_threads(4);

		#pragma omp parallel
		{
			#pragma omp for

			for (int FirstIndex = 0; FirstIndex < Height; FirstIndex++) {

				for (int SecondIndex = 0; SecondIndex < Length; SecondIndex++) {

					Temporary->InnerMatrix[FirstIndex][SecondIndex] = InnerMatrix[FirstIndex][SecondIndex] * SomeMatrix->InnerMatrix[FirstIndex][SecondIndex];
				}
			}
		}

		return Temporary;
	};
	Matrix* P_Add(Matrix* SomeMatrix) {

		CheckSummationAndSubtraction(SomeMatrix);

		Matrix<Type>* Temporary = new Matrix<Type>(Height, Length);

		omp_set_num_threads(4);

		#pragma omp parallel
		{
			#pragma omp for

			for (int FirstIndex = 0; FirstIndex < Height; FirstIndex++) {

				for (int SecondIndex = 0; SecondIndex < Length; SecondIndex++) {

					Temporary->InnerMatrix[FirstIndex][SecondIndex] = InnerMatrix[FirstIndex][SecondIndex] + SomeMatrix->InnerMatrix[FirstIndex][SecondIndex];
				}
			}
		}

		return Temporary;
	};
	Matrix* P_Multiply(Type Coefficient) {

		CheckParameters();

		Matrix<Type>* Temporary = new Matrix<Type>(Height, Length);

		omp_set_num_threads(4);

		#pragma omp parallel
		{
			#pragma omp for

			for (int FirstIndex = 0; FirstIndex < Height; FirstIndex++) {

				for (int SecondIndex = 0; SecondIndex < Length; SecondIndex++) {

					Temporary->InnerMatrix[FirstIndex][SecondIndex] = InnerMatrix[FirstIndex][SecondIndex] * Coefficient;
				}
			}
		}

		return Temporary;
	};
	Matrix* P_Inverse(double Epsilon) {

		CheckSquareness();

		Type SecondRate = 0, FirstRate = 0;

		Matrix<Type>* SecondTemporary = new Matrix<Type>(Height, Length);

		SecondTemporary = SecondTemporary->P_UnitFilling()->P_Multiply(2);

		Matrix<Type>* FirstTemporary = new Matrix<Type>(this);

		for (int FirstIndex = 0; FirstIndex < Height; FirstIndex++) {

			Type ColumnSummation = 0, RowSummation = 0;

			for (int SecondIndex = 0; SecondIndex < Length; SecondIndex++) {

				ColumnSummation += fabs(FirstTemporary->InnerMatrix[SecondIndex][FirstIndex]);

				RowSummation += fabs(FirstTemporary->InnerMatrix[FirstIndex][SecondIndex]);
			}

			SecondRate = fmax(RowSummation, SecondRate);

			FirstRate = fmax(ColumnSummation, FirstRate);
		}

		FirstTemporary = FirstTemporary->P_Transpose()->P_Multiply(1 / (FirstRate * SecondRate));

		Matrix<Type>* InverseMatrix = new Matrix(FirstTemporary);

		delete FirstTemporary;

		if (this->P_Determinant() != 0) {

			for (; fabs((this->P_Multiply(InverseMatrix))->P_Determinant() - 1) >= Epsilon;) {

				Matrix<Type>* PreviousStep = new Matrix<Type>(InverseMatrix);

				InverseMatrix = this->P_Multiply(PreviousStep);

				InverseMatrix = InverseMatrix->P_Multiply(-1)->P_Add(SecondTemporary);

				InverseMatrix = PreviousStep->P_Multiply(InverseMatrix);

				delete PreviousStep;
			}

			delete SecondTemporary;

			return InverseMatrix;
		}

		return nullptr;
	};
	Type P_Determinant() {

		CheckSquareness();

		Matrix<Type>* Temporary = new Matrix(this);

		double Value = 1;

		omp_set_num_threads(4);

		#pragma omp parallel
		{
			#pragma omp for

			for (int FirstIndex = 0; FirstIndex < Height - 1; FirstIndex++) {

				for (int SecondIndex = FirstIndex + 1; SecondIndex < Length; SecondIndex++) {

					double Coefficient = -Temporary->InnerMatrix[SecondIndex][FirstIndex] / Temporary->InnerMatrix[FirstIndex][FirstIndex];

					for (int AuxiliaryIndex = FirstIndex; AuxiliaryIndex < Height; AuxiliaryIndex++) {

						Temporary->InnerMatrix[SecondIndex][AuxiliaryIndex] += Temporary->InnerMatrix[FirstIndex][AuxiliaryIndex] * Coefficient;
					}
				}
			}
		}

		for (int Index = 0; Index < Height; Index++) {

			Value *= Temporary->InnerMatrix[Index][Index];
		}

		delete Temporary;

		return Value;
	};

	/*[Matrix] Outter Methods*/

	Matrix R_LoadFromFIle(string Path, char Delimiter) {

		ifstream Input(Path);

		if (!Input.is_open()) {

			return nullptr;
		}

		int Height = 0, Length = 0;

		string Item, Text;

		getline(Input, Text, '\n');

		Height++;

		stringstream Stream(Text);

		for (; getline(Stream, Item, Delimiter);) {

			Length++;
		}

		for (; getline(Input, Text, '\n');) {

			Height++;
		}

		Input.clear();

		Input.seekg(0, Input.beg);

		Matrix<Type> Temporary(Height, Length);

		Height = 0, Length = 0;

		if (Input.is_open()) {

			for (; getline(Input, Text, '\n');) {

				stringstream FirstStream(Text);

				for (; getline(FirstStream, Item, Delimiter);) {

					stringstream SecondStream(Item);

					SecondStream >> Temporary.InnerMatrix[Height][Length];

					SecondStream.clear();

					Length++;
				}

				Length = 0;

				Height++;

				FirstStream.clear();
			}

			Input.close();

			Height = 0;
		}

		return Temporary;
	};
	Matrix R_UseFunction(Type(*Function)(Type)) {

		CheckParameters();

		Matrix<Type>* Temporary = new Matrix<Type>(Height, Length);

		for (int FirstIndex = 0; FirstIndex < Height; FirstIndex++) {

			for (int SecondIndex = 0; SecondIndex < Length; SecondIndex++) {

				Temporary->InnerMatrix[FirstIndex][SecondIndex] = Function(InnerMatrix[FirstIndex][SecondIndex]);
			}
		}

		return Temporary;
	};
	Matrix R_Reshape(int Height, int Length) {

		CheckParameters();

		Type* DeconvolutedMatrix = new Type[this->Length * this->Height];

		Matrix<Type> Temporary(Height, Length);

		int InnerHeight = 0;

		for (int FirstIndex = 0; FirstIndex < this->Height; FirstIndex++) {

			for (int SecondIndex = 0; SecondIndex < this->Length; SecondIndex++) {

				DeconvolutedMatrix[InnerHeight] = InnerMatrix[FirstIndex][SecondIndex];

				InnerHeight++;
			}
		}

		InnerHeight = 0;

		int Counter = 0;

		for (int Index = 0; Index < this->Length * this->Height; Index++) {

			if (Counter == Length) {

				Counter = 0;

				InnerHeight++;
			}

			Temporary.InnerMatrix[InnerHeight][Counter] = DeconvolutedMatrix[Index];

			Counter++;
		}

		delete[] DeconvolutedMatrix;

		return Temporary;
	};
	Matrix R_Exponentiate(int Power) {

		CheckParameters();

		Matrix<Type> Temporary(Height, Length);

		for (int Index = 0; Index < Power; Index++) {

			Temporary = Temporary.R_Add(this->R_Multiply(this));
		}

		return Temporary;
	};
	Matrix R_Fill(Type Value) {

		CheckParameters();

		for (int FirstIndex = 0; FirstIndex < Height; FirstIndex++) {

			for (int SecondIndex = 0; SecondIndex < Length; SecondIndex++) {

				InnerMatrix[FirstIndex][SecondIndex] = Value;
			}
		}

		return this;
	};
	Matrix R_Randomize() {

		CheckParameters();

		srand(1 + rand() % 100);

		for (int FirstIndex = 0; FirstIndex < Height; FirstIndex++) {

			for (int SecondIndex = 0; SecondIndex < Length; SecondIndex++) {

				InnerMatrix[FirstIndex][SecondIndex] = static_cast <Type> (rand()) / static_cast <Type> (RAND_MAX) - 0.5f;
			}
		}

		return this;
	};
	Matrix R_UnitFilling() {

		CheckSquareness();

		for (int Index = 0; Index < Height; Index++) {

			InnerMatrix[Index][Index] = 1;
		}

		return this;
	};
	Matrix R_Transpose() {

		CheckParameters();

		Matrix<Type> Temporary(Length, Height);

		omp_set_num_threads(4);

		#pragma omp parallel
		{
			#pragma omp for

			for (int FirstIndex = 0; FirstIndex < Length; FirstIndex++) {

				for (int SecondIndex = 0; SecondIndex < Height; SecondIndex++) {

					Temporary.InnerMatrix[FirstIndex][SecondIndex] = InnerMatrix[SecondIndex][FirstIndex];
				}
			}
		}

		return Temporary;
	};
	Type R_Average() {

		CheckParameters();

		Type Value = 0;

		for (int FirstIndex = 0; FirstIndex < Height; FirstIndex++) {

			for (int SecondIndex = 0; SecondIndex < Length; SecondIndex++) {

				Value += InnerMatrix[FirstIndex][SecondIndex];
			}
		}

		return Value / (Height * Length);
	}
	Matrix R_Clear() {

		CheckParameters();

		for (int FirstIndex = 0; FirstIndex < Height; FirstIndex++) {

			for (int SecondIndex = 0; SecondIndex < Length; SecondIndex++) {

				InnerMatrix[FirstIndex][SecondIndex] = 0;
			}
		}

		return this;
	};
	Matrix R_Print() {

		CheckParameters();

		for (int FirstIndex = 0; FirstIndex < Height; FirstIndex++) {

			for (int SecondIndex = 0; SecondIndex < Length; SecondIndex++) {

				cout << " " << InnerMatrix[FirstIndex][SecondIndex];
			}

			cout << endl;
		}

		cout << endl;

		return this;
	};
	Type R_Rate() {

		CheckParameters();

		Type Value = 0;

		for (int FirstIndex = 0; FirstIndex < Height; FirstIndex++) {

			for (int SecondIndex = 0; SecondIndex < Length; SecondIndex++) {

				Value += pow(InnerMatrix[FirstIndex][SecondIndex], 2);
			}
		}

		return sqrt(Value);
	};
	Type R_Max() {

		CheckParameters();

		Type Value = 0;

		for (int FirstIndex = 0; FirstIndex < Height; FirstIndex++) {

			for (int SecondIndex = 0; SecondIndex < Length; SecondIndex++) {

				Value = fmax(Value, InnerMatrix[FirstIndex][SecondIndex]);
			}
		}

		return Value;
	};

	Matrix* P_LoadFromFIle(string Path, char Delimiter) {

		ifstream Input(Path);

		if (!Input.is_open()) {

			return nullptr;
		}

		int Height = 0, Length = 0;

		string Item,Text;

		getline(Input, Text, '\n');

		Height++;

		stringstream Stream(Text);

		for (; getline(Stream, Item, Delimiter);) {

			Length++;
		}

		for (; getline(Input, Text, '\n');) {

			Height++;
		}

		Input.clear();

		Input.seekg(0, Input.beg);

		Matrix<Type>* Temporary = new Matrix<Type>(Height, Length);

		Height = 0, Length = 0;

		if (Input.is_open()) {

			for (; getline(Input, Text, '\n');) {

				stringstream FirstStream(Text);

				for (; getline(FirstStream, Item, Delimiter);) {

					stringstream SecondStream(Item);

					SecondStream >> Temporary->InnerMatrix[Height][Length];

					SecondStream.clear();

					Length++;
				}

				Length = 0;

				Height++;

				FirstStream.clear();
			}

			Input.close();

			Height = 0;
		}

		return Temporary;
	};
	Matrix* P_UseFunction(Type(*Function)(Type)) {

		CheckParameters();

		Matrix<Type>* Temporary = new Matrix<Type>(Height, Length);

		for (int FirstIndex = 0; FirstIndex < Height; FirstIndex++) {

			for (int SecondIndex = 0; SecondIndex < Length; SecondIndex++) {

				Temporary->InnerMatrix[FirstIndex][SecondIndex] = Function(InnerMatrix[FirstIndex][SecondIndex]);
			}
		}

		return Temporary;
	};
	Matrix* P_Reshape(int Height, int Length) {

		CheckParameters();

		Type* DeconvolutedMatrix = new Type[this->Length * this->Height];

		Matrix<Type>* Temporary = new Matrix<Type>(Height, Length);

		int InnerHeight = 0;

		for (int FirstIndex = 0; FirstIndex < this->Height; FirstIndex++) {

			for (int SecondIndex = 0; SecondIndex < this->Length; SecondIndex++) {

				DeconvolutedMatrix[InnerHeight] = InnerMatrix[FirstIndex][SecondIndex];

				InnerHeight++;
			}
		}

		InnerHeight = 0;

		int Counter = 0;

		for (int Index = 0; Index < this->Length * this->Height; Index++) {

			if (Counter == Length) {

				Counter = 0;

				InnerHeight++;
			}

			Temporary->InnerMatrix[InnerHeight][Counter] = DeconvolutedMatrix[Index];

			Counter++;
		}

		delete[] DeconvolutedMatrix;

		return Temporary;
	};
	Matrix* P_Exponentiate(int Power) {

		CheckParameters();

		Matrix<Type>* Temporary = new Matrix<Type>(Height, Length);

		for (int Index = 0; Index < Power; Index++) {

			Temporary = Temporary->P_Add(this->P_Multiply(this));
		}

		return Temporary;
	};
	Matrix* P_Fill(Type Value) {

		CheckParameters();

		for (int FirstIndex = 0; FirstIndex < Height; FirstIndex++) {

			for (int SecondIndex = 0; SecondIndex < Length; SecondIndex++) {

				InnerMatrix[FirstIndex][SecondIndex] = Value;
			}
		}

		return this;
	};
	Matrix* P_Randomize() {

		CheckParameters();

		srand(1 + rand() % 100);

		for (int FirstIndex = 0; FirstIndex < Height; FirstIndex++) {

			for (int SecondIndex = 0; SecondIndex < Length; SecondIndex++) {

				InnerMatrix[FirstIndex][SecondIndex] = static_cast <Type> (rand()) / static_cast <Type> (RAND_MAX) - 0.5f;
			}
		}

		return this;
	};
	Matrix* P_UnitFilling() {

		CheckSquareness();

		for (int Index = 0; Index < Height; Index++) {

			InnerMatrix[Index][Index] = 1;
		}

		return this;
	};
	Matrix* P_Transpose() {

		CheckParameters();

		Matrix<Type>* Temporary = new Matrix<Type>(Length, Height);

		omp_set_num_threads(4);

		#pragma omp parallel
		{
			#pragma omp for

			for (int FirstIndex = 0; FirstIndex < Length; FirstIndex++) {

				for (int SecondIndex = 0; SecondIndex < Height; SecondIndex++) {

					Temporary->InnerMatrix[FirstIndex][SecondIndex] = InnerMatrix[SecondIndex][FirstIndex];
				}
			}
		}

		return Temporary;
	};
	Matrix* P_Clear() {

		CheckParameters();

		for (int FirstIndex = 0; FirstIndex < Height; FirstIndex++) {

			for (int SecondIndex = 0; SecondIndex < Length; SecondIndex++) {

				InnerMatrix[FirstIndex][SecondIndex] = 0;
			}
		}

		return this;
	};
	Matrix* P_Print() {

		CheckParameters();

		for (int FirstIndex = 0; FirstIndex < Height; FirstIndex++) {

			for (int SecondIndex = 0; SecondIndex < Length; SecondIndex++) {

				cout << " " << InnerMatrix[FirstIndex][SecondIndex];
			}

			cout << endl;
		}

		cout << endl;

		return this;
	};
	Type P_Average() {

		CheckParameters();

		Type Value = 0;

		for (int FirstIndex = 0; FirstIndex < Height; FirstIndex++) {

			for (int SecondIndex = 0; SecondIndex < Length; SecondIndex++) {

				Value += InnerMatrix[FirstIndex][SecondIndex];
			}
		}

		return Value / (Height * Length);
	}
	Type P_Rate() {

		CheckParameters();

		Type Value = 0;

		for (int FirstIndex = 0; FirstIndex < Height; FirstIndex++) {

			for (int SecondIndex = 0; SecondIndex < Length; SecondIndex++) {

				Value += pow(InnerMatrix[FirstIndex][SecondIndex], 2);
			}
		}

		return sqrt(Value);
	};
	Type P_Max() {

		CheckParameters();

		Type Value = 0;

		for (int FirstIndex = 0; FirstIndex < Height; FirstIndex++) {

			for (int SecondIndex = 0; SecondIndex < Length; SecondIndex++) {

				Value = fmax(Value, InnerMatrix[FirstIndex][SecondIndex]);
			}
		}

		return Value;
	};

	/*[Matrix] Outter Getters*/

	int GetLength() {

		return Length;
	};
	int GetHeight() {

		return Height;
	};
	Type** GetMatrix() {

		return InnerMatrix;
	};

	/*[Matrix] Properties*/

	Matrix(const Matrix& SomeMatrix) : Length(SomeMatrix.Length), Height(SomeMatrix.Height) {

		InnerMatrix = new Type*[Height];

		for (int FirstIndex = 0; FirstIndex < Height; FirstIndex++) {

			InnerMatrix[FirstIndex] = new Type[Length];

			for (int SecondIndex = 0; SecondIndex < Length; SecondIndex++) {

				InnerMatrix[FirstIndex][SecondIndex] = SomeMatrix.InnerMatrix[FirstIndex][SecondIndex];
			}
		}
	};
	Matrix(int Height, int Length) : Length(Length), Height(Height) {

		InnerMatrix = new Type*[Height];

		for (int FirstIndex = 0; FirstIndex < Height; FirstIndex++) {

			InnerMatrix[FirstIndex] = new Type[Length];

			for (int SecondIndex = 0; SecondIndex < Length; SecondIndex++) {

				InnerMatrix[FirstIndex][SecondIndex] = 0;
			}
		}
	};
	Matrix(Matrix* SomeMatrix) : Length(SomeMatrix->Length), Height(SomeMatrix-> Height) {

		InnerMatrix = new Type*[Height];

		for (int FirstIndex = 0; FirstIndex < Height; FirstIndex++) {

			InnerMatrix[FirstIndex] = new Type[Length];

			for (int SecondIndex = 0; SecondIndex < Length; SecondIndex++) {

				InnerMatrix[FirstIndex][SecondIndex] = SomeMatrix->InnerMatrix[FirstIndex][SecondIndex];
			}
		}
	};
	Matrix(int Height) : Length(1), Height(Height) {

		InnerMatrix = new Type*[Height];

		for (int Index = 0; Index < Height; Index++) {

			InnerMatrix[Index] = new Type[1];

			InnerMatrix[Index][0] = 0;
		}
	};
	~Matrix() {

		for (int Index = 0; Index < Height; Index++) {

			delete[] InnerMatrix[Index];
		}

		delete[] InnerMatrix;
	};
	Matrix() : Length(0), Height(0){

		InnerMatrix = nullptr;
	};

	/*[Matrix] Operators*/

	template<typename Type>
	friend ostream& operator << (ostream& Stream, const Matrix<Type>& SomeMatrix) {

		for (int FirstIndex = 0; FirstIndex < SomeMatrix.Height; FirstIndex++) {

			for (int SecondIndex = 0; SecondIndex < SomeMatrix.Length; SecondIndex++) {

				Stream << setw(to_string(SomeMatrix.InnerMatrix[FirstIndex][SecondIndex]).size()) << SomeMatrix.InnerMatrix[FirstIndex][SecondIndex];
			}

			Stream << endl;
		}

		return Stream;
	};
	template<typename Type>
	friend istream& operator >> (istream& Stream, const Matrix<Type>& SomeMatrix) {

		for (int FirstIndex = 0; FirstIndex < SomeMatrix.Height; FirstIndex++) {

			for (int SecondIndex = 0; SecondIndex < SomeMatrix.Length; SecondIndex++) {

				Stream >> SomeMatrix.InnerMatrix[FirstIndex][SecondIndex];
			}
		}

		return Stream;
	};
	template<typename Type>
	friend ostream& operator << (ostream& Stream, Matrix<Type>* SomeMatrix) {

		if (SomeMatrix != nullptr) {

			for (int FirstIndex = 0; FirstIndex < SomeMatrix->Height; FirstIndex++) {

				for (int SecondIndex = 0; SecondIndex < SomeMatrix->Length; SecondIndex++) {

					Stream << setw(to_string(SomeMatrix->InnerMatrix[FirstIndex][SecondIndex]).size()) << SomeMatrix->InnerMatrix[FirstIndex][SecondIndex];
				}

				Stream << endl;
			}
		}

		return Stream;
	};
	template<typename Type>
	friend istream& operator >> (istream& Stream, Matrix<Type>* SomeMatrix) {

		if (SomeMatrix != nullptr) {

			for (int FirstIndex = 0; FirstIndex < SomeMatrix->Height; FirstIndex++) {

				for (int SecondIndex = 0; SecondIndex < SomeMatrix->Length; SecondIndex++) {

					Stream >> SomeMatrix->InnerMatrix[FirstIndex][SecondIndex];
				}
			}
		}

		return Stream;
	};
	Matrix operator = (const Matrix& SomeMatrix) {

		InnerMatrix = new Type*[Height];

		for (int FirstIndex = 0; FirstIndex < Height; FirstIndex++) {

			InnerMatrix[FirstIndex] = new Type[Length];

			for (int SecondIndex = 0; SecondIndex < Length; SecondIndex++) {

				InnerMatrix[FirstIndex][SecondIndex] = SomeMatrix.InnerMatrix[FirstIndex][SecondIndex];
			}
		}

		return this;
	};
	bool operator == (Matrix* SomeMatrix) {

		for (int FirstIndex = 0; FirstIndex < Height; FirstIndex++) {

			for (int SecondIndex = 0; SecondIndex < Length; SecondIndex++) {

				if (InnerMatrix[FirstIndex][SecondIndex] != SomeMatrix->InnerMatrix[FirstIndex][SecondIndex]) {

					return false;
				}
			}
		}

		return true;
	};
	Type* operator[] (int Index) {

		return InnerMatrix[Index];
	};
#pragma endregion
};

#pragma region
/*[Matrix] Additional Functions*/

template<typename Type>
inline Type** ConvertMatrixToArrayMatrix(Matrix<Type>* SomeMatrix) {

	int Length = SomeMatrix->GetLength();

	int Height = SomeMatrix->GetHeight();

	Type** Temporary = new Type*[Height];

	for (int FirstIndex = 0; FirstIndex < Height; FirstIndex++) {

		Temporary[FirstIndex] = new Type[Length];

		for (int SecondIndex = 0; SecondIndex < Length; SecondIndex++) {

			Temporary[FirstIndex][SecondIndex] = SomeMatrix->GetMatrix()[FirstIndex][SecondIndex];
		}
	}

	return Temporary;
};
template<typename Type>
inline Matrix<Type>* ConvertArrayMatrixToMatrix(Type** SomeMatrix) {

	int Length = _msize(SomeMatrix[0]) / sizeof(SomeMatrix[0][0]);

	int Height = _msize(SomeMatrix) / sizeof(SomeMatrix[0]);

	Matrix<Type>* Temporary = new Matrix<Type>(Height, Length);

	for (int FirstIndex = 0; FirstIndex < Height; FirstIndex++) {

		for (int SecondIndex = 0; SecondIndex < Length; SecondIndex++) {

			Temporary->GetMatrix()[FirstIndex][SecondIndex] = SomeMatrix[FirstIndex][SecondIndex];
		}
	}

	return Temporary;
};
template<typename Type>
inline Type* ConvertVectorToArray(Matrix<Type>* SomeVector) {

	int Length = SomeVector->GetLength();

	Type* Temporary = new Type[Length];

	if (Length == 1) {

		Length = SomeVector->GetHeight();

		Temporary = new Type[Length];

		for (int Index = 0; Index < Length; Index++) {

			Temporary[Index] = SomeVector->GetMatrix()[Index][0];
		}
	}
	else {

		for (int Index = 0; Index < Length; Index++) {

			Temporary[Index] = SomeVector->GetMatrix()[0][Index];
		}
	}

	return Temporary;
};
template<typename Type>
inline Matrix<Type>* ConvertArrayToVector(Type* SomeArray) {

	int Length = _msize(SomeArray) / sizeof(SomeArray[0]);

	Matrix<Type>* Temporary = new Matrix<Type>(Length);

	for (int Index = 0; Index < Length; Index++) {

		Temporary->GetMatrix()[Index][0] = SomeArray[Index];
	}

	return Temporary;
};
#pragma endregion
