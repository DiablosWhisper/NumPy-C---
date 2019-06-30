#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <omp.h>
#include <vector>
#include <cmath>
#include <string>
#include <ctime>
#include <tuple>
#pragma once

const int NUMBER_OF_CORES = 4;

using namespace std;

#pragma region
template<typename Type>
class Matrix;
#pragma endregion

template<typename Type>
class LinearEquationSystem {

#pragma region
private:
	void CheckSquareness(const Matrix<Type>& SomeMatrix) const {

		if (SomeMatrix.GetLength() != SomeMatrix.GetHeight()) {

			throw "The length doesn't equal the height";
		}
	};
	void CheckParameters(const Matrix<Type>& SomeMatrix) const {

		if (SomeMatrix.GetLength() == 0 || SomeMatrix.GetHeight() == 0) {

			throw "The length/height equals zero";
		}
	};
	double CheckDivisionByZero(double Value) const {

		if (Value == -INFINITY || Value == INFINITY || Value == NAN || Value != Value) {

			return 0;
		}

		return Value;
	};
#pragma endregion

#pragma region
public:
	Matrix<Type> JacobiRotation(const Matrix<Type>& SystemOfLinearEquations, const double Precision) const {

		CheckSquareness(SystemOfLinearEquations);

		const int Length = SystemOfLinearEquations.GetLength();

		const int Height = SystemOfLinearEquations.GetHeight();

		Matrix<Type> Temporary(SystemOfLinearEquations);

		Matrix<Type> Diagonal(Height, Length);

		Matrix<Type> Solution(Height);

		Diagonal.UnitFilling();

		double MaxCoefficient = 0, AngleFi = 0, Fault = 0;

		int MaxCoordinateX = 0, MaxCoordinateY = 0;

		#pragma omp parallel
		{
			#pragma omp for

			for (int FirstIndex = 0; FirstIndex < Height; FirstIndex++) {

				for (int SecondIndex = FirstIndex + 1; SecondIndex < Height; SecondIndex++) {

					Fault += Temporary.GetMatrix()[FirstIndex][SecondIndex] * Temporary.GetMatrix()[FirstIndex][SecondIndex];
				}
			}

			Fault = sqrt(Fault * 2);

			do {

				MaxCoefficient = 0;

				for (int FirstIndex = 0; FirstIndex < Height; FirstIndex++) {

					for (int SecondIndex = FirstIndex + 1; SecondIndex < Length; SecondIndex++) {

						if (Temporary.GetMatrix()[FirstIndex][SecondIndex] > 0 && Temporary.GetMatrix()[FirstIndex][SecondIndex] > MaxCoefficient) {

							MaxCoefficient = Temporary.GetMatrix()[FirstIndex][SecondIndex];

							MaxCoordinateY = SecondIndex;

							MaxCoordinateX = FirstIndex;
						}
						else if (Temporary.GetMatrix()[FirstIndex][SecondIndex] < 0 && -Temporary.GetMatrix()[FirstIndex][SecondIndex] > MaxCoefficient) {

							MaxCoefficient = -Temporary.GetMatrix()[FirstIndex][SecondIndex];

							MaxCoordinateY = SecondIndex;

							MaxCoordinateX = FirstIndex;
						}
					}
				}

				Matrix<Type> Rotation(Height, Length);

				Rotation.UnitFilling();

				if (Temporary.GetMatrix()[MaxCoordinateX][MaxCoordinateX] == Temporary.GetMatrix()[MaxCoordinateY][MaxCoordinateY]) {

					Rotation.GetMatrix()[MaxCoordinateX][MaxCoordinateX] = Rotation.GetMatrix()[MaxCoordinateY][MaxCoordinateY] = Rotation.GetMatrix()[MaxCoordinateY][MaxCoordinateX] = sqrt(2.0) / 2.0;

					Rotation.GetMatrix()[MaxCoordinateX][MaxCoordinateY] = -sqrt(2.0) / 2.0;
				}
				else {

					AngleFi = 0.5 * atan((2.0 * Temporary.GetMatrix()[MaxCoordinateX][MaxCoordinateY]) / (Temporary.GetMatrix()[MaxCoordinateX][MaxCoordinateX] - Temporary.GetMatrix()[MaxCoordinateY][MaxCoordinateY]));

					Rotation.GetMatrix()[MaxCoordinateX][MaxCoordinateX] = Rotation.GetMatrix()[MaxCoordinateY][MaxCoordinateY] = cos(AngleFi);

					Rotation.GetMatrix()[MaxCoordinateX][MaxCoordinateY] = -sin(AngleFi);

					Rotation.GetMatrix()[MaxCoordinateY][MaxCoordinateX] = sin(AngleFi);
				}

				Temporary = Rotation.Transpose().Multiply(Temporary).Multiply(Rotation);

				Fault = 0.0;

				for (int FirstIndex = 0; FirstIndex < Height; FirstIndex++) {

					for (int SecondIndex = FirstIndex + 1; SecondIndex < Height; SecondIndex++) {

						Fault += Temporary.GetMatrix()[FirstIndex][SecondIndex] * Temporary.GetMatrix()[FirstIndex][SecondIndex];
					}
				}

				Fault = sqrt(Fault * 2);

				Diagonal *= Rotation;

			} while (Fault > Precision);

			for (int Index = 0; Index < Height; Index++) {

				Solution.GetMatrix()[Index][0] = Temporary.GetMatrix()[Index][Index];
			}
		}

		return Solution;
	};
	Matrix<Type> GaussSeidel(const Matrix<Type>& SystemOfLinearEquations, const double Precision) const {

		CheckParameters(SystemOfLinearEquations);

		const int Length = SystemOfLinearEquations.GetLength();

		const int Height = SystemOfLinearEquations.GetHeight();

		Matrix<Type> Temporary = SystemOfLinearEquations.Cut(0, 0, Height - 1, Length - 2);

		Matrix<Type> Solution(Height);

		double Rate = 1;

		omp_set_num_threads(NUMBER_OF_CORES);

		#pragma omp parallel
		{
			#pragma omp for

			do {

				Matrix<Type> Vector = Solution;

				for (int FirstIndex = 0; FirstIndex < Height; FirstIndex++) {

					double Coefficient = 0;

					for (int SecondIndex = 0; SecondIndex < FirstIndex; SecondIndex++) {

						Coefficient += Temporary.GetMatrix()[FirstIndex][SecondIndex] * Solution.GetMatrix()[SecondIndex][0];
					}

					for (int SecondIndex = FirstIndex + 1; SecondIndex < Length - 1; SecondIndex++) {

						Coefficient += Temporary.GetMatrix()[FirstIndex][SecondIndex] * Vector.GetMatrix()[SecondIndex][0];
					}

					Solution.GetMatrix()[FirstIndex][0] = (SystemOfLinearEquations.GetMatrix()[FirstIndex][Length - 1] - Coefficient) / Temporary.GetMatrix()[FirstIndex][FirstIndex];
				}

				Rate = Solution.Subtract(Vector).Rate();

			} while (Rate > Precision);
		}

		return Solution;
	};
	Matrix<Type> Kaczmarz(const Matrix<Type>& SystemOfLinearEquations, const double Precision) const {

		CheckParameters(SystemOfLinearEquations);

		const int Length = SystemOfLinearEquations.GetLength();

		const int Height = SystemOfLinearEquations.GetHeight();

		Matrix<Type> Temporary = SystemOfLinearEquations.Cut(0, 0, Height - 1, Length - 2);

		Matrix<Type> Solution = Temporary.GetRow(0);

		double Rate = 1;

		int Counter = 0;

		do {

			Matrix<Type> Vector = Temporary.GetRow(Counter);

			double Coefficient = CheckDivisionByZero(SystemOfLinearEquations.GetMatrix()[Counter][Length - 1] - Vector.Scalar(Solution)) / (Vector.Rate() * Vector.Rate());

			Counter = Counter < Height - 1 ? Counter += 1 : 0;

			Rate = Solution.Add(Vector.Multiply(Coefficient)).Subtract(Solution).Rate();

			Solution += (Vector.Multiply(Coefficient));

		} while (Rate > Precision);

		Solution = Solution.Transpose();

		return Solution;
	};
	Matrix<Type> Jacobi(const Matrix<Type>& SystemOfLinearEquations, const double Precision) const {

		CheckParameters(SystemOfLinearEquations);

		const int Length = SystemOfLinearEquations.GetLength();

		const int Height = SystemOfLinearEquations.GetHeight();

		Matrix<Type> Temporary = SystemOfLinearEquations.Cut(0, 0, Height - 1, Length - 2);

		Matrix<Type> Solution(Height);

		double Rate = 1;

		omp_set_num_threads(NUMBER_OF_CORES);

		#pragma omp parallel
		{
			#pragma omp for

			do {

				Matrix<Type> Vector(Height);

				for (int FirstIndex = 0; FirstIndex < Height; FirstIndex++) {

					Vector.GetMatrix()[FirstIndex][0] = SystemOfLinearEquations.GetColumn(Length - 1).GetMatrix()[FirstIndex][0];

					for (int SecondIndex = 0; SecondIndex < Length - 1; SecondIndex++) {

						if (FirstIndex != SecondIndex) {

							Vector.GetMatrix()[FirstIndex][0] -= Temporary.GetMatrix()[FirstIndex][SecondIndex] * Solution.GetMatrix()[SecondIndex][0];
						}
					}

					Vector.GetMatrix()[FirstIndex][0] = CheckDivisionByZero(Vector.GetMatrix()[FirstIndex][0] / Temporary.GetMatrix()[FirstIndex][FirstIndex]);
				}

				Rate = fabs(Solution.GetMatrix()[0][0] - Vector.GetMatrix()[0][0]);

				for (int Index = 0; Index < Height; Index++) {

					if (fabs(Solution.GetMatrix()[Index][0] - Vector.GetMatrix()[Index][0]) > Rate) {

						Rate = fabs(Solution.GetMatrix()[Index][0] - Vector.GetMatrix()[Index][0]);
					}

					Solution.GetMatrix()[Index][0] = Vector.GetMatrix()[Index][0];
				}

			} while (Rate > Precision);
		}

		return Solution;
	};
	Matrix<Type> GaussJordan(const Matrix<Type>& SystemOfLinearEquations) const {

		CheckParameters(SystemOfLinearEquations);

		const int Length = SystemOfLinearEquations.GetLength();

		const int Height = SystemOfLinearEquations.GetHeight();

		Matrix<Type> Temporary(SystemOfLinearEquations);

		Matrix<Type> Solution(Height);

		omp_set_num_threads(NUMBER_OF_CORES);

		#pragma omp parallel
		{
			#pragma omp for

			for (int FirstIndex = 1; FirstIndex < Height; FirstIndex++) {

				for (int SecondIndex = FirstIndex; SecondIndex < Height; SecondIndex++) {

					double Coefficient = CheckDivisionByZero(Temporary.GetMatrix()[SecondIndex][FirstIndex - 1] / Temporary.GetMatrix()[FirstIndex - 1][FirstIndex - 1]);

					if (Coefficient != 0) {

						for (int AuxiliaryIndex = 0; AuxiliaryIndex < Length; AuxiliaryIndex++) {

							Temporary.GetMatrix()[SecondIndex][AuxiliaryIndex] -= Temporary.GetMatrix()[FirstIndex - 1][AuxiliaryIndex] * Coefficient;
						}
					}
				}
			}

			for (int FirstIndex = Height - 1; FirstIndex >= 0; FirstIndex--) {

				for (int SecondIndex = FirstIndex; SecondIndex > 0; SecondIndex--) {

					double Coefficient = CheckDivisionByZero(Temporary.GetMatrix()[SecondIndex - 1][FirstIndex] / Temporary.GetMatrix()[FirstIndex][FirstIndex]);

					if (Coefficient != 0) {

						for (int AuxiliaryIndex = Length - 1; AuxiliaryIndex >= FirstIndex; AuxiliaryIndex--) {

							Temporary.GetMatrix()[SecondIndex - 1][AuxiliaryIndex] -= Temporary.GetMatrix()[FirstIndex][AuxiliaryIndex] * Coefficient;
						}
					}
				}
			}

			for (int Index = 0; Index < Height; Index++) {

				Solution.GetMatrix()[Index][0] = CheckDivisionByZero(Temporary.GetMatrix()[Index][Length - 1] / Temporary.GetMatrix()[Index][Index]);
			}
		}

		return Solution;
	};
	Matrix<Type> Kramer(const Matrix<Type>& SystemOfLinearEquations) const {

		CheckParameters(SystemOfLinearEquations);

		const int Length = SystemOfLinearEquations.GetLength();

		const int Height = SystemOfLinearEquations.GetHeight();

		Matrix<Type> Temporary = SystemOfLinearEquations.Cut(0, 0, Height - 1, Length - 2);

		Matrix<Type> Solution(Height);

		double Delta = (double)Temporary.Determinant();

		for (int Index = 0; Index < Length - 1; Index++) {

			Matrix<Type> Column = SystemOfLinearEquations.GetColumn(Length - 1);

			Matrix<Type> Modified = Temporary.Replace(Index, Column);

			Solution.GetMatrix()[Index][0] = CheckDivisionByZero(Modified.Determinant() / Delta);
		}

		return Solution;
	};
	Matrix<Type> Inverse(const Matrix<Type>& SystemOfLinearEquations) const {

		CheckParameters(SystemOfLinearEquations);

		const int Length = SystemOfLinearEquations.GetLength();

		const int Height = SystemOfLinearEquations.GetHeight();

		Matrix<Type> Temporary = SystemOfLinearEquations.Cut(0, 0, Height - 1, Length - 2);

		Matrix<Type> Column = SystemOfLinearEquations.GetColumn(Length - 1);

		Matrix<Type> Solution(Height);

		Solution = Temporary.Inverse().Multiply(Column);

		return Solution;
	};
	Matrix<Type> Gauss(const Matrix<Type>& SystemOfLinearEquations) const {

		CheckParameters(SystemOfLinearEquations);

		const int Length = SystemOfLinearEquations.GetLength();

		const int Height = SystemOfLinearEquations.GetHeight();

		Matrix<Type> Temporary(SystemOfLinearEquations);

		Matrix<Type> Solution(Height);

		#pragma omp parallel
		{
			#pragma omp for

			for (int FirstIndex = 1; FirstIndex < Height; FirstIndex++) {

				for (int SecondIndex = FirstIndex; SecondIndex < Height; SecondIndex++) {

					double Coefficient = CheckDivisionByZero(Temporary.GetMatrix()[SecondIndex][FirstIndex - 1] / Temporary.GetMatrix()[FirstIndex - 1][FirstIndex - 1]);

					if (Coefficient != 0) {

						for (int AuxiliaryIndex = 0; AuxiliaryIndex < Length; AuxiliaryIndex++) {

							Temporary.GetMatrix()[SecondIndex][AuxiliaryIndex] -= Temporary.GetMatrix()[FirstIndex - 1][AuxiliaryIndex] * Coefficient;
						}
					}
				}
			}

			for (int FirstIndex = Height - 1; FirstIndex >= 0; FirstIndex--) {

				Solution.GetMatrix()[FirstIndex][0] = CheckDivisionByZero(Temporary.GetMatrix()[FirstIndex][Length - 1] / Temporary.GetMatrix()[FirstIndex][FirstIndex]);

				for (int SecondIndex = Height - 1; SecondIndex > FirstIndex; SecondIndex--) {

					Solution.GetMatrix()[FirstIndex][0] -= CheckDivisionByZero(Temporary.GetMatrix()[FirstIndex][SecondIndex] * Solution.GetMatrix()[SecondIndex][0] / Temporary.GetMatrix()[FirstIndex][FirstIndex]);
				}
			}
		}

		return Solution;
	};

	~LinearEquationSystem() = default;
	LinearEquationSystem() = default;
#pragma endregion
};

template<typename Type>
class Matrix {

#pragma region
private:
	void CheckSummationAndSubtraction(const Matrix& SomeMatrix) const {

		if (Length != SomeMatrix.Length || Height != SomeMatrix.Height) {

			throw "The length/height of the first matrix doesn't equal the length/height of the second matrix";
		}
	};
	void CheckMultiplication(const Matrix& SomeMatrix) const {

		if (Length != SomeMatrix.Height) {

			throw "The length of the first matrix doesn't equal the height of the second matrix";
		}
	};
	double CheckDivisionByZero(double Value) const {

		if (Value == -INFINITY || Value == INFINITY || Value == NAN || Value != Value) {

			return 0;
		}

		return Value;
	};
	void CheckSquareness() const {

		if (Length != Height) {

			throw "The length doesn't equal the height";
		}
	};
	void CheckParameters() const {

		if (Length == 0 || Height == 0) {

			throw "The length/height equals zero";
		}
	};

	void SwapRows(int FirstIndex, int SecondIndex) {

		Type* Temporary = InnerMatrix[FirstIndex];

		InnerMatrix[FirstIndex] = InnerMatrix[SecondIndex];

		InnerMatrix[SecondIndex] = Temporary;
	};

	Type** InnerMatrix;
	int Length;
	int Height;
#pragma endregion

#pragma region
public:
	template<typename Type>
	friend ostream& operator << (ostream& Stream, const Matrix<Type>& SomeMatrix) {

		const int Length = SomeMatrix.Length;

		const int Height = SomeMatrix.Height;

		for (int FirstIndex = 0; FirstIndex < Height; FirstIndex++) {

			for (int SecondIndex = 0; SecondIndex < Length; SecondIndex++) {

				Stream << setw(12) << SomeMatrix.InnerMatrix[FirstIndex][SecondIndex];
			}

			Stream << endl;
		}

		return Stream;
	};
	template<typename Type>
	friend istream& operator >> (istream& Stream, const Matrix<Type>& SomeMatrix) {

		const int Length = SomeMatrix.Length;

		const int Height = SomeMatrix.Height;

		for (int FirstIndex = 0; FirstIndex < Height; FirstIndex++) {

			for (int SecondIndex = 0; SecondIndex < Length; SecondIndex++) {

				Stream >> SomeMatrix.InnerMatrix[FirstIndex][SecondIndex];
			}
		}

		return Stream;
	};
	template<typename Type>
	Type** ToPointers(const Matrix<Type>& SomeMatrix) {

		const int Length = SomeMatrix.Length;

		const int Height = SomeMatrix.Height;

		Type** Temporary = new Type*[Height];

		omp_set_num_threads(NUMBER_OF_CORES);

		#pragma omp parallel
		{
			#pragma omp for

			for (int FirstIndex = 0; FirstIndex < Height; FirstIndex++) {

				Temporary[FirstIndex] = new Type[Length];

				for (int SecondIndex = 0; SecondIndex < Length; SecondIndex++) {

					Temporary[FirstIndex][SecondIndex] = SomeMatrix.InnerMatrix[FirstIndex][SecondIndex];
				}
			}
		}

		return Temporary;
	};
	template<typename Type>
	Matrix ToMatrix(Type** SomeMatrix) {

		const int Length = _msize(SomeMatrix[0]) / sizeof(SomeMatrix[0][0]);

		const int Height = _msize(SomeMatrix) / sizeof(SomeMatrix[0]);

		Matrix<Type> Temporary(Height, Length);

		omp_set_num_threads(NUMBER_OF_CORES);

		#pragma omp parallel
		{
			#pragma omp for

			for (int FirstIndex = 0; FirstIndex < Height; FirstIndex++) {

				for (int SecondIndex = 0; SecondIndex < Length; SecondIndex++) {

					Temporary.InnerMatrix[FirstIndex][SecondIndex] = SomeMatrix[FirstIndex][SecondIndex];
				}
			}
		}

		return Temporary;
	};

	Matrix Divide(const Matrix& SomeMatrix, const double Precision) const {

		return this->Multiply(SomeMatrix.Inverse(Precision));
	};
	Matrix Subtract(const Matrix& SomeMatrix) const {

		CheckSummationAndSubtraction(SomeMatrix);

		const int Length = this->Length;

		const int Height = this->Height;

		Matrix<Type> Temporary(Height, Length);

		omp_set_num_threads(NUMBER_OF_CORES);

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
	Matrix Multiply(const Matrix& SomeMatrix) const {

		CheckMultiplication(SomeMatrix);

		const int ThatLength = SomeMatrix.Length;

		const int ThatHeight = SomeMatrix.Height;

		const int ThisLength = Length;

		const int ThisHeight = Height;

		Matrix<Type> Temporary(ThisHeight, ThatLength);

		Type* ThatColumn = new Type[ThatHeight];

		omp_set_num_threads(NUMBER_OF_CORES);

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
	Matrix Product(const Matrix& SomeMatrix) const {

		CheckSummationAndSubtraction(SomeMatrix);

		const int Length = this->Length;

		const int Height = this->Height;

		Matrix<Type> Temporary(Height, Length);

		omp_set_num_threads(NUMBER_OF_CORES);

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
	Matrix Add(const Matrix& SomeMatrix) const {

		CheckSummationAndSubtraction(SomeMatrix);

		const int Length = SomeMatrix.Length;

		const int Height = SomeMatrix.Height;

		Matrix<Type> Temporary(Height, Length);

		omp_set_num_threads(NUMBER_OF_CORES);

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
	Type Scalar(const Matrix& SomeMatrix) {

		CheckSummationAndSubtraction(SomeMatrix);

		const int Length = this->Length;

		const int Height = this->Height;

		Type Scalar = 0;

		omp_set_num_threads(NUMBER_OF_CORES);

		#pragma omp parallel
		{
			#pragma omp for

			for (int FirstIndex = 0; FirstIndex < Height; FirstIndex++) {

				for (int SecondIndex = 0; SecondIndex < Length; SecondIndex++) {

					Scalar += InnerMatrix[FirstIndex][SecondIndex] * SomeMatrix.InnerMatrix[FirstIndex][SecondIndex];
				}
			}
		}

		return Scalar;
	};
	Matrix Multiply(Type Coefficient) const {

		CheckParameters();

		const int Length = this->Length;

		const int Height = this->Height;

		Matrix<Type> Temporary(Height, Length);

		omp_set_num_threads(NUMBER_OF_CORES);

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
	Type Determinant() const {

		CheckSquareness();

		const int Length = this->Length;

		const int Height = this->Height;

		Matrix<Type> Temporary(*this);

		double Value = 1;

		omp_set_num_threads(NUMBER_OF_CORES);

		#pragma omp parallel
		{
			#pragma omp for

			for (int FirstIndex = 1; FirstIndex < Height; FirstIndex++) {

				for (int SecondIndex = FirstIndex; SecondIndex < Height; SecondIndex++) {

					double Coefficient = CheckDivisionByZero(Temporary.InnerMatrix[SecondIndex][FirstIndex - 1] / Temporary.InnerMatrix[FirstIndex - 1][FirstIndex - 1]);

					if (Coefficient != 0) {

						for (int AuxiliaryIndex = 0; AuxiliaryIndex < Length; AuxiliaryIndex++) {

							Temporary.InnerMatrix[SecondIndex][AuxiliaryIndex] -= Temporary.InnerMatrix[FirstIndex - 1][AuxiliaryIndex] * Coefficient;
						}
					}
				}
			}

			for (int Index = 0; Index < Height; Index++) {

				Value *= Temporary.InnerMatrix[Index][Index];
			}
		}

		return Value;
	};
	Matrix Inverse() const {

		CheckSquareness();

		const int Length = this->Length;

		const int Height = this->Height;

		Matrix<Type> Solution(Height, Length);

		Matrix<Type> Temporary(*this);

		Solution.UnitFilling();

		omp_set_num_threads(NUMBER_OF_CORES);

		#pragma omp parallel
		{
			#pragma omp for

			for (int FirstIndex = 1; FirstIndex < Height; FirstIndex++) {

				for (int SecondIndex = FirstIndex; SecondIndex < Height; SecondIndex++) {

					double Coefficient = CheckDivisionByZero(Temporary.InnerMatrix[SecondIndex][FirstIndex - 1] / Temporary.InnerMatrix[FirstIndex - 1][FirstIndex - 1]);

					if (Coefficient != 0) {

						for (int AuxiliaryIndex = 0; AuxiliaryIndex < Length; AuxiliaryIndex++) {

							Temporary.InnerMatrix[SecondIndex][AuxiliaryIndex] -= Temporary.InnerMatrix[FirstIndex - 1][AuxiliaryIndex] * Coefficient;

							Solution.InnerMatrix[SecondIndex][AuxiliaryIndex] -= Solution.InnerMatrix[FirstIndex - 1][AuxiliaryIndex] * Coefficient;
						}
					}
				}
			}

			for (int FirstIndex = 0; FirstIndex < Height; FirstIndex++) {

				double Coefficient = CheckDivisionByZero(1 / Temporary.InnerMatrix[FirstIndex][FirstIndex]);

				for (int SecondIndex = 0; SecondIndex < Length; SecondIndex++) {

					if (Coefficient != 1) {

						Temporary.InnerMatrix[FirstIndex][SecondIndex] *= Coefficient;

						Solution.InnerMatrix[FirstIndex][SecondIndex] *= Coefficient;
					}
				}
			}

			for (int FirstIndex = Height - 1; FirstIndex >= 0; FirstIndex--) {

				for (int SecondIndex = FirstIndex; SecondIndex > 0; SecondIndex--) {

					double Coefficient = CheckDivisionByZero(Temporary.InnerMatrix[SecondIndex - 1][FirstIndex] / Temporary.InnerMatrix[FirstIndex][FirstIndex]);

					if (Coefficient != 0) {

						for (int AuxiliaryIndex = Length - 1; AuxiliaryIndex >= 0; AuxiliaryIndex--) {

							Temporary.InnerMatrix[SecondIndex - 1][AuxiliaryIndex] -= Temporary.InnerMatrix[FirstIndex][AuxiliaryIndex] * Coefficient;

							Solution.InnerMatrix[SecondIndex - 1][AuxiliaryIndex] -= Solution.InnerMatrix[FirstIndex][AuxiliaryIndex] * Coefficient;
						}
					}
				}
			}
		}

		return Solution;
	};

	Matrix Replace(int Position, const Matrix& SomeVector) const {

		CheckParameters();

		const int Length = this->Length;

		const int Height = this->Height;

		Matrix<Type> Temporary(*this);

		omp_set_num_threads(NUMBER_OF_CORES);

		#pragma omp parallel
		{
			#pragma omp for

			if (SomeVector.Length == 1) {

				for (int Index = 0; Index < Height; Index++) {

					Temporary.InnerMatrix[Index][Position] = SomeVector.InnerMatrix[Index][0];
				}

				return Temporary;
			}

			for (int Index = 0; Index < Length; Index++) {

				Temporary.InnerMatrix[Position][Index] = SomeVector.InnerMatrix[0][Index];
			}
		}

		return Temporary;
	};
	Matrix Cut(int StartX, int StartY, int EndX, int EndY) const {

		CheckParameters();

		const int Height = (EndX - StartX) + 1;

		const int Length = (EndY - StartY) + 1;

		Matrix<Type> Temporary(Height, Length);

		omp_set_num_threads(NUMBER_OF_CORES);

		#pragma omp parallel
		{
			#pragma omp for

			for (int FirstIndex = StartY; FirstIndex <= EndY; FirstIndex++) {

				for (int SecondIndex = StartX; SecondIndex <= EndX; SecondIndex++) {

					Temporary.InnerMatrix[FirstIndex - StartY][SecondIndex - StartX] = InnerMatrix[FirstIndex][SecondIndex];
				}
			}
		}

		return Temporary;
	};
	Matrix LoadFromFIle(string Path, char Delimiter) const {

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
	Matrix UseFunction(Type(*Function)(Type)) const {

		CheckParameters();

		const int Length = this->Length;

		const int Height = this->Height;

		Matrix<Type> Temporary(Height, Length);

		omp_set_num_threads(NUMBER_OF_CORES);

		#pragma omp parallel
		{
			#pragma omp for

			for (int FirstIndex = 0; FirstIndex < Height; FirstIndex++) {

				for (int SecondIndex = 0; SecondIndex < Length; SecondIndex++) {

					Temporary.InnerMatrix[FirstIndex][SecondIndex] = Function(InnerMatrix[FirstIndex][SecondIndex]);
				}
			}
		}

		return Temporary;
	};
	Matrix Reshape(int Height, int Length) const {

		CheckParameters();

		const int ThisLength = this->Length;

		const int ThisHeight = this->Height;

		Type* DeconvolutedMatrix = new Type[ThisLength * ThisHeight];

		Matrix<Type> Temporary(Height, Length);

		omp_set_num_threads(NUMBER_OF_CORES);

		#pragma omp parallel
		{
			#pragma omp for

			int InnerHeight = 0;

			for (int FirstIndex = 0; FirstIndex < ThisHeight; FirstIndex++) {

				for (int SecondIndex = 0; SecondIndex < ThisLength; SecondIndex++) {

					DeconvolutedMatrix[InnerHeight] = InnerMatrix[FirstIndex][SecondIndex];

					InnerHeight++;
				}
			}

			InnerHeight = 0;

			int Counter = 0;

			for (int Index = 0; Index < ThisLength * ThisHeight; Index++) {

				if (Counter == Length) {

					Counter = 0;

					InnerHeight++;
				}

				Temporary.InnerMatrix[InnerHeight][Counter] = DeconvolutedMatrix[Index];

				Counter++;
			}
		}

		delete[] DeconvolutedMatrix;

		return Temporary;
	};
	Matrix Exponentiate(int Power) const {

		CheckParameters();

		Matrix<Type> Temporary(Height, Length);

		for (int Index = 0; Index < Power; Index++) {

			Temporary = Temporary.Add(this->Multiply(*this));
		}

		return Temporary;
	};
	Matrix Fill(Type Value) const {

		CheckParameters();

		const int Length = this->Length;

		const int Height = this->Height;

		omp_set_num_threads(NUMBER_OF_CORES);

		#pragma omp parallel
		{
			#pragma omp for

			for (int FirstIndex = 0; FirstIndex < Height; FirstIndex++) {

				for (int SecondIndex = 0; SecondIndex < Length; SecondIndex++) {

					InnerMatrix[FirstIndex][SecondIndex] = Value;
				}
			}
		}

		return *this;
	};
	Matrix Randomize() const {

		CheckParameters();

		const int Length = this->Length;

		const int Height = this->Height;

		srand(1 + rand() % 100);

		omp_set_num_threads(NUMBER_OF_CORES);

		#pragma omp parallel
		{
			#pragma omp for

			for (int FirstIndex = 0; FirstIndex < Height; FirstIndex++) {

				for (int SecondIndex = 0; SecondIndex < Length; SecondIndex++) {

					InnerMatrix[FirstIndex][SecondIndex] = static_cast <float> (rand()) / static_cast <float> (RAND_MAX) - 0.5f;
				}
			}
		}

		return *this;
	};
	Matrix UnitFilling() const {

		CheckSquareness();

		const int Height = this->Height;

		for (int Index = 0; Index < Height; Index++) {

			InnerMatrix[Index][Index] = 1;
		}

		return *this;
	};
	Matrix Transpose() const {

		CheckParameters();

		const int Length = this->Length;

		const int Height = this->Height;

		Matrix<Type> Temporary(Length, Height);

		omp_set_num_threads(NUMBER_OF_CORES);

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
	Matrix Clear() const {

		CheckParameters();

		const int Length = this->Length;

		const int Height = this->Height;

		omp_set_num_threads(NUMBER_OF_CORES);

		#pragma omp parallel
		{
			#pragma omp for

			for (int FirstIndex = 0; FirstIndex < Height; FirstIndex++) {

				for (int SecondIndex = 0; SecondIndex < Length; SecondIndex++) {

					InnerMatrix[FirstIndex][SecondIndex] = 0;
				}
			}
		}

		return *this;
	};
	Matrix Print() const {

		CheckParameters();

		const int Length = this->Length;

		const int Height = this->Height;

		omp_set_num_threads(NUMBER_OF_CORES);

		#pragma omp parallel
		{
			#pragma omp for

			for (int FirstIndex = 0; FirstIndex < Height; FirstIndex++) {

				for (int SecondIndex = 0; SecondIndex < Length; SecondIndex++) {

					cout << setw(12) << InnerMatrix[FirstIndex][SecondIndex];
				}

				cout << endl;
			}
		}

		cout << endl;

		return *this;
	};

	Matrix operator * (const Matrix& SomeMatrix) const {

		return this->Multiply(SomeMatrix);
	};
	Matrix operator + (const Matrix& SomeMatrix) const {

		return this->Add(SomeMatrix);
	};
	Matrix operator / (const Matrix& SomeMatrix) const {

		return this->Divide(SomeMatrix, 0.0001);
	};
	Matrix operator - (const Matrix& SomeMatrix) const {

		return this->Subtract(SomeMatrix);
	};
	Matrix operator * (const Type Coefficient) const {

		return this->Multiply(Coefficient);
	};
	Matrix operator ^ (const int Power) const {

		return this->Exponentiate(Power);
	};

	Matrix& operator *= (const Matrix& SomeMatrix) {

		return *this = this->Multiply(SomeMatrix);
	};
	Matrix& operator += (const Matrix& SomeMatrix) {

		return *this = this->Add(SomeMatrix);
	};
	Matrix& operator /= (const Matrix& SomeMatrix) {

		return *this = this->Divide(SomeMatrix, 0.0001);
	};
	Matrix& operator -= (const Matrix& SomeMatrix) {

		return *this = this->Subtract(SomeMatrix);
	};
	Matrix& operator = (const Matrix& SomeMatrix) {

		if (this == &SomeMatrix) {

			return *this;
		}

		if (Length != SomeMatrix.Length || Height != SomeMatrix.Height) {

			this->~Matrix();

			Length = SomeMatrix.Length;

			Height = SomeMatrix.Height;

			InnerMatrix = new Type*[Height];

			for (int FirstIndex = 0; FirstIndex < Height; FirstIndex++) {

				InnerMatrix[FirstIndex] = new Type[Length];

				for (int SecondIndex = 0; SecondIndex < Length; SecondIndex++) {

					InnerMatrix[FirstIndex][SecondIndex] = SomeMatrix.InnerMatrix[FirstIndex][SecondIndex];
				}
			}
		}

		for (int FirstIndex = 0; FirstIndex < Height; FirstIndex++) {

			for (int SecondIndex = 0; SecondIndex < Length; SecondIndex++) {

				InnerMatrix[FirstIndex][SecondIndex] = SomeMatrix.InnerMatrix[FirstIndex][SecondIndex];
			}
		}

		return *this;
	};
	Matrix& operator *= (const Type Coefficient) {

		return *this = this->Multiply(Coefficient);
	};
	Matrix& operator ^= (const int Power) {

		return *this = this->Exponentiate(Power);
	};

	tuple<Matrix, Matrix> LUPDecomposition() const {

		CheckSquareness();

		const int Length = this->Length;

		const int Height = this->Height;

		Matrix<Type> Diagonal(Height, Length);

		Matrix<Type> Temporary(*this);

		Diagonal.UnitFilling();

		omp_set_num_threads(NUMBER_OF_CORES);

		#pragma omp parallel
		{
			#pragma omp for

			for (int FirstIndex = 0; FirstIndex < Height; FirstIndex++) {

				double PivotValue = 0;

				int Pivot = -1;

				for (int SecondIndex = FirstIndex; SecondIndex < Height; SecondIndex++) {

					if (fabs(Temporary.InnerMatrix[SecondIndex][FirstIndex]) > PivotValue) {

						PivotValue = fabs(Temporary.InnerMatrix[SecondIndex][FirstIndex]);

						Pivot = SecondIndex;
					}
				}

				if (PivotValue != 0) {

					Diagonal.SwapRows(Pivot, FirstIndex);

					Temporary.SwapRows(Pivot, FirstIndex);

					for (int SecondIndex = FirstIndex + 1; SecondIndex < Height; SecondIndex++) {

						Temporary.InnerMatrix[SecondIndex][FirstIndex] = CheckDivisionByZero(Temporary.InnerMatrix[SecondIndex][FirstIndex] / Temporary.InnerMatrix[FirstIndex][FirstIndex]);

						for (int AuxiliaryIndex = FirstIndex + 1; AuxiliaryIndex < Length; AuxiliaryIndex++) {

							Temporary.InnerMatrix[SecondIndex][AuxiliaryIndex] -= Temporary.InnerMatrix[SecondIndex][FirstIndex] * Temporary.InnerMatrix[FirstIndex][AuxiliaryIndex];
						}
					}
				}
			}
		}

		return tuple<Matrix, Matrix>(Temporary, Diagonal);
	};
	tuple<Matrix, Matrix> LUDecomposition() const {

		CheckSquareness();

		const int Length = this->Length;

		const int Height = this->Height;

		Matrix<Type> Modified(Height, Length);

		Matrix<Type> Temporary(*this);

		omp_set_num_threads(NUMBER_OF_CORES);

		#pragma omp parallel
		{
			#pragma omp for

			for (int FirstIndex = 0; FirstIndex < Height; FirstIndex++) {

				for (int SecondIndex = FirstIndex; SecondIndex < Length; SecondIndex++) {

					Modified.InnerMatrix[SecondIndex][FirstIndex] = CheckDivisionByZero(Temporary.InnerMatrix[SecondIndex][FirstIndex] / Temporary.InnerMatrix[FirstIndex][FirstIndex]);
				}
			}

			for (int FirstIndex = 1; FirstIndex < Height; FirstIndex++) {

				for (int SecondIndex = FirstIndex - 1; SecondIndex < Height; SecondIndex++) {

					for (int AuxiliaryIndex = SecondIndex; AuxiliaryIndex < Length; AuxiliaryIndex++) {

						Modified.InnerMatrix[AuxiliaryIndex][SecondIndex] = CheckDivisionByZero(Temporary.InnerMatrix[AuxiliaryIndex][SecondIndex] / Temporary.InnerMatrix[SecondIndex][SecondIndex]);
					}
				}

				for (int SecondIndex = FirstIndex; SecondIndex < Height; SecondIndex++) {

					for (int AuxiliaryIndex = FirstIndex - 1; AuxiliaryIndex < Length; AuxiliaryIndex++) {

						Temporary.InnerMatrix[SecondIndex][AuxiliaryIndex] -= Modified.InnerMatrix[SecondIndex][FirstIndex - 1] * Temporary.InnerMatrix[FirstIndex - 1][AuxiliaryIndex];
					}
				}
			}
		}

		return tuple<Matrix, Matrix>(Temporary, Modified);
	};
	Matrix CholeskyDecomposition() const {

		CheckSquareness();

		const int Length = this->Length;

		const int Height = this->Height;

		Matrix<Type> Temporary(Height, Length);

		omp_set_num_threads(NUMBER_OF_CORES);

		#pragma omp parallel
		{
			#pragma omp for

			for (int FirstIndex = 0; FirstIndex < Height; FirstIndex++) {

				double Coefficient = 0;

				for (int SecondIndex = 0; SecondIndex < FirstIndex; SecondIndex++) {

					Coefficient = 0;

					for (int AuxiliaryIndex = 0; AuxiliaryIndex < SecondIndex; AuxiliaryIndex++) {

						Coefficient += Temporary.InnerMatrix[FirstIndex][AuxiliaryIndex] * Temporary.InnerMatrix[SecondIndex][AuxiliaryIndex];
					}

					Temporary.InnerMatrix[FirstIndex][SecondIndex] = CheckDivisionByZero((InnerMatrix[FirstIndex][SecondIndex] - Coefficient) / Temporary.InnerMatrix[SecondIndex][SecondIndex]);
				}

				Coefficient = InnerMatrix[FirstIndex][FirstIndex];

				for (int AuxiliaryIndex = 0; AuxiliaryIndex < FirstIndex; AuxiliaryIndex++) {

					Coefficient -= Temporary.InnerMatrix[FirstIndex][AuxiliaryIndex] * Temporary.InnerMatrix[FirstIndex][AuxiliaryIndex];
				}

				Temporary.InnerMatrix[FirstIndex][FirstIndex] = sqrt(Coefficient);
			}
		}

		return Temporary;
	};

	Matrix GetColumn(int Position) const {

		CheckParameters();

		const int Height = this->Height;

		Matrix<Type> Temporary(Height);

		for (int Index = 0; Index < Height; Index++) {

			Temporary.InnerMatrix[Index][0] = InnerMatrix[Index][Position];
		}

		return Temporary;
	};
	Matrix GetRow(int Position) const {

		CheckParameters();

		const int Length = this->Length;

		Matrix<Type> Temporary(1, Length);

		for (int Index = 0; Index < Length; Index++) {

			Temporary.InnerMatrix[0][Index] = InnerMatrix[Position][Index];
		}

		return Temporary;
	};
	Matrix GetMainDiagonal() const {

		CheckParameters();

		const int Height = this->Length;

		Matrix<Type> Temporary(Height);

		for (int Index = 0; Index < Height; Index++) {

			Temporary.InnerMatrix[Index][0] = InnerMatrix[Index][Index];
		}

		return Temporary;
	};
	Matrix GetSideDiagonal() const {

		CheckParameters();

		const int Height = this->Height;

		Matrix<Type> Temporary(Height);

		for (int Index = 0; Index < Height; Index++) {

			Temporary.InnerMatrix[Index][0] = InnerMatrix[Index][(Height - Index) - 1];
		}

		return Temporary;
	};
	Type** GetMatrix() const {

		return InnerMatrix;
	};
	int GetLength() const {

		return Length;
	};
	int GetHeight() const {

		return Height;
	};

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

		delete InnerMatrix;
	};
	Matrix() : Length(0), Height(0) {

		InnerMatrix = nullptr;
	};

	Type ColumnSum(int Position) {

		Type Value = 0;

		for (int Index = 0; Index < Height; Index++) {

			Value += InnerMatrix[Index][Position];
		}

		return Value;
	};
	Type RowSum(int Position) {

		Type Value = 0;

		for (int Index = 0; Index < Length; Index++) {

			Value += InnerMatrix[Position][Index];
		}

		return Value;
	};
	Type Average() {

		CheckParameters();

		const int Length = this->Length;

		const int Height = this->Height;

		Type Value = 0;

		omp_set_num_threads(NUMBER_OF_CORES);

		#pragma omp parallel
		{
			#pragma omp for

			for (int FirstIndex = 0; FirstIndex < Height; FirstIndex++) {

				for (int SecondIndex = 0; SecondIndex < Length; SecondIndex++) {

					Value += InnerMatrix[FirstIndex][SecondIndex];
				}
			}
		}

		return Value / (Height * Length);
	}
	Type Rate() {

		CheckParameters();

		const int Length = this->Length;

		const int Height = this->Height;

		Type Value = 0;

		omp_set_num_threads(NUMBER_OF_CORES);

		#pragma omp parallel
		{
			#pragma omp for

			for (int FirstIndex = 0; FirstIndex < Height; FirstIndex++) {

				for (int SecondIndex = 0; SecondIndex < Length; SecondIndex++) {

					Value += pow(InnerMatrix[FirstIndex][SecondIndex], 2);
				}
			}
		}

		return sqrt(Value);
	};
	Type Max() {

		CheckParameters();

		const int Length = this->Length;

		const int Height = this->Height;

		Type Value = 0;

		omp_set_num_threads(NUMBER_OF_CORES);

		#pragma omp parallel
		{
			#pragma omp for

			for (int FirstIndex = 0; FirstIndex < Height; FirstIndex++) {

				for (int SecondIndex = 0; SecondIndex < Length; SecondIndex++) {

					Value = fmax(Value, InnerMatrix[FirstIndex][SecondIndex]);
				}
			}
		}

		return Value;
	};
#pragma endregion
};
