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
#include <list>
#pragma once

using namespace std;

const int MAX_CORES = omp_get_max_threads();

template<typename SomeType>
class OMPMatrix {

#pragma region
private:
	inline void CheckSummationAndSubtraction(const OMPMatrix& SomeMatrix) const {

		if (Length != SomeMatrix.Length || Height != SomeMatrix.Height) {

			throw "The length/height of the first matrix doesn't equal the length/height of the second matrix";
		}
	};
	inline void CheckMultiplication(const OMPMatrix& SomeMatrix) const {

		if (Length != SomeMatrix.Height) {

			throw "The length of the first matrix doesn't equal the height of the second matrix";
		}
	};
	inline double CheckDivisionByZero(double Value) const {

		if (Value == -INFINITY || Value == INFINITY || Value == NAN || Value != Value) {

			return 0;
		}

		return Value;
	};
	inline void CheckSquareness() const {

		if (Length != Height) {

			throw "The length doesn't equal the height";
		}
	};
	inline void CheckParameters() const {

		if (Length == 0 || Height == 0) {

			throw "The length/height equals zero";
		}
	};

	inline void SwapRows(int FirstIndex, int SecondIndex) {

		SomeType* Temporary = InnerMatrix[FirstIndex];

		InnerMatrix[FirstIndex] = InnerMatrix[SecondIndex];

		InnerMatrix[SecondIndex] = Temporary;
	};

	SomeType** InnerMatrix;
	int NUMBER_OF_CORES;
	int Length;
	int Height;
#pragma endregion

#pragma region
public:
	template<typename SomeType>
	friend ostream& operator << (ostream& Stream, const OMPMatrix<SomeType>& SomeMatrix) {

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
	template<typename SomeType>
	friend istream& operator >> (istream& Stream, const OMPMatrix<SomeType>& SomeMatrix) {

		const int Length = SomeMatrix.Length;

		const int Height = SomeMatrix.Height;

		for (int FirstIndex = 0; FirstIndex < Height; FirstIndex++) {

			for (int SecondIndex = 0; SecondIndex < Length; SecondIndex++) {

				Stream >> SomeMatrix.InnerMatrix[FirstIndex][SecondIndex];
			}
		}

		return Stream;
	};

	OMPMatrix Replace(int Position, const OMPMatrix& SomeVector) const {

		CheckParameters();

		const int Length = this->Length;

		const int Height = this->Height;

		OMPMatrix<SomeType> Temporary(*this);

		if (SomeVector.Length == 1) {

			for (int Index = 0; Index < Height; Index++) {

				Temporary.InnerMatrix[Index][Position] = SomeVector.InnerMatrix[Index][0];
			}

			return Temporary;
		}

		for (int Index = 0; Index < Length; Index++) {

			Temporary.InnerMatrix[Position][Index] = SomeVector.InnerMatrix[0][Index];
		}

		return Temporary;
	};
	OMPMatrix UseFunction(SomeType(*Function)(SomeType)) const {

		CheckParameters();

		const int Length = this->Length;

		const int Height = this->Height;

		OMPMatrix<SomeType> Temporary(Height, Length);

		for (int FirstIndex = 0; FirstIndex < Height; FirstIndex++) {

			for (int SecondIndex = 0; SecondIndex < Length; SecondIndex++) {

				Temporary.InnerMatrix[FirstIndex][SecondIndex] = Function(InnerMatrix[FirstIndex][SecondIndex]);
			}
		}

		return Temporary;
	};
	OMPMatrix Cut(int StartX, int StartY, int EndX, int EndY) const {

		CheckParameters();

		const int Height = (EndX - StartX) + 1;

		const int Length = (EndY - StartY) + 1;

		OMPMatrix<SomeType> Temporary(Height, Length);

		for (int FirstIndex = StartY; FirstIndex < Height; FirstIndex++) {

			for (int SecondIndex = StartX; SecondIndex < Length; SecondIndex++) {

				Temporary.InnerMatrix[FirstIndex - StartY][SecondIndex - StartX] = InnerMatrix[FirstIndex][SecondIndex];
			}
		}

		return Temporary;
	};
	OMPMatrix LoadFromFIle(string Path, char Delimiter) const {

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

		OMPMatrix<SomeType> Temporary(Height, Length);

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
	OMPMatrix Reshape(int Height, int Length) const {

		CheckParameters();

		const int ThisLength = this->Length;

		const int ThisHeight = this->Height;

		SomeType* DeconvolutedMatrix = new SomeType[ThisLength * ThisHeight];

		OMPMatrix<SomeType> Temporary(Height, Length);

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

		delete[] DeconvolutedMatrix;

		return Temporary;
	};
	OMPMatrix Exponentiate(int Power) const {

		CheckParameters();

		OMPMatrix<SomeType> Temporary(*this);

		for (int Index = 1; Index < Power; Index++) {

			Temporary = Temporary.Multiply(*this);
		}

		return Temporary;
	};
	OMPMatrix Fill(SomeType Value) const {

		CheckParameters();

		const int Length = this->Length;

		const int Height = this->Height;

		for (int FirstIndex = 0; FirstIndex < Height; FirstIndex++) {

			for (int SecondIndex = 0; SecondIndex < Length; SecondIndex++) {

				InnerMatrix[FirstIndex][SecondIndex] = Value;
			}
		}

		return *this;
	};
	OMPMatrix Randomize() const {

		CheckParameters();

		const int Length = this->Length;

		const int Height = this->Height;

		srand(1 + rand() % 100);

		for (int FirstIndex = 0; FirstIndex < Height; FirstIndex++) {

			for (int SecondIndex = 0; SecondIndex < Length; SecondIndex++) {

				InnerMatrix[FirstIndex][SecondIndex] = static_cast <float> (rand()) / static_cast <float> (RAND_MAX) - 0.5f;
			}
		}

		return *this;
	};
	OMPMatrix UnitFilling() const {

		CheckSquareness();

		const int Height = this->Height;

		for (int Index = 0; Index < Height; Index++) {

			InnerMatrix[Index][Index] = 1;
		}

		return *this;
	};
	OMPMatrix Clear() const {

		CheckParameters();

		const int Length = this->Length;

		const int Height = this->Height;

		for (int FirstIndex = 0; FirstIndex < Height; FirstIndex++) {

			for (int SecondIndex = 0; SecondIndex < Length; SecondIndex++) {

				InnerMatrix[FirstIndex][SecondIndex] = 0;
			}
		}

		return *this;
	};
	OMPMatrix Print() const {

		CheckParameters();

		const int Length = this->Length;

		const int Height = this->Height;

		for (int FirstIndex = 0; FirstIndex < Height; FirstIndex++) {

			for (int SecondIndex = 0; SecondIndex < Length; SecondIndex++) {

				cout << setw(12) << InnerMatrix[FirstIndex][SecondIndex];
			}

			cout << endl;
		}

		cout << endl;

		return *this;
	};

	OMPMatrix operator * (const OMPMatrix& SomeMatrix) const {

		return this->Multiply(SomeMatrix);
	};
	OMPMatrix operator + (const OMPMatrix& SomeMatrix) const {

		return this->Add(SomeMatrix);
	};
	OMPMatrix operator / (const OMPMatrix& SomeMatrix) const {

		return this->Divide(SomeMatrix);
	};
	OMPMatrix operator - (const OMPMatrix& SomeMatrix) const {

		return this->Subtract(SomeMatrix);
	};
	OMPMatrix operator * (const SomeType Coefficient) const {

		return this->Multiply(Coefficient);
	};
	OMPMatrix operator ^ (const int Power) const {

		return this->Exponentiate(Power);
	};
	OMPMatrix operator ~ () {

		return this->Transpose();
	};
	
	OMPMatrix(int Height, int Length, long NUMBER_OF_CORES) : Length(Length), Height(Height), NUMBER_OF_CORES(NUMBER_OF_CORES) {

		const int ThisLength = Length;

		const int ThisHeight = Height;

		InnerMatrix = new SomeType*[ThisHeight];

		omp_set_num_threads(NUMBER_OF_CORES);

		#pragma omp parallel for schedule(static, NUMBER_OF_CORES)

		for (int FirstIndex = 0; FirstIndex < ThisHeight; FirstIndex++) {

			InnerMatrix[FirstIndex] = new SomeType[ThisLength];

			for (int SecondIndex = 0; SecondIndex < ThisLength; SecondIndex++) {

				InnerMatrix[FirstIndex][SecondIndex] = 0;
			}
		}
	};
	OMPMatrix(int Height, long NUMBER_OF_CORES) : Length(1), Height(Height), NUMBER_OF_CORES(NUMBER_OF_CORES) {

		const int ThisHeight = Height;

		InnerMatrix = new SomeType*[ThisHeight];

		omp_set_num_threads(NUMBER_OF_CORES);

		#pragma omp parallel for schedule(static, NUMBER_OF_CORES)

		for (int Index = 0; Index < ThisHeight; Index++) {

			InnerMatrix[Index] = new SomeType[1];

			InnerMatrix[Index][0] = 0;
		}
	};
	OMPMatrix(const OMPMatrix& SomeMatrix) : Length(SomeMatrix.Length), Height(SomeMatrix.Height), NUMBER_OF_CORES(SomeMatrix.NUMBER_OF_CORES) {

		const int ThisLength = Length;

		const int ThisHeight = Height;

		InnerMatrix = new SomeType*[ThisHeight];

		omp_set_num_threads(NUMBER_OF_CORES);

		#pragma omp parallel for schedule(static, NUMBER_OF_CORES)

		for (int FirstIndex = 0; FirstIndex < ThisHeight; FirstIndex++) {

			InnerMatrix[FirstIndex] = new SomeType[ThisLength];

			for (int SecondIndex = 0; SecondIndex < ThisLength; SecondIndex++) {

				InnerMatrix[FirstIndex][SecondIndex] = SomeMatrix.InnerMatrix[FirstIndex][SecondIndex];
			}
		}
	};
	OMPMatrix(int Height, int Length) : Length(Length), Height(Height), NUMBER_OF_CORES(MAX_CORES) {

		const int ThisLength = Length;

		const int ThisHeight = Height;

		InnerMatrix = new SomeType*[ThisHeight];

		omp_set_num_threads(NUMBER_OF_CORES);

		#pragma omp parallel for schedule(static, NUMBER_OF_CORES)

		for (int FirstIndex = 0; FirstIndex < ThisHeight; FirstIndex++) {

			InnerMatrix[FirstIndex] = new SomeType[ThisLength];

			for (int SecondIndex = 0; SecondIndex < ThisLength; SecondIndex++) {

				InnerMatrix[FirstIndex][SecondIndex] = 0;
			}
		}
	};
	OMPMatrix(int Height) : Length(1), Height(Height), NUMBER_OF_CORES(MAX_CORES) {

		const int ThisHeight = Height;

		InnerMatrix = new SomeType*[ThisHeight];

		omp_set_num_threads(NUMBER_OF_CORES);

		#pragma omp parallel for schedule(static, NUMBER_OF_CORES)

		for (int Index = 0; Index < ThisHeight; Index++) {

			InnerMatrix[Index] = new SomeType[1];

			InnerMatrix[Index][0] = 0;
		}
	};
	~OMPMatrix() {

		const int ThisHeight = Height;

		omp_set_num_threads(NUMBER_OF_CORES);

		#pragma omp parallel for schedule(static, NUMBER_OF_CORES)

		for (int Index = 0; Index < ThisHeight; Index++) {

			delete[] InnerMatrix[Index];
		}

		delete[] InnerMatrix;
	};
	OMPMatrix() : Length(0), Height(0) {

		InnerMatrix = nullptr;
	};

	OMPMatrix Subtract(const OMPMatrix& SomeMatrix) const {

		CheckSummationAndSubtraction(SomeMatrix);

		const int Length = this->Length;

		const int Height = this->Height;

		OMPMatrix<SomeType> Temporary(Height, Length);

		omp_set_num_threads(NUMBER_OF_CORES);

		#pragma omp parallel for schedule(static, NUMBER_OF_CORES)

		for (int FirstIndex = 0; FirstIndex < Height; FirstIndex++) {

			for (int SecondIndex = 0; SecondIndex < Length; SecondIndex++) {

				Temporary.InnerMatrix[FirstIndex][SecondIndex] = InnerMatrix[FirstIndex][SecondIndex] - SomeMatrix.InnerMatrix[FirstIndex][SecondIndex];
			}
		}

		return Temporary;
	};
	OMPMatrix Multiply(const OMPMatrix& SomeMatrix) const {

		CheckMultiplication(SomeMatrix);

		const int ThatLength = SomeMatrix.Length;

		const int ThatHeight = SomeMatrix.Height;

		const int ThisLength = Length;

		const int ThisHeight = Height;

		OMPMatrix<SomeType> Temporary(ThisHeight, ThatLength);

		SomeType* ThatColumn = new SomeType[ThatHeight];

		omp_set_num_threads(NUMBER_OF_CORES);

		#pragma omp parallel for schedule(static, NUMBER_OF_CORES)

		for (int SecondIndex = 0; SecondIndex < ThatLength; SecondIndex++) {

			for (int AuxiliaryIndex = 0; AuxiliaryIndex < ThisLength; AuxiliaryIndex++) {

				ThatColumn[AuxiliaryIndex] = SomeMatrix.InnerMatrix[AuxiliaryIndex][SecondIndex];
			}

			for (int FirstIndex = 0; FirstIndex < ThisHeight; FirstIndex++) {

				SomeType* ThisRow = InnerMatrix[FirstIndex];

				SomeType Summation = 0;

				for (int AuxiliaryIndex = 0; AuxiliaryIndex < ThisLength; AuxiliaryIndex++) {

					Summation += ThisRow[AuxiliaryIndex] * ThatColumn[AuxiliaryIndex];
				}

				Temporary.InnerMatrix[FirstIndex][SecondIndex] = Summation;
			}
		}

		delete[] ThatColumn;

		return Temporary;
	};
	OMPMatrix Product(const OMPMatrix& SomeMatrix) const {

		CheckSummationAndSubtraction(SomeMatrix);

		const int Length = this->Length;

		const int Height = this->Height;

		OMPMatrix<SomeType> Temporary(Height, Length);

		omp_set_num_threads(NUMBER_OF_CORES);

		#pragma omp parallel for schedule(static, NUMBER_OF_CORES)

		for (int FirstIndex = 0; FirstIndex < Height; FirstIndex++) {

			for (int SecondIndex = 0; SecondIndex < Length; SecondIndex++) {

				Temporary.InnerMatrix[FirstIndex][SecondIndex] = InnerMatrix[FirstIndex][SecondIndex] * SomeMatrix.InnerMatrix[FirstIndex][SecondIndex];
			}
		}

		return Temporary;
	};
	OMPMatrix Divide(const OMPMatrix& SomeMatrix) const {

		return this->Multiply(SomeMatrix.Inverse());
	};
	SomeType Scalar(const OMPMatrix& SomeMatrix) const {

		CheckSummationAndSubtraction(SomeMatrix);

		const int Length = this->Length;

		const int Height = this->Height;

		SomeType Scalar = 0;

		omp_set_num_threads(NUMBER_OF_CORES);

		#pragma omp parallel for schedule(static, NUMBER_OF_CORES)

		for (int FirstIndex = 0; FirstIndex < Height; FirstIndex++) {

			for (int SecondIndex = 0; SecondIndex < Length; SecondIndex++) {

				Scalar += InnerMatrix[FirstIndex][SecondIndex] * SomeMatrix.InnerMatrix[FirstIndex][SecondIndex];
			}
		}

		return Scalar;
	};
	OMPMatrix Add(const OMPMatrix& SomeMatrix) const {

		CheckSummationAndSubtraction(SomeMatrix);

		const int Length = SomeMatrix.Length;

		const int Height = SomeMatrix.Height;

		OMPMatrix<SomeType> Temporary(Height, Length);

		omp_set_num_threads(NUMBER_OF_CORES);

		#pragma omp parallel for schedule(static, NUMBER_OF_CORES)

		for (int FirstIndex = 0; FirstIndex < Height; FirstIndex++) {

			for (int SecondIndex = 0; SecondIndex < Length; SecondIndex++) {

				Temporary.InnerMatrix[FirstIndex][SecondIndex] = InnerMatrix[FirstIndex][SecondIndex] + SomeMatrix.InnerMatrix[FirstIndex][SecondIndex];
			}
		}

		return Temporary;
	};
	OMPMatrix Multiply(SomeType Coefficient) const {

		CheckParameters();

		const int Length = this->Length;

		const int Height = this->Height;

		OMPMatrix<SomeType> Temporary(Height, Length);

		omp_set_num_threads(NUMBER_OF_CORES);

		#pragma omp parallel for schedule(static, NUMBER_OF_CORES)

		for (int FirstIndex = 0; FirstIndex < Height; FirstIndex++) {

			for (int SecondIndex = 0; SecondIndex < Length; SecondIndex++) {

				Temporary.InnerMatrix[FirstIndex][SecondIndex] = InnerMatrix[FirstIndex][SecondIndex] * Coefficient;
			}
		}

		return Temporary;
	};
	SomeType Determinant() const {

		CheckSquareness();

		const int Length = this->Length;

		const int Height = this->Height;

		OMPMatrix<SomeType> Temporary(*this);

		double Value = 1;

		omp_set_num_threads(NUMBER_OF_CORES);

		#pragma omp parallel for schedule(static, NUMBER_OF_CORES)

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

		return Value;
	};
	OMPMatrix Transpose() const {

		CheckParameters();

		const int Length = this->Length;

		const int Height = this->Height;

		OMPMatrix<SomeType> Temporary(Length, Height);

		omp_set_num_threads(NUMBER_OF_CORES);

#pragma omp parallel for schedule(static, NUMBER_OF_CORES)

		for (int FirstIndex = 0; FirstIndex < Length; FirstIndex++) {

			for (int SecondIndex = 0; SecondIndex < Height; SecondIndex++) {

				Temporary.InnerMatrix[FirstIndex][SecondIndex] = InnerMatrix[SecondIndex][FirstIndex];
			}
		}

		return Temporary;
	};
	OMPMatrix Inverse() const {

		CheckSquareness();

		const int Length = this->Length;

		const int Height = this->Height;

		OMPMatrix<SomeType> Solution(Height, Length);

		OMPMatrix<SomeType> Temporary(*this);

		Solution.UnitFilling();

		omp_set_num_threads(NUMBER_OF_CORES);

#pragma omp parallel for schedule(static, NUMBER_OF_CORES)

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

#pragma omp parallel for schedule(static, NUMBER_OF_CORES)

		for (int FirstIndex = 0; FirstIndex < Height; FirstIndex++) {

			double Coefficient = CheckDivisionByZero(1 / Temporary.InnerMatrix[FirstIndex][FirstIndex]);

			for (int SecondIndex = 0; SecondIndex < Length; SecondIndex++) {

				if (Coefficient != 1) {

					Temporary.InnerMatrix[FirstIndex][SecondIndex] *= Coefficient;

					Solution.InnerMatrix[FirstIndex][SecondIndex] *= Coefficient;
				}
			}
		}

#pragma omp parallel for schedule(static, NUMBER_OF_CORES)

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

		return Solution;
	};

	OMPMatrix& operator *= (const OMPMatrix& SomeMatrix) {

		return *this = this->Multiply(SomeMatrix);
	};
	OMPMatrix& operator += (const OMPMatrix& SomeMatrix) {

		return *this = this->Add(SomeMatrix);
	};
	OMPMatrix& operator /= (const OMPMatrix& SomeMatrix) {

		return *this = this->Divide(SomeMatrix);
	};
	OMPMatrix& operator -= (const OMPMatrix& SomeMatrix) {

		return *this = this->Subtract(SomeMatrix);
	};
	OMPMatrix& operator = (const OMPMatrix& SomeMatrix) {

		if (this == &SomeMatrix) {

			return *this;
		}

		NUMBER_OF_CORES = SomeMatrix.NUMBER_OF_CORES;

		omp_set_num_threads(NUMBER_OF_CORES);

		if (Length != SomeMatrix.Length || Height != SomeMatrix.Height) {

			this->~OMPMatrix();

			Length = SomeMatrix.Length;

			Height = SomeMatrix.Height;

			InnerMatrix = new SomeType*[Height];

			#pragma omp parallel for schedule(static, NUMBER_OF_CORES)

			for (int FirstIndex = 0; FirstIndex < Height; FirstIndex++) {

				InnerMatrix[FirstIndex] = new SomeType[Length];

				for (int SecondIndex = 0; SecondIndex < Length; SecondIndex++) {

					InnerMatrix[FirstIndex][SecondIndex] = SomeMatrix.InnerMatrix[FirstIndex][SecondIndex];
				}
			}
		}

		#pragma omp parallel for schedule(static, NUMBER_OF_CORES)

		for (int FirstIndex = 0; FirstIndex < Height; FirstIndex++) {

			for (int SecondIndex = 0; SecondIndex < Length; SecondIndex++) {

				InnerMatrix[FirstIndex][SecondIndex] = SomeMatrix.InnerMatrix[FirstIndex][SecondIndex];
			}
		}

		return *this;
	};
	OMPMatrix& operator *= (const SomeType Coefficient) {

		return *this = this->Multiply(Coefficient);
	};
	OMPMatrix& operator ^= (const int Power) {

		return *this = this->Exponentiate(Power);
	};

	tuple<OMPMatrix, OMPMatrix> LUPDecomposition() const {

		CheckSquareness();

		const int Length = this->Length;

		const int Height = this->Height;

		OMPMatrix<SomeType> Diagonal(Height, Length);

		OMPMatrix<SomeType> Temporary(*this);

		Diagonal.UnitFilling();

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

		return tuple<OMPMatrix, OMPMatrix>(Temporary, Diagonal);
	};
	tuple<OMPMatrix, OMPMatrix> LUDecomposition() const {

		CheckSquareness();

		const int Length = this->Length;

		const int Height = this->Height;

		OMPMatrix<SomeType> Modified(Height, Length);

		OMPMatrix<SomeType> Temporary(*this);

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

		return tuple<OMPMatrix, OMPMatrix>(Temporary, Modified);
	};
	OMPMatrix CholeskyDecomposition() const {

		CheckSquareness();

		const int Length = this->Length;

		const int Height = this->Height;

		OMPMatrix<SomeType> Temporary(Height, Length);

		for (int FirstIndex = 0; FirstIndex < Height; FirstIndex++) {

			double Coefficient;

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

		return Temporary;
	};

	OMPMatrix GetColumn(int Position) const {

		CheckParameters();

		const int Height = this->Height;

		OMPMatrix<SomeType> Temporary(Height);

		for (int Index = 0; Index < Height; Index++) {

			Temporary.InnerMatrix[Index][0] = InnerMatrix[Index][Position];
		}

		return Temporary;
	};
	OMPMatrix GetRow(int Position) const {

		CheckParameters();

		const int Length = this->Length;

		OMPMatrix<SomeType> Temporary(1, Length);

		for (int Index = 0; Index < Length; Index++) {

			Temporary.InnerMatrix[0][Index] = InnerMatrix[Position][Index];
		}

		return Temporary;
	};
	OMPMatrix GetMainDiagonal() const {

		CheckParameters();

		const int Height = this->Length;

		OMPMatrix<SomeType> Temporary(Height);

		for (int Index = 0; Index < Height; Index++) {

			Temporary.InnerMatrix[Index][0] = InnerMatrix[Index][Index];
		}

		return Temporary;
	};
	OMPMatrix GetSideDiagonal() const {

		CheckParameters();

		const int Height = this->Height;

		OMPMatrix<SomeType> Temporary(Height);

		for (int Index = 0; Index < Height; Index++) {

			Temporary.InnerMatrix[Index][0] = InnerMatrix[Index][Length - Index - 1];
		}

		return Temporary;
	};
	int GetLength() const {

		return Length;
	};
	int GetHeight() const {

		return Height;
	};
	SomeType* operator[] (int Index) const {

		return InnerMatrix[Index];
	};
	
	vector<vector<SomeType>> ToVector() {

		CheckParameters();

		const int Length = this->Length;

		const int Height = this->Height;

		vector<vector<SomeType>> Temporary;

		for (int FirstIndex = 0; FirstIndex < Height; FirstIndex++) {

			Temporary.push_back(vector<SomeType>(Length));

			for (int SecondIndex = 0; SecondIndex < Length; SecondIndex++) {

				Temporary[FirstIndex][SecondIndex] = InnerMatrix[FirstIndex][SecondIndex];
			}
		}

		return Temporary;
	};
	SomeType** ToPointer() {

		CheckParameters();

		const int Length = this->Length;

		const int Height = this->Height;

		SomeType** Temporary = new SomeType*[Height];

		for (int FirstIndex = 0; FirstIndex < Height; FirstIndex++) {

			Temporary[FirstIndex] = new SomeType[Length];

			for (int SecondIndex = 0; SecondIndex < Length; SecondIndex++) {

				Temporary[FirstIndex][SecondIndex] = InnerMatrix[FirstIndex][SecondIndex];
			}
		}

		return Temporary;
	};

	SomeType ColumnSum(int Position) {

		SomeType Value = 0;

		for (int Index = 0; Index < Height; Index++) {

			Value += InnerMatrix[Index][Position];
		}

		return Value;
	};
	SomeType RowSum(int Position) {

		SomeType Value = 0;

		for (int Index = 0; Index < Length; Index++) {

			Value += InnerMatrix[Position][Index];
		}

		return Value;
	};
	SomeType Average() {

		CheckParameters();

		const int Length = this->Length;

		const int Height = this->Height;

		SomeType Value = 0;

		for (int FirstIndex = 0; FirstIndex < Height; FirstIndex++) {

			for (int SecondIndex = 0; SecondIndex < Length; SecondIndex++) {

				Value += InnerMatrix[FirstIndex][SecondIndex];
			}
		}

		return Value / (Height * Length);
	}
	SomeType Rate() {

		CheckParameters();

		const int Length = this->Length;

		const int Height = this->Height;

		SomeType Value = 0;

		for (int FirstIndex = 0; FirstIndex < Height; FirstIndex++) {

			for (int SecondIndex = 0; SecondIndex < Length; SecondIndex++) {

				Value += InnerMatrix[FirstIndex][SecondIndex] * InnerMatrix[FirstIndex][SecondIndex];
			}
		}

		return sqrt(Value);
	};
	SomeType Max() {

		CheckParameters();

		const int Length = this->Length;

		const int Height = this->Height;

		SomeType Value = 0;

		for (int FirstIndex = 0; FirstIndex < Height; FirstIndex++) {

			for (int SecondIndex = 0; SecondIndex < Length; SecondIndex++) {

				Value = InnerMatrix[FirstIndex][SecondIndex] > Value ? InnerMatrix[FirstIndex][SecondIndex] : Value;
			}
		}

		return Value;
	};
#pragma endregion
};
