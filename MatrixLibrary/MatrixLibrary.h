#include <iostream>
#include <complex>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <chrono>
#include <vector>
#include <string>
#include <omp.h>
#include <cmath>
#include <ctime>
#include <tuple>

using namespace std;

template<typename Type>
class OmpMatrix {
#pragma region
protected:
	inline OmpMatrix MultiplicationForLengthDividedBy_5(const OmpMatrix& SomeMatrix) const {

		const int ThatLength = SomeMatrix.Length;
		const int ThatHeight = SomeMatrix.Height;

		const int Length = this->Length;
		const int Height = this->Height;

		omp_set_num_threads(NUMBER_OF_CORES);

		OmpMatrix<Type> Temporary(Height, ThatLength);

		#pragma omp parallel for schedule(static, NUMBER_OF_CORES)

		for (int SecondIndex = 0; SecondIndex < ThatLength; SecondIndex += 5) {

			Type* ThatColumn_4 = new Type[ThatHeight];
			Type* ThatColumn_3 = new Type[ThatHeight];
			Type* ThatColumn_2 = new Type[ThatHeight];
			Type* ThatColumn_1 = new Type[ThatHeight];
			Type* ThatColumn_0 = new Type[ThatHeight];

			for (int AuxiliaryIndex = 0; AuxiliaryIndex < Length; AuxiliaryIndex++) {

				ThatColumn_4[AuxiliaryIndex] = SomeMatrix.InnerMatrix[AuxiliaryIndex][4 + SecondIndex];
				ThatColumn_3[AuxiliaryIndex] = SomeMatrix.InnerMatrix[AuxiliaryIndex][3 + SecondIndex];
				ThatColumn_2[AuxiliaryIndex] = SomeMatrix.InnerMatrix[AuxiliaryIndex][2 + SecondIndex];
				ThatColumn_1[AuxiliaryIndex] = SomeMatrix.InnerMatrix[AuxiliaryIndex][1 + SecondIndex];
				ThatColumn_0[AuxiliaryIndex] = SomeMatrix.InnerMatrix[AuxiliaryIndex][0 + SecondIndex];
			}

			for (int FirstIndex = 0; FirstIndex < Height; FirstIndex++) {

				Type* ThisRow = InnerMatrix[FirstIndex];

				Type Summation_4 = 0;
				Type Summation_3 = 0;
				Type Summation_2 = 0;
				Type Summation_1 = 0;
				Type Summation_0 = 0;

				for (int AuxiliaryIndex = 0; AuxiliaryIndex < Length; AuxiliaryIndex++) {

					Summation_4 += ThisRow[AuxiliaryIndex] * ThatColumn_4[AuxiliaryIndex];
					Summation_3 += ThisRow[AuxiliaryIndex] * ThatColumn_3[AuxiliaryIndex];
					Summation_2 += ThisRow[AuxiliaryIndex] * ThatColumn_2[AuxiliaryIndex];
					Summation_1 += ThisRow[AuxiliaryIndex] * ThatColumn_1[AuxiliaryIndex];
					Summation_0 += ThisRow[AuxiliaryIndex] * ThatColumn_0[AuxiliaryIndex];
				}

				Temporary.InnerMatrix[FirstIndex][4 + SecondIndex] = Summation_4;
				Temporary.InnerMatrix[FirstIndex][3 + SecondIndex] = Summation_3;
				Temporary.InnerMatrix[FirstIndex][2 + SecondIndex] = Summation_2;
				Temporary.InnerMatrix[FirstIndex][1 + SecondIndex] = Summation_1;
				Temporary.InnerMatrix[FirstIndex][0 + SecondIndex] = Summation_0;
			}

			delete[] ThatColumn_4;
			delete[] ThatColumn_3;
			delete[] ThatColumn_2;
			delete[] ThatColumn_1;
			delete[] ThatColumn_0;
		}

		return Temporary;
	};
	inline OmpMatrix MultiplicationForLengthDividedBy_4(const OmpMatrix& SomeMatrix) const {

		const int ThatLength = SomeMatrix.Length;
		const int ThatHeight = SomeMatrix.Height;

		const int Length = this->Length;
		const int Height = this->Height;

		omp_set_num_threads(NUMBER_OF_CORES);

		OmpMatrix<Type> Temporary(Height, ThatLength);

		#pragma omp parallel for schedule(static, NUMBER_OF_CORES)

		for (int SecondIndex = 0; SecondIndex < ThatLength; SecondIndex += 4) {

			Type* ThatColumn_3 = new Type[ThatHeight];
			Type* ThatColumn_2 = new Type[ThatHeight];
			Type* ThatColumn_1 = new Type[ThatHeight];
			Type* ThatColumn_0 = new Type[ThatHeight];

			for (int AuxiliaryIndex = 0; AuxiliaryIndex < Length; AuxiliaryIndex++) {

				ThatColumn_3[AuxiliaryIndex] = SomeMatrix.InnerMatrix[AuxiliaryIndex][3 + SecondIndex];
				ThatColumn_2[AuxiliaryIndex] = SomeMatrix.InnerMatrix[AuxiliaryIndex][2 + SecondIndex];
				ThatColumn_1[AuxiliaryIndex] = SomeMatrix.InnerMatrix[AuxiliaryIndex][1 + SecondIndex];
				ThatColumn_0[AuxiliaryIndex] = SomeMatrix.InnerMatrix[AuxiliaryIndex][0 + SecondIndex];
			}

			for (int FirstIndex = 0; FirstIndex < Height; FirstIndex++) {

				Type* ThisRow = InnerMatrix[FirstIndex];

				Type Summation_3 = 0;
				Type Summation_2 = 0;
				Type Summation_1 = 0;
				Type Summation_0 = 0;

				for (int AuxiliaryIndex = 0; AuxiliaryIndex < Length; AuxiliaryIndex++) {

					Summation_3 += ThisRow[AuxiliaryIndex] * ThatColumn_3[AuxiliaryIndex];
					Summation_2 += ThisRow[AuxiliaryIndex] * ThatColumn_2[AuxiliaryIndex];
					Summation_1 += ThisRow[AuxiliaryIndex] * ThatColumn_1[AuxiliaryIndex];
					Summation_0 += ThisRow[AuxiliaryIndex] * ThatColumn_0[AuxiliaryIndex];
				}

				Temporary.InnerMatrix[FirstIndex][3 + SecondIndex] = Summation_3;
				Temporary.InnerMatrix[FirstIndex][2 + SecondIndex] = Summation_2;
				Temporary.InnerMatrix[FirstIndex][1 + SecondIndex] = Summation_1;
				Temporary.InnerMatrix[FirstIndex][0 + SecondIndex] = Summation_0;
			}

			delete[] ThatColumn_3;
			delete[] ThatColumn_2;
			delete[] ThatColumn_1;
			delete[] ThatColumn_0;
		}

		return Temporary;
	};
	inline OmpMatrix MultiplicationForLengthDividedBy_3(const OmpMatrix& SomeMatrix) const {

		const int ThatLength = SomeMatrix.Length;
		const int ThatHeight = SomeMatrix.Height;

		const int Length = this->Length;
		const int Height = this->Height;

		omp_set_num_threads(NUMBER_OF_CORES);

		OmpMatrix<Type> Temporary(Height, ThatLength);

		#pragma omp parallel for schedule(static, NUMBER_OF_CORES)

		for (int SecondIndex = 0; SecondIndex < ThatLength; SecondIndex += 3) {

			Type* ThatColumn_2 = new Type[ThatHeight];
			Type* ThatColumn_1 = new Type[ThatHeight];
			Type* ThatColumn_0 = new Type[ThatHeight];

			for (int AuxiliaryIndex = 0; AuxiliaryIndex < Length; AuxiliaryIndex++) {

				ThatColumn_2[AuxiliaryIndex] = SomeMatrix.InnerMatrix[AuxiliaryIndex][2 + SecondIndex];
				ThatColumn_1[AuxiliaryIndex] = SomeMatrix.InnerMatrix[AuxiliaryIndex][1 + SecondIndex];
				ThatColumn_0[AuxiliaryIndex] = SomeMatrix.InnerMatrix[AuxiliaryIndex][0 + SecondIndex];
			}

			for (int FirstIndex = 0; FirstIndex < Height; FirstIndex++) {

				Type* ThisRow = InnerMatrix[FirstIndex];

				Type Summation_2 = 0;
				Type Summation_1 = 0;
				Type Summation_0 = 0;

				for (int AuxiliaryIndex = 0; AuxiliaryIndex < Length; AuxiliaryIndex++) {

					Summation_2 += ThisRow[AuxiliaryIndex] * ThatColumn_2[AuxiliaryIndex];
					Summation_1 += ThisRow[AuxiliaryIndex] * ThatColumn_1[AuxiliaryIndex];
					Summation_0 += ThisRow[AuxiliaryIndex] * ThatColumn_0[AuxiliaryIndex];
				}

				Temporary.InnerMatrix[FirstIndex][2 + SecondIndex] = Summation_2;
				Temporary.InnerMatrix[FirstIndex][1 + SecondIndex] = Summation_1;
				Temporary.InnerMatrix[FirstIndex][0 + SecondIndex] = Summation_0;
			}

			delete[] ThatColumn_2;
			delete[] ThatColumn_1;
			delete[] ThatColumn_0;
		}

		return Temporary;
	};
	inline OmpMatrix MultiplicationForLengthDividedBy_2(const OmpMatrix& SomeMatrix) const {

		const int ThatLength = SomeMatrix.Length;
		const int ThatHeight = SomeMatrix.Height;

		const int Length = this->Length;
		const int Height = this->Height;

		omp_set_num_threads(NUMBER_OF_CORES);

		OmpMatrix<Type> Temporary(Height, ThatLength);

		#pragma omp parallel for schedule(static, NUMBER_OF_CORES)

		for (int SecondIndex = 0; SecondIndex < ThatLength; SecondIndex += 2) {

			Type* ThatColumn_1 = new Type[ThatHeight];
			Type* ThatColumn_0 = new Type[ThatHeight];

			for (int AuxiliaryIndex = 0; AuxiliaryIndex < Length; AuxiliaryIndex++) {

				ThatColumn_1[AuxiliaryIndex] = SomeMatrix.InnerMatrix[AuxiliaryIndex][1 + SecondIndex];
				ThatColumn_0[AuxiliaryIndex] = SomeMatrix.InnerMatrix[AuxiliaryIndex][0 + SecondIndex];
			}

			for (int FirstIndex = 0; FirstIndex < Height; FirstIndex++) {

				Type* ThisRow = InnerMatrix[FirstIndex];

				Type Summation_1 = 0;
				Type Summation_0 = 0;

				for (int AuxiliaryIndex = 0; AuxiliaryIndex < Length; AuxiliaryIndex++) {

					Summation_1 += ThisRow[AuxiliaryIndex] * ThatColumn_1[AuxiliaryIndex];
					Summation_0 += ThisRow[AuxiliaryIndex] * ThatColumn_0[AuxiliaryIndex];
				}

				Temporary.InnerMatrix[FirstIndex][1 + SecondIndex] = Summation_1;
				Temporary.InnerMatrix[FirstIndex][0 + SecondIndex] = Summation_0;
			}

			delete[] ThatColumn_1;
			delete[] ThatColumn_0;
		}

		return Temporary;
	};
	inline OmpMatrix MultiplicationForLengthDividedBy_1(const OmpMatrix& SomeMatrix) const {
		
		const int ThatLength = SomeMatrix.Length;
		const int ThatHeight = SomeMatrix.Height;

		const int Length = this->Length;
		const int Height = this->Height;

		omp_set_num_threads(NUMBER_OF_CORES);

		OmpMatrix<Type> Temporary(Height, ThatLength);

		#pragma omp parallel for schedule(static, NUMBER_OF_CORES)

		for (int SecondIndex = 0; SecondIndex < ThatLength; SecondIndex++) {

			Type* ThatColumn = new Type[ThatHeight];

			for (int AuxiliaryIndex = 0; AuxiliaryIndex < Length; AuxiliaryIndex++) {

				ThatColumn[AuxiliaryIndex] = SomeMatrix.InnerMatrix[AuxiliaryIndex][SecondIndex];
			}

			for (int FirstIndex = 0; FirstIndex < Height; FirstIndex++) {

				Type* ThisRow = InnerMatrix[FirstIndex];

				Type Summation = 0;

				for (int AuxiliaryIndex = 0; AuxiliaryIndex < Length; AuxiliaryIndex++) {

					Summation += ThisRow[AuxiliaryIndex] * ThatColumn[AuxiliaryIndex];
				}

				Temporary.InnerMatrix[FirstIndex][SecondIndex] = Summation;
			}

			delete[] ThatColumn;
		}

		return Temporary;
	};

	inline void CheckSummationAndSubtraction(const OmpMatrix& SomeMatrix) const {

		if (Length != SomeMatrix.Length || Height != SomeMatrix.Height) {

			throw "The length/height of the first matrix doesn't equal the length/height of the second matrix";
		}
	};
	inline bool CheckPositiveDefinition(const OmpMatrix& SomeMatrix) const {

		OmpMatrix<Type> Temporary(SomeMatrix.Height);

		Temporary.Randomize();

		if ((Temporary.Transpose().Multiply(SomeMatrix).Multiply(Temporary))[0][0] > 0) {

			return true;
		}

		return false;
	};
	inline void CheckMultiplication(const OmpMatrix& SomeMatrix) const {

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

		Type* Temporary = InnerMatrix[FirstIndex];

		InnerMatrix[FirstIndex] = InnerMatrix[SecondIndex];

		InnerMatrix[SecondIndex] = Temporary;
	};

	const int MAX_NUMBER = omp_get_max_threads();
	int NUMBER_OF_CORES;
	Type** InnerMatrix;
	int Length;
	int Height;
#pragma endregion

#pragma region
public:
    OmpMatrix(Type** SomeMatrix, long NUMBER_OF_CORES) : Length(_msize(SomeMatrix[0]) / sizeof(SomeMatrix[0][0])), Height(_msize(SomeMatrix) / sizeof(SomeMatrix[0])), NUMBER_OF_CORES(NUMBER_OF_CORES) {

		const int ThisLength = Length;
		const int ThisHeight = Height;

		InnerMatrix = new Type*[ThisHeight];

		omp_set_num_threads(NUMBER_OF_CORES);

		#pragma omp parallel for schedule(static, NUMBER_OF_CORES)

		for (int FirstIndex = 0; FirstIndex < ThisHeight; FirstIndex++) {

			InnerMatrix[FirstIndex] = new Type[ThisLength];

			for (int SecondIndex = 0; SecondIndex < ThisLength; SecondIndex++) {

				InnerMatrix[FirstIndex][SecondIndex] = SomeMatrix[FirstIndex][SecondIndex];
			}
		}
	};
    OmpMatrix(Type** SomeMatrix) : Length(_msize(SomeMatrix[0]) / sizeof(SomeMatrix[0][0])), Height(_msize(SomeMatrix) / sizeof(SomeMatrix[0])), NUMBER_OF_CORES(MAX_NUMBER) {

		const int ThisLength = Length;
		const int ThisHeight = Height;

		InnerMatrix = new Type*[ThisHeight];

		omp_set_num_threads(NUMBER_OF_CORES);

		#pragma omp parallel for schedule(static, NUMBER_OF_CORES)

		for (int FirstIndex = 0; FirstIndex < ThisHeight; FirstIndex++) {

			InnerMatrix[FirstIndex] = new Type[ThisLength];

			for (int SecondIndex = 0; SecondIndex < ThisLength; SecondIndex++) {

				InnerMatrix[FirstIndex][SecondIndex] = SomeMatrix[FirstIndex][SecondIndex];
			}
		}
	};
    OmpMatrix(vector<vector<Type>> SomeMatrix, long NUMBER_OF_CORES) : Length(size(SomeMatrix[0])), Height(size(SomeMatrix)), NUMBER_OF_CORES(NUMBER_OF_CORES) {

		const int ThisLength = Length;
		const int ThisHeight = Height;

		InnerMatrix = new Type*[ThisHeight];

		omp_set_num_threads(NUMBER_OF_CORES);

		#pragma omp parallel for schedule(static, NUMBER_OF_CORES)

		for (int FirstIndex = 0; FirstIndex < ThisHeight; FirstIndex++) {

			InnerMatrix[FirstIndex] = new Type[ThisLength];

			for (int SecondIndex = 0; SecondIndex < ThisLength; SecondIndex++) {

				InnerMatrix[FirstIndex][SecondIndex] = SomeMatrix[FirstIndex][SecondIndex];
			}
		}
	};
    OmpMatrix(const OmpMatrix& SomeMatrix) : Length(SomeMatrix.Length), Height(SomeMatrix.Height), NUMBER_OF_CORES(SomeMatrix.NUMBER_OF_CORES) {

		const int ThisLength = Length;
		const int ThisHeight = Height;

		InnerMatrix = new Type*[ThisHeight];

		omp_set_num_threads(NUMBER_OF_CORES);

		#pragma omp parallel for schedule(static, NUMBER_OF_CORES)

		for (int FirstIndex = 0; FirstIndex < ThisHeight; FirstIndex++) {

			InnerMatrix[FirstIndex] = new Type[ThisLength];

			for (int SecondIndex = 0; SecondIndex < ThisLength; SecondIndex++) {

				InnerMatrix[FirstIndex][SecondIndex] = SomeMatrix.InnerMatrix[FirstIndex][SecondIndex];
			}
		}
	};
	OmpMatrix(vector<vector<Type>> SomeMatrix) : Length(size(SomeMatrix[0])), Height(size(SomeMatrix)), NUMBER_OF_CORES(MAX_NUMBER) {

		const int ThisLength = Length;
		const int ThisHeight = Height;

		InnerMatrix = new Type*[ThisHeight];

		omp_set_num_threads(NUMBER_OF_CORES);

		#pragma omp parallel for schedule(static, NUMBER_OF_CORES)

		for (int FirstIndex = 0; FirstIndex < ThisHeight; FirstIndex++) {

			InnerMatrix[FirstIndex] = new Type[ThisLength];

			for (int SecondIndex = 0; SecondIndex < ThisLength; SecondIndex++) {

				InnerMatrix[FirstIndex][SecondIndex] = SomeMatrix[FirstIndex][SecondIndex];
			}
		}
	};
	OmpMatrix(int Height, int Length, long NUMBER_OF_CORES) : Length(Length), Height(Height), NUMBER_OF_CORES(NUMBER_OF_CORES) {

		const int ThisLength = Length;
		const int ThisHeight = Height;

		InnerMatrix = new Type*[ThisHeight];

		omp_set_num_threads(NUMBER_OF_CORES);

		#pragma omp parallel for schedule(static, NUMBER_OF_CORES)

		for (int FirstIndex = 0; FirstIndex < ThisHeight; FirstIndex++) {

			InnerMatrix[FirstIndex] = new Type[ThisLength];

			for (int SecondIndex = 0; SecondIndex < ThisLength; SecondIndex++) {

				InnerMatrix[FirstIndex][SecondIndex] = 0;
			}
		}
	};
	OmpMatrix(int Height, long NUMBER_OF_CORES) : Length(1), Height(Height), NUMBER_OF_CORES(NUMBER_OF_CORES) {

		const int ThisHeight = Height;

		InnerMatrix = new Type*[ThisHeight];

		omp_set_num_threads(NUMBER_OF_CORES);

		#pragma omp parallel for schedule(static, NUMBER_OF_CORES)

		for (int Index = 0; Index < ThisHeight; Index++) {

			InnerMatrix[Index] = new Type[1];

			InnerMatrix[Index][0] = 0;
		}
	};
	OmpMatrix(int Height, int Length) : Length(Length), Height(Height), NUMBER_OF_CORES(MAX_NUMBER) {

		const int ThisLength = Length;
		const int ThisHeight = Height;

		InnerMatrix = new Type*[ThisHeight];

		omp_set_num_threads(NUMBER_OF_CORES);

		#pragma omp parallel for schedule(static, NUMBER_OF_CORES)

		for (int FirstIndex = 0; FirstIndex < ThisHeight; FirstIndex++) {

			InnerMatrix[FirstIndex] = new Type[ThisLength];

			for (int SecondIndex = 0; SecondIndex < ThisLength; SecondIndex++) {

				InnerMatrix[FirstIndex][SecondIndex] = 0;
			}
		}
	};
	OmpMatrix(int Height) : Length(1), Height(Height), NUMBER_OF_CORES(MAX_NUMBER) {

		const int ThisHeight = Height;

		InnerMatrix = new Type*[ThisHeight];

		omp_set_num_threads(NUMBER_OF_CORES);

		#pragma omp parallel for schedule(static, NUMBER_OF_CORES)

		for (int Index = 0; Index < ThisHeight; Index++) {

			InnerMatrix[Index] = new Type[1];

			InnerMatrix[Index][0] = 0;
		}
	};
	~OmpMatrix() {

		const int ThisHeight = Height;
		
		omp_set_num_threads(NUMBER_OF_CORES);

		#pragma omp parallel for schedule(static, NUMBER_OF_CORES)

		for (int Index = 0; Index < ThisHeight; Index++) {

			delete[] InnerMatrix[Index];
		}

		delete[] InnerMatrix;
	};
	OmpMatrix() {

		InnerMatrix = nullptr;

		NUMBER_OF_CORES = 0;

		Length = 0;
		Height = 0;
	};

	friend ostream& operator << (ostream& Stream, const OmpMatrix<Type>& SomeMatrix) {

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

	friend istream& operator >> (istream& Stream, const OmpMatrix<Type>& SomeMatrix) {

		const int Length = SomeMatrix.Length;

		const int Height = SomeMatrix.Height;

		for (int FirstIndex = 0; FirstIndex < Height; FirstIndex++) {

			for (int SecondIndex = 0; SecondIndex < Length; SecondIndex++) {

				Stream >> SomeMatrix.InnerMatrix[FirstIndex][SecondIndex];
			}
		}

		return Stream;
	};


	OmpMatrix Replace(int Position, const OmpMatrix& SomeVector) const {

		CheckParameters();

		const int Length = this->Length;
		const int Height = this->Height;

		OmpMatrix<Type> Temporary(*this);

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
	OmpMatrix Cut(int StartX, int StartY, int EndX, int EndY) const {

		CheckParameters();

		const int Height = (EndX - StartX) + 1;
		const int Length = (EndY - StartY) + 1;

		OmpMatrix<Type> Temporary(Height, Length);

		for (int FirstIndex = StartY; FirstIndex < Height; FirstIndex++) {

			for (int SecondIndex = StartX; SecondIndex < Length; SecondIndex++) {

				Temporary.InnerMatrix[FirstIndex - StartY][SecondIndex - StartX] = InnerMatrix[FirstIndex][SecondIndex];
			}
		}

		return Temporary;
	};
	OmpMatrix LoadFromFIle(string Path, char Delimiter) const {

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

		OmpMatrix<Type> Temporary(Height, Length);

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
	OmpMatrix FillBand(int LowerSize, int UpperSize) const {

		CheckSquareness();

		const int Length = this->Length;
		const int Height = this->Height;

		for (int FirstIndex = 0; FirstIndex < Height; FirstIndex++) {

			for (int SeconIndex = 0; SeconIndex <= UpperSize; SeconIndex++) {

				if (FirstIndex + SeconIndex < Length) {

					InnerMatrix[FirstIndex][FirstIndex + SeconIndex] = static_cast <float> (rand()) / static_cast <float> (RAND_MAX) - 0.5f;
				}
			}

			for (int SeconIndex = 0; SeconIndex <= LowerSize; SeconIndex++) {

				if (FirstIndex - SeconIndex >= 0) {

					InnerMatrix[FirstIndex][FirstIndex - SeconIndex] = static_cast <float> (rand()) / static_cast <float> (RAND_MAX) - 0.5f;
				}
			}
		}

		return *this;
	};
	OmpMatrix FillColumn(int Position, Type Value) const {

		CheckParameters();

		const int Height = this->Height;

		for (int Index = 0; Index < Height; Index++) {

			InnerMatrix[Index][Position] = Value;
		}

		return *this;
	};
    OmpMatrix UseFunction(Type(*Function)(Type)) const {

		CheckParameters();

		const int Length = this->Length;
		const int Height = this->Height;

		OmpMatrix<Type> Temporary(Height, Length);

		for (int FirstIndex = 0; FirstIndex < Height; FirstIndex++) {

			for (int SecondIndex = 0; SecondIndex < Length; SecondIndex++) {

				Temporary.InnerMatrix[FirstIndex][SecondIndex] = Function(InnerMatrix[FirstIndex][SecondIndex]);
			}
		}

		return Temporary;
	};
	OmpMatrix FillRow(int Position, Type Value) const {

		CheckParameters();

		const int Height = this->Height;

		for (int Index = 0; Index < Height; Index++) {

			InnerMatrix[Position][Index] = Value;
		}

		return *this;
	};
	OmpMatrix Reshape(int Height, int Length) const {

		CheckParameters();

		const int ThisLength = this->Length;
		const int ThisHeight = this->Height;

		Type* DeconvolutedMatrix = new Type[ThisLength * ThisHeight];

		OmpMatrix<Type> Temporary(Height, Length);

		int InnerHeight = 0;

		for (int FirstIndex = 0; FirstIndex < ThisHeight; FirstIndex++) {

			for (int SecondIndex = 0; SecondIndex < ThisLength; SecondIndex++) {

				(DeconvolutedMatrix)[InnerHeight] = InnerMatrix[FirstIndex][SecondIndex];

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

			Temporary.InnerMatrix[InnerHeight][Counter] = (DeconvolutedMatrix)[Index];

			Counter++;
		}

		delete[] DeconvolutedMatrix;

		return Temporary;
	};
	OmpMatrix FillMainDiagonal(Type Value) const {

		CheckSquareness();

		const int Height = this->Height;

		for (int Index = 0; Index < Height; Index++) {

			InnerMatrix[Index][Index] = Value;
		}

		return *this;
	};
	OmpMatrix FillSideDiagonal(Type Value) const {

		CheckSquareness();

		const int Height = this->Height;

		for (int Index = 0; Index < Height; Index++) {

			InnerMatrix[Index][Height - Index - 1] = Value;
		}

		return *this;

	};
	OmpMatrix FillSimmetrically() const {

		CheckParameters();

		const int Length = this->Length;
		const int Height = this->Height;

		srand(1 + rand() % 100);

		for (int FirstIndex = 0; FirstIndex < Height; FirstIndex++) {

			for (int SecondIndex = 0; SecondIndex < Length; SecondIndex++) {

				InnerMatrix[FirstIndex][SecondIndex] = InnerMatrix[SecondIndex][FirstIndex] = static_cast <float> (rand()) / static_cast <float> (RAND_MAX) - 0.5f;
			}
		}

		return *this;
	};
	OmpMatrix Fill(Type Value) const {

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
	OmpMatrix Randomize() const {

		CheckParameters();

		const int Length = this->Length;
		const int Height = this->Height;

		srand(1 + rand() % Height * 2);

		for (int FirstIndex = 0; FirstIndex < Height; FirstIndex++) {

			for (int SecondIndex = 0; SecondIndex < Length; SecondIndex++) {

				InnerMatrix[FirstIndex][SecondIndex] = static_cast <Type> (rand());
			}
		}

		return *this;
	};
	OmpMatrix Clear() const {

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
	OmpMatrix Print() const {

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

	OmpMatrix operator * (const OmpMatrix& SomeMatrix) const {

		return this->Multiply(SomeMatrix);
	};
	OmpMatrix operator + (const OmpMatrix& SomeMatrix) const {

		return this->Add(SomeMatrix);
	};
	OmpMatrix operator / (const OmpMatrix& SomeMatrix) const {

		return this->Divide(SomeMatrix);
	};
	OmpMatrix operator - (const OmpMatrix& SomeMatrix) const {

		return this->Subtract(SomeMatrix);
	};
	OmpMatrix operator * (const Type Coefficient) const {

		return this->Multiply(Coefficient);
	};
	OmpMatrix operator ^ (const int Power) const {

		return this->Exponentiate(Power);
	};
	OmpMatrix operator ~ () {

		return this->Transpose();
	};	

	OmpMatrix Subtract(const OmpMatrix& SomeMatrix) const {

		CheckSummationAndSubtraction(SomeMatrix);

		const int Length = this->Length;
		const int Height = this->Height;

		omp_set_num_threads(NUMBER_OF_CORES);

		OmpMatrix<Type> Temporary(Height, Length);

		#pragma omp parallel for schedule(static, NUMBER_OF_CORES)

		for (int FirstIndex = 0; FirstIndex < Height; FirstIndex++) {

			for (int SecondIndex = 0; SecondIndex < Length; SecondIndex++) {

				Temporary.InnerMatrix[FirstIndex][SecondIndex] = InnerMatrix[FirstIndex][SecondIndex] - SomeMatrix.InnerMatrix[FirstIndex][SecondIndex];
			}
		}

		return Temporary;
	};
	OmpMatrix Multiply(const OmpMatrix& SomeMatrix) const {

		CheckMultiplication(SomeMatrix);

		const int Length = SomeMatrix.Length;

		if (Length % 5 == 0) {

			return MultiplicationForLengthDividedBy_5(SomeMatrix);
		}
		if (Length % 4 == 0) {

			return MultiplicationForLengthDividedBy_4(SomeMatrix);
		}
		if (Length % 3 == 0) {

			return MultiplicationForLengthDividedBy_3(SomeMatrix);
		}
		if (Length % 2 == 0) {

			return MultiplicationForLengthDividedBy_2(SomeMatrix);
		}

		return MultiplicationForLengthDividedBy_1(SomeMatrix);
	};
	OmpMatrix Product(const OmpMatrix& SomeMatrix) const {

		CheckSummationAndSubtraction(SomeMatrix);

		const int Length = this->Length;
		const int Height = this->Height;

		omp_set_num_threads(NUMBER_OF_CORES);

		OmpMatrix<Type> Temporary(Height, Length);

		#pragma omp parallel for schedule(static, NUMBER_OF_CORES)

		for (int FirstIndex = 0; FirstIndex < Height; FirstIndex++) {

			for (int SecondIndex = 0; SecondIndex < Length; SecondIndex++) {

				Temporary.InnerMatrix[FirstIndex][SecondIndex] = InnerMatrix[FirstIndex][SecondIndex] * SomeMatrix.InnerMatrix[FirstIndex][SecondIndex];
			}
		}

		return Temporary;
	};
	OmpMatrix Divide(const OmpMatrix& SomeMatrix) const {

		return this->Multiply(SomeMatrix.Inverse());
	};
	OmpMatrix Add(const OmpMatrix& SomeMatrix) const {

		CheckSummationAndSubtraction(SomeMatrix);

		const int Length = SomeMatrix.Length;
		const int Height = SomeMatrix.Height;

		omp_set_num_threads(NUMBER_OF_CORES);

		OmpMatrix<Type> Temporary(Height, Length);

		#pragma omp parallel for schedule(static, NUMBER_OF_CORES)

		for (int FirstIndex = 0; FirstIndex < Height; FirstIndex++) {

			for (int SecondIndex = 0; SecondIndex < Length; SecondIndex++) {

				Temporary.InnerMatrix[FirstIndex][SecondIndex] = InnerMatrix[FirstIndex][SecondIndex] + SomeMatrix.InnerMatrix[FirstIndex][SecondIndex];
			}
		}

		return Temporary;
	};
	Type Scalar(const OmpMatrix& SomeMatrix) const {

		CheckSummationAndSubtraction(SomeMatrix);

		Type Value = 0;

		const int Length = this->Length;
		const int Height = this->Height;

		omp_set_num_threads(NUMBER_OF_CORES);

		#pragma omp parallel for schedule(static, NUMBER_OF_CORES)

		for (int FirstIndex = 0; FirstIndex < Height; FirstIndex++) {

			for (int SecondIndex = 0; SecondIndex < Length; SecondIndex++) {

				Value += InnerMatrix[FirstIndex][SecondIndex] * SomeMatrix.InnerMatrix[FirstIndex][SecondIndex];
			}
		}

		return Value;
	};
	OmpMatrix Multiply(Type Coefficient) const {

		CheckParameters();

		const int Length = this->Length;
		const int Height = this->Height;

		omp_set_num_threads(NUMBER_OF_CORES);

		OmpMatrix<Type> Temporary(Height, Length);

		#pragma omp parallel for schedule(static, NUMBER_OF_CORES)

		for (int FirstIndex = 0; FirstIndex < Height; FirstIndex++) {

			for (int SecondIndex = 0; SecondIndex < Length; SecondIndex++) {

				Temporary.InnerMatrix[FirstIndex][SecondIndex] = InnerMatrix[FirstIndex][SecondIndex] * Coefficient;
			}
		}

		return Temporary;
	};
	OmpMatrix Exponentiate(int Power) const {

		CheckSquareness();

		OmpMatrix<Type> Temporary(*this);

		for (int Index = 1; Index < Power; Index++) {

			Temporary = Temporary.Multiply(*this);
		}

		return Temporary;
	};
	OmpMatrix Transpose() const {

		CheckParameters();

		const int Length = this->Length;
		const int Height = this->Height;

		omp_set_num_threads(NUMBER_OF_CORES);

		OmpMatrix<Type> Temporary(Length, Height);

		#pragma omp parallel for schedule(static, NUMBER_OF_CORES)

		for (int FirstIndex = 0; FirstIndex < Length; FirstIndex++) {

			for (int SecondIndex = 0; SecondIndex < Height; SecondIndex++) {

				Temporary.InnerMatrix[FirstIndex][SecondIndex] = InnerMatrix[SecondIndex][FirstIndex];
			}
		}

		return Temporary;
	};
	OmpMatrix Inverse() const {

		CheckSquareness();

		const int Length = this->Length;
		const int Height = this->Height;

		OmpMatrix<Type> Temporary(*this);

		omp_set_num_threads(NUMBER_OF_CORES);

		OmpMatrix<Type> Solution(Height, Length);

		Solution.FillMainDiagonal(1);

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
	Type Determinant() const {

		CheckSquareness();

		double Value = 1;

		const int Length = this->Length;
		const int Height = this->Height;

		OmpMatrix<Type> Temporary(*this);

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

	OmpMatrix& operator *= (const OmpMatrix& SomeMatrix) {

		return *this = this->Multiply(SomeMatrix);
	};
	OmpMatrix& operator += (const OmpMatrix& SomeMatrix) {

		return *this = this->Add(SomeMatrix);
	};
	OmpMatrix& operator /= (const OmpMatrix& SomeMatrix) {

		return *this = this->Divide(SomeMatrix);
	};
	OmpMatrix& operator -= (const OmpMatrix& SomeMatrix) {

		return *this = this->Subtract(SomeMatrix);
	};
	OmpMatrix& operator = (const OmpMatrix& SomeMatrix) {

		if (this == &SomeMatrix) {

			return *this;
		}

		NUMBER_OF_CORES = SomeMatrix.NUMBER_OF_CORES;

		omp_set_num_threads(NUMBER_OF_CORES);

		if (Length != SomeMatrix.Length || Height != SomeMatrix.Height) {

			this->~OmpMatrix();

			Length = SomeMatrix.Length;
			Height = SomeMatrix.Height;

			InnerMatrix = new Type*[Height];

			#pragma omp parallel for schedule(static, NUMBER_OF_CORES)

			for (int FirstIndex = 0; FirstIndex < Height; FirstIndex++) {

				InnerMatrix[FirstIndex] = new Type[Length];

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
	OmpMatrix& operator *= (const Type Coefficient) {

		return *this = this->Multiply(Coefficient);
	};
	OmpMatrix& operator ^= (const int Power) {

		return *this = this->Exponentiate(Power);
	};

	tuple<OmpMatrix, OmpMatrix> LUFactorization() const {

		CheckSquareness();

		const int Length = this->Length;
		const int Height = this->Height;

		OmpMatrix<Type> Temporary(*this);

		OmpMatrix<Type> Modified(Height, Length);

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

		return tuple<OmpMatrix, OmpMatrix>(Temporary, Modified);
	};
	OmpMatrix CholeskyFactorization() const {

		CheckSquareness();

		const int Length = this->Length;
		const int Height = this->Height;

		OmpMatrix<Type> Temporary(Height, Length);

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
	OmpMatrix UpperFactorization() const {

		CheckParameters();

		const int Length = this->Length;
		const int Height = this->Height;

		OmpMatrix<Type> Temporary(*this);

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

		return Temporary;
	};
	OmpMatrix LowerFactorization() const {

		CheckParameters();

		const int Length = this->Length;
		const int Height = this->Height;

		OmpMatrix<Type> Temporary(*this);

		for (int FirstIndex = Height - 1; FirstIndex >= 0; FirstIndex--) {

			for (int SecondIndex = FirstIndex; SecondIndex > 0; SecondIndex--) {

				double Coefficient = CheckDivisionByZero(Temporary.InnerMatrix[SecondIndex - 1][FirstIndex] / Temporary.InnerMatrix[FirstIndex][FirstIndex]);

				if (Coefficient != 0) {

					for (int AuxiliaryIndex = Length - 1; AuxiliaryIndex >= FirstIndex; AuxiliaryIndex--) {

						Temporary.InnerMatrix[SecondIndex - 1][AuxiliaryIndex] -= Temporary.InnerMatrix[FirstIndex][AuxiliaryIndex] * Coefficient;
					}
				}
			}
		}


		return Temporary;
	};

	OmpMatrix GetColumn(int Position) const {

		CheckParameters();

		const int Height = this->Height;

		OmpMatrix<Type> Temporary(Height);

		for (int Index = 0; Index < Height; Index++) {

			Temporary.InnerMatrix[Index][0] = InnerMatrix[Index][Position];
		}

		return Temporary;
	};
	OmpMatrix GetRow(int Position) const {

		CheckParameters();

		const int Length = this->Length;

		OmpMatrix<Type> Temporary(1, Length);

		for (int Index = 0; Index < Length; Index++) {

			Temporary.InnerMatrix[0][Index] = InnerMatrix[Position][Index];
		}

		return Temporary;
	};
	OmpMatrix GetMainDiagonal() const {

		CheckParameters();

		const int Height = this->Length;

		OmpMatrix<Type> Temporary(Height);

		for (int Index = 0; Index < Height; Index++) {

			Temporary.InnerMatrix[Index][0] = InnerMatrix[Index][Index];
		}

		return Temporary;
	};
	OmpMatrix GetSideDiagonal() const {

		CheckParameters();

		const int Height = this->Height;

		OmpMatrix<Type> Temporary(Height);

		for (int Index = 0; Index < Height; Index++) {

			Temporary.InnerMatrix[Index][0] = InnerMatrix[Index][Length - Index - 1];
		}

		return Temporary;
	};
    Type* operator[] (int Index) const {

		return InnerMatrix[Index];
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

	vector<vector<Type>> ToVector() {

		CheckParameters();

		const int Length = this->Length;
		const int Height = this->Height;

		vector<vector<Type>> Temporary;

		for (int FirstIndex = 0; FirstIndex < Height; FirstIndex++) {

			Temporary.push_back(vector<Type>(Length));

			for (int SecondIndex = 0; SecondIndex < Length; SecondIndex++) {

				Temporary[FirstIndex][SecondIndex] = InnerMatrix[FirstIndex][SecondIndex];
			}
		}

		return Temporary;
	};
	Type** ToPointer() {

		CheckParameters();

		const int Length = this->Length;
		const int Height = this->Height;

		Type** Temporary = new Type*[Height];

		for (int FirstIndex = 0; FirstIndex < Height; FirstIndex++) {

			Temporary[FirstIndex] = new Type[Length];

			for (int SecondIndex = 0; SecondIndex < Length; SecondIndex++) {

				Temporary[FirstIndex][SecondIndex] = InnerMatrix[FirstIndex][SecondIndex];
			}
		}

		return Temporary;
	};

	Type ColumnSum(int Position) {

		Type Value = 0;

		const int Height = this->Height;

		for (int Index = 0; Index < Height; Index++) {

			Value += InnerMatrix[Index][Position];
		}

		return Value;
	};
	Type RowSum(int Position) {

		Type Value = 0;

		const int Length = this->Length;

		for (int Index = 0; Index < Length; Index++) {

			Value += InnerMatrix[Position][Index];
		}

		return Value;
	};
	Type Average() {

		CheckParameters();

		Type Value = 0;

		const int Length = this->Length;
		const int Height = this->Height;

		for (int FirstIndex = 0; FirstIndex < Height; FirstIndex++) {

			for (int SecondIndex = 0; SecondIndex < Length; SecondIndex++) {

				Value += InnerMatrix[FirstIndex][SecondIndex];
			}
		}

		return Value / (Height * Length);
	}
	Type Rate() {

		CheckParameters();

		const int Length = this->Length;
		const int Height = this->Height;

		Type Value = 0;

		for (int FirstIndex = 0; FirstIndex < Height; FirstIndex++) {

			for (int SecondIndex = 0; SecondIndex < Length; SecondIndex++) {

				Value += InnerMatrix[FirstIndex][SecondIndex] * InnerMatrix[FirstIndex][SecondIndex];
			}
		}

		return sqrt(Value);
	};
	Type Max() {

		CheckParameters();

		Type Value = 0;

		const int Length = this->Length;
		const int Height = this->Height;

		for (int FirstIndex = 0; FirstIndex < Height; FirstIndex++) {

			for (int SecondIndex = 0; SecondIndex < Length; SecondIndex++) {

				Value = InnerMatrix[FirstIndex][SecondIndex] > Value ? InnerMatrix[FirstIndex][SecondIndex] : Value;
			}
		}

		return Value;
	};
	Type Min() {

		CheckParameters();

		Type Value = INFINITY;

		const int Length = this->Length;
		const int Height = this->Height;

		for (int FirstIndex = 0; FirstIndex < Height; FirstIndex++) {

			for (int SecondIndex = 0; SecondIndex < Length; SecondIndex++) {

				Value = InnerMatrix[FirstIndex][SecondIndex] < Value ? InnerMatrix[FirstIndex][SecondIndex] : Value;
			}
		}

		return Value;
	};
#pragma endregion
};
