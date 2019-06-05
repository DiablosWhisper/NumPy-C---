%module MatrixLibrary
%{
	#include "MatrixLibrary.h"
%}
using namespace std;
%include "std_string.i"
%include "MatrixLibrary.h"
%template(FMatrix) Matrix<float>;