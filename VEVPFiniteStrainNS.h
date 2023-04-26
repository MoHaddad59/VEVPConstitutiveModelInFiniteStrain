#ifndef VEVPFINITESTRAINNS_H
#define VEVPFINITESTRAINNS_H 1

#include<arrayobject.h>
#include<Python.h>


typedef struct MaterialProps
{
    double E0;
    double nu;
    double yield_strs;

    // VE props
    int nldi;
    double ld0;
    double ldinf;
    double* ldi;
    double* li;

    int ngmi;
    double gm0;
    double gminf;
    double gmi;
    double gi;

    // Yield function 
    int YieldType;
    double keq;
    double kh;
    double q;

    // Hardening Law
    int HardType;
    double HardMod1;
    double HardMod2;
    double HardExp;

    // VP function
    int VPFuncType;
    double pDot0;
    double A;
    double T;
    double alpha;
    double m;

    // VP potential
    int VpPotType;
    double TanPhi;
    double ksi;
    double nup;
    double kappa;

    // Softening law
    int SoftType;
    double s0, s1, s2;
    double h1,h2;
    double ps;
    double k;
    
    // Rehardening law
    double RehardType;
    double CR;
    double N;

} MATERIALPROPS;

typedef struct State
{
    double **strn, **strs;
    double **DefGrad, **PK1st;
    double p,s,pMulti;
    double **Cvp, **Eve;
    double *SphPartSveViscous;
    double ***DevPartSveViscous;
}STATE;

void MallocVect(const int dim, double** Vect);
void FreeVect(double *Vect);
void CopyVect(const int dim, double* Vect, double* VectCopy);
void DisplayVect(const int dim, double* Vect);

void MallocMatrix(const int dim1, const int dim2, double*** Matrix);
void FreeMatrix(const int dim1, double** Matrix);
void CopyMatrix(const int dim1, const int dim2, double** Matrix, double**CopyMatrix);
void DisplayMatrix(const int dim1, const int dim2, double** Matrix);

void MallocMatrix3D(const int dim1, const int dim2, const int dim3, double**** Matrix3D);
void FreeMatrix3D(const int dim1, const int dim2, double*** Matrix3D);
void CopyMatrix3D(const int dim1, const int dim2, const int dim3, double*** Matrix3D, double*** CopyMatrix3D);
void DisplayMatrix3D(const int dim1, const int dim2, const int dim3, double*** Matrix3D);

void MallocState(const int nldi, const int ngmi, double s0, STATE* state);
void FreeState(const int nldi, const int ngmi, STATE* state);
void CopyState(const int nldi, const int ngmi, STATE* st, STATE* stCopy);
void DisplayState(const int nldi, const int ngmi, STATE* st);

void ConvertVectFromPyToCpp(const int nldi,PyObject* PyVect, double* CppVect);
void ConvertMatrixFromPyToCpp(PyObject *PyMatrix, double** CppMatrix);
void Convert3DMatrixFromPyToCpp(const int ngmi, PyObject* Py3DMatrix, double*** Cpp3DMatrix);
void ConvertStateFromPyToCpp(const int nldi, const int ngmi, PyObject* PyState, STATE* CppState);

void ConvertVectFromCppToPy(const int nldi, double* CppVect, PyObject* PyVect);
void ConvertMatrixFromCppToPy(double **CppMatrix, PyObject* PyMatrix);
void Convert3DMatrixFromCppToPy(const int ngmi, double*** Cpp3DMatrix, PyObject* Py3DMatrix);
void ConvertStateFromCppToPy(const int nldi, const int ngmi, STATE* CppState, PyObject* PyState);

int CppPythonInterfaceConstitutiveBox(const char* InputFileName, STATE* stTn, double dt, double** F, STATE *stTau);
void ReadInputFile(const char* InputFileName, double *E0, double *nu);
 

#endif