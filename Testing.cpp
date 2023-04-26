#include <stdio.h>
#include <iostream>
#include "VEVPFiniteStrainNS.h"

using namespace std;

int main()
{
   
    Py_Initialize();
    /*const char* InputFileName="InputFile.txt";
    PyObject* PyStructuresModule=PyImport_ImportModule("Structures");
    PyObject* PyBasicOperationsModule=PyImport_ImportModule("BasicOperations");
    PyObject* PyConstitutiveVEVPModule=PyImport_ImportModule("ConstitutiveVEVP");
    PyObject* PyHomogenizationBox=PyImport_ImportModule("Homogenization");
   
    if(PyBasicOperationsModule==NULL){
        PyErr_Print();
        cout << "I am here \n" << endl;
    }// Read input file
    PyObject* PyReadInputsFunc = PyObject_GetAttrString(PyBasicOperationsModule,"ReadInputs");
    
    PyObject* PyInputFileNameAsArg=PyUnicode_FromString(InputFileName);
    PyObject* PyArgForReadInputsFunc=PyTuple_Pack(1,PyInputFileNameAsArg);

    PyObject* PyResultsOfReadInputsFunc=PyObject_CallObject(PyReadInputsFunc,PyArgForReadInputsFunc);
     
    // call homgenization box
    PyObject* Pymatrix = PyTuple_GetItem(PyResultsOfReadInputsFunc,0);
    PyObject* Pyhomogenization = PyTuple_GetItem(PyResultsOfReadInputsFunc,1);
    PyObject* Pyloading = PyTuple_GetItem(PyResultsOfReadInputsFunc,2);
    PyObject* PystrainCond = PyTuple_GetItem(PyResultsOfReadInputsFunc,3);
    PyObject* Pytime = PyTuple_GetItem(PyResultsOfReadInputsFunc,4);

    // Test different conversion functions 
    PyObject* PyTestingModule=PyImport_ImportModule("Testing");
    PyObject* PyTest1Func = PyObject_GetAttrString(PyTestingModule, "Test1");
    PyObject* PyTest1FuncResults=PyObject_CallObject(PyTest1Func,NULL);
    double **CppMatrix; MallocMatrix(3,3,&CppMatrix);
    ConvertMatrixFromPyToCpp(PyTest1FuncResults,CppMatrix);
    cout << "\n CppMatrix:\n" ;
    DisplayMatrix(3,3,CppMatrix);

    PyObject* PyTest2Func = PyObject_GetAttrString(PyTestingModule, "Test2");
    PyObject* PyTest2FuncResults=PyObject_CallObject(PyTest2Func,NULL);
    double *CppVect; MallocVect(4,&CppVect);
    ConvertVectFromPyToCpp(4,PyTest2FuncResults, CppVect);
    cout << "\n CppVect:\n";
    DisplayVect(4,CppVect);


    PyObject* PyTest3Func = PyObject_GetAttrString(PyTestingModule, "Test3");
    PyObject* PyTest3FuncResults=PyObject_CallObject(PyTest3Func,NULL);
    double ***Cpp3DMatrix; MallocMatrix3D(3,3,4,&Cpp3DMatrix);
    Convert3DMatrixFromPyToCpp(4,PyTest3FuncResults,Cpp3DMatrix);
    cout << "\n Cpp3DMatrix:\n";
    DisplayMatrix3D(3,3,4,Cpp3DMatrix);

    // Create state class
    PyObject* PyArgForState =PyTuple_Pack(1,Pymatrix);
    PyObject* PyStateClass= PyObject_GetAttrString(PyStructuresModule,"State");
    PyObject* PyState= PyObject_CallFunction(PyStateClass,"O",Pymatrix);
    STATE st; MallocState(4,4,&st);
    ConvertStateFromPyToCpp(4,4,PyState,&st);
    DisplayState(4,4,&st);*/
    double **F; MallocMatrix(3,3,&F); F[0][0]=1+1.e-4;F[1][1]=1;F[2][2]=1;
    double dt = 1e-5;
    STATE stTn; MallocState(4,4,107e6,&stTn);
    STATE stTau; MallocState(4,4,107e6,&stTau);
    int ok;
    ok=CppPythonInterfaceConstitutiveBox("../InputFile.txt",&stTn,dt,F,&stTau);
    DisplayState(4,4,&stTau);

    cout << "convergd = " << ok << endl;

    /*double E0, nu;
    ReadInputFile("InputFile.txt",&E0,&nu);
    cout << "E0=" << E0 << "  nu=" << nu << endl;*/

    
    Py_Finalize();
    return 0;
    
}