#ifndef VEVPFiniteStrainNS_C
#define VEVPFiniteStrainNS_C 1

#include <stdio.h>
#include <iostream>
#include "VEVPFiniteStrainNS.h"
//#include "numpy/arrayobject.h"
//#include "/usr/include/python3.8/Python.h"
//#include<arrayobject.h>
//#include<Python.h>

using namespace std;

//Description: Allocate memory and a zero vector of dimesion dim
//Author: Mohamed Haddad 
//Created: 16/02/2023
//Copyright: See COPYING file that comes with this distribution 
void MallocVect(const int dim, double** Vect)
{
    *Vect = (double*)calloc(dim,sizeof(double));
}


//Description: Freeing memory allocated for a vector  of dimesion dim
//Author: Mohamed Haddad 
//Created: 16/02/2023
//Copyright: See COPYING file that comes with this distribution 
void FreeVect(double *Vect)
{
    free(Vect);
}


//Description: Copy a vector  of dimesion dim on an another vector 
//Author: Mohamed Haddad 
//Created: 17/02/2023
//Copyright: See COPYING file that comes with this distribution 
void CopyVect(const int dim, double* Vect, double* VectCopy)
{
    int i;
    for(i=0;i<dim;i++){
        VectCopy[i]=Vect[i];
    }
}


//Description: Display components of a vector 
//Author: Mohamed Haddad 
//Created: 17/02/2023
//Copyright: See COPYING file that comes with this distribution 
void DisplayVect(const int dim, double* Vect)
{
    int i;
    for(i=0;i<dim;i++){
        cout << Vect[i] << " ";
    }
    cout << endl;
}


//Description: Allocate memory and a zero matrix of dimesion dim1*dim2
//Author: Mohamed Haddad 
//Created: 16/02/2023
//Copyright: See COPYING file that comes with this distribution 
void MallocMatrix(const int dim1, const int dim2, double*** Matrix)
{
    int i;

    *Matrix=(double**)calloc(dim1,sizeof(double));
    for(i=0;i<dim1;i++){
        (*Matrix)[i]=(double*)calloc(dim2,sizeof(double));
    }
}


//Description: Freeing memory allocated for a matrix  of dimesion dim1*dim2
//Author: Mohamed Haddad 
//Created: 16/02/2023
//Copyright: See COPYING file that comes with this distribution 
void FreeMatrix(const int dim1, double** Matrix)
{
    int i;
    for(i=0;i<dim1;i++){
        free(Matrix[i]);
    }
}


//Description: Copy a matrix on an another matrix
//Author: Mohamed Haddad 
//Created: 16/02/2023
//Copyright: See COPYING file that comes with this distribution 
void CopyMatrix(const int dim1, const int dim2, double** Matrix, double**CopyMatrix)
{
    int i,j;
    for(i=0;i<dim1;i++){
        for(j=0;j<dim2;j++){
            CopyMatrix[i][j]=Matrix[i][j];
        }
    }
}


//Description: Display matrix components
//Author: Mohamed Haddad 
//Created: 16/02/2023
//Copyright: See COPYING file that comes with this distribution 
void DisplayMatrix(const int dim1, const int dim2, double** Matrix)
{
    int i,j;
    for(i=0;i<dim1;i++){
        for(j=0;j<dim2;j++){
            cout << Matrix[i][j] << " ";
        }
        cout << "\n";
    }
    cout << endl;
}


//Description: Allocate memory and a zero 3D matrix of dimesion dim1*dim2*dim3
//Author: Mohamed Haddad 
//Created: 16/02/2023
//Copyright: See COPYING file that comes with this distribution
void MallocMatrix3D(const int dim1, const int dim2, const int dim3, double**** Matrix3D)
{
    int i,j;
    (*Matrix3D)=(double***)calloc(dim1,sizeof(double));
    for(i=0;i<dim1;i++){
        (*Matrix3D)[i]=(double**)calloc(dim2,sizeof(double));
    }

    for(i=0;i<dim1;i++){
        for(j=0;j<dim2;j++){
            (*Matrix3D)[i][j]=(double*)calloc(dim3,sizeof(double));
        }
    }
}


//Description: Freeing memory allocated for a matrix of dimesion dim1*dim2*dim3
//Author: Mohamed Haddad 
//Created: 16/02/2023
//Copyright: See COPYING file that comes with this distribution 
void FreeMatrix3D(const int dim1, const int dim2, double*** Matrix3D)
{
    int i,j;
    for(i=0;i<dim1;i++){
        for(j=0;j<dim2;j++){
            free(Matrix3D[i][j]);
        }
        free(Matrix3D[i]);
    }
}


//Description: Copy a 3D matrix on an 3D 3D another matrix
//Author: Mohamed Haddad 
//Created: 16/02/2023
//Copyright: See COPYING file that comes with this distribution 
void CopyMatrix3D(const int dim1, const int dim2, const int dim3, double*** Matrix3D, double*** CopyMatrix3D)
{
    int i,j,k;
    for(i=0;i<dim1;i++){
        for(j=0;j<dim2;j++){
            for(k=0;k<dim3;k++){
                CopyMatrix3D[i][j][k]=Matrix3D[i][j][k];
            }
        }
    }
}


//Description: Display a 3D matrix
//Author: Mohamed Haddad 
//Created: 16/02/2023
//Copyright: See COPYING file that comes with this distribution 
void DisplayMatrix3D(const int dim1, const int dim2, const int dim3, double*** Matrix3D)
{
    int i,j,k;
    for(i=0;i<dim1;i++){
        cout << "\n[" << i << "]\n";
        for(j=0;j<dim2;j++){
            cout << "\n";
            for(k=0;k<dim3;k++){
                cout << Matrix3D[i][j][k] << " " ;
            }
        }
    }
    cout << endl;
}


//Description: Allocate memory and a mechanical satate with all variables set to zeros
//Author: Mohamed Haddad 
//Created: 16/02/2023
//Copyright: See COPYING file that comes with this distribution
void MallocState(const int nldi, const int ngmi, double s0, STATE* state)
{
    MallocMatrix(3,3,&(state->strn));
    MallocMatrix(3,3,&(state->strs));
    MallocMatrix(3,3,&(state->DefGrad));
    (state->DefGrad)[0][0]=1;(state->DefGrad)[1][1]=1;(state->DefGrad)[2][2]=1;
    MallocMatrix(3,3,&(state->PK1st));
    state->p =0; state->s=s0; state->pMulti=0;
    MallocMatrix(3,3,&(state->Cvp));
    (state->Cvp)[0][0]=1;(state->Cvp)[1][1]=1;(state->Cvp)[2][2]=1;
    MallocMatrix(3,3,&(state->Eve));

    MallocVect(nldi,&(state->SphPartSveViscous));
    MallocMatrix3D(3,3,ngmi,&(state->DevPartSveViscous));   
}


//Description: Freeing memory allocated for mechanical state
//Author: Mohamed Haddad 
//Created: 16/02/2023
//Copyright: See COPYING file that comes with this distribution
void FreeState(const int nldi, const int ngmi, STATE* state)
{
    FreeMatrix(3,(state->strn));
    FreeMatrix(3,(state->strs));
    FreeMatrix(3,(state->DefGrad));
    FreeMatrix(3,(state->PK1st));

    FreeMatrix(3,(state->Cvp));
    FreeMatrix(3,(state->Eve));

    FreeVect(state->SphPartSveViscous);
    FreeMatrix3D(3,3,(state->DevPartSveViscous));
}


//Description: Copy state on another state
//Author: Mohamed Haddad 
//Created: 16/02/2023
//Copyright: See COPYING file that comes with this distribution
void CopyState(const int nldi, const int ngmi, STATE* st, STATE* stCopy)
{
    CopyMatrix(3,3,st->strn,stCopy->strn);
    CopyMatrix(3,3,st->strs,stCopy->strs);
    CopyMatrix(3,3,st->DefGrad,stCopy->DefGrad);
    CopyMatrix(3,3,st->PK1st,stCopy->PK1st);

    (stCopy->p)=st->p; 
    (stCopy->s)=st->s;
    (stCopy->pMulti)=st->pMulti;

    CopyMatrix(3,3,st->Cvp,stCopy->Cvp);
    CopyMatrix(3,3,st->Eve,stCopy->Eve);

    CopyVect(nldi,st->SphPartSveViscous,stCopy->SphPartSveViscous);
    CopyMatrix3D(3,3,ngmi,st->DevPartSveViscous,stCopy->DevPartSveViscous);
}


//Description: Display state
//Author: Mohamed Haddad 
//Created: 16/02/2023
//Copyright: See COPYING file that comes with this distribution
void DisplayState(const int nldi, const int ngmi, STATE* st)
{
    cout << "True strain:\n"; DisplayMatrix(3,3,st->strn); 
    cout << "Cauchy stress:\n"; DisplayMatrix(3,3,st->strs); 
    cout << "Deformation gradient:\n"; DisplayMatrix(3,3,st->DefGrad); 
    cout << "1st P-K stress:\n"; DisplayMatrix(3,3,st->PK1st); 

    cout << "p=" << st->p << "; s=" << st->s << "; p multiplier =" << st->pMulti <<"\n" << endl;

    cout << "Cvp:\n"; DisplayMatrix(3,3,st->Cvp); 
    cout << "Eve:\n"; DisplayMatrix(3,3,st->Eve); 

    cout << "Viscous stress Vol part:\n"; DisplayVect(nldi,st->SphPartSveViscous);
    cout << "Viscous stress Dev part:\n"; DisplayMatrix3D(3,3,ngmi,st->DevPartSveViscous);
}


//Description: Convert Python vector to Cpp vector
//Author: Mohamed Haddad 
//Created: 20/02/2023
//Copyright: See COPYING file that comes with this distribution
void ConvertVectFromPyToCpp(const int nldi, PyObject* PyVect, double* CppVect)
{
    PyArrayObject* pArray = reinterpret_cast<PyArrayObject*>(PyVect);
    double* pData = reinterpret_cast<double*>(PyArray_DATA(pArray));
    npy_intp* pShape = PyArray_SHAPE(pArray);
    int numRows = pShape[0];

    if(numRows!=nldi){
        cout << "Problem while converting python matrix to c++ array. See function ConvertMatrixFromPyToCpp" << endl;
        exit(0);
    }
    
    for(int i=0;i<numRows;i++){
            CppVect[i]= pData[i];
        
    }

    Py_DECREF(pArray);
    Py_DECREF(pShape);
    delete(pData);
    
}


//Description: Convert Python matrix to Cpp matrix
//Author: Mohamed Haddad 
//Created: 20/02/2023
//Copyright: See COPYING file that comes with this distribution
void ConvertMatrixFromPyToCpp(PyObject *PyMatrix, double** CppMatrix)
{
    PyArrayObject* pArray = reinterpret_cast<PyArrayObject*>(PyMatrix);
    double* pData = reinterpret_cast<double*>(PyArray_DATA(pArray));
    npy_intp* pShape = PyArray_SHAPE(pArray);
    int numRows = pShape[0];
    int numCols = pShape[1];

    if(numRows!=3 || numCols!=3){
        cout << "Problem while converting python matrix to c++ array. See function ConvertMatrixFromPyToCpp" << endl;
        exit(0);
    }
    
    for(int i=0;i<numRows;i++){
        for(int j=0;j<numCols;j++){
            CppMatrix[i][j] = pData[i*numCols + j];
        }
    }

    Py_DECREF(pArray);
    Py_DECREF(pShape);
    delete(pData);

}


//Description: Convert Python 3D matrix to Cpp 3D matrix
//Author: Mohamed Haddad 
//Created: 21/02/2023
//Copyright: See COPYING file that comes with this distribution
void Convert3DMatrixFromPyToCpp(const int ngmi, PyObject* Py3DMatrix, double*** Cpp3DMatrix)
{
    PyArrayObject* pArray = reinterpret_cast<PyArrayObject*>(Py3DMatrix);
    //double* pData = reinterpret_cast<double*>(PyArray_DATA(pArray));
    npy_intp* pShape = PyArray_SHAPE(pArray);
    int numRows = pShape[0];
    int numCols = pShape[1];
    int numLayers = pShape[2];

    if(numRows!=3 || numCols!=3 || numLayers!=ngmi){
        cout << "Problem while converting python matrix to c++ array. See function ConvertMatrixFromPyToCpp" << endl;
        exit(0);
    }
    
    for(int i=0;i<numRows;i++){
        for(int j=0;j<numCols;j++){
            for(int k=0;k<numLayers;k++){
                Cpp3DMatrix[i][j][k] = *reinterpret_cast<double*>(PyArray_GETPTR3(pArray, i, j, k));
            }
            
        }
    }

    Py_DECREF(pArray);
    Py_DECREF(pShape);
    //delete(pData);
}


//Description: Convert Python state to C++ state
//Author: Mohamed Haddad 
//Created: 16/02/2023
//Copyright: See COPYING file that comes with this distribution
void ConvertStateFromPyToCpp(const int nldi, const int ngmi, PyObject* PyState, STATE* CppState)
{
    PyObject* Pystrn=PyObject_GetAttrString(PyState,"strn");
    ConvertMatrixFromPyToCpp(Pystrn,CppState->strn);

    PyObject* Pystrs=PyObject_GetAttrString(PyState,"strs");
    ConvertMatrixFromPyToCpp(Pystrs,CppState->strs);

    //PyObject* PyF=PyObject_GetAttrString(PyState,"F");
    //ConvertMatrixFromPyToCpp(PyF,CppState->DefGrad);

    //PyObject* PyP=PyObject_GetAttrString(PyState,"P");
    //ConvertMatrixFromPyToCpp(PyP,CppState->PK1st);
    


    PyObject* PyEve=PyObject_GetAttrString(PyState,"Eve");
    ConvertMatrixFromPyToCpp(PyEve,CppState->Eve);

    PyObject* PyCvp=PyObject_GetAttrString(PyState,"Cvp");
    ConvertMatrixFromPyToCpp(PyCvp,CppState->Cvp);

    PyObject* PySphPartViscous=PyObject_GetAttrString(PyState,"SphPartSveViscous");
    ConvertVectFromPyToCpp(nldi,PySphPartViscous,CppState->SphPartSveViscous);

    PyObject* PyDevPartViscous=PyObject_GetAttrString(PyState,"DevPartSveViscous");
    Convert3DMatrixFromPyToCpp(nldi,PyDevPartViscous,CppState->DevPartSveViscous);

    
    PyObject* Pyp=PyObject_GetAttrString(PyState,"p");
    CppState->p = PyFloat_AsDouble(Pyp);
    
    PyObject* Pys=PyObject_GetAttrString(PyState,"s");
    CppState->s = PyFloat_AsDouble(Pys);

    Py_DECREF(Pystrn);Py_DECREF(Pystrs);
    //Py_DECREF(PyF);Py_DECREF(PyP);
    Py_DECREF(PyEve);Py_DECREF(PyCvp);
    Py_DECREF(PySphPartViscous);Py_DECREF(PyDevPartViscous);
    Py_DECREF(Pyp);Py_DECREF(Pys);

}


//Description: Copy Cpp vector elements to numpy vector elements
//Author: Mohamed Haddad 
//Created: 24/02/2023
//Copyright: See COPYING file that comes with this distribution
void ConvertVectFromCppToPy(const int nldi, double* CppVect, PyObject* PyVect)
{
    PyArrayObject *pArray = reinterpret_cast<PyArrayObject*>(PyVect);
    double *pStrnArray = reinterpret_cast<double*>(PyArray_DATA(pArray));
    for (int i = 0; i < nldi; i++) {
        pStrnArray[i] = CppVect[i];
    }
    //pStrnArray[0]=1;pStrnArray[1]=2;pStrnArray[2]=3;pStrnArray[3]=4;
}


//Description: Copy Cpp array elements to numpy matrix elements
//Author: Mohamed Haddad 
//Created: 24/02/2023
//Copyright: See COPYING file that comes with this distribution
void ConvertMatrixFromCppToPy(double **CppMatrix, PyObject* PyMatrix)
{
    PyArrayObject *pArray = reinterpret_cast<PyArrayObject*>(PyMatrix);
    double (*pStrnArray)[3] = reinterpret_cast<double(*)[3]>(PyArray_DATA(pArray));
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            pStrnArray[i][j] = CppMatrix[i][j];
        }
    }
     
    /*pStrnArray[0][0] =1;pStrnArray[1][1] =2;pStrnArray[2][2] =3;
    pStrnArray[0][1] =4;pStrnArray[0][2] =5;pStrnArray[1][2] =6;
    pStrnArray[1][0] =7;pStrnArray[2][0] =8;pStrnArray[2][1] =9;*/
}


//Description: Copy Cpp array elements to numpy matrix elements
//Author: Mohamed Haddad 
//Created: 24/02/2023
//Copyright: See COPYING file that comes with this distribution
void Convert3DMatrixFromCppToPy(const int ngmi, double*** Cpp3DMatrix, PyObject* Py3DMatrix)
{
    PyArrayObject *pArray = reinterpret_cast<PyArrayObject*>(Py3DMatrix);
    double* pData = reinterpret_cast<double*>(PyArray_DATA(pArray));
    double (*pStrnArray)[3][ngmi]= reinterpret_cast<double(*)[3][ngmi]>(pData);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 4; k++) {
                pStrnArray[i][j][k] = Cpp3DMatrix[i][j][k];
            }
        }
    }

    /*pStrnArray[0][0][0] =1;pStrnArray[1][1][0] =2;pStrnArray[2][2][0] =3;
    pStrnArray[0][1][0] =4;pStrnArray[0][2][0] =5;pStrnArray[1][2][0] =6;
    pStrnArray[1][0][0] =7;pStrnArray[2][0][0] =8;pStrnArray[2][1][0] =9;

    pStrnArray[0][0][3] =1;pStrnArray[1][1][3] =2;pStrnArray[2][2][3] =3;
    pStrnArray[0][1][3] =4;pStrnArray[0][2][3] =5;pStrnArray[1][2][3] =6;
    pStrnArray[1][0][3] =7;pStrnArray[2][0][3] =8;pStrnArray[2][1][3] =9;*/
}


//Description: Convert Cpp state to python class
//Author: Mohamed Haddad 
//Created: 16/02/2023
//Copyright: See COPYING file that comes with this distribution
void ConvertStateFromCppToPy(const int nldi, const int ngmi, STATE* CppState, PyObject* PyState)
{
    PyObject* Pystrn = PyObject_GetAttrString(PyState,"strn");
    ConvertMatrixFromCppToPy(CppState->strn,Pystrn);  

    PyObject* Pystrs = PyObject_GetAttrString(PyState,"strs");
    ConvertMatrixFromCppToPy(CppState->strs,Pystrs);

    PyObject* PyDefGrad = PyObject_GetAttrString(PyState,"F");
    ConvertMatrixFromCppToPy(CppState->DefGrad,PyDefGrad);

    PyObject* PyPK1st = PyObject_GetAttrString(PyState,"P");
    ConvertMatrixFromCppToPy(CppState->PK1st,PyPK1st);

    PyObject* PyCvp = PyObject_GetAttrString(PyState,"Cvp");
    ConvertMatrixFromCppToPy(CppState->Cvp,PyCvp);

    PyObject* PyEve = PyObject_GetAttrString(PyState,"Eve");
    ConvertMatrixFromCppToPy(CppState->Eve,PyEve);


    PyObject* PySphPartSveViscous= PyObject_GetAttrString(PyState,"SphPartSveViscous");
    ConvertVectFromCppToPy(nldi, CppState->SphPartSveViscous,PySphPartSveViscous);

    PyObject* PyDevPartSveViscous= PyObject_GetAttrString(PyState,"DevPartSveViscous");
    Convert3DMatrixFromCppToPy(ngmi, CppState->DevPartSveViscous,PyDevPartSveViscous);


    PyObject* Pyp= PyFloat_FromDouble(CppState->p);
    PyObject_SetAttrString(PyState, "p", Pyp);

    PyObject* Pys= PyFloat_FromDouble(CppState->s);
    PyObject_SetAttrString(PyState, "s", Pys);

    PyObject* PypMulti= PyFloat_FromDouble(5);
    PyObject_SetAttrString(PyState, "pMulti", PypMulti);

    Py_DECREF(Pystrn);Py_DECREF(Pystrs);
    Py_DECREF(PyDefGrad);Py_DECREF(PyPK1st);
    Py_DECREF(PyCvp);Py_DECREF(PyEve);
    Py_DECREF(PySphPartSveViscous);Py_DECREF(PyDevPartSveViscous);
}


//Description: Interface function between the c++ and python codes
//Author: Mohamed Haddad 
//Created: 16/02/2023
//Copyright: See COPYING file that comes with this distribution
//int CppPythonInterfaceConstitutiveBox(const char* InputFileName, STATE* stateTn, double dt, double **FTau, double **PTau, double ****CalgoTau)
int CppPythonInterfaceConstitutiveBox(const char* InputFileName, STATE* stTn, double dt, double** F, STATE *stTau)
{
    Py_Initialize();
    PyObject* PyBasicOperationsModulde = PyImport_ImportModule("BasicOperations");
    PyObject* PyReadInputsFunction = PyObject_GetAttrString(PyBasicOperationsModulde,"ReadInputs");
    PyObject* PyInputFileName = PyUnicode_FromString(InputFileName);
    PyObject* PyArgsForReadInputsFunction=PyTuple_Pack(1,PyInputFileName);

    PyObject* PyResultsOfReadInputsFunction=PyObject_CallObject(PyReadInputsFunction,PyArgsForReadInputsFunction);
     
    
    PyObject* PyMaterialProps= PyTuple_GetItem(PyResultsOfReadInputsFunction,0);
    PyObject* Pyhomogenization = PyTuple_GetItem(PyResultsOfReadInputsFunction,1);
    PyObject* PyLoading = PyTuple_GetItem(PyResultsOfReadInputsFunction,2);
    PyObject* PyStrainCond = PyTuple_GetItem(PyResultsOfReadInputsFunction,3);
    PyObject* PyTime = PyTuple_GetItem(PyResultsOfReadInputsFunction,4);
    

    PyObject* PyStructuresModule = PyImport_ImportModule("Structures");
    PyObject* PyStateClass = PyObject_GetAttrString(PyStructuresModule,"State");
    PyObject* PyArgsForState =PyTuple_Pack(1,PyMaterialProps);
    PyObject* PyStateInstanceTn= PyObject_CallFunction(PyStateClass,"O",PyMaterialProps);

    ConvertStateFromCppToPy(4,4,stTn,PyStateInstanceTn);

    // Call the python constitutive box
    PyObject* PyConstitutiveBoxModule = PyImport_ImportModule("ConstitutiveVEVP");
    PyObject* PyConstitutiveBoxFunction = PyObject_GetAttrString(PyConstitutiveBoxModule,"constitutiveVEVP");
    PyObject* Pydt =  PyFloat_FromDouble(dt);
    PyObject* PyDefGrad = PyObject_GetAttrString(PyStateInstanceTn,"F");
    ConvertMatrixFromCppToPy(F,PyDefGrad);
    
    PyObject* PyPrintStateFunction = PyObject_GetAttrString(PyStateInstanceTn,"printState");
    PyObject* PyResultsofPrintStateFunction=PyObject_CallObject(PyPrintStateFunction,NULL);

    PyObject* PyArgsForConstitutiveBoxFunction=PyTuple_Pack(5,PyMaterialProps,PyLoading,PyStateInstanceTn,Pydt,PyDefGrad);

    PyObject* PyResultsOfConstitutiveVEVPFunction=PyObject_CallObject(PyConstitutiveBoxFunction,PyArgsForConstitutiveBoxFunction);
    
    PyObject* PyStateInstanceTau=PyTuple_GetItem(PyResultsOfConstitutiveVEVPFunction,0);
    PyObject* PyConverged=PyTuple_GetItem(PyResultsOfConstitutiveVEVPFunction,1);    

    ConvertStateFromPyToCpp(4,4,PyStateInstanceTau,stTau);
    int converged = (int)PyLong_AsLong(PyConverged);

    // Freeing
    Py_DECREF(PyBasicOperationsModulde); Py_DECREF(PyReadInputsFunction);Py_DECREF(PyInputFileName);
    Py_DECREF(PyArgsForReadInputsFunction);
    Py_DECREF(PyResultsOfReadInputsFunction);

    Py_DECREF(PyMaterialProps);Py_DECREF(Pyhomogenization);Py_DECREF(PyLoading);Py_DECREF(PyStrainCond);Py_DECREF(PyTime);
    Py_DECREF(PyStructuresModule);Py_DECREF(PyStateClass);Py_DECREF(PyStateInstanceTn);
    Py_DECREF(PyConstitutiveBoxModule);Py_DECREF(PyConstitutiveBoxFunction);

    Py_Finalize();

    return converged ;
}



void ReadInputFile(const char* InputFileName, double *E0, double *nu)
{
    
    PyObject* PyBasicOperationsModulde = PyImport_ImportModule("BasicOperations");
    PyObject* PyReadInputsFunction = PyObject_GetAttrString(PyBasicOperationsModulde,"ReadInputs");
    PyObject* PyInputFileName = PyUnicode_FromString(InputFileName);
    PyObject* PyArgsForReadInputsFunction=PyTuple_Pack(1,PyInputFileName);
    PyObject* PyResultsOfReadInputsFunction=PyObject_CallObject(PyReadInputsFunction,PyArgsForReadInputsFunction);
    
    PyObject* PyMaterialProps= PyTuple_GetItem(PyResultsOfReadInputsFunction,0);
    PyObject* PyE0 = PyObject_GetAttrString(PyMaterialProps,"E0");
    *E0 = PyFloat_AsDouble(PyE0);
    PyObject* Pynu = PyObject_GetAttrString(PyMaterialProps,"nu");
    *nu = PyFloat_AsDouble(Pynu);


    Py_DECREF(PyBasicOperationsModulde); Py_DECREF(PyReadInputsFunction) ; Py_DECREF(PyInputFileName); 
    Py_DECREF(PyArgsForReadInputsFunction);
    Py_DECREF(PyResultsOfReadInputsFunction);Py_DECREF(PyMaterialProps);Py_DECREF(PyE0); Py_DECREF(Pynu);
   
}

#endif