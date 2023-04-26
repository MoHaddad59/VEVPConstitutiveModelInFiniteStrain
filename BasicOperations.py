import Structures as st
import pandas as pd
import sys
import numpy as np
import math as mt
from inspect import currentframe

def ReadInputs(InputFileName):
    file=open(InputFileName,"r")
    lines=file.readlines()
    inputs=[]
    for ll in lines:
       inputs.append(ll.split(' ')[0])
    file.close()
    
    line=10
    nldi = int(inputs[line-1]) 
    line=20
    ngmi = int(inputs[line-1])
    
    matrix =st.PhaseProp(nldi,ngmi)
    matrix.nldi=nldi
    matrix.ngmi=ngmi
    
    line=6
    matrix.E0 = float(inputs[line-1])
    line=7
    matrix.nu = float(inputs[line-1])
    
    #matrix.gm0=matrix.computeGm0(matrix.E0,matrix.nu)
    #matrix.ld0=matrix.computeLd0(matrix.E0,matrix.nu)
    matrix.computeGm0()
    matrix.computeLd0()
    
    line=8
    matrix.strsy = float(inputs[line-1])
    line=11
    for ildi in range(0,nldi):
        matrix.ldi[ildi] = float(inputs[line-1 + ildi])
        matrix.li[ildi]= float(inputs[line-1 + matrix.nldi + ildi])
    line=21
    for igmmi in range(0,ngmi):
        matrix.gmi[igmmi] = float(inputs[line-1 + igmmi])
        matrix.gi[igmmi]= float(inputs[line-1 + matrix.ngmi + igmmi])
    
    matrix.computeGminf()
    matrix.ComputeLdinf()

    line=30
    matrix.YieldType = float(inputs[line-1])
    line=31
    matrix.keq= float(inputs[line-1])
    line=32
    matrix.kh = float(inputs[line-1])
    line=33
    matrix.q = float(inputs[line-1])
  
    line=35
    matrix.HardType= float(inputs[line-1])
    line=36
    matrix.HardMod1= float(inputs[line-1])
    line=37
    matrix.HardExp = float(inputs[line-1])
    line=38
    matrix.HardMod2 = float(inputs[line-1])  
    
    line=40
    matrix.VpPotType= float(inputs[line-1])
    line=41
    matrix.ksi= float(inputs[line-1])
    line=42
    matrix.nup = float(inputs[line-1])
    matrix.computeTanPhi()
    matrix.computekappa()
    
    line=44
    matrix.VpFuncType= float(inputs[line-1])  
    line=45
    matrix.pDot0= float(inputs[line-1])    
    line=46
    matrix.T= float(inputs[line-1]) 
    line=47
    matrix.A= float(inputs[line-1]) 
    line=48
    matrix.alpha= float(inputs[line-1]) 
    line=49
    matrix.m= float(inputs[line-1])     

    line=51
    matrix.SoftType= float(inputs[line-1]) 
    line=52
    matrix.s0= float(inputs[line-1]) 
    line=53
    matrix.s1= float(inputs[line-1]) 
    line=54
    matrix.s2= float(inputs[line-1]) 
    line=55
    matrix.h1= float(inputs[line-1]) 
    line=56
    matrix.h2= float(inputs[line-1]) 
    line=57
    matrix.ps= float(inputs[line-1]) 
    line=58
    matrix.k= float(inputs[line-1]) 

    line=60
    matrix.RehardType= float(inputs[line-1]) 
    line=61
    matrix.CR= float(inputs[line-1]) 
    line=62
    matrix.N= float(inputs[line-1])     
    
    matrix.printPhaseProp()
    
    homo=st.Homogenization()
    line=66
    homo.activated= str(inputs[line-1])
    if homo.activated not in ["ON","OFF"]:
        print("Invalid input in line %i"%line)
        sys.exit()
        
    line=67
    homo.scheme= str(inputs[line-1])
    if homo.scheme not in ["MORI-TANAKA"]:
        print("Invalid input in line %i"%line)
        sys.exit()
        
    line=68
    homo.isotrop= str(inputs[line-1])
    if homo.isotrop not in ["ESH-ISO"]:
        print("Invalid input in line %i"%line)
        sys.exit()
    
    line=69
    homo.ratiodt= float(inputs[line-1])
    
    line=70
    homo.StatMoment= int(inputs[line-1])
    if homo.StatMoment not in [1,2]:
        print("Invalid input in line %i"%line)
        sys.exit()
    
    homo.printHomo()
    
    loading=st.Loading()
    line=72
    loading.history= str(inputs[line-1])
    if loading.history not in ["MONOTONIC","CYCLIC"]:
        print("Invalid input in line %i"%line)
        sys.exit()
    
    line=73
    loading.testtype = str(inputs[line-1])
    if loading.testtype not in ["COMPRESSION","TENSILE"]:
        print("Invalid input in line %i"%line)
        sys.exit()
    
    line=74
    loading.loading= str(inputs[line-1])
    if loading.loading not in ["UNIAXIAL_1","BIAXIAL_12","SHEAR_12", "TRIAXIAL", "RELAXATION_1","CREEP_1"]:
        print("Invalid input in line %i"%line)
        sys.exit()
        
    line=75
    loading.triaxiality= float(inputs[line-1])
    
    line=76
    loading.nCycles= int(inputs[line-1])
    
    loading.printLoading()
    
    strncond=st.StrainCond()
    line=78
    strncond.strs0=float(inputs[line-1])
    line=79
    strncond.strn0=float(inputs[line-1])
    line=81
    strncond.StrainRate=float(inputs[line-1])
    line=82
    truestrnf=float(inputs[line-1])
    
    strncond.strnf=mt.exp(truestrnf)-1
    
    time=st.Time()
    line=83
    time.nsteps=float(inputs[line-1])
    
    # Final time computed from the applied strain rate and final engineering strain
    time.tf=(strncond.strnf)/(strncond.StrainRate)
    
    strncond.printStrainCond()
    time.printTime()
    
    return matrix, homo, loading, strncond, time

def DevPartMadelStress(Y):
        return(Y-((1/3)*np.trace(Y))*(np.identity(3)))

def EqMandelStress(Y):
    Ydev=DevPartMadelStress(Y)
    return(mt.sqrt((3/2)*np.tensordot(Ydev,Ydev,axes=2)))

def invL(arg):
    return(3*((3-mt.pow(arg,2))/(1-mt.pow(arg,2))))

def computeH1(PhaseProp,ptau):
    angle= (ptau-(PhaseProp.ps))/((PhaseProp.k)*PhaseProp.ps)
    return(-(PhaseProp.h1)*(mt.tanh(angle)-1))    

def computeH2(PhaseProp,ptau):
    angle=(ptau-(PhaseProp.ps))/((PhaseProp.k)*PhaseProp.ps)
    return((PhaseProp.h2)*(mt.tanh(angle)+1))

def SymmetricPart(Y):
    return((1/2)*(Y+Y.transpose()))

def CurrLineNumber():
    cf = currentframe()
    return cf.f_back.f_lineno
    



        

    
    