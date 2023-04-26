import BasicOperations as bo
import VEVPOperations as vevp
import Structures  as st
from numpy.linalg import det,inv
from scipy.linalg import fractional_matrix_power
import numpy as np
import sys
from scipy.optimize import newton


def constitutiveVEVP(PhaseProp,Load,Statetn,dt,F):
    loctol=1.e-8
    locmaxiter=100
    Statetau=st.State(PhaseProp)
    converged=0
    try:
        J=det(F)
    except:
        print("Problems occuured when computing det(F) in constitutiveVEVP line=%i"%(bo.CurrLineNumber()))
        print("F=\n",F,"J=%g"%J)
        sys.exit()
    
    Cvppred=Statetn.Cvp
    Evepred,Svepred,Spred,Ypred,Yeffpred=vevp.StrainStressMeasuresFromCvp(PhaseProp,Statetn,dt,F,Cvppred)
    fpred=vevp.YieldFunction(PhaseProp,Load,Ypred,Yeffpred,Statetn.p)
    if(fpred<=0):
        Statetau.Cvp=Statetn.Cvp
        
        Statetau.p=Statetn.p
        Statetau.s=Statetn.s
        Statetau.pMulti=Statetn.pMulti
        
        Statetau.Sve=Svepred
        Statetau.S=Spred
        Statetau.Y=Ypred
        Statetau.Yeff=Yeffpred
        Statetau.strs=(1/J)*np.dot(F,np.dot(Spred,F.transpose()))
        Statetau.P=np.dot(F,Spred)
        converged=1
    else:
        try:
            Uvpn= fractional_matrix_power(Statetn.Cvp,1/2)
        except:
            print("Problems occured when computing Cvp**(1/2) in constitutiveVEVP line=%i"%(bo.CurrLineNumber()))
            print("Cvpn=\n",Statetn.Cvp,"Uvpn=\n",Uvpn)
            sys.exit()
            
        try:
            invUvpn=inv(Uvpn)
        except:
            print("Problems occured when computing inv(Uvp) in constitutiveVEVP line=%i"%(bo.CurrLineNumber()))
            print("invUvpn=\n",invUvpn)
            sys.exit()
            
        Soliter0=np.zeros(7,dtype=float)
        Soliter0[0]=invUvpn[0,0];Soliter0[1]=invUvpn[0,1];Soliter0[2]=invUvpn[0,2]
        Soliter0[3]=invUvpn[1,1];Soliter0[4]=invUvpn[1,2];Soliter0[5]=invUvpn[2,2];
        #Soliter0[6]=Statetn.p
        Soliter0[6]=Statetn.pMulti
        
        #print(Soliter0)
        
        sol = newton(func=vevp.LocalNR,x0=Soliter0,args=(PhaseProp,Load,Statetn,dt,F),tol=loctol,maxiter=locmaxiter,full_output=True)
        #print("NR solution=\n",sol.root)
        #print("NR convergence=\n",sol.converged)
        flag=1
        for ii in range(0,len(sol.converged)):
            if sol.converged[ii]==False:
                flag=0
                
        if flag==1:
            converged=True
        
        #print("conveged=%i",converged)
        #print("Sol=" ,sol.root)
        Statetau.Cvp, Statetau.pMulti, Statetau.p, Statetau.s = vevp.ExtractFromSol(PhaseProp,Statetn,dt,F,sol.root)
        Statetau.Eve,Statetau.Sve,Statetau.S,Statetau.Y,Statetau.Yeff=vevp.StrainStressMeasuresFromCvp(PhaseProp,Statetn,dt,F,Statetau.Cvp)
        Statetau.strs=(1/J)*np.dot(F,np.dot(Statetau.S,F.transpose()))  
        Statetau.P = np.dot(F,(Statetau.S))
    return Statetau, converged



def GlobalNR(Fopt,PhaseProp,Load,Statetn,dt):
    statepert=st.State(PhaseProp)
    DefGrad=np.zeros((3,3))
    ok=False
    
    DefGrad[0,0]=Fopt[0];DefGrad[0,1]=Fopt[1];DefGrad[0,2]=Fopt[2]
    DefGrad[1,0]=Fopt[3];DefGrad[1,1]=Fopt[4];DefGrad[1,2]=Fopt[5]
    DefGrad[2,0]=Fopt[6];DefGrad[2,1]=Fopt[7];DefGrad[2,2]=Fopt[8]
    
    statepert,ok=constitutiveVEVP(PhaseProp,Load,Statetn,dt,DefGrad)
    
    '''iterRed=1
    while(ok==False):
        print("%i th redcution of the time step" %iterRed)
        dt/=2
        statepert, ok=constitutiveVEVP(PhaseProp,Load,Statetn,dt,F)
        
        if (dt<1.e-6):
            print("Very small time increment attended dt=%g. Exiting..."%dt)
            sys.exit()
        iterRed+=1'''
    
    Fx=np.zeros(9,dtype=float)
    Fx[0]= 0
    Fx[1]=(statepert.strs)[0,1]
    Fx[2]=(statepert.strs)[0,2]
    Fx[3]=(statepert.strs)[1,0]
    Fx[4]=(statepert.strs)[1,1]
    Fx[5]=(statepert.strs)[1,2]
    Fx[6]=(statepert.strs)[2,0]
    Fx[7]=(statepert.strs)[2,1]
    Fx[8]=(statepert.strs)[2,2]
    
    del DefGrad, statepert,ok
    return Fx