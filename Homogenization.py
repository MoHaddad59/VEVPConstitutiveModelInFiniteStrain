import ConstitutiveVEVP as cb
import Structures as st
import VEVPOperations as vevp
import numpy as np
import sys
from csv import DictWriter
import csv  
from scipy.optimize import newton, fsolve

def HomogeneousMaterialVEVPFiniteStrain(PhaseProp,Load,StrainCond,Time,OutputFile):
    
    Statetn=st.State(PhaseProp)
    Statetau=st.State(PhaseProp)
    globtol=1e-3
    globmaxiter=1000
    
    #Initial state to be done
    #????
    
    # t=0
    tn=0.0
    ColNames=['t','logstrn11','logstrn22','logstn33','Cauchystrs11','Cauchystrs22','Cauchystrs33','logstrn12','logstrn13','logstn23','Cauchystrs12','Cauchystrs13','Cauchystrs23','pMultiplier','p','s','F11','F22','F33','Yeff11','Yeff22','Yeff33','F12','F13','F23','F21','F31','F32','Yeff12','Yeff13','Yeff23','Yeff21','Yeff31','Yeff32']    
    with open(OutputFile,mode='w',encoding="utf8") as Outfile:
        writer=csv.writer(Outfile)
        writer.writerow(ColNames)
        Outfile.close()
        
    Statetn.printStateOnTxtFile(OutputFile,tn)
    
    #t>0
    dF=np.zeros((3,3));Ftau=np.zeros((3,3))
    dt=Time.getTimeStep()
    tau=tn+dt
    
    if(Load.loading=="SHEAR_12"):
        while(tau<=Time.tf):
            dt=Time.getTimeStep()
            dF=vevp.DefGradientIncrement(Load,Time,StrainCond)
            Ftau=(Statetn.F)+dF
            
            Statetau, ok=cb.constitutiveVEVP(PhaseProp,Load,Statetn,dt,Ftau)
            
            iterRed=1
            while(ok==False):
                print("%i th redcution of the time step" %iterRed)
                dt/=2
                Statetau, ok=cb.constitutiveVEVP(PhaseProp,Load,Statetn,dt,Ftau)
                
                if (dt<1.e-6):
                    print("Very small time increment attended dt=%g. Exiting..."%dt)
                    sys.exit()
                iterRed+=1
                
            Statetau.DevPartSveViscous,Statetau.SphPartSveViscous,Statetau.strn=vevp.UpdateVariables(PhaseProp,Statetn,Statetau,dt,Ftau)  
            Statetau.F = Ftau
            Statetau.copyStatetIn(Statetn)
            Statetn.printStateOnTxtFile(OutputFile,tau)
            print("[%g]"%tau)
            #Update t
            tau+=dt
    
    elif(Load.loading=="UNIAXIAL_1"):
        while(tau<=Time.tf):
            dt=Time.getTimeStep()
            Ftau=vevp.DeformationGradienttau(Statetn,Load,Time,StrainCond)
            Statetau, ok=cb.constitutiveVEVP(PhaseProp,Load,Statetn,dt,Ftau)
            
            iterRed=1
            while(ok==False):
                print("%i th redcution of the time step" %iterRed)
                dt/=2
                Statetau, ok=cb.constitutiveVEVP(PhaseProp,Load,Statetn,dt,Ftau)
                
                if (dt<1.e-6):
                    print("Very small time increment attended dt=%g. Exiting..."%dt)
                    sys.exit()
                iterRed+=1
                
            '''if(vevp.GlobalResidual(Load,Statetau.P)>globtol):
                Fopt=np.zeros(9,dtype=float)
                Fopt=vevp.F3DTo2D(Ftau)
                print(Fopt) 
                solglob = newton(func=cb.GlobalNR,x0=Fopt,args=(PhaseProp,Load,Statetn,dt),tol=globtol,maxiter=globmaxiter,full_output=True)
                #print(solglob)
                Ftau=vevp.ExtractFromGlobSol(solglob.root)
                
                solglob=fsolve(func=cb.GlobalNR,x0=Fopt,args=(PhaseProp,Load,Statetn,dt),maxfev=globmaxiter,xtol=globtol,full_output=True)
                print((solglob[1])['fvec'])
                Ftau=vevp.ExtractFromGlobSol(solglob[0])
                
                
                Statetau, ok=cb.constitutiveVEVP(PhaseProp,Load,Statetn,dt,Ftau)
                #print("Residual:%g" %vevp.GlobalResidual(Load,Statetau.P)) 
                del Fopt,solglob'''
            
            Statetau.DevPartSveViscous,Statetau.SphPartSveViscous,Statetau.strn=vevp.UpdateVariables(PhaseProp,Statetn,Statetau,dt,Ftau)  
            Statetau.F = Ftau
            Statetau.copyStatetIn(Statetn)
            Statetn.printStateOnTxtFile(OutputFile,tau)
            print("[%g]"%tau)
            #Update t
            tau+=dt
    
    del dF,Ftau,Statetau,Statetn,tn,dt,tau, globtol,globmaxiter
