import math as mt
import sys
import numpy as np
from numpy.linalg import det,inv,eig
from scipy.linalg import expm
import BasicOperations as bo
from scipy.linalg import fractional_matrix_power
from inspect import currentframe

def HardeningLaw(PhaseProp,p):
    # Case 1: R(p) = k1(1-exp(-n*p))
    if(PhaseProp.HardType==1):
        Rp=0
        Rp=(PhaseProp.HardMod1)*(1-mt.exp(-(PhaseProp.HardExp)*p))
        return Rp
    
    else:
        print("Invalid hardening model type in InputFile line 35")
        sys.exit()
    
def YieldFunction(PhaseProp,Load,Y,Yeff,p):
    # Case 1: Dracker-Prager yield function 
    if(PhaseProp.YieldType ==1):
        f=0; Yeffeq=0;Rp=0;Ym=0;StrsPBar=0;StrsTP=0
        Ym = (1/3)*np.trace(Y)
        Yeffeq= bo.EqMandelStress(Yeff)
        Rp=HardeningLaw(PhaseProp,p)
        StrsPBar=(PhaseProp.strsy) + Rp
        
        if(Load.loading=="UNIAXIAL_1"):
            if(Load.testtype=="COMPRESSION"):
                try:
                    StrsTP=(PhaseProp.keq)*mt.pow((StrsPBar/(PhaseProp.strsy)),(PhaseProp.q)) - ((PhaseProp.kh)*(StrsPBar/3))
                except:
                    print("Problems when computing StrsTP in YieldFunction in VEVPOperations line %i\n"%(bo.CurrLineNumber()))
                    print("p=%g, StrsPBar=%g, StrsPBar/(PhaseProp.strsy)=%g, StrsPBar=%g\n"%(p,StrsPBar,StrsPBar/(PhaseProp.strsy),StrsPBar))
            elif(Load.testtype=="TENSILE"):
                try:
                    StrsTP=(PhaseProp.keq)*mt.pow((StrsPBar/(PhaseProp.strsy)),(PhaseProp.q)) + ((PhaseProp.kh)*(StrsPBar/3))
                except:
                    print("Problems when computing StrsTP in YieldFunction in VEVPOperations line %i\n"%(bo.CurrLineNumber()))
                    print("p=%g, StrsPBar=%g, StrsPBar/(PhaseProp.strsy)=%g, StrsPBar=%g\n"%(p,StrsPBar,StrsPBar/(PhaseProp.strsy),StrsPBar))
                
        elif(Load.loading=="SHEAR_12"):
            try:
                StrsTP=(PhaseProp.keq)*mt.pow((StrsPBar/(PhaseProp.strsy)),(PhaseProp.q))- ((PhaseProp.kh)*(StrsPBar/3))
                #StrsTP=(PhaseProp.keq)*mt.pow((StrsPBar/(PhaseProp.strsy)),(PhaseProp.q))
            except:
                print("Problems when computing StrsTP in YieldFunction in VEVPOperations line %i\n"%(bo.CurrLineNumber()))
                print("p=%g, StrsPBar=%g, StrsPBar/(PhaseProp.strsy)=%g, StrsPBar=%g\n"%(p,StrsPBar,StrsPBar/(PhaseProp.strsy),StrsPBar))
        else:
            print("Invalid StrsTP in function YieldFunction in VEVPOperations.py line =%i"%(bo.CurrLineNumber()))
            sys.exit()
            
        f = (PhaseProp.keq)*mt.pow((Yeffeq/(PhaseProp.strsy)),(PhaseProp.q)) - ((PhaseProp.kh)*Ym) - StrsTP
        
        del Ym, Yeffeq, Rp, StrsPBar, StrsTP
        return f
    
    # Case 2: Ngueyen yield function type
    elif(PhaseProp.YieldType ==2):
        f=0; Yeffeq=0;Rp=0;Ym=0;StrsPBar=0
        Ym = (1/3)*np.trace(Y)
        Yeffeq= bo.EqMandelStress(Yeff)
        Rp=HardeningLaw(PhaseProp,p)
        StrsPBar=(PhaseProp.strsy) + Rp
        
        m = PhaseProp.kh  # kh is m and q is alpha for this yield function type
        a2=1/(mt.pow(StrsPBar,PhaseProp.q)); a1=3*((mt.pow(m,PhaseProp.q)-1)/(m+1))*(1/StrsPBar); a0=(mt.pow(m,PhaseProp.q)+m)/(m+1)
        
        f= a2*mt.pow(Yeffeq,PhaseProp.q) - a1*Ym - a0
        
        del Yeffeq,Rp,Ym,StrsPBar,m,a0,a1,a2
        return f
    
    else:
        print("Invalid yield function type in InputFile line 30")
        sys.exit()
    
def VPPotential(PhaseProp,Y,Yeff):
    # Case 1: VP Potential from Murali et Doghri 2018 (Dracker-Prager)
    if(PhaseProp.VpPotType==1):
        Ym=0;Yeffeq=0;Yeffdev=np.zeros((3,3));tmp1=0;g=0;DgDYm=0;DgDYeffeq=0;DgDY=np.zeros((3,3))
        Ym= (1/3)*np.trace(Y)
        Yeffeq=bo.EqMandelStress(Yeff)
        Yeffdev=bo.DevPartMadelStress(Yeff)
    
        tmp1=(PhaseProp.ksi)*(PhaseProp.strsy)*(PhaseProp.TanPhi)
        
        g=mt.sqrt(mt.pow(tmp1,2)+mt.pow(Yeffeq,2)) - (PhaseProp.TanPhi)*Ym
        
        DgDYm=-PhaseProp.TanPhi
        DgDYeffeq=Yeffeq/mt.sqrt(mt.pow(tmp1,2)+mt.pow(Yeffeq,2))
        DgDY=((3/2)*DgDYeffeq*(1/Yeffeq))*Yeffdev + ((1/3)*DgDYm)*np.identity(3)
        
        del Ym,Yeffeq,Yeffdev,tmp1,DgDYm,DgDYeffeq
        return g, DgDY.transpose()
    
    # Case 2: VP potential develpped in Nguyen et al. 2016
    elif(PhaseProp.VpPotType==2):
        Ym=0;Yeffeq=0;Yeffdev=np.zeros((3,3));tmp1=0;g=0;DgDYm=0;DgDYeffeq=0;DgDY=np.zeros((3,3))
        Ym=(1/3)*np.trace(Y)
        Yeffeq=bo.EqMandelStress(Yeff)
        Yeffdev=bo.DevPartMadelStress(Yeff)
        
        #beta= (9- 18*(PhaseProp.nup))/(2*(1+(PhaseProp.nup)))
        
        g= mt.pow(Yeffeq,2) - (PhaseProp.TanPhi)*mt.pow(Ym,2)
        
        DgDYm= - 2*(PhaseProp.TanPhi)*Ym
        DgDYeffeq=2*Yeffeq
        DgDY=((3/2)*DgDYeffeq*(1/Yeffeq))*Yeffdev + ((1/3)*DgDYm)*np.identity(3)
        
        del Ym,Yeffeq,Yeffdev,DgDYm,DgDYeffeq
        return g, DgDY.transpose()
         
    else:
        print("Invalid VP potential type in InputFile line 40")
        sys.exit()
        
def VPFunction(PhaseProp,Load,Y,Yeff,ptau,stau):    
    if(PhaseProp.VpFuncType==1):
        f=0;Rp=0;Ym=0;Yeffeq=0;tmp1=0;prt1=0;prt2=0; gv=0
        f=YieldFunction(PhaseProp,Load,Y,Yeff,ptau)
        Rp=HardeningLaw(PhaseProp,ptau)
        Ym=(1/3)*np.trace(Y)
        Yeffeq=bo.EqMandelStress(Yeff)
        tmp1= stau - (PhaseProp.alpha)*Ym
        
        try:
            prt1=-((PhaseProp.A)/(PhaseProp.T))*tmp1
            prt2=1-mt.pow((Yeffeq/tmp1),(PhaseProp.m))
            gv=(PhaseProp.pDot0)*mt.exp(prt1*prt2)
        except:
            print("Problems when computing gv in VPFunction in VEVPOperations line %i\n"%(bo.CurrLineNumber()))
            print("stau=%g, Ym=%g, A/T*(s-alphaYm)=%g, s-alpha*Ym=%g, Yeffeq=%g, 1-(Yeffeq/(s-alphaYm))**m=%g\n"%(stau,Ym,prt1,(stau-PhaseProp.alpha*Ym),Yeffeq,prt2))
            print("part1=%g,part2=%g"%(prt1,prt2))
            sys.exit()
            
         
        
        del tmp1,prt1,prt2,f,Rp,Ym,Yeffeq
        return gv
    
    else:
        print("Invalid VP function type in InputFile line 44")
        sys.exit()

def Gtld(PhaseProp,dt):
    res=PhaseProp.gminf
    for i in range(0,(PhaseProp.ngmi)):
        res+= (PhaseProp.gmi[i])*(1-mt.exp(-dt/(PhaseProp.gi[i])))*((PhaseProp.gi[i])/dt)
    
    return res

def Ltld(PhaseProp,dt):
    res=PhaseProp.ldinf
    for i in range(0,(PhaseProp.nldi)):
        res+=(PhaseProp.ldi[i])*(1-mt.exp(-dt/(PhaseProp.li[i])))*((PhaseProp.li[i])/dt)
    return res

def DevPartViscousStress(PhaseProp, Statetn, dt, dEve):
    StrsGmi = np.zeros((3,3,PhaseProp.ngmi))
    for ii in range(0,3):
        for jj in range(0,3):
            for kk in range(0,PhaseProp.ngmi):
                StrsGmi[ii,jj,kk] = mt.exp(-dt/(PhaseProp.gi[kk]))*(Statetn.DevPartSveViscous[ii,jj,kk])\
                    +(PhaseProp.gmi[kk])*(1-mt.exp(-dt/(PhaseProp.gi[kk])))*((PhaseProp.gi[kk])/dt)*(dEve[ii,jj])

    return StrsGmi

def SphPartViscousStress(PhaseProp,Statetn,dt,trdEve):
    StrsLdi =np.zeros(PhaseProp.nldi)
    for kk in range(0,(PhaseProp.nldi)):
        StrsLdi[kk]=mt.exp(-dt/PhaseProp.li[kk])*(Statetn.SphPartSveViscous[kk])\
            +(PhaseProp.ldi[kk])*(1-mt.exp(-dt/(PhaseProp.li[kk])))*((PhaseProp.li[kk])/dt)*trdEve
    
    return StrsLdi

def SumDevPartsViscousStresses(PhaseProp,Statetn,dt):
    res=np.zeros((3,3))
    for ii in range(0,3):
        for jj in range(0,3):
            res[ii,jj]=0
            for kk in range(0,(PhaseProp.ngmi)):
                res[ii,jj]+=mt.exp(-dt/(PhaseProp.gi[kk]))*(Statetn.DevPartSveViscous[ii,jj,kk])
                
    return res

def SumSphPartsViscousStress(PhaseProp,Statetn,dt):
    res=0
    for kk in range(0,PhaseProp.nldi):
        res+=mt.exp(-dt/(PhaseProp.li[kk]))*(Statetn.SphPartSveViscous[kk])
    
    return res

def StressSve(PhaseProp,Statetn,dt,dEve,trdEve):
    gltd=Gtld(PhaseProp,dt)
    ltld=Ltld(PhaseProp,dt)
    sumDev=SumDevPartsViscousStresses(PhaseProp,Statetn,dt)
    sumSph=SumSphPartsViscousStress(PhaseProp,Statetn,dt)
    Sve= 2*((PhaseProp.gminf)*Statetn.Eve + gltd*dEve + sumDev \
        + ((1/2)*((PhaseProp.ldinf)*(np.trace(Statetn.Eve))+ltld*trdEve+sumSph))*np.identity(3))
    
    '''Sve= (PhaseProp.gminf)*Statetn.Eve + gltd*dEve + sumDev \
        + ((PhaseProp.ldinf)*(np.trace(Statetn.Eve))+ltld*trdEve+sumSph)*np.identity(3)'''
    
    del gltd, ltld, sumDev, sumSph
    return Sve

'''def StressSinst(PhaseProp,F):
    b=np.dot(F,F.transpose())
    bdev=bo.DevPartMadelStress(b)
    ldp= mt.sqrt(np.trace(b)/3)
    tmp1=ldp/mt.sqrt(PhaseProp.N)
    
    StrsInst=np.zeros((3,3))   
    StrsInst= ((1/3)*(PhaseProp.CR)*(1/tmp1)*(bo.invL(tmp1)))*bdev
    
    try:
        J=det(F)
    except:
        print("Proplems occured when computing det(F) in StressSinst function in VEVPOperations line %i"%(bo.CurrLineNumber()))
        print("F=\n",F,"J=%g"%J)
        sys.exit()
    
    invF=np.zeros((3,3))
    try:
        invF=inv(F)
    except:
        print("Proplems occured when computing inv(F) in StressSinst function in VEVPOperations line %i"%(bo.CurrLineNumber()))
        print("invF=\n",invF)
    
    S=np.zeros((3,3))
    S= J*(np.dot(np.dot(invF,StrsInst),invF.transpose()))
    
    del b,bdev,ldp,tmp1,StrsInst,J,invF
    return S'''

def StressSinst(PhaseProp,Cvp,F):
    bvp=Cvp.transpose()
    bvpdev=bo.DevPartMadelStress(bvp)
    ldp= mt.sqrt(np.trace(bvp)/3)
    tmp1=ldp/mt.sqrt(PhaseProp.N)
    
    StrsInst=np.zeros((3,3))   
    StrsInst= ((1/3)*(PhaseProp.CR)*(1/tmp1)*(bo.invL(tmp1)))*bvpdev
    
    try:
        J=det(F)
    except:
        print("Proplems occured when computing det(F) in StressSinst function in VEVPOperations line %i"%(bo.CurrLineNumber()))
        print("F=\n",F,"J=%g"%J)
        sys.exit()
    
    invF=np.zeros((3,3))
    try:
        invF=inv(F)
    except:
        print("Proplems occured when computing inv(F) in StressSinst function in VEVPOperations line %i"%(bo.CurrLineNumber()))
        print("invF=\n",invF)
    
    S=np.zeros((3,3))
    S= J*(np.dot(np.dot(invF,StrsInst),invF.transpose()))
    
    del bvp,bvpdev,ldp,tmp1,StrsInst,J,invF
    return S

def StrainStressMeasuresFromCvp(PhaseProp,Statetn,dt,F,Cvp):
    C=np.zeros((3,3))
    C=np.dot(F.transpose(),F)
    try:
        invCvp=inv(Cvp)
    except:
        print("Proplems occured when computing inv(Cv) in StrainStressMeasuresFromCvp in VEVPOperations line %i"%(bo.CurrLineNumber()))
        print("Cvp=\n",Cvp,"invCvp=\n",invCvp)
        sys.exit()
    
    
    Eve=(1/2)*(np.dot(invCvp,np.dot(C,invCvp))-np.identity(3))
    dEve=Eve-(Statetn.Eve)
    trdEve=np.trace(dEve)
    
    
    Sve=StressSve(PhaseProp,Statetn,dt,dEve,trdEve)
    Sinst=StressSinst(PhaseProp,Cvp,F)
    
    S = np.dot(invCvp,np.dot(Sve,invCvp)) + Sinst
    
    Y=np.dot(C,S)
    Yeff= Y - np.dot(C,Sinst)
    
    del C,invCvp,dEve,trdEve,Sinst
    return Eve,Sve,S,Y,Yeff

def ComputeSoftFromPlasticity(PhaseProp,Statetn,ptau):
    dp = ptau-(Statetn.p)
    tmp= 1/(1+(((bo.computeH1(PhaseProp,ptau))/(PhaseProp.s1))+((bo.computeH2(PhaseProp,ptau))/(PhaseProp.s2)))*dp)
    stau=(((bo.computeH1(PhaseProp,ptau)+bo.computeH2(PhaseProp,ptau))*dp)+ Statetn.s)*tmp
    
    del dp, tmp
    return stau

def LocalNR(SolIteri,PhaseProp,Load,Statetn,dt,F):
    invUvp=np.zeros((3,3))
    invUvp[0,0]=SolIteri[0];invUvp[0,1]=SolIteri[1];invUvp[0,2]=SolIteri[2]
    invUvp[1,0]=SolIteri[1];invUvp[1,1]=SolIteri[3];invUvp[1,2]=SolIteri[4]
    invUvp[2,0]=SolIteri[2];invUvp[2,1]=SolIteri[4];invUvp[2,2]=SolIteri[5]
    #ptau=SolIteri[6]
    pMultitau=SolIteri[6]
    DpMulti= pMultitau-Statetn.pMulti
    
    Uvp=np.zeros((3,3))
    try:
        Uvp=inv(invUvp)
    except:
        print("Proplems occured when computing inv(invUvp) in LoaclNR function in VEVPOperations line %i"%(bo.CurrLineNumber()))
        print("invUvp=\n",invUvp,"Uvp=\n",Uvp)
        sys.exit()
    
    Cvp=np.dot(Uvp,Uvp)
    
    Eve,Sve,S,Y,Yeff = StrainStressMeasuresFromCvp(PhaseProp,Statetn,dt,F,Cvp)
    g,DgDYT=VPPotential(PhaseProp,Y,Yeff)
    
    DgDyTCvp = np.dot(DgDYT.transpose(),Cvp)
    ptau=(PhaseProp.kappa)*DpMulti*(mt.sqrt(np.trace(np.dot(DgDyTCvp,DgDyTCvp.transpose())))) + Statetn.p
    
    symmetricf=bo.SymmetricPart(np.dot(Cvp,DgDYT))
    dp=ptau-(Statetn.p)
    Tt=DpMulti*(np.dot(invUvp,np.dot(symmetricf,invUvp)))
    
    expTt=np.zeros((3,3))
    try:
        expTt=expm(Tt)
    except:
        print("Proplems occured when computing exp(Tt) in LoaclNR function in VEVPOperations line %i" %(bo.CurrLineNumber()))
        print("Tt=\n",Tt,"expTt= \n",expTt)
        sys.exit()
    
    invUvpexpTtinvUvp=np.dot(invUvp,np.dot(expTt,invUvp))
    
    invCvpn=np.zeros((3,3))
    try:
        invCvpn = inv(Statetn.Cvp)   
    except:
        print("Proplems occured when computing inv(Cvpn) in LoaclNR function in VEVPOperations line %i"%(bo.CurrLineNumber()))
        print("Cvpn=\n",Statetn.Cvp,"invCvpn=\n",invCvpn)
        sys.exit()
    
    Fx=np.zeros(7,dtype=float)
    Fx[0]=invUvpexpTtinvUvp[0,0]-invCvpn[0,0]
    Fx[1]=invUvpexpTtinvUvp[0,1]-invCvpn[0,1]
    Fx[2]=invUvpexpTtinvUvp[0,2]-invCvpn[0,2]
    Fx[3]=invUvpexpTtinvUvp[1,1]-invCvpn[1,1]
    Fx[4]=invUvpexpTtinvUvp[1,2]-invCvpn[1,2]
    Fx[5]=invUvpexpTtinvUvp[2,2]-invCvpn[2,2]
    stau=ComputeSoftFromPlasticity(PhaseProp,Statetn,ptau)
    Fx[6]= (dt*VPFunction(PhaseProp,Load,Y,Yeff,ptau,stau)) - dp 
    
    del invUvp,Uvp,ptau,Cvp,Eve,Y,Yeff,S,Sve,g,DgDYT,symmetricf,dp,Tt,expTt,invUvpexpTtinvUvp,invCvpn,stau,pMultitau,DpMulti,DgDyTCvp
    return Fx

def ExtractFromSol(PhaseProp,Statetn,dt,F,Sol):
    threshold=1.e-12
    invUvp=np.zeros((3,3))
    invUvp[0,0]=Sol[0];invUvp[0,1]=Sol[1];invUvp[0,2]=Sol[2]
    invUvp[1,0]=Sol[1];invUvp[1,1]=Sol[3];invUvp[1,2]=Sol[4]
    invUvp[2,0]=Sol[2];invUvp[2,1]=Sol[4];invUvp[2,2]=Sol[5]
    #ptau=Sol[6]
    pMultitau=Sol[6]
    DpMulti=pMultitau-Statetn.pMulti
    
    # Regularisation of Cvp
    for i in range(3):
        for j in range(3):
            if(abs(invUvp[i,j])<threshold):
                invUvp[i,j]=0
                
    try:
        Uvp=inv(invUvp)
    except:
        print("Proplems occured when computing inv(invUvp) in ExtractFromSol function in VEVPOperations line %i"%(bo.CurrLineNumber()))
        print("invUvp=\n",invUvp,"\nUvp=\n",Uvp)
        sys.exit()
        
    try:
        Cvp=fractional_matrix_power(Uvp,2)
    except:
        print("Proplems occured when computing Cvp=Uvp**2 in ExtractFromSol function in VEVPOperations line %i"%(bo.CurrLineNumber()))
        print("Uvp=\n",Uvp,"\nCvp=\n",Cvp)
        sys.exit()
    
    #
    Eve,Sve,S,Y,Yeff = StrainStressMeasuresFromCvp(PhaseProp,Statetn,dt,F,Cvp)
    g,DgDYT=VPPotential(PhaseProp,Y,Yeff)
    
    ptau=(PhaseProp.kappa)*DpMulti*(mt.sqrt(np.trace(np.dot(DgDYT.transpose(),Cvp)))) + Statetn.p
        
    stau=ComputeSoftFromPlasticity(PhaseProp,Statetn,ptau)
    
    del invUvp,Uvp,DpMulti,Eve,Sve,S,Y,Yeff
    return Cvp,pMultitau,ptau,stau

def DefGradientIncrement(Load,Time,StrainCond):
    dF=np.zeros((3,3))
    if ((Load.loading)=="SHEAR_12"):
        dstrech=(2*(StrainCond.strnf))/Time.nsteps
        dF[0,1]=dstrech
        
    elif((Load.loading)=="UNIAXIAL_1"):
        dstrech=(StrainCond.strnf)/Time.nsteps
        if((Load.testtype)=="COMPRESSION"):
            dF[0,0]=dstrech
        elif((Load.testtype)=="TENSILE"):
            dF[0,0]=dstrech
    
    return dF

def DeformationGradienttau(Statetn,Load,Time,StrainCond):
    Ftau = np.zeros((3,3))
    dF = DefGradientIncrement(Load,Time,StrainCond)
    #Ftau= Statetn.F + dF
    if(Load.loading=="UNIAXIAL_1"):
        if(Load.testtype=="COMPRESSION"):
            Ftau[0,0]= Statetn.F[0,0] + dF[0,0]
            Ftau[1,1]= 1/mt.sqrt(Ftau[0,0])
            Ftau[2,2]= 1/mt.sqrt(Ftau[0,0])
        if(Load.testtype=="TENSION"):
            Ftau[0,0]= Statetn.F[0,0] + dF[0,0]
            Ftau[1,1]= 1/mt.sqrt(Ftau[0,0])
            Ftau[2,2]= 1/mt.sqrt(Ftau[0,0])
    
    del dF
    return Ftau
    
def UpdateVariables(PhaseProp,Statetn,Statetau,dt,Ftau):
    dEve=np.zeros((3,3))
    dEve=(Statetau.Eve)-(Statetn.Eve)
    trdEve=np.trace(dEve)
    
    DevPartSveVisc=DevPartViscousStress(PhaseProp,Statetn,dt,dEve)
    SphPartSveVisc=SphPartViscousStress(PhaseProp,Statetn,dt,trdEve)
    
    b=np.dot(Ftau,Ftau.transpose())
    strnEigenVal, strnEigenVect = eig(b)
    strntau=mt.log(mt.sqrt(strnEigenVal[0]))*np.tensordot(strnEigenVect[:,0],strnEigenVect[:,0],axes=0) + mt.log(mt.sqrt(strnEigenVal[1]))*np.tensordot(strnEigenVect[:,1],strnEigenVect[:,1],axes=0)\
        +mt.log(mt.sqrt(strnEigenVal[2]))*np.tensordot(strnEigenVect[:,2],strnEigenVect[:,2],axes=0)
        
    return DevPartSveVisc,SphPartSveVisc,strntau
    
def GlobalResidual(Load,P):
    if(Load.loading=="UNIAXIAL_1"):
        resd=max(abs(P[1,1]/P[0,0]),abs(P[2,2]/P[0,0]))   
    return resd

def ExtractFromGlobSol(Fopt):
    #print(Fopt)
    F=np.zeros((3,3))
    
    F[0,0]=Fopt[0];F[0,1]=Fopt[1];F[0,2]=Fopt[2]
    F[1,0]=Fopt[3];F[1,1]=Fopt[4];F[1,2]=Fopt[5]
    F[2,0]=Fopt[6];F[2,1]=Fopt[7];F[2,2]=Fopt[8]

    return F
        
def F2DTo3D(F2D):
    F3D = np.zeros((3,3))
    F3D[0,0]=F2D[0]
    F3D[0,1]=F2D[1]
    F3D[0,2]=F2D[2]

    F3D[1,0]=F2D[3]
    F3D[1,1]=F2D[4]
    F3D[1,2]=F2D[5]

    F3D[2,0]=F2D[6]
    F3D[2,1]=F2D[7]
    F3D[2,2]=F2D[8]

    return F3D

def F3DTo2D(F3D):
    F2D = np.zeros(9)
    
    F2D[0]=F3D[0,0]; F2D[1]=F3D[0,1]; F2D[2]=F3D[0,2]
    F2D[3]=F3D[1,0]; F2D[4]=F3D[1,1]; F2D[5]=F3D[1,2]
    F2D[6]=F3D[2,0]; F2D[7]=F3D[2,1]; F2D[8]=F3D[2,2]
    
    return F2D
        