import numpy as np
import math as mt
import csv

class   PhaseProp:
    def __init__(self,nldi,ngmi):
        self.E0=float(0)
        self.nu=float(0)
        self.strsy=float(0)
        self.ar=float(0)
        self.orient= np.zeros(3, dtype=float)
        
        self.nldi=int(0)
        self.ld0=float(0)
        self.ldinf=float(0)
        self.ldi=np.zeros(nldi, dtype=float)
        self.li=np.zeros(nldi,dtype=float)
        
        self.ngmi=int(0)
        self.gm0=float(0)
        self.gminf=float(0)
        self.gmi=np.zeros(ngmi, dtype=float)
        self.gi=np.zeros(ngmi, dtype=float)
        
        self.YieldType = int(1)
        self.keq =float(0)
        self.kh=float(0)
        self.q=float(0)
        
        self.HardType = int(1)
        self.HardMod1 = float(0)
        self.HardExp = float(0)
        self.HardMod2 = float(0)
        
        self.VpFuncType=int(1)
        self.pDot0=float(0)
        self.A=float(0)
        self.T=float(0)
        self.alpha=float(0)
        self.m=float(0)
        
        self.VpPotType=int(1)
        self.TanPhi=float(0)
        self.ksi = float(0)
        self.nup = float(0)
        self.kappa=float(0)
        
        self.SoftType=int(1)
        self.s0=float(0)
        self.s1=float(0)
        self.s2=float(0)
        self.h1=float(0)
        self.h2=float(0)
        self.ps=float(0)
        self.k=float(0)
        
        self.RehardType=int(0)
        self.CR = float(0)
        self.N=float(0)
    
    def printPhaseProp(self):
        print("## VE Properties")
        print("E0=%g, nu=%g, strsy=%g"%(self.E0,self.nu,self.strsy))
        print("ld0=%g, lbdinf=%g, nlbi=%i"%(self.ld0,self.ldinf,self.nldi))
        print("ldi:", self.ldi)
        print("li:",self.li)
        print("gm0=%g, gminf=%g, ngmmi=%i"%(self.gm0,self.gminf,self.ngmi))
        print("gmi:",self.gmi)
        print("gi:",self.gi)
        print("## VP Properties")
        print("YieldType=%i, keq=%g, ks=%g, q=%g"%(self.YieldType,self.keq,self.kh,self.q))
        print("VPPotType=%i, ksi=%g, nup=%g, tanPhi=%g, kappa=%g"%(self.VpPotType,self.ksi,self.nup,self.TanPhi,self.kappa))
        print("VPFuncType=%i, pdot0=%g, T=%g, A=%g, alpha=%g, m=%g"%(self.VpFuncType,self.pDot0,self.T,self.A,self.alpha,self.m))
        print("SoftType=%i, s0=%g, s1=%g, s2=%g, h1=%g, h2=%g, ps=%g, k=%g"%(self.SoftType,self.s0,self.s1,self.s2,self.h1,self.h2,self.ps,self.k))
        print("RehadType=%i, CR=%g, N=%g"%(self.RehardType,self.CR,self.N))
    
    def computeGm0(self):
        #return E0/(2*(1+nu))
        self.gm0 = (self.E0)/(2*(1+(self.nu)))
        
    def computeLd0(self):
        #return (E0*nu)/((1-2*nu)*(1+nu))
        self.ld0 = ((self.E0)*(self.nu))/((1-2*(self.nu))*(1+(self.nu)))
        
    def computeGminf(self):
        sum=0
        for i in self.gmi:
            sum+= i
        return self.gm0-sum
    
    def ComputeLdinf(self):
        sum=0
        for i in self.ldi:
            sum+= i
        return self.ld0-sum
    
    def computeTanPhi(self):
        self.TanPhi=-(9/2)*((1-2*self.nup)/(1+2*self.nup))

    def computekappa(self):
        self.kappa = 1/(mt.sqrt(1+2*mt.pow((self.nup),2)))
    
        
class Time:
    def __init__(self):
        self.tf=float(0)
        self.nsteps=float(0)
        
    def printTime(self):
        print("tf=%g, nsteps=%i"%(self.tf,self.nsteps))
    
    def getTimeStep(self):
        return((self.tf)/(self.nsteps))

        
class StrainCond:
    def __init__(self):
        self.StrainRate=float(0)
        self.strnf=float(0)
        self.strs0=float(0)
        self.strn0=float(0)
        
    def printStrainCond(self):
        print("StrainRate=%g, strnf=%g, strn0=%g, strs0=%g"%(self.StrainRate,self.strnf,self.strn0,self.strs0))

    
class Loading:
    def __init__(self):
        self.history="MONOTONIC"
        self.testtype="COMPRESSION"
        self.loading="UNIAXIAL_1"
        self.nCycles=int(1)
        self.triaxiality=float(1)
        
    def printLoading(self):
        print("History=%s, Test type=%s, Load=%s, triaxiality=%g, ncycles=%i"%(self.history,self.testtype,self.loading,self.triaxiality,self.nCycles))


class Homogenization:
    def __init__(self):
        self.activated="OFF"
        self.scheme="MORI-TANAKA"
        self.isotrop="ESH-ISO"
        self.ratiodt=float(1)
        self.StatMoment=int(1)
        
    def printHomo(self):
        print("Homo=%s, Scheme=%s, Isotrop=%s, ratiodt=%g, StatMoment=%i"%(self.activated,self.scheme,self.isotrop,self.ratiodt,self.StatMoment))


class State:
    def __init__(self,PhaseProp):
        self.F=np.identity(3)
        self.Cvp = np.identity(3)
        self.strn=np.zeros((3,3))
        self.Eve=np.zeros((3,3))
        
        self.pMulti=0
        self.p=0
        self.s=PhaseProp.s0
        
        self.strs=np.zeros((3,3))
        self.P=np.zeros((3,3))
        self.Sve=np.zeros((3,3))
        self.S=np.zeros((3,3))
        self.Y=np.zeros((3,3))
        self.Yeff=np.zeros((3,3))
        
        self.SphPartSveViscous=np.zeros((PhaseProp.nldi))
        self.DevPartSveViscous=np.zeros((3,3,PhaseProp.ngmi))
        
    def printState(self):
        print("\n F:\n",self.F)
        print("logarithmic strain:\n",self.strn)
        print("Eve:\n",self.Eve)
        print("Cvp:\n",self.Cvp)
        print("pMultiplier=%g, p=%g, s=%g \n"%(self.pMulti,self.p,self.s))
        print("Cauchy's stress:\n",self.strs)
        print("1st P-K stress:\n",self.P)
        print("Sve:\n",self.Sve)
        print("Y:\n",self.Y)
        print("Yeff:\n",self.Yeff)
        print("strsLdi:\n",self.SphPartSveViscous)
        print("strsGmi:\n", self.DevPartSveViscous)
    
    def printStateOnTxtFile(self,OutputFile,tn):
        Row=[tn,self.strn[0,0],self.strn[1,1],self.strn[2,2],self.strs[0,0],self.strs[1,1],self.strs[2,2],\
            self.strn[0,1],self.strn[0,2],self.strn[1,2],self.strs[0,1],self.strs[0,2],self.strs[1,2],\
            self.pMulti,self.p,self.s,\
            self.F[0,0],self.F[1,1],self.F[2,2],self.Yeff[0,0],self.Yeff[1,1],self.Yeff[2,2],\
            self.F[0,1],self.F[0,2],self.F[1,2],self.F[1,0],self.F[2,0],self.F[2,1],\
            self.Yeff[0,1],self.Yeff[0,2],self.Yeff[1,2],self.Yeff[1,0],self.Yeff[2,0],self.Yeff[2,1]]
        with open(OutputFile,mode='a',encoding="utf8") as Outfile:
            writer=csv.writer(Outfile)
            writer.writerow(Row)
            Outfile.close()
  
    def copyStatetIn(self,Statetn):
        Statetn.F=self.F
        Statetn.Cvp=self.Cvp
        Statetn.strn=self.strn
        Statetn.Eve=self.Eve
        
        Statetn.pMulti=self.pMulti
        Statetn.p=self.p
        Statetn.s=self.s
        
        Statetn.strs=self.strs
        Statetn.P=self.P
        Statetn.Sve=self.Sve
        Statetn.S=self.S
        Statetn.Y=self.Y
        Statetn.Yeff=self.Yeff

        Statetn.SphPartSveViscous = self.SphPartSveViscous
        Statetn.DevPartSveViscous= self.DevPartSveViscous