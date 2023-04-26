import BasicOperations as bo
import Structures as st
import numpy as np
import VEVPOperations as vevp
import ConstitutiveVEVP as const
matrix, homogenization, loading, strainconditions, time =bo.ReadInputs("InputFile.txt")
'''Y=1.e9*np.ones((3,3))
Yeq=bo.EqMandelStress(Y)
print(Yeq)
print(bo.invL(3))
print(bo.computeH1(matrix,0.001))
print(bo.computeH2(matrix,0.001))
state1 = st.State(matrix)
state2 = st.State(matrix)

#state1.printState()
bo.copyState(state1,state2)

print(vevp.HardeningLaw(matrix,0.01))
print(vevp.Gtld(matrix,0.01))
print(vevp.DevPartViscousStress(matrix,state1,0.01,0.01*np.ones((3,3))))
F=0.8*np.identity(3);F[1,2]=0.2
print("DevVisctau:\n",vevp.DevPartViscousStress(matrix,state1,0.01,0.01*np.ones((3,3))))
print("sumDev:\n",vevp.SumDevPartsViscousStresses(matrix,state1,0.01))
print("sumSph:\n",vevp.SumSphPartsViscousStress(matrix,state1,0.01))
print("F:\n",F)
print("Sve:\n",vevp.StressSve(matrix,state1,0.01,0.01*np.ones((3,3)),0.03))
print("S:\n",vevp.StressSinst(matrix,F))
Y=1.e6*np.identity(3)
Y[0,1]=1.e6
C=0.01*np.identity(3)
C[0,1]=C[1,0]=0.01
Yeff=np.dot(C,Y)
g, DgDYT=vevp.VPPotential(matrix,Y,Y,0.01)
print("g,DgdYT:\n",g,DgDYT)
print("Y:\n",Y)
print("Ydev:\n",bo.DevPartMadelStress(Y))
print("Yeq=%g\n"%Yeq)
#vevp.VPFunction(matrix,loading,Y,Yeff,0.01,matrix.s0)
print("gv:\n",vevp.VPFunction(matrix,loading,Y,Yeff,0.001,matrix.s0))
print("Sinst:\n", vevp.StressSinst(matrix,F))
vevp.StrainStressMeasuresFromCvp(matrix,state1,0.01,F,np.identity(3))
print("s=%g\n"%vevp.ComputeSoftFromPlasticity(matrix,state1,0.001))
state1.s=matrix.s0
sol=np.zeros(7)
sol[0]=sol[3]=sol[5]=1
sol[6]=0
state1.printState()
vevp.LocalNR(sol,matrix,loading,state1,0.01,F)
state2, ok=const.constitutiveVEVP(matrix,loading,state1,0.01,F)
state2.printState()'''
'''def Test1():
    A = 1.e6*np.ones((3,3))
    A[2,2] =0; A[0,1]=5e6; A[2,0]=1e-3
    return A
print("PyMatrix:\n",Test1())

def Test2():
    A=1.e-3*np.ones(4)
    A[3]= 5; A[2]=2
    return A
print("PyVect:\n",Test2())

def Test3():
    A=1e6*np.ones((3,3,4))
    for i in range(3):
        for j in range(3):
            for k in range(4):
                A[i,j,k]= i*1e6 + k*2e6 + j*3e6
    return A

print("Py3DMatrix:\n",Test3())'''

Y=1.e6 * np.ones((3,3)); print("Y=\n",Y);
Y_dev = bo.DevPartMadelStress(Y); print("Y_dev=\n",Y_dev);
print("Y_eq=",bo.EqMandelStress(Y));
print("invL=",bo.invL(0.1));
print("H1(0.1)=",bo.computeH1(matrix,0.1));
print("H2(0.1)=",bo.computeH2(matrix,0.1));
