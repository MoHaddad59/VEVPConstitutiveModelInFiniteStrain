import matplotlib.pyplot as plt
import csv
import numpy as np
import BasicOperations as bo

def DrawPlots(OutputFile,Test):
    if(Test.loading=="SHEAR_12"):
        file = open(OutputFile, mode='r',encoding='utf8')
        dic= csv.DictReader(file)
        
        time= []
        strn12= []
        strs12= []
        Yeff12=[]
        s=[]
        p=[]
        pmulti=[]
        
        for col in dic:
            time.append(float(col['t']))
            strn12.append(float(col['logstrn12']))
            strs12.append(float(col['Cauchystrs12']))
            Yeff12.append(float(col['Yeff12']))
            p.append(float(col['p']))
            s.append(float(col['s']))
            pmulti.append(float(col['pMultiplier']))
        
        fig,axs=plt.subplots(3,1)  
    
        axs[0].plot(strn12,strs12,linewidth = 2, label='$\sigma_{12}$')
        axs[0].plot(strn12,Yeff12,linewidth = 2, color="red", label='$Y^{eff}_{12}$')
        axs[0].set_xlabel('$\epsilon_{12}$[-]')
        axs[0].set_ylabel('$\sigma_{12}$ & $Y^{eff}_{12}$ [MPa]')
        axs[0].legend(loc=4, fontsize=10, numpoints=1, framealpha=0.0)
        
        axs[1].plot(time,pmulti,linewidth = 2, color='red',label='$\gamma$')
        axs[1].plot(time,p,linewidth = 2, label='p')
        axs[1].set_xlabel('t [s]')
        axs[1].set_ylabel('$\gamma$ & p [-]')
        axs[1].legend(loc=4, fontsize=8, numpoints=1, framealpha=0.0)
        
        axs[2].plot(time,s,linewidth = 2, label='Shear Test')
        axs[2].set_xlabel('t [s]')
        axs[2].set_ylabel('s [Pa]')
        
        for ax in axs.flat:
            ax.grid(color = 'grey', linestyle = '--', linewidth = 0.5)
            #ax.legend(loc=4, fontsize=14, numpoints=1, framealpha=0.0)

        #fig.savefig('sig(eps)_EP_uniaxial.png')
        plt.show()
    
    elif(Test.loading=="UNIAXIAL_1"):
        file = open(OutputFile, mode='r',encoding='utf8')
        dic= csv.DictReader(file)
        
        time= []
        F11=[]
        strn11= []
        strs11= []
        Yeff11=[]
        
        strn22= []
        strs22= []
        strn33= []
        strs33= []
        
        pmulti=[]
        s=[]
        p=[]
        
        for col in dic:
            time.append(float(col['t']))
            F11.append(float(np.real(col['F11'])))
            strn11.append(float(np.real(col['logstrn11'])))
            strs11.append(float(col['Cauchystrs11']))
            Yeff11.append(float(col['Yeff11']))
            
            strn22.append(float(col['logstrn22']))
            strs22.append(float(col['Cauchystrs22']))
            strn33.append(float(col['logstn33']))
            strs33.append(float(col['Cauchystrs33']))
            pmulti.append(float(col['pMultiplier']))
            p.append(float(col['p']))
            s.append(float(col['s']))
        
        '''if(Test.testtype=="COMPRESSION"):
            strn11= -1*np.array(strn11); strs11=-1*np.array(strs11)
            strn22= np.array(strn22); strs22=-1*np.array(strs22)
            strn33= np.array(strn33); strs33=-1*np.array(strs33)'''
        
        fig,axs=plt.subplots(3,1)  
    
        axs[0].plot(strn11,strs11,linewidth = 2, label='$\sigma_{11}$')
        axs[0].plot(strn11,Yeff11,linewidth = 2, color='red', label='$Y^{eff}_{11}$')
        #axs[0].plot(strn22,strs22,linewidth = 2, color='green',label='$\sigma_{22}$')
        #axs[0].plot(strn33,strs33,linewidth = 2, color='red',label='$\sigma_{33}$')
        
        axs[0].set_xlabel('$\epsilon_{11}$[-]')
        axs[0].set_ylabel('$\sigma_{11}$  [MPa]')
        axs[0].legend(loc=4, fontsize=10, numpoints=1, framealpha=0.0)
        
        axs[1].plot(time,pmulti,linewidth = 2, color='red', label='$\gamma$')
        axs[1].plot(time,p,linewidth = 2, label='p')
        axs[1].set_xlabel('t [s]')
        axs[1].set_ylabel('$\gamma$ & p [-]')
        axs[1].legend(loc=4, fontsize=10, numpoints=1, framealpha=0.0)
        
        
        axs[2].plot(time,s,linewidth = 2, label='Uniaxial tensile Test')
        axs[2].set_xlabel('t [s]')
        axs[2].set_ylabel('s [Pa]')
        
        for ax in axs.flat:
            ax.grid(color = 'grey', linestyle = '--', linewidth = 0.5)
            #ax.legend(loc=4, fontsize=14, numpoints=1, framealpha=0.0)

        #fig.savefig('sig(eps)_EP_uniaxial.png')
        plt.show()

matrix, homogenization, loading, strainconditions, time =bo.ReadInputs("InputFile.txt")
DrawPlots("results.csv",loading)