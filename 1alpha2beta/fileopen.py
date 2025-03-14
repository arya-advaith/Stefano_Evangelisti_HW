
import pandas as pd
import matplotlib.pyplot as plt


list = [4,8,16,32,64,128,256,512,1024]

m=0

fig, ax = plt.subplots(len(list),2,figsize=(20,35))

for i in list:
    x = pd.read_csv('eigenvalues_'+str(i)+'e_1a2b1close.txt',sep='\t',header=None)
    x = x.sort_values(by=0,ascending=True)
    ax[m,0].plot(range(1,len(x)+1),x)
    ax[m,0].scatter(range(1,len(x)+1),x,alpha=0.3,color='green')
    ax[m,0].set_xlabel(str(i)+' atoms with '+str(i)+' electrons')
    ax[m,0].set_ylabel('Energy Eigenvalue')
    ax[m,1].plot(range(1,int(len(x)/2)+1),x[0:int(len(x)/2)])
    ax[m,1].scatter(range(1,int(len(x)/2)+1),x[0:int(len(x)/2)],alpha=0.3,color='green')
    ax[m,1].set_xlabel(str(i)+' atoms with '+str(int(i/2))+' electrons')
    ax[m,1].set_ylabel('Energy Eigenvalue')
    m+=1

plt.savefig('eigenvalues_1a2b1close.png')
