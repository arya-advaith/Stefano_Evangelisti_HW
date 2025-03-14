import pandas as pd
import matplotlib.pyplot as plt

for i in [4,8,16,32,64,128,256,512,1024]:
    pd.read_csv('eigenvectors_'+str(i)+'e_1a1b1close.txt')
    plt.plot(pd.read_csv('eigenvectors_'+str(i)+'e_1a1b1close.txt'))
    plt.show()
