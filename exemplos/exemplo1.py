import numpy as np
import sys,os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from Auxiliares.pontoFlutuante import PontoFlutuante

def g(x):
    #10-20(exp(-0.2*x)-exp(-0.75*x))
    #queremos f(x) = 5
    #então g(x) = 10-20(exp(-0.2*x)-exp(-0.75*x))-5 = 5-20(exp(-0.2*x)-exp(-0.75*x))
    return 5-20*(np.exp(-0.2*x)-np.exp(-0.75*x))

def dg(x):
    return 4*np.exp(-0.2*x)-15*np.exp(-0.75*x)

def erro(x,xant,tol):
    if np.abs((x-xant)/x) < tol: 
        return True
    return False

x0 = 1
tol = 10E-2
print(PontoFlutuante.newton_rapson(x0,g,dg,3,tol,erro))