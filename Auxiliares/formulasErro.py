import numpy as np

class FormulasErro:

    def errrel(x,xant,tol):
        return np.abs((x-xant)/x) < tol
    def erreelnum(x,xant):
        return np.abs((x-xant)/x)   
    def errabs(x,xant,tol):
        return np.abs(x-xant) < tol