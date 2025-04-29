import numpy as np

class FormulasErro:

    def errRel(x,xant,tol):
        return np.abs((x-xant)/x) < tol
    def errRelNum(x,xant):
        return np.abs((x-xant)/x)   
    def errAbs(x,xant,tol):
        return np.abs(x-xant) < tol