from typing import Callable
from Auxiliares.formulasErro import FormulasErro
import numpy as np

class PontoFlutuante:

    def newtonRapson(x: float, f: Callable[[float], float], df: Callable[[float], float], k: int, tol: float, funcerr: Callable[[float,float,float], bool] = FormulasErro.errRel):
        """
        :x0: ponto inicial
        :f: função a ser testada
        :df: derivada de f
        :k: número de interações
        :funcerr: função de avaliação do erro
        :tol: tolerancia do erro
        """
        if k==0: return x
        
        xant = x
        x -= f(x)/df(x)

        if funcerr(x,xant,tol):
            return x
        return PontoFlutuante.newtonRapson(x,f,df,k-1,tol,funcerr)
            
    def secante(x0: float, x1: float, f: Callable[[float], float], k: int, tol: float, funcerr: Callable[[float,float,float], bool] = FormulasErro.errRel):
        """
        :x0: primeiro ponto inicial
        :x1: segundo ponto inicial
        :f: função a ser testada
        :k: número de interações
        :funcerr: função de avaliação do erro
        :tol: tolerancia do erro
        """
    
        d = (f(x1)-f(x0))/(x1-x0)
        x = x1-(f(x1)/d)
        k -= 1
        if funcerr(x,x1,tol) or k == 0:
            return x
        return PontoFlutuante.secante(x1,x,f,k,tol,funcerr)
    
    def bissecao(x0: float, x1: float, f: Callable[[float], float], k: int, tol: float, funcerr: Callable[[float,float,float], bool] = FormulasErro.errAbs):
        """
        :x0: primeiro ponto inicial
        :x1: segundo ponto inicial
        :f: função a ser testada
        :k: número de interações
        :funcerr: função de avaliação do erro
        :tol: tolerancia do erro
        """
        x = (x0+x1)/2
        k -= 1
        
        if np.sign(f(x0))*np.sign(f(x)) < 0:
            if funcerr(f(x0),f(x),tol) or k == 0:
                return (x+x0)/2
            return PontoFlutuante.bissecao(x0,x,f,k,tol)
        
        if funcerr(f(x),f(x1),tol) or k == 0:
                return (x+x1)/2
        return PontoFlutuante.bissecao(x,x1,f,k,tol)