from typing import Callable
import numpy as np
from sympy import Matrix

class Auxiliares:
    def errrel(x,xant,tol):
        return np.abs((x-xant)/x) < tol
    def erreelnum(x,xant):
        return np.abs((x-xant)/x)
    
    def errabs(x,xant,tol):
        return np.abs(x-xant) < tol
    
    def newton_rapson(x: float, f: Callable[[float], float], df: Callable[[float], float], k: int, tol: float, funcerr: Callable[[float,float,float], bool] = errrel):
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
        return Auxiliares.newton_rapson(x,f,df,k-1,tol,funcerr)
            
    def secante(x0: float, x1: float, f: Callable[[float], float], k: int, tol: float, funcerr: Callable[[float,float,float], bool] = errrel):
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
        return Auxiliares.secante(x1,x,f,k,tol,funcerr)
    
    def bissecao(x0: float, x1: float, f: Callable[[float], float], k: int, tol: float, funcerr: Callable[[float,float,float], bool] = errabs):
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
            return Auxiliares.bissecao(x0,x,f,k,tol)
        
        if funcerr(f(x),f(x1),tol) or k == 0:
                return (x+x1)/2
        return Auxiliares.bissecao(x,x1,f,k,tol)

    def gauss_seidel(vx, matrix, k, tol, funcerr = erreelnum):
        vxnew = vx.copy()
        exnew = []
        for i in range(len(vx)):
            vxnew[i] = Auxiliares.funcpad(vxnew,matrix.row(i),i)
            exnew.append(funcerr(vxnew[i],vx[i]))
        k -= 1
        if (k==0) or (max(exnew) <= tol): return vxnew

        return Auxiliares.gauss_seidel(vxnew,matrix,k,tol,funcerr)        
    
    def funcpad (vx, sist, pos):
        soma = 0
        for i in range(len(sist)):
            if i != pos and i != len(sist)-1:
                soma += vx[i]*sist[i]
            elif i == len(sist)-1:
                soma += sist[i]

        return soma/sist[pos]
    
    def Sassenfeld(matrix):
        for i in range(matrix.rows):
            soma = 0
            divisor = np.abs(matrix[i,i])
            if divisor == 0:
                return False
            
            #usa o numero de linhas para pegar apenas a matriz 
            for j in range(matrix.rows):
                if i != j:
                    soma += np.abs(matrix[i,j])
            if soma/divisor >= 1:
                return False
        return True