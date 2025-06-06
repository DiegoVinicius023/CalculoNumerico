import sys,os, math
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from Auxiliares.integracao import Integracao
# from Auxiliares.interpolacao import Interpolacao
import numpy as np

def extrapeziounico():
    x0 = 3
    x1 = 3.6
    # npontos = 6

    f = lambda x: 1/x

    area = Integracao.trapezios(x0,x1,f,6)

    dfn = lambda x: 2/x**3
    erro = Integracao.erroTrapezios(x0,x1,6,dfn)
    return area,erro

def extercodesimpson():
    df4 = lambda x: (24*((5*(x**4))-(10*(x**2))+1))/((1+(x**2))**5)
    x0,x1 = 0,1
    emax = 1E-4

    n = Integracao.nPontosTercoSimpson(x0,x1,df4,emax)
    erro = Integracao.erroTercoSimpson(x0,x1,n,df4)
    
    f = lambda x: 1/(1+x**2)
    valor = Integracao.tercoSimpson(x0,x1,n,f)
    return erro, round(4*valor,6)

def ex3oitavossimpson():
    x0 = 1
    x1 = 4
    nsubs = 3
    func = lambda x: np.log(x**3+(np.exp(x)+1)**(1/2))
    
    teste3 = Integracao.tresOitavosSimpson(x0,x1,nsubs,func)

    nsubs = 9
    teste9 = Integracao.tresOitavosSimpson(x0,x1,nsubs,func)
    return teste3,teste9

def exNPontosTresOitavosSimpson():
    e = 1E-4
    x0,x1 = 0, np.pi/2

    df4 = lambda x: math.sin(x)

    n = Integracao.nPontosTresOitavosSimpson(x0,x1,df4,e)

    h = (x1-x0)/n

    A = Integracao.tresOitavosSimpson(x0,x1,n,df4)

    return f"N: {n}, H: {h}, Área: {A*2}"

# print(f"Trapezio unico: {extrapeziounico()}")
# print(f"Erro terço de simpson: {extercodesimpson()}")
# print(f"3/8 de simpson: {ex3oitavossimpson()}")
print(f"n pontos 3/8 de simpson: {exNPontosTresOitavosSimpson()}")