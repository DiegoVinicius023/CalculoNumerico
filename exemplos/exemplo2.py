from numpy import pi
from sympy import symbols, Matrix
import sys, os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from Auxiliares.interpolacao import Interpolacao

def exlagrangeslides():
    vx = [-1,0,2]
    vy = [4,1,-1]

    l = Interpolacao.lagrange(vx,vy,2,1)
    return l[0]

def exlinear():
    vx = [0,1]
    vy = [1.35,2.94]
    xval = 0.73
    npol = 1
    narred = 2
    return Interpolacao.interpolacao(vx,vy,npol,xval,narred)

def exquad():
    pontos = [
        (0,0),
        (pi/6,0.328),
        (pi/4,0.560)
    ]
    vx,vy = Interpolacao.gera_vetores(pontos)

    npol = 2
    narred = 3  
    return Interpolacao.interpolacao(vx,vy,npol,None,narred)

def exlagrange():
    vx = [0,0.2,0.4,0.5]
    vy = [0,2.008,4.064,5.125]

    return round(Interpolacao.lagrange(vx,vy,len(vx)-1,0.3)[0],3)

def exdifdiv():
    vx = [0.3,1.5,2.1]
    vy = [3.09,17.25,25.41]
    dict_diferencas_div = {0: vy}
    
    return Interpolacao.gera_diferencas_divididas(vx,dict_diferencas_div)

def exinternewton():
    vx = [0,0.2,0.3,0.5,0.6]
    vy = [1.008,1.064,1.125,1.343,1.512]
    xval = 0.4

    return Interpolacao.internewton(vx,vy,xval)[0]

def exinternewton2():
    pontos = [
        (0,1),
        (0.1,2.001),
        (0.3,4.081),
        (0.6,8.296),
        (1,21)
    ]
    vx,vy = Interpolacao.gera_vetores(pontos)
    xval = 0.2

    return Interpolacao.internewton(vx,vy,xval)[0]

def exerrointerpol():
    vx = [2,2.5,4]
    vy = [0.5,0.4,0.25]
    npol = 2
    xval = 3

    resultado, interpolacao = Interpolacao.interpolacao(vx,vy,npol,3)

    dfn = lambda x: -6/(x**4)

    return Interpolacao.maxerror(vx,dfn,npol,xval)

def exh():
    vx = [0.1,0.2,0.3,0.4,0.5]
    emax = 5E-4
    npol = 1

    dfn = lambda x: -1/(x**2)
    maj = Interpolacao.majorante_n(vx[2:4],dfn)
    h = Interpolacao.h(emax,maj,npol)
    return h

if __name__ == '__main__': 
    print(f"linear: {exlinear()}")
    print(f"quadrático: {exquad()}")
    print(f"lagrange: {exlagrange()}")
    print(f"Diferenças Divididas: {exdifdiv()}")
    print(f"Interpolação de Newton: {exinternewton()}")
    print(f"Interpolação de Newton 2: {exinternewton2()}")
    print(f"Erro de interpolação <= {exerrointerpol()}")
    print(f"h para erro : {exh()}")