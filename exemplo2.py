from numpy import pi
from sympy import symbols, Matrix
from auxiliares import Auxiliares

def exlagrangeslides():
    vx = [-1,0,2]
    vy = [4,1,-1]

    l = Auxiliares.lagrange(vx,vy,2,1)
    return l[0]

def exlinear():
    vx = [0,1]
    vy = [1.35,2.94]
    xval = 0.73
    npol = 1
    narred = 2
    return Auxiliares.interpolacao(vx,vy,npol,narred,xval)

def exquad():
    pontos = [
        (0,0),
        (pi/6,0.328),
        (pi/4,0.560)
    ]
    vx,vy = Auxiliares.gera_vetores(pontos)

    npol = 2
    narred = 3  
    return Auxiliares.interpolacao(vx,vy,npol,narred)

def exlagrange():
    vx = [0,0.2,0.4,0.5]
    vy = [0,2.008,4.064,5.125]

    return round(Auxiliares.lagrange(vx,vy,len(vx)-1,0.3)[0],3)

def exdifdiv():
    vx = [0.3,1.5,2.1]
    vy = [3.09,17.25,25.41]
    dict_diferencas_div = {0: vy}
    
    return Auxiliares.gera_diferencas_divididas(vx,dict_diferencas_div)

def exinternewton():
    vx = [0,0.2,0.3,0.5,0.6]
    vy = [1.008,1.064,1.125,1.343,1.512]
    xval = 0.4

    return Auxiliares.internewton(vx,vy,xval)[0]

def exinternewton2():
    pontos = [
        (0,1),
        (0.1,2.001),
        (0.3,4.081),
        (0.6,8.296),
        (1,21)
    ]
    vx,vy = Auxiliares.gera_vetores(pontos)
    xval = 0.2

    return Auxiliares.internewton(vx,vy,xval)[0]

if __name__ == '__main__': 
    print(f"linear: {exlinear()}")
    print(f"quadrático: {exquad()}")
    print(f"lagrange: {exlagrange()}")
    print(f"Diferenças Divididas: {exdifdiv()}")
    print(f"Interpolação de Newton: {exinternewton()}")
    print(f"Interpolação de Newton 2: {exinternewton2()}")