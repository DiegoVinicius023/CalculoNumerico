from sympy import Matrix, symbols
from auxiliares import Auxiliares
import numpy as np

class Atividade4:
    def ex1():
        vx = [0.1,0.2,0.3,0.4,0.5]
        vy = [-2.303,-1.609,-1.204,-0.916,-0.693]
        xval = 0.32

        vxlinear = vx[2:4]
        vylinear = vy[2:4]
        linear = Auxiliares.interpolacao(vxlinear,vylinear,1,xval)

        vxquad = vx[1:4]
        vyquad = vy[1:4]
        quad = Auxiliares.interpolacao(vxquad,vyquad,2,xval)

        emax = 5E-4
        npol = 1

        dfn = lambda x: -1/(x**2)
        maj = Auxiliares.majorante_n(vx[2:4],dfn)
        h = Auxiliares.h(emax,maj,npol)
        return f"\nLinear: {linear[0]}\nQuadrática: {quad[0]}\nh: {h}\n"
    
    def ex2():
        pontos = [
            (2.25,6.930),
            (2.5,8.726),
            (2.75,10.870)
        ]
        vx,vy = Auxiliares.gera_vetores(pontos)
        xval = 2.4

        res, pol = Auxiliares.internewton(vx,vy,xval)

        dfn = lambda x: ((6+x)/8)*np.exp(x/2)
        maxe = Auxiliares.maxerror(vx,dfn,len(vx)-1,xval)
        return f"\nPolinômio: {pol}\nPn(2.4): {res}\nmaxerror: {maxe}\n"
    
    def ex3():
        narred=5
        vx = [0.5,0.6,0.7,0.8,0.9,1.0]
        vy = [4.12706,3.48692,3.03787,2.70861,2.45959,2.26712]
        vyg = []
        npol = 2
        for i in range(len(vx)):
            vyg.append(vx[i]*vy[i])

        xval = 0.55

        valg = Auxiliares.internewton(vx[:npol+1],vyg[:npol+1],xval,narred)[0]
        pnxval = round(valg/xval,narred)

        maxe = Auxiliares.maxerrorpontos(vx,vyg,npol,xval,10)

        return pnxval, maxe
    
    def ex4():
        vx = [2.4,2.6,2.8,3.0,3.2,3.4,3.6,3.8]
        vy = [11.02,13.46,16.44,20.08,24.53,29.96,36.59,44.70]
        xval = 3.1

        lagrange = Auxiliares.lagrange(vx[3:6],vy[3:6],3,xval,4)
        return lagrange
    
    def ex5():
        vx = [1,3,4,5]
        vy = [0,6,24,60]
        xval = 3.5

        return Auxiliares.lagrange(vx,vy,4,xval) 
    
    def ex6():
        return "não resolvido"


if __name__ == '__main__':
    print(f"1: {Atividade4.ex1()}")
    print(f"2: {Atividade4.ex2()}")
    print(f"3: {Atividade4.ex3()}")
    print(f"4: {Atividade4.ex4()}")
    print(f"5: {Atividade4.ex5()}")
    print(f"6: {Atividade4.ex6()}")