from sympy import Matrix
from auxiliares import Auxiliares
import numpy as np

class Atividade3:
    def ex1():
        vt = [0.2,0.3,0.4,0.5,0.6,0.7,0.8]
        vI = [3.16,2.38,1.75,1.34,1.00,0.74,0.56]

        vlnI = [np.log(i) for i in vI]

        npol = 1
        narred = 4
        minquad = Auxiliares.minimos_quadrados(vt,vlnI,npol,narred)

        I0 = round(np.exp(float(minquad[0])),narred)
        alpha = round(-minquad[1],narred)
        return [I0,alpha]

    def ex2():
        vx = [-2,-1,0,1,2]
        vy = [6,3,-1,2,4]
        npol = 2
        narred = 4
        return Auxiliares.minimos_quadrados(vx,vy,npol,narred)
    
    def ex3():
        vf = [12.211,12.583,11.556,12.792,12.654,11.863,12.082,11.819,12.070,11.543]
        vc = [4.042,4.117,3.911,4.158,4.131,3.973,4.016,3.964,4.014,3.909]

        return Auxiliares.minimos_quadrados(vf,vc,1)
    
    def ex4():
        vx = [1,9,18,28,42,56,63,70,85,94]
        vpt100 = [100.4,102.8,107.2,111.1,116.7,122.0,125.4,128.8,132.6,137.0]
        vntc = [510,450,290,180,100.9,55,40,27.5,16.2,10.3]

        pt100 = Auxiliares.minimos_quadrados(vx,vpt100,1)
        alpha = round(pt100[1]/pt100[0],5)
        ntc = Auxiliares.minimos_quadrados(vx,vntc,1)
        beta = "não feito"
        return [alpha, beta]

    def ex5():
        vx = [-3,-1,1,2,3]
        vy = [-1,0,1,1,-1]
        npol = 2
        narred = 4
        return Auxiliares.minimos_quadrados(vx,vy,npol,narred)

print(f"1: {Atividade3.ex1()}")
print(f"2: {Atividade3.ex2()}")
print(f"3: {Atividade3.ex3()}")
print(f"4: {Atividade3.ex4()}")
print(f"5: {Atividade3.ex5()}")