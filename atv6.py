from Auxiliares.diferenciacao import Diferenciacao
import numpy as np
import math

class Atividade4:
    def ex1():
        ylinha = lambda x,y: y*(x**2) - 3.0*y
        
        x0,x1 = 0,0.8
        h = 0.2

        m = math.ceil((x1-x0)/h)
        y0 = 1

        return round(Diferenciacao.rungeKutta4Ordem(x0,y0,h,m,ylinha),3)
    
    def ex2():

        ylinha = lambda x,y: -3*y + 7*np.exp(4*x)

        x0,x1 = 0,1
        y0 = 2

        h = 0.05

        m = math.ceil((x1-x0)/h)

        return round(Diferenciacao.ralston(x0,y0,h,m,ylinha),3)
    
    def ex3():

        ylinha = lambda m,y: (2E-6)*(m-y)*y 
        m = 100000
        y0 = 1000

        t0,t1 = 0,30

        h = (t1-t0)/m

        res = Diferenciacao.rungeKutta4OrdemIterativo(t0,y0,h,m,ylinha)
        return round(res,1)
    
    def ex4():

        ylinha = lambda x,y: y*(x**2)-3.0*y
        y0 = 0.5

        x0,x1 = 0,1.2

        h = 0.2

        m = math.ceil((x1-x0)/h)

        return round(Diferenciacao.heun(x0,y0,h,m,ylinha),3)
    
    def ex5():
        y0 = 1.80 
        v0 = 12

        h = 0.2

        t0,t1 = 0,1.22

        m = math.ceil((t1-t0)/h)

        vlinha = lambda v: v - 9.8*0.2
        ylinha = lambda y,v: y + v*0.2

        t = t0
        y = y0
        v = v0
        while m > 0:
            t += h

            v = vlinha(v)
            y = ylinha(y,v)
            
            m -=1

        return round(y,3)

    def ex6():
        ylinha = lambda x,y: -1.2*y + 7* np.exp(-0.3*x)

        x0,x1 = 0, 2.5
        h = 0.5

        y0 = 3

        m = math.ceil((x1-x0)/h)

        return round(Diferenciacao.rungeKutta4Ordem(x0,y0,h,m,ylinha),3)
    
    def ex7():
        mlinha = lambda m: m-5E-2*(m**(2/3))

        m0,m1 = 2,1

        t = 0
        m = m0
        while m > m1:
            t+=1
            m = mlinha(m)

        return t

    def ex8():

        ylinha = lambda x,y: -0.09*math.sqrt(y)
        y0 = 10
        h = 4

        t0,t1 = 0,12

        m = math.ceil((t1-t0)/h)

        return Diferenciacao.heun(t0,y0,h,m,ylinha)

if __name__ == '__main__':
    print(f"1: {Atividade4.ex1()}")
    print(f"2: {Atividade4.ex2()}")
    print(f"3: {Atividade4.ex3()}")
    print(f"4: {Atividade4.ex4()}")
    print(f"5: {Atividade4.ex5()}")
    print(f"6: {Atividade4.ex6()}")
    print(f"7: {Atividade4.ex7()}")
    print(f"8: {Atividade4.ex8()}")