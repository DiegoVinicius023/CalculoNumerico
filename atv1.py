import numpy as np
from auxiliares import Auxiliares

class Atividade1:
    def ex1():
        # F(x)= x^2/2+x(ln(x)-1)
        # F(x)= x^2/2+xln(x)-x

        #queremos o ponto critico, então onde a derivada é 0
        #chamaremos a derivada de F de f
        #f(x) = F'(x) = x+ln(x)
        def f1(x):
            return x+np.log(x)

        def df1(x):
            return 1+1/x

        x0 = 2
        tol = 10E-4
        k = 15

        return Auxiliares.newton_rapson(x0,f1,df1,k,tol)
    
    def ex2():
        def f2(x):
            return (np.power(x,3)*np.exp(0.3*x)/6)+ (np.power(x,2)/2)+x+1-np.exp(x)
        
        def df2(x):
            return (np.power(x,2)*np.exp(0.3*x)/2)+(np.power(x,3)*np.exp(0.3*x)/20)+x+1-np.exp(x)
        
        x0=5
        tol = 1E-7
        k=40

        return Auxiliares.newton_rapson(x0,f2,df2,k,tol)
    
    def ex3():
        def f3(x):
            return 2*np.power(x,4)-2*np.power(x,3)-np.power(x,2)-x-3
        
        def df3(x):
            return 8*np.power(x,3)-6*np.power(x,2)-2*x-1
        
        x0= 1
        x1 = 1.5
        tol = 1E-5
        k = 10
        # return Auxiliares.newton_rapson(x0,df3,df3,k,tol)
        return Auxiliares.secante(x0,x1,df3,k,tol)
    
    def ex4():
        def f4(h,l=3,r=0.3,v=0.33):
            return l*(0.5*np.pi*np.power(r,2)-np.power(r,2)*np.arcsin(h/r)-h*np.power(np.power(r,2)-np.power(h,2),0.5))-v
        
        x0= 0.01
        x1= 0.05
        k= 20
        tol = 1E-3
        #faltou fazer que a profundidade é R-h
        # return Auxiliares.secante(x0,x1,f4,k,tol)
        h = Auxiliares.secante(x0,x1,f4,k,tol)
        r = 0.3
        return r-h
    
    def ex5():
        def f5(x):
            return 2*np.sin(x)-x
        x0 = -2
        x1 = 2
        k = 300
        tol = 10E-5

        return Auxiliares.bissecao(x0,x1,f5,k,tol)
    
    def ex6():
        #precisamos achar o 0 da derivada para ter o ponto minimo
        # def f6(x):
        #     return np.power(x,4)+np.power((x-1),2)
        
        def f6(x):
            return 4*np.power(x,3)+2*x-2
        
        def df6(x): return 12*np.power(x,2)+2

        tol = 10E-4
        x0 = 1
        k=20

        return Auxiliares.newton_rapson(x0,f6,df6,k,tol)
    
    def ex7():
        def f7(x):
            return np.power(x,4)-3*np.power(x,2)-3
        
        x0= 1
        x1 = 2
        #a notação 10E-3 dá 10*10^-3, o correto é 1E-3
        tol = 1E-3
        tol = np.power(10.0,-3)
        k=20

        return Auxiliares.secante(x0,x1,f7,k,tol)


print(f"Atv1: {Atividade1.ex1()}\n")
print(f"Atv2: {Atividade1.ex2()}\n")
print(f"Atv3: {Atividade1.ex3()}\n")
print(f"Atv4: {Atividade1.ex4()}\n")
print(f"Atv5: {Atividade1.ex5()}\n")
print(f"Atv6: {Atividade1.ex6()}\n")
print(f"Atv7: {Atividade1.ex7()}\n")