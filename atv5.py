from sympy import Matrix, symbols
from Auxiliares.integracao import Integracao
import numpy as np

class Atividade5:
    def ex1():

        f = lambda x: (x**4)-(4*x**2)
        df4 = lambda x: 24

        # g = lambda x: -f(x)
        # dg4 = lambda x: -24

        x0,x1,x2 = -3,-2,1

        emax = 5e-5

        n1 = Integracao.n_pontos_terco_simpson(x0,x1,df4,emax)
        n2 = Integracao.n_pontos_terco_simpson(x1,x2,df4,emax)

        erro1 = Integracao.erro_terco_de_simpson(x0,x1,n1,df4)
        erro2 = Integracao.erro_terco_de_simpson(x1,x2,n2,df4)
        erro = erro1 + erro2
        
        a1 = Integracao.terco_de_simpson(x0,x1,n1,f)
        a2 = Integracao.terco_de_simpson(x1,x2,n2,f)
        return round(a1+a2,5)
    
    def ex2():
        x0 = 2
        x1 = 18 

        nsubs = 8

        def retorna(x):
            dic_pontos = {
                2:0.5,
                4:0.9,
                6:1.1,
                8:1.3,
                10:1.7,
                12:2.1,
                14:1.5,
                16:1.1,
                18:0.6
            }
            return dic_pontos[int(x)]
        
        primeira_parte = Integracao.tres_oitavos_de_simpson(x0,x1,nsubs,retorna)
        segunda_parte = Integracao.trapezios(14,18,retorna,2)
        return primeira_parte + segunda_parte
    
    def ex3():
        x0 = 0
        x1 = 1.2

        def retorna(x):
            dic_pontos = {
                0.0: (1.0,0.0),
                0.2: (0.8187,0.1987),
                0.4: (0.6703,0.3894),
                0.6: (0.5488,0.5646),
                0.8: (0.4493,0.7174),
                1.0: (0.3679,0.8415),
                1.2: (0.3012,0.9320),
            }
            tupla = dic_pontos[round(x,1)]
            return tupla[0]*tupla[1]
        
        n1,n2 = 3,6

        a1 = Integracao.tres_oitavos_de_simpson(x0,x1,n1,retorna)
        a2 = Integracao.tres_oitavos_de_simpson(x0,x1,n2,retorna)
        return f"{round(a1,4)}; {round(a2,4)}"
    
    def ex4():
        x0 = 0
        x1 = 120
        n = 6

        def retorna(x):
            dic_pontos = {
                0: 0,
                20: 22,
                40: 41,
                60: 53,
                80: 38,
                100: 17,
                120: 0
            }
            return dic_pontos[int(x)]

        return Integracao.tres_oitavos_de_simpson(x0,x1,6,retorna)


    def ex5():
        pass

    def ex6():
        df2 = lambda x: (4*(x**2)*np.exp(-x**2))*(-2*np.exp(-x**2))

        df4 = lambda x: (16*(x**4)-(48*x**2)+(12))*np.exp(-x**2)

        ntrapezio = Integracao.n_pontos_trapezio(0,1,df2,1E-4,2)
        ntercosimpson = Integracao.n_pontos_terco_simpson(0,1,df4,1E-4)
        return ntrapezio, ntercosimpson
    
    def ex7():

        def func(x):
            dic_pontos = {
                -2: -1,
                -1: 5,
                0: 1,
                1: 5,
                2: 35
            }
            return dic_pontos[int(x)]
        
        x0 = -2
        x1 = 2
        n = 4

        teste = Integracao.terco_de_simpson(x0,x1,n,func)
        # erro = Integracao.erro_terco_de_simpson(x0,x1,4,)
        return teste

    def ex8():
        x0 = 1
        x1 = 1.3
        n = 6

        def retorna(x):
            dic_pontos = {
                1.00: 1.0000,
                1.05: 1.0247,
                1.10: 1.0488,
                1.15: 1.0723,
                1.20: 1.0954,
                1.25: 1.1180,
                1.30: 1.1401
            }
            return dic_pontos[round(x,2)]

        return Integracao.trapezios(x0,x1,retorna,6)

    def ex9():
        x0=0
        x1=4
        n=5
        df4 = lambda x: 0
        return Integracao.erro_terco_de_simpson(x0,x1,n,df4)

    def ex10():
        f = lambda x: -x**4+(4*x**2)
        df4 = lambda x: -24

        n = Integracao.n_pontos_terco_simpson(-2,2,df4,1E-4)
        
        Area = Integracao.terco_de_simpson(-2,2,n,f)
        return Area

    def ex11():
        def retorna(x):
            dic_pontos = {
                1.0: 0.0,
                1.1: 0.095,
                1.2: 0.182,
                1.3: 0.262,
                1.4: 0.336,
                1.5: 0.405,
                1.6: 0.470,
            }
            return dic_pontos[round(x,1)]

        x0=1.0
        x1=1.6
        n=6

        return Integracao.terco_de_simpson(x0,x1,n,retorna)
    
    def ex12():
        x0 = 0
        x1 = 0.6

        f = lambda x: 1/(1+x)

        df2 = lambda x: 2 * (1+x)**(-3)

        df4 = lambda x: 24 * (1+x)**(-5)

        n = Integracao.n_pontos_trapezio(x0,x1,df2,5E-4)
        area_trapezio = Integracao.trapezios(x0,x1,f,n)

        n = Integracao.n_pontos_terco_simpson(x0,x1,df4,5E-4)
        area_simpson = Integracao.terco_de_simpson(x0,x1,n,f)
        return area_trapezio, area_simpson

if __name__ == '__main__':
    print(f"1: {Atividade5.ex1()}")
    print(f"2: {Atividade5.ex2()}")
    print(f"3: {Atividade5.ex3()}")
    print(f"4: {Atividade5.ex4()}")
    print(f"5: {Atividade5.ex5()}")
    print(f"6: {Atividade5.ex6()}")
    print(f"7: {Atividade5.ex7()}")
    print(f"8: {Atividade5.ex8()}")
    print(f"9: {Atividade5.ex9()}")
    print(f"10: {Atividade5.ex10()}")
    print(f"11: {Atividade5.ex11()}")
    print(f"12: {Atividade5.ex12()}")