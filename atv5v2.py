from sympy import Matrix, symbols
from Auxiliares.integracao import Integracao
import numpy as np

class Atividade5V2:

    def ex1():
        emax = 1E-4
        x0,x1 = 0,1

        df2 = lambda x: (4*(x**2)-2)*np.exp(-x**2)
        df4 = lambda x: (16*(x**4)-48*(x**2)+12)*np.exp(-x**2)

        n_trapezio = Integracao.n_pontos_trapezio(x0,x1,df2,emax)

        n_terco_simpson = Integracao.n_pontos_terco_simpson(x0,x1,df4,emax)

        return n_trapezio, n_terco_simpson
    
    def ex5():

        f = lambda x: -(x**4)+(4*x**2)
        df4 = lambda x: -24

        # g = lambda x: -f(x)
        # dg4 = lambda x: -24

        x0,x1 = -3,1

        emax = 5e-5

        n_interacoes = Integracao.n_pontos_terco_simpson(x0,x1,df4,emax)

        erro1 = Integracao.erro_terco_de_simpson(x0,x1,n_interacoes,df4)
        
        area = Integracao.terco_de_simpson(x0,x1,n_interacoes,f)
        return round(area,5)


if __name__ == '__main__':
    print(f"1: {Atividade5V2.ex1()}")
    print(f"5: {Atividade5V2.ex5()}")