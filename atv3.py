from sympy import Matrix
from auxiliares import Auxiliares

class Atividade3:
    def ex1():
        pass

    def ex2():
        vx = [-2,-1,0,1,2]
        vy = [6,3,-1,2,4]
        npol = 2

        expandida = Auxiliares.gera_expandida(vx,vy,npol)
        vteste = [0]*(npol+1)
        res = Auxiliares.inicia_gauss_seidel(vteste,expandida,50,1E-3)
        
        pass


print(f"2: {Atividade3.ex2()}")