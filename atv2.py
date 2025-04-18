from Auxiliares.solucaoSistema import SolucaoSistema
from sympy import Matrix

class Atividade2:
    def ex1a():
        matriz = Matrix([
            [2,2,5,1],
            [4,1,-1,5],
            [-1,3,1,-4]
            
        ])
        if not SolucaoSistema.Sassenfeld(matriz):
            return "Matriz não converge"
        vx = [0,0,0]
        return SolucaoSistema.gauss_seidel(vx,matriz,5,1E-3)
    
    def ex1b():
        matriz = Matrix([
            [4,1,-1,1,-2],
            [1,4,-1,-1,-1],
            [-1,-1,-5,1,0],
            [1,-1,2,1,1]
        ])
        if not SolucaoSistema.Sassenfeld(matriz):
            return "Matriz não converge"
        vx = [0,0,0,0]
        return SolucaoSistema.gauss_seidel(vx,matriz,5,1E-3)
    
    def ex1c():
        matriz = Matrix([
            [0.112,0.46,0.24,0.7],
            [0.652,0.36,0.12,0.8],
            [0.147,0.21,0.65,0.9]
        ])
        if not SolucaoSistema.Sassenfeld(matriz):
            return "Matriz não converge"
        vx = [0,0,0]
        return SolucaoSistema.gauss_seidel(vx,matriz,5,1E-3)
    
    def ex2():
        matriz = Matrix([
            
            [4,-1,2,2],
            [-1,5,3,3],
            [4,1,6,1]
        ])
        return SolucaoSistema.Sassenfeld(matriz)
    
    # def ex3():
    #     a=0.5
    #     c1,c2,c3=1

    #     [1,-a,0,c1] => b1 = a
    #     [-a,1,-a,c2] => b2 = 2a
    #     [-a,0,1,c3] => b3 = a
    #     como b1,b2 e b3 são menores que 1, 2a < 1, a < 0.5

    def ex4():
        matriz = Matrix([
            [5,-2,5,2,7],
            [-6,4,-8,1,-9],
            [9,-6,19,1,23],
            [6,-4,-6,15,11]
        ])

        matriz_escada,pivos = matriz.rref()
        return matriz_escada

    def ex5():
        matrizA = Matrix([
            [1,-2,7,2],
            [2,5,-3,1],
            [9,-6,4,1],
            [4,-3,-6,7]
        ])
        matrizB = Matrix([
            [-18],
            [31],
            [35],
            [15]
        ])
        raise NotImplementedError()

if __name__ == '__main__':
    print(f"1a: {Atividade2.ex1a()}")
    print(f"1b: {Atividade2.ex1b()}")
    print(f"1c: {Atividade2.ex1c()}")
    print(f"2: {Atividade2.ex2()}")
    print(f"4:{Atividade2.ex4()}")
    # print(f"5: {Atividade2.ex5()}")