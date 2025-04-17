from sympy import Matrix
from Auxiliares.solucaoSistema import SolucaoSistema

class MinimosQuadrados:

    def vetor_matriz(vx):
        l = []
        for x in vx:
            l.append([x])
        return Matrix(l)
    
    def gera_expandida(vx,vy,npol):
        y = MinimosQuadrados.vetor_matriz(vy)

        vgt = []
        for i in range(npol+1):
            vi = [x**i for x in vx]
            vgt.append(vi)
        Gt = Matrix(vgt)
        G = Gt.transpose()

        A = Gt*G
        B = Gt*y
        return A.col_insert(npol+1,B)
    
    def minimos_quadrados(vx,vy,npol,narred=5):
        expandida = MinimosQuadrados.gera_expandida(vx,vy,npol)
        return SolucaoSistema.eliminacao_gauss(expandida,npol,narred)
    
    def gera_matriz_coef(vx,vy,npol):
        vgt = []
        for i in range(npol+1):
            vi = [x**i for x in vx]
            vgt.append(vi)
        vgt.append(vy)
        return Matrix(vgt).transpose()