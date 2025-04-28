from sympy import Matrix
from Auxiliares.solucaoSistema import SolucaoSistema

class MinimosQuadrados:

    def vetorMatriz(vx):
        l = []
        for x in vx:
            l.append([x])
        return Matrix(l)
    
    def geraExpandida(vx,vy,npol):
        y = MinimosQuadrados.vetorMatriz(vy)

        vgt = []
        for i in range(npol+1):
            vi = [x**i for x in vx]
            vgt.append(vi)
        Gt = Matrix(vgt)
        G = Gt.transpose()

        A = Gt*G
        B = Gt*y
        return A.col_insert(npol+1,B)
    
    def minimosQuadrados(vx,vy,npol,narred=5):
        expandida = MinimosQuadrados.geraExpandida(vx,vy,npol)
        return SolucaoSistema.eliminacaoGauss(expandida,npol,narred)
    
    def geraMatrizCoef(vx,vy,npol):
        vgt = []
        for i in range(npol+1):
            vi = [x**i for x in vx]
            vgt.append(vi)
        vgt.append(vy)
        return Matrix(vgt).transpose()