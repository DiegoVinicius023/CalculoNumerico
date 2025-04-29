from Auxiliares.formulasErro import FormulasErro
import numpy as np

class SolucaoSistema:

    def gaussSeidel(vx, matrix, k, tol,narred = 5, funcerr = FormulasErro.errRelNum):
        vxnew = vx.copy()
        exnew = []
        for i in range(len(vx)):
            vxnew[i] = round(SolucaoSistema.funcpad(vxnew,matrix.row(i),i),narred)
            exnew.append(funcerr(vxnew[i],vx[i]))
        k -= 1
        if (k==0) or (max(exnew) <= tol): return vxnew

        return SolucaoSistema.gauss_seidel(vxnew,matrix,k,tol,narred,funcerr)        
    
    def funcPad (vx, sist, pos):
        soma = 0
        for i in range(len(sist)):
            if i != pos and i != len(sist)-1:
                soma += vx[i]*sist[i]
            elif i == len(sist)-1:
                soma += sist[i]

        return soma/sist[pos]
    
    def Sassenfeld(matrix):
        for i in range(matrix.rows):
            soma = 0
            divisor = np.abs(matrix[i,i])
            if divisor == 0:
                return False
            
            #usa o numero de linhas para pegar apenas a matriz 
            for j in range(matrix.rows):
                if i != j:
                    soma += np.abs(matrix[i,j])
            if soma/divisor >= 1:
                return False
        return True
    
    def iniciaGaussSeidel(vx, matrix, k, tol,narred = 5, funcerr = FormulasErro.errRelNum):
        if not SolucaoSistema.Sassenfeld(matrix):
            return "Matriz não converge"
        return SolucaoSistema.gauss_seidel(vx,matrix,k,tol, narred,funcerr)
    
    def resolveLinha(vx,matriz_linha,pos,narred=5):
        soma = 0
        for i in range(len(vx)+1):
            if i != pos:
                if i == len(vx):
                    soma += matriz_linha[i]
                else:
                    soma -= vx[i]*matriz_linha[i]
        return round(soma/matriz_linha[pos],narred)

    def eliminacaoGauss(matrix,npol,narred):
        escada = matrix.rref()[0]
        vfinal = [0]*(npol+1)

        for i in range(npol+1):
            vfinal[npol-i] = SolucaoSistema.resolveLinha(vfinal,escada.row(npol-i),npol-i,narred)
        return vfinal