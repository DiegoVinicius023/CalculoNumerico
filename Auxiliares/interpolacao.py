from sympy import symbols
from math import factorial
from Auxiliares.minimosQuadrados import MinimosQuadrados
from Auxiliares.solucaoSistema import SolucaoSistema

class Interpolacao:

    def interpolacao(vx,vy,npol,xval=None,narred=4):
        matriz = MinimosQuadrados.geraMatrizCoef(vx,vy,npol)
        gauss = SolucaoSistema.eliminacaoGauss(matriz,npol,narred)

        x = symbols('x')
        equacao = 0
        for i in range(len(gauss)):
            equacao += gauss[i]*(x**i)
        equacao = equacao.simplify()
        if xval is None:
            return equacao
        return round(equacao.subs(symbols('x'),xval),narred),equacao
    
    def lagrange(vx,vy,npol,xval=None,narred=4):
        x = symbols('x')

        soma = 0
        for i in range(npol):
            produto = vy[i]
            for j in range(npol):
                if j != i:
                    produto *= (x-vx[j])/(vx[i]-vx[j])
            soma += produto
        
        soma = soma.simplify()
        if xval is None:
            return soma
        return round(soma.subs(symbols('x'),xval),narred),soma
    
    def geraVetores(pontos):
        return zip(*pontos)
    
    def difDiv(vx,vy,pos,ordem):
        return (vy[pos+1]-vy[pos])/(vx[pos+ordem+1]-vx[pos])

    def geraDiferencasDivididas(vx,dict_diferencas_div,narred=4):
        ordem = max(dict_diferencas_div.keys())

        if len(dict_diferencas_div[ordem]) > 1:
            ordem_n = []
            for i in range(len(dict_diferencas_div[ordem])-1):
                valor_dif_div = Interpolacao.difDiv(vx,dict_diferencas_div[ordem],i,ordem)
                ordem_n.append(round(valor_dif_div,narred))
            dict_diferencas_div[ordem+1] = ordem_n
            return Interpolacao.geraDiferencasDivididas(vx,dict_diferencas_div)
        return dict_diferencas_div
    
    def interNewton(vx,vy,xval=None,narred=4):
        diferencas_divididas = {0:vy}

        diferencas_divididas = Interpolacao.geraDiferencasDivididas(vx,diferencas_divididas)

        x = symbols('x')
        soma = 0
        for i in range(max(diferencas_divididas.keys())+1):
            produto = diferencas_divididas[i][0]
            for j in range(i):
                produto *= (x-vx[j])
            soma += produto
        
        soma = soma.simplify()
        if xval:
            return round(soma.subs(symbols('x'),xval),narred), soma
        return soma
    
    def maxError(vx,dfn,npol,xval,narred=4):
        x = symbols('x')
        maj = 0
        produto = 1

        for xi in vx:
            produto *= (x-xi)
            dxn = dfn(xi)
            if abs(dxn) > maj:
                maj = abs(dxn)
        produto *=(maj/factorial(npol+1))
        emax = round(abs(produto.subs(symbols('x'),xval)),narred)
        return emax
    
    def maxErrorPontos(vx,vy,npol,xval,narred=10):
        difer_div = Interpolacao.geraDiferencasDivididas(vx[:npol+2],{0:vy[:npol+2]},narred)

        x = symbols('x')
        produto = difer_div[max(difer_div.keys())][0]
        for xi in vx[:npol+1]:
            produto *= (x-xi)
        return produto.subs(symbols('x'),xval)
    
    def majoranteN(vx,dfn):
        maj = 0

        for xi in vx:
            dxn = dfn(xi)
            if abs(dxn) > maj:
                maj = abs(dxn)
        return maj

    def h(emax,maj,npol,narred=6):
        return round(pow((emax*4*(npol+1))/maj,1/(npol+1)),narred)