from scipy.optimize import minimize_scalar
from math import ceil

class Integracao:
    def encontra_maximo(func, a, b):
        # Define the negative of the function to find the maximum
        neg_func = lambda x: -(func(x))
        
        # Use minimize_scalar to find the minimum of the negative function
        result = minimize_scalar(neg_func, bounds=(a, b), method='bounded')
        
        # The maximum value is the negative of the minimum of the negative function
        max_value = -result.fun
        return max_value
    
    def encontra_maximo_abs(func, a, b):
        # Define the negative of the function to find the maximum
        neg_func = lambda x: -abs(func(x))
        
        # Use minimize_scalar to find the minimum of the negative function
        result = minimize_scalar(neg_func, bounds=(a, b), method='bounded')
        
        # The maximum value is the negative of the minimum of the negative function
        max_value = -result.fun
        return max_value
    
    def trapezios(x0,x1,f,npontos):
        h = (x1-x0)/npontos
        
        yk = f(x0)
        A = 0
        for k in range (1,npontos+1):
            ykant = yk
            yk = f(x0+h*k)
            A+= (ykant+yk)
        return (h/2)*A
    
    def erro_trapezios(x0,x1,npontos,df2,eps = None):
        if not eps:
            eps = Integracao.encontra_maximo(df2,x0,x1)

        return -((x1-x0)**3/(12*npontos**2))*eps
    
    def n_pontos_trapezio(x0,x1,df2,emax,eps = None):
        if not eps:
            eps = Integracao.encontra_maximo_abs(df2,x0,x1)

        n = ((((x1-x0)**3)*eps)/(12*emax))**(1/2)

        return ceil(n)
    
    def terco_de_simpson(x0,x1,npontos,func):#,df4):
        h = (x1-x0)/npontos

        xkminus2 = x0
        ykminus2 = func(xkminus2)
        xkminus1 = x0+h
        ykminus1 = func(xkminus1)
        soma = 0
        # eps = 0
        for k in range(2,npontos+1,2):
            xk = x0+k*h
            yk = func(xk)

            soma += ykminus2+4*ykminus1+yk

            xkminus2 = xk
            if xkminus2 > x1:
                break
            ykminus2 = yk

            xkminus1 = xk+h
            if xkminus1 > x1:
                break
            ykminus1 = func(xkminus1)
            # eps = max(eps,Integracao.encontra_maximo(df4,xkminus2,xk))
        return (h/3)*soma

    def erro_terco_de_simpson(x0,x1,npontos,df4,eps=None):

        if not eps:
            eps = Integracao.encontra_maximo(df4,x0,x1)
        
        return (((x1-x0)**5)/(180*npontos**4))*eps
    
    def n_pontos_terco_simpson(x0,x1,df4,emax,eps = None):
        if not eps:
            eps = Integracao.encontra_maximo(df4,x0,x1)

        numerador = ((x1-x0)**5)*eps
        denominador = 180*emax

        ceiln = ceil(abs(numerador/denominador)**(1/4))
        if ceiln%2 != 0:
            return ceiln+1
        return ceiln        
    
    def tres_oitavos_de_simpson(x0,x1,nsubintervalos,func):#,df4):
        h = (x1-x0)/nsubintervalos

        xkminus3 = x0
        ykminus3 = func(xkminus3)
        
        xkminus2 = x0+h
        ykminus2 = func(xkminus2)

        xkminus1 = x0+2*h
        ykminus1 = func(xkminus1)

        soma = 0
        # eps = 0
        for k in range(3,nsubintervalos+1,3):
            xk = x0+k*h
            if xk > x1:
                break
            yk = func(xk)

            soma += ykminus3 + 3*ykminus2 + 3*ykminus1 + yk

            xkminus3 = xk
            if xkminus3 > x1:
                break
            ykminus3 = yk

            xkminus2 = xk+h
            if xkminus2 > x1:
                break
            ykminus2 = func(xkminus2)

            xkminus1 = xkminus2+h
            if xkminus1 > x1:
                break
            ykminus1 = func(xkminus1)
            # eps = max(eps,Integracao.encontra_maximo(df4,xkminus2,xk))
        return (3*h/8)*soma