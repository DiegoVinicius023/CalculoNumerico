from Auxiliares.integracao import Integracao

class Diferenciacao:
    def euler(x, y , h, m, f):
        if m == 0: return y

        x_atual = x+h
        y_atual = y+h*f(x,y)

        return Diferenciacao.euler(x_atual,y_atual,h,m-1,f)
    
    def erroPassoEuler(x0,x1,h,df2):

        majorante = Integracao.encontraMaximo(df2,x0,x1)
        return (h**2)*majorante/2

    def rungeKutta4Ordem(x,y,h,m,f):
        if m == 0: return y
        k1 = f(x,y)
        k2 = f(x+(h/2),y+k1*(h/2))
        k3 = f(x+(h/2),y+k2*(h/2))
        k4 = f(x+h,y+ h*k3)

        y_atual = y+ (h/6)*(k1+(2*k2)+(2*k3)+k4)
        return Diferenciacao.rungeKutta4Ordem(x+h,y_atual,h,m-1,f)
    
    def rungeKutta4OrdemIterativo(x,y,h,m,f):

        while m > 0:
        
            k1 = f(x,y)
            k2 = f(x+(h/2),y+k1*(h/2))
            k3 = f(x+(h/2),y+k2*(h/2))
            k4 = f(x+h,y+ h*k3)

            y += (h/6)*(k1+(2*k2)+(2*k3)+k4)
            x += h
            m -= 1

        # return Diferenciacao.rungeKutta4Ordem(x+h,y_atual,h,m-1,f)
    
        return y
    
    def ralston(x,y,h,m,f):
        if m == 0: return y

        k1 = f(x,y)
        k2 = f(x+h*(3/4),y+k1*h*(3/4))

        y_atual = y + h*(k1*(1/3)+k2*(2/3))
        return Diferenciacao.ralston(x+h,y_atual,h,m-1,f)
    
    def heun(x,y,h,m,f):
        if m == 0: return y

        k1 = f(x,y)
        k2 = f(x+h,y+k1*h)

        y_atual = y + (h/2)*(k1+k2)
        return Diferenciacao.heun(x+h,y_atual,h,m-1,f)
