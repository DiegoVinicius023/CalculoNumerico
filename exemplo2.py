from sympy import symbols

def lagrange(vx,vy,npol):
    x = symbols('x')

    soma = 0
    for i in range(npol+1):
        produto = vy[i]
        for j in range(npol+1):
            if j != i:
                produto *= (x-vx[j])/(vx[i]-vx[j])
        soma += produto
    return soma

vx = [-1,0,2]
vy = [4,1,-1]

l = lagrange(vx,vy,2)
print(l.subs(symbols('x'),1))
