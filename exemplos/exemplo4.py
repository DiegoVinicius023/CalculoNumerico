import sys,os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from Auxiliares.diferenciacao import Diferenciacao


ylinha = lambda x,y: x-y+2

x0,x1 = 0,1
y0 = 2

h = 0.1
m = (x1-x0)/ h

print(Diferenciacao.euler(x0,y0,h,m,ylinha))