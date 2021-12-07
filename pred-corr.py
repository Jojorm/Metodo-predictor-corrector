#Adams-Bashforth-Moulton Cuarto-Order Metodo Predictor-Corrector

import math, sys 
import numpy as np
import matplotlib.pyplot as plt

def f(y, t):
    return -math.tan(t)*y+math.cos(t) #<-----------función a calcular

y0 = 0 #y(0) = 5

def RungeKutta4(f, y0, T, n):
    #inicializador: RungeKutta4
    """Resuelve y'= f (t, y), y(0) = y0, con n pasos hasta que t=T. """
    t = np.zeros(n+1)  # variable usada para arreglos
    y = np.zeros(n+1)  # variable usada para arreglos
    y[0] = y0
    t[0] = 0
    dt = T/float(n)

    for k in range(n):
        t[k+1] = t[k] + dt
        K1 = dt * f( y[k], t[k])
        K2 = dt * f(y[k] + 0.5*K1, t[k] + 0.5*dt)
        K3 = dt * f(y[k] + 0.5*K2, t[k] + 0.5*dt)
        K4 = dt * f(y[k] + K3, t[k] + dt)
        y[k+1] = y[k] + (1/6.0)*(K1 + 2*K2 + 2*K3 + K4)
    return y, t

def PredictorCorrector(f, y0, T, n):
    """Resuelve y '= f (t, y), y (0) = y0, con n pasos hasta que t = T. """
    t = np.zeros(n+2)  # variable usada para arreglos
    y = np.zeros(n+2)  # variable usada para arreglos
    y[0] = y0
    t[0] = 0
    dt = T/float(n)
    print("\nValor de h:",dt)
    y, t = RungeKutta4(f, y0, T, n)  #<---- Se utilizan los valores iniciales del metodo anterior
    
    f0 = f(y[0],t[0])
    f1 = f(y[1],t[1])
    f2 = f(y[2],t[2])
    f3 = f(y[3],t[3])
    f4 = f(y[4],t[4])
    
    print('Cálculos inciales - Runge kutta:')
    print('f0:',round(f0,5),'\nf1:',round(f1,5),'\nf2:',round(f2,5),'\nf3:',round(f3,5))
    print('f4 RK4:',round(f4,5))
    
    for k in range(n-1,0,-1):
        #Predictor: El metodo Adams-Bashforth de cuarto orden, un método explícito de cuatro pasos, se define como: 
        y[k+1] = y[k] + (dt/24) *(55*f3 - 59*f2 + 37*f1 - 9*f0)
        #print('\nPredictor y4:', round(y[k+1],5))
        f4_AB = f(y[k+1],t[k+1])
        #print('f4 AB4:',round(f4_AB,5))
        
        #Corrector: La técnica de Adams-Moulton de cuarto orden, un método implícito de tres pasos, se define como: 
        y[k+1] = y[k] + (dt/24) *(9*f(y[k+1],t[k+1]) + 19*f3 - 5*f2 + f1)
        # print('\nCorrector one y4:', round(y[k+1],6))
        
        y[k+1] = y[k] + (dt/24) *(9*f(y[k+1],t[k+1]) + 19*f3 - 5*f2 + f1)
        # print('Corrector two y4:', round(y[k+1],6))
    return y, t

#------------------------------------------------------------------------------------------------------------
#  t = 4 and n = 64 

Kutta4 = y, t = RungeKutta4(f, y0, 5, 30)  #tercer valor representa la x en este caso la tfinal y el cuarto el numero de iteraciones
BashMoulton = y1, t1 = PredictorCorrector(f, y0, 5, 30)  #tercer valor representa la x en este caso la tfinal y el cuarto el numero de iteraciones

print('\nRunge Kutta')
print(' t \t  y')
for i in range(30):
    print(round(t[i+1],1),'\t',round(y[i],8), end ='\n')
    
print('\n\nEl resultado final de y es:',y[i])

print('\nAdams Bashforth Moulton')
print(' t \t  y')
for i in range(30):
    print(round(t[i+1],1),'\t',round(y1[i],8), end ='\n')
    
print('\n\nEl resultado final de y es:',y1[i],'\n')

#-----------------------------Gráfica----------------------------------------
plt.plot(BashMoulton[1],BashMoulton[0])
plt.plot(Kutta4[1],Kutta4[0])

plt.legend(['Adams Bashforth Moulton','RungeKutta4'], loc = 2)
plt.xlabel('t')
plt.ylabel('y')
plt.show()
