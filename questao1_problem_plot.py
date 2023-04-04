## 2023.02.08

# MAP3122

# plots for general implicit one-Step methods.

# problem with unkwnown exact solution 
#              (1) x'= 0.5x - 0.25xy               0<=t<=1.75
#              (2) y'= -0.75y + 0.25xy          x(0) = 1.5; y(0) = 2

import matplotlib.pyplot as plt
import numpy as np

#############################################################################

def phi(t1, y1, t2, y2, f):
    # define discretization function 
    return 0.5*(f(t1, y1)+f(t2, y2))     # euler 

############################################################################

def f(t, y):
   # bidimensional problem
    f0 =  0.5*y[0] - 0.25*y[0]*y[1]
    f1 =  -0.75*y[1] + 0.25*y[0]*y[1]
    
    return np.array([f0, f1])

############################################################################

def implicitMethod(T, n, yn, tn, f):

    dt = (T - tn[-1]) / n

    while tn[-1] < T:
        # initial guess
        ytil = yn[-1] + dt*f(tn[-1], yn[-1])
        diff = 1.0

        # fixed point iteration
        r = 0
        while r<20 and diff > 0.0001:
            ytil0 = ytil
            ytil = yn[-1] + dt*phi(tn[-1], yn[-1], tn[-1] + dt, ytil, f)
            diff = np.linalg.norm(ytil - ytil0)
            r = r + 1
        yn.append(ytil) # y(i+1) = ytil
        tn.append(tn[-1] + dt)
        dt = min(dt, T-tn[-1])
    yn = np.array(yn)

    return yn, tn

############################################################################

# other relevant data
t_n_1 = [0]; t_n_2 = [0]; t_n_3 = [0]; T = 1.75;        # time interval: t in [t0,T]
y_n_1 = [np.array([1.5, 2])]; y_n_2 = [np.array([1.5, 2])];
y_n_3 = [np.array([1.5, 2])]; # initial condition

n_1 = 8                # time interval partition (discretization)
y_n_1, t_n_1 = implicitMethod(T, n_1, y_n_1, t_n_1, f)

n_2 = 128                # time interval partition (discretization)
y_n_2, t_n_2 = implicitMethod(T, n_2, y_n_2, t_n_2, f)

n_3 = 256                # time interval partition (discretization)
y_n_3, t_n_3 = implicitMethod(T, n_3, y_n_3, t_n_3, f)

## plotting the graphic for x
plt.plot(t_n_1, y_n_1[:,0], 'k:', label = 'n = 8')
plt.plot(t_n_2, y_n_2[:,0], 'k--', label = 'n = 128')
plt.plot(t_n_3, y_n_3[:,0], 'k-', label = 'n = 256')


plt.xlabel('t   (em unidade de tempo)')
plt.ylabel('x(t)  (em unidade de x)')
plt.title('Aproximação Numérica da Variável de Estado x')
plt.legend()
plt.show()

## plotting the graphic for y
plt.plot(t_n_1, y_n_1[:,1], 'k:', label = 'n = 8')
plt.plot(t_n_2, y_n_2[:,1], 'k--', label = 'n = 128')
plt.plot(t_n_3, y_n_3[:,1], 'k-', label = 'n = 256')


plt.xlabel('t   (em unidade de tempo)')
plt.ylabel('y(t)  (em unidade de y)')
plt.title('Aproximação Numérica da Variável de Estado y')
plt.legend()
plt.show()
