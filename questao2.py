## 2023.02.18
## Keith Ando Ogawa - keith.ando@usp.br
## Vinicius Viana viniciusviana@usp.br

## MAP3122

## Obtains Runge-Kutta methods' absolute stability regions.
## The curves are obtained using and angle between 0 and 2pi, for which is created
## a straight line crossing (-1 + 0j); then, we deslocate throughout the line, starting
## at (-1 + 0j), until its found a point for which |psi(z)| >= 1 and this point will be
## the initial guess using Newton's method to obtain a value for a certain point in the curve.
## This process is repeated to methods with stages between 1 and 4.

import numpy as np
import matplotlib.pyplot as plt

############################################################

def increase(signalR, signalI, dz, theta): 
    ## Advances throughout the straight line

    if theta == 0 or theta == np.pi: ## Horizontal move only
        p = signalR*dz + 0j
    elif theta == np.pi/2 or theta == 3*np.pi/2: ## Vertical move only
        p = 0 + signalI*dz*1j
    else: 
        p = signalR*dz + abs(np.tan(theta))*signalI*dz*1j ## Moving in the complex plane
    return p

############################################################

def walkOnTheLine(theta, stages):
    dz = 1e-4 ## Spacing between to points analyzed

    ## signalR indicates the advance direction in Re(z) 
    if theta == np.pi/2 or theta == 3*np.pi/2: 
        signalR = 0
    else:
        signalR = np.cos(theta)/np.abs(np.cos(theta))
    
    ## signalI indicates the advance direction in Im(z)
    if theta == 0 or theta == np.pi:
        signalI = 0
    else:
        signalI = np.sin(theta)/np.abs(np.sin(theta))

    abs = absPsi(-1, stages)
    z = -1 + 0j
    while abs < 1: ## Finds a first point to which |psi(z)| >= 1
        z = z + increase(signalR, signalI, dz, theta)
        abs = absPsi(z, stages)
    return z 

############################################################

def absPsi(z, stages):

     if stages == 1:
         abs = np.abs(1 + z)
     elif stages == 2:
         abs = np.abs(1 + z + 1/2*z**2)
     elif stages == 3:
         abs = np.abs(1 + z + 1/2*z**2 + 1/6*z**3)
     else: ## stages == 4
         abs = np.abs(1 + z + 1/2*z**2 + 1/6*z**3 + 1/24*z**4)
     return abs

############################################################

def g(z, stages): ## Function used in Newton's Method
     ## Analyzing cases where |psi(z)|=1
     if stages == 1:
         g = np.abs(1 + z) - 1 
     elif stages == 2:
         g = np.abs(1 + z + 1/2*z**2) - 1
     elif stages == 3:
         g = np.abs(1 + z + 1/2*z**2 + 1/6*z**3) - 1
     else: ## stages == 4
         g = np.abs(1 + z + 1/2*z**2 + 1/6*z**3 + 1/24*z**4) - 1
     return g  ## Defined function

############################################################

def dg(z, stages): ## Derivative of the function used in Newton's Method
     if stages == 1:
         dg = (1 + z)/np.abs(1 + z)
     elif stages == 2:
         dg = (2 + 2*z + z**2)*(z + 1)/np.abs(2 + 2*z + z**2)
     elif stages == 3:
         dg = (z**2 + 2*z + 2)*(6 + 6*z + 3*z**2 + z**3)/(2*np.abs(6 + 6*z + 3*z**2 + z**3)) 
     else: ## stages == 4
         dg = (6 + 6*z + 3*z**2 + z**3)*(24 + 24*z + 12*z**2 + 4*z**3 + z**4)/(6*np.abs(24 + 24*z + 12*z**2 + 4*z**3 + z**4))
     return dg ## Derivative of the function

############################################################

def newton (theta, stages):
    delta = 1e-10 ## Small delta that defines when the algorithm should end

    if theta == 2*np.pi: ## theta is in [0, 2pi]; theta = 2pi is equivalent to theta = 0
        ang = 0 
    else:
        ang = theta

    guess = walkOnTheLine(ang, stages)

    for n in range(10):
        next_guess = guess - g(guess, stages)/dg(guess, stages)
        if abs(next_guess - guess) < delta or g(guess, stages) > 1e-4:
            break
        guess = next_guess

    return next_guess

############################################################

def RK1(theta):
  return newton(theta, 1)

def RK2(theta):
  return newton(theta, 2)

def RK3(theta):
  return newton(theta, 3)

def RK4(theta):
  return newton(theta, 4)

############################################################
## Obtaining values to generate absolute stability regions

total=200
t = np.linspace(0,2*np.pi,total) ## Array with 40 values from 0 to 2pi
y1 = [0] * total; b1 = [0] * total; a1 = [0] * total
y2 = [0] * total; b2 = [0] * total; a2 = [0] * total
y3 = [0] * total; b3 = [0] * total; a3 = [0] * total
y4 = [0] * total; b4 = [0] * total; a4 = [0] * total
for i in range(total):
    y1[i] = RK1(t[i]); a1[i] = y1[i].real; b1[i] = y1[i].imag
    y2[i] = RK2(t[i]); a2[i] = y2[i].real; b2[i] = y2[i].imag
    y3[i] = RK3(t[i]); a3[i] = y3[i].real; b3[i] = y3[i].imag
    y4[i] = RK4(t[i]); a4[i] = y4[i].real; b4[i] = y4[i].imag

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
ax.spines['left'].set_position('center')
ax.spines['bottom'].set_position('center')
ax.spines['right'].set_color('none')
ax.spines['top'].set_color('none')

plt.plot(a1, b1, 'k-')
plt.plot(a2, b2, 'k-')
plt.plot(a3, b3, 'k-')
plt.plot(a4, b4, 'k-')


plt.xlim(-5.5, 5.5)
plt.ylim(-4, 4)
plt.xticks([-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5])

plt.xlabel('Re(z)', loc='right')
plt.ylabel('Im(z)', loc='top')
plt.title('Regiões de estabilidade absoluta para RK')

plt.show()