## 2023.02.18
## Keith Ando Ogawa - keith.ando@usp.br
## Vinicius Viana viniciusviana@usp.br

## MAP3122

import numpy as np
import matplotlib.pyplot as plt

############################################################

def incremento(sinalR, sinalI, dz, theta): 
    ## Avança na reta designada

    if theta == 0 or theta == np.pi: ## Caso em que há movimento apenas na horizontal
        p = sinalR*dz + 0j
    elif theta == np.pi/2 or theta == 3*np.pi/2: ## Caso em que há movimento apenas na vertical
        p = 0 + sinalI*dz*1j
    else: 
        p = sinalR*dz + abs(np.tan(theta))*sinalI*dz*1j ## Há movimento no plano complexo
    return p

############################################################

def percorreReta(theta, estagios):
    dz = 1e-4 ## Espaçamento entre dois pontos escolhidos na reta

    ## SinalR indica o sentido para o qual o avanço ocorre em Re(z)
    if theta == np.pi/2 or theta == 3*np.pi/2: 
        sinalR = 0
    else:
        sinalR = np.cos(theta)/np.abs(np.cos(theta))
    
    ## SinalI indica o sentido para o qual o avanço ocorre em Im(z)
    if theta == 0 or theta == np.pi:
        sinalI = 0
    else:
        sinalI = np.sin(theta)/np.abs(np.sin(theta))

    abs = absPsi(-1, estagios)
    z = -1 + 0j
    while abs < 1: ## Encontra o primeiro ponto para o qual |psi(z)| >= 1
        z = z + incremento(sinalR, sinalI, dz, theta)
        abs = absPsi(z, estagios)
    return z 

############################################################

def absPsi(z, estagios):

     if estagios == 1:
         abs = np.abs(1 + z)
     elif estagios == 2:
         abs = np.abs(1 + z + 1/2*z**2)
     elif estagios == 3:
         abs = np.abs(1 + z + 1/2*z**2 + 1/6*z**3)
     else: ## estagios == 4
         abs = np.abs(1 + z + 1/2*z**2 + 1/6*z**3 + 1/24*z**4)
     return abs

############################################################

def g(z, estagios): ## Função usada no método de Newton
     ## Analisando os casos |psi(z)|=1
     if estagios == 1:
         g = np.abs(1 + z) - 1 
     elif estagios == 2:
         g = np.abs(1 + z + 1/2*z**2) - 1
     elif estagios == 3:
         g = np.abs(1 + z + 1/2*z**2 + 1/6*z**3) - 1
     else: ## estagios == 4
         g = np.abs(1 + z + 1/2*z**2 + 1/6*z**3 + 1/24*z**4) - 1
     return g  ## Função definida

############################################################

def dg(z, estagios): ## Derivada da função usada no método de Newton
     if estagios == 1:
         dg = (1 + z)/np.abs(1 + z)
     elif estagios == 2:
         dg = (2 + 2*z + z**2)*(z + 1)/np.abs(2 + 2*z + z**2)
     elif estagios == 3:
         dg = (z**2 + 2*z + 2)*(6 + 6*z + 3*z**2 + z**3)/(2*np.abs(6 + 6*z + 3*z**2 + z**3)) 
     else: ## estagios == 4
         dg = (6 + 6*z + 3*z**2 + z**3)*(24 + 24*z + 12*z**2 + 4*z**3 + z**4)/(6*np.abs(24 + 24*z + 12*z**2 + 4*z**3 + z**4))
     return dg ## Derivada da  função

############################################################

def newton (theta, estagios):
    delta = 1e-10 ## Delta pequeno que define quando o algoritmo interrompe

    if theta == 2*np.pi: ## intervalo: [0, 2pi]; caso de theta = 2pi é equivalente ao caso de theta = 0
        ang = 0 
    else:
        ang = theta

    chute = percorreReta(ang, estagios)

    for n in range(10):
        prox_chute = chute - g(chute, estagios)/dg(chute, estagios)
        if abs(prox_chute - chute) < delta or g(chute, estagios) > 1e-4:
            break
        chute = prox_chute

    return prox_chute

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
## Obtenção dos valores para a obtenção das regiões de estabilidade absoluta

total=200
t = np.linspace(0,2*np.pi,total) ## Vetor de 40 valores de x de 0 a 2pi
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


plt.xlim(-3, 3)
plt.ylim(-3, 3)

plt.show()