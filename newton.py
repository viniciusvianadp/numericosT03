import numpy as np
import matplotlib.pyplot as plt


def f(x):
    return  x**2#Função definida

def df(x):
    return 2*x#Derivada da  função

xvalues = np.linspace(-4,4,100) # Vetor de 100 valores de x de -4 a 4
yvalues = f(xvalues) #f(x) literalmente dos 100 valores

plt.axhline(0)
plt.plot(xvalues,yvalues) #plot

delta = 1e-10 #pequeno delta que define quando o algoritmo para

rvalues = np.linspace(-5,5,100) #Vetor de valores igualmente espaçacos de -5 a 5. Serão usados como valores iniciais
ivalues = np.linspace(-5,5,100) #Vetor igual acima mas para a parte complexa
#os valores -5,5 se baseam no intervalo que suspeita que está raiz, no exemplo foi -5 a 5.

roots = [] # vetor das raizes

#Feito de um modo que o vetor raizes nao tem valores duplicados

for r in rvalues:
    for i in ivalues:
        guess = r + i*1j
        for n in range(1000):#parada por iterações
            next_guess = guess -f(guess)/df(guess)
            if abs(next_guess - guess) < delta: # abs vai virar modulo
                already_in = False
                for root in roots:#checa se raiz n esta no vetor
                    if abs(next_guess - root) < delta: # se a diferença for menor que delta é a mesma raiz pras limitações do algoritmo
                        already_in = True
                        break

                if not already_in:#se a raiz não estiver
                    roots.append(next_guess)
                break
            guess = next_guess

print(roots)


















#Parte do Vídeo em que ele usa valores iniciais tendenciosos
#guess = 0 #Chute Inicial pode ser Complexo para gerar as raízes complexas
#while abs(f(guess)) > delta: #abs provavelmente vai virar um modulo vetorial
#    guess = guess - f(guess)/df(guess) #literalmente newton
#print(guess)
#
#for n in range(1000): #Parada por iterações 
#    next_guess = guess - f(guess)/df(guess)
#    if abs(next_guess - guess) < delta:
#        print(next_guess)
#        break
#    guess = next_guess 