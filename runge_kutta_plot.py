import numpy as np
import matplotlib.pyplot as plt

# Define a função que calcula o módulo do fator de amplificação
def R2(z):
    return np.abs(1 + z + 1/2*z**2)

def R3(z):
    return np.abs(1 + z + 1/2*z**2 + 1/6*z**3)

def R4(z):
    return np.abs(1 + z + 1/2*z**2 + 1/6*z**3 + 1/24*z**4)
# Define a malha de pontos no plano complexo
x = np.linspace(-4, 4, 1000)
y = np.linspace(-4, 4, 1000)
xx, yy = np.meshgrid(x, y)
z = xx + 1j*yy

# Calcula o valor do módulo do fator de amplificação em cada ponto - Possívelmente Newton entra aqui
Rz2 = R2(z)
Rz3 = R3(z)
Rz4 = R4(z)

# Esboça o contorno de igualdade do módulo do fator de amplificação
plt.contourf(xx, yy, Rz2, levels=[0,1], colors=['white', 'gray'])
plt.contour(xx, yy, Rz2, levels=[1])
plt.contourf(xx, yy, Rz3, levels=[0,1], colors=['white', 'gray'])
plt.contour(xx, yy, Rz3, levels=[1])
plt.contourf(xx, yy, Rz4, levels=[0,1], colors=['white', 'gray'])
plt.contour(xx, yy, Rz4, levels=[1])
plt.xlabel('Parte Real')
plt.ylabel('Parte Imaginária')
plt.title('Região de estabilidade absoluta para RK2')
plt.show()
