import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# 1. Carregar os dados do CSV exportado pelo C++
df = pd.read_csv('trajetoria_iiwa.csv')

# 2. Criar um vetor de tempo (assumindo h = 0.01s)
h = 0.01
tempo = np.arange(len(df)) * h

# 3. Configurar a figura com 3 gráficos empilhados
fig, axs = plt.subplots(3, 1, figsize=(10, 14), sharex=True)

# Paleta com 10 cores distintas nativa do Matplotlib
cores = plt.cm.tab10.colors 

# --- Gráfico 1: Posições Angulares (x1 até x7) ---
for i in range(1, 8):
    axs[0].plot(tempo, df[f'x{i}'], label=f'$\\theta_{i}$', linewidth=2, color=cores[i-1])

axs[0].set_title('Posições Angulares (KUKA iiwa 14)', fontsize=14, fontweight='bold')
axs[0].set_ylabel('Ângulo [rad]', fontsize=12)
axs[0].grid(True, linestyle=':', alpha=0.7)
axs[0].legend(loc='best', ncol=4) 

# --- Gráfico 2: Velocidades Angulares (x8 até x14) ---
for i in range(1, 8):
    axs[1].plot(tempo, df[f'x{i+7}'], label=f'$\\omega_{i}$', linewidth=2, color=cores[i-1])

axs[1].set_title('Velocidades Angulares', fontsize=14, fontweight='bold')
axs[1].set_ylabel('Velocidade [rad/s]', fontsize=12)
axs[1].grid(True, linestyle=':', alpha=0.7)
axs[1].legend(loc='best', ncol=4)

# --- Gráfico 3: Esforço de Controle (Torque u1 até u7) ---
for i in range(1, 8):
    axs[2].plot(tempo, df[f'u{i}'], label=f'Torque {i} ($u_{i}$)', linewidth=2, color=cores[i-1])

axs[2].set_title('Esforço de Controle', fontsize=14, fontweight='bold')
axs[2].set_xlabel('Tempo [s]', fontsize=12)
axs[2].set_ylabel('Torque [Nm]', fontsize=12)
axs[2].grid(True, linestyle=':', alpha=0.7)
axs[2].legend(loc='best', ncol=4)

plt.tight_layout()
plt.show()
