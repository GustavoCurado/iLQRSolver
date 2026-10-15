import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# 1. Carregar os dados do CSV exportado pelo C++
df = pd.read_csv('trajetoria.csv')

# 2. Criar um vetor de tempo (sabemos que h = 0.05s pelo seu código)
h = 0.05
tempo = np.arange(len(df)) * h

# 3. Configurar a figura com 3 gráficos empilhados (Subplots)
fig, axs = plt.subplots(3, 1, figsize=(10, 12), sharex=True)

# --- Gráfico 1: Posições Angulares ---
axs[0].plot(tempo, df['theta1'], label=r'$\theta_1$ (Ombro)', linewidth=2, color='blue')
axs[0].plot(tempo, df['theta2'], label=r'$\theta_2$ (Cotovelo)', linewidth=2, color='orange')
# Linha de referência para a meta (pi/2)
axs[0].axhline(y=np.pi/2, color='red', linestyle='--', alpha=0.7, label='Objetivo $\\theta_1$ ($\pi/2$)')
axs[0].axhline(y=0.0, color='gray', linestyle='--', alpha=0.5)
axs[0].set_title('Posições Angulares (Swing-up do Acrobot)', fontsize=14, fontweight='bold')
axs[0].set_ylabel('Ângulo [rad]', fontsize=12)
axs[0].grid(True, linestyle=':', alpha=0.7)
axs[0].legend(loc='best')

# --- Gráfico 2: Velocidades Angulares ---
axs[1].plot(tempo, df['omega1'], label=r'$\omega_1$', linewidth=2, color='cyan')
axs[1].plot(tempo, df['omega2'], label=r'$\omega_2$', linewidth=2, color='magenta')
axs[1].set_title('Velocidades Angulares', fontsize=14, fontweight='bold')
axs[1].set_ylabel('Velocidade [rad/s]', fontsize=12)
axs[1].grid(True, linestyle=':', alpha=0.7)
axs[1].legend(loc='best')

# --- Gráfico 3: Esforço de Controle (Torque) ---
axs[2].plot(tempo, df['torque'], label='Torque Aplicado ($u$)', color='purple', linewidth=2)
axs[2].set_title('Esforço de Controle', fontsize=14, fontweight='bold')
axs[2].set_xlabel('Tempo [s]', fontsize=12)
axs[2].set_ylabel('Torque [Nm]', fontsize=12)
axs[2].grid(True, linestyle=':', alpha=0.7)
axs[2].legend(loc='best')

# Ajustar o layout para não cortar os textos e mostrar a janela
plt.tight_layout()
plt.show()
