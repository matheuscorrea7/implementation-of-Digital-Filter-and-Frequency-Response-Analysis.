from scipy.fft import rfft, rfftfreq, fft, fftfreq, irfft, ifft
import chardet
import scipy.constants as constants
from scipy import signal
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import math as m
import os
import shutil
from scipy.signal import find_peaks

M =151
N = 1024
string_M = ' ( M = ' + str(M) + ' )'
sm = 0.5 #Hz
S = int(N*sm)
T = float(1/sm)
n = np.arange(N) #Vetor de tamanho 1024 com elementos em sequência

################################################################################################	

def intesity(x, y, N, absolute):
	for i in range (int(N)):
		absolute.append(np.sqrt(x[i]*x[i] + y[i]*y[i]))
		
################################################################################################

def plotting(x, y, x_label, y_label, title_graphic, legend_graphic, name_graphic, start, end):
	font1 = {'family':'serif','color':'blue','size':20}
	font2 = {'family':'serif','color':'darkred','size':15}
	plt.figure(figsize=(12, 8))
	#plt.ticklabel_format(axis='both', style='sci', scilimits=(0,0))
	#plt.text(8, 8, 'ssss', style='italic', bbox={'facecolor': 'green', 'alpha': 0.5, 'pad': 10})
	plt.plot(x, y, color = 'b', label = legend_graphic, linewidth=0.5)
	plt.title(title_graphic + string_M, fontdict = font1)
	plt.xlabel(x_label)
	plt.ylabel(y_label)
	plt.xlim([start,end])
	plt.legend(loc='upper right')
	#if title_graphic == 'Senoide': plt.show()
	plt.savefig(name_graphic)
	plt.close()

################################################################################################

def displacement(Re, Im, N, M, PI,n):
	absolute = []
	T = np.zeros(N)
	for i in range(N):
		index = int(i + M/2)
		if index > N-1: index = index - N
		T[index] = Re[i]
		Re[i] = 0.0
		Im[i] = 0.0

	for i in range(N): Re[i]=T[i]
		
	for i in range(N):
		if i <= M:
			Re[i] = Re[i]*(0.54 - 0.46*(np.cos((2.0*PI*i)/M)))
		else:
			Re[i] = 0.0
			Im[i] = 0.0
			
	intesity(Re, Im, N, absolute)
	plotting(n, Re, 'Amostras', 'Amplitude', 'Filter Kernel', 'Filter-kernel', 'graphics_M'+ str(M) + '/Filter-kernel.png', 0, N)

################################################################################################

def reading_of_data(column1, column2, N):
	dataset = open('pssdata.txt', 'r')
	for line in dataset:
		line = line.split("\t")
		column1.append(float(line[0]))
		column2.append(float(line[1]))

	dataset.close()
	
	#Gráfico
	#plotting(x, y, x_label, y_label, title_graphic, legend_graphic, name_graphic)
	plotting(column1, column2, 'Frequência', 'Amplitude', 'Graphic-Initial', 'signal', 'graphics_M'+ str(M) + '/graphic-initial.png', 0, sm)
	
	#Completando com zeros de N/2 até N	
	for i in range(int(N/2),N,1):
		column1.extend(([0.0]))
		column2.extend(([0.0]))

################################################################################################

def create_folder(dir):
	existe = os.path.exists(dir)
	if existe == 'False':
		os.mkdir(dir) #Criando a pasta
	else:
		shutil.rmtree(dir, ignore_errors = True) #Apagando a pasta
		os.mkdir(dir) #Criando a pasta
		
################################################################################################

def peaks(y, sm):
	for i in range(0,int(N/2),1):
	  	 if y[i+1] > y[i] < y[i-1 or y[i + 1] < y[i] > y[i - 1]] or y[i + 1] == y[i] == y[i - 1]:
	   		y[i] = y[i] - sm
	
################################################################################################

def ajust(y, sm):
	for i in range(0,int(N/2 + 1),1):
	   		y[i] = 2.3*y[i] - float(sm)
	
################################################################################################

def convolution(y_sin, y_2):
	#eixo horizontal
	f_xx = []
	
	for i in range(0,N,1): f_xx.append(float(i/N))
	
	a = np.arange(2*N - 1) #Vetor de tamanho 1024 com elementos em sequência


	#Convolução
	h = signal.fftconvolve(y_sin, y_2, mode='full')
	print(len(h))

	#Plot da Saída
	plotting(a, h, 'Frequência', 'Amplitude', 'Convolution', 'h[i] = f[i]*g[i]', 'graphics_M'+ str(M) + '/convolucao.png', 0, N)
	
	return h
	
################################################################################################

def senoide(fr):
	fase = 0 #Fase para senoide
	fw=(float(fr/(M)))
	freqx = [] 
	y_sen = [] #Vetor para senoide
	x = [] #Vetor para quantidade de amostras 
	A = 1 #Amplitude	
	
	#Quantidade de amostras
	x = np.arange(0, 1024, 1) #não causa divisão por zero

	#Senoide
	y_sen = (A*np.sin(fr*x - fase))
	
	#Plot da saída
	plotting(x, y_sen, 'Amostras', 'Amplitude', 'Senoide', 'sen(freq*x)', 'graphics_M'+ str(M) + '/senoide.png', 0, 1024)

	return y_sen

################################################################################################

#Vetores
column1 = []
column2 = []
aux_x = []
f_y = []
f_x = []

#Variáveis comuns 
PI = np.pi

#Criando uma pasta onde estarão os plots
dir = 'graphics_M'+ str(M)
create_folder(dir)

#Chamada de função leitura
reading_of_data(column2, column1, N)

#Transformada inversa de fourrier do sinal
y = ifft(column1)
y1 = ifft(column2)
yt = y + y1

#Plot da #Tranformada inversa de fourrier
#Instruções: plotting(n, yt, 'Amostras', 'Amplitude', 'Impulse Response', 'Impulse-Response', 'graphics/Impulse-Response.png')
plotting(n, yt.real, 'Amostras', 'Amplitude', 'Impulse Response', 'Impulse-Response', 'graphics_M'+ str(M) + '/Impulse-Response.png', 0, N)

#Aplicando o deslocamento da função
for i in range(0,int(N),1): aux_x.append(float(i/N))

#Deslocamento
displacement(yt.real, aux_x, N, M, PI, n)

#Retornando para o domínio da frequência
yf = abs(rfft(yt.real))
for i in range(0,int(N/2 + 1),1): f_x.append(float(i/N))

if len(yf) == int(N/2 + 1):
	if len(f_x) == int(N/2 + 1):
		for i in range(0,int(N/2 + 1),1):
			f_y.append(np.sqrt(yf[i]*yf[i] + f_x[i]*f_x[i]))
	else:
		print('Parada Forçada')
		exit()

#Ajuste dos dados
ajust(f_y, sm)

#Plot da Saída
plotting(f_x, f_y, 'Frequência', 'Amplitude', 'Actual-Frequency Response', 'actual-frequency-response', 'graphics_M'+ str(M) + '/actual-frequency-response.png', 0, 0.5)

#Criando uma senoide
freq = 0.25 #pow(10,6)
y_sen = senoide(freq)

#Criando convolução entre o sinal obtido é uma senoide
y_conv = convolution(y_sen, yt.real)

print("Programa Finalizado")




















