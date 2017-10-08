from cmath import exp, pi
from operator import mul
import sys

def fft(x,dir=1):
    N = len(x)
    if N <= 1: return x
    even = fft(x[0::2],dir)
    odd =  fft(x[1::2],dir)
    T= [exp((dir*-2j*pi*k)/N)*odd[k] for k in range(N//2)]
    return [even[k] + T[k] for k in range(N//2)] + [even[k] - T[k] for k in range(N//2)]

def ifft(x):
	a=fft(x,-1)
	tam_x=len(x)
	b=[nu/float(tam_x) for nu in a]
	return b

def parte_real(x):
	return [caca.real for caca in x]

def parte_real_redondeada(x):
	return [round(mierda) for mierda in parte_real(x)]


def determina_pot_2_minima(num):
	exp_2=0
	while num>2**exp_2:
		exp_2+=1
	return 2**exp_2

def normalizar_enteros(ent1,ent2):
	max_tam=max(len(ent1),len(ent2))
	max_tam_red=determina_pot_2_minima(max_tam)<<1
	return ent1+[0]*(max_tam_red-len(ent1)),ent2+[0]*(max_tam_red-len(ent2))

def multiplicar_enteros(ent1,ent2):
	ent1,ent2=normalizar_enteros(ent1,ent2)
#	print("ent1 redondeado {}".format(ent1))
#	print("ent2 redondeado {}".format(ent2))
	ent1_t=fft(ent1)
	#print("etn1 t {}".format(ent1_t))
	ent2_t=fft(ent2)
	entr_t=list(map(mul,ent1_t,ent2_t))
#	print("etnr t {}".format(entr_t))
	entr_tmp=parte_real_redondeada(ifft(entr_t))
	entr=entr_tmp[:]
#	print("resultado tmp {}".format(entr_tmp))
	for idx in range(len(entr)-1):
		coef=entr[idx]
		entr[idx]=coef%10
		entr[idx+1]+=coef//10
	return entr

def imprime_entero(ent):
	imprimir=False
	for coef in reversed(ent):
		if(coef):
			imprimir=True
		if(imprimir):
			print("{}".format(coef),end="")
	print("")

def caca():
	lineas_cnt=0
	lineas=[""]*2
	for linea in sys.stdin:
		lineas[lineas_cnt&1]=linea[:].strip()
		if(lineas_cnt and not(lineas_cnt&1)):
			num1=list(reversed(list(map(int,lineas[0]))))
			num2=list(reversed(list(map(int,lineas[1]))))
			numr=multiplicar_enteros(num1,num2)
			imprime_entero(numr)
#		print("linea act {}".format(linea))
		lineas_cnt+=1

	num1=list(reversed(list(map(int,lineas[0]))))
	num2=list(reversed(list(map(int,lineas[1]))))
	numr=multiplicar_enteros(num1,num2)
	imprime_entero(numr)

if __name__=="__main__":
#	a=[3,2,1]
#	b=[6,5,4]
#	c=multiplicar_enteros(a,b)
#	imprime_entero(c)
#	#print("tengo {}".format(c))
	caca()
