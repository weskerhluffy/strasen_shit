from math import log
import sys
class numerin():
	def __init__(self,num, pot_10):
		self.num=num
		self.pot_10=pot_10
	def __repr__(self):
		return "{}".format(self.num)

	def __mul__(self,other):
		return numerin(self.num*other.num,self.pot_10)
	__rmul__=__mul__

def determina_pot_min_10(numero):
	exp_10=0
	while numero>10**exp_10:
		exp_10+=1
	return 10**exp_10

def determina_pot_min_10_polinomio(polinomio):
	max_exp=len(polinomio)
	max_coef=max(polinomio)
	pot_10=determina_pot_min_10(max_exp*(max_coef**2))
	return pot_10

def empaca_polinomio(pol,pot_10):
	exp_pot_10=0
	num=0
	for coef in pol:
		num+=coef*(pot_10**exp_pot_10)
		exp_pot_10+=1
	return numerin(num,pot_10)

def desempaca_polinomio(nume):
	num=nume.num
	pot_10=nume.pot_10
	pol=[]
	while num:
		pol.append(num%pot_10)
		num//=pot_10
	return pol if pol else [0]*int(log(pot_10,10))

def multiplica_polinomios(pol1,pol2):
	pot_10=max(determina_pot_min_10_polinomio(pol1),determina_pot_min_10_polinomio(pol2))
	#print("pot 10 max es {}".format(pot_10))
	pol1=list(pol1)
	pol2=list(pol2)

	num1=empaca_polinomio(pol1,pot_10)
	num2=empaca_polinomio(pol2,pot_10)

	#print("pol 1 {} num 1 {}".format(list(pol1),num1))
	#print("pol 2 {} num 2 {}".format(list(pol2),num2))

	numr=num1*num2

	polr=desempaca_polinomio(numr)

	return polr


def completa_polinomio(pol, tam):
	return pol+[0]*(tam-len(pol))

def homogeiniza_polimonios(pol1,pol2):
	max_exp=max(len(pol1),len(pol2))
	return completa_polinomio(pol1,max_exp), completa_polinomio(pol2,max_exp)

def suma_polinomios(pol1,pol2):
	pol1,pol2=homogeiniza_polimonios(pol1,pol2)
	#print("sumando {} y {}".format(list(reversed(pol1)),list(reversed(pol2))))
	polr=list(map(lambda x,y:x+y,pol1,pol2))
	return polr

def resta_polinomios(pol1,pol2):
	pol1,pol2=homogeiniza_polimonios(pol1,pol2)
	polr=list(map(lambda x,y:x-y,pol1,pol2))
	return polr

def multiplica_polinomios_signados(pol1_p,pol1_n,pol2_p,pol2_n):
	polr_p=suma_polinomios(multiplica_polinomios(pol1_p,pol2_p),multiplica_polinomios(pol1_n,pol2_n))
	polr_n=suma_polinomios(multiplica_polinomios(pol1_p,pol2_n),multiplica_polinomios(pol1_n,pol2_p))
	#print("el pol p {} el n {}".format(list(reversed(polr_p)),list(reversed(polr_n))))
	polr=resta_polinomios(polr_p,polr_n)
	#print("podria ser q al fina {}".format(list(reversed(polr))))
	return polr

def separa_polinomio_por_signo(pol):
	pol_p=[0]*len(pol)
	pol_n=[0]*len(pol)
	for idx,coef in enumerate(pol):
		if(coef<0):
			pol_n[idx]=-coef
		else:
			pol_p[idx]=coef
	return pol_p,pol_n


def multiplica_polinomios_con_signo(pol1,pol2):
	pol1,pol2=homogeiniza_polimonios(pol1,pol2)

	pol1_p,pol1_n=separa_polinomio_por_signo(pol1)
	pol2_p,pol2_n=separa_polinomio_por_signo(pol2)

	polr=multiplica_polinomios_signados(pol1_p,pol1_n,pol2_p,pol2_n)

	return polr

def unados():
	lineas=list(sys.stdin)
	pol1=[int(x) for x in lineas[1].strip().split(" ")]
	pol2=[int(x) for x in lineas[2].strip().split(" ")]
	polr=multiplica_polinomios_con_signo(pol1,pol2)
	#print("{}".format(polr))
	for idx,coef in enumerate(polr):
		print("{}".format(coef),end="")
		if(idx):
			print("x^{}".format(idx),end="")
		if(idx<len(polr)-1):
			print(" + ",end="")

	print("")

if __name__=="__main__":
#	a=[41,49,38,29]
#	an=[0,0,0,0]
#	b=[19,23,46,21]
#	bn=[0,0,0,0]

#	a=[6,0,0,5]
#	an=[0,20,0,0]
#	b=[0,4,0,1]
#	bn=[0,0,2,0]
#	c=multiplica_polinomios_signados(a,an,b,bn)
#	a=[5,0,-20,6]
#	b=[1,-2,4]
#	c=multiplica_polinomios_con_signo(a,b)
#	#print("es un solo {}".format(list(reversed(c))))
	unados()
