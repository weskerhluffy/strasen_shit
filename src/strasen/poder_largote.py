#!/programas/python3.5/bin/python3


import logging
from math import pow
from sys import stdin
from operator import mul
from functools import reduce
from fractions import gcd

nivel_log = logging.ERROR
#nivel_log = logging.DEBUG
#logger_cagada = None

# XXX: https://stackoverflow.com/questions/18114138/computing-eulers-totient-function
def phi(n):
    amount = 0

    for k in range(1, n + 1):
        if gcd(n, k) == 1:
            amount += 1

    return amount

def numero_a_digitos(num,base):
    digitos = []
    #logger_cagada.debug("convirtiendo num {}".format(num))
    while num:
        digitos.append(num % base)
        num //= base
    #logger_cagada.debug("kedo en digitos {}".format(digitos))
    return digitos

def digitos_a_numero(digitos, base):
    factor = 1
    num = 0
    for digito in digitos:
        num += digito * factor
        factor *= base
    return num

def sacacaca(x):
    powers = []
    i = 1
    while i <= x:
        if i & x:
            powers.append(i)
        i <<= 1
    return powers

def mad_power(a,b,m=None):
        res=1
        assert a<m
        assert b>=0
#        pots_2=sacacaca(b)
#        #logger_cagada.debug("de num {} su descom {}".format(b,pots_2))
#        res=reduce(lambda x,y: mul(x%m,y%m)%m,(map(lambda x:pow(a,x)%m,pots_2)),1)%m
        pot=a
        b_tmp=b
        while b_tmp:
            if b_tmp&1:
                res=(res*pot)%m
#            #logger_cagada.debug("la pot {} el exp {} ele res {}".format(pot,b_tmp,res))
            pot=(pot*pot)%m
            b_tmp>>=1
        return res

def quita_sobrantes_coeficientes(coeficientes):
    while coeficientes and coeficientes[-1] == 0:
        coeficientes.pop()  # normalize

def reducir_a_modulo(digitos,modulo):
        digitos_tam=len(digitos)
        idx_act=digitos_tam
        num_final=0
        digitos_max=10
        while idx_act>0:
                idx_act=max(digitos_tam-digitos_max,0)
                #logger_cagada.debug("los digitos ant de red {}".format(digitos))
                porcion_num=digitos[idx_act:]
                #logger_cagada.debug("la porcion tomada {} a partir de {}".format(porcion_num, idx_act))
                num_conv=digitos_a_numero(porcion_num,10)
                #logger_cagada.debug("el num ia conv {}".format(num_conv))
                num_conv=num_conv%modulo
                #logger_cagada.debug("el num ia mod {}".format(num_conv))
                nueva_porcion=numero_a_digitos(num_conv,10)
                #logger_cagada.debug("el num ia mod en digitos {}".format(nueva_porcion))
                digitos_tam=(digitos_tam-min(digitos_max,digitos_tam)+len(nueva_porcion))
                digitos[idx_act:]=nueva_porcion
                #logger_cagada.debug("los digitos reducidos {}".format(digitos))
        num_final=digitos_a_numero(digitos,10)%modulo
        return num_final
        
def puta_mierda(digitos_base,digitos_exp, modulo):
#        base=reducir_a_modulo(digitos_base,320)
        base=reducir_a_modulo(digitos_base,modulo)
        #logger_cagada.debug("ke berg {}".format(base))
        # XXX: http://www.javascripter.net/math/calculators/eulertotientfunction.htm
        totiente=1000000006
        #logger_cagada.debug("el tototiente {}".format(totiente))
        exponente=reducir_a_modulo(digitos_exp,totiente)
        #logger_cagada.debug("el expon {}".format(exponente))
        res=mad_power(base%modulo,exponente,modulo)
        return res
        

def mierda():
        lineas=list(stdin)
        for linea in lineas[1:]:
                par_digitos=linea.strip().split(" ")
                digitos_base=[int(x) for x in reversed(par_digitos[0])]
                digitos_exp=[int(x) for x in reversed(par_digitos[1])]
                res=puta_mierda(digitos_base,digitos_exp,int(1E9+7))
                #logger_cagada.debug("el resultado {}".format(res))
                print("{}".format(res))

if __name__ == "__main__":
    FORMAT = "[%(filename)s:%(lineno)s - %(funcName)20s() ] %(message)s"
    logging.basicConfig(level=nivel_log, format=FORMAT)
    #logger_cagada = logging.getLogger("asa")
    #logger_cagada.setLevel(nivel_log)
    mierda()
