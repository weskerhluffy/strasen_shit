'''
Created on 08/10/2017

@author: ernesto
'''

# XXX: https://www.codechef.com/problems/MULTIPLY
# XXX: http://web.maths.unsw.edu.au/~davidharvey/talks/kronecker-talk.pdf

from math import log
import sys
import logging
from asyncio.log import logger
from itertools import zip_longest

nivel_log = logging.ERROR
# nivel_log = logging.DEBUG
logger_cagada = None

class enterote():
    def __init__(self, representacion):
        self.digitos = []
        self.digitos_tam = 0
        if isinstance(representacion, str):
            self.init_de_cadeana(representacion)
        else:
            self.init_de_digitos(representacion)
    
    def init_de_cadeana(self, cadenota):
        self.digitos = list(reversed(list(map(int, cadenota.strip()))))
        self.digitos_tam = len(self.digitos)

    def init_de_digitos(self, digitos):
        self.digitos = digitos
        self.digitos_tam = len(self.digitos)
    
    @classmethod
    def fft(clazz, x, direccion=1):    
        N = len(x)
        if N <= 1: return x
        even = enterote.fft(x[0::2], direccion)
        odd = enterote.fft(x[1::2], direccion)
        T = [exp((direccion * -2j * pi * k) / N) * odd[k] for k in range(N // 2)]
        return [even[k] + T[k] for k in range(N // 2)] + [even[k] - T[k] for k in range(N // 2)]
    
    @classmethod
    def ifft(clazz, x):
        a = enterote.fft(x, -1)
        tam_x = len(x)
        b = [nu / float(tam_x) for nu in a]
        return b
    
    @classmethod
    def determina_pot_2_minima(clazz, num):
        exp_2 = 0
        while num > 2 ** exp_2:
            exp_2 += 1
        return 2 ** exp_2
    
    def normalizar_a_tam(self, tam):
        tam_normalizado = tam
        if(tam_normalizado <= self.digitos_tam):
            return
        self.digitos += [0] * (tam_normalizado - self.digitos_tam)
        logger_cagada.debug("normalizado de {} a {}".format(self.digitos_tam, tam_normalizado))
        self.digitos_tam = tam_normalizado
    
    @classmethod
    def normalizar_para_fft(clazz, ent1, ent2):
        max_tam = max(ent1.digitos_tam, ent2.digitos_tam)
        tam_normalizado = enterote.determina_pot_2_minima(max_tam) << 1
        ent1.normalizar_a_tam(tam_normalizado)
        ent2.normalizar_a_tam(tam_normalizado)
    
    @classmethod
    def parte_real_de_complejos(clazz, x):
        return [caca.real for caca in x]

    @classmethod
    def parte_real_redondeada_de_complejos(clazz, x):
        return [round(mierda) for mierda in enterote.parte_real_de_complejos(x)]
    
    def __mul__(self, otro):
        enterote.normalizar_para_fft(self, otro)
        logger_cagada.debug("ent 1 normalizado a {}".format(self.digitos))
        logger_cagada.debug("ent 2 normalizado a {}".format(otro.digitos))
#    print("ent1 redondeado {}".format(ent1))
#    print("ent2 redondeado {}".format(ent2))
        ent1_t = enterote.fft(self.digitos)
        logger_cagada.debug("la trans 1 {}".format(ent1_t))
    # print("etn1 t {}".format(ent1_t))
        ent2_t = enterote.fft(otro.digitos)
        entr_t = list(map(mul, ent1_t, ent2_t))
#    print("etnr t {}".format(entr_t))
        entr_tmp = enterote.parte_real_redondeada_de_complejos(enterote.ifft(entr_t))
        logger_cagada.debug("resultado bryto {}".format(entr_tmp))
        entr = entr_tmp[:]
#    print("resultado tmp {}".format(entr_tmp))
        for idx in range(len(entr) - 1):
            coef = entr[idx]
            entr[idx] = coef % 10
            entr[idx + 1] += coef // 10
        logger_cagada.debug("resultado ya arregladito {}".format(entr))
        return enterote(entr)
        
    def __repr__(self):
        imprimir = False
        cadena = ""
        for coef in reversed(self.digitos):
            if(coef):
                imprimir = True
            if(imprimir):
                cadena += "{}".format(coef)
        if(not imprimir):
            cadena += "{}".format(0)
        return cadena
    
    __str__ = __repr__

class poligamio_positivo():
    def __init__(self, representacion):
        self.coeficientes = []
        self.max_exp = 0
        self.max_coef = 0
        self.exp_10 = 0
        self.pot_10 = 0
        if isinstance(representacion, str):
            self.init_cadena(representacion)
        
    def init_cadena(self, cadena):
        self.coeficientes = [int(x) for x in cadena.strip().split(" ")]
        self.max_coef = max(self.coeficientes)
        self.max_exp = len(self.coeficientes)
        self.exp_10, self.pot_10 = poligamio_positivo.determina_pot_min_10(self.max_exp * (self.max_coef ** 2))

    def init_coeficientes(self, coeficientes):
        self.coeficientes = coeficientes
        self.max_coef = max(self.coeficientes)
        self.max_exp = len(self.coeficientes)
        self.exp_10, self.pot_10 = poligamio_positivo.determina_pot_min_10(self.max_exp * (self.max_coef ** 2))
        
    @classmethod
    def determina_pot_exp_min_10(clazz, numero):
        exp_10 = 0
        while numero > 10 ** exp_10:
            exp_10 += 1
        return exp_10, 10 ** exp_10
    
    @staticmethod 
    def numero_a_digitos(clazz, num):
        digitos = []
        while num:
            digitos.append(num % 10)
            num //= 10
        return digitos
    
    @staticmethod
    def numero_a_digitos_padeado(clazz, num, tam):
        digitos = poligamio_positivo.numero_a_digitos(num)
        digitos += digitos + [0] * (tam - len(digitos))
        return digitos
    
    @staticmethod
    def enterote_de_poligamio_positivo(clazz, poligamio, exp_10):
        digitos = []
        pol = poligamio.coeficientes
        for coef in pol:
            digitos += poligamio_positivo.numero_a_digitos_padeado(coef, exp_10)
        logger_cagada.debug("el polinomio {} kedo como digitos {}".format(poligamio, digitos))
        return enterote(digitos)
    
    @staticmethod
    def digitos_a_numero(clazz, digitos):
        factor = 1
        num = 0
        for digito in digitos:
            num += digito * factor
            factor *= 10
        return num
        
# XXX: https://stackoverflow.com/questions/312443/how-do-you-split-a-list-into-evenly-sized-chunks
    @staticmethod
    def poligamio_positivo_de_enterote(clazz, ent, exp_10):
        coefs = []
        for i in range(0, len(ent.digitos), exp_10):
            digitos = ent.digitos[i:i + exp_10]
            coefs.append(poligamio_positivo.digitos_a_numero(digitos))
        logger_cagada.debug("el entero {} se paso a coefs {}".format(ent, coefs))
        return poligamio_positivo(coefs)
        
    
    def __mul__(self, other):
        max_exp = max(self.max_exp, other.max_exp)
        max_coef = max(self.max_coef, other.max_coef)
        exp_10, pot_10 = poligamio_positivo.determina_pot_exp_min_10(max_exp * max_coef ** 2)
        
        ent1 = poligamio_positivo.enterote_de_poligamio_positivo(self, exp_10)
        ent2 = poligamio_positivo.enterote_de_poligamio_positivo(other, exp_10)
        
        logger_cagada.debug("el polinom 1 {} kedo como ent {}".format(self, ent1))
        logger_cagada.debug("el polinom 2 {} kedo como ent {}".format(self, ent2))
        
        entr = ent1 * ent2
        
        pol = poligamio_positivo.poligamio_positivo_de_enterote(entr, exp_10)
        
        logger_cagada.debug("la mult de pol {} y {} resulta {}".format(self, other, pol))
        
        return pol
    
    def __add__(self, orto):
        coefs = [x[0] + x[1] for x in zip_longest(self.digitos, orto.digitos, fillvalue=0)]
        return poligamio_positivo(coefs)
    
    def __sub__(self, orto):
        coefs = [x[0] - x[1] for x in zip_longest(self.digitos, orto.digitos, fillvalue=0)]
        return poligamio_positivo(coefs)
        
        
    def __repr__(self):
        return "{}".format(self.num)

    __rmul__ = __mul__
    

def determina_pot_min_10(numero):
    exp_10 = 0
    while numero > 10 ** exp_10:
        exp_10 += 1
    return 10 ** exp_10

def determina_pot_min_10_polinomio(polinomio):
    max_exp = len(polinomio)
    max_coef = max(polinomio)
    pot_10 = determina_pot_min_10(max_exp * (max_coef ** 2))
    return pot_10

def empaca_polinomio(pol, pot_10):
    exp_pot_10 = 0
    num = 0
    for coef in pol:
        num += coef * (pot_10 ** exp_pot_10)
        exp_pot_10 += 1
    return numerin(num, pot_10)

def desempaca_polinomio(nume):
    num = nume.num
    pot_10 = nume.pot_10
    pol = []
    while num:
        pol.append(num % pot_10)
        num //= pot_10
    return pol if pol else [0] * int(log(pot_10, 10))

def multiplica_polinomios(pol1, pol2):
    pot_10 = max(determina_pot_min_10_polinomio(pol1), determina_pot_min_10_polinomio(pol2))
    # print("pot 10 max es {}".format(pot_10))
    pol1 = list(pol1)
    pol2 = list(pol2)

    num1 = empaca_polinomio(pol1, pot_10)
    num2 = empaca_polinomio(pol2, pot_10)

    # print("pol 1 {} num 1 {}".format(list(pol1),num1))
    # print("pol 2 {} num 2 {}".format(list(pol2),num2))

    numr = num1 * num2

    polr = desempaca_polinomio(numr)

    return polr


def completa_polinomio(pol, tam):
    return pol + [0] * (tam - len(pol))

def homogeiniza_polimonios(pol1, pol2):
    max_exp = max(len(pol1), len(pol2))
    return completa_polinomio(pol1, max_exp), completa_polinomio(pol2, max_exp)

def suma_polinomios(pol1, pol2):
    pol1, pol2 = homogeiniza_polimonios(pol1, pol2)
    # print("sumando {} y {}".format(list(reversed(pol1)),list(reversed(pol2))))
    polr = list(map(lambda x, y:x + y, pol1, pol2))
    return polr

def resta_polinomios(pol1, pol2):
    pol1, pol2 = homogeiniza_polimonios(pol1, pol2)
    polr = list(map(lambda x, y:x - y, pol1, pol2))
    return polr

def multiplica_polinomios_signados(pol1_p, pol1_n, pol2_p, pol2_n):
    polr_p = suma_polinomios(multiplica_polinomios(pol1_p, pol2_p), multiplica_polinomios(pol1_n, pol2_n))
    polr_n = suma_polinomios(multiplica_polinomios(pol1_p, pol2_n), multiplica_polinomios(pol1_n, pol2_p))
    # print("el pol p {} el n {}".format(list(reversed(polr_p)),list(reversed(polr_n))))
    polr = resta_polinomios(polr_p, polr_n)
    # print("podria ser q al fina {}".format(list(reversed(polr))))
    return polr

def separa_polinomio_por_signo(pol):
    pol_p = [0] * len(pol)
    pol_n = [0] * len(pol)
    for idx, coef in enumerate(pol):
        if(coef < 0):
            pol_n[idx] = -coef
        else:
            pol_p[idx] = coef
    return pol_p, pol_n


def multiplica_polinomios_con_signo(pol1, pol2):
    pol1, pol2 = homogeiniza_polimonios(pol1, pol2)

    pol1_p, pol1_n = separa_polinomio_por_signo(pol1)
    pol2_p, pol2_n = separa_polinomio_por_signo(pol2)

    polr = multiplica_polinomios_signados(pol1_p, pol1_n, pol2_p, pol2_n)

    return polr

def unados():
    lineas = list(sys.stdin)
    pol1 = [int(x) for x in lineas[1].strip().split(" ")]
    pol2 = [int(x) for x in lineas[2].strip().split(" ")]
    polr = multiplica_polinomios_con_signo(pol1, pol2)
    # print("{}".format(polr))
    for idx, coef in enumerate(polr):
        print("{}".format(coef), end="")
        if(idx):
            print("x^{}".format(idx), end="")
        if(idx < len(polr) - 1):
            print(" + ", end="")

    print("")

if __name__ == "__main__":
    FORMAT = "[%(filename)s:%(lineno)s - %(funcName)20s() ] %(message)s"
    logging.basicConfig(level=nivel_log, format=FORMAT)
    logger_cagada = logging.getLogger("asa")
    logger_cagada.setLevel(nivel_log)
    unados()

