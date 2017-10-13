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
from _operator import mul
from cmath import exp, pi

nivel_log = logging.ERROR
nivel_log = logging.DEBUG
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
#        logger_cagada.debug("la trans 1 {}".format(ent1_t))
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
        else:
            self.init_coeficientes(representacion)
        
    def init_cadena(self, cadena):
        coeficientes = [int(x) for x in cadena.strip().split(" ")]
        self.init_coeficientes(coeficientes)

    def init_coeficientes(self, coeficientes):
        self.coeficientes = coeficientes
        self.max_coef = max(self.coeficientes)
        self.max_exp = len(self.coeficientes)
        self.exp_10, self.pot_10 = poligamio_positivo.determina_pot_exp_min_10(self.max_exp * (self.max_coef ** 2))
        
    @classmethod
    def determina_pot_exp_min_10(clazz, numero):
        exp_10 = 0
        while numero > 10 ** exp_10:
            exp_10 += 1
        return exp_10, 10 ** exp_10
    
    @classmethod 
    def numero_a_digitos(clazz, num):
        digitos = []
        logger_cagada.debug("convirtiendo num {}".format(num))
        while num:
            digitos.append(num % 10)
            num //= 10
        logger_cagada.debug("kedo en digitos {}".format(digitos))
        return digitos
    
    @classmethod
    def numero_a_digitos_padeado(clazz, num, tam):
        digitos = poligamio_positivo.numero_a_digitos(num)
        logger_cagada.debug("el num {} kedo en digitos {} ".format(num, digitos))
        digitos += [0] * (tam - len(digitos))
        logger_cagada.debug(" padeado a {} kedo {}".format(tam, digitos))
        return digitos
    
    @classmethod
    def enterote_de_poligamio_positivo(clazz, poligamio, exp_10):
        digitos = []
        pol = poligamio.coeficientes
        for coef in pol:
            digitos += poligamio_positivo.numero_a_digitos_padeado(coef, exp_10)
        logger_cagada.debug("el polinomio {} kedo como digitos {}".format(poligamio, digitos))
        return enterote(digitos)
    
    @classmethod
    def digitos_a_numero(clazz, digitos):
        factor = 1
        num = 0
        for digito in digitos:
            num += digito * factor
            factor *= 10
        return num
        
# XXX: https://stackoverflow.com/questions/312443/how-do-you-split-a-list-into-evenly-sized-chunks
    @classmethod
    def poligamio_positivo_de_enterote(clazz, ent, exp_10):
        coefs = []
        if not exp_10:
            return poligamio_positivo([0])
        for i in range(0, len(ent.digitos), exp_10):
            digitos = ent.digitos[i:i + exp_10]
            coefs.append(poligamio_positivo.digitos_a_numero(digitos))
        logger_cagada.debug("el entero {} se paso a coefs {}".format(ent, coefs))
        return poligamio_positivo(coefs)
        
    
    def __mul__(self, other):
        max_exp = max(self.max_exp, other.max_exp)
        max_coef = max(self.max_coef, other.max_coef)
        exp_10, pot_10 = poligamio_positivo.determina_pot_exp_min_10(max_exp * max_coef ** 2)
        logger_cagada.debug("de pol {} y {} el exp 10 {}".format(self, other, exp_10))
        
        ent1 = poligamio_positivo.enterote_de_poligamio_positivo(self, exp_10)
        ent2 = poligamio_positivo.enterote_de_poligamio_positivo(other, exp_10)
        
        logger_cagada.debug("el polinom 1 {} kedo como ent {}".format(self, ent1))
        logger_cagada.debug("el polinom 2 {} kedo como ent {}".format(other, ent2))
        
        entr = ent1 * ent2
        
        pol = poligamio_positivo.poligamio_positivo_de_enterote(entr, exp_10)
        
        logger_cagada.debug("la mult de pol {} y {} resulta {}".format(self, other, pol))
        
        return pol
    
    def __add__(self, orto):
        coefs = [x[0] + x[1] for x in zip_longest(self.coeficientes, orto.coeficientes, fillvalue=0)]
        return poligamio_positivo(coefs)
    
    def __sub__(self, orto):
        coefs = [x[0] - x[1] for x in zip_longest(self.coeficientes, orto.coeficientes, fillvalue=0)]
        return poligamio_positivo(coefs)
        
        
    def __repr__(self):
        return "{}".format(self.coeficientes)

    __rmul__ = __mul__

    

class poligamio():
    def __init__(self, representacion):
        self.coeficientes = []
        self.polinomio_positivo = None
        self.polinomio_negativo = None
        if isinstance(representacion, str):
            self.init_de_cadena(representacion)
        else:
            self.init_de_coeficientes(representacion)
    
    def init_de_cadena(self, cadena):
        coeficientes = [int(x) for x in cadena.strip().split(" ")]
        self.init_de_coeficientes(coeficientes)
        
    def init_de_coeficientes(self, coeficientes):
        self.coeficientes = coeficientes
        poligamio.quita_sobrantes_coeficientes(self.coeficientes)
        coeficientes_positivos = [0] * len(coeficientes)
        coeficientes_negativos = [0] * len(coeficientes)
        for idx, coef in enumerate(coeficientes):
            if(coef < 0):
                coeficientes_negativos[idx] = -coef
            else:
                coeficientes_positivos[idx] = coef
        self.polinomio_positivo = poligamio_positivo(coeficientes_positivos)
        self.polinomio_negativo = poligamio_positivo(coeficientes_negativos)
        logger_cagada.debug("los coefs {} el pol p {} n {}".format(self.coeficientes, self.polinomio_positivo, self.polinomio_negativo))
    
    @classmethod
    def quita_sobrantes_coeficientes(cls, coeficientes):
#        ultimo_coef = 0
#        for idx, coef in enumerate(coeficientes):
#            if coef:
#                ultimo_coef = idx
#        return coeficientes[:ultimo_coef + 1]
        while coeficientes and coeficientes[-1] == 0:
            coeficientes.pop()   # normalize
    
    def __mul__(self, orto):
        polr = self.polinomio_positivo * orto.polinomio_positivo + self.polinomio_negativo * orto.polinomio_negativo
        logger_cagada.debug("el polr solo pos {}".format(polr))
        polr -= (self.polinomio_positivo * orto.polinomio_negativo + self.polinomio_negativo * orto.polinomio_positivo)
        logger_cagada.debug("el polr final {}".format(polr))
        return poligamio(polr.coeficientes)
    
    __rmul__ = __mul__
    
    def __repr__(self):
        cadena = ""
        for idx, coef in enumerate(self.coeficientes):
            cadena += "{}".format(coef)
            if(idx):
                cadena += "x^{}".format(idx)
            if(idx < len(self.coeficientes) - 1):
                cadena += " + "
        return cadena
#        return "{}".format(self.coeficientes)


    def __truediv__(self,orto):
        N=self.coeficientes
        D=orto.coeficientes
    
#    enterote.normalizar_a_tam(N,max_exp)
#    enterote.normalizar_a_tam(D,max_exp)
        poligamio.quita_sobrantes_coeficientes(D)
        dD = len(D)-1
        dN = len(N)-1
        if dD < 0: raise ZeroDivisionError
        if dN >= dD:
            q = [0] * dN
            while dN >= dD:
                d = [0]*(dN - dD) + D
                mult = q[dN - dD] = N[-1] // d[-1]
                d = [coeff*mult for coeff in d]
                N = [ coeffN - coeffd  for coeffN, coeffd in zip(N, d)]
                poligamio.quita_sobrantes_coeficientes(N)
                dN = len(N)-1
            r = N
        else:
            q = [0]
            r = N
        logger_cagada.debug("el q s {} l d {}".format(q,r))
        return poligamio(q),poligamio(r)
    
    __rtruediv__=__truediv__

def laconchadelamadre():
    lineas = list(sys.stdin)
#    polr = poligamio(lineas[1]) * poligamio(lineas[2])
    polq,polr = poligamio(lineas[1]) / poligamio(lineas[2])
    logger_cagada.debug("el p res {},{}".format(polq,polr))
    print("{}".format(polr))

if __name__ == "__main__":
    FORMAT = "[%(filename)s:%(lineno)s - %(funcName)20s() ] %(message)s"
    logging.basicConfig(level=nivel_log, format=FORMAT)
    logger_cagada = logging.getLogger("asa")
    logger_cagada.setLevel(nivel_log)
    laconchadelamadre()

