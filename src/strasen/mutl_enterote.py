'''
Created on 08/10/2017

@author: ernesto
'''
# XXX: https://quickgrid.wordpress.com/2015/11/16/uva-problem-10106-product-solution-lattice-multiplication/
from cmath import exp, pi
from operator import mul
import sys
import logging

nivel_log = logging.ERROR
#nivel_log = logging.DEBUG
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

def caca():
    lineas_cnt = 0
    lineas = [""] * 2
    for linea in sys.stdin:
        logger_cagada.debug("linea act {} lniea cnt {}".format(linea.strip(), lineas_cnt))
        if(lineas_cnt and not(lineas_cnt % 2)):
            logger_cagada.debug("lineas {}".format(lineas))
            num1 = enterote(lineas[0])
            num2 = enterote(lineas[1])
            numr = num1 * num2
            print("{}".format(numr))
        lineas_cnt += 1
        lineas[lineas_cnt % 2] = linea[:].strip()

    num1 = enterote(lineas[0])
    num2 = enterote(lineas[1])
    numr = num1 * num2
    print("{}".format(numr))

if __name__ == "__main__":
    FORMAT = "[%(filename)s:%(lineno)s - %(funcName)20s() ] %(message)s"
    logging.basicConfig(level=nivel_log, format=FORMAT)
    logger_cagada = logging.getLogger("asa")
    logger_cagada.setLevel(nivel_log)
    caca()
