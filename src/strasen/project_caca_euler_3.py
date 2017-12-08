#!/programas/python3.5/bin/python3
'''
Created on 08/10/2017

@author: ernesto
'''
# XXX: 

from math import log,ceil
import sys
import logging
from itertools import zip_longest
from operator import mul
from cmath import exp, pi, acos
from ctypes import c_int
import re
from collections import defaultdict
from functools import partial

nivel_log = logging.ERROR
nivel_log = logging.DEBUG
logger_cagada = None


class enterote():
    def __init__(self, representacion, base):
        self.digitos = []
        self.digitos_tam = 0
        self.base = base
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
    
    dos_pi = 4 * acos(0)
    @classmethod
    def ffft_int(clazz, com_in, com_in_inicio, com_out, com_out_inicio, pasito, tam, direccion, exps):    
#        #logger_cagada.debug("l idx in ini {} l idx out ini {} el pasito {} el tam {}".format(com_in_inicio,com_out_inicio, pasito, tam))
        if tam == 1:
            com_out[com_out_inicio] = com_in[com_in_inicio]
            return
        tam_mitad = tam >> 1
        pasito_doble = pasito << 1
        enterote.ffft_int(com_in, com_in_inicio, com_out, com_out_inicio, pasito_doble, tam_mitad, direccion, exps)
        enterote.ffft_int(com_in, com_in_inicio + pasito, com_out, com_out_inicio + tam_mitad, pasito_doble, tam_mitad, direccion, exps)
#        #logger_cagada.debug("la salida {} el pasito {} el tam {}".format(com_out, pasito, tam))
        for i in range(tam_mitad):
            idx_out_par = i + com_out_inicio
            idx_out_impar = idx_out_par + tam_mitad
#            #logger_cagada.debug("idx out {} com ini {} tam mitad {}".format(idx_out,com_out_inicio,tam_mitad))
            com_par = com_out[idx_out_par]
            com_impar = com_out[idx_out_impar]
            exp1 = exps[i][tam]
            factor_caca = exp1 * com_impar
#            exp1tmp=exp(direccion*enterote.dos_pi*i*1j/tam)
#            assert exp1==exp1tmp, "el exp cache {} el otro {}".format(exp1,exp1tmp)
#            #logger_cagada.debug("el exp1 {} para meirda {} {} ".format(exp1, i,tam))
#            #logger_cagada.debug("el exp2 {} para meirda {}".format(exp2, i+tam_mitad))
            
            com_out[idx_out_par] = com_par + factor_caca
            com_out[idx_out_impar] = com_par - factor_caca
#            #logger_cagada.debug("calculando {} + {} * {} = {} en {}".format(com_par,exp1,com_impar,com_out[idx_out],idx_out))
#            #logger_cagada.debug("calculando {} - {} * {} = {} en {}".format(com_par,exp1,com_impar,com_out[idx_out+tam_mitad],idx_out+tam_mitad))

    @classmethod
    def iffft(clazz, com_in, com_out):    
        tam = len(com_in)
        exps = []
        for i in range(tam):
            mapita = {}
            exp_2_act = 2
            while(exp_2_act <= tam):
                mapita[exp_2_act] = exp(-1 * enterote.dos_pi * i * 1j / exp_2_act)
                exp_2_act <<= 1
            exps.append(mapita)
        enterote.ffft_int(com_in, 0, com_out, 0, 1, len(com_in), -1, exps)
        for i in range(len(com_in)):
            com_out[i] /= len(com_in)

    @classmethod
    def ffft(clazz, com_in, com_out):    
        tam = len(com_in)
        exps = []
        for i in range(tam):
            mapita = {}
            exp_2_act = 2
            while(exp_2_act <= tam):
                mapita[exp_2_act] = exp(enterote.dos_pi * i * 1j / exp_2_act)
                exp_2_act <<= 1
            exps.append(mapita)
        enterote.ffft_int(com_in, 0, com_out, 0, 1, tam, 1, exps)

# XXX: https://rosettacode.org/wiki/Fast_Fourier_transform#Python:_Recursive
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
#        #logger_cagada.debug("normalizado de {} a {}".format(self.digitos_tam, tam_normalizado))
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
        assert self.base == otro.base
        enterote.normalizar_para_fft(self, otro)
#        #logger_cagada.debug("ent 1 normalizado a {}".format(self.digitos))
#        #logger_cagada.debug("ent 2 normalizado a {}".format(otro.digitos))
        ent1_t = [0j] * len(self.digitos)
        ent2_t = [0j] * len(self.digitos)
        # logger_cagada.debug("ent1 redondeado {}".format(ent1_t))
        # logger_cagada.debug("ent2 redondeado {}".format(ent2_t))
        enterote.ffft(self.digitos, ent1_t)
#        #logger_cagada.debug("la trans 1 {}".format(ent1_t))
    # print("etn1 t {}".format(ent1_t))
        enterote.ffft(otro.digitos, ent2_t)
#        #logger_cagada.debug("la trans 2 {}".format(ent2_t))
        entr_t = list(map(mul, ent1_t, ent2_t))
        # logger_cagada.debug("etnr t {}".format(entr_t))
        entr_t_r = [0j] * len(entr_t)
        enterote.iffft(entr_t, entr_t_r)
        entr_tmp = enterote.parte_real_redondeada_de_complejos(entr_t_r)
        # logger_cagada.debug("resultado bryto {}".format(entr_tmp))
        entr = entr_tmp[:]
        # logger_cagada.debug("resultado tmp {}".format(entr_tmp))
        for idx in range(len(entr) - 1):
            coef = entr[idx]
            entr[idx] = coef % self.base
            entr[idx + 1] += coef // self.base
        # logger_cagada.debug("resultado ya arregladito {}".format(entr))
        return enterote(entr, self.base)
        
    def __repr__(self):
        imprimir = False
        cadena = ""
        for coef in reversed(self.digitos):
            if(coef):
                imprimir = True
            if(imprimir):
                cadena += "{},".format(coef)
        if(not imprimir):
            cadena += "{}".format(0)
        cadena += " base {}".format(self.base)
        return cadena
    
    __str__ = __repr__


class poligamio_positivo():
    def __init__(self, representacion, base):
        self.coeficientes = []
        self.max_exp = 0
        self.max_coef = 0
        self.exp_10 = 0
        self.pot_10 = 0
#        self.base=256
        self.base = base
        if isinstance(representacion, str):
            self.init_cadena(representacion)
        else:
            self.init_coeficientes(representacion)
        
    def init_cadena(self, cadena):
        coeficientes = [int(x) for x in cadena.strip().split(" ")]
        self.init_coeficientes(coeficientes)

    def init_coeficientes(self, coeficientes):
        self.coeficientes = coeficientes if coeficientes else [0]
        self.max_coef = max(self.coeficientes)
        self.max_exp = len(self.coeficientes)
#        self.exp_10, self.pot_10 = poligamio_positivo.determina_pot_exp_min_10(self.max_exp * (self.max_coef ** 2, self.base))
        
    @classmethod
    def determina_pot_exp_min_base(clazz, numero, base):
        exp_10 = 0
        while numero > base ** exp_10:
            exp_10 += 1
        return exp_10, base ** exp_10
    
    @classmethod 
    def numero_a_digitos(clazz, num, base):
        digitos = []
        # logger_cagada.debug("convirtiendo num {}".format(num))
        while num:
            digitos.append(num % base)
            num //= base
        # logger_cagada.debug("kedo en digitos {}".format(digitos))
        return digitos
    
    @classmethod
    def numero_a_digitos_padeado(clazz, num, base, pad):
        digitos = poligamio_positivo.numero_a_digitos(num, base)
        # logger_cagada.debug("el num {} kedo en digitos {} ".format(num, digitos))
        digitos += [0] * (pad - len(digitos))
        # logger_cagada.debug(" padeado a {} kedo {}".format(pad, digitos))
        return digitos
    
    @classmethod
    def enterote_de_poligamio_positivo(clazz, poligamio, base, pad):
        digitos = []
        pol = poligamio.coeficientes
        for coef in pol:
            digitos += poligamio_positivo.numero_a_digitos_padeado(coef, base, pad)
        # logger_cagada.debug("el polinomio {} kedo como digitos {}".format(poligamio, digitos))
        return enterote(digitos, base)
    
    @classmethod
    def digitos_a_numero(clazz, digitos, base, modulo):
        factor = 1
        num = 0
        base_mod = base % modulo
        for digito in digitos:
            num += ((digito % modulo) * factor) % modulo
            factor = (factor * base_mod) % modulo
        return num
        
# XXX: https://stackoverflow.com/questions/312443/how-do-you-split-a-list-into-evenly-sized-chunks
    @classmethod
    def poligamio_positivo_de_enterote(clazz, ent, base, exp):
        coefs = []
        if not exp:
            return poligamio_positivo([ent.digitos[0]], base)
        for i in range(0, len(ent.digitos), exp):
            digitos = ent.digitos[i:i + exp]
            coefs.append(poligamio_positivo.digitos_a_numero(digitos, base, base))
        # logger_cagada.debug("el entero {} se paso a coefs {} de base {}".format(ent, coefs, base))
        return poligamio_positivo(coefs, base)
        
    
    # profile
    def __mul__(self, other):
        base = self.base
        max_exp = max(self.max_exp, other.max_exp)
        max_coef = max(self.max_coef, other.max_coef)
        max_coef_esperado = max_exp * max_coef ** 2
        exp_10, pot_10 = poligamio_positivo.determina_pot_exp_min_base(max_coef_esperado, base)
        # logger_cagada.debug("de pol {} y {} el exp base {} es {} para maximo coef sperado {}".format(self, other, base, exp_10,hex(max_coef_esperado)))
        
        ent1 = poligamio_positivo.enterote_de_poligamio_positivo(self, base, exp_10)
        ent2 = poligamio_positivo.enterote_de_poligamio_positivo(other, base, exp_10)
        
        # logger_cagada.debug("el polinom 1 {} kedo como ent {}".format(self, ent1))
        # logger_cagada.debug("el polinom 2 {} kedo como ent {}".format(other, ent2))
        # logger_cagada.debug("el ent1 tam {}".format(len(ent1.digitos)))
        # logger_cagada.debug("el ent2 tam {}".format(len(ent2.digitos)))
        
        entr = ent1 * ent2
        
        pol = poligamio_positivo.poligamio_positivo_de_enterote(entr, base, exp_10)
        
        # logger_cagada.debug("la mult de pol {} y {} resulta {}".format(self, other, pol))
        
        return pol
    
    def __add__(self, orto):
        assert self.base == orto.base
        modulo = self.base
        coefs = [(x[0] % modulo + x[1] % modulo) % modulo for x in zip_longest(self.coeficientes, orto.coeficientes, fillvalue=0)]
        return poligamio_positivo(coefs, self.base)
    
    def __sub__(self, orto):
        assert self.base == orto.base
        modulo = self.base
        coefs = [(x[0] % modulo - x[1] % modulo) % modulo for x in zip_longest(self.coeficientes, orto.coeficientes, fillvalue=0)]
        return poligamio_positivo(coefs, self.base)
        
        
    def __repr__(self):
        return "{}".format(self.coeficientes)

    __rmul__ = __mul__

    

class poligamio():
    def __init__(self, representacion, base):
        self.coeficientes = []
        self.polinomio_positivo = None
        self.polinomio_negativo = None
        self.base = base
        if isinstance(representacion, str):
            self.init_de_cadena(representacion)
        else:
            self.init_de_coeficientes(representacion)
        self.formato_mamalon = True
    
    def init_de_cadena(self, cadena):
        coeficientes = [int(x) for x in cadena.strip().split(" ")]
        assert all(map(lambda x:x < self.base, coeficientes))
        self.init_de_coeficientes(coeficientes)
        
    def init_de_coeficientes(self, coeficientes):
        self.coeficientes = coeficientes
        poligamio.quita_sobrantes_coeficientes(self.coeficientes)
        if not self.coeficientes:
            self.coeficientes = [0]
        coeficientes_positivos = [0] * len(coeficientes)
        coeficientes_negativos = [0] * len(coeficientes)
        for idx, coef in enumerate(coeficientes):
            if(coef < 0):
                coeficientes_negativos[idx] = -coef
            else:
                coeficientes_positivos[idx] = coef
        self.polinomio_positivo = poligamio_positivo(coeficientes_positivos, self.base)
        self.polinomio_negativo = poligamio_positivo(coeficientes_negativos, self.base)
#        #logger_cagada.debug("los coefs {} el pol p {} n {}".format(self.coeficientes, self.polinomio_positivo, self.polinomio_negativo))
    
    @classmethod
    def quita_sobrantes_coeficientes(cls, coeficientes):
#        ultimo_coef = 0
#        for idx, coef in enumerate(coeficientes):
#            if coef:
#                ultimo_coef = idx
#        return coeficientes[:ultimo_coef + 1]
        while coeficientes and coeficientes[-1] == 0:
            coeficientes.pop()  # normalize
    
    # profile
    def __mul__(self, orto):
        polr = self.polinomio_positivo * orto.polinomio_positivo + self.polinomio_negativo * orto.polinomio_negativo
        # logger_cagada.debug("el polr solo pos {}".format(polr))
        polr -= (self.polinomio_positivo * orto.polinomio_negativo + self.polinomio_negativo * orto.polinomio_positivo)
        # logger_cagada.debug("el polr final {}".format(polr))
        return poligamio(polr.coeficientes, self.base)
    
    __rmul__ = __mul__
    
    def __repr__(self):
        cadena = ""
        for idx, coef in enumerate(self.coeficientes):
#            cadena += "{}".format(c_int(coef).value)
            cadena += "{}".format(coef)
            if(idx):
                if self.formato_mamalon:
                    cadena += "x^{}".format(idx)
            if(idx < len(self.coeficientes) - 1):
                if self.formato_mamalon:
                    cadena += " + "
                else:
                    cadena += " "
        return cadena
#        return "{}".format(self.coeficientes)

    def es_divisor_sintetico(self):
        return False
        if self.coeficientes[0] == self.coeficientes[-1] == 1 and not any(self.coeficientes[1:-1]):
            return True
        else:
            return False
        
    # profile
    def __truediv__(self, orto):
        N = self.coeficientes
        D = orto.coeficientes
        assert self.base == orto.base
        modulo = self.base
        # logger_cagada.debug("dividendo {} divisor {}".format(self, orto,))
    
#    enterote.normalizar_a_tam(N,max_exp)
#    enterote.normalizar_a_tam(D,max_exp)
        poligamio.quita_sobrantes_coeficientes(D)
        dN = len(N) - 1
        dD = len(D) - 1
        if(orto.es_divisor_sintetico()):
            return self.div_sintetica(orto)
        # logger_cagada.debug("usando div normal {} entre {}".format(self,orto))
        # logger_cagada.debug("dN {} dD {}".format(dN, dD))
        if dD < 0: raise ZeroDivisionError
        if dN >= dD:
            q = [0] * (dN + 1)
            while dN >= dD:
                dividendo_ant = N[:]
                d = [0] * (dN - dD) + D
                mult = q[dN - dD] = (N[-1] // d[-1]) % modulo
                # logger_cagada.debug("l mult es {}" .format(mult))
                d = [((coeff % modulo) * mult) % modulo for coeff in d]
                N = [ (coeffN - coeffd) % modulo  for coeffN, coeffd in zip(N, d)]
                poligamio.quita_sobrantes_coeficientes(N)
                # logger_cagada.debug("aora N es {}".format(N))
                if(N == dividendo_ant):
                    break
                dN = len(N) - 1
            r = N
        else:
            q = [0]
            r = N
            # logger_cagada.debug("nada q acer r es {}".format(r))
        poligamio.quita_sobrantes_coeficientes(q)
#        if(not q):
#            q=[0]
        # logger_cagada.debug("el q s {} l d {}".format(q, r))
        return poligamio(q, modulo), poligamio(r, modulo)
    
    __rtruediv__ = __truediv__

    def div_sintetica(self, orto):
        # logger_cagada.debug("usando div sintactica {} entre {}".format(self,orto))
        dividendo = self.coeficientes
        divisor = orto.coeficientes
        poligamio.quita_sobrantes_coeficientes(divisor)
        grado_dividendo = len(dividendo) - 1
        grado_divisor = len(divisor) - 1
        r = [0] * (grado_divisor)
        if grado_divisor < 0: raise ZeroDivisionError
        q = [0]
        if grado_dividendo >= grado_divisor:
            grado_act = 0
            while grado_act <= grado_dividendo:
                grado_contraparte = grado_act + grado_divisor
                r[grado_act % grado_divisor] += dividendo[grado_act]
                # logger_cagada.debug("sumando a grado {} {} de {}".format(grado_act % grado_divisor, dividendo[grado_act], grado_act))

                if(grado_contraparte <= grado_dividendo):
                    r[grado_act % grado_divisor] -= dividendo[grado_contraparte]
                    # logger_cagada.debug("restando a grado {} {} de {}".format(grado_act % grado_divisor, dividendo[grado_contraparte], grado_contraparte))

                grado_act += 1
                if(not (grado_act % grado_divisor)):
                    # logger_cagada.debug("brincando de grado {} a {}".format(grado_act, grado_divisor + grado_act))
                    grado_act += grado_divisor
        else:
            r = dividendo
        return poligamio(q), poligamio(r)

    @property
    def grado(self):
        poligamio.quita_sobrantes_coeficientes(self.coeficientes)
        grado = len(self.coeficientes)
        if grado <= 0:
            grado = 1
            self.coeficientes = [0]
        return grado - 1

class nodo_arbol():
    def __init__(self, valor, hijo_izq=None, hijo_der=None):
        self.valor = valor
        self.hijo_izq = hijo_izq
        self.hijo_der = hijo_der
    def __repr__(self):
        return "!!!valor {}\nhijo izq {}\nhijo der {}!!!".format(self.valor, self.hijo_izq, self.hijo_der)
    
class funcionsilla():
    def __init__(self, coeficientes, modulo, base):
        self.polinomio = None
        self.modulo = modulo
        self.inicializa(coeficientes, base)
    
    def inicializa(self, coeficientes, base):
        modulo = self.modulo
        self.polinomio = poligamio(list(map(lambda x:x % modulo, coeficientes)), base)
    
    def evalua(self, putos):
        modulo = self.modulo
        p = self.polinomio
        raiz_arbolin = self.genera_arbolin_producto(putos)
        evaluaciones = defaultdict(lambda:p.coeficientes[0] % modulo)
        logger_cagada.debug("el puto arbol\n{}".format(raiz_arbolin))
        logger_cagada.debug("putos {}".format(putos))
        logger_cagada.debug("el pendejo {}".format(p))
        self.eval_multicaca_traversea(raiz_arbolin, p, evaluaciones)
        # logger_cagada.debug("las evaluaciones {}".format(evaluaciones))
        return evaluaciones

    def eval_multicaca_traversea(self, nodo, residuo_ant, evaluaciones):
        if(not nodo.valor.grado):
            return
        _, mierda = residuo_ant / nodo.valor
        # logger_cagada.debug("para l nodo {} grado {} el res ant {} i el res {}:{}".format(nodo.valor, nodo.valor.grado, residuo_ant, caca, mierda.grado))
        if nodo.valor.grado == 1:
            evaluaciones[-nodo.valor.coeficientes[0]] = mierda
        if nodo.hijo_izq and nodo.hijo_izq.valor.grado and mierda.grado:
            self.eval_multicaca_traversea(nodo.hijo_izq, mierda, evaluaciones)
        if nodo.hijo_der and nodo.hijo_der.valor.grado and mierda.grado:
            self.eval_multicaca_traversea(nodo.hijo_der, mierda, evaluaciones)

# profile

# profile
    def genera_arbolin_producto(self, numeros):
        modulo = self.modulo
        numeros_tam = len(numeros)
        tam_normalizado = enterote.determina_pot_2_minima(numeros_tam)
        numeros_normalizados = numeros + [sys.maxsize] * (tam_normalizado - numeros_tam)
        raiz = self.genera_arbolin_product_recursivo(numeros_normalizados)
        return raiz

# profile
    def genera_arbolin_product_recursivo(self, numeros):
        nodo_act = None
        numeros_tam = len(numeros)
        modulo = self.modulo
        # logger_cagada.debug("los nums {}".format(numeros))
        if(numeros_tam > 1):
            hijo_izq = self.genera_arbolin_product_recursivo(numeros[:numeros_tam // 2])
            hijo_der = self.genera_arbolin_product_recursivo(numeros[numeros_tam // 2:])
            if hijo_izq.valor.grado == 0 and hijo_izq.valor.coeficientes[0] == 1:
                nodo_act = nodo_arbol(hijo_der.valor, hijo_izq, hijo_der)
            else:
                if hijo_der.valor.grado == 0 and hijo_der.valor.coeficientes[0] == 1:
                    nodo_act = nodo_arbol(hijo_izq.valor, hijo_izq, hijo_der)
                else:
                    nodo_act = nodo_arbol(hijo_izq.valor * hijo_der.valor, hijo_izq, hijo_der)
            # logger_cagada.debug("el pol {} viene de {} por {}".format(nodo_act.valor, hijo_izq.valor, hijo_der.valor))
        else:
            if(numeros[0] != sys.maxsize):
                nodo_act = nodo_arbol(poligamio([(-numeros[0]), 1], modulo))
                # logger_cagada.debug("el pol single {} el num {}".format(nodo_act.valor, numeros[0]))
            else:
                nodo_act = nodo_arbol(poligamio([1], modulo))
        return nodo_act

    @classmethod
    def genera_polinomio_recursivo(clazz, numeros,modulo):
        nodo_act = None
        numeros_tam = len(numeros)
        # logger_cagada.debug("los nums {}".format(numeros))
        if(numeros_tam > 1):
            hijo_izq = funcionsilla.genera_polinomio_recursivo(numeros[:numeros_tam // 2],modulo)
            hijo_der = funcionsilla.genera_polinomio_recursivo(numeros[numeros_tam // 2:],modulo)
            if hijo_izq.valor.grado == 0 and hijo_izq.valor.coeficientes[0] == 1:
                nodo_act = nodo_arbol(hijo_der.valor, hijo_izq, hijo_der)
            else:
                if hijo_der.valor.grado == 0 and hijo_der.valor.coeficientes[0] == 1:
                    nodo_act = nodo_arbol(hijo_izq.valor, hijo_izq, hijo_der)
                else:
                    nodo_act = nodo_arbol(hijo_izq.valor * hijo_der.valor, hijo_izq, hijo_der)
            # logger_cagada.debug("el pol {} viene de {} por {}".format(nodo_act.valor, hijo_izq.valor, hijo_der.valor))
        else:
            if(numeros[0] != sys.maxsize):
                nodo_act = nodo_arbol(poligamio([(numeros[0]), 1], modulo))
                # logger_cagada.debug("el pol single {} el num {}".format(nodo_act.valor, numeros[0]))
            else:
                nodo_act = nodo_arbol(poligamio([1], modulo))
        return nodo_act

    @classmethod
    def genera_polinomio(clazz, numeros,modulo):
        numeros_tam = len(numeros)
        tam_normalizado = enterote.determina_pot_2_minima(numeros_tam)
        numeros_normalizados = numeros + [sys.maxsize] * (tam_normalizado - numeros_tam)
        rais=funcionsilla.genera_polinomio_recursivo(numeros,modulo)
        return rais.valor

def mad_power(a,b,m=None):
    res=1
#    logger_cagada.debug("asdd {}".format([a,b,m]))
    assert a<m
    assert b>=0
    pot=a
    b_tmp=b
    while b_tmp:
        if b_tmp&1:
            res=(res*pot)%m
        pot=(pot*pot)%m
        b_tmp>>=1
#    logger_cagada.debug("pero q mierda {}".format(res))
    return res

def mult_con_mod(a,b,m):
#    logger_cagada.debug("asdd {}".format([a,b,m]))
    res=((a%m)*(b%m))%m
#    logger_cagada.debug("pero q mierda {}".format(res))
    return res

def sumar_con_mod(a,b,m):
#    logger_cagada.debug("asdd {}".format([a,b,m]))
    res=((a%m)+(b%m))%m
#    logger_cagada.debug("pero q mierda {}".format(res))
    return res

def cagar_mierda(b, c, d, e, m,k):
#    logger_cagada.debug("ass {}".format([b, c, d, e, m,k]))
    pote=partial(mad_power,m=m)
    suma=partial(sumar_con_mod,m=m)
    multi=partial(mult_con_mod,m=m)
#    res=sumar_con_mod(sumar_con_mod(mult_con_mod(b,mad_power(c,4*k,m),m),mult_con_mod(d,mad_power(c,2*k,m),m),m),e,m)
    res=suma(suma(multi(b,pote(c,4*k)),multi(d,pote(c,2*k))),e)
#    logger_cagada.debug("m corto los webs {}".format(res))
    return res
    

def pce_genera_putos(n):
    d=ceil(n**(1/4))
    return list(range(1,d+1))

def pce_core(n):
    numeros=pce_genera_putos(n)
    modulo=n
    putos=[0]+list(numeros[:-1])
    poli=funcionsilla.genera_polinomio(numeros,modulo)
    logger_cagada.debug("n {} los nums {} los puts {} el pol {}" .format(n,list(numeros),list(putos),poli))
    funcion_caca = funcionsilla(poli.coeficientes, modulo, modulo)
    evaluaciones = funcion_caca.evalua(putos)
    logger_cagada.debug("las evals {}".format(evaluaciones))


def pce_main():
    t = int(input().strip())
    for a0 in range(t):
        n = int(input().strip())
        pce_core(n)


if __name__ == "__main__":
    FORMAT = "[%(filename)s:%(lineno)s - %(funcName)20s() ] %(message)s"
    logging.basicConfig(level=nivel_log, format=FORMAT)
    logger_cagada = logging.getLogger("asa")
    logger_cagada.setLevel(nivel_log)
    pce_main()
