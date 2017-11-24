#!/programas/python3.5/bin/python3

# XXX: https://www.hackerrank.com/challenges/python-power-mod-power/problem
# XXX: https://www.khanacademy.org/computing/computer-science/cryptography/modarithmetic/a/fast-modular-exponentiation

import logging
from math import pow
from sys import stdin
from operator import mul
from functools import reduce

nivel_log = logging.ERROR
nivel_log = logging.DEBUG
logger_cagada = None

def sacacaca(x):
    powers = []
    i = 1
    while i <= x:
        if i & x:
            powers.append(i)
        i <<= 1
    return powers

def mad_power(a,b,m=None):
	res=0
	assert a<m
	assert b>=0
	pots_2=sacacaca(b)
	logger_cagada.debug("de num {} su descom {}".format(b,pots_2))
	res=reduce(lambda x,y: mul(x,y)%m,(map(lambda x:pow(a,x)%m,pots_2)),1)%m
	return res

def mierda():
	lineas=list(stdin)
	a=int(lineas[0].strip())
	b=int(lineas[1].strip())
	m=int(lineas[2].strip())
	res =mad_power(a%m,b,m)
	print("{}".format(int(pow(a,b))))
	print("{}".format(int(res)))

if __name__ == "__main__":
    FORMAT = "[%(filename)s:%(lineno)s - %(funcName)20s() ] %(message)s"
    logging.basicConfig(level=nivel_log, format=FORMAT)
    logger_cagada = logging.getLogger("asa")
    logger_cagada.setLevel(nivel_log)
    mierda()
