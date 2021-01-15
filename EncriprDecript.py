import numpy as np
import time
import sys
import random 
import itertools
def random_nonsingular_matrix(size,base_ring=GF(2)): 
    V = base_ring**size
    vectors = []
    for i in range(size):
        v = V.random_element()
        while v in V.span(vectors):
            v = V.random_element()
        vectors.append(v)
    return(matrix(vectors))


def Gsplit(p):

        Phi = p.parent()
        p0 = Phi([sqrt(c) for c in p.list()[0::2]])
        p1 = Phi([sqrt(c) for c in p.list()[1::2]])
        return (p0,p1)
def Gnorm(a,b):

    X = g.parent().gen()
    return 2^((a^2+X*b^2).degree())

def Glattice_basis_reduce(s):

    t = g.degree()
    a = [] 
    a.append(0)
    b = []
    b.append(0);
    (q,r) = g.quo_rem(s)
    (a[0],b[0]) = simplify((g - q*s, 0 - q))

    if Gnorm(a[0],b[0]) > 2^t:
        a.append(0)
        b.append(0)
        (q,r) = s.quo_rem(a[0]);
        (a[1],b[1]) = (r, 1 - q*b[0])
        if a[1] == 0:
                return (s,1)

    else:
        return (a[0], b[0])
    i = 1
    while Gnorm(a[i],b[i]) > 2^t:
        a.append(0) 
        b.append(0)
        (q,r) = a[i-1].quo_rem(a[i])
        (a[i+1],b[i+1]) = (r, b[i-1] - q*b[i])
        i+=1
    return (a[i],b[i]);




def Patterson(g,goppa_code,word):
    t = g.degree()
    F2 = GF(2)
    F_2m = g.base_ring()
    Z = F_2m.gen()
    PR_F_2m = g.parent() 
    X = PR_F_2m.gen()

    factor_list = list(factor(2^m-1))
    final_factor_list = []
    for i in range(len(factor_list)):
        for j in range(factor_list[i][1]):
            final_factor_list.append(factor_list[i][0])

    while 1:
        primitive_root = F_2m.random_element()
        if primitive_root == 0:
            continue;

        for i in range(len(final_factor_list)):
            for j in itertools.combinations(final_factor_list,i):
                exponent = 1
                for _ in range(len(j)):
                    exponent *= j[_]
                if primitive_root^(exponent) == 1:
                    output = False
                    break
                else:
                    output = True
                    continue
            if output == False:
                break
        if output == True:
            break


    codelocators = []
    for i in range(2^m-1):
        codelocators.append(primitive_root^(i+1));
    codelocators.append(F_2m(0))

    SyndromeCalculator = matrix(PR_F_2m, 1, len(codelocators))
    for i in range(len(codelocators)):
        SyndromeCalculator[0,i] = inverse_mod((X - codelocators[i]),g)
    word = matrix(GF(2),word )


    X = g.parent().gen()
    synd = SyndromeCalculator*word.transpose()

    syndrome_poly = 0;
    for i in range (synd.nrows()):
        syndrome_poly += synd[i,0]*X^i
        
    error = matrix(GF(2),1,goppa_code.parity_check_matrix().ncols())
    (g0,g1) = Gsplit(g)

    (d,u,v) = xgcd(g1,g);
    sqrt_X = g0* u.mod(g)
    T = syndrome_poly.inverse_mod(g)

    (T0,T1) = Gsplit(T - X)
    R = (T0+ sqrt_X*T1).mod(g)
    (alpha, beta) = Glattice_basis_reduce(R)
    sigma = (alpha*alpha) + (beta*beta)*X
    if (X^(2^m)).mod(sigma) != X:
        error = error

    for i in range(len(codelocators)):
        if sigma(codelocators[i]) == 0:
            error[0,i] = 1
    return error


m = 5
F = GF(2^m)
R.<x> = F[]
g = x^3 + x^2 + 1
L = [a for a in F.list() if g(a) != 0]
C = codes.GoppaCode(g, L)
E = codes.encoders.GoppaCodeEncoder(C)
n = 32
k = 17
P = matrix(GF(2),n,[random.random()<0.5 for _ in range(n^2)])
while (rank(P) < n) :
    P[floor(n*random.random()),floor(n*random.random())] +=1 
S = random_nonsingular_matrix(k)
G = C.generator_matrix()
Gfin = S * G *P
Gfin = matrix(GF(2),Gfin)
worde = [0]*k
worde[6] = 1
worde = vector(GF(2),worde)
e = [0]*n
e[20] = 1 
e = vector(GF(2),e)
cip = ( worde * Gfin) + e
inS = S.inverse()
inP = P.inverse() 
cdec = cip * inP

err = Patterson(g,C,cdec)
mS=[0]*n
for i in range(0,n):
    mS[i] = (cdec[i] + err[0][i] ) % 2
    
mS = vector(GF(2),mS)
message = (S*G).solve_left(mS)
