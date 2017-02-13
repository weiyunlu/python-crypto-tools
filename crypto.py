from math import sqrt

def discrete_log(g, h, p):
    x = 0
    s = 1
    for i in range(1,p):
        x += 1
        s *= g
        s = s % p
        if s == h:
            return x
    return -1
        
def powermod(a, b, n):
    return (a**b) % n
    
def multimod(a, b, n):
    return (a*b) % n
    
def findordermod(x, n):
    for i in range(1, n):
        if powermod(x, i, n) == 1:
            return i
    return -1
    
def factor(n):     # attempt to factor a product of primes
    for i in range(2, int(sqrt(n))):
        if n % i == 0:
            return (i, n/i)
    return -1
    
def gcd(a,b):           # extended Euclidean algorithm
    u = 1
    g = a
    x = 0
    y = b
    q = 0
    s = 0
    t = 0
    
    while y != 0:
        q = g / y
        t = g % y
        s = u - (q * x)
        u = x
        g = y
        x = s
        y = t
        
    v = (g - a*u) / b
    
    if u > 0:
        while u - b/g > 0:
            u = u - b/g
            v = v - a/g
    else:
        while u < 0:
            u = u + b/g
            v = v - a/g
    
    return (g,u,v)          # g is the gcd, u and v satisfy au + bv = 1.  
                            # can use to find inverses in Z_p with gcd(a,p)[1]

def shanks(g, h, p):
    N = findordermod(g, p)
    n = 1 + int(sqrt(N))
    List_1 = []
    List_2 = []
    
    for i in range(n):
        List_1.append(powermod(g, i, p))

    u = powermod(gcd(g,p)[1], n, p)
    
    for i in range(n):
        List_2.append(multimod(h, powermod(u, i, p), p))
    
    collide_i = 0
    collide_j = 0
    
    for i in range(n):
        for j in range(n):
            if List_1[i] == List_2[j]:
                collide_i = i
                collide_j = j
                break
            
    return collide_i + collide_j * n
    
def miller_rabin(n, a):                     # Miller-Rabin test for primality of n
    if n % 2 == 0 or 1 < gcd(a, n)[0] < n:
        return "Composite."
    k = 0
    q = n - 1
    while q % 2 == 0:
        q = q/2
        k += 1
    a = powermod(a,q,n)
    if a % n == 1:
        return "Test fails."
    for i in range(k):
        if a % n == n-1:
            return "Test fails."
        a = powermod(a,2,n)
    return "Composite."
    
def pollard(N, a=2, b=100):
    for j in range(2, b):
        a = a**j % N
        d = gcd(a - 1, N)[0]
        if 1 < d < N:
            return d
    return -1
    
def quadres(a, p):
    a_p = a % p
    for i in range(p):
        if i**2 % p == a_p:
            return 1
    return -1
    
def find_roots(a, p):
    a_p = a % p
    roots = []
    for i in range(p):
        value = i**2 % p
        if value == a_p:
            roots.append(i)
    return roots
    
def bday_same_me(n):    # probability that at least one person in room of n people has same birthday as me
    return 1 - (float(364)/365)**n
    
def bday_same_anytwo(n):    # probability that -any- two people in room of n share a birthday
    product = 1
    for i in range(1, n+1):
        factor = float(365 - (i - 1)) / 365
        product *= factor
    return 1 - product
    
def same_card_twolists(N, n, m): # have N cards, mark n, pick m at random; chance of at least 1 match?
    product = 1
    for i in range(1, m+1):
        factor = float(N-n-i+1)/(N-i+1)
        product *= factor
    return 1 - product
    
class elliptic_curve:
    def __init__(self, A, B):
        eps = 1E-12
        valid = abs(4*(A**3) + 27*(B**2)) > eps
        msg = "Invalid curve!  This one doesn't have three distinct roots!"
        assert valid, msg
        self.A = A
        self.B = B
        
    def plus(self, P, Q):
        A = self.A
        B = self.B
        eps = 1E-12
        if P == "O":
            return Q
        if Q == "O":
            return P
        x1 = P[0]; y1 = P[1]
        x2 = Q[0]; y2 = Q[1]
        P_on_E = abs(y1**2 - x1**3 - A*x1 - B) < eps
        Q_on_E = abs(y2**2 - x2**3 - A*x2 - B) < eps
        msg_P = "Invalid point!  P does not lie on E!"
        msg_Q = "Invalid point!  Q does not lie on E!"
        assert P_on_E, msg_P
        assert Q_on_E, msg_Q
        if (x1 == x2 and y1 == -y2):
            return "O"
        else:
            if P == Q:
                slope = float(3*x1**2 + A)/(2*y1)
            else:
                slope = float(y2 - y1)/(x2 - x1)
        x3 = slope**2 - x1 - x2
        y3 = slope*(x1 - x3) - y1
        return (x3, y3)
        
class elliptic_curve_ff:
    def __init__(self, A, B, p):
        eps = 1E-12
        valid = abs(4*(A**3) + 27*(B**2)) > eps
        msg = "Invalid curve!  This one doesn't have three distinct roots!"
        assert valid, msg
        self.A = A
        self.B = B
        self.p = p
        
    def find_elements(self):
        A = self.A
        B = self.B
        p = self.p
        elements = ['O']
        for x in range(p):
            y_squared = x**3 + A*x + B
            roots = find_roots(y_squared, p)
            for y in roots:
                elements.append((x,y))
        return elements
        
    def count(self):
        elements = self.find_elements()
        return len(elements)
        
    def check_point(self,P):
        A, B, p = self.A, self.B, self.p
        x = P[0]
        y = P[1]
        success = (y**2 - x**3 - A*x - B) % p == 0
        msg = "Invalid input!  You entered a point not on the curve!"
        assert success, msg
        
    def plus(self, P, Q):
        A = self.A
        p = self.p
        if P == "O":
            return Q
        if Q == "O":
            return P
        x1 = P[0]; y1 = P[1]
        x2 = Q[0]; y2 = Q[1]
        self.check_point(P)
        self.check_point(Q)
        if (x1 == x2 and y1 % p == -y2 % p):
            return "O"
        else:
            if P == Q:
                numerator = (3*x1**2 + A)
                denominator_inv = gcd(2*y1, p)[1]
            else:
                numerator = y2 - y1
                denominator_inv = gcd(x2 - x1, p)[1]
        slope = numerator * denominator_inv
        x3 = (slope**2 - x1 - x2) % p
        y3 = (slope*(x1 - x3) - y1) % p
        return (x3, y3)
        
    def mult(self, n, P):
        self.check_point(P)
        Q = P
        R = 'O'
        while n > 0:
            if n % 2 == 1:
                R = self.plus(R,Q)
            Q = self.plus(Q,Q)
            n = int(n/2)
        return R
        
    def log(self, P, Q):
        self.check_point(P)
        self.check_point(Q)
        if Q == 'O':
            return 0
        order = self.count()
        n = 0
        result = 'O'
        while n < order:
            result = self.plus(result, P)
            n += 1
            if result == Q:
                return n
        return -1
        
    def table(self):
        elements = self.find_elements()
        print ' '*10 + '|',
        for x in elements:
            print '%10s' %(x,),
        print ''
        print '-'*10 + '|',
        for x in elements:
            print '-'*10,
        print ''
        for y in elements:
            print '%10s' %(y,) + '|',
            for x in elements:
                print '%10s' %(self.plus(x,y),),
            print ''