import NISTP256 as EC

import numpy as np
import random as rand
import os

leakage = []

""" useful functions
    some are based on the field F_p
"""

def HW(val):
    return bin(val).count('1')

def random(size):
    return rand.randint(0, 2**size-1)

def sqrt(a):
    b = pow(a, EC.p_m_3_d_4, EC.p)        
    x = (b * a) % EC.p     
    if(((b * x) % EC.p) == 1):
        return x
    else:
        return False      
        
def LIM( a, b ):
    
    res = [ 0 for i in range(16)]
    cst = 2**32
    
    for i in range(8):
        
        tmp = [0,0]
        c = 0 
        
        for j in range(8):
            
            m = a[7-j] * b[7-i]
            tmp[1] = (m % cst) + res[15 - i - j]
            tmp[0] = (m >> 32) % cst
            
            # WRITE MSB LEAKAGE INTO FILE WHICH IS THE HAMMING WEIGHT OF TMP[0]
            global leakage 
            leakage.append(HW(tmp[0]))
            
            if( tmp[1] >= cst ):
                tmp[0] += 1
                tmp[1] %= cst
            
            tmp[1] += c
            
            if( tmp[1] >= cst ):
                tmp[0] += 1
                tmp[1] %= cst
            
            res[15 - i - j] = tmp[1]
            c = tmp[0]
        
        res[7 - i] = c
    
    return res
            
def write_leakage(file_name):
    np.array(leakage).tofile(file_name)
        
def read_leakage(file_name):
    return np.fromfile(file_name, dtype=int)
    
#def read_leakage(file_name):
    #return np.fromfile(file_name, dtype=np.uint64)
        
""" class to handle field F_p element
    A field element is a long integer 
    If the field element is handled in binary or hex,
    it is in BIG ENDIAN -> msb to the left and lsb to the right. 
"""
class Field_element:
        
    def __init__(self, v = 0):
        self.value = v % EC.p
   
    def random():
        return Field_element( rand.randint(0, EC.p-1))
        
    def random_msb_set():
        ks = True
        value = 0
        while(ks):
            value = rand.randint(0, EC.p-1)
            value += 2**255
            if( value < EC.p):
                ks = False
        return Field_element(value)
    
    def switch_bit(self, idx):
        mask = 1 << (255 - idx) 
        return Field_element(self.value ^ mask)
        
    def to_regrep(self):
        tbin = bin(self.value)[2:]
        tbin = tbin.rjust(256, '0')
        tab = [ int( tbin[32*i : 32*(i+1)] , 2 ) for i in range(8)]
        return tab
    
    def from_regrep(tab):
        sum = 0
        l = len(tab)
        c = 2**32
        for i in range(l):
            sum += tab[l - i - 1]*(c**i)
        return Field_element(sum)
        

    """ print the field element
        base : 'b' for binary 
               'd' for decimal
               'x' for hexadecimal (default) 
               'o' for octal
    """    
    def print(self, base = 'x'):
        print(format(self.value, base))
    
    def size_bit(self):
        return (self.value.bit_length())
    
    def size_byte(self):
        return (math.ceil(self.value.bit_length()/8))

    def sqrt(self):
        return(sqrt(self.value))
    
    def __neg__(self):
        return Field_element(-self.value % EC.p)
    
    def __add__(self, other):
        return Field_element(( self.value + other.value) % EC.p)
    
    def __radd__(self, other):
        return Field_element( (self.value + other) % EC.p)
    
    def __sub__(self, other):
        return Field_element(( self.value - other.value) % EC.p)
    
    def __mul__(self, other):
        return Field_element(( self.value * other.value) % EC.p)
            
    def __pow__(self, other):
        return Field_element( pow(self.value, other.value , EC.p))
    
    def __eq__(self, other):
        return (self.value == other.value)
    
    def LIM(self, other):
        return Field_element.from_regrep(LIM( self.to_regrep(), other.to_regrep()) )
    
    
    
""" class to handle a point on E(F_p)
"""    
    
class Point:
    
    def __init__(self, x , y , z):
        self.x = Field_element(x) 
        self.y = Field_element(y) 
        self.z = Field_element(z)
        
    def infinity():
        return Point(1,1,0)
    
    def is_infinity(self):
        return ( self.z.value == 0 )
    
    def opposite(self):
        return( Point( self.x.value , -self.y.value % EC.p , self.z.value))
    
    def equal(self, Q):
        return( (Q.x * self.z**Field_element(2) == self.x * Q.z**Field_element(2))
                &(Q.y * self.z**Field_element(3)  == self.y * Q.z**Field_element(3)))
    
    def is_opposite(self, Q):
        return( (Q.x * self.z**Field_element(2) == self.x * Q.z**Field_element(2))
                &(Q.y * self.z**Field_element(3) == -self.y * Q.z**Field_element(3)) )
        
    def print(self):
        print(str(self.x.value) + ' , ' + str(self.y.value) + ' , ' + str(self.z.value) )
        
    def write(self, file_name):
        with open(file_name, "w") as f:
            f.write(hex(self.x.value)+'\n')
            f.write(hex(self.y.value)+'\n')
            f.write(hex(self.z.value)+'\n')
            f.close()
            
    def read(file_name):
        with open(file_name, "r") as f:
            coo = f.readlines()
            f.close()
        return Point( int(coo[0], 16) , int(coo[1], 16), int(coo[2], 16))
        
    
    def generate(randz = False):
        while(1):
            x = rand.randint(0, EC.p-1)
            y2 = ((x**3) + (EC.a*x) + EC.b) % EC.p 
            y = sqrt(y2)
            if y != False:
                z = 1 
                if randz:
                    z = rand.randint(0, EC.p-1)
                    x = (z**2) * x % EC.p
                    y = (z**3) * y % EC.p
                return Point(x,y,z)
    
    def is_on_curve(self):
        eq = self.x.value**3 
        eq += EC.a * self.x.value * (self.z.value**4)
        eq += EC.b * (self.z.value**6) 
        eq %= EC.p 
        Y2 = self.y.value**2 
        Y2 %= EC.p 
        return (eq == Y2)
    
            
    def add(self, other):
        
        if (self.is_infinity()) : return other 
        if (other.is_infinity()) : return self 
        if (other.is_opposite(self)) : return Point.infinity()
        
        A = self.z.LIM(self.z)      # self.z * self.z 
        B = other.z.LIM(other.z)    # other.z * other.z
        C = self.x.LIM(B)           # self.x * B 
        D = other.x.LIM(A)          # other.x * A
        E = C - D
        F = self.y.LIM(B)           # self.y * B 
        F = F.LIM(other.z)          # F * other.z
        G = other.y.LIM(A)          # other.y * A 
        G = G.LIM(self.z)           # G * self.z
        H = F - G 
        I = E.LIM(E)                # E * E  
        J = I.LIM(E)                # I * E
        K = C.LIM(I)                # C * I 
        
        X = H.LIM(H) + J - K - K    # H * H
        Y = (K - X)
        Y = H.LIM(Y)                # H * Y
        tmp = F.LIM(J)              # F * J
        Y = Y - tmp 
        Z = self.z.LIM(other.z)     # self.z * other.z 
        Z = Z.LIM(E)                # Z * E 
        
        return Point(X.value ,Y.value, Z.value)
        
    
    def doubl(self):
        
        if (self.is_infinity()) : return Point.infinity()
        if (self.y.value == 0) : return Point.infinity()
        
        A = self.x.LIM(self.x)           # self.x * self.x
        C = self.z.LIM(self.z)           # self.z * self.z 
        D = A + A + A
        tmp = C.LIM(C)                   # C * C
        tmp = Field_element(-3*tmp.value)
        D = D + tmp
        B = self.y.LIM(self.y)           # self.y * self.y
        E = B.LIM(B)                     # B * B 
        F = self.x.LIM(B)                # self.x * B 
        F = F + F
        F = F + F
        
        X = (D.LIM(D)) - F - F           # D * D
        Z = self.y.LIM(self.z)           # self.y * self.z
        Z = Z + Z
        Y = F - X 
        Y = D.LIM(Y)                     # D * Y
        tmp = E + E
        tmp = tmp + tmp 
        tmp = tmp + tmp
        Y = Y - tmp  
        
        return Point(X.value, Y.value, Z.value)
    
    def montgomery_ladder(self, k):
        
        kbin = bin(k.value)[2:]
        
        R0 = self
        R1 = self.doubl()
        
        global leakage 
        leakage = []
        
        for i in range(1, len(kbin)):
            if (kbin[i] == str(0)) :
                R1 = R1.add(R0)
                R0 = R0.doubl()
            elif (kbin[i] == str(1)) :
                R0 = R0.add(R1)
                R1 = R1.doubl()
                
        return R0
   

dir = './test'
os.makedirs(dir)

for i in range(1000):
    curr_dir = dir + '/point' + str(i)
    os.makedirs(curr_dir)
    
    P = Point.generate(True)
    P.write(curr_dir + '/point')
    
    leakage = []
    Point.montgomery_ladder(P, Field_element(2))
    write_leakage( curr_dir + '/trc0')
    
    leakage = []
    Point.montgomery_ladder(P, Field_element(3))
    write_leakage( curr_dir + '/trc1')
      

l = read_leakage("./test/point0/trc0")

print(l[:100])

