import sys
import os
import re
import math
from optparse import OptionParser
from math import ceil, log
import fileinput

from mrjob.job import MRJob
from mrjob.job import MRStep 

def getDimensions(fileName):
    with open(fileName, 'r+') as f: 
        line = f.readline() 
        data = f.read()
        f.seek(0)
        f.write(data) 
        f.truncate() 
        f.close()
    r, c = line.split()
    
    return int(r), int(c)
def writeBackDimensions(fileName, a,b):
    with open(fileName, 'r+') as f: 
        line = str(a)+" "+ str(b)+ "\r\n"
        data = f.read()
        f.seek(0)
        f.write(line)
        f.write(data) 
        f.truncate()
        f.close() 
    
def matrix_product(A, B):
    n = len(A)
    C = [[0 for i in range(n)] for j in range(n)]
    for i in range(n):
        for k in range(n):
            for j in range(n):
                C[i][j] += A[i][k] * B[k][j]
    return C

def add(A, B):
    n = len(A)
    C = [[0 for j in range(0, n)] for i in range(0, n)]
    for i in range(0, n):
        for j in range(0, n):
            C[i][j] = A[i][j] + B[i][j]
    return C


def subtract(A, B):
    n = len(A)
    C = [[0 for j in range(0, n)] for i in range(0, n)]
    for i in range(0, n):
        for j in range(0, n):
            C[i][j] = A[i][j] - B[i][j]
    return C


def recursiveStrassen(A, B):
    N = len(A)
    N = int(N)
   
    if N == 1:
        return matrix_product(A, B)
        
    else:
        # Create sub-matrices of zeros
       
        n = math.floor(N/2)
        
        a11 = [[0 for j in range(0, n)] for i in range(0, n)]
        a12 = [[0 for j in range(0, n)] for i in range(0, n)]
        a21 = [[0 for j in range(0, n)] for i in range(0, n)]
        a22 = [[0 for j in range(0, n)] for i in range(0, n)]

        b11 = [[0 for j in range(0, n)] for i in range(0, n)]
        b12 = [[0 for j in range(0, n)] for i in range(0, n)]
        b21 = [[0 for j in range(0, n)] for i in range(0, n)]
        b22 = [[0 for j in range(0, n)] for i in range(0, n)]

        aResult = [[0 for j in range(0, n)] for i in range(0, n)]
        bResult = [[0 for j in range(0, n)] for i in range(0, n)]

        # Dividing the matrices in 4 sub-matrices:
        #---------
        # Q1 | Q2
        #---------
        # Q3 | Q4
        #---------
        for i in range(0, n):
            for j in range(0, n):
                a11[i][j] = A[i][j]         # Q1
                a12[i][j] = A[i][j + n]     # Q2
                a21[i][j] = A[i + n][j]     # Q3
                a22[i][j] = A[i + n][j + n] # Q4

                b11[i][j] = B[i][j]         # Q1
                b12[i][j] = B[i][j + n]     # Q2
                b21[i][j] = B[i + n][j]     # Q3
                b22[i][j] = B[i + n][j + n] # Q4

        # Calculating p1 to p7:
        '''
        p1 = (a11+a22) * (b11+b22)
        p2 = (a21+a22) * (b11)
        p3 = (a11) * (b12 - b22)
        p4 = (a22) * (b21 - b11)
        p5 = (a11+a12) * (b22)
        p6 = (a21-a11) * (b11+b12)
        p7 = (a12-a22) * (b21+b22)
        '''
        aResult = add(a11, a22)
        bResult = add(b11, b22)
        p1 = recursiveStrassen(aResult, bResult)

        aResult = add(a21, a22)              
        p2 = recursiveStrassen(aResult, b11)

        bResult = subtract(b12, b22) 
        p3 = recursiveStrassen(a11, bResult)  

        bResult = subtract(b21, b11) 
        p4 =recursiveStrassen(a22, bResult)   

        aResult = add(a11, a12)     
        p5 = recursiveStrassen(aResult, b22)  

        aResult = subtract(a21, a11) 
        bResult = add(b11, b12)     
        p6 = recursiveStrassen(aResult, bResult) 

        aResult = subtract(a12, a22) 
        bResult = add(b21, b22)      
        p7 = recursiveStrassen(aResult, bResult)

        c12 = add(p3, p5)
        c21 = add(p2, p4)  

        aResult = add(p1, p4) 
        bResult = add(aResult, p7) 
        c11 = subtract(bResult, p5) 

        aResult = add(p1, p3) 
        bResult = add(aResult, p6)
        c22 = subtract(bResult, p2) 

        # C = C11+C12+C21+C22
        C = [[0 for j in range(0, N)] for i in range(0, N)]
        for i in range(0, n):
            for j in range(0, n):
                C[i][j] = c11[i][j]
                C[i][j + n] = c12[i][j]
                C[i + n][j] = c21[i][j]
                C[i + n][j + n] = c22[i][j]
        return C
        
def strassen(A, B, l):
    assert type(A) == list and type(B) == list
    assert len(A) == len(A[0]) == len(B) == len(B[0])

    # Make the matrices bigger so that you can apply the strassen
    # algorithm recursively without having to deal with odd
    # matrix sizes
    n = l
    m = PowerOf2(n)
    APrep = [[0 for i in range(m)] for j in range(m)]
    BPrep = [[0 for i in range(m)] for j in range(m)]
    for i in range(l):
        for j in range(l):
            APrep[i][j] = A[i][j]
            BPrep[i][j] = B[i][j]
    CPrep = recursiveStrassen(APrep, BPrep)
    C = [[0 for i in range(n)] for j in range(n)]
    for i in range(n):
        for j in range(n):
            C[i][j] = CPrep[i][j]
    return C

def PowerOf2(n):
    #print (2**int(ceil(log(n,2))))
    return int(2**int(ceil(log(n,2))))


class MR_Matrix_Multiplication_B(MRJob):
    
    def mapperOne(self,_,line): # id
        #print ("ola1")
        i, j, value = line.split()	
       	filename = os.environ['map_input_file']
        l = max(m,n,p)
        L = int(PowerOf2(int(l)))-1
        i = int(i)
        j = int(j)
        value =float(value)
      
        # splitting up the matrices into 4 sub-matrices
        if filename == sys.argv[1]:
            #quad 1 A11
            if i < L/2 and j <L/2:
                #print ("A1")
                yield 'A11B11', ('a', int(i),int(j), float(value)) # C
                yield 'A11B12', ('a', int(i),int(j), float(value))
            #quad 3 A21
            elif i > L/2 and j <L/2:
                #print ("A2")
                yield 'A21B11', ('a', int(i),int(j), float(value))
                yield 'A21B12', ('a', int(i),int(j), float(value))
            #quad 2 A12
            elif i < L/2 and j >L/2:
                #print ("A3")
                yield 'A12B21', ('a', int(i),int(j), float(value)) # T
                yield 'A12B22', ('a', int(i),int(j), float(value))
            #quad 4 A22
            elif i > L/2 and j >L/2:
                #print ("A4")
                yield 'A22B21', ('a', int(i),int(j), float(value))
                yield 'A22B22', ('a', int(i),int(j), float(value))
                
        elif filename == sys.argv[2]:
             #quad 1 B11
             if i < L/2 and j <L/2:
                 yield 'A11B11', ('b', int(i),int(j), float(value)) 
                 yield 'A21B11', ('b', int(i),int(j), float(value))
             #quad 2 B12
             elif i < L/2 and j >L/2:
                 yield 'A11B12', ('b', int(i),int(j), float(value))
                 yield 'A21B12', ('b', int(i),int(j), float(value))
             #quad 3 B21
             elif i > L/2 and j <L/2:
                 yield 'A12B21', ('b', int(i),int(j), float(value)) 
                 yield 'A22B21', ('b', int(i),int(j), float(value))
             #quad 4 B22
             elif i > L/2 and j >L/2:
                 yield 'A12B22', ('b', int(i),int(j), float(value))
                 yield 'A22B22', ('b', int(i),int(j), float(value))
           
    def reducerOne(self, key, values): # Use Strassen Algorithm and output C and T matrices
        l = max(m,n,p)
        l = int(PowerOf2(int(l))/2)
        A = [[0 for col in range(l)] for row in range(l)]
        B = [[0 for col in range(l)] for row in range(l)]
        
        for v in values:
            if v[0] == 'a':
                A[v[1]%l][v[2]%l] = v[3]
               
            elif v[0] =='b':
                B[v[1]%l][v[2]%l] = v[3]

        C = strassen(A,B,l)
     
        identity = ""
        if key == 'A11B11': 
            identity = "C11"
        elif key == 'A11B12': 
            identity = "C12"
        elif key == 'A21B11': 
            identity = "C21"
        elif key == 'A21B12': 
            identity = "C22"
        
        if key == 'A12B21':
            identity = "T11"
        elif key == 'A12B22': 
            identity = "T12"
        elif key == 'A22B21':
            identity = "T21"
        elif key == 'A22B22': 
            identity = "T22"
        for a in range(len(C)):
            for b in range(len(C[a])):
                yield identity[0], [(identity,C[a][b],a,b)]
        
    def mapperTwo(self, key, values): # Map matrices C and T to indices (i,k)
        #print ("ola3")
        l = max(m,n,p)
        l = PowerOf2(l)/2
        i = 0
        k = 0
    
        for v in values:
           
            if v[0] == "C11" or v[0] == "T11" :
                i = v[2]
                k = v[3]
            elif v[0] == "C12" or v[0] == "T12":
                i = v[2]
                k = v[3]+l
            elif v[0] == "C21" or v[0] == "T21":
                i = v[2]+l
                k = v[3]+0
            elif v[0] == "C22" or v[0] == "T22":
                i = v[2]+l
                k = v[3]+l
            i= int(i)
            k= int(k)
            #print (v)
            
            yield (i,k), (float(v[1]))
    
    def reducerTwo(self, key, values): # sum matrices C and T
        #print ("ola4")
        ret = sum(values)
        if ret > 0:
            yield (key), ret
        
    def steps(self):
        return  [MRStep(mapper=self.mapperOne,reducer=self.reducerOne),
        MRStep(mapper=self.mapperTwo,reducer=self.reducerTwo)]
n=0
m=0
p=0

if __name__ == '__main__':
 
    A = sys.argv[1]
    B = sys.argv[2]
    n, m = getDimensions(A)
    m, p = getDimensions(B)
    MR_Matrix_Multiplication_B.run()
    writeBackDimensions(A,n,m)
    writeBackDimensions(B,m,p)
    writeBackDimensions("outC.list", n,p)
    