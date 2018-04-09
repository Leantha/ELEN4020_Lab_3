import sys
import os
import re

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
        
# ----- Algorithm A for matrix multiplication ----- ---- ---- 
# Algorithm A executes by mapping and reducing twice
# ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- 
# The first map: maps each element in matrix A(I,J,V) and B(J,K,W) to their key value pairs
# A(I,J,V) -> (J, (A,I,Aij)) and N(J,K,W) -> (J, (B,K, Bjk))
# The first reduce: For each key j, it examines its list of associated values
# and produces a tuple: (I,K,(V = Aij*Bjk))
# The output is: ( J, [ (I1,K1,V1),(I2,K2,V2), ... , (Ip,Kp,Vp) ] ) or a list of all the tuples
# The second map: produces p key-value pairs: ((I1,K1), V1), ((I2,K2), V2), ..., ((Ip,Kp), Vp))
# The second reduce: for each key (i,k), produce the sum of the list values associated with this key.
# It outputs ((I,K),V)
# V is the element in the row I and column K of the matrix C = AB
# An extra mapreduce job has to be initially run, in order to retrieve the values i j k.
# ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- 
class MR_Matrix_Multiplication_A(MRJob):
    #OUTPUT_PROTOCOL = MRJob.protocol.RawValueProtocol
    #def configure_options(self):
        #super(MR_Matrix_Multiplication_A, self).configure_options()
        #self.add_passthrough_option('--matrix1', help="File1ForLab3.txt")
        #self.add_passthrough_option('--matrix2', help="File3ForLab3.txt")
        #MRJob.mrjob_conf_option_name(output_dir(os.getcwd()))
        #self.add_file_option('--scoring-db', help = "output.txt")
    #OUTPUT_PROTOCOL= JSONProtocol
    def mapperOne(self,_,line): #emit values
        i, j, value = line.split()
        filename = os.environ['map_input_file']

        if filename == sys.argv[1]:
            yield j, ('a', i, value)
        elif filename == sys.argv[2]:
            yield i, ('b', j, value)

    def mapperTwo(self, j, value): # identify
        yield (j), value

    def reducerOne(self, j, values): #multiply values
        values_A = []
        values_B = []

        for key in values:
            if key[0] == 'a':
                values_A.append(key)
            elif key[0] =='b':
                values_B.append(key)

        for (a, i, v1) in values_A:
            for (b, k, v2) in values_B:
                yield (i,k), int(v1)*int(v2)

    def reducerTwo(self, key, value): #add values
        yield key, sum(value)

    def steps(self):
        return  [MRStep(mapper=self.mapperOne,reducer=self.reducerOne),
                 MRStep(mapper=self.mapperTwo,reducer=self.reducerTwo)]

if __name__ == '__main__':
    A = sys.argv[1]
    B = sys.argv[2]
    n, m = getDimensions(A)
    m, p = getDimensions(B)
    MR_Matrix_Multiplication_A.run()
    writeBackDimensions(A,n,m)
    writeBackDimensions(B,m,p)
    writeBackDimensions("outC.list", n,p)
    
