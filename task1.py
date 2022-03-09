import numpy as np
import time


class OnionMatrix:
    def __init__(self, blocks):
        self.n = len(blocks)
        self.m = blocks[0].shape[0]
        self.blocks = blocks
 
    # A1 cebulowa dolna trójkątna bez przekątnej taka że A1 + A2 = A
    def multiplyA1(self, vArray):
        resultArray = np.zeros((self.n, self.m))
        
        for index in range (1, self.n):
            resultArray[index] = resultArray[index-1] + self.blocks[index-1]@vArray[index-1]
          
        return np.hstack(resultArray)
        
    
    # A2 cebulowa górna trójkątna z przekątną taka że A1 + A2 = A
    def multiplyA2(self, vArray):
        vSum = np.zeros(self.m)
        resultArray = np.zeros((self.n, self.m))
        
        for index in range (self.n-1, -1, -1):
            vSum += vArray[index]
            resultArray[index] = self.blocks[index]@vSum
            
        return np.hstack(resultArray)

    
    def multiply(self, v):
        vArray = np.split(v, self.n)
        return self.multiplyA1(vArray) + self.multiplyA2(vArray)

    
    def solve(self, b):
        bArray = np.split(b, self.n)
        
        diagonal = np.zeros((self.n, self.m, self.m))
        # odejmujemy i-1 cebulowy wiersz od i-tego, żeby uzyskać macierz cebulową górną trójkątną
        for index in range (self.n-1, 0, -1):
            diagonal[index] = self.blocks[index] - self.blocks[index-1]
            bArray[index] = bArray[index] - bArray[index-1]
        diagonal[0] = self.blocks[0]
        
        xArray = np.zeros((self.n, self.m))
        xSum = np.zeros(self.m)
        for index in range (self.n-1, -1, -1):
            x_i = np.linalg.solve(diagonal[index], bArray[index] - diagonal[index]@xSum)
            xArray[index] = x_i
            
            xSum += xArray[index]
            
        return np.hstack(xArray)

