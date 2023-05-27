import numpy as np

#Biharmonicを用いたRBF補間
class Biharmonic:
    def __init__(self, points, cs, lambdas):
        self.N = len(points)
        self.points = points
        self.cs = cs
        self.lambdas = lambdas

    def p1(self, x, y, z):
        return self.cs[0] + self.cs[1]* x + self.cs[2]*y + self.cs[3]*z

    def p2(self, x, y, z):
        a = 0
        for i in range(self.N):
            a += self.lambdas[i]*np.linalg.norm([x, y, z] - self.points[i, 1:4])
        return a

    def cal(self, x, y, z):
        return self.p1(x, y, z) + self.p2(x, y, z)
