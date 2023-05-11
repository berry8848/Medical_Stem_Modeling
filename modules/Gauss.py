#モジュールインポート
import numpy as np

class Gauss:
    def __init__(self, points):
        #前提条件の入力
        self.N = len(points)
        self.A = []
        print('N : ' , self.N, 'A: ',self.A)
        for i in range(self.N):
            A_ = []
            for j in range(self.N):
                phi = np.linalg.norm(points[i, 1:4] - points[j, 1:4])
                A_.append(phi)
            A_.extend([1, points[i, 1], points[i, 2], points[i, 3]])
            #print('A_i: ', i)
            self.A.append(A_)

        # [1,...,1, 0, 0, 0, 0] の作成
        A_ = [1]*self.N
        A2_ = [0]*4
        A_.extend(A2_)
        self.A.append(A_)
        # P, 0行列の作成
        x_ = []
        y_ = []
        z_ = []
        for k in range(self.N):
            x_.append(points[k, 1])
            y_.append(points[k, 2])
            z_.append(points[k, 3])
        x2_ = [0]*4
        x_.extend(x2_)
        self.A.append(x_)
        y2_ = [0]*4
        y_.extend(y2_)
        self.A.append(y_)
        z2_ = [0]*4
        z_.extend(z2_)
        self.A.append(z_)
        print('A: ', np.shape(self.A))

        self.B=np.append(points[:,4], [0,0,0,0])
        print('B = ', len(self.B))
    
    #部分ピポット選択
    def select_pivot(self, xs, i):
        k = np.abs(xs[i:,i]).argmax() + i
        if k != i:
            temp = xs[i].copy()
            xs[i] = xs[k]
            xs[k] = temp

    #Gaussの消去法
    def gauss1(self):
        self.A = np.array(self.A)
        # 拡大係数行列の生成
        n = len(self.A)
        zs = np.c_[self.A.astype(np.float_), self.B.astype(np.float_)]
        # 前進消去
        for i in range(n - 1):
            self.select_pivot(zs, i)
            for j in range(i + 1, n):
                temp = zs[j, i] / zs[i, i]
                zs[j, i:] -= temp * zs[i, i:]
        print(zs)
        # 後退代入
        for i in range(n - 1, -1, -1):
            zs[i, n] -= zs[i, i+1:n] @ zs[i+1:, n]
            zs[i, n] /= zs[i, i]
        print('zs: ', len(zs))
        return zs[:, n], zs[:self.N, n], zs[self.N:, n]


