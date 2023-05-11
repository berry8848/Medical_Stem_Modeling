# 未完成

from modules import Point
from modules import Biharmonic
from modules import CrossingNumberAlgorithm
import numpy as np
import math
import random

class PDS(Biharmonic.Biharmonic):
    def __init__(self, points, cs, lambdas):
        #Biharmonic.__init__(self, points, cs, lambdas)
        self.biharmonic = Biharmonic.Biharmonic(points, cs, lambdas)
        # PDSでの点の生成範囲の設定
        self.x_max = points[0][1]
        self.x_min = points[0][1]
        self.y_max = points[0][2]
        self.y_min = points[0][2]
        self.z_max = points[0][3]
        self.z_min = points[0][3]

        for point in points:
            if point[1] > self.x_max:
                self.x_max = point[1]
            if point[1] < self.x_min:
                self.x_min = point[1]
            
            if point[2] > self.y_max:
                self.y_max = point[2]
            if point[2] < self.y_min:
                self.y_min = point[2]

            if point[3] > self.z_max:
                self.z_max = point[3]
            if point[3] < self.z_min:
                self.z_min = point[3]

        # PDSでの点の生成範囲の表示
        print("x_max = ", self.x_max, "x_min = ", self.x_min)
        print("y_max = ", self.y_max, "y_min = ", self.y_min)
        print("z_max = ", self.z_max, "z_min = ", self.z_min)
        # pds_point が含まれる四面体を divide_tetrahedron_list から捜索．点が連続でN回生成できなかったら終了
        self.N = 80
        self.num = 0
        self.fixed_points = [] # 確定点
        # CrossNumberAlgorithm
        self.CNA = CrossingNumberAlgorithm.CrossingNumberAlgorithm()
        

    def generate(self):
        while self.num < self.N:
            if len(self.fixed_points) >= 450:
                break

            flg_P = False
            while not flg_P: # 物体内部の点を生成するまでループ。内部であればflg_PはTrueとなる。
                # ランダムな点を生成
                pds_x = random.uniform(self.x_min, self.x_max)
                pds_y = random.uniform(self.y_min, self.y_max)
                pds_z = random.uniform(self.z_min, self.z_max)
                pds_point = [pds_x, pds_y, pds_z]
                # 物体内部の点か判定
                flg_P = self.CNA.cramer(pds_point)

            # 生成点の座標情報をPointに格納し候補点にする
            # candidate_point = Point(pds_point)
            # candidate_point.pds_coordinate()
            candidate_point = pds_point
            # Biharmonicより候補点の応力を導出
            # new_stress = self.biharmonic.cal(candidate_point.x, candidate_point.y, candidate_point.z)
            new_stress = self.biharmonic.cal(candidate_point[0], candidate_point[1], candidate_point[2])
            # 導出した応力より候補点の満たすべき密度を導出
            density = stress_to_density(new_stress)
            # 導出した密度より候補点の満たすべき点間距離を導出
            long = density_to_long(density)

            # 点間距離内に他の点が存在するかを確認
            flg = check_distance(self.fixed_points, candidate_point, long)

            # 点間距離内に他の点が存在しないとき候補点を確定点に追加
            if flg:
                self.fixed_points.append(pds_point)
                num = 0
                print('num : ', num, 'fixed_points : ', len(self.fixed_points))
            # 点間距離内に他の点が存在するとき
            else :
                num = num + 1
                print('num : ', num)
        return self.fixed_points



def stress_to_density(stress):
    if stress >= 0: #生データが正の場合，0を返す．（通常はマイナスの値をとる）
        density = 1.448
    else:
        #density = -stress/15
        density = 3*0.00001*(stress)*(stress) + 0.01*(stress) + 1.448
    return density

def density_to_long(density):
    long = 10 * 0.1 * math.sqrt(math.sqrt(2)*math.pi/ density)
    return long

def check_distance(fixed_points, candidate_point, long):
    # b_x = candidate_point.x 
    # b_y = candidate_point.y
    # b_z = candidate_point.z
    b_x = candidate_point[0] 
    b_y = candidate_point[1]
    b_z = candidate_point[2]

    check = True
    for a in fixed_points:
        a_x = a[0]
        a_y = a[1]
        a_z = a[2]
        distance = (a_x- b_x)**2 + (a_y- b_y)**2 + (a_z- b_z)**2
        if long**2 > distance:
            check = False
            break
    return check
    

