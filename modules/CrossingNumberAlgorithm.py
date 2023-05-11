import numpy as np
import random
import math
from modules import Point
from modules import Biharmonic

class CrossingNumberAlgorithm:
    def __init__(self):
       

        self.lefts = np.array([[None, None, None]])
        #self.origin = np.array([(x_max+x_min)/2, (y_max+y_min)/2, (z_max+z_min)/2])
        # 表面形状データの読み込み
        print("読み込み中...")
        f = open("Input/SYS.txt")
        self.list1 = []
        self.list2 = []
        count = 1
        for line in f.readlines():
            # vertexとfaceに分けてlist化
            if count < 3313: # vertex 形状データによって値を変更
                self.list1.append(line.split())
            else: # list
                self.list2.append(line.split())
            count+=1
            # strをfloatに変換
            self.list1 = [[float(x) for x in y] for y in self.list1]
            self.list2 = [[int(x) for x in y] for y in self.list2]
        print('list1 = ',self.list1[:5])
        print('list2 = ',self.list2[:5])
        print("読み込み終了")

        #surface用にnp変換
        self.np_list1 = np.array(self.list1)
        self.np_list2 = np.array(self.list2)
        self.np_list1 = self.np_list1[:,0:3] #[x, y, z]座標
        self.np_list2 = self.np_list2[:,1:4] #[点番号１, 点番号2, 点番号3] ex[0, 5, 2]
    
    def biharmonic(self, points, cs, lambdas):
        # Biharmonic
        self.biharmonic = Biharmonic.Biharmonic(points, cs, lambdas)
        print(self.biharmonic.cal(15.78, 4.09, 52.35))

    # PDS点群に表面点群を連結
    def outline(self, outline_points):
        for list in self.list1:
            array = [list[0], list[1], list[2]]
            outline_points.append(array)

    #物体表面上に点を2個ずつ配置     
    def surface(self, fixed_points):
        for list2 in self.np_list2:
            # 三角形メッシュの中心に点を生成(1つ目の点)
            pds_point = (self.np_list1[list2[0]]+self.np_list1[list2[1]]+self.np_list1[list2[2]])/3
            fixed_points.append(pds_point)

            #三角形メッシュ内に任意に点を生成（2つ目の点）
            # # ベクトルの生成
            # a = self.np_list1[list2[1]] - self.np_list1[list2[0]] 
            # b = self.np_list1[list2[2]] - self.np_list1[list2[0]]
            # s = random.random()
            # t = random.random()
            # while s+t>1:
            #     s = random.random()
            #     t = random.random()
            # pds_point = self.np_list1[list2[0]] + s*a + t*b #新点の生成
            # fixed_points.append(pds_point)
    

    #物体表面上でPDS       
    def surface_pds(self, fixed_points):
        N = 800
        num = 0
        #表面形状のデータ整理
        # self.np_list1 = np.array(self.list1)
        # self.np_list2 = np.array(self.list2)
        # self.np_list1 = self.np_list1[:,0:3] #[x, y, z]座標
        # self.np_list2 = self.np_list2[:,1:4] #[点番号１, 点番号2, 点番号3] ex[0, 5, 2]

        # print('np_list1 = ', self.np_list1)
        # print('np_list2 = ', self.np_list2)

        while num < N:
            if len(fixed_points) >= 1000: #4to6のPDSの点の数も含む
                break 

            # 候補点の生成
            i = random.randint(0, len(self.np_list2)-1)
            list2 = self.np_list2[i]
            # ベクトルの生成
            a = self.np_list1[list2[1]] - self.np_list1[list2[0]] 
            b = self.np_list1[list2[2]] - self.np_list1[list2[0]]
            #print(self.np_list1[list2[1]])
            s = random.random()
            t = random.random()
            while s+t>1:
                s = random.random()
                t = random.random()
            pds_point = self.np_list1[list2[0]] + s*a + t*b #新点の生成
            
            # 生成点の座標情報をPointに格納し候補点にする
            candidate_point = Point.Point(pds_point)
            candidate_point.pds_coordinate()
            # Biharmonicより候補点の応力を導出
            new_stress = self.biharmonic.cal(candidate_point.x, candidate_point.y, candidate_point.z)
            # 導出した応力より候補点の満たすべき密度を導出
            density = stress_to_density(new_stress)
            # 導出した密度より候補点の満たすべき点間距離を導出
            long = density_to_long(density)
            # 点間距離内に他の点が存在するかを確認
            flg = check_distance(fixed_points, candidate_point, long)

            # 点間距離内に他の点が存在しないとき候補点を確定点に追加
            if flg:
                fixed_points.append(pds_point)
                num = 0
                print('num : ', num, 'fixed_points : ', len(fixed_points))
            # 点間距離内に他の点が存在するとき
            else :
                num = num + 1
                print('num : ', num)

            

    def cramer(self, pds_point): #クラメルの公式を用いて内外判定
        #ray = self.origin - pds_point
        ray = np.array([1, 0.1, 0.1]) - pds_point
        cross_num = 0
        for l in self.list2:
            v0 = np.array([self.list1[l[1]][0], self.list1[l[1]][1], self.list1[l[1]][2]]) # [v0_x, v0_y, v0_z]の順
            v1 = np.array([self.list1[l[2]][0], self.list1[l[2]][1], self.list1[l[2]][2]])
            v2 = np.array([self.list1[l[3]][0], self.list1[l[3]][1], self.list1[l[3]][2]])

            OA = v1 - v0
            OB = v2 - v0

            left = [[OA[0], OB[0], ray[0]], 
                    [OA[1], OB[1], ray[1]], 
                    [OA[2], OB[2], ray[2]]]
            
            right = pds_point - v0

            u, v, t = np.linalg.solve(left, right)

            if t>=0 and u>=0 and u<=1 and v>=0 and v<=1 and u+v>=0 and u+v<=1:
                cross_num += 1 #三角形の平面内で交点を持つ → カウント
                #print("u, v, t = ", u, v, t)
            else:
                continue #三角形の平面外で交点を持つ
        
        #print("cross_num", cross_num)
        if cross_num % 2 == 0:
            flg = False #外側
        else:
            flg = True #内側
        return flg


# 応力と密度の関係式．PDSのみに使用　※関係式が微妙なため，臨時で別の関数
def stress_to_density(stress):
    if stress >= 0: #生データが正の場合，0を返す．（通常はマイナスの値をとる）
        density = 1.448
    else:
        #density = -stress/15
        density = 3*0.00001*(stress)*(stress) + 0.01*(stress) + 1.448
    return density

def density_to_long(density):
    long = 6 * 0.1 * math.sqrt(math.sqrt(2)*math.pi/ density)
    return 1*long

def check_distance(fixed_points, candidate_point, long):
    b_x = candidate_point.x 
    b_y = candidate_point.y
    b_z = candidate_point.z
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
