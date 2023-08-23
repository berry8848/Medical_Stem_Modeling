import numpy as np
import random
import math

#交差数判定法
class CrossingNumberAlgorithm:
    def __init__(self, SPLIT, mesh_data):
        #self.origin = np.array([(x_max+x_min)/2, (y_max+y_min)/2, (z_max+z_min)/2])
        # 表面形状データの読み込み
        print("読み込み中...")
        f = open(mesh_data)
        self.list1 = []
        self.list2 = []
        count = 1
        for line in f.readlines():
            # vertexとfaceに分けてlist化
            if count < SPLIT: # vertex 形状データによって値を変更
                self.list1.append(line.split())
            else: # list
                self.list2.append(line.split())
            count+=1
        # strをfloatに変換．list1はメートル表記からmm表記に変換
        self.list1 = [[float(x) * 1000 for x in y] for y in self.list1]
        self.list2 = [[int(x) for x in y] for y in self.list2]
        print('list1 = ',self.list1[:5])
        print('list2 = ',self.list2[:5])
        print("読み込み終了")

        #surface用にnp変換
        self.np_list1 = np.array(self.list1)
        self.np_list2 = np.array(self.list2)
        self.np_list1 = self.np_list1[:,0:3] #[x, y, z]座標
        self.np_list2 = self.np_list2[:,1:4] #[点番号１, 点番号2, 点番号3] ex[0, 5, 2]
        # print('np_list1 = ',self.np_list1[:5])
        # print('np_list2 = ',self.np_list2[:5])    


    #物体表面上に点を生成（キカラボさん手法）       
    def surface_kikalab(self, fixed_points, PITCH, PDS_PITCH):
        tentative_points = [] # 物体表面上の点群用
        for list2 in self.np_list2:
            #値の設定
            p0 = self.np_list1[list2[0]]
            p1 = self.np_list1[list2[1]]
            p2 = self.np_list1[list2[2]]
            base0 = np.linalg.norm(p1-p2) #三角形の底辺
            base1 = np.linalg.norm(p2-p0) #三角形の底辺
            base2 = np.linalg.norm(p0-p1) #三角形の底辺
            area = calculate_triangle_area(p0, p1, p2) #三角形の面積
            h0 = 2*area/base0
            h1 = 2*area/base1
            h2 = 2*area/base2
            hight = max([h0, h1, h2]) # 高さを比較し，最大の高さを適用
            #print('hight比較：',h0,h1,h2,hight)
            ndiv = int(hight/PITCH) + 1.0 #int()+1とすることでndiv>=1となるため，hがpitchより小さいとき生成される点は三頂点のみになる
            l0_points = [] #l0上の区分点
            l1_points = [] #l1上の区分点

            #l0, l1上をndiv等分した点の座標を求める
            i = 0
            while i <= ndiv:
                l0_points.append(p0 + i*(p1-p0)/ndiv)
                l1_points.append(p0 + i*(p2-p0)/ndiv)
                i+=1
            
            #di上にndiv_i等分した点の座標を求める
            i = 0
            while i < len(l0_points):
                ndiv_i = int(np.linalg.norm(l0_points[i]-l1_points[i])/PITCH)+1.0
                j = 0
                while j <= ndiv_i:
                    tentative_points.append(l0_points[i] + j * (l1_points[i] - l0_points[i]) / ndiv_i)
                    j+=1
                i+=1
        
        # post_tentative_points = thinning(tentative_points, RATE_OF_THINNINGS) # 間引き
        post_tentative_points = thinning_pds(tentative_points, PDS_PITCH) # 間引き
        print('間引き前len(tentative_points) = ', len(tentative_points))
        print('間引き後len(tentative_points) = ', len(post_tentative_points))
        fixed_points.extend(post_tentative_points)
        print('len(fixed_points) = ', len(fixed_points))


    def cramer(self, pds_point): #クラメルの公式を用いて内外判定
        #ray = self.origin - pds_point
        ray = np.array([10, 10, 10]) - pds_point
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
            #print('u, v, t = ', u, v, t)

            if t>=0 and u>=0 and u<=1 and v>=0 and v<=1 and u+v>=0 and u+v<=1:
                cross_num += 1 #三角形の平面内で交点を持つ → カウント
                #print('u, v, t = ', u, v, t)
                #print('{:.1000f}'.format(u))

            else:
                continue #三角形の平面外で交点を持つ
        
        #print('cross_num', cross_num)
        if cross_num % 2 == 0:
            flg = False #外側
        else:
            flg = True #内側
        #print('内外判定結果：', flg)
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

# ３点の座標を引数とし、その3点で構成される三角形の面積を返す
def calculate_triangle_area(p0, p1, p2): 
    # ベクトルの外積を計算
    cross_product = np.cross(p1 - p0, p2 - p0)
    # 三角形の面積はベクトルの長さの半分
    area = 0.5 * np.linalg.norm(cross_product)
    return area

# 直線p1-p2上にl/ndiv間隔で点列を取得
def point_sequence(p1, p2, ndiv):
    #単位ベクトルの生成
    e = (p2-p1)/np.linalg.norm(p2-p1)
    points=[]
    i=0
    while i < ndiv:
        points.append(p1 + i*e)
    return points


# 間引き
def thinning(points, RATE_OF_THINNINGS):
    points, _ = np.unique(points, return_index=True, axis=0) #重複した座標を削除．return_index=Trueとすることでsortなしになる
    print(points)
    points = points.tolist() # ndarrayをlistに変換
    print('重複削除len(tentative) = ' , len(points))
    num_trials = math.floor(len(points)*(1.0 - RATE_OF_THINNINGS))
    print('num_trials = ' , num_trials)

    for _ in range(num_trials):
        num_points = len(points)
        target = random.randint(0, num_points-1)
        del points[target]

    return points



# PDSを用いた間引き
def thinning_pds(points, PDS_PITCH):
    points, _ = np.unique(points, return_index=True, axis=0) #重複した座標を削除．return_index=Trueとすることでsortなしになる
    print(points)
    points = points.tolist() # ndarrayをlistに変換
    print('重複削除len(tentative) = ' , len(points))
    PDS_PITCH_SQUARE = PDS_PITCH**2
    fixed_points = [] # 確定点
    fixed_points.append(points[0])
    i = 0
    for attention_point in points:
        flg = True
        for fixed_point in fixed_points:
            x1, y1, z1 = attention_point
            x2, y2, z2 = fixed_point
            distance_square = (x2 - x1)**2 + (y2 - y1)**2 + (z2 - z1)**2

            if distance_square < PDS_PITCH_SQUARE:
                # PDSピッチ内に他の点が存在したとき、注目点と確定点の比較を終了。
                flg = False
                pass
        
        if flg:
            fixed_points.append(attention_point)
        i+=1
        print('i : ', i)
    print('間引き前表面生成点(重複削除済み): ', i)

    return fixed_points