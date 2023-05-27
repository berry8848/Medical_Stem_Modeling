# 目的：ドロネー分割を節点番号表記に変更．その後，PDSによる点群生成．
# Input：ドロネー分割の結果（座標値のみ：③の結果）＆　②の結果
# Output：PDSの点群座標.plyファイル，PDSの点群座標.csvファイル

from modules import Point
from modules import Gauss
from modules import Biharmonic
from modules import CrossingNumberAlgorithm

import numpy as np
import csv
import math
import time
import random

start = time.time()  # 時間計測用

file_name_2 = './Input/Column10_0615.csv' # ANSYSのデータファイル
# file_name_2 = './Input/square1220.csv' # ANSYSのデータファイル\
# file_name_2 = './Input/joint.csv' # ANSYSのデータファイル

# Outputファイル
file_resultPly = 'Output/result_main/result.ply'
file_resultCsv = 'Output/result_main/result.csv'

a2 = [] #file2の入力用
points_obj_list = []  #Pointのオブジェクトを保持。

#ファイル2読み込み
with open(file_name_2) as f2:
    reader2 = csv.reader(f2)
    for row in reader2:
        a2.append(row) 

# ANSYS上の点群を取得し座標値を取得 ※取得方法不明のため，csvからの読み取りに臨時変更
points = np.loadtxt(file_name_2, delimiter=',')
print('points = ',points[0:10])


# FEMの節点の読み込み
for i in range(len(points)):
    point = Point.Point(points[i]) # クラスに格納
    point.system_guid_obj_to_coordinate()
    points_obj_list.append(point)

# Gauss
gauss = Gauss.Gauss(points)
zs, lambdas , cs= gauss.gauss1()
print('lambdas: ', lambdas)
print('cs: ', cs)

# Biharmonic
biharmonic = Biharmonic.Biharmonic(points, cs, lambdas)
print(biharmonic.cal(15.78, 4.09, 52.35))



#start = time.time()  # 時間計測用

# 以下，PDS
# PDSでの点の生成範囲の設定
x_max = points[0][1]
x_min = points[0][1]
y_max = points[0][2]
y_min = points[0][2]
z_max = points[0][3]
z_min = points[0][3]

for point in points:
    if point[1] > x_max:
        x_max = point[1]
    if point[1] < x_min:
        x_min = point[1]
    
    if point[2] > y_max:
        y_max = point[2]
    if point[2] < y_min:
        y_min = point[2]

    if point[3] > z_max:
        z_max = point[3]
    if point[3] < z_min:
        z_min = point[3]

# PDSでの点の生成範囲の表示
print("x_max = ", x_max, "x_min = ", x_min)
print("y_max = ", y_max, "y_min = ", y_min)
print("z_max = ", z_max, "z_min = ", z_min)

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
    return 10*long

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

# 点が連続でN回生成できなかったら終了
N = 800
num = 0
fixed_points = [] # 確定点

# CrossNumberAlgorithm
CNA = CrossingNumberAlgorithm.CrossingNumberAlgorithm()

while num < N:
    if len(fixed_points) >= 0:
        break

    flg_P = False
    while not flg_P: # 物体内部の点を生成するまでループ。内部であればflg_PはTrueとなる。
        # ランダムな点を生成
        pds_x = random.uniform(x_min, x_max)
        pds_y = random.uniform(y_min, y_max)
        pds_z = random.uniform(z_min, z_max)
        pds_point = [pds_x, pds_y, pds_z]
        #pds_point = [0, 0, 5]
        # 物体内部の点か判定
        flg_P = CNA.cramer(pds_point)
        #print(flg_P)


    # 生成点の座標情報をPointに格納し候補点にする
    candidate_point = Point.Point(pds_point)
    candidate_point.pds_coordinate()
    # Biharmonicより候補点の応力を導出
    new_stress = biharmonic.cal(candidate_point.x, candidate_point.y, candidate_point.z)
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


#物体表面上でPDS
# CNA.biharmonic(points, cs, lambdas)
# CNA.surface_pds(fixed_points)
# CNA.surface(fixed_points)
CNA.surface_kikalab(fixed_points)

#重複した座標を削除
fixed_points = np.unique(fixed_points, axis=0)



# PDS点群に表面点群を連結
#CNA.outline(fixed_points)

# ply にPDSの結果出力
print("fixed_points  = ", len(fixed_points), "個")
with open(file_resultPly, 'w', newline="") as f:
    f.write('ply\n')
    f.write('format ascii 1.0\n')
    f.write('element vertex '+str(len(fixed_points))+'\n')
    f.write('property double x\n')
    f.write('property double y\n')
    f.write('property double z\n')
    f.write('end_header\n')
    for ele in fixed_points:
        f.write(str(ele[0])+' '+str(ele[1])+' '+str(ele[2])+' '+'\n')

# csv にPDSの結果出力
with open(file_resultCsv, 'w', newline="") as f:
    writer = csv.writer(f)
    writer.writerows(fixed_points)


# 計測結果
elapsed_time = time.time() - start
print ("elapsed_time:{0}".format(elapsed_time) + "[sec]")
