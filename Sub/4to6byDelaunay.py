# 目的：ドロネー分割を節点番号表記に変更．その後，PDSによる点群生成．
# Input：ドロネー分割の結果（座標値のみ：③の結果）＆　②の結果
# Output：PDSの点群座標.plyファイル，PDSの点群座標.csvファイル

from platform import node
from tkinter import NE
from modules import Tetrahedron
from modules import Point

import numpy as np
import csv
import math
import time
import random

start = time.time()  # 時間計測用

file_name_1 = './Input/Column10_tri1220_delaunay.csv' # C++でのドロネー分割の結果ファイル
file_name_2 = './Input/Column10_0615.csv' # ANSYSのデータファイル
# file_name_1 = './Input/joint_delaunay.csv' # C++でのドロネー分割の結果ファイル
# file_name_2 = './Input/joint.csv' # ANSYSのデータファイル


# Outputファイル
#file_name_3 = './Output/joint_delaunay_x.csv' # ドロネー分割の四面体のノード番号を出力する．ドロネー分割の節点番号表記を見たい場合#を外す．
file_resultPly = 'Output/result_4to6.ply'
file_resultCsv = 'Output/result_4to6.csv'

a1 = [] # file1の入力用
a2 = [] # file2の入力用
delaunay_list = [] # ドロネー分割の節点番号表記
points_obj_list = []  # Pointのオブジェクトを保持。
divide_tetrahedron_list = []  # このリストに四面体分割を保持していく。[p1, p2, p3, p4]．※ p1 = [x1, y1, z1]

# ファイル1読み込み
with open(file_name_1) as f1:
    reader1 = csv.reader(f1)
    for row in reader1:
        a1.append(row)

# ファイル2読み込み
with open(file_name_2) as f2:
    reader2 = csv.reader(f2)
    for row in reader2:
        a2.append(row)

# ドロネー分割を節点番号での表記に変換
for b1 in a1:
    c = [None, None, None, None]
    flg0 = False
    flg1 = False
    flg2 = False
    flg3 = False
    b1 = [float(n) for n in b1]

    for b2 in a2:
        b2 = [float(n) for n in b2]
        if b1[0] == b2[1] and b1[1] == b2[2] and b1[2] == b2[3]:
            c[0] = b2[0]
            flg0 = True
        
        if b1[3] == b2[1] and b1[4] == b2[2] and b1[5] == b2[3]:
            c[1] = b2[0]
            flg1 = True

        if b1[6] == b2[1] and b1[7] == b2[2] and b1[8] == b2[3]:
            c[2] = b2[0]
            flg2 = True
        
        if b1[9] == b2[1] and b1[10] == b2[2] and b1[11] == b2[3]:
            c[3] = b2[0]
            flg3 = True
        
        if flg0 and flg1 and flg2 and flg3:
            delaunay_list.append(c)
            break
    #print(delaunay_list[-1])
print(delaunay_list[0:10])

# # ファイル書き込み．ドロネー分割の節点番号表記を見たい場合実行．それ以外は必要なし．
# with open(file_name_3, 'w', newline="") as f:
#     writer = csv.writer(f)
#     writer.writerows(delaunay_list)

# ANSYS上の点群を取得し座標値を取得 ※取得方法不明のため，csvからの読み取りに臨時変更
points = np.loadtxt(file_name_2, delimiter=',')
print(points[0:10])

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

# FEMの節点の読み込み
for i in range(len(points)):
    point = Point.Point(points[i]) # クラスに格納
    point.system_guid_obj_to_coordinate()
    points_obj_list.append(point)

# divide_tetrahedron_listの作成
# delaunay_listに格納された四面体の節点番号から，点の座標値と応力値を紐づけ．※ 節点番号p1→[座標値，応力値]
flgnum = 0
for tri in delaunay_list:
    c = [None, None, None, None]
    flg = False
    flg0 = False
    flg1 = False
    flg2 = False
    flg3 = False

    for select_point in points_obj_list:
        if int(select_point.node) == int(tri[0]):
            c[0] = select_point
            flg0 = True
        
        if int(select_point.node) == int(tri[1]):
            c[1] = select_point
            flg1 = True

        if int(select_point.node) == int(tri[2]):
            c[2] = select_point
            flg2 = True
        
        if int(select_point.node) == int(tri[3]):
            c[3] = select_point
            flg3 = True
        
        if flg0 and flg1 and flg2 and flg3:
            new_tetrahedron = Tetrahedron.Tetrahedron(c[0], c[1], c[2], c[3])
            new_tetrahedron.set()
            divide_tetrahedron_list.append(new_tetrahedron)
            flg1 = new_tetrahedron.cul_center_p_and_radius()
            if flg1:
                del divide_tetrahedron_list[-1]
                print("Error",flgnum)
            flgnum += 1
            flg = True
            break

divide_tetrahedron_list_count = len(divide_tetrahedron_list)
print(len(divide_tetrahedron_list))

# 以下，PDS
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
    return 24*long

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

# pds_point が含まれる四面体を divide_tetrahedron_list から捜索．点が連続でN回生成できなかったら終了
N = 80
num = 0
fixed_points = [] # 確定点


while num < N:
    if len(fixed_points) >= 450:
        break
    count = 0
    # ランダムな点を生成
    pds_x = random.uniform(x_min, x_max)
    pds_y = random.uniform(y_min, y_max)
    pds_z = random.uniform(z_min, z_max)
    pds_point = [pds_x, pds_y, pds_z]
    # 新点の座標情報をPointに格納し候補点にする
    candidate_point = Point.Point(pds_point)
    candidate_point.pds_coordinate()
    # 最初に探索するドロネー四面体をランダムに選択
    i = random.randint(0, len(divide_tetrahedron_list)-1)
    NextTri = divide_tetrahedron_list[i]

    # 新点を内包するドロネー四面体を探索
    for _ in range(len(divide_tetrahedron_list)):
        # 内包するかを確認
        check, node_list = NextTri.check(candidate_point)

        # 内包する四面体のとき
        if check:
            if not fixed_points:
                fixed_points.append(pds_point)
                break
            new_stress = NextTri.output_stress(pds_point)
            density = stress_to_density(new_stress)
            long = density_to_long(density)
            flg = check_distance(fixed_points, candidate_point, long)
            # 点間距離内に他の点が存在しないとき
            if flg:
                fixed_points.append(pds_point)
                num = 0
                print('num : ', num, 'fixed_points : ', len(fixed_points))
                #time.sleep(0.5)
                break
            else :
                num = num + 1
                print('num : ', num)
                break

        # 四面体の外側に点があるとき
        else:
            for tetrahedron in divide_tetrahedron_list:
                list_and = set(tetrahedron.node) & set(node_list)
                if len(list_and) == 3:
                    list_and = set(tetrahedron.node) & set(NextTri.node)
                    if len(list_and) == 4:
                        count += 1
                        #print("pass", count)
                        pass
                    else:
                        #print('NextTri = ',NextTri.node)
                        NextTri = tetrahedron
                        count += 1
                        break
            if count > divide_tetrahedron_list_count/50:
                print("count ", divide_tetrahedron_list_count/50, "!!")
                break


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