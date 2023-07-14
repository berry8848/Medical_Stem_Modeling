# DrawingFaceと結果は同じ。numpy-stlを使っているかどうか。
# 内容が分かりやすいため、基本こっちを使う。

import numpy as np
import time
from modules import CrossingNumberAlgorithm

# define
SPLIT = 703 # CNAで用いる．vertexとfaceの分け目．

start = time.time()  # 時間計測用

Input_file = 'Input/Delaunay_C++/check2_delaunay.csv' # Inputファイル
Output_file = 'Output/Drawing_Mesh/Drawing_Mesh_ver1.ply' # Outputファイル

vertices = [] # Input_file用
faces = [] # 面法線ベクトル用

# ファイルの読み込み。物体の頂点を定義する
vertices = np.loadtxt(Input_file, delimiter=',')

# CrossNumberAlgorithm
CNA = CrossingNumberAlgorithm.CrossingNumberAlgorithm(SPLIT)

print(len(vertices))

# ノードの作成
for i in range(0, len(vertices), 4): 
    # 各稜線の中点が物体の内側か外側かを判断する。内側だとTrue、外側だとFalseを返す。
    flg01 = CNA.cramer((vertices[i]+vertices[i+1])/2)
    flg02 = CNA.cramer((vertices[i]+vertices[i+2])/2)
    flg03 = CNA.cramer((vertices[i]+vertices[i+3])/2)
    flg12 = CNA.cramer((vertices[i+1]+vertices[i+2])/2)
    flg13 = CNA.cramer((vertices[i+1]+vertices[i+3])/2)
    flg23 = CNA.cramer((vertices[i+2]+vertices[i+3])/2)

    # 全てTrueのときのみfaceを生成。
    if flg01 and flg02 and flg12:
        faces.append([i, i+1, i+2])

    if flg02 and flg03 and flg23:
        faces.append([i, i+2, i+3])

    if flg03 and flg01 and flg13:
        faces.append([i, i+3, i+1])

    if flg12 and flg13 and flg23:
        faces.append([i+1, i+2, i+3])
    
    print(i)

# listをnumpyに変換
faces = np.array(faces)

# plyで保存
with open(Output_file, 'w', newline="") as f:
    f.write('ply\n')
    f.write('format ascii 1.0\n')
    f.write('comment VCGLIB generated\n')
    f.write('element vertex '+str(len(vertices))+'\n')
    f.write('property float x\n')
    f.write('property float y\n')
    f.write('property float z\n')
    f.write('element face '+str(len(faces))+'\n')
    f.write('property list uchar int vertex_indices\n')
    f.write('end_header\n')

    for ele in vertices: # 点群の座標値入力
        f.write(str(ele[0])+' '+str(ele[1])+' '+str(ele[2])+' '+'\n')
    
    for ele in faces: # 三角形を構成する点群のノード入力
        f.write('3'+' '+str(ele[0])+' '+str(ele[1])+' '+str(ele[2])+' '+'\n')

# 計測結果
elapsed_time = time.time() - start
print ("elapsed_time:{0}".format(elapsed_time) + "[sec]")
