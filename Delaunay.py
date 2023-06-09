import numpy as np
from scipy.spatial import Delaunay
import time
from modules import CrossingNumberAlgorithm

# define
SPLIT = 703 # CNAで用いる．vertexとfaceの分け目．

start = time.time()  # 時間計測用
edges = [] # PLYファイルのedge用


Input_file = 'Output/result_main/result_112.csv' # Inputファイル
Output_file = 'Output/Delaunay_Python/delaunay_mesh.ply' # Outputファイル

# ファイルの読み込み。物体の頂点を定義する
vertices = np.loadtxt(Input_file, delimiter=',')
print('vertices', vertices[:3])

# Delaunay分割の作成
tri = Delaunay(vertices)
print('simplices：', tri.simplices)

# CrossNumberAlgorithm
CNA = CrossingNumberAlgorithm.CrossingNumberAlgorithm(SPLIT)

print('len(tri.simplice)：', len(tri.simplices))
i = 0

# edgeノードの作成
for simplice in tri.simplices:
    print(i)
    i+=1

    # 各稜線の中点が物体の内側か外側かを判断する。内側だとTrue、外側だとFalseを返す。
    flg01 = CNA.cramer((vertices[simplice[0]]+vertices[simplice[1]])/2)
    flg02 = CNA.cramer((vertices[simplice[0]]+vertices[simplice[2]])/2)
    flg03 = CNA.cramer((vertices[simplice[0]]+vertices[simplice[3]])/2)
    flg12 = CNA.cramer((vertices[simplice[1]]+vertices[simplice[2]])/2)
    flg13 = CNA.cramer((vertices[simplice[1]]+vertices[simplice[3]])/2)
    flg23 = CNA.cramer((vertices[simplice[2]]+vertices[simplice[3]])/2)

    # Trueのとき，edge生成
    if flg01: edges.append([simplice[0] + simplice[1]])
    if flg02: edges.append([simplice[0] + simplice[2]])
    if flg03: edges.append([simplice[0] + simplice[3]])
    if flg12: edges.append([simplice[1] + simplice[2]])
    if flg13: edges.append([simplice[1] + simplice[3]])
    if flg23: edges.append([simplice[2] + simplice[3]])


# plyで保存
with open(Output_file, 'w', newline="") as f:
    f.write('ply\n')
    f.write('format ascii 1.0\n')
    f.write('comment VCGLIB generated\n')
    f.write('element vertex '+str(len(vertices))+'\n')
    f.write('property float x\n')
    f.write('property float y\n')
    f.write('property float z\n')
    f.write('element edge '+str(len(edges))+'\n')
    f.write('property int vertex1\n')
    f.write('property int vertex2\n')
    f.write('end_header\n')

    for ele in vertices: # 点群の座標値入力
        f.write(str(ele[0])+' '+str(ele[1])+' '+str(ele[2])+' '+'\n')
    
    for ele in edges: # 三角形を構成する点群のノード入力
        f.write(str(ele[0])+' '+str(ele[1])+'\n')

# 計測結果
elapsed_time = time.time() - start
print ("elapsed_time:{0}".format(elapsed_time) + "[sec]")

