# Delaunay分割後、メッシュ生成

import numpy as np
from scipy.spatial import Delaunay
import time
from modules import CrossingNumberAlgorithm

start = time.time()  # 時間計測用
edges = [] # PLYファイルのedge用


# Input_file = 'Output/result_main/check/check3_delaunay.csv' # Inputファイル
Input_file = 'Output/result_main/result_144.csv' # Inputファイル
mesh_data = 'Input/Mesh_Data/cube_50x50mm_mesh.ply' # 物体の表面形状データ。
# mesh_data = 'Input/Mesh_Data/cube_50x50_check.ply' # 物体の表面形状データ。
Output_file = 'Output/Delaunay_Python/delaunay_mesh.ply' # Outputファイル

# ファイルの読み込み。物体の頂点を定義する
vertices = np.loadtxt(Input_file, delimiter=',')
print('vertices', vertices[:3])

# Delaunay分割の作成
tri = Delaunay(vertices)
print('simplices：', tri.simplices)

# CrossNumberAlgorithm
CNA = CrossingNumberAlgorithm.CrossingNumberAlgorithm(mesh_data)

print('len(tri.simplice)：', len(tri.simplices))
i = 0

# # edgeノードの作成
# for simplice in tri.simplices:
#     # 各稜線の中点が物体の内側か外側かを判断する。内側だとTrue、外側だとFalseを返す。
#     flg01 = CNA.cramer((vertices[simplice[0]]+vertices[simplice[1]])/2)
#     flg02 = CNA.cramer((vertices[simplice[0]]+vertices[simplice[2]])/2)
#     flg03 = CNA.cramer((vertices[simplice[0]]+vertices[simplice[3]])/2)
#     flg12 = CNA.cramer((vertices[simplice[1]]+vertices[simplice[2]])/2)
#     flg13 = CNA.cramer((vertices[simplice[1]]+vertices[simplice[3]])/2)
#     flg23 = CNA.cramer((vertices[simplice[2]]+vertices[simplice[3]])/2)

#     # Trueのとき，edge生成
#     if flg01: edges.append([simplice[0], simplice[1]])
#     if flg02: edges.append([simplice[0], simplice[2]])
#     if flg03: edges.append([simplice[0], simplice[3]])
#     if flg12: edges.append([simplice[1], simplice[2]])
#     if flg13: edges.append([simplice[1], simplice[3]])
#     if flg23: edges.append([simplice[2], simplice[3]])

#     print(i)
#     i+=1

# edgeノードの作成
for simplice in tri.simplices:
#     print(f"""
#           -------------------- simplice: {simplice} --------------------
#           """)
#     # 各稜線の中点が物体の内側か外側かを判断する。内側だとTrue、外側だとFalseを返す。
#     print(f"""  flg01 :  
#           始点座標： {vertices[simplice[0]]}
#           終点座標： {vertices[simplice[1]]}
#           中間点座標： {(vertices[simplice[0]]+vertices[simplice[1]])/2}
# """)
    flg01 = CNA.majority_vote((vertices[simplice[0]]+vertices[simplice[1]])/2)
    # print('flg01結果：', flg01)
    # print('              ')
    

    # print(f"""  flg02 :  
    #       始点座標： {vertices[simplice[0]]}
    #       終点座標： {vertices[simplice[2]]}
    #       中間点座標： {(vertices[simplice[0]]+vertices[simplice[2]])/2}
    #         """)
    flg02 = CNA.majority_vote((vertices[simplice[0]]+vertices[simplice[2]])/2)
    # print('flg02結果：', flg02)
    # print('              ')

    
    # print(f"""  flg03 :  
    #       始点座標： {vertices[simplice[0]]}
    #       終点座標： {vertices[simplice[3]]}
    #       中間点座標： {(vertices[simplice[0]]+vertices[simplice[3]])/2}
    #         """)
    flg03 = CNA.majority_vote((vertices[simplice[0]]+vertices[simplice[3]])/2)
    # print('flg03結果：', flg03)
    # print('              ')
    

    # print(f"""  flg12 :  
    #       始点座標： {vertices[simplice[1]]}
    #       終点座標： {vertices[simplice[2]]}
    #       中間点座標： {(vertices[simplice[1]]+vertices[simplice[2]])/2}
    #         """)
    flg12 = CNA.majority_vote((vertices[simplice[1]]+vertices[simplice[2]])/2)
    # print('flg12結果：', flg12)
    # print('              ')


    # print(f"""  flg13 :  
    #       始点座標： {vertices[simplice[1]]}
    #       終点座標： {vertices[simplice[3]]}
    #       中間点座標： {(vertices[simplice[1]]+vertices[simplice[3]])/2}
    #         """)
    flg13 = CNA.majority_vote((vertices[simplice[1]]+vertices[simplice[3]])/2)
    # print('flg13結果：', flg13)
    # print('              ')


    # print(f"""  flg23 :  
    #       始点座標： {vertices[simplice[2]]}
    #       終点座標： {vertices[simplice[3]]}
    #       中間点座標： {(vertices[simplice[2]]+vertices[simplice[3]])/2}
    #         """)
    flg23 = CNA.majority_vote((vertices[simplice[2]]+vertices[simplice[3]])/2)
    # print('flg23結果：', flg23)
    # print('              ')


    # flg01 = True
    # flg02 = True
    # flg03 = True
    # flg12 = True
    # flg13 = True
    # flg23 = True



    # Trueのとき，edge生成
    if flg01: edges.append([simplice[0], simplice[1]])
    if flg02: edges.append([simplice[0], simplice[2]])
    if flg03: edges.append([simplice[0], simplice[3]])
    if flg12: edges.append([simplice[1], simplice[2]])
    if flg13: edges.append([simplice[1], simplice[3]])
    if flg23: edges.append([simplice[2], simplice[3]])

    print(i)
    i+=1


# # edgeノードの作成
# for simplice in tri.simplices:
#     # 各稜線の中点が物体の内側か外側かを判断する。内側だとTrue、外側だとFalseを返す。
#     flg01 = CNA.ray_triangle_intersection((vertices[simplice[0]]+vertices[simplice[1]])/2)
#     flg02 = CNA.ray_triangle_intersection((vertices[simplice[0]]+vertices[simplice[2]])/2)
#     flg03 = CNA.ray_triangle_intersection((vertices[simplice[0]]+vertices[simplice[3]])/2)
#     flg12 = CNA.ray_triangle_intersection((vertices[simplice[1]]+vertices[simplice[2]])/2)
#     flg13 = CNA.ray_triangle_intersection((vertices[simplice[1]]+vertices[simplice[3]])/2)
#     flg23 = CNA.ray_triangle_intersection((vertices[simplice[2]]+vertices[simplice[3]])/2)

#     # Trueのとき，edge生成
#     if flg01: edges.append([simplice[0], simplice[1]])
#     if flg02: edges.append([simplice[0], simplice[2]])
#     if flg03: edges.append([simplice[0], simplice[3]])
#     if flg12: edges.append([simplice[1], simplice[2]])
#     if flg13: edges.append([simplice[1], simplice[3]])
#     if flg23: edges.append([simplice[2], simplice[3]])

#     print(i)
#     i+=1


print('edge 重複削除前個数：', len(edges))
#重複した座標を削除
edges = np.unique(edges, axis=0)
print('edge 重複削除後個数：', len(edges))

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

