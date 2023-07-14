# DrawingFaceと結果は同じ。numpy-stlを使っているかどうか。
# 内容が分かりやすいため、基本こっちを使う。

import numpy as np
import math
import csv
import time
from modules import CrossingNumberAlgorithm

# define
SPLIT = 3313 # CNAで用いる．vertexとfaceの分け目．

start = time.time()  # 時間計測用

Input_file = 'Input/Delaunay_C++/result_2002_delaunay.csv' # Inputファイル
Output_file = 'Output/Drawing_Mesh/Drawing_Mesh_ver2.ply' # Outputファイル

vertices = [] # Input_file用
edges = [] # ply edge用

# ファイルの読み込み。物体の頂点を定義する
vertices = np.loadtxt(Input_file, delimiter=',')

# CrossNumberAlgorithm
CNA = CrossingNumberAlgorithm.CrossingNumberAlgorithm(SPLIT)

print(len(vertices))
print(vertices)


# voxelの範囲指定
voxel_max = np.max(vertices, axis=0) #(x_max, y_max, z_max)
voxel_max = [math.ceil(voxel_max[0]+1), math.ceil(voxel_max[1]+1),math.ceil(voxel_max[2]+1)] #max値に1を足し、小数点切り上げ。1を足さないとスーパーボックスの表面上に点を含んでしまう
voxel_min = np.min(vertices, axis=0)
print('voxel_max: ', voxel_max)
# voxelの初期値      
x0 = math.floor(voxel_min[0]-1)
y0 = math.floor(voxel_min[1]-1)
z0 = math.floor(voxel_min[2]-1)
print('x0: ',x0,'y0: ',y0,'z0: ',z0)
# voxelのピッチ指定
px = 2
py = 2
pz = 4


# voxel生成
indx_x = int(math.ceil((voxel_max[0]-x0)/px))
indx_y = int(math.ceil((voxel_max[1]-y0)/py))
indx_z = int(math.ceil((voxel_max[2]-z0)/pz))
print('indx = ' ,indx_x, indx_y, indx_z)
voxel = [[[0 for k in range(indx_x)] for j in range(indx_y)] for i in range(indx_z)]
#print('voxel = ', voxel)
print('shape of voxel = ', np.shape(voxel))

k = 0
while k<indx_z:
    j = 0
    print(k)
    while j<indx_y:
        i = 0
        while i<indx_x:
            #voxelの頂点が物体の内部にあるかの判定。内側だとTrue、外側だとFalseを返す。
            flg1 = CNA.cramer([x0+i*px, y0+j*py, z0+k*pz])
            flg2 = CNA.cramer([x0+i*px, y0+j*py, z0+(k+1)*pz])
            flg3 = CNA.cramer([x0+i*px, y0+(j+1)*py, z0+k*pz])
            flg4 = CNA.cramer([x0+i*px, y0+(j+1)*py, z0+(k+1)*pz])
            flg5 = CNA.cramer([x0+(i+1)*px, y0+j*py, z0+k*pz])
            flg6 = CNA.cramer([x0+(i+1)*px, y0+j*py, z0+(k+1)*pz])
            flg7 = CNA.cramer([x0+(i+1)*px, y0+(j+1)*py, z0+k*pz])
            flg8 = CNA.cramer([x0+(i+1)*px, y0+(j+1)*py, z0+(k+1)*pz])

            if flg1 and flg2 and flg3 and flg4 and flg5 and flg6 and flg7 and flg8: #全て内側
                voxel[i][j][k] = 0
            elif (not flg1) and (not flg2) and (not flg3) and (not flg4) and (not flg5) and (not flg6) and (not flg7) and (not flg8):#全て外側
                voxel[i][j][k] = 2
            else: #物体表面を含むvoxel
                voxel[i][j][k] = 1
            print(i, j, k, 'voxel = ',voxel[i][j][k])
            i+=1
        j+=1
    k+=1
    
elapsed_time_1 = time.time() - start
print('finish')


for i in range(0, len(vertices),4):
    #頂点01
    voxel_indx_x = (((vertices[i][0] + vertices[i+1][0])/2)-x0) // px
    voxel_indx_y = (((vertices[i][1] + vertices[i+1][1])/2)-y0) // py
    voxel_indx_z = (((vertices[i][2] + vertices[i+1][2])/2)-z0) // pz
    #print('0 : ',int(voxel_indx_x),int(voxel_indx_y),int(voxel_indx_z))
    flg01 = voxel[int(voxel_indx_x)][int(voxel_indx_y)][int(voxel_indx_z)]
    if flg01==0:
        edges.append([i, i+1])
    elif flg01==1:
        flg = CNA.cramer((vertices[i]+vertices[i+1])/2)
        if flg:
            edges.append([i, i+1])
    #頂点12
    voxel_indx_x = (((vertices[i+1][0] + vertices[i+2][0])/2)-x0) // px
    voxel_indx_y = (((vertices[i+1][1] + vertices[i+2][1])/2)-y0) // py
    voxel_indx_z = (((vertices[i+1][2] + vertices[i+2][2])/2)-z0) // pz
    #print('1 : ',int(voxel_indx_x),int(voxel_indx_y),int(voxel_indx_z))
    flg12 = voxel[int(voxel_indx_x)][int(voxel_indx_y)][int(voxel_indx_z)]

    if flg12==0:
        edges.append([i+1, i+2])
    elif flg12==1:
        flg = CNA.cramer((vertices[i+1]+vertices[i+2])/2)
        if flg:
            edges.append([i+1, i+2])
    #頂点23
    voxel_indx_x = (((vertices[i+2][0] + vertices[i+3][0])/2)-x0) // px
    voxel_indx_y = (((vertices[i+2][1] + vertices[i+3][1])/2)-y0) // py
    voxel_indx_z = (((vertices[i+2][2] + vertices[i+3][2])/2)-z0) // pz
    #print('2 : ',int(voxel_indx_x),int(voxel_indx_y),int(voxel_indx_z))
    flg23 = voxel[int(voxel_indx_x)][int(voxel_indx_y)][int(voxel_indx_z)]
    if flg23==0:
        edges.append([i+2, i+3])
    elif flg23==1:
        flg = CNA.cramer((vertices[i+2]+vertices[i+3])/2)
        if flg:
            edges.append([i+2, i+3])
    #頂点03
    voxel_indx_x = (((vertices[i+3][0] + vertices[i][0])/2)-x0) // px
    voxel_indx_y = (((vertices[i+3][1] + vertices[i][1])/2)-y0) // py
    voxel_indx_z = (((vertices[i+3][2] + vertices[i][2])/2)-z0) // pz
    #print('z', (vertices[i+3][2] + vertices[i][2])/2)
    #print('3 : ',int(voxel_indx_x),int(voxel_indx_y),int(voxel_indx_z))
    flg03 = voxel[int(voxel_indx_x)][int(voxel_indx_y)][int(voxel_indx_z)]
    if flg03==0:
        edges.append([i+3, i])
    elif flg03==1:
        flg = CNA.cramer((vertices[i+3]+vertices[i])/2)
        if flg:
            edges.append([i+3, i])
    #頂点02
    voxel_indx_x = (((vertices[i+2][0] + vertices[i][0])/2)-x0) // px
    voxel_indx_y = (((vertices[i+2][1] + vertices[i][1])/2)-y0) // py
    voxel_indx_z = (((vertices[i+2][2] + vertices[i][2])/2)-z0) // pz
    #print('4 : ',int(voxel_indx_x),int(voxel_indx_y),int(voxel_indx_z))
    flg02 = voxel[int(voxel_indx_x)][int(voxel_indx_y)][int(voxel_indx_z)]
    if flg02==0:
        edges.append([i+2, i])
    elif flg02==1:
        flg = CNA.cramer((vertices[i+2]+vertices[i])/2)
        if flg:
            edges.append([i+2, i])
    #頂点13
    voxel_indx_x = (((vertices[i+1][0] + vertices[i+3][0])/2)-x0) // px
    voxel_indx_y = (((vertices[i+1][1] + vertices[i+3][1])/2)-y0) // py
    voxel_indx_z = (((vertices[i+1][2] + vertices[i+3][2])/2)-z0) // pz
    #print('5 : ',int(voxel_indx_x),int(voxel_indx_y),int(voxel_indx_z))
    flg13 = voxel[int(voxel_indx_x)][int(voxel_indx_y)][int(voxel_indx_z)]
    if flg13==0:
        edges.append([i+1, i+3])
    elif flg13==1:
        flg = CNA.cramer((vertices[i+1]+vertices[i+3])/2)
        if flg:
            edges.append([i+1, i+3])
    print(i)

# listをnumpyに変換
edges = np.array(edges)



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
elapsed_time_2 = time.time() - start
print ("elapsed_time1:{0}".format(elapsed_time_1) + "[sec]")
print ("elapsed_time2:{0}".format(elapsed_time_2) + "[sec]")
