# DrawingMeshと結果は同じ。numpy-stlを使っているかどうか。
# 基本こっちは使わない
# 最終更新2022/10/23
# ラティス中点除去を行っていない

import numpy as np
from stl import mesh

Input_file = './Input/delaunay/result_4to6_1150_delaunay.csv' # Inputファイル
Output_file = 'Output/Drawing_Face.stl' # Outputファイル

vertices = [] # Input_file用
faces = [] # 面法線ベクトル用

# ファイルの読み込み。物体の頂点を定義する
vertices = np.loadtxt(Input_file, delimiter=',')

# ノードの作成
for i in range(0, len(vertices), 4): 
    faces.append([i, i+1, i+2])
    faces.append([i, i+1, i+3])
    faces.append([i, i+2, i+3])
    faces.append([i+1, i+2, i+3])

# listをnumpyに変換
faces = np.array(faces)

# メッシュ（物体）作成
obj= mesh.Mesh(np.zeros((faces.shape[0]), dtype=mesh.Mesh.dtype))
for i, f in enumerate(faces):
    for j in range(3):
        obj.vectors[i][j] = vertices[f[j],:]

# 保存
obj.save(Output_file)