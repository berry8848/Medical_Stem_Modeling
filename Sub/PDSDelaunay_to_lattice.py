# delaunay.cpp で出力した結果を、四面体を構成する4点のうち、
# 2点(辺)ごとに整理するプログラム
# 2022/10/23時点で、DrawingMeshやDrawingFaceを用いる場合似たような内容のため、このプログラムは使わない。
# こっちは線分で出力、DrawinMeshらは面で出力

import numpy as np
import csv

# ","で区切る
arr1 = np.loadtxt('Input/result_4to6_1150_450_delaunay.csv', delimiter=',')
lattices = []

# for arr in arr1:
#     arr2.append([arr[0], arr[1], arr[2], arr[3], arr[4], arr[5]])
#     arr2.append([arr[0], arr[1], arr[2], arr[6], arr[7], arr[8]])
#     arr2.append([arr[0], arr[1], arr[2], arr[9], arr[10], arr[11]])
#     arr2.append([arr[3], arr[4], arr[5], arr[6], arr[7], arr[8]])
#     arr2.append([arr[3], arr[4], arr[5], arr[9], arr[10], arr[11]])
#     arr2.append([arr[6], arr[7], arr[8], arr[9], arr[10], arr[11]])

for i in range(0, len(arr1), 4): 
    lattices.append([arr1[i][0],arr1[i][1],arr1[i][2], arr1[i+1][0], arr1[i+1][1], arr1[i+1][2]])
    lattices.append([arr1[i][0],arr1[i][1],arr1[i][2], arr1[i+2][0], arr1[i+2][1], arr1[i+2][2]])
    lattices.append([arr1[i][0],arr1[i][1],arr1[i][2], arr1[i+3][0], arr1[i+3][1], arr1[i+3][2]])
    lattices.append([arr1[i+1][0], arr1[i+1][1], arr1[i+1][2], arr1[i+2][0], arr1[i+2][1], arr1[i+2][2]])
    lattices.append([arr1[i+1][0], arr1[i+1][1], arr1[i+1][2], arr1[i+3][0], arr1[i+3][1], arr1[i+3][2]])
    lattices.append([arr1[i+2][0], arr1[i+2][1], arr1[i+2][2], arr1[i+3][0], arr1[i+3][1], arr1[i+3][2]])

print(lattices[:9])
print(len(lattices))

# ファイル書き込み
with open('Output/result_4to6_1150_450_delaunay_lattice.csv', 'w', newline="") as f:
    writer = csv.writer(f)
    writer.writerows(lattices)
