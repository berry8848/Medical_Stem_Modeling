import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.spatial import Delaunay

# 3次元データの生成
points = np.random.rand(20, 3)  # 20個のランダムな3次元点を生成

# Delaunay分割の作成
tri = Delaunay(points)

# 結果の可視化
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# 3次元点をプロット
ax.scatter(points[:, 0], points[:, 1], points[:, 2])

# 分割三角形をプロット
ax.plot_trisurf(points[:, 0], points[:, 1], points[:, 2], triangles=tri.simplices)

# 軸ラベルの設定
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')

# 表示
plt.show()