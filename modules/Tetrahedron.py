
# coding: utf-8
import math
from operator import ne
from platform import node
import numpy as np
from numpy import linalg as LA
import time


class Tetrahedron:

    def __init__(self, p1, p2, p3, p4):
        self.p1 = p1
        self.p2 = p2
        self.p3 = p3
        self.p4 = p4
        self.face1 = None
        self.face2 = None
        self.face3 = None
        self.face4 = None
        self.radius = None
        self.center_p = None
        self.node = None
    
    def set(self): # 四面体の再番号付け　頂点を1する．頂点から見て底面を反時計回りで2, 3, 4
        self.p3, self.p4 = re_number([self.p1, self.p2, self.p3, self.p4])
        self.node = [self.p1.node, self.p2.node, self.p3.node, self.p4.node]
        self.face1 = [self.p1, self.p3, self.p2]
        self.face2 = [self.p1, self.p4, self.p3]
        self.face3 = [self.p1, self.p2, self.p4]
        self.face4 = [self.p2, self.p3, self.p4]
        #print('p3 reset ', self.p3.node, self.p4.node, "node is ", self.node)
    
    def set_print(self):
        return print(self.node)


    def cul_center_p_and_radius(self):
        '''外接球の中心点と半径を計算'''

        a = [
            [self.p1.x, self.p1.y, self.p1.z, 1],
            [self.p2.x, self.p2.y, self.p2.z, 1],
            [self.p3.x, self.p3.y, self.p3.z, 1],
            [self.p4.x, self.p4.y, self.p4.z, 1]
        ]

        d_x = [
            [pow(self.p1.x, 2) + pow(self.p1.y, 2) + pow(self.p1.z, 2), self.p1.y, self.p1.z, 1],
            [pow(self.p2.x, 2) + pow(self.p2.y, 2) + pow(self.p2.z, 2), self.p2.y, self.p2.z, 1],
            [pow(self.p3.x, 2) + pow(self.p3.y, 2) + pow(self.p3.z, 2), self.p3.y, self.p3.z, 1],
            [pow(self.p4.x, 2) + pow(self.p4.y, 2) + pow(self.p4.z, 2), self.p4.y, self.p4.z, 1]
        ]

        d_y = [
            [pow(self.p1.x, 2) + pow(self.p1.y, 2) + pow(self.p1.z, 2), self.p1.x, self.p1.z, 1],
            [pow(self.p2.x, 2) + pow(self.p2.y, 2) + pow(self.p2.z, 2), self.p2.x, self.p2.z, 1],
            [pow(self.p3.x, 2) + pow(self.p3.y, 2) + pow(self.p3.z, 2), self.p3.x, self.p3.z, 1],
            [pow(self.p4.x, 2) + pow(self.p4.y, 2) + pow(self.p4.z, 2), self.p4.x, self.p4.z, 1]
        ]

        d_z = [
            [pow(self.p1.x, 2) + pow(self.p1.y, 2) + pow(self.p1.z, 2), self.p1.x, self.p1.y, 1],
            [pow(self.p2.x, 2) + pow(self.p2.y, 2) + pow(self.p2.z, 2), self.p2.x, self.p2.y, 1],
            [pow(self.p3.x, 2) + pow(self.p3.y, 2) + pow(self.p3.z, 2), self.p3.x, self.p3.y, 1],
            [pow(self.p4.x, 2) + pow(self.p4.y, 2) + pow(self.p4.z, 2), self.p4.x, self.p4.y, 1]
        ]

        a = np.linalg.det(a)
        if a ==0: # a=0になる場合がある．おそらく四面体が潰れた形？
            flg = True
            return flg
        d_x = np.linalg.det(d_x)
        d_y = np.linalg.det(d_y)
        d_z = np.linalg.det(d_z)

        center_p_x = np.nan_to_num(d_x / (2 * a))
        center_p_y = np.nan_to_num(d_y / (2 * a))
        center_p_z = np.nan_to_num(d_z / (2 * a))

        self.center_p = [float(center_p_x), float(center_p_y), float(center_p_z)]

        distance = math.sqrt(
            pow((self.p1.x - self.center_p[0]), 2) + pow((self.p1.y - self.center_p[1]), 2) + pow((self.p1.z - self.center_p[2]), 2))
        
        # c = np.linalg.det([
        #     [pow(self.p1.x, 2) + pow(self.p1.y, 2) + pow(self.p1.z, 2), self.p1.x, self.p1.y, self.p1.z],
        #     [pow(self.p2.x, 2) + pow(self.p2.y, 2) + pow(self.p2.z, 2), self.p2.x, self.p2.y, self.p2.z],
        #     [pow(self.p3.x, 2) + pow(self.p3.y, 2) + pow(self.p3.z, 2), self.p3.x, self.p3.y, self.p3.z],
        #     [pow(self.p4.x, 2) + pow(self.p4.y, 2) + pow(self.p4.z, 2), self.p4.x, self.p4.y, self.p4.z]
        #     ])
        
        #r2 = math.sqrt(pow(d_x, 2) + pow(d_y, 2) + pow(d_z, 2) - 4 *a * c)/ (2*np.abs(a))
        
        #print('r1 = ', distance, "r2 = ", r2)

        self.radius = float(distance)
        flg = False
        return flg



    def check_point_include_circumsphere(self, check_p):
        '''外接球に任意の点が内包されているかどうか判定'''

        distance = math.sqrt(pow((check_p.x - self.center_p[0]), 2) + pow((check_p.y - self.center_p[1]), 2) + pow((check_p.z - self.center_p[2]), 2))

        if distance < self.radius:
            return True
        else:
            return False


    def check(self, select_point):
        node_list = []
        near90theta = 10000
        p0 = np.array([select_point.x, select_point.y, select_point.z])
        #print('p0 = ' + str(p0))
        p1 = np.array([self.p1.x, self.p1.y, self.p1.z])
        # print('p1 = ' + str(p1))
        p2 = np.array([self.p2.x, self.p2.y, self.p2.z])
        # print('p2 = ' + str(p2))
        p3 = np.array([self.p3.x, self.p3.y, self.p3.z])
        p4 = np.array([self.p4.x, self.p4.y, self.p4.z])
        loc1 = (p1+p2+p3)/3
        loc2 = (p1+p3+p4)/3
        loc3 = (p1+p4+p2)/3
        loc4 = (p2+p4+p3)/3
        #print('loc1 = ' + str(loc1))

        flg1, theta = check_in_out(p0, p1, p3, p2, loc1)
        if flg1:
            if theta < near90theta:
                near90theta = theta
                node_list = [self.p1.node, self.p3.node, self.p2.node]
        
        flg2, theta = check_in_out(p0, p1, p4, p3, loc2)
        if flg2:
            if theta < near90theta:
                near90theta = theta
                node_list = [self.p1.node, self.p4.node, self.p3.node]
        
        flg3, theta = check_in_out(p0, p1, p2, p4, loc3)
        if flg3:
            if theta < near90theta:
                near90theta = theta
                node_list = [self.p1.node, self.p2.node, self.p4.node]

        flg4, theta = check_in_out(p0, p2, p3, p4, loc4)
        if flg4:
            if theta < near90theta:
                near90theta = theta
                node_list = [self.p2.node, self.p3.node, self.p4.node]
        
        if any([flg1, flg2, flg3, flg4]) == True: # 1つでも四面体の外側のときFalseを返す
            check = False
            #print(node_list, near90theta)
            #time.sleep(3)            

        else: # 全て四面体の内側出会った場合，Trueを返す
            check = True
            #print('True')

        return check, node_list
    
    def output_stress(self, point):
        l1 = math.sqrt(pow(self.p1.x-point[0], 2) + pow(self.p1.y-point[1], 2) + pow(self.p1.z-point[2], 2))
        l2 = math.sqrt(pow(self.p2.x-point[0], 2) + pow(self.p2.y-point[1], 2) + pow(self.p2.z-point[2], 2))
        l3 = math.sqrt(pow(self.p3.x-point[0], 2) + pow(self.p3.y-point[1], 2) + pow(self.p3.z-point[2], 2))
        l4 = math.sqrt(pow(self.p4.x-point[0], 2) + pow(self.p4.y-point[1], 2) + pow(self.p4.z-point[2], 2))
        
        # 長さの割合
        u1 = l1/ (l1+l2+l3+l4)
        u2 = l2/ (l1+l2+l3+l4)
        u3 = l3/ (l1+l2+l3+l4)
        u4 = l4/ (l1+l2+l3+l4)        

        output_stress = (u2 + u3 + u4) * self.p1.stress + (u1 + u3 + u4) * self.p2.stress + (u2 + u1 + u4) * self.p3.stress + (u2 + u3 + u1) * self.p4.stress

        return output_stress
        
def check_in_out(t, pa, pb, pc, loc):
    u = pb - pa
    v = pc - pa
    #法線ベクトルw
    w = np.cross(u, v)
    #任意の点へのベクトルh
    h = t - loc
    #法線ベクトルwと平面から任意点のベクトルhの内積
    i = np.inner(w, h)
    n = LA.norm(w) * LA.norm(h)

    #内積から角度thetaを求める
    c = i / n
    theta = np.rad2deg(np.arccos(np.clip(c, -1.0, 1.0)))

    # 座標点が四面体の外側であればTrueを返す
    if theta > 90:
        flg = True
        #print('座標点は立体の外側')
    else:
        flg = False
        theta = -1
        #print('座標点は立体の内側')
    
    return flg, theta



def re_number(list):
    pa = np.array([list[0].x, list[0].y, list[0].z])
    pb = np.array([list[1].x, list[1].y, list[1].z])
    pc = np.array([list[2].x, list[2].y, list[2].z])
    pd = np.array([list[3].x, list[3].y, list[3].z])

    loc = (pb + pc + pd) / 3
    u = pc - pb
    v = pd - pb
    w = np.cross(u, v)
    h = pa - loc
    i = np.inner(w, h)
    n = LA.norm(w) * LA.norm(h)
    # if LA.norm(w)==0:
    #     print("w  [", pa, pb, pc, pd, "]")
    # if LA.norm(h)==0:
    #     print("h  [", pa, pb, pc, pd, "]")
    # if LA.norm(n)==0:
    #     print("n  [", pa, pb, pc, pd, "]")
    c = i / n
    a = np.rad2deg(np.arccos(np.clip(c, -1.0, 1.0)))
    #print("a = ", a)
    p_a = list[0]
    p_b = list[1]
    p_c = list[2]
    p_d = list[3]

    if a >= 90:
        p_c = list[3]
        p_d = list[2]
        #print("node change = ", pc.node, pd.node)
    
    return p_c, p_d

    

    
