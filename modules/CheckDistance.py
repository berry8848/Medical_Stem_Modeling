import math

class CheckDistance:
    '''応力値から点間距離を求め，点間距離内に他の点が含まれるか否かを判断する'''
    def __init__(self, COEFFICIENT_OF_LONG):
        self.COEFFICIENT_OF_LONG = COEFFICIENT_OF_LONG

    def check_distance(self, fixed_points, candidate_point, stress):
        density = stress_to_density(stress) # 応力から密度の変換
        long = self.COEFFICIENT_OF_LONG * density_to_long(density) # 密度から点間距離の変換

        # 点間距離内に他の点が含まれているか否かを判定．
        # 点間距離内に他の点が含まれていたらFalseを返す
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
    
# 応力と密度の関係式．※関係式が微妙なため，臨時で別の関数
def stress_to_density(stress):
    if stress >= 0: #生データが正の場合，0を返す．（通常はマイナスの値をとる）
        density = 1.448
    else:
        #density = -stress/15
        density = 3*0.00001*(stress)*(stress) + 0.01*(stress) + 1.448
    return density

# 密度と点間距離の関係式
def density_to_long(density):
    long = 6 * 0.1 * math.sqrt(math.sqrt(2)*math.pi/ density)
    return long
