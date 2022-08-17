'''
计算高熵合金的Warren-Cowley参数
输入：POSCAR、目标原子对
'''

import numpy as np

# 读取POSCAR文件并处理数据
def read_POSCAR(POSCAR):
    with open(POSCAR, 'r', encoding='utf-8') as f:
        POSCAR = f.readlines()
        elements = POSCAR[5].strip().split()                    # 元素名称
        element_number = POSCAR[6].strip().split()              #
        element_number = [int(x) for x in element_number]       # 每种元素个数
        elements = dict(zip(elements, element_number))
        # element_list = []                                      # 储存每个原子
        # for i in range(len(elements)):
        #     for j in range(element_number[i]):
        #         element_list.append(elements[i] + str(j+1))
        # 读取坐标信息
        if POSCAR[7].strip() == "Selective dynamics":
            coordiate = POSCAR[9:]
        else:
            coordiate = POSCAR[8:]
        coor = []                                               # 储存每个原子的坐标
        for i in coordiate:
            one_coor = i.strip().split()[:3]
            if one_coor:
                one_coor = [float(x) for x in one_coor]
                coor.append(one_coor)
        # 生成每种元素对应的坐标字典
        elements_coor = []          # 每种元素的坐标为一个小列表，储存在一个大列表中，最后合并为一个大字典
        index = 0
        for i in element_number:
            one_element_coor = coor[index:index+i]
            index += i
            elements_coor.append(one_element_coor)
        atoms = dict(zip(elements, elements_coor))

    return atoms, elements, coor


# 计算两点之间的距离，考虑周期性边界条件
def distance(i, j, scalor):
    dx = (i[0] - j[0])*scalor
    if (abs(dx) > scalor*0.5):
        dx = scalor - abs(dx)
    dy = (i[1] - j[1])*scalor
    if (abs(dy) > scalor*0.5):
        dy = scalor - abs(dy)
    dz = (i[2] - j[2])*scalor
    if (abs(dz) > scalor*0.5):
        dz = scalor - abs(dz)
    dis = np.sqrt(dx**2 + dy**2 + dz**2)
    return dis

# 计算第一近邻距离，只需将第一个原子与剩余原子的距离算一遍，取最小值即可
def first_neigh(coor, scalor):
    first = coor[0]
    distance_list = []
    for i in coor[1:]:
        dis = distance(i, first, scalor=12)
        distance_list.append(dis)
    first_neighboor = min(distance_list)
    return first_neighboor

# 计算warren-cowley参数，一次性计算出所有原子对的SRO参数
def cacu_SRO(elements, first_neighboor):
    SRO = []                            # 储存计算出的SRO参数，是按顺序排列的，之后转换为n*n矩阵即可
    for element_i in elements:
        Zi = 0
        for element_k in elements:
            for i in atoms[element_i]:
                for j in atoms[element_k]:
                    if distance(i, j, scalor=12) == first_neighboor:
                        Zi += 1
        for element_j in elements:
            Zij = 0
            for i in atoms[element_i]:
                for j in atoms[element_j]:
                    if distance(i, j, scalor=12) == first_neighboor:
                        # i原子对应的j近邻原子数加1
                        Zij += 1
            # 计算该原子对的SRO，需要该原子数，原子分数
            fraction_j = elements[element_j] / sum(elements.values())
            WCP_ij = 1 - Zij / (Zi * fraction_j)
            SRO.append(WCP_ij)
            # print(Zij, Zi, fraction_j)
    return np.array(SRO).reshape(4, 4)


atoms, elements, coor = read_POSCAR('POSCAR-5')
first_neighboor = first_neigh(coor, scalor=12)
SRO = cacu_SRO(elements, first_neighboor)
print(SRO)

# print(first_neighboor)
# print(atoms)
# print(elements)

atoms, elements, coor = read_POSCAR('best1-POSCAR')
first_neighboor = first_neigh(coor, scalor=12)
SRO = cacu_SRO(elements, first_neighboor)
print(SRO)

atoms, elements, coor = read_POSCAR('best2-POSCAR')
first_neighboor = first_neigh(coor, scalor=12)
SRO = cacu_SRO(elements, first_neighboor)
print(SRO)