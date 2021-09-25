from Node2d import Node2d
from d2cps3 import d2cps3
from Boundary2d import Boundary2d
from FEM2d import FEM2d
import numpy as np

if __name__ == '__main__':
    # 節点のリストを作成する
    node1 = Node2d(1, 0.0, 0.0)
    node2 = Node2d(2, 1.0, 0.0)
    node3 = Node2d(3, 0.0, 1.0)
    node4 = Node2d(4, 1.0, 1.0)
    nodes1 = [node1, node2, node4]
    nodes2 = [node1, node4, node3]
    nodes = [node1, node2, node3, node4]

    # 要素のリストを作成する
    thickness = 1
    young = 210000
    poisson = 0.3
    elem1 = d2cps3(1, nodes1, thickness, young, poisson)
    elem2 = d2cps3(2, nodes2, thickness, young, poisson)
    elems = [elem1, elem2]

    # 境界条件を設定する
    theta = 45
    mpc1 = np.array([0, 0, -np.sin(np.deg2rad(theta)), np.cos(np.deg2rad(theta)), 0, 0, 0, 0])
    d1 = 0
    bound = Boundary2d(len(nodes))
    bound.addSPC(1, 0.0, 0.0)
    bound.addSPC(3, 0.0, 0.0)
    bound.addMPC(mpc1, d1)
    bound.addForce(4, 0.0, 1000.0)

    # 有限要素法で変位を計算する
    fem2d = FEM2d(nodes, elems, bound)
    vecDisp = fem2d.analysis()
    print("変位ベクトル : ", vecDisp)