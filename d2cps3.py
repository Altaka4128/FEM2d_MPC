import numpy as np
import numpy.linalg as LA
from Dmatrix2d import Dmatrix2d

class d2cps3:

    def __init__(self, no, nodes, thickness, young, poisson):
        self.no = no                 # 要素番号
        self.nodes = nodes           # nodesは反時計回りの順番になっている前提(Node2d型のリスト形式)
        self.thickness = thickness   # 厚さ
        self.young = young           # ヤング率
        self.poisson = poisson       # ポアソン比
        self.ipNum = 1               # 積分点数

        # Dマトリクスを計算する
        self.matD = self.makeDmatrix()

        # ヤコビ行列を計算する
        self.matJ = self.makeJmatrix()

        # Bマトリクスを計算する
        self.matB = self.makeBmatrix()
    
    # Keマトリクスを作成する
    def makeKematrix(self):
        
        w = 0.5   # 積分点の重み係数
        matKe = w * self.thickness * self.matB.T * self.matD * self.matB * LA.det(self.matJ)

        return matKe
    
    # ヤコビ行列を作成する
    def makeJmatrix(self):

        dxda = self.nodes[0].x - self.nodes[2].x
        dxdb = self.nodes[1].x - self.nodes[2].x
        dyda = self.nodes[0].y - self.nodes[2].y
        dydb = self.nodes[1].y - self.nodes[2].y
        matJ = np.matrix([[dxda, dyda],
                          [dxdb, dydb]])
        
        return matJ

    # Dマトリクスを作成する
    def makeDmatrix(self):
        d2dmat = Dmatrix2d(self.young, self.poisson)
        matD = d2dmat.makeDcpsmatrix()

        return matD

    # Bマトリクスを作成する
    def makeBmatrix(self):

        # dNi/da, dNi/dbを計算する
        dN1da = 1.0
        dN1db = 0.0
        dN2da = 0.0
        dN2db = 1.0
        dN3da = -1.0
        dN3db = -1.0

        # dNi/dx, dNi/dyを計算する
        matdNdab = np.matrix([[dN1da, dN2da, dN3da],
                              [dN1db, dN2db, dN3db]])
        #dNdxy = matJinv * matdNdab
        dNdxy = LA.solve(self.matJ, matdNdab)
        
        # Bマトリクスを計算する
        matB = np.matrix([[dNdxy[0, 0], 0.0, dNdxy[0, 1], 0.0, dNdxy[0, 2], 0.0],
                          [0.0, dNdxy[1, 0], 0.0, dNdxy[1, 1], 0.0, dNdxy[1, 2]],
                          [dNdxy[1, 0], dNdxy[0, 0], dNdxy[1, 1], dNdxy[0, 1], dNdxy[1, 2], dNdxy[0, 2]]])
        
        return matB
