import numpy as np
import numpy.linalg as LA

class FEM2d:

    def __init__(self, nodes, elements, bound):
        self.nodes = nodes         # 節点は1から始まる順番で並んでいる前提(Node2d型のリスト)
        self.elements = elements   # 要素のリスト(d2cps型のリスト)
        self.bound = bound         # 境界条件(Boundary2d型)

    def analysis(self):
        # 境界条件を考慮しないKマトリクスを作成する
        matK = self.makeKmatrix()

        # 荷重ベクトルを作成する
        vecf = self.bound.makeForceVector()

        # 境界条件を考慮したKマトリクス、荷重ベクトルを作成する
        matKc, vecfc, mpcNum = self.setBoundCondition(matK, vecf)

        # 変位ベクトルを計算する
        vecTmp = LA.solve(matKc, vecfc)
        vecDisp = np.delete(vecTmp, slice(len(vecTmp) - mpcNum, len(vecTmp)))
        self.vecDisp = vecDisp

        return vecDisp

    # 境界条件を考慮しないKマトリクスを作成する
    def makeKmatrix(self):

        matK = np.matrix(np.zeros((len(self.nodes) * 2, len(self.nodes) * 2)))
        for elem in self.elements:
            # keマトリクスを計算する
            matKe = elem.makeKematrix()

            # Kマトリクスに代入する
            for c in range(len(elem.nodes) * 2):
                ct = (elem.nodes[c // 2].no - 1) * 2 + c % 2
                for r in range(len(elem.nodes) * 2):
                    rt = (elem.nodes[r // 2].no - 1) * 2 + r % 2
                    matK[ct, rt] += matKe[c, r]
        return matK

    # Kマトリクス、荷重ベクトルに境界条件を考慮する
    def setBoundCondition(self, matK, vecf):
        matKc = np.copy(matK)
        vecfc = np.copy(vecf)

        # 単点拘束条件を考慮したKマトリクス、荷重ベクトルを作成する
        vecDisp = self.bound.makeDispVector()
        for i in range(len(vecDisp)):
            if not vecDisp[i] == None:
                # Kマトリクスからi列を抽出する
                vecx = np.array(matK[:, i]).flatten()

                # 変位ベクトルi列の影響を荷重ベクトルに適用する
                vecfc = vecfc - vecDisp[i] * vecx

                # Kマトリクスのi行、i列を全て0にし、i行i列の値を1にする
                matKc[:, i] = 0.0
                matKc[i, :] = 0.0
                matKc[i, i] = 1.0
        for i in range(len(vecDisp)):
            if not vecDisp[i] == None:
                vecfc[i] = vecDisp[i]

        # 多節点拘束条件を考慮したKマトリクス、荷重ベクトルを作成する
        matC, vecd = self.bound.makeMPCmatrixes()
        mpcNum = len(vecd)   # 多節点拘束の条件式の数を求める
        if mpcNum > 0:
            matKc = np.hstack((matKc, matC.T))
            tmpMat = np.hstack((matC, np.zeros((matC.shape[0], matC.shape[0]))))
            matKc = np.vstack((matKc, tmpMat))
            vecfc = np.hstack((vecf, vecd))

        return matKc, vecfc, mpcNum