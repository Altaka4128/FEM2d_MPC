import numpy as np

class Dmatrix2d:

    def __init__(self, young, poisson):
        self.young = young       # ヤング率
        self.poisson = poisson   # ポアソン比
    
    # 平面応力状態のDマトリクスを作成する
    def makeDcpsmatrix(self):
        coef = self.young / (1 - self.poisson * self.poisson)
        matD = np.matrix([[1.0, self.poisson, 0.0],
                          [self.poisson, 1.0, 0.0],
                          [0.0, 0.0, 0.5 * (1 - self.poisson)]])
        matD *= coef

        return matD