

class Solver:

    def __init__(self, a, b, g):
        assert 0 < a < 1 and 0 < b < 1 and 0 < g < 1
        self.a = a
        self.b = b
        self.g = g
        self.inA = 1 - a
        self.inB = 1 - b
        self.inG = 1 - g
        self.q0 = [1]
        self.q1 = [a / (self.inA * b)]

    def calc_q0(self, i, flag=False):
        assert i >= 0
        if not flag:
            for k in range(i):
                self.q0.append(None)
        if self.q0[i]:
            return self.q0[i]
        else:
            self.q0[i] = self.a * self.inB * self.inG**i * self.calc_q0(i - 1, True) / (1 - self.inG**i)
            return self.q0[i]

    def calc_q1(self, i):
        assert i >= 0
        if i == 0:
            return self.a * self.calc_q0(0) / (self.inA * self.b)
        else:
            return self.a * self.inB * (1 - self.inA * self.inG**i) * self.calc_q1(i-1) / (self.inA * self.inB * (1 - self.inG**i))

    def calc_p00(self):
        summ1, summ2 = 0, 0
        for i in range(1000):
            summ1 += self.calc_q0(i)
        for i in range(1000):
            summ2 += self.calc_q1(i)
        return 1 / (summ1+summ2)
