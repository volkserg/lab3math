

class Solver:

    def __init__(self, a, b, g, infin=100):
        assert 0 < a < 1 and 0 < b < 1 and 0 < g < 1
        self.a = a
        self.b = b
        self.g = g
        self.inA = 1 - a
        self.inB = 1 - b
        self.inG = 1 - g
        self.q0 = [1]
        self.q1 = [a / (self.inA * b)]
        self.rho = a / b
        self.P1 = self.rho
        self.P0 = 1 - self.rho
        self.inf = infin

    # -------------------------1st-----------------

    def calc_q0(self, i):
        assert i >= 0
        if len(self.q0) <= i:
            # for k in range(len(self.q0), i):
            while len(self.q0) < i+1:
                self.q0.append(None)
        if self.q0[i]:
            return self.q0[i]
        else:
            self.q0[i] = self.a * self.inB * self.inG**i * self.calc_q1(i - 1) / (1 - self.inG**i)
            return self.q0[i]

    def calc_q1(self, i):
        assert i >= 0
        if len(self.q1) <= i:
            # for k in range(len(self.q1), i):
            while len(self.q1) < i+1:
                self.q1.append(None)
        if self.q1[i]:
            return self.q1[i]
        else:
            self.q1[i] = self.a * self.inB * (1 - self.inA * self.inG**i) * self.calc_q1(i-1) / (self.inA * self.b * (1 - self.inG**i))
            return self.q1[i]

    def calc_p00(self):
        summ1, summ2 = 0, 0
        for i in range(self.inf):
            summ1 += self.calc_q0(i)
        for i in range(self.inf):
            summ2 += self.calc_q1(i)
        return 1 / (summ1+summ2)

    def calc_p0i(self, i):
        return self.calc_p00() * self.calc_q0(i)

    def calc_p1i(self, i):
        return self.calc_p00() * self.calc_q1(i)

    def calc_p00_check(self):
        mult = 1
        for i in range(1, self.inf):
            mult *= self.inA * (self.b - self.a*self.inB*self.inG**i) / (self.inA*self.b - self.a*self.inB*self.inG**i)
        return (1 - self.rho) * (1/mult)

    def calc_pxi(self, i):
        if i == 0:
            return self.P0 + self.P1 * self.b
        else:
            return self.calc_p1i(i-1) * self.inB

    #--------------------2nd-------------------------

