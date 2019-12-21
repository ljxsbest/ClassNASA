import numpy as np
import matplotlib.pyplot as plt
from sympy import *

class NASA:

    # [Name, tlow, thigh, tmid, a[1-14]]
    spec1 = ["", 0, 0, 0, []]
    spec2 = ["", 0, 0, 0, []]
    spec3 = ["", 0, 0, 0, []]
    Tlow = 0
    Thigh = 0
    fraction = []
    R = 8.314
    
    def nasa_coff(self, species):  # 此函数用于查找文件中的species
        with open("thermo30.dat", 'r') as filein:
            text = filein.readlines()

        data = []
        spe_flag = 0
        # search the given species
        for index, line in enumerate(text):
            if line.startswith(species + " "):
            # found
                spe_flag = 1
                data.append(text[index])
                data.append(text[index + 1])
                data.append(text[index + 2])
                data.append(text[index + 3])  # 将species对应的数据交给 data
                break

        # if not found
        if spe_flag == 0:
            print("Species %s not found." % (species))
            return False
        # 读取每个species后面的温度，因为有的species不到5000k就分解完了。
        Tlow, Thi, Tbk = float(data[0].split()[-4]), float(data[0].split()[-3]), float(
        data[0].split()[-2])
        # 读species的系数a1至a14
        a = []
        for i in range(14):
            a.append(float(data[int(i / 5) + 1][15 * (i % 5):15 * (i % 5 + 1)]))
        # 在这一行中，[int(i/5) + 1]用来找到需要的数字在的行；[15*(i%5):15*(i%5 + 1)] 用来找到数字所在的列区间。
        return [Tlow, Thi, Tbk, a]

    # 计算熵
    def nasa_S(self, data, T):
        R = self.R
        Tdata = [data[0], data[1], data[2]]
        coff = data[3]
        print(Tdata, T)
        if T >= Tdata[2] and T <= Tdata[1]:
            S = R * (coff[0] * np.log(T) + coff[1] * T + (coff[2] * T ** 2) / 2 + (coff[3] * T ** 3) / 3 + (
                    coff[4] * T ** 4) / 4 + coff[6])
        elif T < Tdata[2] and T >= Tdata[0]:
            S = R * (coff[7] * np.log(T) + coff[8] * T + (coff[9] * T ** 2) / 2 + (coff[10] * T ** 3) / 3 + (
                    coff[11] * T ** 4) / 4 + coff[13])

        else:
            print("Temperature %s K is not supported." % (T))
        return False
        return S

    # 计算焓
    def nasa_H(self, data, T):
        R = self.R
        Tdata = [data[0], data[1], data[2]]
        coff = data[3]
        if T >= Tdata[2] and T <= Tdata[1]:
            H = R * T * (coff[0] + (coff[1] * T) / 2 + (coff[2] * T ** 2) / 3 + (coff[3] * T ** 3) / 4 + (coff[4] * T ** 4) / 5 + coff[5] / T)
        elif T < Tdata[2] and T >= Tdata[0]:
            H = R * T * (coff[7] + (coff[8] * T) / 2 + (coff[9] * T ** 2) / 3 + (coff[10] * T ** 3) / 4 + (coff[11] * T ** 4) / 5 + coff[12] / T)
        else:
            print("Temperature %s K is not supported." % (T))
            return False
        return H

    # 计算摩尔分数
    def fraction(self, Kp):
        n = symbols('n')

        Nmole = solve(2 * n * n / Kp + n - 1, n)
        if Nmole[0] > 0:
            N = Nmole[0]
        else:
            N = Nmole[1]
            N2 = (1 - N) / 2
            return N2 / (N2 + N)
    
    def set_temp_range(self, l1, l2, l3, h1, h2, h3):
        self.Tlow = max(l1, l2, l3)
        self.Thigh = min(h1, h2, h3)

    def __init__(self, sp1name, sp2name, sp3name):
    
        self.spec1[0] = sp1name
        print(sp1name)
        self.spec1[1] = self.nasa_coff(sp1name)[0]
        self.spec1[2] = self.nasa_coff(sp1name)[1]
        self.spec1[3] = self.nasa_coff(sp1name)[2]
        self.spec1[4] = self.nasa_coff(sp1name)[3]

        self.spec2[0] = sp2name
        self.spec2[1] = self.nasa_coff(sp2name)[0]
        self.spec2[2] = self.nasa_coff(sp2name)[1]
        self.spec2[3] = self.nasa_coff(sp2name)[2]
        self.spec2[4] = self.nasa_coff(sp2name)[3]

        self.spec3[0] = sp3name
        self.spec3[1] = self.nasa_coff(sp3name)[0]
        self.spec3[2] = self.nasa_coff(sp3name)[1]
        self.spec3[3] = self.nasa_coff(sp3name)[2]
        self.spec3[4] = self.nasa_coff(sp3name)[3]

        print(self.spec1[1], self.spec2[1], self.spec3[1], self.spec1[2], self.spec2[2], self.spec3[2])
        # Determine the computing range of temperature
        self.set_temp_range(self.spec1[1], self.spec2[1], self.spec3[1], self.spec1[2], self.spec2[2], self.spec3[2])
        print(self.Tlow, self.Thigh)

    def compute(self):
        R = self.R
        self.fraction = []
        for t in range(int(self.Tlow), int(self.Thigh)+1):
            deltaS = self.nasa_S(self.spec2[4], t) + self.nasa_S(self.spec3[4], t) - self.nasa_S(self.spec1[4], t)
            deltaH = self.nasa_H(self.spec2[4], t) + self.nasa_H(self.spec3[4], t) - self.nasa_H(self.spec1[4], t)
            Kp = np.exp(-deltaH / (R * t)) * np.exp(deltaS / R)
            (self.fraction).append(self.fraction(Kp))
            print(t, frac[-1])

    def plotim(self):
        plt.plot(range(self.Tlow, self.Thigh+1), frac)
        plt.savefig("molefraction.png", dpi=800)

test = NASA("CH4", "CH3", "H")
test.compute()
test.plotim()
