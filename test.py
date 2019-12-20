import numpy as np
import matplotlib.pyplot as plt
from sympy import *


class NASA:
    def nasa_coff(species):  # 此函数用于查找文件中的species
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
            Tlow(species), Thi(species), Tbk(species) = float(data[0].split()[-4]), float(data[0].split()[-3]), float(
            data[0].split()[-2])
    # 读species的系数a1至a14
        a = []
        for i in range(14):
            a.append(float(data[int(i / 5) + 1][15 * (i % 5):15 * (i % 5 + 1)]))
    # 在这一行中，[int(i/5) + 1]用来找到需要的数字在的行；[15*(i%5):15*(i%5 + 1)] 用来找到数字所在的列区间。
        return [Tlow.species, Thi.species, Tbk.species, a]


    def nasa_S(data, T):
        Tdata, coff = data
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


    def nasa_H(data, T):
        Tdata, coff = data
        if T >= Tdata[2] and T <= Tdata[1]:
            H = R * T * (coff[0] + (coff[1] * T) / 2 + (coff[2] * T ** 2) / 3 + (coff[3] * T ** 3) / 4 + (
                coff[4] * T ** 4) / 5 + coff[5] / T)
        elif T < Tdata[2] and T >= Tdata[0]:
            H = R * T * (coff[7] + (coff[8] * T) / 2 + (coff[9] * T ** 2) / 3 + (coff[10] * T ** 3) / 4 + (
                coff[11] * T ** 4) / 5 + coff[12] / T)
        else:
            print("Temperature %s K is not supported." % (T))
            return False
        return H


    def fraction(Kp):
        n = symbols('n')


        Nmole = solve(2 * n * n / Kp + n - 1, n)
        if Nmole[0] > 0:
            N = Nmole[0]
        else:
            N = Nmole[1]
            N2 = (1 - N) / 2
            return N2 / (N2 + N)


    def input():
        species1 = input()
        species2 = input()
        species3 = input()
        if min(Thi.species1, Thi.species2, Thi.species3) == 3500:
            Temp = range(340, 3500)
        else:
            Temp = range(340, 3600)


    def output():
        data1 = nasa_coff(species1)
        data2 = nasa_coff(species2)
        data3 = nasa_coff(species3)
        frac = []
        for t in Temp:
            deltaS = nasa_S(data2, t) + nasa_S(data3, t) - nasa_S(data1, t)
        deltaH = nasa_H(data2, t) + nasa_H(data3, t) - nasa_H(data1, t)
        Kp = np.exp(-deltaH / (R * t)) * np.exp(deltaS / R)
        frac.append(fraction(Kp))
        print(t, frac[-1])

    def plot():
        plt.plot(Temp, frac)
        plt.savefig("molefraction.png", dpi=800)
        plt.show()


NASA.input()
NASA.output()
NASA.plot()