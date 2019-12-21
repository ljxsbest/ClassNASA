import numpy as np
import matplotlib.pyplot as plt
from sympy import *

class NASA():

    ####################### Overall definitions ######################
    R = 8.314
    # Molefraction
    def frac(Kp):
        n = symbols('n')

        Pmole = solve(n * n + Kp * n - Kp, n)
        if Pmole[0] > 0:
            Product = Pmole[0]
        else:
            Product = Pmole[1]
            Origin = 1 - Product
            return Origin / (Origin + 2*Product)
    # Entropy
    def nasa_S(data, T):
        R = NASA.R
        Tdata = data[1:4]
        coff = data[4]
        if T >= Tdata[2] and T <= Tdata[1]:
            S = R * (coff[0] * np.log(T) + coff[1] * T + (coff[2] * T ** 2) / 2 + (coff[3] * T ** 3) / 3 + (coff[4] * T ** 4) / 4 + coff[6])
        elif T < Tdata[2] and T >= Tdata[0]:
            S = R * (coff[7] * np.log(T) + coff[8] * T + (coff[9] * T ** 2) / 2 + (coff[10] * T ** 3) / 3 + (coff[11] * T ** 4) / 4 + coff[13])
        else:
            print("Temperature %s K is not supported." % (T))
            return False
        return S
    # Enthapy
    def nasa_H(data, T):
        R = NASA.R
        Tdata = data[1:4]
        coff = data[4]
        if T >= Tdata[2] and T <= Tdata[1]:
            H = R * T * (coff[0] + (coff[1] * T) / 2 + (coff[2] * T ** 2) / 3 + (coff[3] * T ** 3) / 4 + (coff[4] * T ** 4) / 5 + coff[5] / T)
        elif T < Tdata[2] and T >= Tdata[0]:
            H = R * T * (coff[7] + (coff[8] * T) / 2 + (coff[9] * T ** 2) / 3 + (coff[10] * T ** 3) / 4 + (coff[11] * T ** 4) / 5 + coff[12] / T)
        else:
            print("Temperature %s K is not supported." % (T))
            return False
        return H

    #################### Object Definitions ############################

    # Tool functionality
    def set_temp_range(self, l1, l2, l3, h1, h2, h3):
        self.Tlow = max(l1, l2, l3)
        self.Thigh = min(h1, h2, h3)
    
    # spec1 => spec2 + spec3
    def __init__(self, spec1, spec2, spec3):

        self.tempout = []
        self.fraction = []
        
        # [Name, tlow, thigh, tmid, a[1-14]]
        self.spec1 = ["", 0, 0, 0, []]
        self.spec2 = ["", 0, 0, 0, []]
        self.spec3 = ["", 0, 0, 0, []]
        self.spec1[0] = spec1
        self.spec2[0] = spec2
        self.spec3[0] = spec3

        with open("thermo30.dat", 'r') as filein:
            self.thermo30 = filein.readlines()

        self.spec1[1] = self.get_coff(spec1)[0][0]
        self.spec1[2] = self.get_coff(spec1)[0][1]
        self.spec1[3] = self.get_coff(spec1)[0][2]
        self.spec1[4] = self.get_coff(spec1)[1][:]

        self.spec2[1] = self.get_coff(spec2)[0][0]
        self.spec2[2] = self.get_coff(spec2)[0][1]
        self.spec2[3] = self.get_coff(spec2)[0][2]
        self.spec2[4] = self.get_coff(spec2)[1][:]

        self.spec3[1] = self.get_coff(spec3)[0][0]
        self.spec3[2] = self.get_coff(spec3)[0][1]
        self.spec3[3] = self.get_coff(spec3)[0][2]
        self.spec3[4] = self.get_coff(spec3)[1][:]

        self.set_temp_range(self.spec1[1], self.spec2[1], self.spec3[1], self.spec1[2], self.spec2[2], self.spec3[2])

    def get_coff(self, species):
        spe_flag = 0
        data = []
        # search the given species
        for index, line in enumerate(self.thermo30):
            if line.startswith(species + " "):
            # found
                spe_flag = 1
                data.append(self.thermo30[index])
                data.append(self.thermo30[index + 1])
                data.append(self.thermo30[index + 2])
                data.append(self.thermo30[index + 3])  # 将species对应的数据交给 data
                break

        # if not found
        if spe_flag == 0:
            print("Species %s not found." % (species))
            return False
        # 读取每个species后面的温度，因为有的species不到5000k就分解完了。
        Tlow, Thi, Tbk = float(data[0].split()[-4]), float(data[0].split()[-3]), float(data[0].split()[-2])
        # 读species的系数a1至a14
        a = []
        for i in range(14):
            a.append(float(data[int(i / 5) + 1][15 * (i % 5):15 * (i % 5 + 1)]))
        # 在这一行中，[int(i/5) + 1]用来找到需要的数字在的行；[15*(i%5):15*(i%5 + 1)] 用来找到数字所在的列区间。
        return [Tlow, Thi, Tbk], a

    def print_rxn(self):
        print("%s ==> %s + %s" %(self.spec1[0], self.spec2[0], self.spec3[0]))
    
    def compute(self):

        self.print_rxn()
        
        R = NASA.R

        tempout = []
        fraction = []

        for t in range(int(self.Tlow), int(self.Thigh)+1):
            if(t % 100 == 0):
                tempout.append(t)
                deltaS = NASA.nasa_S(self.spec2, t) + NASA.nasa_S(self.spec3, t) - NASA.nasa_S(self.spec1, t)
                deltaH = NASA.nasa_H(self.spec2, t) + NASA.nasa_H(self.spec3, t) - NASA.nasa_H(self.spec1, t)
                Kp = np.exp(-deltaH / (R * t)) * np.exp(deltaS / R)
                fraction.append(NASA.frac(Kp))
                if(t % 500 == 0):
                    print("T = %s K,\tMole fraction of %s = %s." %(t, self.spec1[0], fraction[-1]))
        self.tempout = np.array(tempout[:])
        self.fraction = np.array(fraction[:])

    def plotimg(self, filename, dpi):
        plt.plot(self.tempout, self.fraction, "+", c="r")
        plt.savefig(filename, dpi=dpi)
        plt.show()
        
rxn1 = NASA("O2", "O", "O")
rxn1.compute()
rxn2 = NASA("N2", "N", "N")
rxn2.compute()
rxn3 = NASA("CH4", "CH3", "H")
rxn3.compute()
rxn4 = NASA("H2O", "OH", "H")
rxn4.compute()
rxn5 = NASA("C2H2", "C2H", "H")
rxn5.compute()
plt.plot(rxn1.tempout, rxn1.fraction, "*")
plt.plot(rxn2.tempout, rxn2.fraction, "*")
plt.plot(rxn3.tempout, rxn3.fraction, "*")
plt.plot(rxn4.tempout, rxn4.fraction, "*")
plt.plot(rxn5.tempout, rxn5.fraction, "*")
plt.show()
