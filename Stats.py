import numpy as np
from matplotlib import pyplot as plt
# from math import gamma as gamma_func
# Because of:
    # gamma_func dies with "OverflowError: math range error" which is because of:
        #64 bit float used in math library
from mpmath import gamma as gamma_func
from mpmath import mp, isnan
def retriveData(data: list):
    Dpts = []
    Dips = []
    Azi = []
    for item in data:
        Dpts.append(item[0])
        Dips.append(item[1])
        Azi.append(item[2])
    return {'Dpts': Dpts, 'Dips': Dips, 'Azi': Azi}
class CalcStats():
    def __init__(self, dataArr: list, stepH: int, showPrint: bool = False):
        """'dataArr' - list of floats. 'stepH' - number of histogram binss. 'showPrint' - show the calculation step-by-step."""
        self.maxVal = max(dataArr)
        self.minVal = min(dataArr)
        if self.maxVal == self.minVal:
            #Skips over data of same values and draws a rectangle.
            self.mean = self.maxVal
            self.mode = self.maxVal
            self.curveType = 'Skipped due to same values'
            self.x = self.maxVal
            self.y = 100
            self.bins2 = [self.maxVal]
            self.freq = [50]
            self.binWidth = 10
        else:
            #Normal calculation
            self.step = (self.maxVal - self.minVal)/stepH
            self.binsold = np.linspace(self.minVal, self.maxVal, stepH+1)
            self.bins2 = np.linspace(self.minVal, self.maxVal, stepH) #+1 # fucking madness
            self.binWidth = self.bins2[1] - self.bins2[0]
            self.bins = self.bins2 + self.binWidth / 2
            self.freq, _ = np.histogram(dataArr, bins=self.binsold)
            self.binMids = np.ones((stepH), dtype=np.float64)
            i = 0
            while i < len(self.bins)-1:
                self.binMids[i] = self.minVal+self.step*(i-0.5)
                i+=1
            # self.moments = self.binMids*self.freq
            self.moments = self.bins*self.freq
            self.mean = sum(self.moments)/sum(self.freq)
            self.zetta = self.bins/self.step
            if showPrint:
                print(f'Max value = {self.maxVal}')
                print(f'Min value = {self.minVal}')
                print(f'Step = {self.step}')
                print(f'Bins = {self.bins}')
                # print(f'middle of bins = {self.binMids}')
                print(f'Frequency = {self.freq}')
                print(f'Moments = {self.moments}')
                print(f'Zetta = {self.zetta}')
            self.calcNumWhatever(showPrint=showPrint)
            self.calcMu(showPrint=showPrint)
            self.calcKappa(showPrint=showPrint)
            _ = self.determineCurveType()
            print(_)
            if _ != 0:
                self.delta = self.step
                self.x = np.arange(self.minVal, self.maxVal, ((self.maxVal-self.minVal)/1000))
                self.y = self.frequencyCurveCalculation(showPrint=showPrint)
            else: pass
  
    def calcNumWhatever(self, showPrint: bool = False):
        self.calc1 = self.freq*self.zetta**1
        self.calc2 = self.freq*self.zetta**2
        self.calc3 = self.freq*self.zetta**3
        self.calc4 = self.freq*self.zetta**4
        if showPrint:
            print(f"frequency*x':\n{self.calc1}")
            print(f"frequency*x'^2:\n{self.calc2}")
            print(f"frequency*x'^3:\n{self.calc3}")
            print(f"frequency*x'^4:\n{self.calc4}")
            print(f"sum(frequency*x) = {sum(self.calc1)/sum(self.freq)}")
            print(f"sum(frequency*x'^2) = {sum(self.calc2)/sum(self.freq)}")
            print(f"sum(frequency*x'^3) = {sum(self.calc3)/sum(self.freq)}")
            print(f"sum(frequency*x'^4) = {sum(self.calc4)/sum(self.freq)}")
    def calcMu(self,showPrint: bool =False):
        self.sumfreq = sum(self.freq)
        self.mu1 = 0
        self.mu2 = (sum(self.calc2)/self.sumfreq) - (sum(self.calc1)/self.sumfreq)**2
        self.mu3 = (sum(self.calc3)/self.sumfreq) - 3*(sum(self.calc2)/self.sumfreq)*(sum(self.calc1)/self.sumfreq) + 2*(sum(self.calc1)/self.sumfreq)**3 
        self.mu4 = (sum(self.calc4)/self.sumfreq) - 4*(sum(self.calc1)/self.sumfreq)*(sum(self.calc3)/self.sumfreq) + 6*(sum(self.calc1)/self.sumfreq)**2 * (sum(self.calc2)/self.sumfreq) - 3*(sum(self.calc1)/self.sumfreq)**4
        # tmp WHY IS IT ZERO?
        # self.mu4 = 0
        if showPrint:
            print(f'Mu1 = {self.mu1}')
            print(f'Mu2 = {self.mu2}')
            print(f'Mu3 = {self.mu3}')
            print(f'Mu4 = {self.mu4}')
    def calcKappa(self, showPrint: bool =False):
        self.beta1 = (self.mu3**2) / (self.mu2**3)
        self.beta2 = (self.mu4) / (self.mu2**2)
        self.kappa = (self.beta1*(self.beta2+3)**2) / (4*(2*self.beta2 - 3*self.beta1 - 6)*(4*self.beta2 - 3*self.beta1))
        if showPrint:
            print(f'beta1 = {self.beta1}')
            print(f'beta2 = {self.beta2}')
            print(f'kappa = {self.kappa}')
    def determineCurveType(self):
        if self.kappa < 0: self.curveType = 'I'
        elif self.kappa > 0 and self.kappa < 1: self.curveType = 'IV' 
        elif self.kappa > 1: self.curveType = 'VI'
        elif self.kappa == 0 or self.kappa == 1: 
            print(f"Kappa = {self.kappa}. Curve type not implemented.")
            pass #TODO implement other curve types !!!
        #I dont remember why this try-except part is here, I must have been testing something
        try: 
            # print(f'Curve Type = {self.curveType}')
            pass
        except: return 0
    def frequencyCurveCalculation(self, showPrint: bool =False):
        if self.curveType == 'I':
            self.r = (6*(self.beta2-self.beta1-1)) / (3*self.beta1 - 2*self.beta2 + 6)
            self.a1a2 = 0.5*(self.mu2*(self.beta1*((self.r+2)**2)+16*(self.r+1)))**(1/2)
            self.root1 = 0.5*(self.r-2)+0.5*self.r*(self.r+2)*(self.beta1/(self.beta1*((self.r+2)**2)+16*(self.r+1)))**(1/2)
            self.root2 = 0.5*(self.r-2)-0.5*self.r*(self.r+2)*(self.beta1/(self.beta1*((self.r+2)**2)+16*(self.r+1)))**(1/2)
            if self.mu3 < 0:
                self.m1 = self.root1
                self.m2 = self.root2
            else: 
                self.m1 = self.root2
                self.m2 = self.root1
            self.a1 = self.a1a2*self.m1/(self.m1+self.m2)
            self.a2 = self.a1a2*self.m2/(self.m1+self.m2)
            self.A1 = self.a1a2*(self.m1+1)/(self.m1+self.m2+2)
            self.A2 = self.a1a2*(self.m2+1)/(self.m1+self.m2+2)
            tmparr = [self.m1, self.m2]
            self.ye = (self.sumfreq/(self.A1+self.A2)) * ((self.m1+1)**self.m1) * ((self.m2+1)**self.m2) * gamma_func(self.m1+self.m2+2) / (((self.m1+self.m2+2)**(self.m1+self.m2)*gamma_func(self.m1+1)*gamma_func(self.m2+1)))
            self.y0 = self.sumfreq*(np.sign(self.m1) * (np.abs(self.m1)) ** self.m1)*(np.sign(self.m2) * (np.abs(self.m2)) ** self.m2)*gamma_func(self.m1+self.m2+2)/(self.a1a2*(np.sign(self.m1 + self.m2) * (np.abs(self.m1 + self.m2)) ** (self.m1 + self.m2))*gamma_func(self.m1+1)*gamma_func(self.m2+1))
            self.mean = sum(self.moments)/self.sumfreq
            self.mode = self.mean - 0.5*self.mu3*(self.r + 2) / (self.mu2*(self.r-2))*self.step #TODO cumulative
            y = self.ye * ((1 + ((self.x - self.mean)/self.step)/self.A1)**(self.m1)) * ((1 - ((self.x - self.mean)/self.step)/self.A2)**(self.m2))
            # Find the index of the maximum value in y
            max_index = np.argmax(y)
            # Get the corresponding value from x
            self.mode = self.x[max_index]
            if showPrint:
                print(f'r = {self.r}')
                print(f'a1+a2 = {self.a1a2}')
                print(f'root1 = {self.root1}')
                print(f'root2 = {self.root2}')
                print(f'm1 = {self.m1}')
                print(f'm2 = {self.m2}')
                print(f'A1 = {self.A1}')
                print(f'A2 = {self.A2}')
                print(f'ye = {self.ye}')
                print(f'y0 = {self.y0}')
                print(f'mean = {self.mean}')
                print(f'mode = {self.mode}')
                # print(f'y = {y}')       
        elif self.curveType == 'IV':
            self.r = 6*(self.beta2-self.beta1-1)/(2*self.beta2-3*self.beta1-6)
            self.m = (self.r + 2)*0.5
            self.Z = (self.r**2) / (1-self.beta1*((self.r-2)**2)/(16*(self.r-1)))
            if self.mu3 < 0:
                self.v = abs((self.r*(self.r-2)*(self.beta1)**(1/2))/((16*(self.r-1)-self.beta1*(self.r-2)**2)**(1/2)))
            else: 
                self.v = -1 * abs((self.r*(self.r-2)*(self.beta1)**(1/2))/((16*(self.r-1)-self.beta1*(self.r-2)**2)**(1/2)))
            self.a = self.r*((self.mu2*(self.r-1)/(self.Z))**(1/2))
            self.aNew = self.a*self.step
            self.origin = self.v*self.a/self.r
            self.d = -1*self.mu3*(self.r-2)/(2*self.mu2*(self.r+2))
            self.pi = np.pi
            self.z = np.arange(0,1,(1/1000))
            self.teta = self.pi*self.z - self.pi/2 + 0.01
            self.X = self.a*np.tan(self.teta)-self.origin
            self.P = (np.cos(self.pi*self.z-self.pi/2)**(2*self.m))*np.exp(-self.v*(self.pi*self.z-self.pi/2))
            self.Pint = [0]
            i = 1
            while i < len(self.z):
                self.Pint.append(self.Pint[i-1]+0.5*(self.P[i]+self.P[i-1])*(self.X[i]-self.X[i-1]))
                i+=1
            self.y0 = sum(self.freq)/self.Pint[-1]
            # Origin fucks up! Of course it would fuck up if you put a '+' instead of a '*'. Degenerate!
            self.Origin = self.mean+self.v*self.aNew/self.r
            self.mode = self.mean-0.5*self.mu3*(self.r-2)/(self.mu2*(self.r+2)) * self.step
            y = self.y0*((1+(((self.x-self.Origin)/self.aNew)**2))**(-self.m))*np.exp(-self.v*np.arctan((self.x-self.Origin)/self.aNew))
            # Find the index of the maximum value in y
            max_index = np.argmax(y)
            # Get the corresponding value from x
            self.mode = self.x[max_index]
            if showPrint:
                print(f'r = {self.r}')
                print(f'm = {self.m}')
                print(f'Z = {self.Z}')
                print(f'v = {self.v}')
                print(f'a = {self.a}')
                print(f'aNew = {self.aNew}')
                print(f'origin = {self.origin}')
                print(f'd = {self.d}')
                # large arrays, tmp disabled output
                # print(f'z = {self.z}')
                # print(f'teta = {self.teta}')
                # print(f'X = {self.X}')
                # print(f'P = {self.P}')
                # print(f'Pint = {self.Pint}')
                print(f'y0 = {self.y0}')
                # quit()
                print(f'Origin = {self.Origin}')
                print(f'Mode = {self.mode}')
                # print(f'y = {y}')
        elif self.curveType == 'VI':
            self.r = 6*(self.beta2-self.beta1-1)/(-2*self.beta2+3*self.beta1+6)
            if self.mu3 >= 0 :
                self.a = (1/2)*(self.mu2*(self.beta1*(self.r+2)**2 + 16*(self.r+1)))**(1/2)
            else:
                self.a = -1*(1/2)*(self.mu2*(self.beta1*(self.r+2)**2 + 16*(self.r+1)))**(1/2)
            self.qroot1 = 0.5*(self.r - 2) + 0.5*self.r*(self.r+2)*(self.beta1/(self.beta1*((self.r+1)**2)+16*(self.r+1)))**(1/2)
            self.qroot2 = 0.5*(self.r - 2) - 0.5*self.r*(self.r+2)*(self.beta1/(self.beta1*((self.r+1)**2)+16*(self.r+1)))**(1/2)
            self.roots = [self.qroot1, self.qroot2]
            self.q1 = min(self.roots)
            if self.q1 < 0: self.q1 = self.q1*(-1)
            self.q2 = max(self.roots)
            # Thank you stackoverflow. The problem was something with negative floating points and numpy fucks up with them somehow
            # gamma_func dies with "OverflowError: math range error"
            self.y0 = sum(self.freq) * (np.sign(self.a) * (np.abs(self.a))**(self.q1-self.q2-1)) * gamma_func(self.q1)/(gamma_func(self.q1-self.q2-1)*gamma_func(self.q2+1))
            self.anew = self.step*self.a
            self.origin_mean = -1*self.a*(self.q1-1)/(self.q1-self.q2-2)
            self.mode_mean = -1*0.5*self.mu3*(self.r+2)/(self.mu2*(self.r-2))
            self.A1 = self.a*(self.q1-1)/((self.q1-1)-(self.q2+1))
            self.A2 = self.a*(self.q2+1)/(self.q1-1-(self.q2+1))

            tmparr = [self.q1, self.q2]
            if max(tmparr) < 50:
                self.ye = - sum(self.freq) * (np.sign(self.q2 + 1)*np.abs(self.q2 + 1)**self.q2)*(np.sign(self.q1 - self.q2 - 2)*np.abs(self.q1 - self.q2 - 2)**(self.q1 - self.q2)) * gamma_func(self.q1) / (self.a*((self.q1 - 1)**self.q1)*gamma_func(self.q1 - self.q2 - 1)*gamma_func(self.q2 + 1))
            else:
                print("Not implemented yet for max root >=50!")
                self.ye = - sum(self.freq) * (np.sign(self.q2 + 1)*np.abs(self.q2 + 1)**self.q2)*(np.sign(self.q1 - self.q2 - 2)*np.abs(self.q1 - self.q2 - 2)**(self.q1 - self.q2)) * gamma_func(self.q1) / (self.a*((self.q1 - 1)**self.q1)*gamma_func(self.q1 - self.q2 - 1)*gamma_func(self.q2 + 1))
            if self.a >= 0:
                y = - self.ye * ((1 + ((self.x - self.mean)/self.step)/self.A1)**(-self.q1)) * ((1 + ((self.x - self.mean)/self.step)/self.A2)**(self.q2))
                # y = - self.ye * (1 + ((self.x - self.mean)/(self.step*self.A1))**(-self.q1)) * (1 + ((self.x - self.mean)/(self.step*self.A2))**(self.q2))
            else:
                y = self.ye * ((1 + ((self.x - self.mean)/self.step)/self.A1)**(-self.q1)) * ((1 + ((self.x - self.mean)/self.step)/self.A2)**(self.q2))
                # y = self.ye * (1 + ((self.x - self.mean)/(self.step*self.A1))**(-self.q1)) * (1 + ((self.x - self.mean)/(self.step*self.A2))**(self.q2))
            # tmp and bad, do not like TODO (it will obviously stay)
            tmpcounter = 0
            while tmpcounter < len(y):
                if isnan(y[tmpcounter]):
                    y[tmpcounter] = 0
                tmpcounter+=1
            self.mode = self.mean+self.mode_mean*self.step
            # Find the index of the maximum value in y
            max_index = np.argmax(y)
            # Get the corresponding value from x
            self.mode = self.x[max_index]
            if showPrint:
                print(f'r = {self.r}')
                print(f'a = {self.a}')
                print(f'q1 = {self.q1}')
                print(f'q2 = {self.q2}')
                print(f'y0 = {self.y0}')
                print(f'anew = {self.anew}')
                print(f'origin_mean = {self.origin_mean}')
                print(f'mode_mean = {self.mode_mean}')
                print(f'A1 = {self.A1}')
                print(f'A2 = {self.A2}')
                print(f'mode = {self.mode}')
                print(f'ye = {self.ye}')
                # print(f'y = {y}')
                # print(f'x = {self.x}')
        else: quit() #TODO implement other types !!!
        return y

