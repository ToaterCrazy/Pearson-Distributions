import Stats as stats
import json
import matplotlib.pyplot as plt
import numpy as np

#Benchmarking via decorators
def benchmark(func):
    import time
    def wrapper(*args, **kwargs):
        start = time.time()
        #Execution of the function that is benchmarked
        return_value = func(*args, **kwargs)
        end = time.time()
        print(f'[*] Calculation time: {round(end-start, 4)} seconds.')
        return return_value
    return wrapper
#Some questionable color choices
BarColor = '#00000050'
BarEdgeColor = '#000000'
CurveColor = 'blue'
MeanColor = '#00ff00'
ModeColor = '#ff0000'

@benchmark
def RollStats(data: list, cutoff: int = 20, stepH: int = 20, showGraphs: bool = False, details: bool = False, showPrint: bool = False):
    """'data' - loaded JSON. 'cutoff' - minimum number of points in a dataset. 'stepH' - number of histogram bins. 
    'showGraphs' - show resulting graphs. 'details' - show details on graphs. 'showPrint' - show the calculation step-by-step."""
    counter = 0
    DepthArr = []
    DipArr = []
    AziArr = []
    CurveTypes = [0,0,0]
    SineKeyArray = []
    for sine_key, sine_values in data.items():
        array = data[sine_key]
        if len(array) > cutoff:
            SineKeyArray.append(sine_key)
            counter +=1
            print(sine_key)
            splitData = stats.retriveData(array)
            Dpts = splitData['Dpts']
            DepthsAuto = stats.CalcStats(Dpts, stepH, showPrint=showPrint)
            Dips = splitData['Dips']
            DipsAuto = stats.CalcStats(Dips, stepH, showPrint=showPrint)
            Azi = splitData['Azi']
            AziAuto = stats.CalcStats(Azi, stepH, showPrint=showPrint)
            DepthArr.append(DepthsAuto)
            DipArr.append(DipsAuto)
            AziArr.append(AziAuto)
            # Count types of curves
            # An extremely degenerate way to do this :)
            if DepthsAuto.curveType == 'I': CurveTypes[0] +=1
            elif DepthsAuto.curveType == 'IV': CurveTypes[1] +=1
            elif DepthsAuto.curveType == 'VI': CurveTypes[2] +=1
            if DipsAuto.curveType == 'I': CurveTypes[0] +=1
            elif DipsAuto.curveType == 'IV': CurveTypes[1] +=1
            elif DipsAuto.curveType == 'VI': CurveTypes[2] +=1
            if AziAuto.curveType == 'I': CurveTypes[0] +=1
            elif AziAuto.curveType == 'IV': CurveTypes[1] +=1
            elif AziAuto.curveType == 'VI': CurveTypes[2] +=1
            if showGraphs:
                fig, (ax1, ax2, ax3) = plt.subplots(1,3)
                fig.suptitle(f'{sine_key}')
                # Graph 1
                ax1.set_title(f'Curve Type {DepthsAuto.curveType}')
                ax1.plot(DepthsAuto.x, DepthsAuto.y, color=CurveColor, label='Auto')
                ax1.bar(DepthsAuto.bins2, DepthsAuto.freq, width=DepthsAuto.binWidth, align='center', color=BarColor, edgecolor=BarEdgeColor)
                ax1.set_ylim(0, max(DepthsAuto.freq)*1.2) #limit the graph at 120% of the highest hits value
                ax1.set_xlabel('Depth, m.')
                ax1.set_ylabel('Frequency')
                if details:
                    pixel_width = 20
                    pixel_per_data_unit = (ax1.get_xlim()[1] - ax1.get_xlim()[0]) / (fig.get_size_inches()[0] * fig.dpi)
                    bar_width = pixel_width * pixel_per_data_unit   
                    ax1.bar(DepthsAuto.mode, max(DepthsAuto.freq), width=bar_width, align='center', color=ModeColor, edgecolor=BarEdgeColor, label='Mode')
                    ax1.bar(DepthsAuto.mean, max(DepthsAuto.freq), width=bar_width, align='center', color=MeanColor, edgecolor=BarEdgeColor, label='Mean')
                    ax1.text(DepthsAuto.mode, max(DepthsAuto.freq)+max(DepthsAuto.freq)/10+0.2, str(round(DepthsAuto.mode, 2)))
                    ax1.text(DepthsAuto.mean, max(DepthsAuto.freq)+0.2, str(round(DepthsAuto.mean, 2)))
                # Graph 2
                ax2.set_title(f'Curve Type {DipsAuto.curveType}')
                ax2.plot(DipsAuto.x, DipsAuto.y, color=CurveColor, label='Auto')
                ax2.bar(DipsAuto.bins2, DipsAuto.freq, width=DipsAuto.binWidth, align='center', color=BarColor, edgecolor=BarEdgeColor)
                ax2.set_ylim(0, max(DipsAuto.freq)*1.2) #limit the graph at 120% of the highest hits value
                ax2.set_xlabel(f'Dip, deg.')
                ax2.set_ylabel('Frequency')
                if details:   
                    pixel_width = 20
                    pixel_per_data_unit = (ax2.get_xlim()[1] - ax2.get_xlim()[0]) / (fig.get_size_inches()[0] * fig.dpi)
                    bar_width = pixel_width * pixel_per_data_unit   

                    ax2.bar(DipsAuto.mode, max(DipsAuto.freq), width=bar_width, align='center', color=ModeColor, edgecolor=BarEdgeColor, label='Mode')
                    ax2.bar(DipsAuto.mean, max(DipsAuto.freq), width=bar_width, align='center', color=MeanColor, edgecolor=BarEdgeColor, label='Mean')
                    ax2.text(DipsAuto.mode, max(DipsAuto.freq)+max(DipsAuto.freq)/10+0.2, str(round(DipsAuto.mode, 2)))
                    ax2.text(DipsAuto.mean, max(DipsAuto.freq)+0.2, str(round(DipsAuto.mean, 2)))
                # Graph 3
                ax3.set_title(f'Curve Type {AziAuto.curveType}')
                ax3.plot(AziAuto.x, AziAuto.y, color=CurveColor, label='Auto')
                ax3.bar(AziAuto.bins2, AziAuto.freq, width=AziAuto.binWidth, align='center', color=BarColor, edgecolor=BarEdgeColor)
                ax3.set_ylim(0, max(AziAuto.freq)*1.2) #limit the graph at 120% of the highest hits value
                ax3.set_xlabel(f'Azimuth, deg.')
                ax3.set_ylabel('Frequency')
                if details:                     
                    pixel_width = 20
                    pixel_per_data_unit = (ax3.get_xlim()[1] - ax3.get_xlim()[0]) / (fig.get_size_inches()[0] * fig.dpi)
                    bar_width = pixel_width * pixel_per_data_unit   

                    ax3.bar(AziAuto.mode, max(AziAuto.freq), width=bar_width, align='center', color=ModeColor, edgecolor=BarEdgeColor, label='Mode')
                    ax3.bar(AziAuto.mean, max(AziAuto.freq), width=bar_width, align='center', color=MeanColor, edgecolor=BarEdgeColor, label='Mean')
                    ax3.text(AziAuto.mode, max(AziAuto.freq)+max(AziAuto.freq)/10+0.2, str(round(AziAuto.mode, 2)))
                    ax3.text(AziAuto.mean, max(AziAuto.freq)+0.2, str(round(AziAuto.mean, 2)))
                plt.show()
        else: pass
    ReturnList = [DepthArr, DipArr, AziArr] 
    ReturnArr = np.array(ReturnList).T # TODO rework the thingie in the top to np.arrays
    return counter, ReturnArr, CurveTypes, SineKeyArray

if __name__ == '__main__':
    #Load data
    with open(f"DATA.json", "r") as file:
        data = json.load(file)
    #Roll stats over data
    counter, ResultArray, CurveTypes, SineKeyArray = RollStats(data=data, cutoff = 100, stepH = 20, ShowGraphs=True, Details=False, showPrint=False)
    print(f'[*] Datasets calculated: {counter}.')
    #Print resulting mode values
    print('Mode values for each dataset:')
    for i in ResultArray:
        print(f'Depth = {round(i[0].mode, 3)}, Dip = {round(i[1].mode, 3)}, Azi = {round(i[2].mode, 3)}')
    print(f'Total number of distributions of:\nType I = {CurveTypes[0]}, Type IV = {CurveTypes[1]}, Type VI = {CurveTypes[2]}')