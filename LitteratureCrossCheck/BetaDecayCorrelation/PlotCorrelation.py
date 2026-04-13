from pathlib import Path
import sys
import matplotlib.pyplot as plt
from ROOT import TFile
import numpy as np
from matplotlib.colors import LinearSegmentedColormap

# Make local utilities importable regardless of current working directory.
PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT / "include"))

import root2mpl

## main function
if __name__ == "__main__":
    ## get args
    if len(sys.argv) != 2:
        print("Usage: python3 PlotCorrelation.py <input_root_file>")
        sys.exit(1)


    readfilename = sys.argv[1]
    # Open the ROOT file
    file = TFile(readfilename, "READ")
    if not file.IsOpen():
        print(f"Error: Could not open file {readfilename}")
        sys.exit(1)
    

    # Get the canvas 
    for key in file.GetListOfKeys():
        if key.GetClassName() == "TCanvas":
            canvas = key.ReadObj()
            break

    # 2x2 grid of subplots
    tpadcounter = 0
    fig, ax = plt.subplots(2, 2, figsize=(16, 10), constrained_layout=True)
        
    for primitive in canvas.GetListOfPrimitives():
        if primitive.ClassName() == "TPad":
            axi = ax[tpadcounter//2][tpadcounter%2]
            for element in primitive.GetListOfPrimitives():
                if element.ClassName() == "TH1D":
                    root2mpl.DisplayTH1D(element, ax=axi, color="black", draw="E", scaler=2, capsize = 2, rebin = 2)
                    
                if (element.ClassName() == "TF1"):
                    root2mpl.DisplayTF1(element, ax=axi, color="red", zorder = 10, linewidth = 2)

            ### 5 ticks over x axis 
            axi.set_xticks(np.linspace(axi.get_xlim()[0], axi.get_xlim()[1], 5))
            tpadcounter += 1    
                

    plt.show()



