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

# Load the ROOT file
file = TFile("Cascades.root")

counter = 0
# Load Canvas
for key in file.GetListOfKeys():
    if (key.GetName() != "Analytical"):   
        fig, ax = plt.subplots(figsize = (10, 6))
        obj = key.ReadObj()

        for i in range(obj.GetListOfPrimitives().GetEntries()):
            primitive = obj.GetListOfPrimitives().At(i)
            if primitive.InheritsFrom("TH1D"):
                root2mpl.DisplayTH1D(primitive, ax=ax, draw="E", color="black", capsize = 2, rebin = 2, scaler = 2, label = "CRADLE++")

            if primitive.InheritsFrom("TF1"):
                root2mpl.DisplayTF1(primitive, ax=ax, color="red", label = r"$W(\theta)$", lw = 2, zorder = 10)           

        ax.set_xlim(0, 180)
        ax.set_xlabel(r"$\theta (deg)$", fontsize = 20)
        ax.set_ylabel(r"$W(\theta)$", fontsize = 20)
        
        ax.set_ylim(primitive.GetMinimum()*0.8, primitive.GetMaximum()*1.2)
        title = "$" + key.GetName().replace("->", r" \rightarrow ") + "$"
        ax.set_title(r"{}".format(title))
        plt.legend(loc="upper right")
        plt.savefig(f"{key.GetName()}.pdf")

        counter+=1
    
    else:
        cmap = plt.get_cmap("mako")
        colors = [cmap(i) for i in np.linspace(0.1, 0.8, counter)]
        fig, ax = plt.subplots(figsize = (10, 6))
        obj = key.ReadObj()
        for i in range(obj.GetListOfPrimitives().GetEntries()):
            primitive = obj.GetListOfPrimitives().At(i)
            if primitive.InheritsFrom("TF1"):
                label = "$" + primitive.GetName().replace("->", r" \rightarrow ") + "$"
                root2mpl.DisplayTF1(primitive, ax=ax, lw = 2, title = r"{}".format(label), label = r"{}".format(label), color=colors[i])

        ax.set_xlim(0, 180)
        ax.set_ylim(0.4, 2.4)
        ax.set_xlabel(r"$\theta (deg)$", fontsize = 20)
        ax.set_ylabel(r"$W(\theta)$", fontsize = 20)
        ax.set_title("Analytical", fontsize = 20)
        plt.legend(loc="center right", bbox_to_anchor=(1.3, 0.5))
        plt.tight_layout()
        plt.savefig(f"{key.GetName()}.pdf")