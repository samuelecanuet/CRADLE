from ROOT import *
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib as mpl
from matplotlib.lines import Line2D
from matplotlib.colors import Normalize
import re

mpl.rcParams["savefig.dpi"] = 300
# mpl.rcParams['savefig.bbox'] = 'tight'
# mpl.rcParams['savefig.pad_inches'] = 0.05

custom_params = {
        "xtick.direction" : "out",
        "ytick.direction" : "out",
        "lines.markeredgecolor" : "k",
        "lines.markeredgewidth" : 0.5,
        "lines.linewidth" : 1,
        "lines.markersize" : 5,
        "figure.figsize" : (16,9),
        "font.family" : "serif",
        "ytick.labelsize" : 16,
        "xtick.labelsize" : 16,
        "axes.labelsize" : 40,
        "axes.titlesize" : 30,
        "legend.fontsize" : 15,
        "text.usetex" : True,
        'figure.subplot.left': 0.20, 
        'figure.subplot.bottom': 0.15, 
        # 'figure.subplot.right': 0.95, 
        # 'figure.subplot.top': 0.90
        }
sns.set_theme(style = "ticks", rc=custom_params)

def DisplayTH1D(Hist, ax=None, color=None, fill_color=None, label=None, title=None, xlabel=None, ylabel=None, xlim=None, ylim=None, 
                xtick=None, ytick=None, ylog=None, xlog=None, rebin=None, normalized=None, scaler=None, visible=True, lw = None, 
                titlesize = None, labelsize=None, ticksize=None, draw=None, capsize=None, **kwargs):
        
        if ax == None: fig, ax = plt.subplots()
        if rebin : Hist.Rebin(rebin)
        if normalized : integral = Hist.Integral("width")
        else : integral = 1.
        if scaler : integral = scaler
        nbins_x = Hist.GetNbinsX()
        hist_data = np.zeros(nbins_x)
        bin_centers_x = np.zeros(nbins_x)

        if xlim != None: 
            binstart = Hist.GetXaxis().FindBin(xlim[0])
            binend = Hist.GetXaxis().FindBin(xlim[1])
        else: 
            binstart = 1
            binend = nbins_x
        for i in range(binstart, binend + 1):
            hist_data[i - 1] = Hist.GetBinContent(i)/integral
            bin_centers_x[i - 1] = Hist.GetXaxis().GetBinCenter(i)
        
        hist_color = gROOT.GetColor(Hist.GetLineColor())
        hist_fillcolor = gROOT.GetColor(Hist.GetFillColor())
        if color    == None: color = (hist_color.GetRed(), hist_color.GetGreen(), hist_color.GetBlue())
        if color    == "auto": color = None
        if fill_color != None : ax.fill_between(bin_centers_x, 0, hist_data, step='pre', color=fill_color, alpha = kwargs["alpha"], linewidth = 0, label = label)
        if label    == None: label = Hist.GetTitle()    
        if visible  == False: label = '_nolegend_'          
        if title    == None: title = Hist.GetTitle()
        if xlabel   == None: 
            xlabel = Hist.GetXaxis().GetTitle().replace("#theta", r"\theta").replace("#beta", r"\beta").replace("#nu", r"\nu").replace("#phi", r"\phi") 
            xlabel = "$" + xlabel + "$"   
        if ylabel   == None: ylabel = Hist.GetYaxis().GetTitle()
        if xlim     == None: xlim = ( bin_centers_x.min(), bin_centers_x.max() )
        if ylim     == None and hist_data.max()*1.1 > ax.get_ylim()[1] : ylim = ( 0, hist_data.max()*1.1 )
        if xtick    != None: ax.set_xticks(np.linspace(xlim[0], xlim[1], xtick))
        if ytick    != None: ax.set_yticks(np.linspace(ylim[0], ylim[1], ytick))
        if ticksize != None: ax.tick_params(labelsize=ticksize)
        if xlog     != None: ax.set_xscale('log')
        if ylog     != None: 
            ax.set_yscale('log')
            # ylim = ( 1, hist_data.max()*1.1 )
 
        
        if visible != False:
            ax.set_xlabel(xlabel, fontsize = labelsize)
            ax.set_ylabel(ylabel, fontsize = labelsize)
            ax.set_title(title.replace("_", "\_"), fontsize = titlesize)
            ax.set_xlim(xlim)
            ax.set_ylim(ylim)

        if draw == "E":
            # xerr = np.array([Hist.GetBinWidth(i)/2 for i in range(1, nbins_x + 1)])
            yerr = np.array([Hist.GetBinError(i)/scaler for i in range(1, nbins_x + 1)])
            step = ax.errorbar(bin_centers_x, hist_data, yerr=yerr, label = label, fmt='o', color=color, visible = visible, capsize = capsize, **kwargs)
        else:
            step, = ax.step(bin_centers_x, hist_data, label=label, color=color, visible = visible, linewidth=lw)

        return step
    
def DisplayTH2D(Hist, ax=None, color='plasma', label=None, title=None, xlabel=None, ylabel=None, xlim=None, ylim=None, xtick=None, ytick=None, ylog=None, xlog=None, zlog=None, rebinx=None, rebiny=None, vmax = None, vmin = None, view=None, scale = None, **kwargs):
        if rebinx   != None: Hist.RebinX(rebinx)
        if rebiny   != None: Hist.RebinY(rebiny)

    
        nbins_x = Hist.GetNbinsX()
        nbins_y = Hist.GetNbinsY()

        hist_data = np.zeros((nbins_x, nbins_y))
        bin_centers_x = np.zeros(nbins_x)
        bin_centers_y = np.zeros(nbins_y)

        BD, BG, HG, HD = [0, 0], [0, 0], [0, 0], [0, 0]
        for i in range(1, nbins_x + 1):
            for j in range(1, nbins_y + 1):
                hist_data[i - 1, j - 1] = Hist.GetBinContent(i, j)
                bin_centers_x[i - 1] = Hist.GetXaxis().GetBinCenter(i)
                bin_centers_y[j - 1] = Hist.GetYaxis().GetBinCenter(j)

                if BD[0] < Hist.GetBinContent(i, j) and np.sqrt(Hist.GetXaxis().GetBinCenter(i)**2+Hist.GetYaxis().GetBinCenter(j)**2) > BD[1] and Hist.GetXaxis().GetBinCenter(i) > 0. and Hist.GetYaxis().GetBinCenter(j) < 0.:
                    BD = [Hist.GetBinContent(i, j), np.sqrt(Hist.GetXaxis().GetBinCenter(i)**2+Hist.GetYaxis().GetBinCenter(j)**2)]
                
                if BG[0] < Hist.GetBinContent(i, j) and np.sqrt(Hist.GetXaxis().GetBinCenter(i)**2+Hist.GetYaxis().GetBinCenter(j)**2) > BG[1] and Hist.GetXaxis().GetBinCenter(i) < 0. and Hist.GetYaxis().GetBinCenter(j) < 0.:
                    BG = [Hist.GetBinContent(i, j), np.sqrt(Hist.GetXaxis().GetBinCenter(i)**2+Hist.GetYaxis().GetBinCenter(j)**2)]
                
                if HG[0] < Hist.GetBinContent(i, j) and np.sqrt(Hist.GetXaxis().GetBinCenter(i)**2+Hist.GetYaxis().GetBinCenter(j)**2) > HG[1] and Hist.GetXaxis().GetBinCenter(i) < 0. and Hist.GetYaxis().GetBinCenter(j) > 0.:
                    HG = [Hist.GetBinContent(i, j), np.sqrt(Hist.GetXaxis().GetBinCenter(i)**2+Hist.GetYaxis().GetBinCenter(j)**2)]
                
                if HD[0] < Hist.GetBinContent(i, j) and np.sqrt(Hist.GetXaxis().GetBinCenter(i)**2+Hist.GetYaxis().GetBinCenter(j)**2) > HD[1] and Hist.GetXaxis().GetBinCenter(i) > 0. and Hist.GetYaxis().GetBinCenter(j) > 0.:
                    HD = [Hist.GetBinContent(i, j), np.sqrt(Hist.GetXaxis().GetBinCenter(i)**2+Hist.GetYaxis().GetBinCenter(j)**2)]

        lim = max([i[1] for i in [BD, BG, HG, HD]])/np.sqrt(2)*1.1


        if label  == None: label = Hist.GetTitle()
        if title  == None: title = Hist.GetTitle()
        if xlabel == None: xlabel = Hist.GetXaxis().GetTitle()
        if ylabel == None: ylabel = Hist.GetYaxis().GetTitle()
        if xlim   == None: xlim = ( bin_centers_x.min(), bin_centers_x.max() )
        if ylim   == None: ylim = (bin_centers_y.min(), bin_centers_y.max())
        if xlim   == 'auto': xlim = (-lim, lim)
        if ylim   == 'auto': ylim = (-lim, lim)
        if xtick  != None: ax.set_xticks(np.linspace(xlim[0], xlim[1], xtick))
        if ytick  != None: ax.set_yticks(np.linspace(ylim[0], ylim[1], ytick))
        if xlog   != None: ax.set_xscale('log')
        if ylog   != None: ax.set_yscale('log')
        if scale  != None: hist_data = hist_data*scale

        cax = ax.imshow(hist_data.T, extent=(bin_centers_x.min(), bin_centers_x.max(), bin_centers_y.min(), bin_centers_y.max()), origin='lower', aspect='equal', cmap=color, norm=Normalize(vmin=0, vmax=hist_data.max()))    

        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title(title.replace("_", "\_"))
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)

        if vmax      != None: cax.set_clim(vmax = vmax)
        if vmin      != None: cax.set_clim(vmin = vmin)

        return [cax, ax]

def DisplayCanvas(canvas, **kwargs):
    if "ax" not in list(kwargs.keys()): 
        fig, kwargs["ax"]  = plt.subplots(figsize = (canvas.GetXsizeReal()/2, canvas.GetYsizeReal ()/2))

    primitives = canvas.GetListOfPrimitives()

    for i, primitive in enumerate(primitives):
        if isinstance(primitive, TPad) or isinstance(primitive, TFrame):
            ax = plt.axes([primitive.GetAbsXlowNDC(), primitive.GetAbsYlowNDC(), primitive.GetAbsWNDC (), primitive.GetAbsHNDC ()])
            for content in primitive.GetListOfPrimitives():
                Display(content, ax=ax)
        else:
            fig = Display(primitive, **kwargs)
    
    if canvas.GetTitle() != "Canvas": 
        if ("title" in list(kwargs.keys())): kwargs["ax"].set_title(kwargs["title"])
        else: kwargs["ax"].set_title(canvas.GetTitle().replace("_", "\_"))
    # return fig

def GetColor(object):
    color = gROOT.GetColor(object.GetLineColor())
    return (color.GetRed(), color.GetGreen(), color.GetBlue())

def DisplayTF1(object, ax=None, color=None, label=None, title=None, nolimits=None, **kwargs):
    if label    == None : label = object.GetTitle()
    if color    == None : color = GetColor(object)
    # if title    == None : title = object.GetTitle()
    if nolimits == True : x_values = np.linspace(-1e6, 1e6, 2000)
    else : x_values = np.linspace(object.GetXmin(), object.GetXmax(), int(object.GetXmax()-object.GetXmin())*100) 
    ax.plot(x_values, [object.Eval(x) for x in x_values], color=color, label=label, **kwargs)

    return ax

def DisplayGraphError(graph, ax=None, color=None, label=None, title=None, xlabel=None, ylabel=None, xlim=None, ylim=None, ylog=None, **kwargs):
    if not ax:
        fig, ax = plt.subplots()
    
    n_points = graph.GetN()
    x_data = np.zeros(n_points)
    y_data = np.zeros(n_points)
    x_errors = np.zeros(n_points)
    y_errors = np.zeros(n_points)
    
    for i in range(n_points):
        x_data[i] = graph.GetPointX(i)
        y_data[i] = graph.GetPointY(i)
        x_errors[i] = graph.GetErrorX(i)
        y_errors[i] = graph.GetErrorY(i)
    
    if color is None:
        graph_color = gROOT.GetColor(graph.GetMarkerColor())
        color = (graph_color.GetRed(), graph_color.GetGreen(), graph_color.GetBlue())

   
    obj = ax.errorbar(x_data, y_data, xerr=x_errors, yerr=y_errors, color=color, label=label, markeredgecolor = color, ecolor=color, **kwargs)
    
    if xlabel is None:
        xlabel = graph.GetXaxis().GetTitle().replace("#beta", r"$\beta$")
    if ylabel is None:
        ylabel = graph.GetYaxis().GetTitle().replace("#beta", r"$\beta$")
    if title is None:
        title = graph.GetTitle().replace("#beta", r"$\beta$")
    # if xlim is None:
    #     xlim = (x_data.min() - x_errors.max(), x_data.max() + x_errors.max())
    # if ylim is None:
    #     ylim = (y_data.min() - y_errors.max(), y_data.max() + y_errors.max())
    
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title.replace("_", "\_"))
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    if ylog != None:
        ax.set_yscale('log')
    
    if label is not None:
        ax.legend()
    
    return obj


def DisplayTLegend(object, ax=None, **kwargs):
    labels = []
    handles = []
    for entry in object.GetListOfPrimitives():
        labels.append(entry.GetLabel().replace("#beta", r"$\beta$"))
        handles.append(Line2D([0], [0], color=GetColor(entry.GetObject()), linewidth=2))
    ax.legend(handles=handles, labels=labels, fontsize = "medium")

    return ax

def DisplayTLatex(object, ax=None, **kwargs):
    try: 
        x = kwargs.pop("x") * (ax.get_xlim()[1] - ax.get_xlim()[0]) + ax.get_xlim()[0]
    except KeyError:
        x = object.GetX()
    try:
        y = kwargs.pop("y") * (ax.get_ylim()[1] - ax.get_ylim()[0]) + ax.get_ylim()[0]
    except KeyError:
        y = object.GetY()
    text = ConvertRootLatex(object.GetTitle())
    ax.text(x, y, text, **kwargs)

    return ax

def Display(object, **kwargs):
    figure = None
    # print(object.IsA().GetName())
    if isinstance(object, TCanvas):
        figure = DisplayCanvas(object, **kwargs)
    elif isinstance(object, TH2):
        figure = DisplayTH2D(object, **kwargs)
    elif isinstance(object, TH1):
        figure = DisplayTH1D(object, **kwargs)
    elif isinstance(object, TF1):
        figure = DisplayTF1(object, **kwargs)
    elif isinstance(object, TGraphErrors):
        figure = DisplayGraphError(object, **kwargs)
    elif isinstance(object, TLegend):
        if kwargs["legend"] != None:
            figure = DisplayTLegend(object, **kwargs)
            
    else:
        print(f"Sorry, {object.IsA().GetName()} is not yet supported by root2mpl")
        return None
        

    return figure

def ConvertRootLatex(root_string):
    converted = re.sub(r"#\^\{([^}]+)\}", r"^{\1}", root_string)
    
    # Handle subscripts, e.g., #_{text} becomes _{text}
    converted = re.sub(r"#_\{([^}]+)\}", r"_{\1}", converted)
    
    # Replace ROOT-style Greek symbols, e.g., #beta becomes \beta
    converted = re.sub(r"#([a-zA-Z]+)", r"\\\1", converted)
    
    # Wrap the final string in $ symbols for LaTeX rendering in Matplotlib
    return f"${converted}$".replace("#", "")