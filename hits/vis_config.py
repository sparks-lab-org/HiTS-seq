import numpy as np
from matplotlib import cm
from matplotlib.colors import ListedColormap
import matplotlib.pyplot as plt
from matplotlib import ticker

AP_COLOR = [39./255, 174./255, 96./255, 1.0]
P_COLOR = [142./255, 68./255, 173./255, 1.0]
SIG_THRESH = -3
SIG_RANGE = np.linspace(-52, -2, 50)
#RGB = cm.get_cmap(plt.get_cmap('Greys'))(numpy.linspace(0.5, 1.0, 50))[numpy.newaxis,:,:3][0][::-1]
RGB0 = cm.get_cmap(plt.get_cmap('Reds'))(np.linspace(0.5, 1.0, 50))[np.newaxis,:,:3][0][::-1]
RGB1 = cm.get_cmap(plt.get_cmap('Blues'))(np.linspace(0.5, 1.0, 50))[np.newaxis,:,:3][0][::-1]
RGB2 = cm.get_cmap(plt.get_cmap('Greens'))(np.linspace(0.5, 1.0, 50))[np.newaxis,:,:3][0][::-1]
CMAP = ListedColormap(RGB0)
NONSIG_COLOR = [0.5, 0.5, 0.5, 0.5]
#NONSIG_COLOR = [0.5,0.5,0.5]

#RGB_dict = {'Self (AP)': RGB, 'Truncated AP1': RGB, 'Self (P)': RGB1, 'Assisted': RGB2 }
RGB_dict = {'Self': RGB1, 'Assisted': RGB0}

def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return array[idx]

def generate_colors(stats_corrected, RGB_label):
    RGB = RGB_dict.get(RGB_label, RGB0)
    stat_colors = []
    #line_widths = []
    #edge_colors = []
    #sig_count = 0
    for stat_sample in stats_corrected:
        if stat_sample > SIG_THRESH:
            stat_colors.append(NONSIG_COLOR)
        if stat_sample <= SIG_THRESH:
            stat_colors.append(RGB[int(np.where(SIG_RANGE == find_nearest(SIG_RANGE, stat_sample))[0])])
    #print(stat_colors)
    return np.array(stat_colors, dtype='object')

def my_cbar(cbar_ax):
    #legend,swap range
    sm = plt.cm.ScalarMappable(cmap=CMAP,
            norm=plt.Normalize(vmin=int(min(SIG_RANGE)), vmax=int(max(SIG_RANGE))))
    # fake up the array of the scalar mappable. Urgh...
    sm._A = []
    #formatter = ticker.LogFormatter(10, labelOnlyBase=False)
    cb = plt.colorbar(sm, cax=cbar_ax)#, format = formatter)
    font_size = 20 # Adjust as appropriate.
    cb.ax.tick_params(labelsize=font_size)
    tick_locator = ticker.MaxNLocator(nbins=7)
    cb.locator = tick_locator
    cb.update_ticks()