import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import ticker
import scipy.stats as linstats
import numpy as np
from ..vis_config import AP_COLOR, P_COLOR, generate_colors

font = {
        'size'   : 16
        }
matplotlib.rc('font', **font)

def equal_fn(dat):
    x = [1, 2, 3, 4]
    y = [1, 2, 3, 4]
    equality = linstats.linregress(x,y)
    return np.poly1d([equality[0], equality[1]])(dat)

# Not used but retained for legacy purposes
def sel_unsel_plot(unsel, sel, unsel_ap, sel_ap, unsel_null, sel_null, stats_corrected,
    xlab='Unselected', ylab='Selected'):
    plt.clf()
    fig, ax = plt.subplots()
    stat_colors = generate_colors(stats_corrected, ylab)
    edge_colors = stat_colors

    ap_fit = linstats.linregress(unsel_null, sel_null)
    ap_fit_fn = np.poly1d([ap_fit[0], ap_fit[1]])

    ax.scatter(unsel, sel, color=stat_colors, edgecolors=edge_colors, 
            label = 'Parallel: %d'%len(sel))
    ax.scatter(unsel_ap, sel_ap, facecolors='None', edgecolors='k', 
    #ax.scatter(unsel_null, sel_null, facecolors='None', edgecolors='k', 
            label = 'AntiParallel: %d'%len(sel_ap))

    max_s_P = max(sel)
    max_us_P = max(unsel)

    xmax = int(max([max_s_P, max_us_P])*2)
    plot_range = range(-1, xmax)
    ax.plot(plot_range, equal_fn(plot_range), '-', color = '#000000', ms=0)
    ax.plot(plot_range, ap_fit_fn(plot_range), '--', color = AP_COLOR, ms=0)

    ax.set_xscale('symlog')
    ax.set_yscale('symlog')
    ax.set_xlim([-0.5, xmax])
    ax.set_ylim([-0.5, xmax])
    ax.set_xlabel('%s Counts'%xlab)
    ax.set_ylabel('%s Counts'%ylab)
    ax.tick_params(direction='out', top=False, right=False)
    ax.legend(loc='upper left', fontsize=7)

    plt.tight_layout()
    plt.savefig('sel_unsel_plot.png')


def position_plot(fold_frame_parallel, stats_corrected, observed, unobserved, ax, ylab=''):
    gene_length = len(observed) + len(unobserved) + 1

    if ylab=='Assisted':
        max_idx = np.argmin(stats_corrected)
        max_pos = observed[max_idx]
        if stats_corrected[max_idx]<-3:
            ax.arrow(max_pos, fold_frame_parallel[max_idx]+1.75, 0, -1.5, width=(0.001*gene_length/100),\
                 fc='k', overhang=1.0, head_width=(0.5*gene_length/100), head_length=0.3, length_includes_head=True)
    
    stat_colors = generate_colors(stats_corrected, ylab)

    #Visualize examples where counts are equal
    jitter = np.where(fold_frame_parallel == 0)[0]
    fold_frame_parallel[jitter] += 0.2
    fold_frame_parallel = np.clip(fold_frame_parallel, -6, 6)

    #plot fold-change profile
    ax.bar(observed, fold_frame_parallel, 
            color=stat_colors, width=1, linewidth=0, snap=False)

    #Adjust plot configuration
    ax.axhline(y=0, xmin=0, xmax=gene_length+1, linewidth=0.25, color = [0, 0, 0, 1])
    ax.set_xlim([0, gene_length])
    ax.set_ylim([-6., 6.])
    ax.tick_params(direction='out', top=False, right=False)
    ax.yaxis.set_major_locator(ticker.MultipleLocator(2)) #tickspacing
    ax.xaxis.set_major_locator(ticker.MultipleLocator(20))
    ax.xaxis.set_minor_locator(ticker.MultipleLocator(1))
    ax.yaxis.set_ticklabels([])
    ax.set_ylabel(ylab)
    
# Not used but retained for legacy purposes
def fold_abundance_plot(fold_frame_antiparallel, fold_frame_parallel, fit_abund, fit_mean,
        fit_sd, stat_colors, edge_colors, unsel, unsel_ap, xlab='Unselected'):
    plt.clf()
    fig, ax = plt.subplots()

    #print(fold_frame_parallel)
    ymin = 1.1*min((np.nanmin(fold_frame_antiparallel),np.nanmin(fold_frame_parallel)))
    #fold_frame_parallel = [ymin if numpy.isnan(x) else x for x in fold_frame_parallel]
    #fold_frame_antiparallel = [ymin if numpy.isnan(x) else x for x in fold_frame_antiparallel]
    fold_frame_parallel = [(0.9*(ymin*1.1)) if np.isnan(x) else x for x in fold_frame_parallel]
    fold_frame_antiparallel = [(0.9*(ymin*1.1)) if np.isnan(x) else x for x in fold_frame_antiparallel]
    #print len(fold_frame_parallel), len(fold_frame_antiparallel), len(unsel), len(unsel_ap), len(stat_colors), len(edge_colors)

    #background = [x for x in fold_frame_antiparallel if x > -9]
    #background_mean = numpy.average(background)
    #background_mean = fit_mean[-1]
    #print('mean: ' + str(background_mean))
    #print len(unsel_ap)
    #print zip(unsel_ap, fold_frame_antiparallel)

    ax.scatter(unsel, fold_frame_parallel, color=stat_colors,
            edgecolors=edge_colors, label='Parallel')
    ax.scatter(unsel_ap, fold_frame_antiparallel, facecolors='none',
            edgecolors='k', label='Antiparallel')

    ax.plot(fit_abund, fit_mean, '-', color=(0,0,0,1), label="Mean")
    mean_plusSD = [fit_mean[j] + 2*fit_sd[j] for j in range(0, len(fit_mean))]
    mean_minusSD = [fit_mean[j] - 2*fit_sd[j] for j in range(0, len(fit_mean))]

    ax.plot(fit_abund, mean_plusSD, '--', color=(0,0,0,1), label="Mean + 2$\sigma$")
    #ax.plot(fit_abund, mean_minusSD, '--', color=(0,0,0,1), label="Mean - 2$\sigma$")

    xmin, xmax = ax.get_xlim()
    #ax.axhline(y=background_mean, xmin=xmin, xmax=xmax,
    #        linewidth=1, color='k', label='mean')
    ax.axhline(y=ymin, xmin=xmin, xmax=xmax,
            linewidth=8, color=[0, 0, 0, 0.4], label='No Longer Observed')
    ax.legend(fontsize=7)
    #ax.legend(loc='center left', bbox_to_anchor=(1.05, 0.5) , fontsize=7)
    _, ymax = ax.get_ylim()
    ax.set_ylim([ymin*1.05, ymax])
    #ax.set_ylim([ymin*1.05, 6])
    ax.tick_params(direction='out', top=False, right=False)
    ax.yaxis.set_major_locator(ticker.MultipleLocator(2)) #tickspacing
    ax.set_xlim([-0.5, max([max(unsel)]) * 1.5])
    ax.set_xscale('symlog')
    ax.set_ylabel("Apparent Fold Change (log2)")
    ax.set_xlabel('%s Counts'%xlab)
    plt.tight_layout()
    plt.savefig('fchange.png')

