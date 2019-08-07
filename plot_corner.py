from imports import *
from priors import *
from astropy.io import fits


def get_result_array():
    # get MCMC results
    samples = fits.open('output_data/3planets_TOI175_H164_HalphaGP_samplesv1')[0].data
    #samples = fits.open('output_data/3planets_TOI175_H113_HalphaGP_samplesv2_reducedmemory')[0].data
    
    results = get_results(samples)
    # fix Kb
    v = np.percentile(samples[:,18], (16,50,84))
    results[:,18] = v[1], v[2]-v[1], v[1]-v[0]
    hdu = fits.PrimaryHDU(results)
    hdu.writeto('output_data/3planets_TOI175_H164_HalphaGP_resultsv1', overwrite=True)

    #results = fits.open('output_data/3planets_TOI175_H113_HalphaGP_resultsv2')[0].data
    
    return samples, results


def get_xticks(j):
    if j == 0: # lna
        ticks = np.arange(0,6,1)
        ticklabels = np.array(['%i'%i for i in ticks])
        ticklabels[np.arange(1,6,2)] = ''
    elif j == 1: # lnl
        ticks = np.arange(0,12,2)
        ticklabels = np.array(['%i'%i for i in ticks])
        ticklabels[np.arange(1,6,2)] = ''
    elif j == 2: # lnG
        ticks = np.arange(-4,2,2)
        ticklabels = np.array(['%i'%i for i in ticks])
        #ticklabels[np.arange(1,,2)] = ''
    elif j == 3: # lnP
        ticks = np.arange(2,4,.5)
        ticklabels = np.array(['%i'%i for i in ticks])
        ticklabels[np.arange(1,4,2)] = ''
    elif j == 4: # s
        ticks = np.arange(0,1.2,.2)
        ticklabels = np.array(['%.1f'%i for i in ticks])
        ticklabels[np.arange(1,6,2)] = ''
    elif j == 5: # V
        ticks = np.arange(-2,3,2)
        ticklabels = np.array(['%i'%i for i in ticks])
        #ticklabels[np.arange(1,3,2)] = ''
    elif j == 6: # Pb
        ticks = np.arange(2.252,2.2545,1e-3)
        ticklabels = np.array(['%.3f'%i for i in ticks])
        ticklabels[np.arange(0,4,2)] = ''
    elif j == 7: # T0-T0true
        ticks = np.arange(-4,5,2)*1e-4
        ticklabels = np.array(['%.4f'%i for i in ticks])
        ticklabels[np.arange(0,5,2)] = ''
    elif j == 8: # Kb
        ticks = np.arange(0,2.5,.5)
        ticklabels = np.array(['%i'%i for i in ticks])
        ticklabels[np.arange(1,5,2)] = ''
    elif j == 9: # hb
        ticks = np.arange(-1,1.5,.5)
        ticklabels = np.array(['%i'%i for i in ticks])
        ticklabels[np.arange(1,5,2)] = ''
    elif j == 10: # kb
        ticks = np.arange(-1,1.5,.5)
        ticklabels = np.array(['%i'%i for i in ticks])
        ticklabels[np.arange(1,5,2)] = ''
    elif j == 11: # Pc
        ticks = np.arange(3.6895,3.69151,1e-3)
        ticklabels = np.array(['%.4f'%i for i in ticks])
        ticklabels[np.arange(0,4,2)] = ''
    elif j == 12: # T0-T0true
        ticks = np.arange(-2,2.1,1)*1e-3
        ticklabels = np.array(['%.3f'%i for i in ticks])
        ticklabels[np.arange(0,5,2)] = ''
    elif j == 13: # Kc
        ticks = np.arange(1,3.6,1)
        ticklabels = np.array(['%i'%i for i in ticks])
        #ticklabels[np.arange(1,3,2)] = ''
    elif j == 14: # hc
        ticks = np.arange(-.5,.6,.25)
        ticklabels = np.array(['%.1f'%i for i in ticks])
        ticklabels[np.arange(1,5,2)] = ''
    elif j == 15: # kc
        ticks = np.arange(-.75,.8,.25)
        ticklabels = np.array(['%.1f'%i for i in ticks])
        ticklabels[np.arange(0,7,2)] = ''
    elif j == 16: # Pd
        ticks = np.arange(7.448,7.455,2e-3)
        ticklabels = np.array(['%.3f'%i for i in ticks])
        ticklabels[np.arange(0,4,2)] = ''
    elif j == 17: # T0-T0true
        ticks = np.arange(-2,3,1)*1e-3
        ticklabels = np.array(['%.3f'%i for i in ticks])
        ticklabels[np.arange(0,5,2)] = ''
    elif j == 18: # Kd
        ticks = np.arange(0,4,1)
        ticklabels = np.array(['%i'%i for i in ticks])
        ticklabels[np.arange(1,5,2)] = ''
    elif j == 19: # hd
        ticks = np.arange(-.6,.61,.3)
        ticklabels = np.array(['%.1f'%i for i in ticks])
        ticklabels[np.arange(1,5,2)] = ''
    elif j == 20: # kd
        ticks = np.arange(-1,1.1,.5)
        ticklabels = np.array(['%i'%i for i in ticks])
        ticklabels[np.arange(1,5,2)] = ''
        
    return ticks, ticklabels



def get_xlim(j):
    if j == 0:
        xlim = (0,5)
    elif j == 1:
        xlim = (0,12)
    elif j == 2:
        xlim = (-4,2)
    elif j == 3:
        xlim = (2,4)
    elif j == 4:
        xlim = (0,1)
    elif j == 5:
        xlim = (-2,2)
    elif j == 6:
        xlim = (2.252,2.2545)
    elif j == 7:
        xlim = (-4e-4,4e-4)
    elif j == 8:
        xlim = (0,2)
    elif j == 9:
        xlim = (-1,1)
    elif j == 10:
        xlim = (-1,1)
    elif j == 11:
        xlim = (3.6895,3.6915)
    elif j == 12:
        xlim = (-2e-3,3e-3)
    elif j == 13:
        xlim = (1,3.5)
    elif j == 14:
        xlim = (-.5,.5)
    elif j == 15:
        xlim = (-.75,.75)
    elif j == 16:
        xlim = (7.448,7.454)
    elif j == 17:
        xlim = (-3e-3,3e-3)
    elif j == 18:
        xlim = (0,3.5)
    elif j == 19:
        xlim = (-.6,.6)
    elif j == 20:
        xlim = (-1,1)
    return xlim


        
def plot_corner(samples, results, pltt=False, label=True):
    bins = 30
    order = np.array([0,1,2,3,4,5,16,17,18,19,20,6,7,8,9,10,11,12,13,14,15])
    samples2 = samples[:,order]
    results2 = results[:,order][0]
    labels = ['$\ln{a}$','$\ln{\lambda}$','$\ln{\Gamma}$','$\ln{P_{GP}}$',
              '$s$\n[m/s]','$\gamma_0$\n[m/s-5678.16]',
	      '$P_b$\n[days]','$T_{0,b}$ [BJD\n-2458366.1708]','$K_b$\n[m/s]','$h_b$','$k_b$',
	      '$P_c$\n[days]','$T_{0,c}$ [BJD\n-2458367.2752]','$K_c$\n[m/s]','$h_c$','$k_c$',
	      '$P_d$\n[days]','$T_{0,d}$ [BJD\n-2458362.7376]','$K_d$\n[m/s]','$h_d$','$k_d$']
    nparam = samples.shape[1]
    #nparam = 8
    
    # offset T0 pdfs
    T0b, T0c, T0d = 2458366.1708, 2458367.2752, 2458362.7376
    samples2[:,7] -= T0b
    samples2[:,12] -= T0c
    samples2[:,17] -= T0d
    
    fig = plt.figure(figsize=(nparam*.9, nparam*.9))
    ind = 1
    for i in range(nparam):
        
        print i/float(nparam)

        for j in range(nparam):
            ax = fig.add_subplot(nparam, nparam, ind)
	    ind += 1

            # 2d histograms
            if i > j:
                z,x,y = plt.histogram2d(samples2[:,j], samples2[:,i], bins=bins)
                ax.pcolormesh(x, y, z.T, cmap='Greys')
                dx, dy = np.diff(x).mean()/2, np.diff(y).mean()/2
                levels = np.percentile(z[z != 0], (16,84))
                ax.contour(x[:-1]+dx, y[:-1]+dy, z.T, levels=levels, colors='k')
                #if len(results2) > 0:
                #    ax.plot(results2[j], results2[i], 'o', ms=6,
                #            color='#31a354') #36e2ac
                
            # Plot histogram
            elif i == j:
                y,x,_ = ax.hist(samples2[:,i], bins=bins, histtype='step',
                                color='k', lw=2, normed=True)
                #if i in [0,1,2,3,4]:
                    #kernel = gaussian_kde(samples2[:,i])
                    #xarr = np.linspace(samples2[:,i].min(), samples2[:,i].max(),500)
                    #probs = kernel.pdf(xarr) / kernel.pdf(xarr).sum()
                    #result = float(xarr[probs==probs.max()])
                    #ax.plot(xarr, kernel.pdf(xarr), '-', c='#31a354',lw=1.5)
                #if len(results2) > 0:
                    #ylims = ax.get_ylim()
                    #ax.plot(np.repeat(results2[i],2), list(ylims), 'k-', lw=.9)
                    #if len(sigs) > 0:
                    #ax.fill_between([results2[i]-sigs[i,0],
                    #                 results2[i]+sigs[i,1]], 0, ylims[1],
                    #                color='k', alpha=.2)
                ylims = ax.get_ylim()
                ax.set_ylim(ylims)

            # No axis
            else:
                ax.set_frame_on(False)
                ax.set_xticks([])
                ax.set_yticks([])

            # General customization
            ax.set_yticklabels('')
            if i != nparam-1:
                ax.set_xticklabels('')
            else:
                ticks, ticklabels = get_xticks(j)
                ax.set_xticks(ticks)
                ax.set_xticklabels(ticklabels, fontsize=8)
                if labels != []:
                    ax.set_xlabel(labels[j], fontsize=10)
            ax.set_xlim(get_xlim(j))
            ax.tick_params('both', length=4, width=1, which='major')
            
            # Add y-labels to 2d histograms
            if j == 0 and i > 0:
                ticks, ticklabels = get_xticks(i)
                ax.set_yticks(ticks)
                ax.set_yticklabels(ticklabels, fontsize=8)
                if labels != []:
                    ax.set_ylabel(labels[i], fontsize=10)
            
   
    # save samples
    hdr = 'ln a,ln lambda, ln Gamma, ln P_GP, s, gamma_0, P_b, T0_b, K_b, h_b, k_b, P_c, T0_c, K_c, h_c, k_c, P_d, T0_d, K_d, h_d, k_d'
    samples2[:,7] += T0b-2457e3
    samples2[:,12] += T0c-2457e3
    samples2[:,17] += T0d-2457e3
    inds = np.arange(samples2.shape[0])
    np.random.shuffle(inds)
    inds = inds[:10000]
    np.savetxt('paper/toi175_samples.csv', samples2[inds], fmt='%.6f', delimiter=',', header=hdr)
 
    plt.subplots_adjust(bottom=.04, left=.05, right=.99, top=.99,
                        hspace=.09, wspace=.09)
    if label:
        plt.savefig('/Users/ryancloutier/Research/TOI_175/plots/corner.png',
                    dpi=100)
    if pltt:
        plt.show()
    plt.close('all')


if __name__ == '__main__':
    #print 'hi'
    samples, results = get_result_array()
    #samples = samples[:10000]
    print samples.shape
    plot_corner(samples, results)
