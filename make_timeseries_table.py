from read_data import *

# get data
p = read_HARPS()
bjd, rv, erv, program, bjdshort, fwhm, bis, Halpha, eHalpha, Hbeta, eHbeta, Hgamma, eHgamma, NaD, eNaD, Sindex, eSindex = p

# add back systemic velocity
rv -= 5678.38

# create header
hdr = '\\begin{table*}[t]\n'
hdr += '\caption{HARPS spectroscopic time series.}\n'
hdr += '\label{tab:timeseries}\n'
hdr += '\centering\n'
hdr += '\small\n'
hdr += '\\begin{tabular}{ccccccccccccccc}\n'
hdr += '\hline\\noalign{\smallskip}\n'
hdr += 'Time & RV & $\sigma_{\\text{RV}}$ & $H\\alpha$ & $\sigma_{H\\alpha}$ & $H\\beta$ & $\sigma_{H\\beta}$ & $H\\gamma$ & $\sigma_{H\\gamma}$ & NaD & $\sigma_{NaD}$ & S-index & $\sigma_{S}$ & FWHM & BIS \\\ \n'
hdr += '$[$BJD - & [\mps{]} & [\mps{]} & $\\times 10^2$ & $\\times 10^2$ & $\\times 10^2$ & $\\times 10^2$ & $\\times 10^2$ & $\\times 10^2$ & $\\times 10^2$ & $\\times 10^2$ & - & - & $[$km s$^{-1}]$ & $[$km s$^{-1}]$ \\\ \n'
hdr += '2,457,000$]$ &&&&&&&&&&&&&& \\\ \n'
hdr += '\hline\\noalign{\smallskip}\n'

for i in range(bjd.size):
    fwhm_entry = '%.4f'%fwhm[bjdshort==bjd[i]] if bjd[i] in bjdshort else '-' 
    bis_entry = '%.4f'%bis[bjdshort==bjd[i]] if bjd[i] in bjdshort else '-'
    p = '' if i < 5 else '%'
    p += '%.6f,%.1f,%.1f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%s,%s'%(bjd[i]-2457e3,rv[i],erv[i],Halpha[i]*1e2,eHalpha[i]*1e2,Hbeta[i]*1e2,eHbeta[i]*1e2,Hgamma[i]*1e2,eHgamma[i]*1e2,NaD[i]*1e2,eNaD[i]*1e2,Sindex[i],eSindex[i],fwhm_entry,bis_entry)
    p = p.replace(',',' & ')
    p += ' \\\ \n'
    hdr += p

# create footer
ftr = '\hline\\noalign{\smallskip}\n'
ftr += '\end{tabular}\n'
ftr += '\end{table*}'

f = open('toi175_timeseries_table.tex','w')
f.write(hdr+ftr)
f.close()
