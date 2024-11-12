"""
OrCAS Program Target Selection Script.

2024/11/12 IJMC: Cleaned up for survey definition paper release.

Licensed for release under the GNU GPL v3.0
"""

from pylab import *
from scipy import stats
import pandas as pd
from scipy import interpolate

# https://github.com/California-Planet-Search/KPF-etc :
from kpf_etc.etc import kpf_photon_noise_estimate, kpf_etc_rv, kpf_etc_snr



plotalot = True # whether to generate plots!


########################################
# First generate 2D selection metric
########################################

# Load the NEXSCI archive:
ncat = pd.read_csv('./tables/PS_2023.08.22_10.46.49_CLEANED_CONFIRMED_OR_VALIDATED_WITH_TSM.csv')
ncat['tsm'] = ncat.pl_tsm
tran = np.logical_and(np.isfinite(ncat.pl_radjerr1 * ncat.pl_radjerr2), ncat.pl_orbsmax<5)
values = np.vstack([np.log10(ncat.pl_orbper)[tran], np.log10(ncat.pl_rade[tran])])
masserr = np.abs(ncat.pl_bmasseerr1 * ncat.pl_bmasseerr2)**0.5
sig = np.logical_and((ncat.pl_bmasse / masserr) > 5, tran)

# Load JWST Targets:
jwst = pd.read_csv('./tables/jwst_cycle_1-and-2_targets.csv')


# Make a Fulton-esque KDE map:
cat = pd.read_csv('./tables/cks_physical_merged.csv')
values = np.vstack([np.log10(cat.koi_period), np.log10(cat.iso_prad)])
valid = np.logical_and(np.logical_and(np.logical_and(np.logical_and(np.logical_and(cat.koi_impact < 0.7, cat.koi_period < 100), cat.iso_srad < 2.1), cat.iso_steff > 4700), cat.iso_steff < 6500), np.isfinite(cat.iso_prad))
values = values[:,valid]

rp_kernel = stats.gaussian_kde(values)
xmin, xmax = values[0].min(), values[0].max()
ymin, ymax = values[1].min(), values[1].max()
rp_X, rp_Y = np.mgrid[xmin:xmax:100j, ymin:ymax:100j]
positions = np.vstack([rp_X.ravel(), rp_Y.ravel()])
rp_Z = np.reshape(rp_kernel(positions).T, rp_X.shape)

values2 = np.vstack([(cat.koi_teq), (cat.iso_prad)])[:,valid]
rp_kernel2 = stats.gaussian_kde(values2)
xmin2, xmax2 = values2[0].min(), values2[0].max()
ymin2, ymax2 = values2[1].min(), values2[1].max()
rp_X2, rp_Y2 = np.mgrid[xmin2:xmax2:100j, ymin2:ymax2:100j]
positions2 = np.vstack([rp_X2.ravel(), rp_Y2.ravel()])
rp_Z2 = np.reshape(rp_kernel2(positions2).T, rp_X2.shape)


xt = [.75, 1, 2, 4, 8, 10, 20]
smallind = np.logical_and(sig, ncat.pl_rade < 4.217)
nepind = np.logical_and(smallind, ncat.pl_rade > 1.8)
nsmall_all = (ncat.pl_rade < 4.217).sum()
nsmall = smallind.sum()
nnep = nepind.sum()
smallgoodind = np.logical_and(smallind, ncat.tsm>30)
nsmallgood = smallgoodind.sum()



def jwst_sample_map(per, rade):
    values = np.vstack([np.log10(per), np.log10(rade)])
    valid = np.isfinite(values.prod(0))
    values = values[:,valid]
    jw_kernel = stats.gaussian_kde(values)
    positions = np.vstack([rp_X.ravel(), rp_Y.ravel()])
    jw_Z = np.reshape(jw_kernel(positions).T, rp_X.shape)
    return(jw_Z)
    
jw_Z =  jwst_sample_map(jwst.pl_orbper, jwst.pl_rade)
diffscale = 1
diffmap = rp_Z - diffscale * jw_Z
neg = diffmap<0
pos = diffmap>0


if plotalot:
    figure(figsize=[13, 4])
    ax1=subplot(131)
    contourf(rp_X, rp_Y, rp_Z, 15, cmap=plt.cm.YlOrRd)
    cb1=colorbar()
    cb1.set_label('Planet Occurrence', fontsize=12)
    ax1.set_ylabel('Planet Radius [$R_\oplus$]', fontsize=14)
    ax1.set_yticklabels( [1, 2, 4, 8, 16])

    ax2=subplot(132)
    contourf(rp_X, rp_Y, jw_Z, 20, cmap=plt.cm.gist_earth_r)
    ax2.axis([0, 2, log10(.75), log10(16)])
    cb2 = colorbar()
    cb2.set_label('JWST Sampling Density', fontsize=12)

    ax3=subplot(133)
    contourf(rp_X, rp_Y, diffmap, 15, cmap=plt.cm.bwr)
    cb3=colorbar()
    cb3.set_label('Over (blue)/Under (red) Representation by JWST', fontsize=10)

    [ax.plot(log10(jwst.pl_orbper), np.log10(jwst.pl_rade), 'h', mew=2, mec='k', mfc='gold', ms=9) for ax in [ax2]]

    for ax in [ax1, ax2, ax3]:
        ax.axis([0, 2, log10(2), log10(16)])
        ax.set_xticks([0,0.5, 1, 1.5, 2])
        ax.set_xticklabels( [1, 3.2, 10, 31.6, 100])
        ax.set_yticks([0,0.3, 0.6, 0.903, 1.204])
        ax.set_xlabel('Orbital Period [days]', fontsize=14)
        ax.axis([0, 2, log10(1), log10(24)])

    [ax.set_yticklabels([]) for ax in [ax2, ax3]]
    [cb.set_ticklabels([]) for cb in [cb1, cb2, cb3]]
    tight_layout()

    ax3.text(1.15, log10(1.2), 'Underrepresented by JWST;\n need more targets here!', horizontalalignment='center', rotation=20)



########################################
# Now consider TOI selection:
########################################


# Now load in TOI catalog:
toi = pd.read_csv('./tables/TOI2023-08-18_aggregate.csv')
valid_toi = np.logical_and(toi.tfopwg_disp=='PC', np.isfinite(toi.TSM))
toi['RA_deg'] = toi.ra
toi['Dec_deg'] = toi.dec



spl = interpolate.RectBivariateSpline(rp_X[:,0], rp_Y[0], diffmap)

# Estimate TOI catalog properties:
toi['mass'] = toi.pl_rade**2.06
toi['mass'][toi.mass>317] = 317.
toi['k_rv'] = 28.4 * toi.mass/317. * toi.st_rad**-0.666667 * (toi.pl_orbper/365.)**-0.3333333
toi['diffmap'] = spl(log10(toi.pl_orbper), log10(toi.pl_rade), grid=False)
if diffscale != 1:
    toi['diffmap'] += 0.5
    toi['diffmap'][toi['diffmap']<0] = 0.
    toi['diffmap'] /= toi['diffmap'].max()


toi['weight'] = toi.TSM * toi.diffmap
toi['weight'] *= np.logical_and(np.isfinite(toi.diffmap), toi.diffmap>0).astype(int)
region_radper = np.logical_and(np.logical_or(toi.pl_rade < 3, toi.pl_orbper>10), toi.pl_rade < 20)
toi.weight[np.logical_not(region_radper)] = 0.
toi.weight[np.logical_not(valid_toi)] = 0.

# Now calculate KPF integration times for all targets:
ai = np.argsort(toi.weight).values
kpf_rms = 0.5
itimes = np.zeros(len(toi)) + 9e9
toi['kpf_itime'] = 9e9

for jj,ii in enumerate(ai[-300:-1]):
    teff = toi.st_teff.values[ii]
    mag = toi.Vmag.values[ii]
    try:
        val = kpf_etc_rv(teff, mag, kpf_rms)
        if not np.isfinite(val): val = 9e9
    except:
        val = 9e9
    itimes[ii] = val
    toi['kpf_itime'][ii] = val
    print('%i' % (jj+1)),
    
toi['kpf_itime'] = itimes

# Now estimate all noise sources:

# Following eq. 7 of Yu et al. 2018,
#   https://ui.adsabs.harvard.edu/abs/2018MNRAS.480L..48Y/abstract
toi['jit_gran_osc'] = 2.0 * 0.63 * (teff/5777.)**2.46 * ((10**toi.st_logg)/27440)**(-0.45)
toi['jit_vsini'] = 0.16*toi.vsini**1.54
toi['total_rv_noise'] = np.sqrt(kpf_rms**2 + toi.jit_gran_osc**2 + toi.jit_vsini**2)
toi['nrv_to_5sig'] = np.vstack((((5/(toi.k_rv / toi.total_rv_noise))**2).values, 30*np.ones(len(toi)))).max(0)
toi['kpf_tottime'] = toi.nrv_to_5sig * (toi.kpf_itime + 120.)
toi['SUR'] = 1000*toi.weight / toi.kpf_tottime


# Identify TOIs we manually decided to drop:
droplist = [1154.01, 1262.01, 1277.02, 1287.01, 1404.01, 1430.01,1434.01, 1437.01, 1438.02, 1440.01, 1443.01, 1445.01, 1451.01,1453.01, 1473.01, 1648.01, 1730.01, 1730.03, 1735.01, 1736.01,  1742.01, 1751.01, 1778.01, 1801.01, 1802.01, 1824.01, 1835.01,  2280.01, 5076.01, 5117.01, 5156.01, 1643.01, 1432.01, 1441.01, 1740.01, 5539.01, 5515.01, 1464.01, 4576.01, 2287.01]


include = np.ones(len(toi), dtype=bool)
for thistoi in droplist:
    include[(toi.toi.astype(int)==int(thistoi)).values] = False


# Select out targets satisfying our criteria:    
potential = np.logical_and(include, np.logical_and(np.isfinite(toi.kpf_tottime), np.logical_and(region_radper, np.logical_and(valid_toi, np.logical_and(np.logical_and(toi.Dec_deg > -30, toi.k_rv < 10), toi.TSM > 10)))))


    
# Rank and print the final sample:
ai2 = np.argsort(toi.SUR[potential])
print(toi.iloc[potential.values.nonzero()[0]][['toi', 'pl_rade', 'pl_orbper', 'Vmag', 'k_rv', 'TSM', 'SUR', 'kpf_tottime']].iloc[ai2[::-1]])


