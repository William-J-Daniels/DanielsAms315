# use cmd argument to choose the dataset
import sys
datafile = sys.argv[1]

# silence console output so I can see gnu parallel progress bar
import os
sys.stdout = open(os.devnull, "w")
sys.stderr = open(os.devnull, "w")

# create a subdir of figures in which to store results
import os
path = "figures/"+sys.argv[1][:6]+"/"
if not os.path.exists(path):
  os.mkdir(path)

# read data in
import pandas as pd
rawdata = pd.read_csv("data/Project1_PartB_data/"+datafile)
# sort for visualization purposes
rawdata.sort_values("x", inplace=True)

# create an initial visualization
import matplotlib.pyplot as plt
fig, ax = plt.subplots(layout="tight")
ax.set_xlabel("x")
ax.set_ylabel("y")

ax.scatter(rawdata["x"], rawdata["y"],
           marker=".", color="black")

fig.savefig(path+"unbinned.jpg")

# create a dataframe with binned values and visualize
import numpy as np
def bin_nearly_repeated(data, col, factor=0.5, dist_factor=1, diagnostics=False):
    """
    A function to bin nearly repeated values, as in
    x    | y         x    | y
    -----+-- becomes -----+--
    1.01 | 2         1.02 | 2
    1.02 | 3         1.02 | 3
    1.03 | 4         1.02 | 4
    Required input:
    - A pandas dataframe
    - The name of the column being binned
    Optional input:
    - The factor by which to scale the threshold; this modifies the mean of the diff and is set to 0.5 by default
    - The factor by which to scale the maximum bin width; this modifies the greatest distance in the selected column and is set to one be default
      Since it's default is one, the default behavior is that there is no maximum
    - A bool controlling whether or not debugging plots are produced
    Returns:
    - A sorted pandas dataframe containing a column of binned values "x" with their corresponding values "y", equal in length to the input dataframe
    """
    data.sort_values(col, inplace=True)

    # break the relevant column into sections where consecutive points are within the threshold from each other
    dx = np.diff(data[col])

    threshold = factor*dx.mean()
    dist_threshold = dist_factor*max(dx)
    
    break_idx = np.argwhere(dx > threshold)[:,0] + 1
    chunks = np.split(np.asarray(data[col]), break_idx)

    #threshold = 0.015
    #cx = np.cumsum(dx) % threshold
    #break_idx = np.argwhere(np.diff(cx) < 0)[:,0] + 1
    #chunks = np.split(np.asarray(data[col]), break_idx)
    
    # make each element of every chunk equal to the average of the chunk
    # take advantage of the fact that np.split returns views
    for chunk in chunks:
      chunk[:] = chunk.mean()
      
    return data

data = rawdata.copy()
data = bin_nearly_repeated(data, "x", 2.5)
data2 = rawdata.copy()

fig, ax = plt.subplots(layout="tight")
ax.set_xlabel("x")
ax.set_ylabel("y")

ax.scatter(data["x"], data["y"],
           marker=".", color="black")

fig.savefig(path+"binned.jpg")

# define an analysis function and a set of transformations
# then check every possible permutation w/ replacement
import scipy.stats as stats
def fit_and_plot(data, xcol, ycol):
    """
    I'll be doing this many times in the analysis so I made it a function
    required input:
    - a pandas dataframe containing the x and y data
    - the string for the dataframe column with independent values
    - the string for the dataframe column with dependent values
    returns:
    - the slope
    - the intercept
    - the R^2 value
    - pvalue of lack of fit test
    - the figure object that the fit and residuals are drawn on
    """
    # useful quantities
    xbar = data[xcol].mean()
    ybar = data[ycol].mean()
    Sxx = np.sum((data[xcol] - xbar)**2)
    Sxy = np.sum((data[xcol] - xbar) * (data[ycol] - ybar))
    Syy = np.sum((data[ycol] - ybar)**2)

    # estimate fit parameters and add fit + residuals to dataframe
    B1 = Sxy / Sxx
    B0 = ybar - B1*xbar
    fitcol = xcol + ycol + "_fit"
    rescol = xcol + ycol + "_res"
    data[fitcol] = np.polyval([B1, B0], data[xcol])
    data[rescol] = data[ycol] - data[fitcol]
    SSres = np.sum(data[rescol]**2)

    # compute r, r squared, degrees of freedom, and mean squared error
    r = Sxy/np.sqrt(Sxx*Syy)
    rsq = (Syy - SSres) / Syy
    n = len(data)
    dof = n - 2
    MSE = np.sqrt(np.sum((data[ycol] - data[fitcol])**2) / dof)

    # compute .95 confidence and prediction intervals and add to dataframe
    t = stats.t.ppf(0.975, dof)
    CI = t * MSE * np.sqrt(1/n + ((data[xcol]-xbar)**2)/Sxx)
    PI = t * MSE * np.sqrt(1 + 1/n + ((data[xcol]-xbar)**2)/Sxx)

    upper_ci_col = xcol + ycol + "_upperCI"
    lower_ci_col = xcol + ycol + "_lowerCI"
    upper_pi_col = xcol + ycol + "_upperPI"
    lower_pi_col = xcol + ycol + "_lowerPI"

    data[upper_ci_col] = data[fitcol] + CI
    data[lower_ci_col] = data[fitcol] - CI
    data[upper_pi_col] = data[fitcol] + PI
    data[lower_pi_col] = data[fitcol] - PI

    # break data into levels to prepare for lack of fit test
    break_idx = np.nonzero(np.diff(data[xcol]))[0] + 1
    levels_x = np.split(np.asarray(data[xcol]), break_idx)
    levels_y = np.split(np.asarray(data[ycol]), break_idx)

    # find pure experimental error and error due to lack of fit
    SSPexp = 0
    SSPexp_dof = 0
    for ly in levels_y:
        lybar = ly.mean()
        SSPexp += np.sum((ly - lybar)**2)
        SSPexp_dof += len(ly) - 1
    SSlack = SSres - SSPexp
    SSlack_dof = n - 2 - SSPexp_dof

    # calculate the statistic and result of the lack of fit test
    F0 = (SSPexp*SSlack_dof) / (SSlack*SSPexp_dof)
    pLOF = stats.f.cdf(F0, SSPexp_dof, SSlack_dof)

    # visualize
    fig, (fax, rax) = plt.subplots(2, layout="tight",
                                   sharex=True)

    # plot the fit
    fax.set_ylabel("$"+ycol+"$")
    fax.scatter(data[xcol], data[ycol],
                marker=".", color="black")
    l = "$\\beta_1 = %.2f$\n$\\beta_0 = %.2f$\n$r^2 = %.2f$\n$r_{s} = %.2f$"\
         %(B1, B0, rsq, r)
    fax.plot(data[xcol], data[fitcol],
             color="red", label=l)
    
    # plot .95 CI and PI
    fax.fill_between(data[xcol], data[upper_ci_col], data[lower_ci_col],
                     color="pink", alpha=0.5,
                     edgecolor="red", linestyle="dashed",
                     label="$95\%$ CI")
    fax.fill_between(data[xcol], data[upper_pi_col], data[lower_pi_col],
                     color="lightpink", alpha=0.5,
                     edgecolor="red", linestyle="dotted",
                     label="$95\%$ PI")

    fax.legend(loc='center left', bbox_to_anchor=(1, 0.5), frameon=False)

    # residuals with mean and 1 sigma error band
    rax.set_ylabel("$"+ycol+"-\hat{y}$")
    rax.set_xlabel("$"+xcol+"$")
    rax.scatter(data[xcol], data[rescol],
                marker=".", color="black")
    res_mean = data[rescol].mean()
    res_std_dev = np.std(data[rescol])
    l = "$\\overline{res} = %.2f$\n$SS_{res} = %.2f$\n$SSP_{exp} = %.2f$\n$SS_{lack} = %.2f$\n$p_{LOF} = %.4f$"\
        %(res_mean, SSres, SSPexp, SSlack, pLOF)
    rax.hlines(res_mean, min(data[xcol]), max(data[xcol]),
               color="red", label=l)
    l = "$\\sigma_{res} = %.2f$" %res_std_dev
    rax.fill_between(data[xcol], res_mean+res_std_dev, res_mean-res_std_dev,
                     color="pink", alpha=0.5,
                     edgecolor="red", linestyle="dashed",
                     label=l)
    rax.legend(loc='center left', bbox_to_anchor=(1, 0.5), frameon=False)

    return B1, B0, rsq, pLOF, fig

def root(data, col, name):
    data[name+col] = np.sqrt(data[col])
def log(data, col, name):
    data[name+col] = np.log(data[col])
def inverse(data, col, name):
    data[name+col] = 1.0/data[col]
def square(data, col, name):
    data[name+col] = data[col]**2
def no_trans(data, col, name):
    pass

trans_names = [
    "Root",
    "Log",
    "Inverse",
    "Square"
]
trans_funcs = [
    root,
    log,
    inverse,
    square
]

for name, func in zip(trans_names, trans_funcs):
    func(data, "x", name)
    func(data, "y", name)
    func(rawdata, "x", name)
    func(rawdata, "y", name)
    func(data2, "x", name)
    func(data2, "y", name)

for name in trans_names:
    data2 = bin_nearly_repeated(data2, name+"x", 2.5)

trans_names.append("")
binned_results = []
raw_results = []
alt_results = []
from itertools import product
for p in product(trans_names, repeat=2):
    thisx = p[0]+"x"
    thisy = p[1]+"y"
    thisp = thisx + "_" + thisy
    binned_results.append([thisp, *fit_and_plot(data, thisx, thisy)])
    raw_results.append([thisp, *fit_and_plot(rawdata, thisx, thisy)])
    alt_results.append([thisp, *fit_and_plot(data2, thisx, thisy)])

# save all plots
for b, r, a in zip(binned_results, raw_results, alt_results):
    b[5].savefig(path+"%s_binned.jpg" %b[0])
    r[5].savefig(path+"%s_raw.jpg" %r[0])
    a[5].savefig(path+"%s_alt.jpg" %a[0])

# check for relationship between r^2 and pLOF
by_rsq = sorted(binned_results, key=lambda binned_results: binned_results[3], reverse=True)
rsq_vals = [x[3] for x in by_rsq]
pLOF_vals = [y[4] for y in by_rsq]
trans = [l[0] for l in by_rsq]

fig, ax = plt.subplots(layout="tight")
ax.set_xlabel("$r^2$")
ax.set_ylabel("$p_{lof}$")
ax.plot(rsq_vals, pLOF_vals,
        marker="x", linestyle="dotted",
        color="black")
for l, x, y in zip(trans, rsq_vals, pLOF_vals):
    ax.annotate(l, (x,y), size=7)#, rotation=15)
ax.hlines(0.05, min(rsq_vals), max(rsq_vals),
          color="yellow", linestyle="dashed")
ax.hlines(0.10, min(rsq_vals), max(rsq_vals),
          color="green", linestyle="dashed")
fig.savefig(path+"rsq_vs_pLOF.jpg")

# check for relationship between r^2 and pLOF
by_rsq2 = sorted(alt_results, key=lambda alt_results: alt_results[3], reverse=True)
rsq_vals2 = [x[3] for x in by_rsq2]
pLOF_vals2 = [y[4] for y in by_rsq2]
trans2 = [l[0] for l in by_rsq2]

fig, ax = plt.subplots(layout="tight")
ax.set_xlabel("$r^2$")
ax.set_ylabel("$p_{lof}$")
ax.plot(rsq_vals2, pLOF_vals2,
        marker="x", linestyle="dotted",
        color="black")
for l, x, y in zip(trans2, rsq_vals2, pLOF_vals2):
    ax.annotate(l, (x,y), size=7)#, rotation=15)
ax.hlines(0.05, min(rsq_vals2), max(rsq_vals2),
          color="yellow", linestyle="dashed")
ax.hlines(0.10, min(rsq_vals2), max(rsq_vals2),
          color="green", linestyle="dashed")
fig.savefig(path+"rsq_vs_pLOF_alt.jpg")