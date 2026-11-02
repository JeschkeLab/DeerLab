import numpy as np 
import warnings 
import matplotlib.pyplot as plt
from numpy import cos, sin, sqrt, arccos
from deerlab.bg_models import _strexp,_hom3d
from deerlab.dd_models import dd_gauss, dd_wormgauss
from deerlab.dipolarmodel import _triangleweight
from deerlab import FitResult
import numexpr as ne
import matplotlib.gridspec as gridspec
from scipy.io import loadmat
from scipy import optimize
from scipy.interpolate import interp1d
import scipy.stats as st

from itertools import permutations,combinations



def _kde2d(x,y,color,bins,ax=None):
    """
    Plot a 2D kernel density estimate as filled contours at the 90%, 50%, and 25% highest-density levels.

    Parameters
    ----------
    x, y : array_like
        Sample coordinates of the two variables.
    color : str or RGB tuple
        Fill color for the contour levels.
    bins : array_like
        Values whose min/max set the axis limits and KDE evaluation range.
    ax : matplotlib.axes.Axes, optional
        Axes to draw into. Defaults to the current axes (``plt.gca()``).

    Returns
    -------
    ax : matplotlib.axes.Axes
        The axes the contours were drawn into.
    """
    if ax is None:
        ax = plt.gca()

    xmin, xmax = min(bins),max(bins)
    ymin, ymax = min(bins),max(bins)
    xx, yy = np.mgrid[xmin:xmax:100j, ymin:ymax:100j]
    positions = np.vstack([xx.ravel(), yy.ravel()])
    values = np.vstack([x, y])
    kernel = st.gaussian_kde(values,bw_method=0.6)
    f = np.reshape(kernel(positions).T, xx.shape)

    f = f/np.sum(f)
    n = 2000
    t = np.linspace(0, f.max(), n)
    integral = ((f >= t[:, None, None]) * f).sum(axis=(1,2))

    from scipy import interpolate
    interp = interpolate.interp1d(integral, t, fill_value='extrapolate')
    t_contours = np.sort(interp(np.array([0.9,0.5,0.25])))
    ax.contourf(xx,yy, f, np.array([t_contours[0],t_contours[1]]), extent=[x.min(),x.max(),x.min(),x.max()], colors=color,alpha=0.2)
    ax.contourf(xx,yy, f, np.array([t_contours[1],t_contours[2]]), extent=[x.min(),x.max(),x.min(),x.max()], colors=color,alpha=0.4)
    ax.contourf(xx,yy, f, np.array([t_contours[2],np.max(f)]), extent=[x.min(),x.max(),x.min(),x.max()], colors=color,alpha=0.7)

    ax.set_xlim([xmin,xmax])
    ax.set_ylim([xmin,xmax])
    ax.grid(alpha=0.3)

    return ax

def _contour2d(x,y,P,color,ax=None):
    """
    Plot a gridded 2D distribution as filled contours at the 90%, 50%, and 25% highest-density levels.

    Gridded counterpart of ``_kde2d``, for distributions that are already evaluated on a grid
    (e.g. the 2D projections of an MMMx simulation), so no sampling or kernel density estimate is needed.

    Parameters
    ----------
    x, y : array_like
        Distance axes of the two variables.
    P : array_like
        2D distribution with rows along ``x`` and columns along ``y``, shape ``(len(x), len(y))``.
    color : str or RGB tuple
        Fill color for the contour levels.
    ax : matplotlib.axes.Axes, optional
        Axes to draw into. Defaults to the current axes (``plt.gca()``).

    Returns
    -------
    ax : matplotlib.axes.Axes
        The axes the contours were drawn into.
    """
    from matplotlib.colors import to_rgba
    if ax is None:
        ax = plt.gca()

    P = np.asarray(P,dtype=float)
    P = P/np.sum(P)

    # Find the density thresholds enclosing 90%, 50% and 25% of the probability mass
    Psorted = np.sort(P.ravel())[::-1]
    mass = np.cumsum(Psorted)
    t_contours = [Psorted[np.searchsorted(mass,level)] for level in [0.9,0.5,0.25]]

    # contourf expects Z with shape (len(y), len(x)), hence P.T
    ax.contourf(x, y, P.T, levels=[*t_contours, P.max()],
                colors=[to_rgba(color,0.2),to_rgba(color,0.4),to_rgba(color,0.7)])

    return ax



def plot_multivariate_old(means, cholfactors, model='gauss', xlim=[0,10], MDfile=None, HSCfile=None, 
                      MMMfile=None, MMMparam=None, graphic=None, skiprows=None, saveas=None):
    """
    Generates a multi-variate plot
    
    
    Parameters
    ----------


    
    """

    means,cholfactors = np.atleast_1d(means,cholfactors)

    MDbins = np.linspace(*xlim,150)

    # Monte-Carlo multivariate distance integration
    Nsamples = 10000
    np.random.seed(seed=1)
    if model=='gauss':
        rsamples = Pmultivar(means,cholfactors.tolist()).rvs(Nsamples)
        rsamples = np.maximum(rsamples,1e-16) # Avoid values exactly at zero 
    elif model=='wormgauss':
        # Monte-Carlo multivariate distance integration
        Nsamples = 500000
        rsamples = []
        for contour,persistence,std in zip(means, cholfactors):
            P = lambda r: dl.dd_wormgauss(r,contour,persistence,std) 
            rsamples.append( inverse_sample_function(P, Nsamples, x_min=0, x_max=8, n=300) )
        rsamples = np.array(rsamples).T
        rsamples = np.maximum(rsamples,1e-16) # Avoid values exactly at zero 
        rsamples = np.maximum(rsamples,1e-16) # Avoid values exactly at zero 

    # Triangle inequalities
    Ndistances = np.shape(rsamples)[1]
    idx = np.arange(Ndistances)
    triangle_condition = np.full_like(rsamples,False)
    for n in range(Ndistances):
        idx = np.roll(idx,-1,axis=0)
        triangle_condition[:,n] = np.sum(rsamples[:,idx[:-1]],axis=1) > rsamples[:,idx[-1]]
    triangle_condition = np.all(triangle_condition,axis=1)
    # Discard samples that do not satisfy triangle inequalities
    rsamples = rsamples[triangle_condition,:]

    if MDfile is not None:
        MDdata = np.loadtxt(MDfile,skiprows=skiprows,usecols=[11,12,13])
        MDdata = MDdata/10
        MDmeans = np.median(MDdata,axis=0)
        sorting = match(means,MDmeans,lambda x1, x2: abs(x1 - x2))
        MDdata = MDdata[:,sorting]

    if HSCfile is not None: 
        HSCdata = np.loadtxt(HSCfile)
        HSCmeans = np.median(HSCdata,axis=0)
        sorting = match(means,HSCmeans,lambda x1, x2: abs(x1 - x2))
        HSCdata = HSCdata[:,sorting]

    if MMMfile is not None: 
        MMMdata = loadmat(MMMfile)['rlist']
        MMMmeans = np.median(MMMdata,axis=0)
        sorting = match(means,MMMmeans,lambda x1, x2: abs(x1 - x2))
        MMMdata = MMMdata[:,sorting]

    if MMMparam is not None:
        truth_means,truth_std = MMMparam
        truth_means,truth_std =np.atleast_1d(truth_means,truth_std)
        sorting = match(means,truth_means,lambda x1, x2: abs(x1 - x2))
        truth_means = truth_means[sorting]
        truth_std = truth_std[sorting]
        truth_corr = [truth_std[0], truth_std[1], truth_std[2], 0, 0, 0]

    # Define custom colors
    violet = '#4550e6'

    # Get dimensionality of the distribution
    Ndims = len(means)

    # Construct grids for numerical evaluation
    x = np.linspace(*xlim,200)
    grid1d = np.stack(np.meshgrid(x, indexing='ij'),axis=-1)
    grid2d = np.stack(np.meshgrid(x, x, indexing='ij'),axis=-1)

    nrow,ncol = Ndims,Ndims
    fig, axes = plt.subplots(
        nrow, ncol,
        gridspec_kw=dict(wspace=0.05, hspace=0.05,
                        top=1. - 0.5 / (nrow + 1), bottom=0.5 / (nrow + 1),
                        left=0.5 / (ncol + 1), right=1 - 0.5 / (ncol + 1)),
        figsize=(9, 9),
    )
    plt.rcParams['font.size'] = '13'

    if graphic is not None:
        import matplotlib.image as mpimg
        ax = plt.subplot(Ndims,Ndims,Ndims-1)
        img =  mpimg.imread(graphic)
        plt.imshow(img)
        plt.gca().get_xaxis().set_visible(False)
        plt.gca().get_yaxis().set_visible(False)
        plt.gca().spines['top'].set_visible(False)
        plt.gca().spines['bottom'].set_visible(False)
        plt.gca().spines['right'].set_visible(False)
        plt.gca().spines['left'].set_visible(False)

    # Loop over all subplots
    for dim1 in range(Ndims):
        for dim2 in range(Ndims):
                
            plt.sca(axes[dim1][dim2])

            if dim2>dim1:
                
                plt.gca().get_xaxis().set_visible(False)
                plt.gca().get_yaxis().set_visible(False)
                plt.gca().spines['top'].set_visible(False)
                plt.gca().spines['bottom'].set_visible(False)
                plt.gca().spines['right'].set_visible(False)
                plt.gca().spines['left'].set_visible(False)
                continue 

            # Univariate marginal distribution
            if dim1==dim2:
                #----------------------------------------------------------------------------------
                # Plot simulations
                if MDfile is not None: 
                    kde = st.gaussian_kde(MDdata[:,dim1])
                    plt.fill(MDbins,kde(MDbins),color='grey',alpha=0.5,density=True,label='MD simulation')
                if HSCfile is not None: 
                    kde = st.gaussian_kde(HSCdata[:,dim1])
                    plt.fill(MDbins,kde(MDbins),color='grey',alpha=0.5,label='HSC simulation')
                if MMMfile is not None: 
                    kde = st.gaussian_kde(MMMdata[:,dim1])
                    plt.fill(MDbins,kde(MDbins),color='grey',alpha=0.5,label='MMM simulation')
                if MMMparam is not None: 
                    Σ = cholesky_covmat(truth_corr)
                    truth_pdf = Punivar(truth_means[dim1],0,Σ[dim1,dim1]).pdf(grid1d)
                    plt.fill(x,truth_pdf,color='grey',alpha=0.3,label='MMM simulation',linewidth=0)
                #----------------------------------------------------------------------------------
                kde = st.gaussian_kde(rsamples[:,dim1])
                plt.fill(MDbins,kde(MDbins),color=violet,alpha=0.7,label='Experiment')

                plt.gca().set_xlim(xlim)
                plt.gca().get_xaxis().set_ticks(np.arange(min(xlim)+1,max(xlim)+1,2).tolist())
                plt.gca().get_xaxis().set_ticklabels(np.arange(min(xlim)+1,max(xlim)+1,2).tolist())
                plt.gca().get_yaxis().set_visible(False)
                plt.gca().spines['top'].set_visible(False)
                plt.gca().spines['right'].set_visible(False)
                plt.gca().spines['left'].set_visible(False)
                if dim1==0:
                    plt.legend(frameon=False)
                if dim1<Ndims-1:
                    plt.gca().get_xaxis().set_ticklabels([])
                else:
                    plt.xlabel(f'$r_{dim2+1}$ (nm)')
                plt.rcParams['axes.ymargin'] = 0

            if dim2<dim1:
                #----------------------------------------------------------------------------------
                # Plot simulations
                if MDfile is not None: 
                    subsampling = np.arange(0,MDdata.shape[0],50)
                    _kde2d(MDdata[subsampling,dim2], MDdata[subsampling,dim1],color='grey', bins=MDbins)
                if HSCfile is not None: 
                    subsampling = np.arange(0,HSCdata.shape[0],50)
                    _kde2d(HSCdata[subsampling,dim2], HSCdata[subsampling,dim1],color='grey', bins=MDbins)
                if MMMfile is not None: 
                    subsampling = np.arange(0,MMMdata.shape[0],50)
                    _kde2d(MMMdata[subsampling,dim2], MMMdata[subsampling,dim1],color='grey', bins=MDbins)
                #----------------------------------------------------------------------------------
                plt.grid(alpha=0.2)
                subsampling = np.arange(0,rsamples.shape[0],50)
                _kde2d(rsamples[subsampling,dim2], rsamples[subsampling,dim1],color=violet, bins=MDbins)

                #cont3 = plt.contour(x,x,pdf,alpha=0.7,cmap='binary')
                plt.xlim(xlim)
                plt.ylim(xlim)
                
                plt.gca().get_xaxis().set_ticks(np.arange(min(xlim)+1,max(xlim)+1,2).tolist())
                plt.gca().get_yaxis().set_ticks(np.arange(min(xlim)+1,max(xlim)+1,2).tolist())
                plt.gca().get_xaxis().set_ticklabels(np.arange(min(xlim)+1,max(xlim)+1,2).tolist())
                plt.gca().get_yaxis().set_ticklabels(np.arange(min(xlim)+1,max(xlim)+1,2).tolist())
                if dim2>0:
                    plt.gca().get_yaxis().set_ticklabels([])
                else:
                    plt.ylabel(f'$r_{dim1+1}$ (nm)')
                if dim1<Ndims-1:
                    plt.gca().get_xaxis().set_ticklabels([])
                else:
                    plt.xlabel(f'$r_{dim2+1}$ (nm)')
    
    if saveas is not None:
        plt.savefig(saveas)
    plt.show()


# -----------------------------------------------------------------------------
# Multivariate plot functions
# -----------------------------------------------------------------------------

def _add_param_table(fig,axs,sources,Ndims,fontsize=8):
    """
    Add a table of the means and Cholesky factors in the empty upper-right corner of the mosaic.

    Parameters
    ----------
    sources : dict
        Maps a label (e.g. 'MMM', 'Fit') to a dict with keys 'means' (Ndims,), 'L' (Ndims,Ndims)
        lower Cholesky factor, both in nm, and optionally 'color' for the text of that block.
    """
    # Span the free area: right of the first column, above the last diagonal panel
    pos_r1, pos_r2, pos_rN = axs['r1'].get_position(), axs['r2'].get_position(), axs[f'r{Ndims}'].get_position()
    ax_table = fig.add_axes([pos_r2.x0, axs[f'r{Ndims-1}'].get_position().y0+0.02, pos_rN.x1-pos_r2.x0, pos_r1.y1-axs[f'r{Ndims-1}'].get_position().y0-0.02])
    ax_table.axis('off')

    cell_text,cell_colors = [],[]
    for label,src in sources.items():
        means,L = np.asarray(src['means']),np.asarray(src['L'])
        for i in range(Ndims):
            cell_text.append([label if i == 0 else '', f'$r_{{{i+1}}}$', f'{means[i]:.2f}'] + [f'{L[i,j]:.3f}' if j <= i else '' for j in range(Ndims)])
            cell_colors.append(src.get('color','k'))
    col_labels = ['', '', r'$\langle r \rangle$ (nm)'] + [f'$L_{{i{j+1}}}$ (nm)' for j in range(Ndims)]

    table = ax_table.table(cellText=cell_text,colLabels=col_labels,loc='upper right',cellLoc='center',edges='horizontal')
    table.auto_set_font_size(False)
    table.set_fontsize(fontsize)
    table.scale(1,1.4)
    table.auto_set_column_width(list(range(len(col_labels))))
    for (row,col),cell in table.get_celld().items():
        if row > 0:
            cell.get_text().set_color(cell_colors[row-1])
            # Only draw a separating line between sources
            cell.visible_edges = 'T' if (row-1) % Ndims == 0 else ''
    return ax_table

def _build_mosaic(Ndims): # for Ndim >=3
    """ Builds the mosaic layout for the multivariate plot. Returns a list of lists of axes labels, where '.' indicates an empty space. """
    if Ndims == 3:
        return [['r1','.','.'],
                ['r12','r2','.'],
                ['r13','r23','r3']]
    elif Ndims == 4:
        return [['r1','.','.','.'],
                ['r12','r2','.','.'],
                ['r13','r23','r3','.'],
                ['r14','r24','r34','r4']]
    elif Ndims == 5:
        return [['r1','.','.','.','.'],
                ['r12','r2','.','.','.'],
                ['r13','r23','r3','.','.'],
                ['r14','r24','r34','r4','.'],
                ['r15','r25','r35','r45','r5']]

def prepare_MMMx(filepath) -> dict:
    """
    Loads a trivariate MMMx simulation from a .mat file and computes the mean distances and
    covariance matrix of its moment-matched multivariate Gaussian.

    The file is expected as written by MMMx ``trivariate_distribution``: a 3D array ``trivariate``
    whose axes are, in order, the site1-site3, site1-site2 and site2-site3 distances (r13, r12, r23),
    with the corresponding distance axes ``r_axis_1``, ``r_axis_2`` and ``r_axis_3`` in Ångström.
    The axes are reordered here so that all outputs follow the order (r12, r13, r23).
    The site order is the one given to MMMx (see the .mcx script or the file name), and must be
    matched to the order of the distances in the fitted model before comparison.

    Note that MMMx smooths the distribution of each conformer with a Gaussian kernel, which slightly
    broadens the standard deviations (by about 0.6 Å in quadrature, at the default 0.5 Å resolution)
    and therefore slightly reduces the correlations.

    Parameters
    ----------
    filepath : str, path-like
        Path to the .mat file containing the MMMx simulation.

    Returns
    -------
    dict
        A dictionary containing the following keys:
        - P_2D (dict): 2D probability distributions for each pair of distances, keyed as e.g. 'r12_r13',
          with rows along the first and columns along the second distance.
        - P_1D (dict): 1D marginal distributions for each distance, keyed 'r12', 'r13', 'r23'.
        - means (array): Mean distances (r12, r13, r23), in nm.
        - cov (array): Covariance matrix of (r12, r13, r23), in nm².
        - L (array): Lower Cholesky factor of the covariance matrix (``cov = L@L.T``), in nm, matching
          the ``chol`` parameters of multi-spin ``dipolarmodel`` models (see ``choleskycovmat``).
        - r (dict): Distance axes for each distance, keyed 'r12', 'r13', 'r23', in nm.
    """

    MMM_data = loadmat(filepath,simplify_cells=True)

    if 'trivariate' not in MMM_data:
        raise ValueError("The provided .mat file does not contain a trivariate MMMx simulation. Please provide a valid file.")

    P = np.asarray(MMM_data['trivariate'], dtype=float)
    P /= np.sum(P)  # Normalize the distribution (MMMx stores a density in 1/Å³)

    # MMMx axes of P are (r13, r12, r23), reorder them to (r12, r13, r23); MMMx distance axes are in Å
    P = np.transpose(P, (1,0,2))
    labels = ['r12','r13','r23']
    r = {label: np.ravel(MMM_data[f'r_axis_{k+1}'])/10 for k,label in zip([1,0,2],labels)}

    P_1D = {
        'r12': np.sum(P, axis=(1,2)),
        'r13': np.sum(P, axis=(0,2)),
        'r23': np.sum(P, axis=(0,1))
    }

    P_2D = {
        'r12_r13': np.sum(P, axis=2),
        'r12_r23': np.sum(P, axis=1),
        'r13_r23': np.sum(P, axis=0)
    }

    means = np.array([np.sum(r[label]*P_1D[label]) for label in labels])

    # Calculate the covariance matrix
    cov_matrix = np.diag([np.sum((r[label] - mean)**2*P_1D[label]) for label,mean in zip(labels,means)])
    for i, j in combinations(range(3), 2):
        Pij = P_2D[f'{labels[i]}_{labels[j]}']
        cov_ij = r[labels[i]] @ Pij @ r[labels[j]] - means[i]*means[j]
        cov_matrix[i, j] = cov_ij
        cov_matrix[j, i] = cov_ij  # Symmetric

    L = np.linalg.cholesky(cov_matrix)

    return {
        'Type': 'MMMx',
        'P_2D': P_2D,
        'P_1D': P_1D,
        'means': means,
        'cov': cov_matrix,
        'L': L,
        'r': r
    }

def prepare_HSC(filepath) -> dict:

    raise NotImplementedError("The HSC simulation preparation function is not yet implemented. Please use the MMMx simulation preparation function instead.")



def _extract_means_and_cholesky(fit:FitResult) -> dict:
    """
    Extracts the means and Cholesky factors from a FitResult object for dd_Gauss or dd_WormGauss models.

    Parameters
    ----------
    fit : FitResult
        The fit result object containing the fitted parameters.
    """

    import re
    param_list = fit.paramlist   
    output = {}
    if any(re.match(r'rmean\d+', p) for p in param_list): # model is dd_Gauss
        output['model'] = 'dd_Gauss'
        means = np.array([getattr(fit, p) for p in param_list if re.match(r'rmean\d+', p)])
        L = np.zeros((len(means), len(means)))
        for i in range(len(means)):
            for j in range(i+1):
                L[i,j] = getattr(fit,f'chol{i+1}{j+1}')
        output.update({'means': means, 'L': L})
        return output
    elif any(re.match(r'contour\d+', p) for p in param_list): # model is dd.WormGauss
        output['model'] = 'dd_WormGauss'
        means = np.array([getattr(fit, p) for p in param_list if re.match(r'contour\d+', p)])
        L = np.zeros((len(means), len(means)))
        for i in range(len(means)):
            for j in range(i+1):
                L[i,j] = getattr(fit,f'chol{i+1}{j+1}')
        output.update({'means': means, 'L': L})
        return output
    else:
        raise ValueError("The fit result does not contain the expected parameters for a dd.Gauss or dd.WormGauss model.")

def _inverse_sample_function(dist, pnts, x_min=-100, x_max=100, n=1e5,**kwargs):
    x = np.linspace(x_min, x_max, int(n))
    cumulative = np.cumsum(dist(x, **kwargs))
    cumulative -= cumulative.min()
    f = interp1d(cumulative/cumulative.max(), x)
    return f(np.random.random(int(pnts)))

def _gaussian_marginals(means, L, Npoints=121, nsigma=5):
    """
    Computes the 1D and 2D marginals of a trivariate normal distance distribution truncated by the triangle
    inequality, as in the forward model of multi-spin ``dipolarmodel`` models (see ``_triangleweight``).

    The distribution is evaluated on a grid spanning ``nsigma`` standard deviations around the mean of each
    distance, and returned in the same format as ``prepare_MMMx``, keyed 'r1', 'r2', 'r3' and 'r1_r2' etc.

    Parameters
    ----------
    means : array_like
        Mean distances, in nm.
    L : array_like
        Lower Cholesky factor of the covariance matrix, in nm.
    Npoints : int, optional
        Number of grid points along each distance, by default 121.
    nsigma : float, optional
        Half-width of the grid in standard deviations, by default 5.
    """
    means,L = np.asarray(means,dtype=float),np.asarray(L,dtype=float)
    Ndims = len(means)
    if Ndims != 3:
        raise NotImplementedError("Only trivariate (three-spin) distance distributions are currently supported.")

    Σ = L@L.T
    sigmas = np.maximum(np.sqrt(np.diag(Σ)),1e-3)
    labels = [f'r{q+1}' for q in range(Ndims)]
    r = {label: np.linspace(max(mean-nsigma*sigma,0),mean+nsigma*sigma,Npoints) for label,mean,sigma in zip(labels,means,sigmas)}

    grid = np.stack(np.meshgrid(*r.values(),indexing='ij'),axis=-1)
    P = st.multivariate_normal(mean=means,cov=Σ,allow_singular=True).pdf(grid)
    P *= _triangleweight(*np.moveaxis(grid,-1,0))
    P /= np.sum(P)

    P_1D = {label: np.sum(P,axis=tuple(k for k in range(Ndims) if k != q)) for q,label in enumerate(labels)}
    P_2D = {f'{labels[i]}_{labels[j]}': np.sum(P,axis=tuple(k for k in range(Ndims) if k not in (i,j))) for i,j in combinations(range(Ndims),2)}

    return {'P_1D': P_1D, 'P_2D': P_2D, 'r': r}

def _plot_marginals(axs, marginals, color, label=None, alpha=0.3):
    """
    Plots gridded 1D marginals on the diagonal panels and 2D marginals as contours on the lower-triangle panels.

    Parameters
    ----------
    marginals : dict
        Dictionary with keys 'r', 'P_1D' and 'P_2D' as returned by ``prepare_MMMx`` or ``_gaussian_marginals``.
        The k-th distance is plotted on diagonal panel r{k+1}, and the joint of the i-th and j-th distance on panel r{i+1}{j+1}.
    """
    labels = list(marginals['P_1D'])
    for k,key in enumerate(labels):
        r = marginals['r'][key]
        P1D = marginals['P_1D'][key]
        P1D = P1D/np.trapezoid(P1D,r) # Normalize as a probability density (1/nm), such that different grids can be compared
        axs[f'r{k+1}'].fill_between(r, P1D, color=color, alpha=alpha, linewidth=0)
        axs[f'r{k+1}'].plot(r, P1D, color=color, label=label)

    for i,j in combinations(range(len(labels)),2):
        _contour2d(marginals['r'][labels[i]], marginals['r'][labels[j]], marginals['P_2D'][f'{labels[i]}_{labels[j]}'], color=color, ax=axs[f'r{i+1}{j+1}'])


def plot_multivariate(*,Nspins, fit:FitResult=None, r_lim=(1,8), sim:dict = None, table=True, fig=None):

    """"
    Generates the plot for multivariate distributions of more than two distances. 
    The plot is formed as a triangle of 2D distrubutions showing the correlation between any two distances, and the diagonal shows the 1D marginal distributions of each distance.

    A simulated distribution from either MMMx or HSC can also be added to the plot for comparison, or plotted on its own.

    A table of the mean distances and the cholesky factors of the covariance matrix can also be added to the plot.
    Supports 3 and 4 spin-systems.

    Parameters
    ----------

    nspins: int
        Number of spins in the system. Can be 3 or 4.

    fit: dl.FitResult, optional
        The fit result to be plotted, by default None.

    r_lim: tuple, optional
        The limits of the distance axis, by default (1,8).

    sim: dict, optional
        The MMMx or HSC simulation to be added. This should be created using the
        `prepare_MMMx` or `prepare_HSC` functions, by default None.
    table: bool, optional
        If True, a table of the mean distances and the cholesky factors of the covariance matrix will be added to the plot, by default True. 
    fig: matplotlib.figure.Figure, optional
        If provided, the plot will be added to this figure, by default None.

    Returns
    -------
    fig : matplotlib.figure.Figure
        The figure containing the plot.
    """

    if Nspins < 3 or Nspins > 4:
        raise ValueError("Only spins systems between 3 and 4 are supported.")


    # Build the Figure Layout
    if fig is None:
        fig,axs = plt.subplot_mosaic(
            _build_mosaic(Nspins),
            figsize=(8,8),
            # sharex=True, sharey=True,
            gridspec_kw={'wspace':0,'hspace':0},
        )
    else:
        axs = fig.subplot_mosaic(
            _build_mosaic(Nspins),
            # sharex=True, sharey=True,
            gridspec_kw={'wspace':0,'hspace':0},
        )

    for i in range(Nspins):
        for j in range(i,Nspins):   
            if i == j: # Format the diagonal subplots
                axs[f'r{i+1}'].set_xlim(*r_lim)
                if i<Nspins-1:
                    axs[f'r{i+1}'].axis('off')
                else:
                    axs[f'r{i+1}'].spines['top'].set_visible(False)
                    axs[f'r{i+1}'].spines['right'].set_visible(False)
                    axs[f'r{i+1}'].yaxis.set_visible(False)
                    axs[f'r{i+1}'].set_xlabel(f'$r_{{{i+1}}}$ (nm)')
            if i < j: # Format the lower triangle subplots
                axs[f'r{i+1}{j+1}'].set_xlim(*r_lim)
                axs[f'r{i+1}{j+1}'].set_ylim(*r_lim)
                axs[f'r{i+1}{j+1}'].grid('both')

                if i == 0: #Format the y-axis of the first row of subplots
                    axs[f'r{i+1}{j+1}'].set_ylabel(f'$r_{{{j+1}}}$ (nm)')
                else:
                    axs[f'r{i+1}{j+1}'].set_yticklabels([])

                if j == Nspins-1: #Format the x-axis of the last column of subplots
                    axs[f'r{i+1}{j+1}'].set_xlabel(f'$r_{{{i+1}}}$ (nm)')
                else:
                    axs[f'r{i+1}{j+1}'].set_xticklabels([])

    # Add Simulation if provided
    table_params = {}
    
    if sim is not None:
        # Plot the 1D and 2D marginals of the simulation, the k-th distance (e.g. r12, r13, r23) on diagonal panel r{k+1}
        _plot_marginals(axs, sim, color='grey', label=sim['Type'])
        table_params[sim['Type']] = {'means': sim['means'], 'L': sim['L'], 'color': 'grey'}

    if fit is not None:
        # Extract the means and cholesky factors for dd.Gauss or dd.WormGauss models, using the parameter to idenify the type or model. dd.Gauss has rmean1, rmean2 etc..
        fit_values = _extract_means_and_cholesky(fit)

        if fit_values['model'] == 'dd_Gauss':
            # Multivariate normal distribution truncated by the triangle inequality, as in the dipolarmodel forward model
            fit_marginals = _gaussian_marginals(fit_values['means'], fit_values['L'])
        else:
            raise NotImplementedError("Plotting of multivariate dd_WormGauss fits is not yet supported.")

        _plot_marginals(axs, fit_marginals, color='blue', label='Fit')
        table_params['Fit'] = {'means': fit_values['means'], 'L': fit_values['L'], 'color': 'blue'}

    # The 1D marginals are densities, start their axes at zero
    for k in range(Nspins):
        axs[f'r{k+1}'].set_ylim(bottom=0)

    if table and table_params:
        _add_param_table(fig,axs,table_params,Nspins)

    return fig

