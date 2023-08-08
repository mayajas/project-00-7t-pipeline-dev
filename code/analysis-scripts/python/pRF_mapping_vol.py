#!/usr/bin/env python
# coding: utf-8

# Imports
import numpy as np
import os
import sys
from os.path import join as opj
import scipy.io
import nibabel as nib
import matplotlib.pyplot as plt
from nilearn.maskers import NiftiMasker
from nilearn import image
import pickle
import math 
from scipy.io import savemat
import nipype.interfaces.freesurfer as fs


fit_hrf        = True  # boolean, optional
                              # Whether or not to fit two extra parameters for hrf derivative and
                              # dispersion. The default is False.
model_name     = 'prfpy_vol_fit_hrf_'+str(fit_hrf)

## Set data directories
local = False
if not local:
    slurm_run = True

# get current sub and hem ids
if local or (not local and not slurm_run):
    sub_id        = 0 
    hem_id        = 1
elif not local and slurm_run:
    sub_id        = int(sys.argv[1])
    hem_id        = int(sys.argv[2])

# list of all sub and hem ids
subject_list  = ['sub-01','sub-02','sub-03','sub-04']
hem_list      = ['lh','rh']
hem_text_list = ['left','right']

# directories
proj_id       = 'project-00-7t-pipeline-dev'
if local:
    proj_dir      = '/home/mayajas/scratch/'+proj_id+'/'
    home_dir      = '/home/mayajas/Documents/'+proj_id+'/'
    programs_dir  = '/home/mayajas/Documents/programs/'
else:
    proj_dir      = '/scratch/mayaaj90/'+proj_id+'/'
    home_dir      = '/home/mayaaj90/projects/'+proj_id+'/'
    programs_dir  = '/home/mayaaj90/programs/'
    
print(model_name)
out_dir     = opj(proj_dir,'output',model_name,subject_list[sub_id])

if subject_list[sub_id] == 'sub-04':
    FS_dir       = opj(proj_dir,'derivatives','wf_advanced_skullstrip_sub-04',
                       '_subject_id_'+subject_list[sub_id],'autorecon_pial')
else:
    FS_dir       = opj(proj_dir,'derivatives','wf_advanced_skullstrip',
                       '_subject_id_'+subject_list[sub_id],'autorecon_pial')

# set FS subjects dir
os.environ["SUBJECTS_DIR"] = FS_dir

# add path to prfpy
prfpy_dir = opj(programs_dir,'prfpy-main')
sys.path.append(prfpy_dir)

from prfpy.stimulus import PRFStimulus2D
from prfpy.model import Iso2DGaussianModel
from prfpy.fit import Iso2DGaussianFitter, Extend_Iso2DGaussianFitter

# set working dir
if os.getcwd() != opj(home_dir,'code','analysis-scripts','python'):
    os.chdir(opj(home_dir,'code','analysis-scripts','python'))

# number of cores to use: either set explicitly or base on settings in slurm job file
import os
if local:
    n_procs = 1
else:
    n_procs = int(os.getenv('OMP_NUM_THREADS'))   
print(n_procs)


###########################################################################################
### Input data filenames
# Image filenames
# pRF mapping runs 
n_runs           = 2
bar1_nii_fn      = opj(out_dir,'reg_bar1.nii')
bar2_nii_fn      = opj(out_dir,'reg_bar2.nii')

# mean functional
meanFunc_nii_fn  = opj(out_dir,'reg_meanFunc.nii')

GM_fn            = opj(out_dir,hem_list[hem_id]+'_GM_inflated_ds.nii')
occ_fn           = opj(out_dir,'funcSpaceOccipitalMask.nii')
    
# Output files
bar_mat_fn           = opj(out_dir,hem_list[hem_id]+'_bar.mat')
grid_fit_avg_fn      = opj(out_dir,hem_list[hem_id]+'_grid_fit_avg.pckl')
iterative_fit_avg_fn = opj(out_dir,hem_list[hem_id]+'_iterative_fit_avg.pckl')
pRF_param_avg_fn     = opj(out_dir,hem_list[hem_id]+'_pRF_params_avg.pckl')

polar_map_nii_fn     = opj(out_dir,hem_list[hem_id]+'.pol.nii')
ecc_map_nii_fn       = opj(out_dir,hem_list[hem_id]+'.ecc.nii')
sigma_map_nii_fn     = opj(out_dir,hem_list[hem_id]+'.sigma.nii')

polar_map_mgh_fn     = opj(out_dir,hem_list[hem_id]+'.pol.mgh')
ecc_map_mgh_fn       = opj(out_dir,hem_list[hem_id]+'.ecc.mgh')



###########################################################################################
### Load preprocessed data
# bar data
bar1     = nib.load(bar1_nii_fn)
bar2     = nib.load(bar2_nii_fn)

# Get mean functional
meanFunc = nib.load(meanFunc_nii_fn)

# Get occipital and gray matter masks, combine
occ      = nib.load(occ_fn)
GM       = nib.load(GM_fn)

mask = image.math_img("np.logical_and(img1, img2)", img1=occ, img2=GM)

###########################################################################################
### Make simple masker for 3D data (no time dim)
# Apply simple masker to mean functional
# This will be used later for "unmasking" pRF parameter maps
masker = NiftiMasker(mask_img=mask)
masked_meanFunc = masker.fit_transform(meanFunc)
masked_meanFunc = np.squeeze(masked_meanFunc)

###########################################################################################
## Clean input data
# - apply occipital mask to constrain analysis to occipital pole
# - detrend, standardize, and bandpass filter each functional pRF run
# - average pRF runs
detrend     = True
standardize = 'zscore'
low_pass    = 0.08       # Low pass filters out high frequency signals from our data: 
                         # fMRI signals are slow evolving processes, any high frequency signals 
                         # are likely due to noise 
high_pass   = 0.009      # High pass filters out any very low frequency signals (below 0.009Hz), 
                         # which may be due to intrinsic scanner instabilities
TR          = 3.0        # repetition time (s)

confounds   = None       # could add motion regressors here

# Detrend, standardize, and bandpass filter each functional pRF run, apply occipital mask
bar_masker = NiftiMasker(mask_img=mask, detrend=detrend, standardize=standardize, t_r=TR,
                    low_pass=low_pass, high_pass=high_pass, verbose=True)

masked_bar1_raw = bar_masker.fit_transform(bar1)
masked_bar2_raw = bar_masker.fit_transform(bar2)

masked_bar1 = masked_bar1_raw.T
masked_bar2 = masked_bar2_raw.T

# Average bar runs
avg_bar = (masked_bar1 + masked_bar2)/2

# Save mat file
if not os.path.exists(bar_mat_fn):
    mymat={'bar_data':avg_bar}
    savemat(bar_mat_fn, mymat)

###########################################################################################
# ### Creating stimulus object
# Get pRF stimulus aperture file
# set design mat from aperture file
Ap_file            = os.path.join(home_dir,'code','stim-scripts','apertures','stimulus_bar.mat')
mat                = scipy.io.loadmat(Ap_file)
design_matrix      = mat["stim"]

# Set max eccentricity
# screen size parameters
screen_height_cm   = 12.00
screen_size_cm     = screen_height_cm/2 
screen_distance_cm = 52.0

# calculate max stim ecc
max_ecc            = math.atan(screen_size_cm/screen_distance_cm)
max_ecc_deg        = math.degrees(max_ecc)
max_ecc_deg


# Define stimulus object
prf_stim = PRFStimulus2D(screen_size_cm=screen_size_cm,
                             screen_distance_cm=screen_distance_cm,
                             design_matrix=design_matrix,
                             TR=TR)

###########################################################################################
## PRF fitting (on data averaged across depths)
## Creating Gaussian model and fitter objects
# Iteratively adjust pRF model parameters (x, y position, pRF size), minimizing the residual sum of squared errors between the prediction and data.
## Define two-dimensional isotropic Gaussian pRF model and model fitter

# Input parameters of Iso2DGaussianModel
hrf                = None     # string, list or numpy.ndarray, optional
                              # HRF shape for this Model.
                              # Can be 'direct', which implements nothing (for eCoG or later convolution),
                              # a list or array of 3, which are multiplied with the three spm HRF basis functions,
                              # and an array already sampled on the TR by the user.
                              # (the default is None, which implements standard spm HRF)
filter_predictions = False    # boolean, optional
                              # whether to high-pass filter the predictions, default False
filter_type        = 'sg'

sg_filter_window_length = 201
sg_filter_polyorder     = 3

filter_params      = {'window_length':sg_filter_window_length, 
                      'polyorder':sg_filter_polyorder}
normalize_RFs      = False    # whether or not to normalize the RF volumes (generally not needed).

# Input parameters of Iso2DGaussianFitter
n_jobs             = n_procs  # int, optional
                              # number of jobs to use in parallelization (iterative search), by default 1

# Define 2D iso Gaussian model
gg = Iso2DGaussianModel(stimulus=prf_stim,
                          filter_predictions=filter_predictions,
                          filter_type=filter_type,
                          filter_params=filter_params,
                          normalize_RFs=normalize_RFs)
# Define 2D iso Gaussian model fitter
gf = Iso2DGaussianFitter(data=avg_bar, model=gg, n_jobs=n_jobs, fit_hrf=fit_hrf)


## Grid fit
# First, conduct a quick, coarse model fitting using provided grids and predictor definitions
# Grid fit parameters
grid_nr       = 30
max_ecc_size  = round(max_ecc_deg,2)

size_grid, ecc_grid, polar_grid = max_ecc_size * np.linspace(0.25,1,grid_nr)**2, \
                    max_ecc_size * np.linspace(0.1,1,grid_nr)**2, \
                        np.linspace(0, 2*np.pi, grid_nr)
verbose       = True        # boolean, optional
                            # Whether to print output. The default is False.


# Run grid fit
if not local and slurm_run:
    if not os.path.exists(grid_fit_avg_fn):
        gf.grid_fit(ecc_grid=ecc_grid,
                    polar_grid=polar_grid,
                    size_grid=size_grid,
                    verbose=verbose,
                    n_batches=n_procs)
    else:
        f = open(grid_fit_avg_fn,'rb')
        gf = pickle.load(f)


# Save grid fit result
if not local and slurm_run:
    if not os.path.exists(grid_fit_avg_fn):
        f = open(grid_fit_avg_fn, 'wb')
        pickle.dump(gf, f)
        f.close()


## Iterative fit
# Next, run fine, iterative fit 

# Iterative fit parameters (2D iso Gaussian model)
rsq_thresh_itfit = 0.0005   # float
                            # Rsq threshold for iterative fitting. Must be between 0 and 1.
verbose          = True     # boolean, optional
                            # Whether to print output. The default is False.


# Run iterative fit
if not local and slurm_run:
    if not os.path.exists(iterative_fit_avg_fn):
        gf.iterative_fit(rsq_threshold=rsq_thresh_itfit, verbose=verbose)
    else:
        f = open(iterative_fit_avg_fn,'rb')
        gf = pickle.load(f)


# Save iterative fit result
if not local and slurm_run:
    if not os.path.exists(iterative_fit_avg_fn):
        f = open(iterative_fit_avg_fn, 'wb')
        pickle.dump(gf, f)
        f.close()


## PRF parameter estimates
# Extract pRF parameter estimates from iterative fit result
if not local and slurm_run and not os.path.exists(pRF_param_avg_fn):
    x=gf.iterative_search_params[:,0]
    y=gf.iterative_search_params[:,1]
    sigma=gf.iterative_search_params[:,2]
    if fit_hrf:
        hrf_1=gf.iterative_search_params[:,5]
        hrf_2=gf.iterative_search_params[:,6]
    total_rsq = gf.iterative_search_params[:,-1]

    #Calculate polar angle and eccentricity maps
    polar = np.angle(x + 1j*y)
    ecc = np.abs(x + 1j*y)

    # Save pRF parameters
    f = open(pRF_param_avg_fn, 'wb')
    if fit_hrf:
        pickle.dump([x, y, sigma, total_rsq, polar, ecc, hrf_1, hrf_2], f)
    else:
        pickle.dump([x, y, sigma, total_rsq, polar, ecc], f)
    f.close()

elif os.path.exists(pRF_param_avg_fn):
    f = open(pRF_param_avg_fn,'rb')
    if fit_hrf:
        x, y, sigma, total_rsq, polar, ecc, hrf_1, hrf_2 = pickle.load(f)
    else:
        x, y, sigma, total_rsq, polar, ecc = pickle.load(f)
    f.close()

###########################################################################################
## Save polar angle and eccentricity maps for delineation purposes

# Threshold pRF maps by rsq, constrain to realistic eccentricities & pRF sizes
rsq_thresh = 0.1
pRF_thresh = max_ecc_deg   

# remove bad fits -avg
x[total_rsq<rsq_thresh]     = np.nan
y[total_rsq<rsq_thresh]     = np.nan
sigma[total_rsq<rsq_thresh] = np.nan
polar[total_rsq<rsq_thresh] = np.nan
ecc[total_rsq<rsq_thresh]   = np.nan
if fit_hrf:
    hrf_1[total_rsq<rsq_thresh] = np.nan
    hrf_2[total_rsq<rsq_thresh] = np.nan

# remove vertices where eccentricity is larger than max stimulus ecc
x[ecc>max_ecc_deg]     = np.nan
y[ecc>max_ecc_deg]     = np.nan
polar[ecc>max_ecc_deg] = np.nan
sigma[ecc>max_ecc_deg] = np.nan
ecc[ecc>max_ecc_deg]   = np.nan
if fit_hrf:
    hrf_1[ecc>max_ecc_deg] = np.nan
    hrf_2[ecc>max_ecc_deg] = np.nan

# remove vertices where pRF size is negative
x[sigma<0]     = np.nan
y[sigma<0]     = np.nan
polar[sigma<0] = np.nan
ecc[sigma<0]   = np.nan
sigma[sigma<0] = np.nan
if fit_hrf:
    hrf_1[sigma<0] = np.nan
    hrf_2[sigma<0] = np.nan

# set max pRF size to max stimulus eccentricity
x[sigma>pRF_thresh]     = np.nan
y[sigma>pRF_thresh]     = np.nan
polar[sigma>pRF_thresh] = np.nan
ecc[sigma>pRF_thresh]   = np.nan
sigma[sigma>pRF_thresh] = np.nan
if fit_hrf:
    hrf_1[sigma>pRF_thresh] = np.nan
    hrf_2[sigma>pRF_thresh] = np.nan

# set nans to 0
x[np.isnan(x)]         = 0.
y[np.isnan(y)]         = 0.
polar[np.isnan(polar)] = 0.
ecc[np.isnan(ecc)]     = 0.
sigma[np.isnan(sigma)] = 0.
if fit_hrf:
    hrf_1[np.isnan(sigma)] = np.nan
    hrf_2[np.isnan(sigma)] = np.nan 

# Unmask avg pRF parameters
unmask_x      = masker.inverse_transform(x) 
unmask_y      = masker.inverse_transform(y) 
unmask_sigma  = masker.inverse_transform(sigma) 
unmask_polar    = masker.inverse_transform(polar) 
unmask_ecc    = masker.inverse_transform(ecc) 
unmask_rsq    = masker.inverse_transform(total_rsq) 
if fit_hrf:
    unmask_hrf_1 = masker.inverse_transform(hrf_1) 
    unmask_hrf_2 = masker.inverse_transform(hrf_2) 


# Save maps to .nii files for manual delineations
if not os.path.exists(sigma_map_nii_fn):
    nib.save(unmask_sigma, sigma_map_nii_fn)  
if not os.path.exists(polar_map_nii_fn):
    nib.save(unmask_polar, polar_map_nii_fn) 
if not os.path.exists(ecc_map_nii_fn):
    nib.save(unmask_ecc, ecc_map_nii_fn) 

# Save maps to .mgh files for manual delineations
if not os.path.exists(polar_map_mgh_fn):
    sampler = fs.SampleToSurface(hemi=hem_list[hem_id])
    sampler.inputs.source_file = polar_map_nii_fn
    sampler.inputs.reg_header = True
    sampler.inputs.subjects_dir = FS_dir
    sampler.inputs.subject_id = subject_list[sub_id]
    sampler.inputs.sampling_method = "point"
    sampler.inputs.sampling_range = 0.0
    sampler.inputs.sampling_units = "mm"
    sampler.inputs.surface = 'pial'
    sampler.inputs.out_file = polar_map_mgh_fn
    sampler.run()

if not os.path.exists(ecc_map_mgh_fn):
    sampler = fs.SampleToSurface(hemi=hem_list[hem_id])
    sampler.inputs.source_file = ecc_map_nii_fn
    sampler.inputs.reg_header = True
    sampler.inputs.subjects_dir = FS_dir
    sampler.inputs.subject_id = subject_list[sub_id]
    sampler.inputs.sampling_method = "point"
    sampler.inputs.sampling_range = 0.0
    sampler.inputs.sampling_units = "mm"
    sampler.inputs.surface = 'pial'
    sampler.inputs.out_file = ecc_map_mgh_fn
    sampler.run()



# meanFunc_mgh_nib = nib.freesurfer.mghformat.load(meanFunc_mgh_fn)
# affine = meanFunc_mgh_nib.affine

# if not os.path.exists(polar_map_mgh):
#     nib.save(nib.freesurfer.mghformat.MGHImage(unmask_polar.astype(np.float32, order = "C"),affine=affine),polar_map_mgh)
# if not os.path.exists(ecc_map_mgh):
#     nib.save(nib.freesurfer.mghformat.MGHImage(unmask_ecc.astype(np.float32, order = "C"),affine=affine),ecc_map_mgh)