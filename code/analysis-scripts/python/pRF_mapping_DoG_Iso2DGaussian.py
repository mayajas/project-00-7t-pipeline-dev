#!/usr/bin/env python
# coding: utf-8

# Imports
import numpy as np
from prfpy.stimulus import PRFStimulus2D
from prfpy.model import Iso2DGaussianModel, DoG_Iso2DGaussianModel
from prfpy.fit import Iso2DGaussianFitter, DoG_Iso2DGaussianFitter
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

# get nr of processing threads (based on slurm script)
n_procs = int(os.getenv('OMP_NUM_THREADS'))   
print(n_procs)

# get current sub and hem ids from input
sub_id       = int(sys.argv[1])
hem_id       = int(sys.argv[2])

# list of all sub and hem ids
subject_list = ['sub-01','sub-02','sub-03','sub-04']
hem_list     = ['lh','rh']

# directories
proj_dir     = '/scratch/mayaaj90/project-00-7t-pipeline-dev/'
data_dir     = opj(proj_dir,'output','func','sliceTimeCorr',
                   '_subject_id_'+subject_list[sub_id])
home_dir  = '/home/mayaaj90/projects/project-00-7t-pipeline-dev/'
prfpy_output_dir = opj(proj_dir,'output','prfpy_DoG_Iso2DGaussianModel',subject_list[sub_id])

# check if prfpy dir exists 
# (it should, and it should contain an occipital mask and GM mask)
if not os.path.isdir(prfpy_output_dir):
    os.makedirs(prfpy_output_dir)

# image files
bar1_file    = opj(data_dir,'_sess_id_task-bar_run-01_sess_nr_0_sess_nvol_124','atask-bar_run-01_roi_warp4D.nii')
bar2_file    = opj(data_dir,'_sess_id_task-bar_run-02_sess_nr_1_sess_nvol_124','atask-bar_run-02_roi_warp4D.nii')
meanFunc_file= opj(proj_dir,'output','func','meanFunc',
                   '_subject_id_'+subject_list[sub_id],'merged_func_mcf.nii_mean_reg.nii')
GM_file      = opj(prfpy_output_dir,
                   hem_list[hem_id]+'_GM_funcSpace.nii')
occ_file     = opj(prfpy_output_dir,
                  'funcSpaceOccipitalMask.nii')
    
# set names of output files
grid_fit_file          = opj(prfpy_output_dir,hem_list[hem_id]+'_grid_fit.pckl')
dog_grid_fit_file      = opj(prfpy_output_dir,hem_list[hem_id]+'_grid_fit_dog.pckl')
iterative_fit_file     = opj(prfpy_output_dir,hem_list[hem_id]+'_iterative_fit.pckl')
dog_iterative_fit_file = opj(prfpy_output_dir,hem_list[hem_id]+'_iterative_fit_dog.pckl')
pRF_param_file         = opj(prfpy_output_dir,hem_list[hem_id]+'_pRF_params.pckl')
bar_mat_file           = opj(prfpy_output_dir,hem_list[hem_id]+'_bar.mat')

# set working dir
if os.getcwd() != opj(home_dir,'code','analysis-scripts','python'):
    os.chdir(opj(home_dir,'code','analysis-scripts','python'))

# set design mat from aperture file
Ap_file            = os.path.join(home_dir,'code','stim-scripts','apertures','stimulus_bar.mat')
mat                = scipy.io.loadmat(Ap_file)
design_matrix      = mat["stim"]

# screen size parameters, TR
screen_height_cm   = 12.0
screen_size_cm     = screen_height_cm/2 
screen_distance_cm = 52.0
TR                 = 3.0

# calculate max stim ecc
max_ecc = math.atan(screen_size_cm/screen_distance_cm)
max_ecc_deg = math.degrees(max_ecc)

# Define stimulus object
prf_stim = PRFStimulus2D(screen_size_cm=screen_size_cm,
                             screen_distance_cm=screen_distance_cm,
                             design_matrix=design_matrix,
                             TR=TR)

# Get func data: 2 runs of sweeping bars
bar1 = nib.load(bar1_file)
bar2 = nib.load(bar2_file)

# Get mean functional
meanFunc = nib.load(meanFunc_file)

# Get occipital and gray matter masks, combine
occ = nib.load(occ_file)
GM = nib.load(GM_file)

mask = image.math_img("np.logical_and(img1, img2)", img1=occ, img2=GM)

# Apply mask to bar data (include detrend, standardization, bandpass filtering)
detrend     = True
standardize = True
low_pass    = 0.1
high_pass   = 0.01
verbose     = True
t_r         = TR

bar_masker = NiftiMasker(mask_img=mask, detrend=detrend, standardize=standardize, t_r=t_r,
                    low_pass=low_pass, high_pass=high_pass, verbose=verbose)

masked_bar1_raw = bar_masker.fit_transform(bar1)
masked_bar2_raw = bar_masker.fit_transform(bar2)

masked_bar1 = masked_bar1_raw.T
masked_bar2 = masked_bar2_raw.T

# Average bar runs
bar_data = (masked_bar1 + masked_bar2)/2

# Save mat file
mymat={'bar_data':bar_data}
savemat(bar_mat_file, mymat)

## Simple masker for 3D data (no time dim)
# Apply simple masker to mean functional
# This will be used later for "unmasking" pRF parameter maps
masker = NiftiMasker(mask_img=mask)
masked_meanFunc = masker.fit_transform(meanFunc)
masked_meanFunc = np.squeeze(masked_meanFunc)


################################################################################################################
## Iso 2D Gaussian model fitting
## Creating Gaussian model and fitter objects

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
fit_hrf            = False    # boolean, optional
                              # Whether or not to fit two extra parameters for hrf derivative and
                              # dispersion. The default is False.

# Define 2D iso Gaussian model
gg = Iso2DGaussianModel(stimulus=prf_stim,
                          filter_predictions=filter_predictions,
                          filter_type=filter_type,
                          filter_params=filter_params,
                          normalize_RFs=normalize_RFs)
# Define 2D iso Gaussian model fitter
gf = Iso2DGaussianFitter(data=bar_data, model=gg, n_jobs=n_jobs, fit_css=False)

# Grid fit parameters (2D iso Gaussian model)
grid_nr      = 30
max_ecc_size = round(max_ecc_deg,2)

size_grid, ecc_grid, polar_grid = max_ecc_size * np.linspace(0.25,1,grid_nr)**2, \
                    max_ecc_size * np.linspace(0.1,1,grid_nr)**2, \
                        np.linspace(0, 2*np.pi, grid_nr)
pos_prfs_only=False

# Iterative fit parameters
rsq_thresh_itfit = 0.0005      # float
                            # Rsq threshold for iterative fitting. Must be between 0 and 1.
verbose       = True        # boolean, optional
                            # Whether to print output. The default is False.

# Run grid fit then iterative fit
try:
    try:
        # if iterative search has already been run
        print("Checking if grid and iterative fits have already been run...")
        f = open(iterative_fit_file,'rb')
        gf = pickle.load(f)

        print("Grid and iterative fits have already been run.")
    except IOError:
        print("Iterative fit not yet run.")
        print("Checking if grid fit has been run...")
        # if grid search has already been run
        f = open(grid_fit_file,'rb')
        gf = pickle.load(f)

        print("Grid fit has already been run")

        print("Now running iterative fit")

        gf.iterative_fit(rsq_threshold=rsq_thresh_itfit, verbose=verbose)

        f = open(iterative_fit_file, 'wb')
        pickle.dump(gf, f)

except IOError:
    print("Neither grid nor iterative fit yet run...")
    print("Now running grid fit.")
    
    gf.grid_fit(ecc_grid=ecc_grid,
                    polar_grid=polar_grid,
                    size_grid=size_grid,
                    verbose=verbose,
                    n_batches=n_procs)
    f = open(grid_fit_file, 'wb')
    pickle.dump(gf, f)
    f.close()
          
    print("Finished running grid fit.")
    print("Now running iterative fit.")
          
    gf.iterative_fit(rsq_threshold=rsq_thresh_itfit, verbose=verbose)

    f = open(iterative_fit_file, 'wb')
    pickle.dump(gf, f)
finally:
    f.close()


################################################################################################################
## DoG Iso 2D Gaussian model fitting
## Creating Gaussian model and fitter objects

# Input parameters of DoG Iso2DGaussianModel
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
fit_hrf            = False    # boolean, optional
                              # Whether or not to fit two extra parameters for hrf derivative and
                              # dispersion. The default is False.

# Define 2D iso Gaussian model
gg_dog = DoG_Iso2DGaussianModel(stimulus=prf_stim,
                                filter_predictions=filter_predictions,
                                filter_type=filter_type,
                                filter_params=filter_params,
                                normalize_RFs=normalize_RFs)

# Define 2D iso Gaussian model fitter
gf_dog = DoG_Iso2DGaussianFitter(data=bar_data, model=gg_dog, n_jobs=n_jobs, fit_css=False,
                                previous_gaussian_fitter=gf)

# Grid fit parameters (DoG 2D iso Gaussian model)
grid_nr = 30
max_ecc_size = round(max_ecc_deg,2)
surround_amplitude_grid, surround_size_grid = np.linspace(0.25,1,grid_nr)**2, \
                    3*max_ecc_size * np.linspace(0.25,1,grid_nr)**2,
pos_prfs_only=False

# Iterative fit parameters
rsq_thresh_itfit = 0.0005      # float
                            # Rsq threshold for iterative fitting. Must be between 0 and 1.
verbose       = True        # boolean, optional
                            # Whether to print output. The default is False.

# Run grid fit then iterative fit
try:
    try:
        # if iterative search has already been run
        print("Checking if DoG grid and iterative fits have already been run...")
        f = open(dog_iterative_fit_file,'rb')
        gf_dog = pickle.load(f)

        print("DoG grid and iterative fits have already been run.")
    except IOError:
        print("DoG iterative fit not yet run.")
        print("Checking if DoG grid fit has been run...")
        # if grid search has already been run
        f = open(dog_grid_fit_file,'rb')
        gf_dog = pickle.load(f)

        print("DoG grid fit has already been run")

        print("Now running DoG iterative fit")

        gf_dog.iterative_fit(rsq_threshold=rsq_thresh_itfit, verbose=verbose)

        f = open(dog_iterative_fit_file, 'wb')
        pickle.dump(gf_dog, f)

except IOError:
    print("Neither DoG grid nor iterative fit yet run...")
    print("Now running DoG grid fit.")
    
    gf_dog.grid_fit(surround_amplitude_grid = surround_amplitude_grid,
                surround_size_grid = surround_size_grid,
                rsq_threshold=rsq_thresh_itfit, verbose=True, gaussian_params=gf.gridsearch_params)
    f = open(dog_grid_fit_file, 'wb')
    pickle.dump(gf_dog, f)
    f.close()
          
    print("Finished running DoG grid fit.")
    print("Now running DoG iterative fit.")
          
    gf_dog.iterative_fit(rsq_threshold=rsq_thresh_itfit, verbose=verbose)

    f = open(dog_iterative_fit_file, 'wb')
    pickle.dump(gf_dog, f)
finally:
    f.close()


# Extract parameters from iterative fit result
try:
    # if pRF parameters have already been extracted
    f = open(pRF_param_file,'rb')
    x, y, prf_size, prf_amp, bold_baseline, srf_amp, srf_size, total_rsq, polar, ecc = pickle.load(f)
except IOError:
    x=gf_dog.iterative_search_params[:,0]
    y=gf_dog.iterative_search_params[:,1]
    prf_size=gf_dog.iterative_search_params[:,2]
    prf_amp=gf_dog.iterative_search_params[:,3]
    bold_baseline=gf_dog.iterative_search_params[:,4]
    srf_amp=gf_dog.iterative_search_params[:,5]
    srf_size=gf_dog.iterative_search_params[:,6]
    total_rsq = gf_dog.iterative_search_params[:,-1]
    
    #Calculate polar angle and eccentricity maps
    polar = np.angle(x + 1j*y)
    ecc = np.abs(x + 1j*y)
    
    f = open(pRF_param_file, 'wb')
    pickle.dump([x, y, prf_size, prf_amp, bold_baseline, srf_amp, srf_size, total_rsq, polar, ecc], f)
finally:
    f.close()


# Threshold by rsq
rsq_thresh_viz = 0.1

x[total_rsq<rsq_thresh_viz] = float('nan')
y[total_rsq<rsq_thresh_viz] = float('nan')
prf_size[total_rsq<rsq_thresh_viz] = float('nan')
prf_amp[total_rsq<rsq_thresh_viz] = float('nan')
srf_size[total_rsq<rsq_thresh_viz] = float('nan')
srf_amp[total_rsq<rsq_thresh_viz] = float('nan')
polar[total_rsq<rsq_thresh_viz] = float('nan')
ecc[total_rsq<rsq_thresh_viz] = float('nan')


# Transform back to brain space
unmasked_x        = masker.inverse_transform(x) 
unmasked_y        = masker.inverse_transform(y) 
unmasked_prf_size = masker.inverse_transform(prf_size) 
unmasked_prf_amp  = masker.inverse_transform(prf_amp) 
unmasked_srf_size = masker.inverse_transform(srf_size) 
unmasked_srf_amp  = masker.inverse_transform(srf_amp) 
unmasked_pol      = masker.inverse_transform(polar) 
unmasked_ecc      = masker.inverse_transform(ecc) 
unmasked_rsq      = masker.inverse_transform(total_rsq) 
#unmasked_meanFunc     = masker.inverse_transform(masked_meanFunc) 
unmasked_barAvg = masker.inverse_transform(bar_data.T) 

# Save Nifti images
nib.save(unmasked_x, opj(prfpy_output_dir, hem_list[hem_id]+'_x.nii'))  
nib.save(unmasked_y, opj(prfpy_output_dir, hem_list[hem_id]+'_y.nii'))  
nib.save(unmasked_prf_size, opj(prfpy_output_dir, hem_list[hem_id]+'_prf_size.nii'))  
nib.save(unmasked_prf_amp, opj(prfpy_output_dir, hem_list[hem_id]+'_prf_amp.nii'))  
nib.save(unmasked_srf_size, opj(prfpy_output_dir, hem_list[hem_id]+'_srf_size.nii'))  
nib.save(unmasked_srf_amp, opj(prfpy_output_dir, hem_list[hem_id]+'_srf_amp.nii')) 
nib.save(unmasked_rsq, opj(prfpy_output_dir, hem_list[hem_id]+'_rsq.nii'))  
nib.save(unmasked_pol, opj(prfpy_output_dir, hem_list[hem_id]+'_pol.nii'))  
nib.save(unmasked_ecc, opj(prfpy_output_dir, hem_list[hem_id]+'_ecc.nii'))  

nib.save(unmasked_barAvg, opj(prfpy_output_dir, hem_list[hem_id]+'_bar.nii'))  



