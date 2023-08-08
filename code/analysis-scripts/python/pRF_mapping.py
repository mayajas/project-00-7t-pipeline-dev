#!/usr/bin/env python
# coding: utf-8

# Population receptive field mapping workflow
import os
import sys
import subprocess
from os.path import join as opj

import numpy as np
import nibabel as nib
import math 
import scipy

import pickle

from fsl.data.freesurfer import loadVertexDataFile
from nilearn import image, surface, plotting, signal
import nipype.interfaces.freesurfer as fs

import matplotlib.pyplot as plt
from matplotlib import animation


fit_hrf        = True  # boolean, optional
                              # Whether or not to fit two extra parameters for hrf derivative and
                              # dispersion. The default is False.
start_from_avg = False   # whether to use avg across depths as starting point for layer fits
model_name     = 'prfpy_fit_hrf_'+str(fit_hrf)+'_start_from_avg_'+str(start_from_avg)

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

# add path to surface tools
surface_tools_dir = opj(programs_dir,'surface_tools','equivolumetric_surfaces')
sys.path.append(surface_tools_dir)

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

# anatomical image
T1_nii_fn        = opj(out_dir,'T1_out.nii')

# Freesurfer mesh filenames
gm_surf_fn        = opj(FS_dir,subject_list[sub_id],'surf',hem_list[hem_id]+'.pial')
wm_surf_fn        = opj(FS_dir,subject_list[sub_id],'surf',hem_list[hem_id]+'.white')

inflated_surf_fn  = opj(FS_dir,subject_list[sub_id],'surf',hem_list[hem_id]+'.inflated')
sulc_surf_fn      = opj(FS_dir,subject_list[sub_id],'surf',hem_list[hem_id]+'.sulc')

# equivolumetric surface output filenames
n_surfs           = 6  # number of equivolumetric surfaces (including pial and white)
n_surfs_str       = str(n_surfs)

equi_surf_fn_list = ['equi0.0.pial',
                     'equi0.2.pial',
                     'equi0.4.pial',
                     'equi0.6.pial',
                     'equi0.8.pial',
                     'equi1.0.pial']

# surface-projected functional runs
meanFunc_mgh_fn  = opj(out_dir,hem_list[hem_id]+'.meanFunc.mgh')
bar_mgh_fn_list   = [[opj(out_dir,hem_list[hem_id]+'.equi'+str(depth)+'.bar'+str(run+1)+'.mgh') 
                      for run in range(0,n_runs)]
                     for depth in range(0,n_surfs)]


# PRF output files
grid_fit_avg_fn      = opj(out_dir,hem_list[hem_id]+'_grid_fit_avg.pckl')
iterative_fit_avg_fn = opj(out_dir,hem_list[hem_id]+'_iterative_fit_avg.pckl')
pRF_param_avg_fn     = opj(out_dir,hem_list[hem_id]+'_pRF_params_avg.pckl')

occ_mask_fn          = opj(out_dir,hem_list[hem_id]+'_occ_mask.pckl')

grid_fit_per_depth_fn      = opj(out_dir,hem_list[hem_id]+'_grid_fit_per_depth.pckl')
iterative_fit_per_depth_fn = opj(out_dir,hem_list[hem_id]+'_iterative_fit_per_depth.pckl')
pRF_param_per_depth_fn     = opj(out_dir,hem_list[hem_id]+'_pRF_params_per_depth.pckl')

polar_map_mgh        = opj(out_dir,hem_list[hem_id]+'.pol.mgh')
ecc_map_mgh          = opj(out_dir,hem_list[hem_id]+'.ecc.mgh')
x_map_mgh            = opj(out_dir,hem_list[hem_id]+'.x.mgh')
y_map_mgh            = opj(out_dir,hem_list[hem_id]+'.y.mgh')
hrf_1_map_mgh        = opj(out_dir,hem_list[hem_id]+'.hrf_1.mgh')
hrf_2_map_mgh        = opj(out_dir,hem_list[hem_id]+'.hrf_2.mgh')


polar_map_per_depth_mgh = [opj(out_dir,hem_list[hem_id]+'.pol.'+str(depth)+'.mgh')
                               for depth in range(0,n_surfs)]

ecc_map_per_depth_mgh   = [opj(out_dir,hem_list[hem_id]+'.ecc.'+str(depth)+'.mgh')
                               for depth in range(0,n_surfs)]

###########################################################################################
### Generate equivolumetric surfaces
output      = hem_list[hem_id]+'.equi'
equivol_path= opj(FS_dir,subject_list[sub_id],'surf',output+'0.0.pial')

if not os.path.exists(equivol_path):
    command = ['python', surface_tools_dir+'/generate_equivolumetric_surfaces.py',gm_surf_fn, wm_surf_fn, n_surfs_str, output,
           '--software', 'freesurfer','--smoothing','0','--subject_id',subject_list[sub_id]]

    # pass the current environment to the subprocess
    env = os.environ.copy()

    # run the command
    result = subprocess.run(command, env=env, stdout=subprocess.PIPE, shell = True)

    # print the output of the command
    print(result.stdout.decode('utf-8'))

    #os.system('python '+surface_tools_dir+'/generate_equivolumetric_surfaces.py --smoothing 0 --software ''freesurfer'' --subject_id ' + subject_list[sub_id] + ' ' + gm_surf_fn  + ' ' + wm_surf_fn  + ' ' + n_surfs_str  + ' ' + output)

###########################################################################################
### Surface-project functional data
# Mean functional
depth = 0

if not os.path.exists(meanFunc_mgh_fn):
    sampler = fs.SampleToSurface(hemi=hem_list[hem_id])
    sampler.inputs.source_file = meanFunc_nii_fn
    sampler.inputs.reg_header = True
    sampler.inputs.subjects_dir = FS_dir
    sampler.inputs.subject_id = subject_list[sub_id]
    sampler.inputs.sampling_method = "point"
    sampler.inputs.sampling_range = 0.0
    sampler.inputs.sampling_units = "mm"
    sampler.inputs.surface = equi_surf_fn_list[depth]
    sampler.inputs.out_file = meanFunc_mgh_fn
    sampler.run()


# Bar runs (iterating over bar run and equivolumetric surface depth)
for depth in range(0,n_surfs):
    for run in range(0,n_runs):
        if not os.path.exists(bar_mgh_fn_list[depth][run]):
            print(bar_mgh_fn_list[depth][run])
            sampler = fs.SampleToSurface(hemi=hem_list[hem_id])
            if run == 0:
                sampler.inputs.source_file = bar1_nii_fn
            elif run == 1:
                sampler.inputs.source_file = bar2_nii_fn
            sampler.inputs.reg_header = True
            sampler.inputs.subjects_dir = FS_dir
            sampler.inputs.subject_id = subject_list[sub_id]
            sampler.inputs.sampling_method = "point"
            sampler.inputs.sampling_range = 0.0
            sampler.inputs.sampling_units = "mm"
            sampler.inputs.surface = equi_surf_fn_list[depth]
            sampler.inputs.out_file = bar_mgh_fn_list[depth][run]
            sampler.run()

###########################################################################################
### Load preprocessed data
# Freesurfer meshes
gm_mesh       = surface.load_surf_mesh(gm_surf_fn) 
wm_mesh       = surface.load_surf_mesh(wm_surf_fn) 
inflated_mesh = surface.load_surf_mesh(inflated_surf_fn) 

# Surface-projected bar data
meanFunc_mgh      = loadVertexDataFile(meanFunc_mgh_fn)
bar_mgh_list   = [[loadVertexDataFile(bar_mgh_fn_list[depth][run]) 
                   for run in range(0,n_runs)] for depth in range(0,n_surfs)]

###########################################################################################
### Make occipital mask
# (based on surface vertex y-coordinate cut-off, including only posterior vertices)
y_coord_cutoff = -25

n_vtx = len(meanFunc_mgh[:])
n_vtx

occ     = np.zeros(n_vtx)
occ[gm_mesh.coordinates[:,1]<y_coord_cutoff]=1.

occ_mask = np.nonzero(occ)[0]

# Save occipital mask coordinates
if not local:
    if not os.path.exists(occ_mask_fn):
        f = open(occ_mask_fn, 'wb')
        pickle.dump([occ_mask,n_vtx], f)
        f.close()

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

# for details, see: https://nilearn.github.io/dev/modules/generated/nilearn.signal.clean.html


# Apply occipital mask
masked_bar_mgh_list   = [[bar_mgh_list[depth][run][occ_mask].T
                          for run in range(0,n_runs)] for depth in range(0,n_surfs)]

bar_mgh_list[0][0].shape

masked_bar_mgh_list[0][0].shape


# Detrend, standardize, and bandpass filter each functional pRF run
filtered_bar_mgh_list  = [[signal.clean(masked_bar_mgh_list[depth][run],
                           confounds=confounds,
                           detrend=detrend, standardize=standardize, 
                           filter='butterworth', low_pass=low_pass, high_pass=high_pass, 
                           t_r=TR)
                           for run in range(0,n_runs)] for depth in range(0,n_surfs)]


# Average over runs
avg_bar_list = [(sum(filtered_bar_mgh_list[depth][:])/len(filtered_bar_mgh_list[depth][:])).T 
                for depth in range(0,n_surfs)]

# Average over depths
avg_bar = sum(avg_bar_list[:])/len(avg_bar_list[:])

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
# PRF fitting (per depth)
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
if start_from_avg:
    gf_per_depth = [Extend_Iso2DGaussianFitter(data=avg_bar_list[depth], 
                                        model=gg, n_jobs=n_jobs, 
                                        previous_gaussian_fitter=gf,
                                        fit_hrf=fit_hrf) for depth in range(0,n_surfs)]
else:
    gf_per_depth = [Iso2DGaussianFitter(data=avg_bar_list[depth], 
                                        model=gg, n_jobs=n_jobs, 
                                        fit_hrf=fit_hrf) for depth in range(0,n_surfs)]

print("Now running grid fit per depth.")

if not local and slurm_run:
    if not os.path.exists(grid_fit_per_depth_fn):
        print('Running grid fit per depth.')
        for depth in range(0,n_surfs):
            print(depth)
            gf_per_depth[depth].grid_fit(ecc_grid=ecc_grid,
                                         polar_grid=polar_grid,
                                         size_grid=size_grid,
                                         verbose=verbose,
                                         n_batches=n_procs)
    else:
        print('Grid fit per depth already run. Loading.')
        f = open(grid_fit_per_depth_fn,'rb')
        gf_per_depth = pickle.load(f)

if not local and slurm_run:
    if not os.path.exists(grid_fit_per_depth_fn):
        print('Saving grid fit per depth.')
        f = open(grid_fit_per_depth_fn, 'wb')
        pickle.dump(gf_per_depth, f)
        f.close()

print("Now running iterative fit per depth.")

if not local and slurm_run:
    if not os.path.exists(iterative_fit_per_depth_fn):
        print('Running iterative fit per depth.')
        for depth in range(0,n_surfs):
            print(depth)
            gf_per_depth[depth].iterative_fit(rsq_threshold=rsq_thresh_itfit, verbose=verbose)
    else:
        print('Iterative fit per depth already run. Loading.')
        f = open(iterative_fit_per_depth_fn,'rb')
        gf_per_depth = pickle.load(f)

if not local and slurm_run:
    if not os.path.exists(iterative_fit_per_depth_fn):
        print('Saving iterative fit per depth.')
        f = open(iterative_fit_per_depth_fn, 'wb')
        pickle.dump(gf_per_depth, f)


## Extract PRF parameter estimates
if not local and slurm_run and not os.path.exists(pRF_param_per_depth_fn):
    x_per_depth = [gf_per_depth[depth].iterative_search_params[:,0] for depth in range(0,n_surfs)]
    y_per_depth = [gf_per_depth[depth].iterative_search_params[:,1] for depth in range(0,n_surfs)]
    sigma_per_depth = [gf_per_depth[depth].iterative_search_params[:,2] for depth in range(0,n_surfs)]
    if fit_hrf:
        hrf_1_per_depth=[gf_per_depth[depth].iterative_search_params[:,5] for depth in range(0,n_surfs)]
        hrf_2_per_depth=[gf_per_depth[depth].iterative_search_params[:,6] for depth in range(0,n_surfs)]
    total_rsq_per_depth = [gf_per_depth[depth].iterative_search_params[:,-1] for depth in range(0,n_surfs)]

    #Calculate polar angle and eccentricity maps
    polar_per_depth = [np.angle(x_per_depth[depth] + 1j*y_per_depth[depth]) for depth in range(0,n_surfs)]
    ecc_per_depth = [np.abs(x_per_depth[depth] + 1j*y_per_depth[depth]) for depth in range(0,n_surfs)]

    # Save pRF parameters
    f = open(pRF_param_per_depth_fn, 'wb')
    if fit_hrf:
        pickle.dump([x_per_depth, y_per_depth, sigma_per_depth, total_rsq_per_depth, polar_per_depth, ecc_per_depth,
                     hrf_1_per_depth,hrf_2_per_depth], f)
    else:
        pickle.dump([x_per_depth, y_per_depth, sigma_per_depth, total_rsq_per_depth, polar_per_depth, ecc_per_depth], f)
    f.close()

elif os.path.exists(pRF_param_per_depth_fn):
    f = open(pRF_param_per_depth_fn,'rb')
    if fit_hrf:
        x_per_depth, y_per_depth, sigma_per_depth, total_rsq_per_depth, polar_per_depth, ecc_per_depth, hrf_1_per_depth,hrf_2_per_depth = pickle.load(f)
    else:
        x_per_depth, y_per_depth, sigma_per_depth, total_rsq_per_depth, polar_per_depth, ecc_per_depth = pickle.load(f)
    f.close()


###########################################################################################
## Save polar angle and eccentricity maps for delineation purposes
f = open(occ_mask_fn,'rb')
occ_mask,n_vtx = pickle.load(f)
f.close()

# Unmask avg pRF parameters
unmask_x               = np.zeros(n_vtx)
unmask_y               = np.zeros(n_vtx)
unmask_sigma           = np.zeros(n_vtx)
unmask_rsq             = np.zeros(n_vtx)
unmask_polar           = np.zeros(n_vtx)
unmask_ecc             = np.zeros(n_vtx)
if fit_hrf:
    unmask_hrf_1       = np.zeros(n_vtx)
    unmask_hrf_2       = np.zeros(n_vtx)

unmask_x[occ_mask]     = x
unmask_y[occ_mask]     = y
unmask_sigma[occ_mask] = sigma
unmask_rsq[occ_mask]   = total_rsq
unmask_polar[occ_mask] = polar
unmask_ecc[occ_mask]   = ecc
if fit_hrf:
    unmask_hrf_1[occ_mask] = hrf_1
    unmask_hrf_2[occ_mask] = hrf_2

# Threshold pRF maps by rsq, constrain to realistic eccentricities & pRF sizes
rsq_thresh = 0.1
pRF_thresh = max_ecc_deg   

# remove bad fits -avg
unmask_x[unmask_rsq<rsq_thresh]     = np.nan
unmask_y[unmask_rsq<rsq_thresh]     = np.nan
unmask_sigma[unmask_rsq<rsq_thresh] = np.nan
unmask_polar[unmask_rsq<rsq_thresh] = np.nan
unmask_ecc[unmask_rsq<rsq_thresh]   = np.nan
if fit_hrf:
    unmask_hrf_1[unmask_rsq<rsq_thresh] = np.nan
    unmask_hrf_2[unmask_rsq<rsq_thresh] = np.nan

# remove vertices where eccentricity is larger than max stimulus ecc
unmask_x[unmask_ecc>max_ecc_deg]     = np.nan
unmask_y[unmask_ecc>max_ecc_deg]     = np.nan
unmask_polar[unmask_ecc>max_ecc_deg] = np.nan
unmask_sigma[unmask_ecc>max_ecc_deg] = np.nan
unmask_ecc[unmask_ecc>max_ecc_deg]   = np.nan
if fit_hrf:
    unmask_hrf_1[unmask_ecc>max_ecc_deg] = np.nan
    unmask_hrf_2[unmask_ecc>max_ecc_deg] = np.nan

# remove vertices where pRF size is negative
unmask_x[unmask_sigma<0]     = np.nan
unmask_y[unmask_sigma<0]     = np.nan
unmask_polar[unmask_sigma<0] = np.nan
unmask_ecc[unmask_sigma<0]   = np.nan
unmask_sigma[unmask_sigma<0] = np.nan
if fit_hrf:
    unmask_hrf_1[unmask_sigma<0] = np.nan
    unmask_hrf_2[unmask_sigma<0] = np.nan

# set max pRF size to max stimulus eccentricity
unmask_x[unmask_sigma>pRF_thresh]     = np.nan
unmask_y[unmask_sigma>pRF_thresh]     = np.nan
unmask_polar[unmask_sigma>pRF_thresh] = np.nan
unmask_ecc[unmask_sigma>pRF_thresh]   = np.nan
unmask_sigma[unmask_sigma>pRF_thresh] = np.nan
if fit_hrf:
    unmask_hrf_1[unmask_sigma>pRF_thresh] = np.nan
    unmask_hrf_2[unmask_sigma>pRF_thresh] = np.nan

# set nans to 0
unmask_x[np.isnan(unmask_x)]         = 0.
unmask_y[np.isnan(unmask_y)]         = 0.
unmask_polar[np.isnan(unmask_polar)] = 0.
unmask_ecc[np.isnan(unmask_ecc)]     = 0.
unmask_sigma[np.isnan(unmask_sigma)] = 0.
if fit_hrf:
    unmask_hrf_1[np.isnan(unmask_sigma)] = np.nan
    unmask_hrf_2[np.isnan(unmask_sigma)] = np.nan

# Save maps to .mgh files for manual delineations
meanFunc_mgh_nib = nib.freesurfer.mghformat.load(meanFunc_mgh_fn)
affine = meanFunc_mgh_nib.affine

if not os.path.exists(polar_map_mgh):
    nib.save(nib.freesurfer.mghformat.MGHImage(unmask_polar.astype(np.float32, order = "C"),affine=affine),polar_map_mgh)
if not os.path.exists(ecc_map_mgh):
    nib.save(nib.freesurfer.mghformat.MGHImage(unmask_ecc.astype(np.float32, order = "C"),affine=affine),ecc_map_mgh)
if not os.path.exists(x_map_mgh):
    nib.save(nib.freesurfer.mghformat.MGHImage(unmask_x.astype(np.float32, order = "C"),affine=affine),x_map_mgh)
if not os.path.exists(y_map_mgh):
    nib.save(nib.freesurfer.mghformat.MGHImage(unmask_y.astype(np.float32, order = "C"),affine=affine),y_map_mgh)
if not os.path.exists(hrf_1_map_mgh):
    nib.save(nib.freesurfer.mghformat.MGHImage(unmask_hrf_1.astype(np.float32, order = "C"),affine=affine),hrf_1_map_mgh)
if not os.path.exists(hrf_2_map_mgh):
    nib.save(nib.freesurfer.mghformat.MGHImage(unmask_hrf_2.astype(np.float32, order = "C"),affine=affine),hrf_2_map_mgh)



###########################################################################################
## Save depth-dependent polar angle and eccentricity maps for vizualization purposes

# per depth
unmask_x_per_depth     = [np.zeros(n_vtx) for depth in range(0,n_surfs)]
unmask_y_per_depth     = [np.zeros(n_vtx) for depth in range(0,n_surfs)]
unmask_sigma_per_depth = [np.zeros(n_vtx) for depth in range(0,n_surfs)]
unmask_rsq_per_depth   = [np.zeros(n_vtx) for depth in range(0,n_surfs)]
unmask_polar_per_depth = [np.zeros(n_vtx) for depth in range(0,n_surfs)]
unmask_ecc_per_depth   = [np.zeros(n_vtx) for depth in range(0,n_surfs)]

for depth in range(0,n_surfs):
    unmask_x_per_depth[depth][occ_mask]     = x_per_depth[depth]
    unmask_y_per_depth[depth][occ_mask]     = y_per_depth[depth]
    unmask_sigma_per_depth[depth][occ_mask] = sigma_per_depth[depth]
    unmask_rsq_per_depth[depth][occ_mask]   = total_rsq_per_depth[depth]
    unmask_polar_per_depth[depth][occ_mask] = polar_per_depth[depth]
    unmask_ecc_per_depth[depth][occ_mask]   = ecc_per_depth[depth]

    # remove bad fits - per depth
    unmask_x_per_depth[depth][unmask_rsq_per_depth[depth]<rsq_thresh]     = np.nan
    unmask_y_per_depth[depth][unmask_rsq_per_depth[depth]<rsq_thresh]     = np.nan
    unmask_sigma_per_depth[depth][unmask_rsq_per_depth[depth]<rsq_thresh] = np.nan
    unmask_polar_per_depth[depth][unmask_rsq_per_depth[depth]<rsq_thresh] = np.nan
    unmask_ecc_per_depth[depth][unmask_rsq_per_depth[depth]<rsq_thresh]   = np.nan

    # remove vertices where eccentricity is larger than max stimulus ecc
    unmask_x_per_depth[depth][unmask_ecc_per_depth[depth]>max_ecc_deg]     = np.nan
    unmask_y_per_depth[depth][unmask_ecc_per_depth[depth]>max_ecc_deg]     = np.nan
    unmask_polar_per_depth[depth][unmask_ecc_per_depth[depth]>max_ecc_deg] = np.nan
    unmask_sigma_per_depth[depth][unmask_ecc_per_depth[depth]>max_ecc_deg] = np.nan
    unmask_ecc_per_depth[depth][unmask_ecc_per_depth[depth]>max_ecc_deg]   = np.nan

    # remove vertices where pRF size is negative
    unmask_x_per_depth[depth][unmask_sigma_per_depth[depth]<0]     = np.nan
    unmask_y_per_depth[depth][unmask_sigma_per_depth[depth]<0]     = np.nan
    unmask_polar_per_depth[depth][unmask_sigma_per_depth[depth]<0] = np.nan
    unmask_ecc_per_depth[depth][unmask_sigma_per_depth[depth]<0]   = np.nan
    unmask_sigma_per_depth[depth][unmask_sigma_per_depth[depth]<0] = np.nan

    # set max pRF size to max stimulus eccentricity
    unmask_x_per_depth[depth][unmask_sigma_per_depth[depth]>pRF_thresh]     = np.nan
    unmask_y_per_depth[depth][unmask_sigma_per_depth[depth]>pRF_thresh]     = np.nan
    unmask_polar_per_depth[depth][unmask_sigma_per_depth[depth]>pRF_thresh] = np.nan
    unmask_ecc_per_depth[depth][unmask_sigma_per_depth[depth]>pRF_thresh]   = np.nan
    unmask_sigma_per_depth[depth][unmask_sigma_per_depth[depth]>pRF_thresh] = np.nan


    if not os.path.exists(polar_map_per_depth_mgh[depth]):
        nib.save(nib.freesurfer.mghformat.MGHImage(unmask_polar_per_depth[depth].astype(np.float32, order = "C"),affine=affine),polar_map_per_depth_mgh[depth])
    if not os.path.exists(ecc_map_per_depth_mgh[depth]):
        nib.save(nib.freesurfer.mghformat.MGHImage(unmask_ecc_per_depth[depth].astype(np.float32, order = "C"),affine=affine),ecc_map_per_depth_mgh[depth])