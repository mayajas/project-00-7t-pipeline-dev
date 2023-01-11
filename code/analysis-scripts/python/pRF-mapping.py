#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
from prfpy.stimulus import PRFStimulus2D
from prfpy.model import Iso2DGaussianModel, Norm_Iso2DGaussianModel
from prfpy.fit import Iso2DGaussianFitter, Norm_Iso2DGaussianFitter
from prfpy.rf import gauss2D_iso_cart

import os

from os.path import join as opj

import scipy.io

import nibabel as nib

import matplotlib.pyplot as plt

from nilearn.maskers import NiftiMasker
from nilearn import image

import pickle


#from nipype.interfaces.fsl.utils import ImageMeants


# In[2]:


if os.getcwd() != '/home/mayaaj90/projects/project-00-7t-pipeline-dev/code/analysis-scripts/python':
    os.chdir('/home/mayaaj90/projects/project-00-7t-pipeline-dev/code/analysis-scripts/python')


# In[3]:


os.getcwd()


# In[4]:


n_procs = int(os.getenv('OMP_NUM_THREADS'))   
print(n_procs)


# In[5]:

import sys
sub_id       = int(sys.argv[1])
hem_id       = int(sys.argv[2])

subject_list = ['sub-01','sub-02','sub-03','sub-04']
hem_list     = ['lh','rh']

proj_dir     = '/scratch/mayaaj90/project-00-7t-pipeline-dev/'
data_dir     = opj(proj_dir,'output','func','sliceTimeCorr',
                   '_subject_id_'+subject_list[sub_id])
bar1_file    = opj(data_dir,'_sess_id_task-bar_run-01_sess_nr_0_sess_nvol_124','atask-bar_run-01_roi_warp4D.nii')
bar2_file    = opj(data_dir,'_sess_id_task-bar_run-02_sess_nr_1_sess_nvol_124','atask-bar_run-02_roi_warp4D.nii')

meanFunc_file= opj(proj_dir,'output','func','meanFunc',
                   '_subject_id_'+subject_list[sub_id],'merged_func_mcf.nii_mean_reg.nii')


prfpy_output_dir = opj(proj_dir,'output','prfpy',subject_list[sub_id])
if not os.path.isdir(prfpy_output_dir):
    os.makedirs(prfpy_output_dir)
    
GM_file      = opj(prfpy_output_dir,
                   hem_list[hem_id]+'_GM_funcSpace.nii')
occ_file     = opj(prfpy_output_dir,
                  'funcSpaceOccipitalMask.nii')
    
grid_fit_file      = opj(prfpy_output_dir,hem_list[hem_id]+'_grid_fit.pckl')
iterative_fit_file = opj(prfpy_output_dir,hem_list[hem_id]+'_iterative_fit.pckl')
pRF_param_file     = opj(prfpy_output_dir,hem_list[hem_id]+'_pRF_params.pckl')

#/home/mayajas/scratch/project-00-7t-pipeline-dev/manualcorr/func/_subject_id_sub-01
#output/func/sliceTimeCorr/_subject_id_sub-01/_sess_id_task-bar_run-01_sess_nr_0_sess_nvol_12


# In[6]:


GM_file


# In[7]:


occ_file


# # Creating stimulus object

# Get stimulus aperture

# In[8]:


project_dir  = '/home/mayaaj90/projects/project-00-7t-pipeline-dev/'
Ap_file      = os.path.join(project_dir,'code','stim-scripts','apertures','stimulus_bar.mat')

mat          = scipy.io.loadmat(Ap_file)


# In[9]:


np.shape(mat["stim"])


# Define screen size, distance, TR and aperture matrix

# In[10]:


screen_height_cm   = 12.0
screen_size_cm     = screen_height_cm/2 
screen_distance_cm = 52.0
TR                 = 3.0
design_matrix      = mat["stim"]


# In[11]:


import math 
max_ecc = math.atan(screen_size_cm/screen_distance_cm)
max_ecc_deg = math.degrees(max_ecc)

print("Max eccentricity of stimulus is "+str(round(max_ecc_deg,2)))


# Define stimulus object

# In[12]:


prf_stim = PRFStimulus2D(screen_size_cm=screen_size_cm,
                             screen_distance_cm=screen_distance_cm,
                             design_matrix=design_matrix,
                             TR=TR)


# # Gaussian model fit

# ## Creating Gaussian model and fitter objects

# In[13]:


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


# In[14]:


help(Iso2DGaussianModel)


# In[15]:


#gg = Iso2DGaussianModel(stimulus=prf_stim)

gg = Iso2DGaussianModel(stimulus=prf_stim,
                          filter_predictions=filter_predictions,
                          filter_type=filter_type,
                          filter_params=filter_params,
                          normalize_RFs=normalize_RFs)

#gf = Iso2DGaussianFitter(data=data, model=gg, n_jobs=n_jobs, fit_hrf=fit_hrf)


# Get func data

# In[16]:


bar1 = nib.load(bar1_file)
bar2 = nib.load(bar2_file)

# bar1_data = bar1.get_fdata()
# bar2_data = bar2.get_fdata()


# In[17]:


occ = nib.load(occ_file)


# In[18]:


GM = nib.load(GM_file)


# In[19]:


meanFunc = nib.load(meanFunc_file)


# Do masking

# In[20]:


mask = image.math_img("np.logical_and(img1, img2)", img1=occ, img2=GM)


# In[21]:


# Visualize it as an ROI
from nilearn.plotting import plot_roi
plot_roi(mask,meanFunc)


# Masker parameters

# In[110]:


detrend     = True
standardize = True
low_pass    = 0.08
high_pass   = 0.02
verbose     = True
t_r         = TR


# In[111]:


masker = NiftiMasker(mask_img=mask, detrend=detrend, standardize=standardize, t_r=t_r,
                    low_pass=low_pass, high_pass=high_pass, verbose=verbose)


# In[112]:


masked_bar1_raw = masker.fit_transform(bar1)
masked_bar2_raw = masker.fit_transform(bar2)


# In[113]:


masked_bar1_raw.shape


# In[137]:


a = np.ones((10, 2))

a[3,0]=7
a[4,1]=5

a


# In[138]:


# a=np.reshape(a,(2,10))
a=a.T
#a=np.rot90(a)
#a =  a.reshape((a.shape[1], a.shape[0]))

a


# In[134]:


masked_bar1_rot = np.rot90(masked_bar1_raw)
masked_bar2_rot = np.rot90(masked_bar2_raw)

masked_bar1_reshape =  masked_bar1_raw.reshape((masked_bar1_raw.shape[1], masked_bar1_raw.shape[0]))
masked_bar2_reshape =  masked_bar2_raw.reshape((masked_bar2_raw.shape[1], masked_bar2_raw.shape[0]))


masked_bar1 = masked_bar1_raw.T
masked_bar2 = masked_bar2_raw.T


# In[115]:


masker = NiftiMasker(mask_img=mask)


# In[116]:


masked_meanFunc = masker.fit_transform(meanFunc)

masked_meanFunc = np.squeeze(masked_meanFunc)


# In[124]:


np.shape(masked_bar1_rot)


# In[118]:


#np.shape(masked_meanFunc)


# In[119]:


np.shape(masked_bar1) == np.shape(masked_bar2)


# Plot timecourses

# In[127]:


idx =  100

# And now plot a few of these
import matplotlib.pyplot as plt
plt.figure(figsize=(7, 5))
plt.plot(masked_bar1_rot[idx, :])
plt.plot(masked_bar2_rot[idx, :])
plt.xlabel('Time [TRs]', fontsize=16)
plt.ylabel('Intensity', fontsize=16)
#plt.xlim(0, 150)
plt.subplots_adjust(bottom=.12, top=.95, right=.95, left=.12)


# In[136]:


idx =  100
# And now plot a few of these
import matplotlib.pyplot as plt
plt.figure(figsize=(7, 5))
plt.plot(masked_bar1_raw[:,idx].T)
plt.plot(masked_bar2_raw[:,idx].T)
# plt.plot(masked_bar1[idx,:].T)
# plt.plot(masked_bar2[idx,:].T)
plt.xlabel('Time [TRs]', fontsize=16)
plt.ylabel('Intensity', fontsize=16)
#plt.xlim(0, 150)
plt.subplots_adjust(bottom=.12, top=.95, right=.95, left=.12)


# Average data from 2 sessions

# In[139]:


bar_data = (masked_bar1 + masked_bar2)/2


# Plot averaged data

# In[140]:


# And now plot a few of these
import matplotlib.pyplot as plt
plt.figure(figsize=(7, 5))
plt.plot(np.rot90(bar_data[:1, :]))
plt.xlabel('Time [TRs]', fontsize=16)
plt.ylabel('Intensity', fontsize=16)
#plt.xlim(0, 150)
plt.subplots_adjust(bottom=.12, top=.95, right=.95, left=.12)


# In[141]:


gf = Iso2DGaussianFitter(data=bar_data, model=gg, n_jobs=n_jobs, fit_css=False)


# ## Gaussian grid fit & Gaussian Iterative Fit

# In[54]:


help(Iso2DGaussianFitter.grid_fit)


# In[55]:


help(Iso2DGaussianFitter.iterative_fit)


# Grid fit parameters

# In[56]:


grid_nr = 30
max_ecc_size = round(max_ecc_deg,2)
size_grid, ecc_grid, polar_grid = max_ecc_size * np.linspace(0.25,1,grid_nr)**2, \
                    max_ecc_size * np.linspace(0.1,1,grid_nr)**2, \
                        np.linspace(0, 2*np.pi, grid_nr)
pos_prfs_only=False


# In[57]:


size_grid


# Iterative fit parameters

# In[34]:


rsq_threshold = 0.0005      # float
                            # Rsq threshold for iterative fitting. Must be between 0 and 1.
verbose       = True        # boolean, optional
                            # Whether to print output. The default is False.
# gauss_bounds  = 
# constraints   =
# xtol          =
# ftol          =


# Check if grid/iterative fit were run; if not, run them

# In[75]:


try:
    try:
        # if pRF parameters have already been extracted
        print("Checking if pRF parameters have already been extracted...")
        f = open(pRF_param_file,'rb')
        print("PRF parameters good to go!")
    except IOError:
        try:
            # if iterative search has already been run
            print("pRF parameters not yet extracted")
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

            gf.iterative_fit(rsq_threshold=rsq_threshold, verbose=verbose)

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
          
    gf.iterative_fit(rsq_threshold=rsq_threshold, verbose=verbose)

    f = open(iterative_fit_file, 'wb')
    pickle.dump(gf, f)
finally:
    f.close()


# ## Plot results

# Extract parameters from iterative fit result

# In[26]:


try:
    # if pRF parameters have already been extracted
    f = open(pRF_param_file,'rb')
    x, y, sigma, total_rsq, polar, ecc = pickle.load(f)
except IOError:
    x=gf.iterative_search_params[:,0]
    y=gf.iterative_search_params[:,1]
    sigma=gf.iterative_search_params[:,2]
    total_rsq = gf.iterative_search_params[:,-1]
    
    #Calculate polar angle and eccentricity maps
    polar = np.angle(x + 1j*y)
    ecc = np.abs(x + 1j*y)
    
    f = open(pRF_param_file, 'wb')
    pickle.dump([x, y, sigma, total_rsq, polar, ecc], f)
finally:
    f.close()


# In[27]:


max(ecc)


# In[28]:


x.shape


# In[32]:


type(x)


# In[33]:


#np.squeeze(masked_meanFunc).shape


# In[32]:


type(masked_meanFunc)


# # In[ ]:


# f = open(iterative_fit_file,'rb')
# gf = pickle.load(f)


# # In[40]:


# from prfpy.rf import gauss2D_iso_cart


# # In[42]:


# vf_extent = [-8, 8]
# nr_vf_pix = 200

# prf_space_x, prf_space_y = np.meshgrid(np.linspace(vf_extent[0], vf_extent[1], nr_vf_pix, endpoint=True),
#                                        np.linspace(vf_extent[0], vf_extent[1], nr_vf_pix, endpoint=True))


# In[43]:


# prf_par_names = ['x', 'y', 'size', 'beta', 'baseline', 'rsq']


# In[76]:


# gf.data_var


# In[ ]:





# In[ ]:

# Threshold by rsq

rsq_thresh = 0.1

x[total_rsq<rsq_thresh] = float('nan')
y[total_rsq<rsq_thresh] = float('nan')
sigma[total_rsq<rsq_thresh] = float('nan')
polar[total_rsq<rsq_thresh] = float('nan')
ecc[total_rsq<rsq_thresh] = float('nan')


# Transform back to brain space

# In[39]:


unmasked_x     = masker.inverse_transform(x) 
unmasked_y     = masker.inverse_transform(y) 
unmasked_sigma = masker.inverse_transform(sigma) 
unmasked_pol   = masker.inverse_transform(polar) 
unmasked_ecc   = masker.inverse_transform(ecc) 
unmasked_rsq   = masker.inverse_transform(total_rsq) 


# In[40]:


#unmasked_meanFunc     = masker.inverse_transform(masked_meanFunc) 
unmasked_barAvg     = masker.inverse_transform(bar_data.T) 


# In[41]:



#plt.stairs(counts, bins)


# Visualize x

# In[42]:



# Save Nifti images

# In[56]:


nib.save(unmasked_x, opj(prfpy_output_dir, hem_list[hem_id]+'_x.nii'))  
nib.save(unmasked_y, opj(prfpy_output_dir, hem_list[hem_id]+'_y.nii'))  
nib.save(unmasked_sigma, opj(prfpy_output_dir, hem_list[hem_id]+'_sigma.nii'))  
nib.save(unmasked_rsq, opj(prfpy_output_dir, hem_list[hem_id]+'_rsq.nii'))  
nib.save(unmasked_pol, opj(prfpy_output_dir, hem_list[hem_id]+'_pol.nii'))  
nib.save(unmasked_ecc, opj(prfpy_output_dir, hem_list[hem_id]+'_ecc.nii'))  

nib.save(unmasked_barAvg, opj(prfpy_output_dir, hem_list[hem_id]+'_bar.nii'))  

