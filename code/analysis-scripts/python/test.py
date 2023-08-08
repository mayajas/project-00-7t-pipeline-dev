#!/usr/bin/env python

import os
import sys
from os.path import join as opj

import numpy as np
import nibabel
import math 
import scipy

import pickle

from fsl.data.freesurfer import loadVertexDataFile
from nilearn import image, surface, plotting, signal
import nipype.interfaces.freesurfer as fs

import matplotlib.pyplot as plt
from matplotlib import animation

import statsmodels.api as sm
from sklearn.model_selection import train_test_split

local = True
if not local:
    slurm_run = False


# In[5]:


# get current sub and hem ids
if local or (not local and not slurm_run):
    sub_id        = 0 
    hem_id        = 0
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
    
prfpy_dir     = opj(proj_dir,'output','prfpy',subject_list[sub_id])

if subject_list[sub_id] == 'sub-04':
    FS_dir       = opj(proj_dir,'derivatives','wf_advanced_skullstrip_sub-04',
                       '_subject_id_'+subject_list[sub_id],'autorecon_pial')
else:
    FS_dir       = opj(proj_dir,'derivatives','wf_advanced_skullstrip',
                       '_subject_id_'+subject_list[sub_id],'autorecon_pial')

# set FS subjects dir
os.environ["SUBJECTS_DIR"] = FS_dir

prfpy_dir = opj(programs_dir,'prfpy-main')
sys.path.append(prfpy_dir)
from prfpy.stimulus import PRFStimulus2D
from prfpy.model import Iso2DGaussianModel
from prfpy.fit import Iso2DGaussianFitter

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


# ### Input data filenames

# Image filenames

# In[50]:


# pRF mapping runs 
n_runs           = 2
bar1_nii_fn      = opj(prfpy_dir,'reg_bar1.nii')
bar2_nii_fn      = opj(prfpy_dir,'reg_bar2.nii')

# mean functional
meanFunc_nii_fn  = opj(prfpy_dir,'reg_meanFunc.nii')

# anatomical image
T1_nii_fn        = opj(prfpy_dir,'T1_out.nii')

# Freesurfer mesh filenames
gm_surf_fn        = opj(FS_dir,subject_list[sub_id],'surf',hem_list[hem_id]+'.pial')
wm_surf_fn        = opj(FS_dir,subject_list[sub_id],'surf',hem_list[hem_id]+'.white')

inflated_surf_fn  = opj(FS_dir,subject_list[sub_id],'surf',hem_list[hem_id]+'.inflated')
sulc_surf_fn      = opj(FS_dir,subject_list[sub_id],'surf',hem_list[hem_id]+'.sulc')

n_surfs_str = str(6)

output      = hem_list[hem_id]+'.equi'
import subprocess

# set the command to be run
# command = ['python', surface_tools_dir+'/generate_equivolumetric_surfaces.py', '--smoothing 0', '--software ''freesurfer''', '--subject_id ''sub-01''',
#            gm_surf_fn, wm_surf_fn, n_surfs_str, output]

command = ['python', surface_tools_dir+'/generate_equivolumetric_surfaces.py',gm_surf_fn, wm_surf_fn, n_surfs_str, output,
           '--software', 'freesurfer','--smoothing','0','--subject_id',subject_list[sub_id]]

# pass the current environment to the subprocess
env = os.environ.copy()

# run the command
result = subprocess.run(command, env=env, stdout=subprocess.PIPE)

# print the output of the command
print(result.stdout.decode('utf-8'))
#subprocess.run(['python',surface_tools_dir+'/generate_equivolumetric_surfaces.py --smoothing 0 --software ''freesurfer'' --subject_id ' + 'sub-01' + ' ' + gm_surf_fn  + ' ' + wm_surf_fn  + ' ' + n_surfs_str  + ' ' + output])

