import numpy as np

import sys
sys.path.append("/home/mayajas/Documents/programs/prfpy-main")

import os
from os.path import join as opj

import nibabel as nib

import matplotlib.pyplot as plt

import nilearn
from nilearn import plotting, datasets, image
from nilearn.input_data import NiftiMasker
from nilearn.plotting import plot_roi

import pickle
import pandas as pd

import itertools

import math

#################################################################################################################
# Important variables 
# layer space (either anat or func): space in which laynii layers are computed
# prf model (either prfpy_vol_fit_hrf_False or prfpy_vol_fit_hrf_True
layer_space    = 'func'
fit_hrf        = True  # boolean, optional
                              # Whether or not to fit two extra parameters for hrf derivative and
                              # dispersion. The default is False.
prf_model      = 'prfpy_vol_fit_hrf_'+str(fit_hrf)

#################################################################################################################
# Set max eccentricity
# screen size parameters
screen_height_cm   = 12.00
screen_size_cm     = screen_height_cm/2 
screen_distance_cm = 52.0

# calculate max stim ecc
max_ecc            = math.atan(screen_size_cm/screen_distance_cm)
max_ecc_deg        = math.degrees(max_ecc)

#################################################################################################################
# Define important variables
subject_list = ['sub-01','sub-02','sub-03','sub-04']
hem_list     = ['lh','rh']
roi_list     = ['V1','V2','V3','V4','V1d','V1v','V2d','V2v','V3d','V3v']

rsq_thresh = 0.1

n_layers   = 6
layers     = np.arange(1,n_layers+1)
n_sub      = len(subject_list)
n_hem      = len(hem_list)
n_roi      = len(roi_list)


#################################################################################################################
## Define data frame to store pRF values (pRF size x depth)

depth = np.linspace(1,n_layers,n_layers).tolist()*n_hem*n_sub*n_roi

rois = []
for roi in roi_list:
    rois.append([roi]*n_layers*n_hem)
rois = list(itertools.chain(*rois))*n_sub

hems = [['lh']*n_layers,['rh']*n_layers]
hems = list(itertools.chain(*hems))*n_sub*n_roi

subs = [['sub-01']*n_layers*n_hem*n_roi,['sub-02']*n_layers*n_hem*n_roi,['sub-03']*n_layers*n_hem*n_roi,['sub-04']*n_layers*n_hem*n_roi]
subs = list(itertools.chain(*subs))

pRF_size = np.empty((1,len(depth),)).tolist()
pRF_size = list(itertools.chain(*pRF_size))

eccs = np.empty((1,len(depth),)).tolist()
eccs = list(itertools.chain(*eccs))

hrf_1s = np.empty((1,len(depth),)).tolist()
hrf_1s = list(itertools.chain(*hrf_1s))

hrf_2s = np.empty((1,len(depth),)).tolist()
hrf_2s = list(itertools.chain(*hrf_2s))

rsq_threshs = np.empty((1,len(depth),)).tolist()
rsq_threshs = list(itertools.chain(*rsq_threshs))

# equivol x depth
df_equivol_per_depth = pd.DataFrame({
    'sub id' : subs,
    'roi' : rois,
    'hem' : hems,
    'layer' : depth,
    'pRF size' : pRF_size,
    'ecc' : eccs,
    'hrf1' : hrf_1s,
    'hrf2' : hrf_2s,
    'rsq' : rsq_threshs
})

#################################################################################################################
## Define data frame to store pRF values (pRF size x depth) @ 2deg iso-ecc

depth = np.linspace(1,n_layers,n_layers).tolist()*n_hem*n_sub*n_roi

rois = []
for roi in roi_list:
    rois.append([roi]*n_layers*n_hem)
rois = list(itertools.chain(*rois))*n_sub

hems = [['lh']*n_layers,['rh']*n_layers]
hems = list(itertools.chain(*hems))*n_sub*n_roi

subs = [['sub-01']*n_layers*n_hem*n_roi,['sub-02']*n_layers*n_hem*n_roi,['sub-03']*n_layers*n_hem*n_roi,['sub-04']*n_layers*n_hem*n_roi]
subs = list(itertools.chain(*subs))

pRF_size = np.empty((1,len(depth),)).tolist()
pRF_size = list(itertools.chain(*pRF_size))

eccs = np.empty((1,len(depth),)).tolist()
eccs = list(itertools.chain(*eccs))

hrf_1s = np.empty((1,len(depth),)).tolist()
hrf_1s = list(itertools.chain(*hrf_1s))

hrf_2s = np.empty((1,len(depth),)).tolist()
hrf_2s = list(itertools.chain(*hrf_2s))

rsq_threshs = np.empty((1,len(depth),)).tolist()
rsq_threshs = list(itertools.chain(*rsq_threshs))

# equivol x depth
df_equivol_per_depth_ecc2 = pd.DataFrame({
    'sub id' : subs,
    'roi' : rois,
    'hem' : hems,
    'layer' : depth,
    'pRF size' : pRF_size,
    'ecc' : eccs,
    'hrf1' : hrf_1s,
    'hrf2' : hrf_2s,
    'rsq' : rsq_threshs
})

#################################################################################################################
## Loop over subjects and hemispheres
for sub_id in range(0,len(subject_list)):
    for hem_id in range(0,len(hem_list)):
        ###########################################################################################
        ## Define input data directories and filenames
        # data directories
        prfpy_dir    = '/home/mayajas/scratch/project-00-7t-pipeline-dev/output/'+prf_model+'/'+subject_list[sub_id]
        proj_dir     = '/home/mayajas/scratch/project-00-7t-pipeline-dev/'

        if layer_space == 'func':
            lay_dir      = opj(prfpy_dir,'layerification_func')
        else:
            lay_dir      = opj(prfpy_dir,'layerification')

        # image files
        meanFunc_fn  = opj(prfpy_dir,'meanFunc.nii')
        GM_fn        = opj(prfpy_dir,
                        hem_list[hem_id]+'_GM_inflated_ds.nii')
        occ_fn       = opj(prfpy_dir,
                        'funcSpaceOccipitalMask.nii')

        # prfpy output files
        pRF_param_fn     = opj(prfpy_dir,hem_list[hem_id]+'_pRF_params_avg.pckl')

        # manual ROI delineations
        V1_fn        = opj(prfpy_dir,'func_'+hem_list[hem_id]+'.V1.nii')
        V2_fn        = opj(prfpy_dir,'func_'+hem_list[hem_id]+'.V2.nii')
        V3_fn        = opj(prfpy_dir,'func_'+hem_list[hem_id]+'.V3.nii')
        V4_fn        = opj(prfpy_dir,'func_'+hem_list[hem_id]+'.V4.nii')
        
        ecc2_fn      = opj(prfpy_dir,'func_'+hem_list[hem_id]+'.ecc2.nii')

        # layers
        lay_equivol_fn = opj(lay_dir,'func_ribbon_rim_layers_equivol.nii')

        ###########################################################################################
        ## Load image files
        occ = nib.load(occ_fn)
        GM = nib.load(GM_fn)
        meanFunc = nib.load(meanFunc_fn)

        V1=nib.load(V1_fn)
        V2=nib.load(V2_fn)
        V3=nib.load(V3_fn)
        V4=nib.load(V4_fn)

        ecc2=nib.load(ecc2_fn)

        lay_equivol=nib.load(lay_equivol_fn)        

        # define mask (GM+occipital) and process bar data and masker
        mask = image.math_img("np.logical_and(img1, img2)", img1=occ, img2=GM)
        masker = NiftiMasker(mask_img=mask)

        # mask input files
        masked_meanFunc = masker.fit_transform(meanFunc)
        masked_meanFunc = np.squeeze(masked_meanFunc)

        masked_V1 = masker.fit_transform(V1)
        masked_V1 = np.squeeze(masked_V1)

        masked_V2 = masker.fit_transform(V2)
        masked_V2 = np.squeeze(masked_V2)

        masked_V3 = masker.fit_transform(V3)
        masked_V3 = np.squeeze(masked_V3)

        masked_V4 = masker.fit_transform(V4)
        masked_V4 = np.squeeze(masked_V4)

        masked_ecc2 = masker.fit_transform(ecc2)
        masked_ecc2 = np.squeeze(masked_ecc2)

        masked_lay_equivol = masker.fit_transform(lay_equivol)
        masked_lay_equivol = np.squeeze(masked_lay_equivol)

        ###########################################################################################
        ## Load & filter pRF results
        f = open(pRF_param_fn,'rb')
        if fit_hrf:
            x, y, sigma, total_rsq, polar, ecc, hrf_1, hrf_2 = pickle.load(f)
        else:
            x, y, sigma, total_rsq, polar, ecc = pickle.load(f)
        f.close()

        # Threshold pRF maps by rsq, constrain to realistic eccentricities & pRF sizes
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
        x[sigma>max_ecc_deg]     = np.nan
        y[sigma>max_ecc_deg]     = np.nan
        polar[sigma>max_ecc_deg] = np.nan
        ecc[sigma>max_ecc_deg]   = np.nan
        sigma[sigma>max_ecc_deg] = np.nan
        if fit_hrf:
            hrf_1[sigma>max_ecc_deg] = np.nan
            hrf_2[sigma>max_ecc_deg] = np.nan


        ###########################################################################################
        ## Extract sigma across eccentricities
        depth_id=0
        for layer in layers:
            for roi_id in range(0,len(roi_list)):
                if roi_id == 0:
                    masked_ROI = masked_V1
                elif roi_id == 1:
                    masked_ROI = masked_V2
                elif roi_id == 2:
                    masked_ROI = masked_V3
                elif roi_id == 3:
                    masked_ROI = masked_V4
                sigmaxlayer = sigma[(masked_ROI==1) & (masked_lay_equivol == layer)]
                eccxlayer = ecc[(masked_ROI==1) & (masked_lay_equivol == layer)]
                rsqxlayer = total_rsq[(masked_ROI==1) & (masked_lay_equivol == layer)]
                if fit_hrf:
                    hrf1xlayer  = hrf_1[(masked_ROI==1) & (masked_lay_equivol == layer)]
                    hrf2xlayer  = hrf_2[(masked_ROI==1) & (masked_lay_equivol == layer)]

                
                idx=df_equivol_per_depth.loc[(df_equivol_per_depth['sub id'] == subject_list[sub_id]) & 
                                        (df_equivol_per_depth['hem'] == hem_list[hem_id]) & 
                                        (df_equivol_per_depth['layer'] == depth[depth_id]) &
                                        (df_equivol_per_depth['roi'] == roi_list[roi_id])].index.tolist()
                df_equivol_per_depth.loc[idx, 'pRF size'] = pd.Series([sigmaxlayer]*len(idx), index=df_equivol_per_depth.index[idx])
                df_equivol_per_depth.loc[idx, 'ecc'] = pd.Series([eccxlayer]*len(idx), index=df_equivol_per_depth.index[idx])
                if fit_hrf:
                    df_equivol_per_depth.loc[idx, 'hrf1'] = pd.Series([hrf1xlayer]*len(idx), index=df_equivol_per_depth.index[idx])
                    df_equivol_per_depth.loc[idx, 'hrf2'] = pd.Series([hrf2xlayer]*len(idx), index=df_equivol_per_depth.index[idx]) 
                df_equivol_per_depth.loc[idx, 'rsq'] = pd.Series([rsqxlayer]*len(idx), index=df_equivol_per_depth.index[idx])
                
            depth_id += 1
        
        ###########################################################################################
        ## Extract sigma at manually-defined iso-ecc
        depth_id=0
        for layer in layers:
            for roi_id in range(0,len(roi_list)):
                if roi_id == 0:
                    masked_ROI = masked_V1
                elif roi_id == 1:
                    masked_ROI = masked_V2
                elif roi_id == 2:
                    masked_ROI = masked_V3
                elif roi_id == 3:
                    masked_ROI = masked_V4
                sigmaxlayer = sigma[(masked_ROI==1) & (masked_lay_equivol == layer) & (masked_ecc2==1)]
                eccxlayer = ecc[(masked_ROI==1) & (masked_lay_equivol == layer) & (masked_ecc2==1)]
                rsqxlayer = total_rsq[(masked_ROI==1) & (masked_lay_equivol == layer) & (masked_ecc2==1)]
                if fit_hrf:
                    hrf1xlayer  = hrf_1[(masked_ROI==1) & (masked_lay_equivol == layer) & (masked_ecc2==1)]
                    hrf2xlayer  = hrf_2[(masked_ROI==1) & (masked_lay_equivol == layer) & (masked_ecc2==1)]

                idx=df_equivol_per_depth_ecc2.loc[(df_equivol_per_depth_ecc2['sub id'] == subject_list[sub_id]) & 
                                        (df_equivol_per_depth_ecc2['hem'] == hem_list[hem_id]) & 
                                        (df_equivol_per_depth_ecc2['layer'] == depth[depth_id]) &
                                        (df_equivol_per_depth_ecc2['roi'] == roi_list[roi_id])].index.tolist()
                df_equivol_per_depth_ecc2.loc[idx, 'pRF size'] = pd.Series([sigmaxlayer]*len(idx), index=df_equivol_per_depth_ecc2.index[idx])
                df_equivol_per_depth_ecc2.loc[idx, 'ecc'] = pd.Series([eccxlayer]*len(idx), index=df_equivol_per_depth_ecc2.index[idx])
                if fit_hrf:
                    df_equivol_per_depth_ecc2.loc[idx, 'hrf1'] = pd.Series([hrf1xlayer]*len(idx), index=df_equivol_per_depth_ecc2.index[idx])
                    df_equivol_per_depth_ecc2.loc[idx, 'hrf2'] = pd.Series([hrf2xlayer]*len(idx), index=df_equivol_per_depth_ecc2.index[idx])
                df_equivol_per_depth_ecc2.loc[idx, 'rsq'] = pd.Series([rsqxlayer]*len(idx), index=df_equivol_per_depth_ecc2.index[idx])
            
            depth_id += 1
        
# save df
if layer_space == 'func':
    f = open(opj(prfpy_dir,'..','df_prf_param_func'), 'wb')
else:
    f = open(opj(prfpy_dir,'..','df_prf_param'), 'wb')
pickle.dump([df_equivol_per_depth,df_equivol_per_depth_ecc2], f)
f.close()
