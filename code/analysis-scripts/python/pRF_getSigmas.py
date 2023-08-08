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
# prf model (either prfpy_fit_hrf_False_start_from_avg_False, 
# prfpy_fit_hrf_True_start_from_avg_False,
# prfpy_fit_hrf_False_start_from_avg_True or
# prfpy_fit_hrf_True_start_from_avg_True)
fit_hrf        = True  # boolean, optional
                              # Whether or not to fit two extra parameters for hrf derivative and
                              # dispersion. The default is False.
start_from_avg = True   # whether to use avg across depths as starting point for layer fits
prf_model      = 'prfpy_fit_hrf_'+str(fit_hrf)+'_start_from_avg_'+str(start_from_avg)

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

n_layers   = 7 # there are 6 layers; layer 7 corresponds to avg across layers
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
    'hrf1': hrf_1s,
    'hrf2': hrf_2s,
    'rsq' : rsq_threshs
})
df_equivol_per_depth.loc[df_equivol_per_depth.layer == 7.0,'layer'] = 'avg'


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
    'hrf1': hrf_1s,
    'hrf2': hrf_2s,
    'rsq' : rsq_threshs
})
df_equivol_per_depth_ecc2.loc[df_equivol_per_depth_ecc2.layer == 7.0,'layer'] = 'avg'

#################################################################################################################
## Loop over subjects and hemispheres
for sub_id in range(0,len(subject_list)):
    for hem_id in range(0,len(hem_list)):
        ## Define input data directories and filenames
        # data directories
        prfpy_dir    = '/home/mayajas/scratch/project-00-7t-pipeline-dev/output/'+prf_model+'/'+subject_list[sub_id]
        proj_dir     = '/home/mayajas/scratch/project-00-7t-pipeline-dev/'

        # prfpy output files
        pRF_param_avg_fn       = opj(prfpy_dir,hem_list[hem_id]+'_pRF_params_avg.pckl')
        pRF_param_per_depth_fn = opj(prfpy_dir,hem_list[hem_id]+'_pRF_params_per_depth.pckl')

        # occipital mask (for unmasking prf outputs)
        occ_mask_fn          = opj(prfpy_dir,hem_list[hem_id]+'_occ_mask.pckl')
        ecc2_fn      = opj(prfpy_dir,hem_list[hem_id]+'.ecc2.label')

        ###########################################################################################
        ## Load occipital mask
        ecc2=nib.freesurfer.io.read_label(ecc2_fn)

        ## Load occipital mask
        f = open(occ_mask_fn,'rb')
        occ_mask,n_vtx = pickle.load(f)
        f.close()

        ###########################################################################################
        ## Load & unmask avg pRF results
        f = open(pRF_param_avg_fn,'rb')
        if fit_hrf:
            x, y, sigma, total_rsq, polar, ecc, hrf_1, hrf_2 = pickle.load(f)
        else:
            x, y, sigma, total_rsq, polar, ecc = pickle.load(f)
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
        # remove bad fits 
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
        unmask_x[unmask_sigma>max_ecc_deg]     = np.nan
        unmask_y[unmask_sigma>max_ecc_deg]     = np.nan
        unmask_polar[unmask_sigma>max_ecc_deg] = np.nan
        unmask_ecc[unmask_sigma>max_ecc_deg]   = np.nan
        unmask_sigma[unmask_sigma>max_ecc_deg] = np.nan
        if fit_hrf:
            unmask_hrf_1[unmask_sigma>max_ecc_deg] = np.nan
            unmask_hrf_2[unmask_sigma>max_ecc_deg] = np.nan

        ###########################################################################################
        ## Load & unmask per-depth pRF results
        n_surfs = int(max(df_equivol_per_depth[df_equivol_per_depth["layer"].apply(lambda x: isinstance(x, float))].layer.values))

        f = open(pRF_param_per_depth_fn,'rb')
        if fit_hrf:
            x_per_depth, y_per_depth, sigma_per_depth, total_rsq_per_depth, polar_per_depth, ecc_per_depth, hrf_1_per_depth,hrf_2_per_depth = pickle.load(f)
        else:
            x_per_depth, y_per_depth, sigma_per_depth, total_rsq_per_depth, polar_per_depth, ecc_per_depth = pickle.load(f)
        f.close()

        # unmask per depth pRF parameters
        unmask_x_per_depth     = [np.zeros(n_vtx) for depth in range(0,n_surfs)]
        unmask_y_per_depth     = [np.zeros(n_vtx) for depth in range(0,n_surfs)]
        unmask_sigma_per_depth = [np.zeros(n_vtx) for depth in range(0,n_surfs)]
        unmask_rsq_per_depth   = [np.zeros(n_vtx) for depth in range(0,n_surfs)]
        unmask_polar_per_depth = [np.zeros(n_vtx) for depth in range(0,n_surfs)]
        unmask_ecc_per_depth   = [np.zeros(n_vtx) for depth in range(0,n_surfs)]
        if fit_hrf:
            unmask_hrf_1_per_depth    = [np.zeros(n_vtx) for depth in range(0,n_surfs)]
            unmask_hrf_2_per_depth    = [np.zeros(n_vtx) for depth in range(0,n_surfs)]

        for depth in range(0,n_surfs):
            unmask_x_per_depth[depth][occ_mask]     = x_per_depth[depth]
            unmask_y_per_depth[depth][occ_mask]     = y_per_depth[depth]
            unmask_sigma_per_depth[depth][occ_mask] = sigma_per_depth[depth]
            unmask_rsq_per_depth[depth][occ_mask]   = total_rsq_per_depth[depth]
            unmask_polar_per_depth[depth][occ_mask] = polar_per_depth[depth]
            unmask_ecc_per_depth[depth][occ_mask]   = ecc_per_depth[depth]
            if fit_hrf:
                unmask_hrf_1_per_depth[depth][occ_mask]    = hrf_1_per_depth[depth]
                unmask_hrf_2_per_depth[depth][occ_mask]    = hrf_2_per_depth[depth]
            
            # remove bad fits - per depth
            unmask_x_per_depth[depth][unmask_rsq_per_depth[depth]<rsq_thresh]         = np.nan
            unmask_y_per_depth[depth][unmask_rsq_per_depth[depth]<rsq_thresh]         = np.nan
            unmask_sigma_per_depth[depth][unmask_rsq_per_depth[depth]<rsq_thresh]     = np.nan
            unmask_polar_per_depth[depth][unmask_rsq_per_depth[depth]<rsq_thresh]     = np.nan
            unmask_ecc_per_depth[depth][unmask_rsq_per_depth[depth]<rsq_thresh]       = np.nan
            if fit_hrf:
                unmask_hrf_1_per_depth[depth][unmask_rsq_per_depth[depth]<rsq_thresh] = np.nan
                unmask_hrf_2_per_depth[depth][unmask_rsq_per_depth[depth]<rsq_thresh] = np.nan

            
            # remove vertices where eccentricity is larger than max stimulus ecc
            unmask_x_per_depth[depth][unmask_ecc_per_depth[depth]>max_ecc_deg]         = np.nan
            unmask_y_per_depth[depth][unmask_ecc_per_depth[depth]>max_ecc_deg]         = np.nan
            unmask_polar_per_depth[depth][unmask_ecc_per_depth[depth]>max_ecc_deg]     = np.nan
            unmask_sigma_per_depth[depth][unmask_ecc_per_depth[depth]>max_ecc_deg]     = np.nan
            if fit_hrf:
                unmask_hrf_1_per_depth[depth][unmask_ecc_per_depth[depth]>max_ecc_deg] = np.nan
                unmask_hrf_2_per_depth[depth][unmask_ecc_per_depth[depth]>max_ecc_deg] = np.nan
            unmask_ecc_per_depth[depth][unmask_ecc_per_depth[depth]>max_ecc_deg]       = np.nan
            
            # remove vertices where pRF size is negative
            unmask_x_per_depth[depth][unmask_sigma_per_depth[depth]<0]         = np.nan
            unmask_y_per_depth[depth][unmask_sigma_per_depth[depth]<0]         = np.nan
            unmask_polar_per_depth[depth][unmask_sigma_per_depth[depth]<0]     = np.nan
            unmask_ecc_per_depth[depth][unmask_sigma_per_depth[depth]<0]       = np.nan
            if fit_hrf:
                unmask_hrf_1_per_depth[depth][unmask_sigma_per_depth[depth]<0] = np.nan
                unmask_hrf_2_per_depth[depth][unmask_sigma_per_depth[depth]<0] = np.nan
            unmask_sigma_per_depth[depth][unmask_sigma_per_depth[depth]<0]     = np.nan

            # set max pRF size to max stimulus eccentricity
            unmask_x_per_depth[depth][unmask_sigma_per_depth[depth]>max_ecc_deg]         = np.nan
            unmask_y_per_depth[depth][unmask_sigma_per_depth[depth]>max_ecc_deg]         = np.nan
            unmask_polar_per_depth[depth][unmask_sigma_per_depth[depth]>max_ecc_deg]     = np.nan
            unmask_ecc_per_depth[depth][unmask_sigma_per_depth[depth]>max_ecc_deg]       = np.nan
            if fit_hrf:
                unmask_hrf_1_per_depth[depth][unmask_sigma_per_depth[depth]>max_ecc_deg] = np.nan
                unmask_hrf_2_per_depth[depth][unmask_sigma_per_depth[depth]>max_ecc_deg] = np.nan
            unmask_sigma_per_depth[depth][unmask_sigma_per_depth[depth]>max_ecc_deg]     = np.nan

        ###########################################################################################
        ## Extract sigma across eccentricities
        for roi_id in range(0,len(roi_list)):
            # manual ROI delineation
            ROI_fn       = opj(prfpy_dir,hem_list[hem_id]+'.'+roi_list[roi_id]+'.label')
            ROI=nib.freesurfer.io.read_label(ROI_fn)
            
            # across depth
            for depth in range(0,n_surfs):      
                sigmaxlayer = unmask_sigma_per_depth[depth][ROI]
                eccxlayer   = unmask_ecc_per_depth[depth][ROI]
                rsqxlayer   = unmask_rsq_per_depth[depth][ROI]
                if fit_hrf:
                    hrf1xlayer  = unmask_hrf_1_per_depth[depth][ROI]
                    hrf2xlayer  = unmask_hrf_2_per_depth[depth][ROI]

                idx=df_equivol_per_depth.loc[(df_equivol_per_depth['sub id'] == subject_list[sub_id]) & 
                                        (df_equivol_per_depth['hem'] == hem_list[hem_id]) & 
                                        (df_equivol_per_depth['layer'] == layers[depth]) &
                                        (df_equivol_per_depth['roi'] == roi_list[roi_id])].index.tolist()
                df_equivol_per_depth.loc[idx, 'pRF size'] = pd.Series([sigmaxlayer]*len(idx), index=df_equivol_per_depth.index[idx])
                df_equivol_per_depth.loc[idx, 'ecc'] = pd.Series([eccxlayer]*len(idx), index=df_equivol_per_depth.index[idx])
                if fit_hrf:
                    df_equivol_per_depth.loc[idx, 'hrf1'] = pd.Series([hrf1xlayer]*len(idx), index=df_equivol_per_depth.index[idx])
                    df_equivol_per_depth.loc[idx, 'hrf2'] = pd.Series([hrf2xlayer]*len(idx), index=df_equivol_per_depth.index[idx]) 
                df_equivol_per_depth.loc[idx, 'rsq'] = pd.Series([rsqxlayer]*len(idx), index=df_equivol_per_depth.index[idx]) 
                
            # avg
            sigmaxlayer = unmask_sigma[ROI]
            eccxlayer   = unmask_ecc[ROI]
            rsqxlayer   = unmask_rsq[ROI]
            if fit_hrf:
                    hrf1xlayer  = unmask_hrf_1[ROI]
                    hrf2xlayer  = unmask_hrf_2[ROI]

            idx=df_equivol_per_depth.loc[(df_equivol_per_depth['sub id'] == subject_list[sub_id]) & 
                                    (df_equivol_per_depth['hem'] == hem_list[hem_id]) & 
                                    (df_equivol_per_depth['layer'] == 'avg') &
                                    (df_equivol_per_depth['roi'] == roi_list[roi_id])].index.tolist()
            df_equivol_per_depth.loc[idx, 'pRF size'] = pd.Series([sigmaxlayer]*len(idx), index=df_equivol_per_depth.index[idx])
            df_equivol_per_depth.loc[idx, 'ecc'] = pd.Series([eccxlayer]*len(idx), index=df_equivol_per_depth.index[idx])
            if fit_hrf:
                    df_equivol_per_depth.loc[idx, 'hrf1'] = pd.Series([hrf1xlayer]*len(idx), index=df_equivol_per_depth.index[idx])
                    df_equivol_per_depth.loc[idx, 'hrf2'] = pd.Series([hrf2xlayer]*len(idx), index=df_equivol_per_depth.index[idx]) 
            df_equivol_per_depth.loc[idx, 'rsq'] = pd.Series([rsqxlayer]*len(idx), index=df_equivol_per_depth.index[idx])


        ###########################################################################################
        ## Extract sigma at manually-defined iso-ecc
        for roi_id in range(0,len(roi_list)):
            # manual ROI delineation
            ROI_fn       = opj(prfpy_dir,hem_list[hem_id]+'.'+roi_list[roi_id]+'.label')
            ROI=nib.freesurfer.io.read_label(ROI_fn)

            # across depth
            for depth in range(0,n_surfs):      
                sigmaxlayer = unmask_sigma_per_depth[depth][np.intersect1d(ROI,ecc2)]
                eccxlayer   = unmask_ecc_per_depth[depth][np.intersect1d(ROI,ecc2)]
                rsqxlayer   = unmask_rsq_per_depth[depth][np.intersect1d(ROI,ecc2)]
                if fit_hrf:
                    hrf1xlayer  = unmask_hrf_1_per_depth[depth][np.intersect1d(ROI,ecc2)]
                    hrf2xlayer  = unmask_hrf_2_per_depth[depth][np.intersect1d(ROI,ecc2)]

                idx=df_equivol_per_depth_ecc2.loc[(df_equivol_per_depth_ecc2['sub id'] == subject_list[sub_id]) & 
                                        (df_equivol_per_depth_ecc2['hem'] == hem_list[hem_id]) & 
                                        (df_equivol_per_depth_ecc2['layer'] == layers[depth]) &
                                        (df_equivol_per_depth_ecc2['roi'] == roi_list[roi_id])].index.tolist()
                df_equivol_per_depth_ecc2.loc[idx, 'pRF size'] = pd.Series([sigmaxlayer]*len(idx), index=df_equivol_per_depth_ecc2.index[idx])
                df_equivol_per_depth_ecc2.loc[idx, 'ecc'] = pd.Series([eccxlayer]*len(idx), index=df_equivol_per_depth_ecc2.index[idx])
                if fit_hrf:
                    df_equivol_per_depth_ecc2.loc[idx, 'hrf1'] = pd.Series([hrf1xlayer]*len(idx), index=df_equivol_per_depth_ecc2.index[idx])
                    df_equivol_per_depth_ecc2.loc[idx, 'hrf2'] = pd.Series([hrf2xlayer]*len(idx), index=df_equivol_per_depth_ecc2.index[idx])
                df_equivol_per_depth_ecc2.loc[idx, 'rsq'] = pd.Series([rsqxlayer]*len(idx), index=df_equivol_per_depth_ecc2.index[idx])
                
            # avg
            sigmaxlayer = unmask_sigma[np.intersect1d(ROI,ecc2)]
            eccxlayer   = unmask_ecc[np.intersect1d(ROI,ecc2)]
            rsqxlayer   = unmask_rsq[np.intersect1d(ROI,ecc2)]
            if fit_hrf:
                hrf1xlayer  = unmask_hrf_1[np.intersect1d(ROI,ecc2)]
                hrf2xlayer  = unmask_hrf_2[np.intersect1d(ROI,ecc2)]

            idx=df_equivol_per_depth_ecc2.loc[(df_equivol_per_depth_ecc2['sub id'] == subject_list[sub_id]) & 
                                    (df_equivol_per_depth_ecc2['hem'] == hem_list[hem_id]) & 
                                    (df_equivol_per_depth_ecc2['layer'] == 'avg') &
                                    (df_equivol_per_depth_ecc2['roi'] == roi_list[roi_id])].index.tolist()
            df_equivol_per_depth_ecc2.loc[idx, 'pRF size'] = pd.Series([sigmaxlayer]*len(idx), index=df_equivol_per_depth_ecc2.index[idx])
            df_equivol_per_depth_ecc2.loc[idx, 'ecc'] = pd.Series([eccxlayer]*len(idx), index=df_equivol_per_depth_ecc2.index[idx])
            if fit_hrf:
                df_equivol_per_depth_ecc2.loc[idx, 'hrf1'] = pd.Series([hrf1xlayer]*len(idx), index=df_equivol_per_depth_ecc2.index[idx])
                df_equivol_per_depth_ecc2.loc[idx, 'hrf2'] = pd.Series([hrf2xlayer]*len(idx), index=df_equivol_per_depth_ecc2.index[idx])
            df_equivol_per_depth_ecc2.loc[idx, 'rsq'] = pd.Series([rsqxlayer]*len(idx), index=df_equivol_per_depth_ecc2.index[idx])


# save df
f = open(opj(prfpy_dir,'..','df_prf_param'), 'wb')
pickle.dump([df_equivol_per_depth,df_equivol_per_depth_ecc2], f)
f.close()
