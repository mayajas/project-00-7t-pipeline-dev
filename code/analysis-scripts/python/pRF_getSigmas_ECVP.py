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
# prf model (either prfpy_Iso2DGaussianModel or prfpy_Iso2DGaussianModel_LinearDeveining)
layer_space        = 'anat'
prf_model          = 'prfpy_Iso2DGaussianModel'
#_LinearDeveining_4layers

#################################################################################################################
# Stimulus size (same for all subjects)
screen_size_cm     = 12.0
screen_distance_cm = 52.0
TR                 = 3.0
#design_matrix      = mat["stim"]
# import math 
# max_ecc = math.atan(screen_size_cm/screen_distance_cm)/2
# max_ecc_deg = math.degrees(max_ecc)



#################################################################################################################
# Define important variables
subject_list = ['sub-01','sub-02','sub-03','sub-04']
hem_list     = ['lh','rh']
roi_list     = ['V1','V2','V3']

rsq_thresh = 0.1

n_layers   = 6
layers     = np.arange(1,n_layers+1)
n_sub      = len(subject_list)
n_hem      = len(hem_list)
n_roi      = len(roi_list)


#################################################################################################################
## Define data frame to store pRF values (pRF size x depth)

depth = np.linspace(1,n_layers,n_layers).tolist()*n_hem*n_sub*n_roi

rois = [['V1']*n_layers*n_hem,['V2']*n_layers*n_hem, ['V3']*n_layers*n_hem]
rois = list(itertools.chain(*rois))*n_sub

hems = [['lh']*n_layers,['rh']*n_layers]
hems = list(itertools.chain(*hems))*n_sub*n_roi

subs = [['sub-01']*n_layers*n_hem*n_roi,['sub-02']*n_layers*n_hem*n_roi,['sub-03']*n_layers*n_hem*n_roi,['sub-04']*n_layers*n_hem*n_roi]
subs = list(itertools.chain(*subs))

pRF_size = np.empty((1,len(depth),)).tolist()
pRF_size = list(itertools.chain(*pRF_size))

eccs = np.empty((1,len(depth),)).tolist()
eccs = list(itertools.chain(*eccs))

rsq_threshs = np.empty((1,len(depth),)).tolist()
rsq_threshs = list(itertools.chain(*rsq_threshs))

# equivol x depth
df_equivol_per_depth = pd.DataFrame({
    'sub id' : subs,
    'roi' : rois,
    'hem' : hems,
    'depth' : depth,
    'pRF size' : pRF_size,
    'ecc' : eccs,
    'rsq' : rsq_threshs
})

# equidist x depth
df_equidist_per_depth = pd.DataFrame({
    'sub id' : subs,
    'roi' : rois,
    'hem' : hems,
    'depth' : depth,
    'pRF size' : pRF_size,
    'ecc' : eccs,
    'rsq' : rsq_threshs
})

#################################################################################################################
## Define data frame to store pRF values (pRF size x depth) @ 2deg iso-ecc

depth = np.linspace(1,n_layers,n_layers).tolist()*n_hem*n_sub*n_roi

rois = [['V1']*n_layers*n_hem,['V2']*n_layers*n_hem, ['V3']*n_layers*n_hem]
rois = list(itertools.chain(*rois))*n_sub

hems = [['lh']*n_layers,['rh']*n_layers]
hems = list(itertools.chain(*hems))*n_sub*n_roi

subs = [['sub-01']*n_layers*n_hem*n_roi,['sub-02']*n_layers*n_hem*n_roi,['sub-03']*n_layers*n_hem*n_roi,['sub-04']*n_layers*n_hem*n_roi]
subs = list(itertools.chain(*subs))

pRF_size = np.empty((1,len(depth),)).tolist()
pRF_size = list(itertools.chain(*pRF_size))

eccs = np.empty((1,len(depth),)).tolist()
eccs = list(itertools.chain(*eccs))

rsq_threshs = np.empty((1,len(depth),)).tolist()
rsq_threshs = list(itertools.chain(*rsq_threshs))

cortdist = np.empty((1,len(depth),)).tolist()
cortdist = list(itertools.chain(*cortdist))

# equivol x depth
df_equivol_per_depth_ecc2 = pd.DataFrame({
    'sub id' : subs,
    'roi' : rois,
    'hem' : hems,
    'depth' : depth,
    'pRF size' : pRF_size,
    'ecc' : eccs,
    'rsq' : rsq_threshs,
    'cortdist' : cortdist
})

# equidist x depth
df_equidist_per_depth_ecc2 = pd.DataFrame({
    'sub id' : subs,
    'roi' : rois,
    'hem' : hems,
    'depth' : depth,
    'pRF size' : pRF_size,
    'ecc' : eccs,
    'rsq' : rsq_threshs,
    'cortdist' : cortdist
})

#################################################################################################################
## Define data frame to store pRF values (pRF size x cortical distance)

rois = [['V1'],['V2'], ['V3']]
rois = list(itertools.chain(*rois))*n_sub*n_hem

hems = [['lh']*n_roi,['rh']*n_roi]
hems = list(itertools.chain(*hems))*n_sub

subs = [['sub-01']*n_hem*n_roi,['sub-02']*n_hem*n_roi,['sub-03']*n_hem*n_roi,['sub-04']*n_hem*n_roi]
subs = list(itertools.chain(*subs))

pRF_size = np.empty((1,len(hems),)).tolist()
pRF_size = list(itertools.chain(*pRF_size))

eccs = np.empty((1,len(hems),)).tolist()
eccs = list(itertools.chain(*eccs))

rsq_threshs = np.empty((1,len(hems),)).tolist()
rsq_threshs = list(itertools.chain(*rsq_threshs))

cort_dist = np.empty((1,len(hems),)).tolist()
cort_dist = list(itertools.chain(*cort_dist))

# equivol x depth
df_pRF_per_cortdist = pd.DataFrame({
    'sub id' : subs,
    'roi' : rois,
    'hem' : hems,
    'pRF size' : pRF_size,
    'ecc' : eccs,
    'rsq' : rsq_threshs,
    'cort dist': cort_dist
})

#################################################################################################################
## Define data frame to store pRF values (pRF size x cortical distance) at manually-defined iso-eccentricity band

rois = [['V1'],['V2'], ['V3']]
rois = list(itertools.chain(*rois))*n_sub*n_hem

hems = [['lh']*n_roi,['rh']*n_roi]
hems = list(itertools.chain(*hems))*n_sub

subs = [['sub-01']*n_hem*n_roi,['sub-02']*n_hem*n_roi,['sub-03']*n_hem*n_roi,['sub-04']*n_hem*n_roi]
subs = list(itertools.chain(*subs))

pRF_size = np.empty((1,len(hems),)).tolist()
pRF_size = list(itertools.chain(*pRF_size))

eccs = np.empty((1,len(hems),)).tolist()
eccs = list(itertools.chain(*eccs))

rsq_threshs = np.empty((1,len(hems),)).tolist()
rsq_threshs = list(itertools.chain(*rsq_threshs))

cort_dist = np.empty((1,len(hems),)).tolist()
cort_dist = list(itertools.chain(*cort_dist))

# equivol x depth
df_pRF_per_cortdist_ecc2 = pd.DataFrame({
    'sub id' : subs,
    'roi' : rois,
    'hem' : hems,
    'pRF size' : pRF_size,
    'ecc' : eccs,
    'rsq' : rsq_threshs,
    'cort dist': cort_dist
})

#################################################################################################################
## Loop over subjects and hemispheres
for sub_id in range(0,len(subject_list)):
    #print('#################################################################################################################')
    #print(subject_list[sub_id])
    for hem_id in range(0,len(hem_list)):
        #print('')
        #print('Hemisphere: '+hem_list[hem_id])
        ###########################################################################################
        ## Define input data directories and filenames
        # data directories
        prfpy_dir    = '/home/mayajas/scratch/project-00-7t-pipeline-dev/output/'+prf_model+'/'+subject_list[sub_id]
        proj_dir     = '/home/mayajas/scratch/project-00-7t-pipeline-dev/'
        data_dir     = opj(proj_dir,'output','func','sliceTimeCorr',
                        '_subject_id_'+subject_list[sub_id])
        ROI_dir      = opj(prfpy_dir,'ROIs')
        if layer_space == 'func':
            lay_dir      = opj(prfpy_dir,'layerification_func')
        else:
            lay_dir      = opj(prfpy_dir,'layerification')
        lay_dir_func     = opj(prfpy_dir,'layerification_func')

        # image files
        meanFunc_fn  = opj(prfpy_dir,'meanFunc.nii')
        UNI_func_fn  = opj(prfpy_dir,'UNI_funcSpace.nii')

        GM_fn        = opj(prfpy_dir,
                        hem_list[hem_id]+'_GM_funcSpace.nii')
        occ_fn       = opj(prfpy_dir,
                        'funcSpaceOccipitalMask.nii')
        bar1_fn      = opj(data_dir,'_sess_id_task-bar_run-01_sess_nr_0_sess_nvol_124','atask-bar_run-01_roi_warp4D.nii')
        bar2_fn      = opj(data_dir,'_sess_id_task-bar_run-02_sess_nr_1_sess_nvol_124','atask-bar_run-02_roi_warp4D.nii')


        # prfpy output files
        grid_fit_fn      = opj(prfpy_dir,hem_list[hem_id]+'_grid_fit.pckl')
        iterative_fit_fn = opj(prfpy_dir,hem_list[hem_id]+'_iterative_fit.pckl')
        pRF_param_fn     = opj(prfpy_dir,hem_list[hem_id]+'_pRF_params.pckl')
        equivol_param_fn = opj(prfpy_dir,hem_list[hem_id]+'_pRF_params_equivol.pckl')
        equidist_param_fn = opj(prfpy_dir,hem_list[hem_id]+'_pRF_params_equidist.pckl')

        # manual ROI delineations
        V1_fn        = opj(ROI_dir,'func_'+hem_list[hem_id]+'_V1.nii')
        V2_fn        = opj(ROI_dir,'func_'+hem_list[hem_id]+'_V2.nii')
        V3_fn        = opj(ROI_dir,'func_'+hem_list[hem_id]+'_V3.nii')
        ecc2_fn      = opj(ROI_dir,'func_'+hem_list[hem_id]+'_ecc2.nii')

        # layers
        lay_equidist_fn = opj(lay_dir,'func_ribbon_rim_layers_equidist.nii')
        lay_equivol_fn = opj(lay_dir,'func_ribbon_rim_layers_equivol.nii')

        # cortical distances (mm)
        cortdist_fn = opj(lay_dir_func,'cortical_distances.nii')

        ###########################################################################################
        ## Load image files
        occ = nib.load(occ_fn)
        GM = nib.load(GM_fn)
        meanFunc = nib.load(meanFunc_fn)
        UNI = nib.load(UNI_func_fn)

        V1=nib.load(V1_fn)
        V2=nib.load(V2_fn)
        V3=nib.load(V3_fn)

        ecc2=nib.load(ecc2_fn)

        lay_equidist=nib.load(lay_equidist_fn)
        lay_equivol=nib.load(lay_equivol_fn)
        cort_distances=nib.load(cortdist_fn)
        

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

        masked_ecc2 = masker.fit_transform(ecc2)
        masked_ecc2 = np.squeeze(masked_ecc2)

        masked_lay_equidist = masker.fit_transform(lay_equidist)
        masked_lay_equidist = np.squeeze(masked_lay_equidist)

        masked_lay_equivol = masker.fit_transform(lay_equivol)
        masked_lay_equivol = np.squeeze(masked_lay_equivol)

        masked_cort_distances = masker.fit_transform(cort_distances)
        masked_cort_distances = np.squeeze(masked_cort_distances)

        ###########################################################################################
        ## Load pRF results
        f = open(pRF_param_fn,'rb')
        x, y, sigma, total_rsq, polar, ecc = pickle.load(f)
        f.close()

        ###########################################################################################
        ## Extract sigma per depth: equivolumetric
        depth_id=0
        for layer in layers:
            for roi_id in range(0,len(roi_list)):
                if roi_id == 0:
                    masked_ROI = masked_V1
                elif roi_id == 1:
                    masked_ROI = masked_V2
                elif roi_id == 2:
                    masked_ROI = masked_V3
                sigmaxlayer = sigma[(masked_ROI==1) & (masked_lay_equivol == layer)]
                eccxlayer = ecc[(masked_ROI==1) & (masked_lay_equivol == layer)]
                rsqxlayer = total_rsq[(masked_ROI==1) & (masked_lay_equivol == layer)]
                
                idx=df_equivol_per_depth.loc[(df_equivol_per_depth['sub id'] == subject_list[sub_id]) & 
                                        (df_equivol_per_depth['hem'] == hem_list[hem_id]) & 
                                        (df_equivol_per_depth['depth'] == depth[depth_id]) &
                                        (df_equivol_per_depth['roi'] == roi_list[roi_id])].index.tolist()
                df_equivol_per_depth.loc[idx, 'pRF size'] = pd.Series([sigmaxlayer]*len(idx), index=df_equivol_per_depth.index[idx])
                df_equivol_per_depth.loc[idx, 'ecc'] = pd.Series([eccxlayer]*len(idx), index=df_equivol_per_depth.index[idx])
                df_equivol_per_depth.loc[idx, 'rsq'] = pd.Series([rsqxlayer]*len(idx), index=df_equivol_per_depth.index[idx])
                
            depth_id += 1
        ###########################################################################################
        ## Extract sigma per depth: equidistant
        depth_id=0
        for layer in layers:
            for roi_id in range(0,len(roi_list)):
                if roi_id == 0:
                    masked_ROI = masked_V1
                elif roi_id == 1:
                    masked_ROI = masked_V2
                elif roi_id == 2:
                    masked_ROI = masked_V3
                sigmaxlayer = sigma[(masked_ROI==1) & (masked_lay_equidist == layer)]
                eccxlayer = ecc[(masked_ROI==1) & (masked_lay_equidist == layer)]
                rsqxlayer = total_rsq[(masked_ROI==1) & (masked_lay_equidist == layer)]

                idx=df_equidist_per_depth.loc[(df_equidist_per_depth['sub id'] == subject_list[sub_id]) & 
                                        (df_equidist_per_depth['hem'] == hem_list[hem_id]) & 
                                        (df_equidist_per_depth['depth'] == depth[depth_id]) &
                                        (df_equidist_per_depth['roi'] == roi_list[roi_id])].index.tolist()
                df_equidist_per_depth.loc[idx, 'pRF size'] = pd.Series([sigmaxlayer]*len(idx), index=df_equidist_per_depth.index[idx])
                df_equidist_per_depth.loc[idx, 'ecc'] = pd.Series([eccxlayer]*len(idx), index=df_equidist_per_depth.index[idx])
                df_equidist_per_depth.loc[idx, 'rsq'] = pd.Series([rsqxlayer]*len(idx), index=df_equidist_per_depth.index[idx])
            depth_id += 1
        
        ###########################################################################################
        ## Extract sigma per depth: equivolumetric (manually-defined iso-ecc)
        depth_id=0
        for layer in layers:
            for roi_id in range(0,len(roi_list)):
                if roi_id == 0:
                    masked_ROI = masked_V1
                elif roi_id == 1:
                    masked_ROI = masked_V2
                elif roi_id == 2:
                    masked_ROI = masked_V3
                sigmaxlayer = sigma[(masked_ROI==1) & (masked_lay_equivol == layer) & (masked_ecc2==1)]
                eccxlayer = ecc[(masked_ROI==1) & (masked_lay_equivol == layer) & (masked_ecc2==1)]
                rsqxlayer = total_rsq[(masked_ROI==1) & (masked_lay_equivol == layer) & (masked_ecc2==1)]
                cortdistxlayer = masked_cort_distances[(masked_ROI==1) & (masked_lay_equivol == layer) & (masked_ecc2==1)]

                
                idx=df_equivol_per_depth_ecc2.loc[(df_equivol_per_depth_ecc2['sub id'] == subject_list[sub_id]) & 
                                        (df_equivol_per_depth_ecc2['hem'] == hem_list[hem_id]) & 
                                        (df_equivol_per_depth_ecc2['depth'] == depth[depth_id]) &
                                        (df_equivol_per_depth_ecc2['roi'] == roi_list[roi_id])].index.tolist()
                df_equivol_per_depth_ecc2.loc[idx, 'pRF size'] = pd.Series([sigmaxlayer]*len(idx), index=df_equivol_per_depth_ecc2.index[idx])
                df_equivol_per_depth_ecc2.loc[idx, 'ecc'] = pd.Series([eccxlayer]*len(idx), index=df_equivol_per_depth_ecc2.index[idx])
                df_equivol_per_depth_ecc2.loc[idx, 'rsq'] = pd.Series([rsqxlayer]*len(idx), index=df_equivol_per_depth_ecc2.index[idx])
                df_equivol_per_depth_ecc2.loc[idx, 'cortdist'] = pd.Series([cortdistxlayer]*len(idx), index=df_equivol_per_depth_ecc2.index[idx])
            
            depth_id += 1
        ###########################################################################################
        ## Extract sigma per depth: equidistant (manually-defined iso-ecc)
        depth_id=0
        for layer in layers:
            for roi_id in range(0,len(roi_list)):
                if roi_id == 0:
                    masked_ROI = masked_V1
                elif roi_id == 1:
                    masked_ROI = masked_V2
                elif roi_id == 2:
                    masked_ROI = masked_V3
                sigmaxlayer = sigma[(masked_ROI==1) & (masked_lay_equidist == layer) & (masked_ecc2==1)]
                eccxlayer = ecc[(masked_ROI==1) & (masked_lay_equidist == layer) & (masked_ecc2==1)]
                rsqxlayer = total_rsq[(masked_ROI==1) & (masked_lay_equidist == layer) & (masked_ecc2==1)]
                cortdistxlayer = masked_cort_distances[(masked_ROI==1) & (masked_lay_equidist == layer) & (masked_ecc2==1)]

                idx=df_equidist_per_depth_ecc2.loc[(df_equidist_per_depth_ecc2['sub id'] == subject_list[sub_id]) & 
                                        (df_equidist_per_depth_ecc2['hem'] == hem_list[hem_id]) & 
                                        (df_equidist_per_depth_ecc2['depth'] == depth[depth_id])  &
                                        (df_equidist_per_depth_ecc2['roi'] == roi_list[roi_id])].index.tolist()
                df_equidist_per_depth_ecc2.loc[idx, 'pRF size'] = pd.Series([sigmaxlayer]*len(idx), index=df_equidist_per_depth_ecc2.index[idx])
                df_equidist_per_depth_ecc2.loc[idx, 'ecc'] = pd.Series([eccxlayer]*len(idx), index=df_equidist_per_depth_ecc2.index[idx])
                df_equidist_per_depth_ecc2.loc[idx, 'rsq'] = pd.Series([rsqxlayer]*len(idx), index=df_equidist_per_depth_ecc2.index[idx])
                df_equidist_per_depth_ecc2.loc[idx, 'cortdist'] = pd.Series([cortdistxlayer]*len(idx), index=df_equidist_per_depth_ecc2.index[idx])
            depth_id += 1
        
        ###########################################################################################
        ## Extract sigma per cortical distance
        for roi_id in range(0,len(roi_list)):
            if roi_id == 0:
                masked_ROI = masked_V1
            elif roi_id == 1:
                masked_ROI = masked_V2
            elif roi_id == 2:
                masked_ROI = masked_V3
            sigmaxcortdist = sigma[(masked_ROI==1)]
            eccxcortdist = ecc[(masked_ROI==1)]
            rsqxcortdist = total_rsq[(masked_ROI==1)]
            cortdist     = masked_cort_distances[(masked_ROI==1)]

            idx=df_pRF_per_cortdist.loc[(df_pRF_per_cortdist['sub id'] == subject_list[sub_id]) & 
                                        (df_pRF_per_cortdist['hem'] == hem_list[hem_id]) &
                                        (df_pRF_per_cortdist['roi'] == roi_list[roi_id])].index.tolist()
            df_pRF_per_cortdist.loc[idx, 'pRF size'] = pd.Series([sigmaxcortdist]*len(idx), index=df_pRF_per_cortdist.index[idx])
            df_pRF_per_cortdist.loc[idx, 'ecc'] = pd.Series([eccxcortdist]*len(idx), index=df_pRF_per_cortdist.index[idx])
            df_pRF_per_cortdist.loc[idx, 'rsq'] = pd.Series([rsqxcortdist]*len(idx), index=df_pRF_per_cortdist.index[idx])
            df_pRF_per_cortdist.loc[idx, 'cort dist'] = pd.Series([cortdist]*len(idx), index=df_pRF_per_cortdist.index[idx])
        

        ###########################################################################################
        ## Extract sigma per cortical distance @ manually-defined iso-eccentricity band (1-3 deg)
        for roi_id in range(0,len(roi_list)):
            if roi_id == 0:
                masked_ROI = masked_V1
            elif roi_id == 1:
                masked_ROI = masked_V2
            elif roi_id == 2:
                masked_ROI = masked_V3
            sigmaxcortdist = sigma[(masked_ROI==1) & (masked_ecc2==1)]
            eccxcortdist = ecc[(masked_ROI==1) & (masked_ecc2==1)]
            rsqxcortdist = total_rsq[(masked_ROI==1) & (masked_ecc2==1)]
            cortdist     = masked_cort_distances[(masked_ROI==1) & (masked_ecc2==1)]

            idx=df_pRF_per_cortdist_ecc2.loc[(df_pRF_per_cortdist_ecc2['sub id'] == subject_list[sub_id]) & 
                                        (df_pRF_per_cortdist_ecc2['hem'] == hem_list[hem_id]) &
                                        (df_pRF_per_cortdist_ecc2['roi'] == roi_list[roi_id])].index.tolist()
            df_pRF_per_cortdist_ecc2.loc[idx, 'pRF size'] = pd.Series([sigmaxcortdist]*len(idx), index=df_pRF_per_cortdist_ecc2.index[idx])
            df_pRF_per_cortdist_ecc2.loc[idx, 'ecc'] = pd.Series([eccxcortdist]*len(idx), index=df_pRF_per_cortdist_ecc2.index[idx])
            df_pRF_per_cortdist_ecc2.loc[idx, 'rsq'] = pd.Series([rsqxcortdist]*len(idx), index=df_pRF_per_cortdist_ecc2.index[idx])
            df_pRF_per_cortdist_ecc2.loc[idx, 'cort dist'] = pd.Series([cortdist]*len(idx), index=df_pRF_per_cortdist_ecc2.index[idx])


# save df
if layer_space == 'func':
    f = open(opj(prfpy_dir,'..','df_layers_allROI_func'), 'wb')
else:
    f = open(opj(prfpy_dir,'..','df_layers_allROI'), 'wb')
pickle.dump([df_equivol_per_depth,df_equivol_per_depth_ecc2,df_equidist_per_depth,df_equidist_per_depth_ecc2,df_pRF_per_cortdist,df_pRF_per_cortdist_ecc2], f)
f.close()
