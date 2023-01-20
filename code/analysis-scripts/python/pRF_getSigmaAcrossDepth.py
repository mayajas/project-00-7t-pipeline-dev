import numpy as np

import sys
sys.path.append("/home/mayajas/Documents/programs/prfpy-main")

# from prfpy.stimulus import PRFStimulus2D
# from prfpy.model import Iso2DGaussianModel, Norm_Iso2DGaussianModel
# from prfpy.fit import Iso2DGaussianFitter, Norm_Iso2DGaussianFitter

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


#################################################################################################################
# Stimulus size (same for all subjects)
screen_size_cm     = 12.0
screen_distance_cm = 52.0
TR                 = 3.0
#design_matrix      = mat["stim"]

import math 
max_ecc = math.atan(screen_size_cm/screen_distance_cm)/2
max_ecc_deg = math.degrees(max_ecc)



#################################################################################################################
# Define important variables
subject_list = ['sub-01','sub-02','sub-03','sub-04']
hem_list     = ['lh','rh']

rsq_thresh = 0.1

#################################################################################################################
## Define data frame to store pRF values (pRF size x ecc)

eccs = np.arange(2, 6, 1).tolist()*8

hems = [['lh']*4,['rh']*4]
hems = list(itertools.chain(*hems))*4

subs = [['sub-01']*8,['sub-02']*8,['sub-03']*8,['sub-04']*8]
subs = list(itertools.chain(*subs))

pRF_size = np.empty((1,len(eccs),)).tolist()
pRF_size = list(itertools.chain(*pRF_size))

# equivol x ecc
depth = np.empty((1,len(eccs),)).tolist()
depth = list(itertools.chain(*depth))

df_equivol_per_ecc = pd.DataFrame({
    'sub id' : subs,
    'hem' : hems,
    'ecc' : eccs,
    'pRF size' : pRF_size,
    'depth' : depth
})
depth = [np.arange(0.2, 1.0, 0.1).tolist()]
idx = np.arange(0,len(df_equivol_per_ecc.index)).tolist()
df_equivol_per_ecc.loc[idx, 'depth'] = pd.Series(depth*len(df_equivol_per_ecc.index), index=df_equivol_per_ecc.index[idx])

# equidist
depth = np.empty((1,len(eccs),)).tolist()
depth = list(itertools.chain(*depth))

df_equidist_per_ecc = pd.DataFrame({
    'sub id' : subs,
    'hem' : hems,
    'ecc' : eccs,
    'pRF size' : pRF_size,
    'depth' : depth
})
depth = [np.arange(0.2, 1.0, 0.1).tolist()]
idx = np.arange(0,len(df_equidist_per_ecc.index)).tolist()
df_equidist_per_ecc.loc[idx, 'depth'] = pd.Series(depth*len(df_equidist_per_ecc.index), index=df_equidist_per_ecc.index[idx])




#################################################################################################################
## Define data frame to store pRF values (pRF size x depth)

depth = np.arange(0.2, 1.0, 0.1).tolist()*8

hems = [['lh']*8,['rh']*8]
hems = list(itertools.chain(*hems))*4

subs = [['sub-01']*16,['sub-02']*16,['sub-03']*16,['sub-04']*16]
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
    'hem' : hems,
    'depth' : depth,
    'pRF size' : pRF_size,
    'ecc' : eccs,
    'rsq' : rsq_threshs
})

# equidist x depth
df_equidist_per_depth = pd.DataFrame({
    'sub id' : subs,
    'hem' : hems,
    'depth' : depth,
    'pRF size' : pRF_size,
    'ecc' : eccs,
    'rsq' : rsq_threshs
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
        prfpy_dir    = '/home/mayajas/scratch/project-00-7t-pipeline-dev/output/prfpy/'+subject_list[sub_id]
        proj_dir     = '/home/mayajas/scratch/project-00-7t-pipeline-dev/'
        data_dir     = opj(proj_dir,'output','func','sliceTimeCorr',
                        '_subject_id_'+subject_list[sub_id])
        ROI_dir      = opj(prfpy_dir,'ROIs')
        lay_dir      = opj(prfpy_dir,'layerification')

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

        # layers
        lay_equidist_fn = opj(lay_dir,'func_ribbon_rim_layers_equidist.nii')
        lay_equivol_fn = opj(lay_dir,'func_ribbon_rim_layers_equivol.nii')

        ###########################################################################################
        ## Load image files
        occ = nib.load(occ_fn)
        GM = nib.load(GM_fn)
        meanFunc = nib.load(meanFunc_fn)
        UNI = nib.load(UNI_func_fn)
        V1=nib.load(V1_fn)
        lay_equidist=nib.load(lay_equidist_fn)
        lay_equivol=nib.load(lay_equivol_fn)

        # define mask (GM+occipital) and process bar data and masker
        mask = image.math_img("np.logical_and(img1, img2)", img1=occ, img2=GM)
        masker = NiftiMasker(mask_img=mask)

        # mask input files
        masked_meanFunc = masker.fit_transform(meanFunc)
        masked_meanFunc = np.squeeze(masked_meanFunc)

        masked_V1 = masker.fit_transform(V1)
        masked_V1 = np.squeeze(masked_V1)

        masked_lay_equidist = masker.fit_transform(lay_equidist)
        masked_lay_equidist = np.squeeze(masked_lay_equidist)

        masked_lay_equivol = masker.fit_transform(lay_equivol)
        masked_lay_equivol = np.squeeze(masked_lay_equivol)

        ###########################################################################################
        ## Load pRF results
        f = open(pRF_param_fn,'rb')
        x, y, sigma, total_rsq, polar, ecc = pickle.load(f)
        f.close()

        ###########################################################################################
        ## Extract sigma per ecc: equivolumetric
        sigma_ecc_2 = []
        sigma_ecc_3 = []
        sigma_ecc_4 = []
        sigma_ecc_5 = []
        #print('')
        #print('Equivolumetric layering:')
        for layer in range(2,np.size(np.unique(masked_lay_equivol))):
            #print('')
            #print('Layer '+str(layer))
            ecc_idx = (ecc >= 1.5) & (ecc < 2.5)
            s_val=np.nanmean(sigma[(ecc_idx) & (masked_V1==1) & (total_rsq>rsq_thresh)  & (masked_lay_equivol == layer)])
            #print("Ecc 2: "+str(s_val))
            sigma_ecc_2.append(s_val)
            
            ecc_idx = (ecc >= 2.5) & (ecc < 3.5)
            s_val=np.nanmean(sigma[(ecc_idx) & (masked_V1==1) & (total_rsq>rsq_thresh)  & (masked_lay_equivol == layer)])
            #print("Ecc 3: "+str(s_val))
            sigma_ecc_3.append(s_val)
            
            ecc_idx = (ecc >= 3.5) & (ecc < 4.5)
            s_val=np.nanmean(sigma[(ecc_idx) & (masked_V1==1) & (total_rsq>rsq_thresh)  & (masked_lay_equivol == layer)])
            #print("Ecc 4: "+str(s_val))
            sigma_ecc_4.append(s_val)
            
            ecc_idx = (ecc >= 4.5) & (ecc < 5.5)
            s_val=np.nanmean(sigma[(ecc_idx) & (masked_V1==1) & (total_rsq>rsq_thresh)  & (masked_lay_equivol == layer)])
            #print("Ecc 5: "+str(s_val))
            sigma_ecc_5.append(s_val)
        
        # save
        f = open(equivol_param_fn, 'wb')
        pickle.dump([sigma_ecc_2,sigma_ecc_3,sigma_ecc_4,sigma_ecc_5], f)
        f.close()

        # write to df
        idx=df_equivol_per_ecc.loc[(df_equivol_per_ecc['sub id'] == subject_list[sub_id]) & 
                           (df_equivol_per_ecc['hem'] == hem_list[hem_id]) & 
                           (df_equivol_per_ecc['ecc'] == 2)].index.tolist()
        df_equivol_per_ecc.loc[idx, 'pRF size'] = pd.Series([sigma_ecc_2]*len(idx), index=df_equivol_per_ecc.index[idx])

        idx=df_equivol_per_ecc.loc[(df_equivol_per_ecc['sub id'] == subject_list[sub_id]) & 
                           (df_equivol_per_ecc['hem'] == hem_list[hem_id]) & 
                           (df_equivol_per_ecc['ecc'] == 3)].index.tolist()
        df_equivol_per_ecc.loc[idx, 'pRF size'] = pd.Series([sigma_ecc_3]*len(idx), index=df_equivol_per_ecc.index[idx])

        idx=df_equivol_per_ecc.loc[(df_equivol_per_ecc['sub id'] == subject_list[sub_id]) & 
                           (df_equivol_per_ecc['hem'] == hem_list[hem_id]) & 
                           (df_equivol_per_ecc['ecc'] == 4)].index.tolist()
        df_equivol_per_ecc.loc[idx, 'pRF size'] = pd.Series([sigma_ecc_4]*len(idx), index=df_equivol_per_ecc.index[idx])

        idx=df_equivol_per_ecc.loc[(df_equivol_per_ecc['sub id'] == subject_list[sub_id]) & 
                           (df_equivol_per_ecc['hem'] == hem_list[hem_id]) & 
                           (df_equivol_per_ecc['ecc'] == 5)].index.tolist()
        df_equivol_per_ecc.loc[idx, 'pRF size'] = pd.Series([sigma_ecc_5]*len(idx), index=df_equivol_per_ecc.index[idx])



        ###########################################################################################
        ## Extract sigma per ecc: equidistant
        sigma_ecc_2 = []
        sigma_ecc_3 = []
        sigma_ecc_4 = []
        sigma_ecc_5 = []
        #print('')
        #print('Equidistant layering:')
        for layer in range(2,np.size(np.unique(masked_lay_equidist))):
            #print('')
            #print('Layer '+str(layer))
            ecc_idx = (ecc >= 1.5) & (ecc < 2.5)
            s_val=np.nanmean(sigma[(ecc_idx) & (masked_V1==1) & (total_rsq>rsq_thresh)  & (masked_lay_equidist == layer)])
            #print("Ecc 2: "+str(s_val))
            sigma_ecc_2.append(s_val)
            
            ecc_idx = (ecc >= 2.5) & (ecc < 3.5)
            s_val=np.nanmean(sigma[(ecc_idx) & (masked_V1==1) & (total_rsq>rsq_thresh)  & (masked_lay_equidist == layer)])
            #print("Ecc 3: "+str(s_val))
            sigma_ecc_3.append(s_val)
            
            ecc_idx = (ecc >= 3.5) & (ecc < 4.5)
            s_val=np.nanmean(sigma[(ecc_idx) & (masked_V1==1) & (total_rsq>rsq_thresh)  & (masked_lay_equidist == layer)])
            #print("Ecc 4: "+str(s_val))
            sigma_ecc_4.append(s_val)
            
            ecc_idx = (ecc >= 4.5) & (ecc < 5.5)
            s_val=np.nanmean(sigma[(ecc_idx) & (masked_V1==1) & (total_rsq>rsq_thresh)  & (masked_lay_equidist == layer)])
            #print("Ecc 5: "+str(s_val))
            sigma_ecc_5.append(s_val)
        
        # save
        f = open(equidist_param_fn, 'wb')
        pickle.dump([sigma_ecc_2,sigma_ecc_3,sigma_ecc_4,sigma_ecc_5], f)
        f.close()

        # write to df
        idx=df_equidist_per_ecc.loc[(df_equidist_per_ecc['sub id'] == subject_list[sub_id]) & 
                           (df_equidist_per_ecc['hem'] == hem_list[hem_id]) & 
                           (df_equidist_per_ecc['ecc'] == 2)].index.tolist()
        df_equidist_per_ecc.loc[idx, 'pRF size'] = pd.Series([sigma_ecc_2]*len(idx), index=df_equidist_per_ecc.index[idx])

        idx=df_equidist_per_ecc.loc[(df_equidist_per_ecc['sub id'] == subject_list[sub_id]) & 
                           (df_equidist_per_ecc['hem'] == hem_list[hem_id]) & 
                           (df_equidist_per_ecc['ecc'] == 3)].index.tolist()
        df_equidist_per_ecc.loc[idx, 'pRF size'] = pd.Series([sigma_ecc_3]*len(idx), index=df_equidist_per_ecc.index[idx])

        idx=df_equidist_per_ecc.loc[(df_equidist_per_ecc['sub id'] == subject_list[sub_id]) & 
                           (df_equidist_per_ecc['hem'] == hem_list[hem_id]) & 
                           (df_equidist_per_ecc['ecc'] == 4)].index.tolist()
        df_equidist_per_ecc.loc[idx, 'pRF size'] = pd.Series([sigma_ecc_4]*len(idx), index=df_equidist_per_ecc.index[idx])

        idx=df_equidist_per_ecc.loc[(df_equidist_per_ecc['sub id'] == subject_list[sub_id]) & 
                           (df_equidist_per_ecc['hem'] == hem_list[hem_id]) & 
                           (df_equidist_per_ecc['ecc'] == 5)].index.tolist()
        df_equidist_per_ecc.loc[idx, 'pRF size'] = pd.Series([sigma_ecc_5]*len(idx), index=df_equidist_per_ecc.index[idx])

        ###########################################################################################
        ## Extract sigma per depth: equivolumetric
        depth_idx=0
        for layer in range(2,np.size(np.unique(masked_lay_equivol))):
            sigmaxlayer = sigma[(masked_V1==1) & (masked_lay_equivol == layer)]
            eccxlayer = ecc[(masked_V1==1) & (masked_lay_equivol == layer)]
            rsqxlayer = total_rsq[(masked_V1==1) & (masked_lay_equivol == layer)]
            
            idx=df_equivol_per_depth.loc[(df_equivol_per_depth['sub id'] == subject_list[sub_id]) & 
                                       (df_equivol_per_depth['hem'] == hem_list[hem_id]) & 
                                       (df_equivol_per_depth['depth'] == depth[depth_idx])].index.tolist()
            df_equivol_per_depth.loc[idx, 'pRF size'] = pd.Series([sigmaxlayer]*len(idx), index=df_equivol_per_depth.index[idx])
            df_equivol_per_depth.loc[idx, 'ecc'] = pd.Series([eccxlayer]*len(idx), index=df_equivol_per_depth.index[idx])
            df_equivol_per_depth.loc[idx, 'rsq'] = pd.Series([rsqxlayer]*len(idx), index=df_equivol_per_depth.index[idx])
            
            depth_idx += 1
        ###########################################################################################
        ## Extract sigma per depth: equidistant
        depth_idx=0
        for layer in range(2,np.size(np.unique(masked_lay_equidist))):
            sigmaxlayer = sigma[(masked_V1==1) & (masked_lay_equidist == layer)]
            eccxlayer = ecc[(masked_V1==1) & (masked_lay_equidist == layer)]
            rsqxlayer = total_rsq[(masked_V1==1) & (masked_lay_equidist == layer)]

            idx=df_equidist_per_depth.loc[(df_equidist_per_depth['sub id'] == subject_list[sub_id]) & 
                                       (df_equidist_per_depth['hem'] == hem_list[hem_id]) & 
                                       (df_equidist_per_depth['depth'] == depth[depth_idx])].index.tolist()
            df_equidist_per_depth.loc[idx, 'pRF size'] = pd.Series([sigmaxlayer]*len(idx), index=df_equidist_per_depth.index[idx])
            df_equidist_per_depth.loc[idx, 'ecc'] = pd.Series([eccxlayer]*len(idx), index=df_equidist_per_depth.index[idx])
            df_equidist_per_depth.loc[idx, 'rsq'] = pd.Series([rsqxlayer]*len(idx), index=df_equidist_per_depth.index[idx])
            depth_idx += 1



# save df
f = open(opj(prfpy_dir,'..','df_layers'), 'wb')
pickle.dump([df_equivol_per_ecc,df_equidist_per_ecc,df_equivol_per_depth,df_equidist_per_depth], f)
f.close()