import nibabel as nib
import numpy as np
from scipy.spatial import KDTree
from os.path import join as opj


#################################################################################################################
# Define important variables
subject_list = ['sub-01','sub-02','sub-03','sub-04']
prf_model    = 'prfpy_Iso2DGaussianModel'#_LinearDeveining'


#################################################################################################################
## Loop over subjects
for sub_id in range(0,len(subject_list)):
    prfpy_dir    = '/home/mayajas/scratch/project-00-7t-pipeline-dev/output/'+prf_model+'/'+subject_list[sub_id]
    
    lay_dir      = opj(prfpy_dir,'layerification_func')

    # Load the rim label mask as a Nifti file
    mask_img = nib.load(opj(lay_dir,'func_ribbon_rim.nii'))
    mask = mask_img.get_fdata()

    # Get the coordinates of all white matter voxels
    white_matter_coords = np.array(np.where(mask == 2)).T

    # Get the coordinates of all gray matter voxels
    gray_matter_coords = np.array(np.where(mask == 3)).T

    ## Distance from WM
    # Create a KDTree to efficiently search for the nearest white matter voxel
    WM_tree = KDTree(white_matter_coords)

    # For each gray matter voxel, find the distance to the nearest white matter voxel
    result = WM_tree.query(gray_matter_coords)
    indices = result[1]
    distances = result[0]

    distances_mm = distances*0.8

    # # Set all gm and wm distances to 0 a priori
    # mask[mask == 3] = np.zeros(mask[mask == 3].shape)
    # mask[mask == 2] = np.zeros(mask[mask == 2].shape)

    # Store the distances as values in the mask for each gray matter voxel
    mask[mask == 3] = distances_mm

    ## Distance from GM
    # Create a KDTree to efficiently search for the nearest white matter voxel
    GM_tree = KDTree(gray_matter_coords)

    # For each gray matter voxel, find the distance to the nearest white matter voxel
    result = GM_tree.query(white_matter_coords)
    indices = result[1]
    distances = result[0]

    distances_mm = -distances*0.8

    # Store the distances as values in the mask for each gray matter voxel
    mask[mask == 2] = distances_mm

    # Create a new Nifti file with the updated mask
    mask_with_distances_img = nib.Nifti1Image(mask, mask_img.affine, header=mask_img.header)

    # Save the updated mask
    nib.save(mask_with_distances_img, opj(lay_dir,'cortical_distances.nii'))