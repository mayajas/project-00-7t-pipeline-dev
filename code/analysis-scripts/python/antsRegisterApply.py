from os.path import join as opj

from nipype.interfaces.ants.registration import Registration
from nipype.interfaces.ants import N4BiasFieldCorrection
from nipype.interfaces.ants import ApplyTransforms


from nipype.interfaces.fsl.maths import BinaryMaths

sub_id='sub-01'
data_dir='/home/mayajas/Documents/project-00-7t-pipeline-dev/tests/'+sub_id

reg = Registration()
reg.inputs.verbose = True
reg.inputs.dimension = 3
reg.inputs.float = True
reg.inputs.output_transform_prefix = 'registered_'
reg.inputs.output_warped_image = data_dir+'/registered_Warped.nii.gz'
reg.inputs.output_inverse_warped_image = data_dir+'/registered_InverseWarped.nii.gz'
reg.inputs.interpolation = 'Linear'
reg.inputs.use_histogram_matching = [False, False, False] # This is the default
reg.inputs.winsorize_lower_quantile = 0.005
reg.inputs.winsorize_upper_quantile = 0.995
reg.inputs.initial_moving_transform = data_dir+'/coreg_itksnap_func2struct.txt'

reg.inputs.fixed_image = data_dir+'/UNI_corrected.nii'
reg.inputs.moving_image = data_dir+'/merged_func_mcf.nii_mean_reg.nii'

reg.inputs.transforms = ['Rigid','Affine','SyN']
reg.inputs.transform_parameters = [(0.1,), (0.1,), (0.1, 3.0, 0.0)]
reg.inputs.number_of_iterations = [[1000,500,250,100], [1000,500,250,100],[100,70,50,20]]

reg.inputs.metric = ['MI','MI','CC']
reg.inputs.metric_weight = [1,1,1] # Default (value ignored currently by ANTs)
reg.inputs.radius_or_number_of_bins = [32,32,4]
reg.inputs.sampling_strategy = ['Regular', 'Regular','None']
reg.inputs.sampling_percentage = [0.25, 0.25, None]

reg.inputs.convergence_threshold = [1.e-6, 1.e-6, 1.e-6]
reg.inputs.convergence_window_size = [10,10,10]
reg.inputs.smoothing_sigmas = [[3,2,1,0], [3,2,1,0], [3,2,1,0]]
reg.inputs.sigma_units = ['vox'] * 3
reg.inputs.shrink_factors = [[8,4,2,1], [8,4,2,1], [8,4,2,1]]

reg.inputs.fixed_image_masks = [data_dir+'/mask_midocc.nii']
reg.run()
