SUBJECT_ID="sub-0$1"

# define directories
LAYNII_DIR=/home/mayajas/Documents/programs/laynii
RIM_SCRIPT_DIR=/home/mayajas/Documents/project-00-7t-pipeline-dev/code/analysis-scripts/python
PRFPY_DIR=/home/mayajas/scratch/project-00-7t-pipeline-dev/output/prfpy/$SUBJECT_ID
LAYER_DIR=$PRFPY_DIR/layerification
COREG_DIR=/home/mayajas/scratch/project-00-7t-pipeline-dev/manualcorr/coreg/_subject_id_$SUBJECT_ID

# define anatomical and func images
FUNC_IMG=/home/mayajas/scratch/project-00-7t-pipeline-dev/output/prfpy/$SUBJECT_ID/meanFunc.nii

# make layer output directory
if [ ! -d "$LAYER_DIR" ]; then
  mkdir $LAYER_DIR
  echo "Making directory ${LAYER_DIR}..."
else
  echo "Layerification outputs will be saved to ${LAYER_DIR}..."
fi

# convert ribbon mgz to nii
mri_convert $PRFPY_DIR/ribbon.mgz $LAYER_DIR/ribbon.nii

# make laynii required rim file
python $RIM_SCRIPT_DIR/makeLayniiRim.py $LAYER_DIR/ribbon.nii

# laynii layerification
$LAYNII_DIR/LN2_LAYERS -rim $LAYER_DIR/ribbon_rim.nii -nr_layers 10 -equivol

# coreg to func space
# equidistant image:
INPUT=$LAYER_DIR/ribbon_rim_layers_equidist.nii
OUTPUT=$LAYER_DIR/func_ribbon_rim_layers_equidist.nii
${ANTSPATH}antsApplyTransforms \
  -d 3 \
  -i $INPUT \
  -r $FUNC_IMG \
  -t [ $COREG_DIR/registered_0GenericAffine.mat, 1 ] \
  -o $OUTPUT \
  --interpolation NearestNeighbor

# equivol image:
INPUT=$LAYER_DIR/ribbon_rim_layers_equivol.nii
OUTPUT=$LAYER_DIR/func_ribbon_rim_layers_equivol.nii
${ANTSPATH}antsApplyTransforms \
  -d 3 \
  -i $INPUT \
  -r $FUNC_IMG \
  -t [ $COREG_DIR/registered_0GenericAffine.mat, 1 ] \
  -o $OUTPUT \
  --interpolation NearestNeighbor














# # GM/WM values in ribbon.mgz
# rh_GM_val=42
# rh_WM_val=41
# lh_GM_val=3
# lh_WM_val=2

# # get nii of ribbon image
# mri_convert ribbon.mgz ribbon.nii

# # extract white matter, borderize
# mri_binarize --i ribbon.nii --o lh_WM.nii --match $lh_WM_val
# $LAYNII_DIR/LN2_BORDERIZE -input lh_WM.nii -jumps 1

# # extract gray matter
# mri_binarize --i ribbon.nii --o lh_GM.nii --match $lh_GM_val

# # get sum of WM borders and GM, borderize
# 3dcalc -a lh_WM_borders.nii -b lh_GM.nii -prefix WMGM_sum.nii -overwrite -expr 'a+b'
# $LAYNII_DIR/LN2_BORDERIZE -input WMGM_sum.nii -jumps 1

# # separately label external and internal borders
# 3dcalc -a WMGM_sum_borders.nii -b lh_WM_borders.nii -prefix lh_borders.nii -overwrite -expr 'a+b'
# mri_binarize --i lh_borders.nii --match 1 --o lh_GM_external.nii
# mri_binarize --i lh_borders.nii --match 2 --o lh_GM_internal.nii

# # make rim image
# 3dcalc -a lh_GM_external.nii -b lh_GM_internal.nii -c lh_GM.nii -prefix lh_GM_rim.nii -overwrite -expr 'a+2*b+3*c-1'

# # compare with rimify output
# $LAYNII_DIR/LN2_RIMIFY -input lh_GM.nii -innergm 0 -outergm 0 -gm 1 -output lh_GM_layniirim.nii

# # run layerification
# $LAYNII_DIR/LN2_LAYERS -rim lh_GM_layniirim.nii -nr_layers 10 -equivol