SUBJECT_ID="sub-0$1"
HEMI="$2h"

PRF_MODEL=prfpy_fit_hrf_True_start_from_avg_True

echo "Merging ventral and dorsal manually-defined ROI labels for subject $SUBJECT_ID"

# define directories
PROJECT_ID=project-00-7t-pipeline-dev
PROJECT_DIR=/home/mayajas/scratch/$PROJECT_ID
SUB_DIR=$PROJECT_DIR/output/${PRF_MODEL}/$SUBJECT_ID

# merge ventral and dorsal labels (V1, V2, V3)
mri_mergelabels -i $SUB_DIR/${HEMI}.V1d.label -i $SUB_DIR/${HEMI}.V1v.label -o $SUB_DIR/${HEMI}.V1.label
mri_mergelabels -i $SUB_DIR/${HEMI}.V2d.label -i $SUB_DIR/${HEMI}.V2v.label -o $SUB_DIR/${HEMI}.V2.label
mri_mergelabels -i $SUB_DIR/${HEMI}.V3d.label -i $SUB_DIR/${HEMI}.V3v.label -o $SUB_DIR/${HEMI}.V3.label