#!/bin/sh

#CompositeTransformUtil --disassemble abi2dsurqec_Composite.h5 Dissasembled

#WarpImageMultiTransform 3 energy.nii.gz warped_withref200.nii.gz 01_Dissasembled_DisplacementFieldTransform.nii.gz 00_Dissasembled_AffineTransform.mat -R dsurqec_200micron_masked.nii.gz
#fslorient -copyqform2sform warped_withref200.nii.gz


antsApplyTransforms --d 3 -i energy.nii.gz -r dsurqec_200micron_masked.nii -o warped_applyTransform.nii.gz -t abi2dsurqec_Composite.h5
fslorient -copyqform2sform warped_withref200.nii.gz

