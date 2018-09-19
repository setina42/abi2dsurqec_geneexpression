#!/bin/sh


WarpImageMultiTransform 3 energy.nii.gz warped_direct.nii.gz abi2dsurqec_0GenericAffine.mat -R dsurqec_200micron.nii.gz

fslorient -copyqform2sform warped_direct.nii.gz

