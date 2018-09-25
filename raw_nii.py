import numpy, sys
import mhd_utils_3d
import nibabel
from mhd_utils_3d import *
import matplotlib.pyplot as plt

#Read raw data and meta  
image_array, meta_header = mhd_utils_3d.load_raw_data_with_mhd("energy.mhd")


#Read header infomormation and create affine matrix    
dims = numpy.array(meta_header["ElementSpacing"].split(" "),dtype=numpy.float)
affine_matrix = numpy.zeros((4,4),dtype=numpy.float)
affine_matrix[0,0] = dims[0]
affine_matrix[1,1] = dims[1]
affine_matrix[2,2] = dims[2]


#Bring to same orientation as the converted atlas.nii
#Super unsafe? :(
if meta_header['AnatomicalOrientation'] == 'RAI':
    #affine_matrix[0,0] = affine_matrix[0,0] * (-1)
    #affine_matrix[1,1] = affine_matrix[1,1] * (-1)
    #affine_matrix[2,2] = affine_matrix[2,2] * (-1)
    zu = 3
    
#image_array = numpy.swapaxes(image_array,0,2)
image_array = numpy.swapaxes(image_array,1,2)

#image_array = image_array[:,:,::-1]
#image_array = image_array[:,::-1,:]
image_array = image_array[::-1,:,:]
image_array = image_array[:,:,::-1]
image_array = image_array[:,::-1,:]

#Bring to the right units?
affine_matrix = affine_matrix*0.001

affine_matrix[3,3] = 1  #?????


img = nibabel.Nifti1Image(image_array,affine_matrix)
nibabel.save(img,"energy.nii.gz")


