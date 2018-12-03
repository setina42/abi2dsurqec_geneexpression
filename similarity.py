import numpy as np
import nibabel
import os
import glob
from nipype.interfaces import  ants as ants
import nipype.interfaces.fsl.maths as fsl
import csv

#TODO: correct???
def transform(x,y,z,affine):
  M = affine[:4, :4]
  A = np.linalg.inv(M)
  A = np.round(A,decimals=0) #really?
  return np.round(A.dot([x,y,z,1]),decimals=0) #?


#TODO:Either already do that when downloading and distribute mirrored images, or else keep them for later use or delete them right after similarity is measured??? Or even distribute raw and mirrored data?
#depending on size /calculation time
def mirror_sagittal(image):
    img = nibabel.load(image)
    img_data = img.get_fdata()
    
    #find coordinates of origin (Bregma?) in data matrix to determine the left-right midline
    origin = transform(0,0,0,img.affine)
    mid = int(origin[0])
    #Now find the first LR slice that has only -1 in it
#    start = origin[0] + 2
#    done = False
#    while not done:
#      if  np.max(img_data[start,:,:]) > 0:
#        start = start - 1   #check what side i want to go beforehand just to be safe
#      else:
#        done = True
#        midline = start
        
#TODO: Make a distinction between even and odd number of x-slices


    #Now we have the first slice where no data was collected. Now, copy and flip???
    left_side = np.copy(img_data[(mid + 1) : np.shape(img_data)[0],:,: ])
    left_side = np.flip(left_side,0)
    #replace
    if np.shape(left_side)[0] >  np.shape(img_data[0:mid-1,:,:])[0]:
      #case 1: origin slightly to the left (or right??), need to trim left_side to the size of the right side
      replace_value = np.shape(left_side)[0] - np.shape(img_data[0:mid,:,:])[0]
      img_data[0:mid,:,:] = left_side[replace_value:np.shape(left_side)[0],:,:]
    elif np.shape(left_side)[0] >  np.shape(img_data[0:mid,:,:])[0]:
      #case 2 : origin slightly to the right (or left??), need to 
      replace_start = np.shape(img_data[0:(mid-1),:,:])[0] - np.shape(left_side)[0]
      img_data[replacevalue:(mid-1),:,:] = left_side
    else:
      #case 3: same size 
      img_data[0:mid-1,:,:] = left_side
    
    img_average = nibabel.Nifti1Image(img_data,img.affine)
    
    filename = str.split(os.path.basename(image),'.nii')[0] + '_mirrored.nii.gz'
    path_to_mirrored = os.path.join(os.path.dirname(image),filename)
    nibabel.save(img_average,path_to_mirrored)
    print("saved")


    z = 0
    return path_to_mirrored

def create_mask(img):
  print("reached")
  mask = fsl.Threshold()
  mask.inputs.thresh = 0
  mask.inputs.args = '-bin'  
  mask.inputs.in_file = img
  img_out = str.split(img,'.nii')[0] + 'mask.nii.gz'
  mask.inputs.out_file = img_out
  mask.run()

  return img_out


def create_experiment_average(imgs,strategy='max'):
    img_data = []
    img = []
    for image in imgs:
        img_2 = nibabel.load(image)
        img_data.append(img_2.get_fdata())
        img.append(img_2)

    if strategy == 'max':
        print(len(img_data))
        print(np.shape(img_data[0]))
        average_img = np.maximum(img_data[0],img_data[1])
        if len(img_data)>2:
          for i in range(2,(len(img_data)-1)):
            average_img = np.maximum(average_img,img_data[i])

    elif strategy == 'min':
        print(len(img_data))
        print(np.shape(img_data[0]))
        average_img = np.minimum(img_data[0],img_data[1])
        if len(img_data)>2:
          for i in range(2, (len(img_data)-1)):
            average_img = np.minimum(average_img,img_data[i])

    elif strategy == 'mean':
        average_img = (np.add(*img_data))/float(len(img_data))
    
    filename = str.split(os.path.basename(imgs[0]),"_")[0] + "_experiments_average.nii.gz"
    path_to_exp_average = os.path.join(os.path.dirname(imgs[0]),"..")
    path_to_exp_average = os.path.join(path_to_exp_average,filename)
    img_average = nibabel.Nifti1Image(average_img,img[0].affine)
    nibabel.save(img_average,path_to_exp_average)
    return path_to_exp_average

#seems to work as intended,tested only without mask :)
def ants_measure_similarity(fixed_image,moving_image,mask_gene = None,mask_map = None,metric = 'GC',metric_weight = 1.0,radius_or_number_of_bins = 64,sampling_strategy='Regular',sampling_percentage=1.0):
    sim = ants.MeasureImageSimilarity()
    sim.inputs.dimension = 3
    sim.inputs.metric = metric
    sim.inputs.fixed_image = fixed_image
    sim.inputs.moving_image = moving_image
    sim.inputs.metric_weight = metric_weight
    sim.inputs.radius_or_number_of_bins = radius_or_number_of_bins
    sim.inputs.sampling_strategy = sampling_strategy
    sim.inputs.sampling_percentage = sampling_percentage
    sim.inputs.fixed_image_mask = mask_map
    sim.inputs.moving_image_mask = mask_gene

    sim_res = sim.run()
    return sim_res.outputs.similarity


def measure_similarity_geneexpression(stat_map,path_to_genes="/home/gentoo/src/abi2dsurqec_geneexpression/ABI_geneexpression_data",metric = 'MeanSquares',radius_or_number_of_bins = 64,create_exp_average = True,strategy='max'):
    
    #TODO: create a mask for the stat map? or userprovided? or both possible?
    #TODO case: no exp average is wanted
    #img = nibabel.load(stat_map)
    #img_data = img.get_fdata()
    mask_map = create_mask(stat_map)
    results = dict()
    #loop through all gene folders, either get data form single experiment or get combined data.
    for dir in os.listdir(path_to_genes):
        path = os.path.join(path_to_genes,dir)
        print(dir)
        #multiple experiment data available per gene
        if len(os.listdir(path)) > 1:
            img = []
            imgs = glob.glob(path + '/*/*_2dsurqec.nii.gz')
            for img_gene in imgs:
              if "sagittal" in img:
                 img.append(mirror_sagittal(img_gene))
              else:
                img.append(img_gene)
            img_gene = create_experiment_average(img,strategy=strategy)
        elif len(os.listdir(path)) == 1:
            img_gene = glob.glob(path + '/*/*_2dsurqec.nii.gz')[0]
            if "sagittal" in img_gene:
                img_gene = mirror_sagittal(img_gene)
        else:
          print("Folder empty or no registered niftifile found. Or weird folder name :)")
          print("Skipping " + dir)
          break
        #TODO: catch unexpected errors as to not interrupt program, print genename
        mask_gene = create_mask(img_gene)
        similarity = ants_measure_similarity(stat_map,img_gene,mask_gene = mask_gene,mask_map=mask_map,metric=metric,radius_or_number_of_bins=radius_or_number_of_bins)
        results[dir] = similarity

    name = metric + "_" + str(radius_or_number_of_bins) +".csv"

    with open(name,'w') as f:
      for key in results.keys():
        f.write("%s,%s\n"%(key,results[key]))


def main():
  img = "/home/gentoo/src/abi2dsurqec_geneexpression/ABI_geneexpression_data/Mef2c/Mef2c_P56_sagittal_79677145_200um/Mef2c_P56_sagittal_79677145_200um_2dsurqec.nii.gz"
 #res = ants_measure_similarity("/home/gentoo/src/abi2dsurqec_geneexpression/ABI_geneexpression_data/Mef2c/Mef2c_P56_sagittal_79677145_200um/Mef2c_P56_sagittal_79677145_200um_2dsurqec.nii.gz","/home/gentoo/src/abi2dsurqec_geneexpression/ABI_geneexpression_data/Mef2c/Mef2c_P56_coronal_79567505_200um//Mef2c_P56_coronal_79567505_200um_2dsurqec.nii.gz")
  #  print(res)
  measure_similarity_geneexpression("/home/gentoo/src/abi2dsurqec_geneexpression/ABI_geneexpression_data/Mef2c/Mef2c_P56_coronal_79567505_200um/Mef2c_P56_coronal_79567505_200um_2dsurqec.nii.gz",metric = 'GC',radius_or_number_of_bins = 64)
    
#    measure_similarity_geneexpression("/home/gentoo/src/abi2dsurqec_geneexpression/ABI_geneexpression_data/Mef2c/Mef2c_P56_coronal_79567505_200um/Mef2c_P56_coronal_79567505_200um_2dsurqec.nii.gz",metric = 'CC',radius_or_number_of_bins = 4)

#    measure_similarity_geneexpression("/home/gentoo/src/abi2dsurqec_geneexpression/ABI_geneexpression_data/Mef2c/Mef2c_P56_coronal_79567505_200um/Mef2c_P56_coronal_79567505_$200um_2dsurqec.nii.gz",metric = 'Mattes',radius_or_number_of_bins = 32)

#    measure_similarity_geneexpression("/home/gentoo/src/abi2dsurqec_geneexpression/ABI_geneexpression_data/Mef2c/Mef2c_P56_coronal_79567505_200um/Mef2c_P56_coronal_79567505_200um_2dsurqec.nii.gz",metric = 'MeanSquares',radius_or_number_of_bins = 64)

#    measure_similarity_geneexpression("/home/gentoo/src/abi2dsurqec_geneexpression/ABI_geneexpression_data/Mef2c/Mef2c_P56_coronal_79567505_200um/Mef2c_P56_coronal_79567505_200um_2dsurqec.nii.gz")
#  create_mask(img)

if __name__ == "__main__":
        main()

