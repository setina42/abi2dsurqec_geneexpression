import numpy as np
import nibabel
import os
import glob
from nipype.interfaces import  ants as ants
from nipype.interfaces import fsl
import csv


def z(x,y,z,affine):
  M = affine[:4, :4]
  A = numpy.linalg.inv(M)
  A = numpy.round(A,decimals=0)
  return numpy.round(A.dot([x,y,z,1]),decimals=0)


#TODO:Either already do that when downloading and distribute mirrored images, or else keep them for later use or delete them right after similarity is measured??? Or even distribute raw and mirrored data?
#depending on size /calculatio time
def mirror_sagittal(imge):
    img = nibabel.load(image)
    img_data = img_get_fdata()
    
    #find coordinates of origin (Bregma?) in data matrix to determine the left-right midline
    origin = z(0,0,0,img.affine)
    origin[0] = int(origin[0])
    origin[1] = int(origin[1]
    origin[2] = int(origin[2])

    #Now find the first LR slice that has only -1 in it
    start = origin[0] + 2
    done = False
    while not done:
      if  numpy.max(img_data[start,:,:]) > 0:
        start = start - 1   #check what side i want to go beforehand just to be safe
      else:
        done = True
        midline = start
        
    #Now we have the first slice where no data was collected. Now, copy and flip???
    left_side = numpy.copy(img_data[(start + 1) : numpy.shape(img_data)[0],:,: ])
    left_side = numpy.flip(left_side,0)
    #replace
    if numpy.shape(left_side)[0] >  numpy.shape(img_data[0:start,:,:])[0]:
      #case 1: origin slightly to the left (or right??), need to trim left_side to the size of the right side
      replace_value = numpy.shape(left_side)[0] - numpy.shape(img_data[0:start,:,:])[0]
      img_data[0:start,:,:] = left_side[replace_value:numpy.shape(left_side)[0],:,:]
    elif numpy.shape(left_side)[0] >  numpy.shape(img_data[0:start,:,:])[0]:
      replace_start = numpy.shape(img_data[0:start,:,:])[0] - numpy.shape(left_side)[0]
    else:
      #case 3: same size 
      img_data[0:start,:,:] = left_side
      img_average = nibabel.Nifti1Image(img_data,img.affine)
    
    nibabel.save(img_average,"C:/Users/tinas/Desktop/Sync/Rec/Stx7_P56_sagittal_69887415_200um/Stx7_P56_sagittal_69887415_200um_2dsurqec_mirrored.nii.gz")
    print("saved")


    z = 0

def create_mask(img):
    z = 0


def create_experiment_average(imgs,strategy='max'):
    img_data = []
    img = []

    for image in imgs:
        img_2 = nibabel.load(image)
        img_data.append(img_2.get_fdata())
        img.append(img_2)

    if strategy == 'max':
        average_img = np.maximum(*img_data)

    elif strategy == 'min':
        average_img = np.minimum(*img_data)

    elif strategy == 'mean':
        average_img = (np.add(*img_data))/float(len(img_data))
    filename = str.split(os.path.basename(imgs[0]),"_")[0] + "_experiments_average.nii.gz"
    path_to_exp_average = os.path.join(os.path.dirname(imgs[0]),filename)
    img_average = nibabel.Nifti1Image(average_img,img[0].affine)
    nibabel.save(img_average,path_to_exp_average)
    return path_to_exp_average

#seems to work as intended :)
def ants_measure_similarity(fixed_image,moving_image,metric = 'GC',metric_weight = 1.0,radius_or_number_of_bins = 64,sampling_strategy='Regular',sampling_percentage=1.0):
    sim = ants.MeasureImageSimilarity()
    sim.inputs.dimension = 3
    sim.inputs.metric = metric
    sim.inputs.fixed_image = fixed_image
    sim.inputs.moving_image = moving_image
    sim.inputs.metric_weight = metric_weight
    sim.inputs.radius_or_number_of_bins = radius_or_number_of_bins
    sim.inputs.sampling_strategy = sampling_strategy
    sim.inputs.sampling_percentage = sampling_percentage
    #sim.inputs.fixed_image_mask = 
    #sim.inputs.moving_image_mask = 

    sim_res = sim.run()
    return sim_res.outputs.similarity


def measure_similarity_geneexpression(stat_map,path_to_genes="/home/gentoo/src/abi2dsurqec_geneexpression/ABI_geneexpression_data",metric = 'MeanSquares',radius_or_number_of_bins = 64,create_exp_average = True,strategy='max'):
    #img = nibabel.load(stat_map)
    #img_data = img.get_fdata()
    results = dict()
    #loop through all gene folders, either get data form single experiment or get combined data.
    for dir in os.listdir(path_to_genes):
        path = os.path.join(path_to_genes,dir)
        if len(os.listdir(path)) > 1:
            imgs = glob.glob(path + '/*/*_2dsurqec.nii.gz')
            img_gene = create_experiment_average(imgs,strategy=strategy)
        elif len(os.listdir(path)) == 1:
            img_gene = glob.glob(path + '/*/*_2dsurqec.nii.gz')[0]
        similarity = ants_measure_similarity(stat_map,img_gene,metric=metric,radius_or_number_of_bins=radius_or_number_of_bins)
        results[dir] = similarity

    name = metric + "_" + str(radius_or_number_of_bins) +".csv"

    with open(name,'w') as f:
      for key in results.keys():
        f.write("%s,%s\n"%(key,results[key]))

    #calculate experiment averages


    #separate function

def main():
    #res = ants_measure_similarity("/home/gentoo/src/abi2dsurqec_geneexpression/ABI_geneexpression_data/Mef2c/Mef2c_P56_sagittal_79677145_200um/Mef2c_P56_sagittal_79677145_200um_2dsurqec.nii.gz","/home/gentoo/src/abi2dsurqec_geneexpression/ABI_geneexpression_data/Mef2c/Mef2c_P56_coronal_79567505_200um//Mef2c_P56_coronal_79567505_200um_2dsurqec.nii.gz")
  #  print(res)
#    measure_similarity_geneexpression("/home/gentoo/src/abi2dsurqec_geneexpression/ABI_geneexpression_data/Mef2c/Mef2c_P56_coronal_79567505_200um/Mef2c_P56_coronal_79567505_200um_2dsurqec.nii.gz",metric = 'Mattes',radius_or_number_of_bins = 128)
    
    measure_similarity_geneexpression("/home/gentoo/src/abi2dsurqec_geneexpression/ABI_geneexpression_data/Mef2c/Mef2c_P56_coronal_79567505_200um/Mef2c_P56_coronal_79567505_200um_2dsurqec.nii.gz",metric = 'CC',radius_or_number_of_bins = 4)

#    measure_similarity_geneexpression("/home/gentoo/src/abi2dsurqec_geneexpression/ABI_geneexpression_data/Mef2c/Mef2c_P56_coronal_79567505_200um/Mef2c_P56_coronal_79567505_$200um_2dsurqec.nii.gz",metric = 'Mattes',radius_or_number_of_bins = 32)

#    measure_similarity_geneexpression("/home/gentoo/src/abi2dsurqec_geneexpression/ABI_geneexpression_data/Mef2c/Mef2c_P56_coronal_79567505_200um/Mef2c_P56_coronal_79567505_200um_2dsurqec.nii.gz",metric = 'MeanSquares',radius_or_number_of_bins = 64)

#    measure_similarity_geneexpression("/home/gentoo/src/abi2dsurqec_geneexpression/ABI_geneexpression_data/Mef2c/Mef2c_P56_coronal_79567505_200um/Mef2c_P56_coronal_79567505_200um_2dsurqec.nii.gz")


if __name__ == "__main__":
        main()

