import numpy
import nibabel
import os
import glob
from nipype.interfaces import  ants as ants
from nipype.interfaces import fsl

def mirror_sagittal(img):
    z = 0

def create_mask(img):
    z = 0
def create_experiment_average(imgs,strategy='max'):
    z = 0

    return path_to_exp_average

#seems to work as intended :)
def ants_measure_similarity(fixed_image,moving_image,metric = 'MI',metric_weight = 1.0,radius_or_number_of_bins = 64,sampling_strategy='Regular',sampling_percentage=1.0):
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


def measure_similarity_geneexpression(stat_map,path_to_genes="/home/gentoo/src/abi2dsurqec_geneexpression/ABI_geneexpression_data",create_exp_average = True,strategy='max'):
    #img = nibabel.load(stat_map)
    #img_data = img.get_fdata()

    #loop through all gene folders, either get data form single experiment or get combined data.
    for dir in os.listdir(path_to_genes):
        path = os.path.join(path_to_genes,dir)
        if len(os.listdir(path)) > 1:
            imgs = glob.glob(path + '/*/*_2dsurqec.nii.gz')
            img_gene = create_experiment_average(imgs,strategy=strategy)
        elif len(os.listdir(path)) == 1:
            img_gene = glob.glob(path + '/*/*_2dsurqec.nii.gz')
        similarity = ants_measure_similarity(stat_map,img_gene)



    #calculate experiment averages


    #separate function

def main():
    res = ants_measure_similarity("/home/gentoo/src/abi2dsurqec_geneexpression/ABI_geneexpression_data/Mef2c/Mef2c_P56_sagittal_79677145_200um/Mef2c_P56_sagittal_79677145_200um_2dsurqec.nii.gz","/home/gentoo/src/abi2dsurqec_geneexpression/ABI_geneexpression_data/Mef2c/Mef2c_P56_coronal_79567505_200um//Mef2c_P56_coronal_79567505_200um_2dsurqec.nii.gz")
    print(res)
    #measure_similarity_geneexpression("bla")


if __name__ == "__main__":
        main()

