import numpy as np
import nibabel
import os
import glob
from nipype.interfaces import ants as ants
import nipype.interfaces.fsl.maths as fsl
import csv

def transform(x,y,z,affine):
	M = affine[:4, :4]
	A = np.linalg.inv(M)
	print("unrounded")
	print(A.dot([x,y,z,1]))
	return np.round(A.dot([x,y,z,1]),decimals=0)

def mirror_sagittal(image):
	"""
	Sagittal datasets form Allen mouse brain are only collected for the right/left? hemisphere.
	Function to mirror feature map at midline and saving image as nifti-file.
	"""
	img = nibabel.load(image)
	img_data = img.get_fdata()

	#find coordinates of origin (Bregma?) in data matrix to determine the left-right midline
	origin = transform(0,0,0,img.affine)
	mid = int(origin[0])

	#TODO:midline point at 31.34, how to mirror properly

	#Copy image right/left(?) from the mid, but keep mid unchanged

	left_side = np.copy(img_data[(mid + 1):,:,: ])
	left_side = np.flip(left_side,0)

	right_side = np.copy(img_data[0:mid,:,:])
	right_side = np.flip(right_side,0)

	#replace
	if np.shape(left_side)[0] > np.shape(img_data[0:mid,:,:])[0]:
		print("case 1")
		#case 1: origin slightly to the left (or right??), need to trim left_side to the size of the right side
		replace_value = np.shape(left_side)[0] - np.shape(img_data[0:mid,:,:])[0]
		img_data[0:mid,:,:][img_data[0:mid,:,:]  == -1] = left_side[(replace_value-1):,:,:][img_data[0:mid,:,:]  == -1]
		np.savetxt("Slice_case1.txt",img_data[:,(int(origin[1]) -5):(int(origin[1])) -2, (int(origin[2]) -5)])
		img_data[mid:np.shape(right_side[0]),:,:][img_data[mid:,:,:]  == -1] = right_side[:,:,:][img_data[mid:,:,:]  == -1]

	elif np.shape(left_side)[0] < np.shape(img_data[0:mid,:,:])[0]:
		print("case 2")
		#case 2 : origin slightly to the right (or left??), need to
		replace_value = np.shape(img_data[0:mid,:,:])[0] - np.shape(left_side)[0]
		img_data[replace_value:mid,:,:] = left_side
		np.savetxt("Slice_case2.txt",img_data[:,(int(origin[1]) -5):(int(origin[1])) -2, (int(origin[2]) -5)])

		img_data[mid:,:,:][img_data[mid:,:,:]  == -1] = right_side[0:np.shape(img_data[mid:,:,:]),:,:][img_data[mid:,:,:]  == -1]
	else:
		print("case 3")
		#case 3: same size
		#TODO: -1 or 0??
		#TODO: There's got to be a better way to write this....
		img_data[0:mid,:,:][img_data[0:mid,:,:]  == -1] = left_side[img_data[0:mid,:,:]  == -1]
		np.savetxt("Slice_case3.txt",img_data[:,(int(origin[1]) -5):(int(origin[1])) -2, (int(origin[2]) -5)])
		img_data[mid+1:,:,:][img_data[mid+1:,:,:]  == -1] = right_side[[img_data[mid+1:,:,:]  == -1]]

	img_average = nibabel.Nifti1Image(img_data,img.affine)
	filename = str.split(os.path.basename(image),'.nii')[0] + '_mirrored.nii.gz'
	path_to_mirrored = os.path.join(os.path.dirname(image),filename)
	nibabel.save(img_average,path_to_mirrored)

	return path_to_mirrored

def create_mask(image):
	#I think i need 0 to be in my mask. This seems not to be possible using fslmaths, so maybe do directly with numpy? thr sets all to zero below the value and bin uses image>0 to binarise.
#	mask = fsl.Threshold()
#	mask.inputs.thresh = 0
#	mask.inputs.args = '-bin'
#	mask.inputs.in_file = image
#	img_out = str.split(image,'.nii')[0] + '_mask.nii.gz'
#	mask.inputs.out_file = img_out
#	mask.run()

	mask_img = nibabel.load("/home/gentoo/src/abi2dsurqec_geneexpression/dsurqec_200micron_mask.nii")
	atlas_mask = mask_img.get_fdata()

	#using numpy instead of fslmaths
	img = nibabel.load(image)
	img_data = img.get_fdata()
	#img_data[img_data >= 0 and atlas_mask == 1] = 1
	img_data[np.logical_and(img_data >= 0,atlas_mask == 1)] = 1
	img_data[img_data < 0] = 0
	img_out = str.split(image,'.nii')[0] + '_mask.nii.gz'
	img_mask = nibabel.Nifti1Image(img_data,img.affine)
	nibabel.save(img_mask,img_out)

	return img_out

def create_experiment_average(imgs,strategy='max'):
	"""
	In case of several datasets present, experiment average is calculated.
	"""

	img_data = []
	img = []
	for image in imgs:
		img_2 = nibabel.load(image)
		img_data.append(img_2.get_fdata())
		img.append(img_2)

	if strategy == 'max':
		average_img = np.maximum(img_data[0],img_data[1])
		if len(img_data)>2:
			for i in range(2,(len(img_data)-1)):
				average_img = np.maximum(average_img,img_data[i])

	elif strategy == 'mean':
		#TODO Ignore values with -1 for averaging, take only other values
		average_img = (np.add(*img_data))/float(len(img_data))

	filename = str.split(os.path.basename(imgs[0]),"_")[0] + "_experiments_average.nii.gz"
	path_to_exp_average = os.path.join(os.path.dirname(imgs[0]),"..")
	path_to_exp_average = os.path.join(path_to_exp_average,filename)
	img_average = nibabel.Nifti1Image(average_img,img[0].affine)
	nibabel.save(img_average,path_to_exp_average)
	return path_to_exp_average

#seems to work as intended,tested only without mask :)
#TODO: masking working correctly?
def ants_measure_similarity(fixed_image,moving_image,mask_gene = None,mask_map = None,metric = 'MI',metric_weight = 1.0,radius_or_number_of_bins = 64,sampling_strategy='Regular',sampling_percentage=1.0):
	"""
	Nipype ants
	"""
	sim = ants.MeasureImageSimilarity()
	sim.inputs.dimension = 3
	sim.inputs.metric = metric
	sim.inputs.fixed_image = fixed_image
	sim.inputs.moving_image = moving_image
	sim.inputs.metric_weight = metric_weight
	sim.inputs.radius_or_number_of_bins = radius_or_number_of_bins
	sim.inputs.sampling_strategy = sampling_strategy
	sim.inputs.sampling_percentage = sampling_percentage
	if not mask_map is None: sim.inputs.fixed_image_mask = mask_map
	if not mask_gene is None: sim.inputs.moving_image_mask = mask_gene
	sim_res = sim.run()
	return sim_res.outputs.similarity

def measure_similarity_geneexpression(stat_map,path_to_genes="/home/gentoo/src/abi2dsurqec_geneexpression/ABI_geneexpression_data",metric = 'MI',radius_or_number_of_bins = 64,comparison = 'gene',strategy='max'):
	"""
	master blabla
	"""
	#TODO: create a mask for the stat map? or userprovided? or both possible? Or use a single mask
	#TODO case: no exp average is wanted
	mask_map = create_mask(stat_map)
	results = dict()
	#loop through all gene folders, either get data form single experiment or get combined data.
	if comparison == 'gene':
		for dir in os.listdir(path_to_genes):
			path = os.path.join(path_to_genes,dir)
			print(dir)
			if not os.path.isdir(path):continue
			#print(path)
			#multiple experiment data available per gene
			if len(os.listdir(path)) > 1:
					img = []
					imgs = glob.glob(path + '/*/*_2dsurqec.nii.gz')
					for img_gene in imgs:
						if "sagittal" in img_gene:
							 img.append(mirror_sagittal(img_gene))
						else:
							img.append(img_gene)
					img_gene = create_experiment_average(img,strategy=strategy)
			elif len(os.listdir(path)) == 1:
					img_gene = glob.glob(path + '/*/*_2dsurqec.nii.gz')[0]
					if "sagittal" in img_gene:
							img_gene = mirror_sagittal(img_gene)
			else:
				#TODO: wrong place. Insert above in globglob. Generally bad idea with len(list.dir) == 1 If user creates a folder inside, program will crash... Used criteria after globglob
				print("Folder empty or no registered niftifile found. Or weird folder name :)")
				print("Skipping " + dir)
				continue
			#TODO: catch unexpected errors as to not interrupt program, print genename
			mask_gene = create_mask(img_gene)
			similarity = ants_measure_similarity(stat_map,img_gene,mask_gene = mask_gene,mask_map=mask_map,metric=metric,radius_or_number_of_bins=radius_or_number_of_bins)
			results[dir] = similarity

	elif comparison == 'experiment':
		for dir in os.listdir(path_to_genes):
			path = os.path.join(path_to_genes,dir)
			if not os.path.isdir(path):continue

			img = []
			imgs = glob.glob(path + '/*/*_2dsurqec.nii.gz')
			for img_gene in imgs:
				if "sagittal" in img_gene: img_gene = mirror_sagittal(img_gene)
				mask_gene = create_mask(img_gene)
				experiment_id = os.path.basename(img_gene).split("_")[3]
				print(experiment_id)
				id = dir + "_" + experiment_id
				similarity = ants_measure_similarity(stat_map,img_gene,mask_gene = mask_gene,mask_map=mask_map,metric=metric,radius_or_number_of_bins=radius_or_number_of_bins)
				results[id] = similarity

	name = metric + "_" + str(radius_or_number_of_bins) +".csv"
	with open(name,'w') as f:
		for key in results.keys():
			f.write("%s,%s\n"%(key,results[key]))



def main():


	img = "/home/gentoo/src/abi2dsurqec_geneexpression/ABI_geneexpression_data/Mef2c/Mef2c_P56_sagittal_79677145_200um/Mef2c_P56_sagittal_79677145_200um_2dsurqec.nii.gz"
# res = ants_measure_similarity("/home/gentoo/src/abi2dsurqec_geneexpression/ABI_geneexpression_data/Mef2c/Mef2c_P56_sagittal_79677145_200um/Mef2c_P56_sagittal_79677145_200um_2dsurqec.nii.gz","/home/gentoo/src/abi2dsurqec_geneexpression/ABI_geneexpression_data/Mef2c/Mef2c_P56_coronal_79567505_200um//Mef2c_P56_coronal_79567505_200um_2dsurqec.nii.gz")
	#	print(res)
	measure_similarity_geneexpression("/home/gentoo/src/abi2dsurqec_geneexpression/ABI_geneexpression_data/Mef2c/Mef2c_P56_coronal_79567505_200um/Mef2c_P56_coronal_79567505_200um_2dsurqec.nii.gz",metric = 'MI',radius_or_number_of_bins = 64,comparison = 'experiment')

#	measure_similarity_geneexpression("/home/gentoo/src/abi2dsurqec_geneexpression/ABI_geneexpression_data/Mef2c/Mef2c_P56_coronal_79567505_200um/Mef2c_P56_coronal_79567505_200um_2dsurqec.nii.gz",metric = 'CC',radius_or_number_of_bins = 4)

#	measure_similarity_geneexpression("/home/gentoo/src/abi2dsurqec_geneexpression/ABI_geneexpression_data/Mef2c/Mef2c_P56_coronal_79567505_200um/Mef2c_P56_coronal_79567505_200um_2dsurqec.nii.gz",metric = 'Mattes',radius_or_number_of_bins = 32)

#	measure_similarity_geneexpression("/home/gentoo/src/abi2dsurqec_geneexpression/ABI_geneexpression_data/Mef2c/Mef2c_P56_coronal_79567505_200um/Mef2c_P56_coronal_79567505_200um_2dsurqec.nii.gz",metric = 'MeanSquares',radius_or_number_of_bins = 64)

#	measure_similarity_geneexpression("/home/gentoo/src/abi2dsurqec_geneexpression/ABI_geneexpression_data/Mef2c/Mef2c_P56_coronal_79567505_200um/Mef2c_P56_coronal_79567505_200um_2dsurqec.nii.gz")
#	create_mask(img)

if __name__ == "__main__":
				main()

