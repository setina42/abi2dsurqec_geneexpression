import numpy as np
import nibabel
import os
import glob
from nipype.interfaces import ants as ants
import nipype.interfaces.fsl.maths as fsl
import csv
import argparse
#TODO:do I need nilearn??
import nilearn
import matplotlib.pyplot as plt
import samri.plotting.maps as maps
from collections import defaultdict
from nilearn._utils.extmath import fast_abs_percentile

#TODO: Currently loading data 3 times with nibabel.laod and get_fdata(), maybe load once and pass to func (if its the same data!)

def transform(x,y,z,affine):
	M = affine[:4, :4]
	A = np.linalg.inv(M)
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
#TODO: checkec for case 3, test for other cases as well
	#replace
	if np.shape(left_side)[0] > np.shape(img_data[0:mid,:,:])[0]:
#		print("case 1")
		#case 1: origin slightly to the left (or right??), need to trim left_side to the size of the right side
		replace_value = np.shape(left_side)[0] - np.shape(img_data[0:mid,:,:])[0]
		img_data[0:mid,:,:][img_data[0:mid,:,:]  == -1] = left_side[(replace_value-1):,:,:][img_data[0:mid,:,:]  == -1]
		#np.savetxt("Slice_case1.txt",img_data[:,(int(origin[1]) -5):(int(origin[1])) -2, (int(origin[2]) -5)])
		img_data[mid:np.shape(right_side[0]),:,:][img_data[mid:,:,:]  == -1] = right_side[:,:,:][img_data[mid:,:,:]  == -1]

	elif np.shape(left_side)[0] < np.shape(img_data[0:mid,:,:])[0]:
#		print("case 2")
		#case 2 : origin slightly to the right (or left??), need to
		replace_value = np.shape(img_data[0:mid,:,:])[0] - np.shape(left_side)[0]
		img_data[replace_value:mid,:,:] = left_side
		#np.savetxt("Slice_case2.txt",img_data[:,(int(origin[1]) -5):(int(origin[1])) -2, (int(origin[2]) -5)])

		img_data[mid:,:,:][img_data[mid:,:,:]  == -1] = right_side[0:np.shape(img_data[mid:,:,:]),:,:][img_data[mid:,:,:]  == -1]
	else:
#		print("case 3")
		#case 3: same size
		#TODO: -1 or 0??
		#TODO: There's got to be a better way to write this....  -> a_slice notation:: a[10:16] is a reference, not a copy!
		img_data[0:mid,:,:][img_data[0:mid,:,:]  == -1] = left_side[img_data[0:mid,:,:]  == -1]
		#np.savetxt("Slice_case3.txt",img_data[:,(int(origin[1]) -5):(int(origin[1])) -2, (int(origin[2]) -5)])
		img_data[mid+1:,:,:][img_data[mid+1:,:,:]  == -1] = right_side[[img_data[mid+1:,:,:]  == -1]]

	img_average = nibabel.Nifti1Image(img_data,img.affine)
	filename = str.split(os.path.basename(image),'.nii')[0] + '_mirrored.nii.gz'
	path_to_mirrored = os.path.join(os.path.dirname(image),filename)
	nibabel.save(img_average,path_to_mirrored)

	return path_to_mirrored

def create_mask(image,threshold):
	#I think i need 0 to be in my mask. This seems not to be possible using fslmaths, so maybe do directly with numpy? thr sets all to zero below the value and bin uses image>0 to binarise.
	#mask = fsl.Threshold()
	#mask.inputs.thresh = 0
	#mask.inputs.args = '-bin'
	#mask.inputs.in_file = image
	#img_out = str.split(image,'.nii')[0] + '_mask_fsl.nii.gz'
	#mask.inputs.out_file = img_out
	#mask.run()
	
	mask_img = nibabel.load("/home/gentoo/src/abi2dsurqec_geneexpression/dsurqec_200micron_mask.nii")
	atlas_mask = mask_img.get_fdata()

	#using numpy instead of fslmaths
	img = nibabel.load(image)
	img_data = img.get_fdata()
	#apparently ambibuous:img_data[img_data >= 0 and atlas_mask == 1] = 1
#	print("mask")
#	print(np.min(img_data))
	#TODO: Atlas mask needs to have the right resolution, load different one for dsurqec_40micron_mask.nii. Do I always compare between the same resolution? Otherwise I need two brain masks...
	#TODO: ensure shape mathces
	img_data[np.logical_and(img_data > threshold,atlas_mask == 1)] = 1
	img_data[np.logical_or(img_data <= threshold,atlas_mask==0)] = 0
	img_out = str.split(image,'.nii')[0] + '_mask.nii.gz'
	img_mask = nibabel.Nifti1Image(img_data,img.affine)
	nibabel.save(img_mask,img_out)

	return img_out

def nan_if(arr,value):
	return np.where(arr == value, np.nan,arr)

#TODO:evaluate mean average, diff to and between single experiments, especially now that we include expression = false as well, with small expr, patterns. Also , habe a loo
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
		for i in range(0,len(img_data)):
			#ignore values of -1 for mean
			img_data[i] = nan_if(img_data[i],-1)

		average_img = np.nanmean([*img_data],axis=0)
		average_img[np.isnan(average_img)] = -1

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
	try:
		sim_res = sim.run()
		res = sim_res.outputs.similarity
	except:
		print("something happened?")
		res = 0
	return res


def _plot(dis):
	fig = plt.figure()
	main = fig.add_axes([0,0,0.6,1.0])
	ax_1=fig.add_axes([0.6,0,0.4,0.3])
	ax_2=fig.add_axes([0.6,0.3,0.4,0.3])
	ax_3=fig.add_axes([0.6,0.6,0.4,0.3])
	
	for ax in [main,ax_1,ax_2,ax_3]:
		ax.get_xaxis().set_visible(False)
		ax.get_yaxis().set_visible(False)
		ax.patch.set_alpha(0.1)
		ax.spines["top"].set_visible(False)
		ax.spines["bottom"].set_visible(False)
		ax.spines["right"].set_visible(False)
		ax.spines["left"].set_visible(False)

	main.imshow(plt.imread("_stat.png"))
	ax_1.imshow(plt.imread("0.png"))
	ax_2.imshow(plt.imread("1.png"))
	ax_3.imshow(plt.imread("2.png"))
	plt.savefig("all.png",dpi=600)


#TODO: maybe use the same cut coords for all plots, may prove difficult bc not possbile to use nilearns func directly
def plot_results(stat_map,results,hits = 3, template = "/usr/share/mouse-brain-atlases/ambmc2dsurqec_15micron_masked.obj",comparison='gene',path_to_genes="/home/gentoo/src/abi2dsurqec_geneexpression/ABI_geneexpression_data"):
	# TODO: put into stat3D or stat, to avoid loading the data twice threshold = fast_abs_percentile(stat_map)
	dis = dict()
	img_s = nibabel.load(stat_map)
	img_data_s = img_s.get_fdata()
	tresh_s = fast_abs_percentile(img_data_s[img_data_s >0],percentile=98)
	print(tresh_s)
	display_stat = maps.stat3D(stat_map,template="/home/gentoo/src/abi2dsurqec_geneexpression/dsurqec_200micron_masked.nii",save_as= '_stat.png',threshold=tresh_s,pos_values=True,figure_title="stat")
	dis["main"] = display_stat
	for i in range(0,hits):
		gene_name = results[i][0].split("_")[0] #TODO:this should work in both cases (also for connectivity??)
		full_path_to_gene = results[i][1][1]
		print("now plotting: ")
		print(full_path_to_gene)
		img = nibabel.load(full_path_to_gene)
		img_data = img.get_fdata()
		tresh = fast_abs_percentile(img_data[img_data > 0],percentile=98)
		display = maps.stat3D(full_path_to_gene,template="/home/gentoo/src/abi2dsurqec_geneexpression/dsurqec_200micron_masked.nii",save_as=str(i) + '.png',threshold=tresh,pos_values=True,figure_title=gene_name)
		dis[str(i)] = display
	_plot(dis)
	print(tresh_s)
	print(tresh)
	#TODO:sep.function?



	return

def output_results(results,hits = 3,output_name=None):

	print("Top " + str(hits) + " hits: ")
	for i in range(0,hits):
		print(str(results[i][1][0]) + " " + str(results[i][1][1]))
	#TODO: some smart name...
	#save to csv
	if output_name is None:
		output_name =  "output_results.csv"
	else:
		output_name = output_name + ".csv"

	with open(output_name,'w') as f:
		for i in range(0,len(results)):
			f.write("%s,%s\n"%(results[i][0],results[i][1]))


	return
#TODO: sorted results: same score, sorting?? warning
def measure_similarity_geneexpression(stat_map,path_to_genes="/home/gentoo/src/abi2dsurqec_geneexpression/ABI_geneexpression_data",metric = 'MI',radius_or_number_of_bins = 64,comparison = 'gene',strategy='mean'):
	"""
	master blabla
	"""

	#TODO: if mirrored or mask files are already present, don't make them again
	#TODO:
	#TODO: create a mask for the stat map? or userprovided? or both possible? Or use a single mask. Also, if yes, include threshold in mask func
	mask_map = create_mask(stat_map,0)
	#results = dict()
	results = defaultdict(list)
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
			mask_gene = create_mask(img_gene,-1)
			similarity = ants_measure_similarity(stat_map,img_gene,mask_gene = mask_gene,mask_map=mask_map,metric=metric,radius_or_number_of_bins=radius_or_number_of_bins)
			print(dir)
			print(similarity)
			results[dir].append(similarity)
			results[dir].append(img_gene)

	elif comparison == 'experiment':
		for dir in os.listdir(path_to_genes):
			path = os.path.join(path_to_genes,dir)
			if not os.path.isdir(path):continue

			img = []
			imgs = glob.glob(path + '/*/*_2dsurqec.nii.gz')
			for img_gene in imgs:
				if "sagittal" in img_gene: img_gene = mirror_sagittal(img_gene)
				mask_gene = create_mask(img_gene,-1)
				experiment_id = os.path.basename(img_gene).split("_")[3]
				id = dir + "_" + experiment_id
				similarity = ants_measure_similarity(stat_map,img_gene,mask_gene = mask_gene,mask_map=mask_map,metric=metric,radius_or_number_of_bins=radius_or_number_of_bins)
				results[id].append(similarity)
				results[id].append(img_gene)
	#TODO: sort, or use sorted dict form beginning, or only keep top scores anyway?
	#print(comparison)
	#view = [ (v,k) for k,v in results.items() ]
	#view.sort(reverse=True)
	#for v,k in view:
	#	print( "%s: %f" % (k,v))

	#TODO: if metric = MSE, sort other way round
	#sort results:
	sorted_results = sorted(results.items(),key=lambda x: x[1][0])
#	print(sorted_results[4]) #5th highest result
#	print(sorted_results[4][0]) #key = gene or exp_id
#	print(sorted_results[4][1][0]) #similarity score
#	print(sorted_results[4][1][1]) #path

	output_results(sorted_results, hits = 3)
	plot_results(stat_map,sorted_results,hits=3)

	return results


def measure_similarity_connectivity(stat_map,path_to_exp="/home/gentoo/src/abi2dsurqec_geneexpression/ABI_connectivity_data",metric = 'MI',radius_or_number_of_bins = 64,resolution=200):
	#TODO: mirror sagittal for connectivity??
	mask_map = create_mask(stat_map,0)
	results = defaultdict(list)
	for dir in os.listdir(path_to_exp):
		path = os.path.join(path_to_exp,dir)
		print(path)
		img = glob.glob(path + '/*' + str(resolution) + '*_2dsurqec.nii*')[0]
		print(img)
		print(len(glob.glob(path + '/*' + str(resolution) + '*_2dsurqec.nii*')))
		mask_gene = create_mask(img,0)
		similarity = ants_measure_similarity(stat_map,img,mask_gene = mask_gene,mask_map =mask_map,metric=metric,radius_or_number_of_bins=radius_or_number_of_bins)
		results[dir].append(similarity)
		results[dir].append(img)  #path for plotting
	sorted_results = sorted(results.items(),key=lambda x: x[1][0])
	print(sorted_results)
	output_results(sorted_results,hits = 3)
	plot_results(stat_map,sorted_results,hits=3)

def main():

	parser = argparse.ArgumentParser(description="Similarity",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('--stat_map','-s',type=str)
	parser.add_argument('--path_genes','-p',type=str)
	parser.add_argument('--comparison','-c',type=str, default='gene')
	parser.add_argument('--radius_or_number_of_bins','-r',type=int,default = 64)
	parser.add_argument('--metric','-m',type=str,default='MI')
	parser.add_argument('--strategy','-y',type=str, default='max')
	args=parser.parse_args()
	img = "/home/gentoo/src/abi2dsurqec_geneexpression/ABI_geneexpression_data/Mef2c/Mef2c_P56_sagittal_79677145_200um/Mef2c_P56_sagittal_79677145_200um_2dsurqec.nii.gz"
# res = ants_measure_similarity("/home/gentoo/src/abi2dsurqec_geneexpression/ABI_geneexpression_data/Mef2c/Mef2c_P56_sagittal_79677145_200um/Mef2c_P56_sagittal_79677145_200um_2dsurqec.nii.gz","/home/gentoo/src/abi2dsurqec_geneexpression/ABI_geneexpression_data/Mef2c/Mef2c_P56_coronal_79567505_200um//Mef2c_P56_coronal_79567505_200um_2dsurqec.nii.gz")
	#	print(res)
#	measure_similarity_geneexpression("/home/gentoo/src/abi2dsurqec_geneexpression/save_small_dataset/Mef2c/Mef2c_P56_coronal_79567505_200um/Mef2c_P56_coronal_79567505_200um_2dsurqec.nii.gz",metric = 'MI',radius_or_number_of_bins = 64,comparison = 'gene')
#	measure_similarity_geneexpression("/home/gentoo/ABI_data_full/data/Tlx2/Tlx2_P56_sagittal_81655554_200um/Tlx2_P56_sagittal_81655554_200um_2dsurqec_mirrored.nii.gz",metric = 'MI',path_to_genes="/home/gentoo/ABI_data_full/data",radius_or_number_of_bins = 64,comparison = 'gene')

	measure_similarity_connectivity("/home/gentoo/src/abi2dsurqec_geneexpression/ABI_geneexpression_data/Kat6a/Kat6a_P56_sagittal_71764326_200um/Kat6a_P56_sagittal_71764326_200um_2dsurqec_mirrored.nii.gz",metric = 'MI',radius_or_number_of_bins = 64)



#	measure_similarity_geneexpression("/home/gentoo/src/abi2dsurqec_geneexpression/ABI_geneexpression_data/Kat6a/Kat6a_P56_sagittal_71764326_200um/Kat6a_P56_sagittal_71764326_200um_2dsurqec_mirrored.nii.gz",metric = 'MI',path_to_genes="/home/gentoo/ABI_data_full/data",radius_or_number_of_bins = 64,comparison = 'gene')



#	measure_similarity_geneexpression("/home/gentoo/src/abi2dsurqec_geneexpression/ABI_geneexpression_data/Mef2c/Mef2c_P56_coronal_79567505_200um/Mef2c_P56_coronal_79567505_200um_2dsurqec.nii.gz",metric = 'CC',radius_or_number_of_bins = 4)

#	measure_similarity_geneexpression("/home/gentoo/src/abi2dsurqec_geneexpression/ABI_geneexpression_data/Mef2c/Mef2c_P56_coronal_79567505_200um/Mef2c_P56_coronal_79567505_200um_2dsurqec.nii.gz",metric = 'Mattes',radius_or_number_of_bins = 32)

#	measure_similarity_geneexpression("/home/gentoo/src/abi2dsurqec_geneexpression/ABI_geneexpression_data/Mef2c/Mef2c_P56_coronal_79567505_200um/Mef2c_P56_coronal_79567505_200um_2dsurqec.nii.gz",metric = 'MeanSquares',radius_or_number_of_bins = 64)

#	measure_similarity_geneexpression("/home/gentoo/src/abi2dsurqec_geneexpression/ABI_geneexpression_data/Mef2c/Mef2c_P56_coronal_79567505_200um/Mef2c_P56_coronal_79567505_200um_2dsurqec.nii.gz")
#	create_mask(img)

if __name__ == "__main__":
				main()

