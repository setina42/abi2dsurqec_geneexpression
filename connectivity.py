import argparse
import copy
import json
import os
import sys
import urllib
import urllib.request
import zipfile
import numpy
import nibabel
import re
import nrrd
from collections import defaultdict
from mhd_utils_3d import *
from nipype.interfaces.ants import ApplyTransforms
from nipype.interfaces.ants.base import ANTSCommand, ANTSCommandInputSpec
from nipype.interfaces.base import BaseInterface, BaseInterfaceInputSpec, traits, File, Str, TraitedSpec, Directory, CommandLineInputSpec, CommandLine, InputMultiPath, isdefined, Bunch, OutputMultiPath

#from nipype.interfaces.fsl import fslorient

API_SERVER = "http://api.brain-map.org/"
API_DATA_PATH = API_SERVER + "api/v2/data/"

#TODO: maybe merge into a single file that downloads either connectivtiy or gene expression data

def GetExpID():
    """
    Queries the Allen Mouse Brain Institute website for all gene expression data available for download.

    Returns:
    --------
    GeneNames: list[dict()]
        list of all genes where expression data is available for download. Dict contains experiment/gene metadata.

    SectionDataSetID : list(int)
        corresponding SectionDataSetID (SectionDataSet: see "http://help.brain-map.org/display/api/Data+Model")
        ID needed to specify download target.

    """

    startRow = 0
    numRows = 3
    totalRows = 3
    rows = []
    GeneNames = []
    SectionDataSetID = []
#    info = defaultdict(list)
    info = list()
#    GeneNames_cor = []
#    SectionDataSetID_cor = []

#    GeneNames_sag = []
#    SectionDataSetID_sag = []
    done = False

    while not done:
        #pagedUrl = API_DATA_PATH +"query.json?criteria=model::SectionDataSet,rma::criteria,products[abbreviation$eq'Mouse'],rma::include,specimen(stereotaxic_injections(primary_injection_structure,structures))startRow=%d&numRows=%d" % (startRow,numRows)
        r = "&start_row=%d&num_rows=%d" % (startRow,numRows)
        pagedUrl = API_DATA_PATH + "query.json?criteria=model::SectionDataSet,rma::criteria,products%5Bid$eq5%5D,rma::include,specimen(stereotaxic_injections(primary_injection_structure,structures))" + r
        print(pagedUrl)
        source = urllib.request.urlopen(pagedUrl).read()
        response = json.loads(source)
        rows += response['msg']
        for x in response['msg']:

            if x['failed'] == False :
                print(x['id'])
                info.append(x['id'])
        if totalRows < 0:
            totalRows = int(response['total_rows'])

        startRow += len(response['msg'])

        if startRow >= totalRows:
            done = True

    return info


def nrrd_to_nifti(file):
    print("Reading " + file)
    readnrrd = nrrd.read(file)
    data = readnrrd[0]
    header = readnrrd[1]
    print("Converting " + file)

    affine_matrix = numpy.array(header["space directions"],dtype=numpy.float)
    affine_matrix = affine_matrix*0.001
    affine_matrix = numpy.insert(affine_matrix,3,[0,0,0], axis=1)
    affine_matrix = numpy.insert(affine_matrix,3,[0,0,0,1], axis=0)

    #Change Orientation from PIR to RAS. Steps: PIR -> RIP -> RPI -> RPS -> RAS
    data.setflags(write=1)
    data = numpy.swapaxes(data,0,2)
    data = numpy.swapaxes(data,1,2)
    data = data[:,:,::-1]
    data = data[:,::-1,:]
    data = data[::-1,:,:]  #TODO: Check for Atlas files!!!!!!
    img = nibabel.Nifti1Image(data,affine_matrix)
    nii_path = os.path.join(os.path.dirname(file), os.path.basename(file).split(".")[0] + '.nii')
    nibabel.save(img,nii_path)

    return nii_path
#TODO: separate python script??
class ResampleImageSpec(ANTSCommandInputSpec):
    dimension = traits.Enum(2,3,4,
            argstr='%d',
            position=1)

    input_image = File(
            argstr='%s',
            mandatory=True,
            position=2)

    output_image = File(
            argstr='%s',
            mandatory = True,
            position=3)
    
    M = traits.Str(
            argstr ='%s',
            mandatory = True,
            position=4)
    size = traits.Int(
            argstr = 'size=%d',
            position=5)

    spacing = traits.Int(
            argstr = 'spacing=%d',
            position=6)
    
    interpolation = traits.Enum(0,1,2,3,4,
            argstr='%d',
            position = 7)
    #TODO:optional parameters size, spacing, interpolation, blablabla...

def ants_int_resample(dim,input_image,resolution,interpolation,output_image = None):
    print(str(dim),input_image,output_image,str(resolution),str(interpolation))
    ri = ResampleImage()
    ri.inputs.dimension = dim
    ri.inputs.input_image = input_image
    if output_image is None:
        #TODO: theres got to be an easier way...
        output_image = ""
        for s in input_image.split("_"):
            if "um" not in s :  
                output_image += s + "_"
            else:
                output_image += (str(resolution) + "um_")
        output_image = output_image[:-1]
    print(output_image)
    ri.inputs.output_image = output_image
    resolution = resolution / float(1000)
    resolution_cmd = str(resolution) + 'x' + str(resolution) + 'x' + str(resolution)
    print(resolution_cmd)
    ri.inputs.M = resolution_cmd
    ri.inputs.size = 1
    ri.inputs.spacing = 0
    ri.inputs.interpolation = interpolation
    print(ri.cmdline)
    ri.run()


class ResampleImageOutputSpec(TraitedSpec):
    output_image = File(exists=True,desc='Resampled Image')

class ResampleImage(ANTSCommand):
    _cmd = 'ResampleImage'
    input_spec = ResampleImageSpec
    output_spec = ResampleImageOutputSpec

def ants_resampleImage(img,resolution):
    at = ApplyTransforms()
    at.inputs.input_image = img
    at.inputs.reference_image = "bl"
    at.inputs.output_image="bl"
    at.inputs.transforms = 'identity'

#TODO:remove unnecessary files
def download_all_connectivity(info):
    """
    Download all given genes corresponding to SectionDataSetID given in 100um and 25um resolution, converts nrrd to nii, registers to dsurqec... and resamples files to 40 and 200 respectively.

    Parameters:
    -----------
        SectionDataSetID : list(int)
            o=[0.200000002980232 0 0 -6.26999998092651; 0 0.200000002980232 0 -10.6000003814697; 0 0 0.200000002980232 -7.88000011444092; 0 0 0 1]list of SectionDataSetID to download.
    """
    os.mkdir("/home/gentoo/src/abi2dsurqec_geneexpression/ABI_connectivity_data")
    download_url = "http://api.brain-map.org/grid_data/download_file/"
    for exp in info:
        #replace brackets with '_' and remove all other special characters
        path_to_exp = os.path.join("/home/gentoo/src/abi2dsurqec_geneexpression/ABI_connectivity_data",str(exp))
        os.mkdir(path_to_exp)
        for resolution in [100,25]:
            resolution_url = "?image=projection_density&resolution=" + str(resolution)
            url = download_url + str(exp) + resolution_url
            print(url)
            fh = urllib.request.urlretrieve(url)
            filename = str.split((fh[1]._headers[6][1]),'filename=')[1]  #TODO: Consistent??
            #TODO: do that differenttly ...
            filename = str.split(filename,";")[0]
            file_path = os.path.join(path_to_exp,filename)
            print(filename)
            print(file_path)
            os.rename(fh[0],file_path)
            #file_path_res=ants_resampleImage(file_path,resolution)
            file_path = nrrd_to_nifti(file_path)
            file_path = apply_composite(file_path)
            if resolution == 25 :
                target_resolution = 40
            elif resolution == 100:
                target_resolution = 200
            print("to resample")
            print(file_path)
            file_path = ants_int_resample(3,file_path,target_resolution,4)
def apply_composite(file):
    """
    Uses ANTS ApplyTransforms to register image to

    Parameters :
    ------------

    file : str
        path to image

    """
    at = ApplyTransforms()
    at.inputs.dimension = 3
    at.inputs.input_image = file
    at.inputs.reference_image = 'dsurqec_200micron_masked.nii'
    name = str.split(os.path.basename(file),'.nii')[0] + '_2dsurqec.nii.gz'
    at.inputs.interpolation = 'NearestNeighbor' #TODO: Sure??
    output_image = os.path.join(os.path.dirname(file),name)
    at.inputs.output_image = output_image
    at.inputs.transforms = 'abi2dsurqec_Composite.h5'
    at.run()

    #TODO sform to qform
    return output_image


def main():

    #ants_int_resample(3,"/home/gentoo/src/abi2dsurqec_geneexpression/dsurqec_200micron_masked.nii", "/home/gentoo/src/abi2dsurqec_geneexpression/dsurqec_100micron_masked.nii",100,4)
    info=GetExpID()
    download_all_connectivity(info)


if __name__ == "__main__":
    main()
