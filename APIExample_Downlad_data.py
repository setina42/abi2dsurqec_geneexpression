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
from collections import defaultdict
from mhd_utils_3d import *
from nipype.interfaces.ants import ApplyTransforms

API_SERVER = "http://api.brain-map.org/"
API_DATA_PATH = API_SERVER + "api/v2/data/"

def GetGeneNames():
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
    numRows = 2000
    totalRows = -1
    rows = []
    GeneNames = []
    SectionDataSetID = []
    info = defaultdict(list)

#    GeneNames_cor = []
#    SectionDataSetID_cor = []

#    GeneNames_sag = []
#    SectionDataSetID_sag = []
    done = False

    while not done:
        pagedUrl = API_DATA_PATH +"query.json?criteria=model::SectionDataSet,rma::criteria,[failed$eqfalse],products[abbreviation$eq'Mouse'],treatments[name$eq'ISH'],rma::include,genes,specimen(donor(age)),plane_of_section" + '&startRow=%d&numRows=%d' % (startRow,numRows)

        print(pagedUrl)
        source = urllib.request.urlopen(pagedUrl).read()
        response = json.loads(source)
        rows += response['msg']
        for x in response['msg']:

            if x['failed'] == False and x['expression'] == True :
                info[x['genes'][0]['acronym']].append(x['id'])
        if totalRows < 0:
            totalRows = int(response['total_rows'])

        startRow += len(response['msg'])

        if startRow >= totalRows:
            done = True

    return info

def download_all_ISH(info):
    """
    Download all given genes corresponding to SectionDataSetID given, converts data format from mhd/raw to nii and registers data to dsurqec template.

    Parameters:
    -----------
        SectionDataSetID : list(int)
            list of SectionDataSetID to download.
    """
    os.mkdir("/home/gentoo/src/abi2dsurqec_geneexpression/ABI_geneexpression_data")
    download_url = "http://api.brain-map.org/grid_data/download/"
    for gene in info:
        path_to_gene = os.path.join("/home/gentoo/src/abi2dsurqec_geneexpression/ABI_geneexpression_data",gene)
        os.mkdir(path_to_gene)
        for id in info[gene]:
            url = download_url + str(id)
            fh = urllib.request.urlretrieve(url)
            zf = zipfile.ZipFile(fh[0])
            filename = str.split((fh[1]._headers[6][1]),'filename=')[1]  #TODO: Consistent???
            filename = str.split(filename,'.zip')[0]
            path_to_folder = os.path.join(path_to_gene,filename)
            zf.extractall(os.path.join(path_to_gene,filename))
            zf.close()
            path_to_mhd = os.path.join(path_to_folder,"energy.mhd")
            path_to_nifti = convert_raw_to_nii(path_to_mhd,filename)
            apply_composite(path_to_nifti)

def convert_raw_to_nii(input_file,output_file):
    """
    Converts mhd/raw format to nifti and orient data matrix in RAS-space

    Parameters:
    -----------
        input_file : str
            path to .mhd file
        output_file : str
            filename prefix (?)
    """
    path = os.path.abspath('.')
    image_array, meta_header = load_raw_data_with_mhd(input_file)

    #Read header infomormation and create affine matrix
    dims = numpy.array(meta_header["ElementSpacing"].split(" "),dtype=numpy.float)
    affine_matrix = numpy.zeros((4,4),dtype=numpy.float)
    affine_matrix[0,0] = dims[0]
    affine_matrix[1,1] = dims[1]
    affine_matrix[2,2] = dims[2]

    #Orient in RAS
    image_array = numpy.swapaxes(image_array,1,2)

    image_array = image_array[::-1,:,:]
    image_array = image_array[:,:,::-1]
    image_array = image_array[:,::-1,:]

    #Bring to the right units
    affine_matrix = affine_matrix*0.001

    affine_matrix[3,3] = 1

    img = nibabel.Nifti1Image(image_array,affine_matrix)
    name = output_file + '.nii.gz'
    nibabel.save(img,os.path.join(os.path.dirname(input_file),name))
    return os.path.join(os.path.dirname(input_file),name)

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
    at.inputs.output_image = os.path.join(os.path.dirname(file),name)
    at.inputs.transforms = 'abi2dsurqec_Composite.h5'
    at.run()


def main():

    info=GetGeneNames()
    download_all_ISH(info)


if __name__ == "__main__":
    main()
