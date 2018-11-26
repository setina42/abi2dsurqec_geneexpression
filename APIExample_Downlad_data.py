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
from mhd_utils_3d import *
from nipype.interfaces.ants import ApplyTransforms

API_SERVER = "http://api.brain-map.org/"
API_DATA_PATH = API_SERVER + "api/v2/data/"

def GetGeneNames():   
    startRow = 0
    numRows = 2000
    totalRows = 2000  #was -1
    rows = []
    GeneNames = []
    SectionDataSetID = []

    GeneNames_cor = []
    SectionDataSetID_cor = []

    GeneNames_sag = []
    SectionDataSetID_sag = []
    done = False

    while not done:
        pagedUrl = API_DATA_PATH +"query.json?criteria=model::SectionDataSet,rma::criteria,[failed$eqfalse],products[abbreviation$eq'Mouse'],treatments[name$eq'ISH'],rma::include,genes,specimen(donor(age)),plane_of_section" + '&startRow=%d&numRows=%d' % (startRow,numRows)

        print(pagedUrl)
        source = urllib.request.urlopen(pagedUrl).read()
        response = json.loads(source)
        rows += response['msg']
        for x in response['msg']:

            if x['failed'] == False and x['expression'] == True :
                GeneNames.append(x['genes'])
                SectionDataSetID.append(x['id'])
                if x['reference_space_id'] == 10:
                    GeneNames_sag.append(x['genes'])
                    SectionDataSetID_sag.append(x['id'])
                if x['reference_space_id'] == 9:
                    GeneNames_cor.append(x['genes'])
                    SectionDataSetID_cor.append(x['id'])
        if totalRows < 0:
            totalRows = int(response['total_rows'])

        startRow += len(response['msg'])

        if startRow >= totalRows:
            done = True

    return GeneNames,SectionDataSetID

def download_all_ISH(GeneNames):
    download_url = "http://api.brain-map.org/grid_data/download/"
    for gene in GeneNames:
            url = download_url + str(gene)
            fh = urllib.request.urlretrieve(url)
            zf = zipfile.ZipFile(fh[0]) 
            zf.extractall(os.path.join("/home/gentoo/src/abi2dsurqec_geneexpression",os.path.basename(fh[0])))   #Somewhere along the line the file name is lost... :/
            zf.close()
            p = os.path.join("/home/gentoo/src/abi2dsurqec_geneexpression",os.path.basename(fh[0]))
            pe = os.path.join(p,"energy.mhd")
            convert_raw_to_nii(pe)
            pn = os.path.join(p,"energy.nii.gz")
            apply_composite(pn)

def convert_raw_to_nii(file):
    path = os.path.abspath('.')
    image_array, meta_header = load_raw_data_with_mhd(file)

    #Read header infomormation and create affine matrix
    dims = numpy.array(meta_header["ElementSpacing"].split(" "),dtype=numpy.float)
    affine_matrix = numpy.zeros((4,4),dtype=numpy.float)
    affine_matrix[0,0] = dims[0]
    affine_matrix[1,1] = dims[1]
    affine_matrix[2,2] = dims[2]

    image_array = numpy.swapaxes(image_array,1,2)

    image_array = image_array[::-1,:,:]
    image_array = image_array[:,:,::-1]
    image_array = image_array[:,::-1,:]

    #Bring to the right units?
    affine_matrix = affine_matrix*0.001

    affine_matrix[3,3] = 1  #?????

    img = nibabel.Nifti1Image(image_array,affine_matrix)
    nibabel.save(img,"energy.nii.gz")

def apply_composite(file):
    at = AntsApplyTransform()
    at.inputs.dimension = 3
    at.inputs.input_image = file
    at.inputs.reference_image = 'dsurqec_200micron_masked.nii'
    at.inputs.output_image = os.path.join(os.path.dirname(file),'energy2dsurqec.nii.gz')
    at.inputs.transforms = 'abi2dsurqec_Composite.h5'
    at.run()


def main():

    GeneNames=GetGeneNames()
    download_all_ISH(GeneNames[1][1:100])


if __name__ == "__main__":
    main()
