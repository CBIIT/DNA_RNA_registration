from xml.dom import minidom
import os
import pandas as pd
from skimage import io, feature, img_as_ubyte, img_as_uint
from scipy import ndimage
import numpy as np
import shutil
from matplotlib import pyplot as plt
import sys
def rna_dna_registration(moving_path, template_path, output_folder):
    
    if os.path.isdir(output_folder) == False:
        os.mkdir(output_folder)
    reg_dataset_path = os.path.join(output_folder, os.path.split(os.path.split(template_path)[0])[-1])

    shutil.copytree(template_path, reg_dataset_path )
    os.remove(os.path.join(reg_dataset_path, 'MeasurementData.mlf'))
    # parse an xml file by name
    metadatafilename = os.path.join(template_path, 'MeasurementData.mlf')
    mydoc = minidom.parse(metadatafilename)
    
    
    ### append moving image metadata
    aaa = mydoc.getElementsByTagName('bts:MeasurementData')[0]
    rna_metadatafilename = os.path.join(moving_path, 'MeasurementData.mlf')
    rna_mydoc = minidom.parse(rna_metadatafilename)
    rna_items = rna_mydoc.getElementsByTagName('bts:MeasurementRecord')
    for i in range(rna_items.length):
        # print(rna_items[i].attributes['bts:Ch'].value)
        if (rna_items[i].attributes['bts:Ch'].value==str(2)) | (rna_items[i].attributes['bts:Ch'].value==str(3)):
            aaa.appendChild(rna_items[i])


    PATH_TO_FILES = os.path.split(metadatafilename)[0]
    items = mydoc.getElementsByTagName('bts:MeasurementRecord')

    df_cols = ["ImageName", "column", "row", "time_point", "field_index", "z_slice", "channel", 
               "x_coordinates", "y_coordinates","z_coordinate", "action_index", "action", "Type", "Time"]
    rows = []

    for i in range(items.length):

        fullstring = items[i].firstChild.data
        substring = "Error"


        if fullstring.find(substring) == -1:
            if items[i].attributes['bts:Type'].value=='IMG':
                rows.append({

                     "ImageName": os.path.join(PATH_TO_FILES, items[i].firstChild.data), 
                     "column": items[i].attributes['bts:Column'].value, 
                     "row": items[i].attributes['bts:Row'].value, 
                     "time_point": items[i].attributes['bts:TimePoint'].value, 
                     "field_index": items[i].attributes['bts:FieldIndex'].value, 
                     "z_slice": items[i].attributes['bts:ZIndex'].value, 
                     "channel": items[i].attributes['bts:Ch'].value,
                     "x_coordinates": items[i].attributes['bts:X'].value,
                     "y_coordinates": items[i].attributes['bts:Y'].value,
                     "z_coordinate": items[i].attributes['bts:Z'].value,
                     "action_index": items[i].attributes['bts:ActionIndex'].value,
                     "action": items[i].attributes['bts:Action'].value, 
                     "Type": items[i].attributes['bts:Type'].value, 
                     "Time": items[i].attributes['bts:Time'].value,

                })

    out_df = pd.DataFrame(rows, columns = df_cols)
    with open(os.path.join(reg_dataset_path,'MeasurementData.mlf'),"w") as file_handle:
        mydoc.writexml(file_handle)
        file_handle.close()



    nuclei_df = out_df.loc[out_df["channel"]==str(1)]
    ch2_action_ind = "A0" + str(out_df.loc[out_df["channel"]==str(2)]["action_index"].iloc[0])
    ch3_action_ind = "A0" + str(out_df.loc[out_df["channel"]==str(3)]["action_index"].iloc[0])
    for i in range(len(nuclei_df)):
        dna_nuclei_img = nuclei_df["ImageName"].iloc[i]
        rna_nuclei_img = os.path.join(moving_path, os.path.split(dna_nuclei_img)[-1])
        rna_spot_name_ch2 = os.path.split(rna_nuclei_img)[-1].replace("C01","C02").replace("A01",ch2_action_ind)
        rna_spot_name_ch3 = os.path.split(rna_nuclei_img)[-1].replace("C01","C03").replace("A01",ch3_action_ind)
        rna_spot_img_ch2 = os.path.join(moving_path, rna_spot_name_ch2)
        rna_spot_img_ch3 = os.path.join(moving_path, rna_spot_name_ch3)
        if os.path.isfile(dna_nuclei_img):
            im1 = io.imread(dna_nuclei_img, as_gray=True)
            if os.path.isfile(rna_nuclei_img):
                im2 = io.imread(rna_nuclei_img, as_gray=True)
                rna_spot_ch2 = io.imread(rna_spot_img_ch2, as_gray=True)
                rna_spot_ch3 = io.imread(rna_spot_img_ch3, as_gray=True)
                sh_row, sh_col = im1.shape

                # Registration of the two images
                translation = feature.register_translation(im1, im2, upsample_factor=1)[0]
                im2_register = ndimage.shift(im2, translation)
                reg_rna_spot_ch2 = ndimage.shift(rna_spot_ch2, translation)
                reg_rna_spot_ch3 = ndimage.shift(rna_spot_ch3, translation)
                
                new_rna_spot_name_ch2 = os.path.split(rna_nuclei_img)[-1].replace("C01","C02").replace("A01",ch2_action_ind)
                new_rna_spot_name_ch3 = os.path.split(rna_nuclei_img)[-1].replace("C01","C03").replace("A01",ch3_action_ind)
                
                io.imsave(os.path.join(reg_dataset_path,new_rna_spot_name_ch2), reg_rna_spot_ch2, check_contrast=False)
                io.imsave(os.path.join(reg_dataset_path,new_rna_spot_name_ch3), reg_rna_spot_ch3, check_contrast=False)

                all1=np.zeros((sh_row, sh_col,3), dtype="uint16")
                all1[:,:,0]=im1
                all1[:,:,2]=im2_register

                all2=np.zeros((sh_row, sh_col,3), dtype="uint16")
                all2[:,:,0]=im1
                all2[:,:,2]=im2


                fig = plt.figure()
                fig.set_size_inches(28, 14)
                fig.subplots_adjust(hspace=0.4, wspace=0.4)
                ax = fig.add_subplot(1, 2, 1)
                ax.imshow(all1)
                offsets_val= "Row_offset=" + str(translation[0]) + "    Column_offset=" + str(translation[1])
                ax.set_title("Registered Nuclei\n" + offsets_val, fontsize=16)

                ax = fig.add_subplot(1, 2, 2)
                ax.imshow(all2)
                ax.set_title("Un-Registered Nuclei", fontsize=16)
                reg_nuclei_path = os.path.join(output_folder, "registered_nuclei")
                if os.path.isdir(reg_nuclei_path) == False:
                    os.mkdir(reg_nuclei_path)
                comaprison_img = os.path.join(reg_nuclei_path,os.path.split(dna_nuclei_img)[-1])
                fig.savefig(comaprison_img) 
                fig.clear()
    
    print("Registratoin of RNA channels has finished.")

moving_path = str(sys.argv[1])
template_path = str(sys.argv[2])
output_folder = str(sys.argv[3])

rna_dna_registration(moving_path, template_path, output_folder)

