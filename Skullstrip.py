#overlay 1 1 /data2/Projects/Maryam/QC/CPAC_run/output/pipeline_cpac_default__freq-filter__scrub/sub-5359869_ses-1/anatomical_reorient/sub-5359869_acq-HCP_T1w_resample.nii.gz -A /data2/Projects/Maryam/QC/CPAC_run/output/pipeline_cpac_default__freq-filter__scrub/sub-5359869_ses-1/anatomical_brain/sub-5359869_acq-HCP_T1w_resample_calc.nii.gz 1 10 /data2/Projects/Maryam/QC/CPAC_run/output/pipeline_cpac_default__freq-filter__scrub/sub-5359869_ses-1/anatomical_brain/mask111111.nii.gz

#slices ./data2/Projects/Maryam/QC/CPAC_run/output/pipeline_cpac_default__freq-filter__scrub/sub-5359869_ses-1/anatomical_brain/mask.nii.gz -o ../data2/Projects/Maryam/QC/CPAC_run/output/pipeline_cpac_default__freq-filter__scrub/sub-5359869_ses-1/anatomical_brain/sub-5359869_acq-HCP_T1w_resample_calc_mask_2_slices_4.png

######################Note############################
# before run this program you need to go inside Docker and change the permission of the C-PAC output:
#docker run -ti -v /data2/Projects/Maryam/QC/CPAC_run/:/bids_dataset --entrypoint=/bin/bash childmind/c-pac:latest
# cd bids_dataset/output
#chmod -R 777 /bids_dataset/output
######################Note############################



import os
import pandas as pd

template='/data2/Projects/Maryam/QC/QC_dashboard/template.html'
#this the root directory path to the output of your C-PAC run
datain = '/data2/Projects/Maryam/QC/CPAC_run/output/pipeline_cpac_default/'
#output_dir = '/data2/Projects/Maryam/QC/QC_dashboard/Pages/2/'
csvf = '/data2/Projects/Maryam/QC/QC_dashboard/sublist.csv'

sub_df = pd.read_csv(csvf)


#f =  open( '/data2/Projects/Maryam/QC/QC_dashboard/overley.sh','w')
for row in range(len(sub_df)):
	try:
		sub_id = sub_df.iloc[row,0]
		command =  'cd ' +datain + sub_id + '_ses-1' 
		(os.system(command))

		command2 = 'overlay 1 1 /data2/Projects/Maryam/QC/CPAC_run/output/pipeline_cpac_default__freq-filter__scrub/'+sub_id+'_ses-1/anatomical_reorient/'+sub_id+ '_acq-HCP_T1w_resample.nii.gz -A /data2/Projects/Maryam/QC/CPAC_run/output/pipeline_cpac_default__freq-filter__scrub/' +sub_id+ '_ses-1/anatomical_brain/' + sub_id+'_acq-HCP_T1w_resample_calc.nii.gz 1 10 /data2/Projects/Maryam/QC/CPAC_run/output/pipeline_cpac_default__freq-filter__scrub/'+sub_id+'_ses-1/anatomical_brain/mask001.nii.gz'
 		command3 = 'slices /data2/Projects/Maryam/QC/CPAC_run/output/pipeline_cpac_default__freq-filter__scrub/'+sub_id+'_ses-1/anatomical_brain/mask001.nii.gz -o /data2/Projects/Maryam/QC/CPAC_run/output/pipeline_cpac_default__freq-filter__scrub/'+sub_id+'_ses-1/anatomical_brain/'+sub_id+'_acq-HCP_T1w_resample_calc_mask_2_slices_001.png'
		os.system(command2)
		os.system(command3)

 		command5 = 'slices /data2/Projects/Maryam/QC/CPAC_run/output/pipeline_cpac_default__freq-filter__scrub/'+sub_id+'_ses-1/anatomical_reorient/'+sub_id+ '_acq-HCP_T1w_resample.nii.gz -o /data2/Projects/Maryam/QC/CPAC_run/output/pipeline_cpac_default__freq-filter__scrub/'+sub_id+'_ses-1/anatomical_brain/'+sub_id+'_acq-HCP_T1w_resample_calc_slices_002.png'
		os.system(command5)

		#f.write("%s\n" %command2) 
		#f.write("%s\n" %command3) 

	except:
		continue


#f.close
