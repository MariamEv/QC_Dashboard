import os
import glob
import csv
import pandas as pd

def create_htmls(sub_id,next_sub,prev_sub):	 

		    #if  os.path.isdir(path): 
			#head,sub=os.path.split(path)
			#copy template into a new html file
			output_file = output_dir+sub_id+'.html'
			a= os.system('cp -r '+template+' '+output_file)
			print sub_id,next_sub,prev_sub	
			f =open(output_file,'r')
			filedata = f.read()
			f.close

			before = '<div class="row"><div class="column"><a id="prev" href="/data2/Projects/Maryam/QC/QC_dashboard/a/sub-12348888.html" class="button" >&laquo;Prev</a></div><div class="column"><h2 id="sub_id" style="text-align: center;" >Subject ID:</h2></div><div class="column_right"><a id="next" href="/data2/Projects/Maryam/QC/QC_dashboard/a/sub-12348888.html" class="button" >Next&raquo;</a></div></div>'

			after = '<div class="row"><div class="column"><a id="prev" href='+output_dir+prev_sub+'.html class="button" >&laquo;Prev</a></div><div class="column"><h2 id="sub_id" style="text-align: center;" >Subject ID:'+sub_id+'</h2></div><div class="column_right"><a id="next" href='+output_dir+next_sub+'.html class="button" >Next&raquo;</a></div></div>'

			skl_before = '<div class="column2"> <b>SkullStrip</b><img src="/data2/Projects/Maryam/QC/CPAC_run/output/pipeline_cpac_default__freq-filter__scrub/sub-5359869_ses-1/anatomical_brain/sub-5359869_acq-HCP_T1w_resample_calc_mask_slices.png" alt="" style="width:100%">'
			skl_after = '<img src="'+SkullStrip_path_to_general_pipline+sub_id+'_ses-1/anatomical_brain/'+sub_id+'_acq-HCP_T1w_resample_calc_mask_2_slices_001.png" alt="" style="width:100%">'
			
			print skl_after 
			skl_before_2 = '<div class="column2"> <b>SkullStrip</b><img src="/data2/Projects/Maryam/QC/CPAC_run/output/pipeline_cpac_default__freq-filter__scrub/sub-5359869_ses-1/anatomical_brain/sub-5359869_acq-HCP_T1w_resample_calc_slices.png" alt="" style="width:100%">'
			skl_after_2 = '<img src="'+SkullStrip_path_to_general_pipline+sub_id+'_ses-1/anatomical_brain/'+sub_id+'_acq-HCP_T1w_resample_calc_slices_002.png" alt="" style="width:100%">'

			fd_before = '/data2/Projects/Maryam/QC/FD_images'
			fd_after = FDs_Path + sub_id
			

			repls = ('/data2/Projects/Maryam/QC/QC_imgs/',datain) , (before,after) , ('sub-5359869',sub_id), (fd_before, fd_after) , (skl_before ,skl_after), (skl_before_2, skl_after_2)
			new_data = reduce(lambda a, kv: a.replace(*kv), repls, filedata)
			#new_data = filedata.replace(b , a)

			f = open(output_file,'w')
			f.write(new_data)
			f.close



template='/data2/Projects/Maryam/QC/QC_dashboard/template.html'
FDs_Path = '/data2/Projects/Maryam/QC/QC_dashboard/Pages/2/FDs/FD/'
SkullStrip_path_to_general_pipline  = '/data2/Projects/Maryam/QC/CPAC_run/output/pipeline_cpac_default__freq-filter__scrub/'
#this the root directory path to the output of your C-PAC run
datain='/data2/Projects/Maryam/QC/CPAC_run/output/'
output_dir = '/data2/Projects/Maryam/QC/QC_dashboard/Pages/2/'
csvf = '/data2/Projects/Maryam/QC/QC_dashboard/sublist.csv'

sub_df = pd.read_csv(csvf)

for row in range(len(sub_df)):
		sub_id = sub_df.iloc[row,0]
		print sub_id
		if(row == 0):
			next_sub=sub_df.iloc[row+1,0]
			prev_sub=sub_df.iloc[len(sub_df)-1,0]	
		elif( row == len(sub_df)-1):
			next_sub=sub_df.iloc[0,0]
			prev_sub=sub_df.iloc[row-1,0]
		else:
			next_sub=sub_df.iloc[row+1,0]
			prev_sub=sub_df.iloc[row-1,0]	
		create_htmls(sub_id,next_sub,prev_sub)




