

# simple script that, when run, grabs each figure from the top level directory
# (where it was presumably just made by running a script) and then pdf_crops
# it and puts it here.

# RUN THIS FROM WITHIN THE IMAGES DIRECTORY

## cwt_ppn_estimation.tex:
imgs="
	recovery_trees_divided.pdf
	recovery_histo_panel_1.pdf
	recovery_histo_panel_2.pdf
	actual_vs_pred_cwt_recoveries_p_marked.pdf
	post_mean_theta_v_counts_f_marked_times_p_marked.pdf
	post_mean_theta_v_counts_release_state.pdf
"

for i in $imgs; do 
	echo $i;
	pdfcrop ../../../$i  ./$i 
done