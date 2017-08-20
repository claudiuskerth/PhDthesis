# # Step 1: create .addTheta files
# seq 0 199 | perl -ne'printf("%03d\n", $_)' | \
# 	parallel -j 12 "java -Xmx1g -cp stairway_plot_es/:stairway_plot_es/swarmops.jar Stairway_fold_random_break5 out_PAR_200_corr/input/{}.unfolded.sfs.dadi.corr.stair 34"

# # Step 2: determine number of break points
mv -f out_PAR_200_corr/input/*.addTheta out_PAR_200_corr/rand34/
java -Xmx1g -cp stairway_plot_es/ Stairpainter PAR_200_corr.blueprint


# bash PAR_200_corr.blueprint.plot.sh
