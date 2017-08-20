# # Step 1: create .addTheta files
# seq 0 99 | perl -ne'printf("%03d\n", $_)' | \
# 	parallel -j 24 "java -Xmx1g -cp stairway_plot_es/:stairway_plot_es/swarmops.jar  Stairway_fold_random_break5  out_ERY_100ninput/input/{}.unfolded.sfs.folded.stair 34"

# Step 2: determine number of break points
# mv -f out_ERY_100ninput/input/*addTheta out_ERY_100ninput/rand34/
java -Xmx1g -cp stairway_plot_es/ Stairpainter ERY.blueprint
bash ERY.blueprint.plot.sh
