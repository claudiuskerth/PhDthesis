# # Step 1: create .addTheta files
# seq 0 199 | perl -ne 'printf("%03d\n", $_)' | \
# 	parallel -j 12 "java -Xmx1g -cp stairway_plot_es/:stairway_plot_es/swarmops.jar Stairway_fold_random_break5 out_ERY_200/input/{}.unfolded.sfs.folded.stair 34"

# # Step 2: determine number of break points
mv -f out_ERY_200/input/*.addTheta out_ERY_200/rand34/
java -Xmx1g -cp stairway_plot_es/ Stairpainter ERY_200.blueprint

## Step 3: plot
# bash ERY_200.blueprint.plot.sh
