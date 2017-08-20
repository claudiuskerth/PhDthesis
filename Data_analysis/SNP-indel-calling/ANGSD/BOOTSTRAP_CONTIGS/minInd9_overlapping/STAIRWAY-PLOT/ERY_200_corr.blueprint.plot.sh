# rand34: average logLikelihood=-1182.31325473289 median logLikelihood=-1154.1075757853687
# best number of break points (average logLikelihood):34
# best number of break points (median logLikelihood):34
# Step 1: estimations for the "final" result
cp -f out_ERY_200_boot_corr/rand34/*.addTheta out_ERY_200_boot_corr/final/
# Step 2: create summaries and plots
java -Xmx4g -cp stairway_plot_es/:stairway_plot_es/gral-core-0.11.jar:stairway_plot_es/VectorGraphics2D-0.9.3.jar Stairway_output_summary_plot ERY_200_corr.blueprint
java -Xmx4g -cp stairway_plot_es/:stairway_plot_es/gral-core-0.11.jar:stairway_plot_es/VectorGraphics2D-0.9.3.jar Stairway_output_summary_plot ERY_200_corr.blueprint rand34
