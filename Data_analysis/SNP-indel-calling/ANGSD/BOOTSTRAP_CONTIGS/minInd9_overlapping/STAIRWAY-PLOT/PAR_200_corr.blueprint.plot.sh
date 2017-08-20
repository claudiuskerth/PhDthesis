# rand34: average logLikelihood=-5376.0078385510205 median logLikelihood=-5312.361039714888
# best number of break points (average logLikelihood):34
# best number of break points (median logLikelihood):34
# Step 1: estimations for the "final" result
cp -f out_PAR_200_corr/rand34/*.addTheta out_PAR_200_corr/final/
# Step 2: create summaries and plots
java -Xmx4g -cp stairway_plot_es/:stairway_plot_es/gral-core-0.11.jar:stairway_plot_es/VectorGraphics2D-0.9.3.jar Stairway_output_summary_plot PAR_200_corr.blueprint
java -Xmx4g -cp stairway_plot_es/:stairway_plot_es/gral-core-0.11.jar:stairway_plot_es/VectorGraphics2D-0.9.3.jar Stairway_output_summary_plot PAR_200_corr.blueprint rand34
