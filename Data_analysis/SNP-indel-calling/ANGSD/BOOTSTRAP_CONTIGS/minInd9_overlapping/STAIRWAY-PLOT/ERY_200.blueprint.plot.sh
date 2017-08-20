# rand34: average logLikelihood=-1341.643679098878 median logLikelihood=-1311.051066402346
# best number of break points (average logLikelihood):34
# best number of break points (median logLikelihood):34
# Step 1: estimations for the "final" result
cp -f out_ERY_200/rand34/*addTheta out_ERY_200/final/
# Step 2: create summaries and plots
java -Xmx4g -cp stairway_plot_es/:stairway_plot_es/gral-core-0.11.jar:stairway_plot_es/VectorGraphics2D-0.9.3.jar Stairway_output_summary_plot ERY_200.blueprint
java -Xmx4g -cp stairway_plot_es/:stairway_plot_es/gral-core-0.11.jar:stairway_plot_es/VectorGraphics2D-0.9.3.jar Stairway_output_summary_plot ERY_200.blueprint rand34
