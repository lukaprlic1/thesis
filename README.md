OECDGas is the main file that implements the method developed in my Master Thesis on One-Step Additive Non-parametric IV Estimator.
OECDGas uses bits_comparison.R file shich contains the functions necessary to run the method.

Some of the parameters can be changed in the main file (OECDGas), like the grid on which lambda (regularization parameter) is estimated, or the number of fold used in cross-validation etc.
Number of repetitions can also be changed, however at the expense of computational cost.
