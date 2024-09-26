OECDGas is the main file that implements the method developed in my Master’s thesis on the One-Step Additive Non-Parametric IV Estimator. It uses the bits_comparison.R file, which contains the necessary functions to run the method.

Some parameters can be adjusted in the main file (OECDGas), such as the grid over which lambda (the regularization parameter) is estimated, the number of folds used in cross-validation, and more. The number of repetitions can also be changed, although this may increase computational cost.
