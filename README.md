OECD Gas Demand 
- One-step additive Non-parametric IV used to obtain the effect of prices on gasoline demand in OECD
- The final graph represents the findings instrumenting prices vs. assuming the prices are exogenous
- Accounting for endogeneity produces more sensible results: the obtained elasticities correspond to the ones found in literature
- Nonparametric framework gives rise to 2SLS estimator for this application

OECDGas is the main file that implements the method developed in my Master’s thesis on the One-Step Additive Non-Parametric IV Estimator. 
It uses the bits_comparison.R file, which contains the necessary functions to run the method.

Some parameters can be adjusted in the main file (OECDGas), such as the grid over which lambda (the regularization parameter) is estimated, 
the number of folds used in cross-validation, and more. The number of repetitions can also be changed, although this may increase computational cost.

Demand function specification:

ln(gas_it/car_it) = η_i + g(ln(price_it)) + ε_it 

where  η_i are country-fixed effects, g some function in prices, and ε_it an error term.

First-differencing this panel gets rid of fixed effects and we end up in an additive non-parametric framework to which I 
can apply the model:

∆ln(gas_it/car_it) = g_t(ln(price_it)) - g_(t-1)(ln(price_it)) + ε_it - ε_i(t-1)
