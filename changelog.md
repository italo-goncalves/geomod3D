## Version 0.0.0.9003

* GP training now uses a dedicated genetic algorithm
* Bug fixes for Cauchy covariance
* Covariance matrix calculations implemented as proper internal methods (no more
`covfun` slots)
* Removed `weights` argument in `GP` and `SPGP`
* Added `data` as default value for `pseudo_inputs` in `SPGP`
* Constant and linear covariances now supported
* Testing a more stable method for matrix inversions and Cholesky decomposition

## Version 0.0.0.9002

* Corrected a bug in `SPGP` object initialization the prevented the use of
`points3DDataFrame` objects as pseudo-inputs
* Corrected a bug that caused the `GetContacts()` method to convert the labels
to factors
* Corrected a bug in the `Predict()` method that caused an error for 
`SPGP_geomod` objects
* Added a default for the `midrange` parameter in `anisotropy3d()` method
* Corrected a bug in the `Simulate()` method for the standard GP
* Added regularization to the smoothing covariance matrix for SPGP simulations
* Corrected a bug in the `Simulate()` method for the SPGP, where the mean was not
added back when making smoothed simulations
* Corrected a bug in the `Simulate()` method for the SPGP, where the mean was 
added improperly during the simulation, generating inflated values
* SPGP simulation now uses only the approximated variance as prior to avoid
value inflation
* Creation of the `covarianceModel3D` class to handle the covariances and the
nugget effect on the same object
* `points3DDataFrame` can now be initialized with no arguments
* Corrected the computations of the log-likelihood
* `GP_geomod` can now handle missing data labels
* `SPGP` method `Predict` now outputs a value that measures the quality of the
sparse approximation
* `GP_geomod` now uses the same covariance model for all classes
* Improved documentation

## Version 0.0.0.9001

* Formalization of structural data as a `directions3DDataFrame` object
* Support for block models as `blocks3DDataFrame` objects
* Inclusion of Sparse Gaussian Processes for regression and classification (implicit modeling)
* Support for simulations
* Removed support for covariance matrices in variogram form
* Support of cross-validation for Sparse Gaussian Processes
* Implementation of `as.data.frame()` method
* Unification of covariance methods for all kinds of spatial objects
* Added a parameter for nugget effect of structural data in the `GP` object
* Calculation of the implicit model's probabilities is now done through sampling
