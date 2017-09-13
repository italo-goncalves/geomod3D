## Version 0.0.0.9002

* Corrected a bug in `SPGP` object initialization the prevented the use of
`points3DDataFrame` objects as pseudo-inputs
* Corrected a bug that caused the `GetContacts()` method to covert the labels
to factors
* Corrected a bug in the `Predict()` method that caused an error for 
`SPGP_geomod` objects

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
