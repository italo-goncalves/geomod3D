# geomod3D
Implicit 3D geological modeling and geostatistics

Welcome to the geomod3D package. This is a project to develop a R package with 
tools for 3D implicit modeling of geological bodies and general manipulation of 
3D spatial data and drillholes. See the vignette for instructions and a 
demonstration (if you have trouble viewing it, check the online version on
<http://rpubs.com/italogoncalves/geomod3D-intro>).

Open 3D geological data is *very* welcome. Send me an e-mail if you can provide 
some or if you know of a repository.

You can install the package with the command 
`devtools::install_github("italo-goncalves/geomod3D")`. The 
`build_vignette = T` option builds the vignette locally, but it can take a 
few hours.

TODO list:
* Add data and a vignette containing structural measurements
* Add examples to help files
* Add more open data and modeling examples
* Improve documentation
* Add data and examples of modeling with continuous variables
* Speed up the `GP` class predictions
* Add infrastructure for multivariate modeling
* Implement GPU computations
* Implement package testing
