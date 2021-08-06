####################################
Include the model sources to the CIF
####################################

Before starting to implement a model to the CIF, please separate elements
that are directly involved in running your numerical model from other elements that produce inputs and read outputs.
The later ones will be included later into the ''model'' pyCIF class, while the first ones are called ''model sources''
below.

Create a new directory in the the folder ''model_sources/'' including all sources for your numerical model.
It includes only codes in Fortran or C, or any other language that can be executed to run the model and produce outputs.

It is highly recommended that the model is built in a way that allows users to compile it once for all
and then pyCIF will only call the corresponding executable anytime necessary.

If necessary, please include here any information useful for compiling the model
(third-party libraries, compiling options, etc.).

