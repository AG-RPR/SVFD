#### Function to "decompose" visual field measurements (total deviations) into the 17 archetypes

# requires the R library 'archetypes' which can be installed by
# install.packages("archetypes")
require("archetypes")

# the 17 archetypes are stored in the R data structure "at17"
# which is loaded from file "archetypes.RData";
# make sure the file is in the same directory, or modify the following line accordingly:
if(!exists("at17")) load("archetypes.RData")


# calculate the 17 archetype coefficients of TD measurements td.measurement,
# either given as a vector of 52 elements or, to facilitate the calculation 
# of more than one measurement, as a matrix with 52 columns:
decompose.vf <- function(td.measurement)
{
	if(is.null(dim(td.measurement)))
	{
		# the following line is a workaround for a bug
		# (if only a single row is given as argument, R crashes!!!):
		x = rbind(td.measurement, numeric(52))
		predict(at17, x)[1,]
	}
	else
		predict(at17, td.measurement)
}


### Example:
# taken from the figure "Illustrative examples of visual field measurements and their
# decompositions into archetypes", subfigure A

# tds = c(-20, -24, -29, -13, -16, -6, -17, -21, -19, -18, -9, -11, -8, -7, -5, -25, -14, -2, -3, -20, -4, -5, -3, -1, -3, -3, -2, -2, -1, -3, -2, 0, 0, -3, 0, -1, 0, -1, -1, -4, -1, -1, -3, -2, -2, 0, -1, -4, -7, -2, -2, -2)
# decompose.vf(tds)
