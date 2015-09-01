eval.FEM <- function(FEM, locations, CPP_CODE = FALSE)
{
  # EVAL_FEM_FD evaluates the FEM fd object at points (X,Y)
  #
  #        arguments:
  # X         an array of x-coordinates.
  # Y         an array of y-coordinates.
  # fit.FEM a FELspline object
  # FAST      a boolean indicating if the walking algorithm should be apply 
  #        output:
  # EVALMAT   an array of the same size as X and Y containing the value of 
  #           fit.FEM at (X,Y).
  
  res = NULL
  if(CPP_CODE == FALSE)
  {
    res = R_eval.FEM(FEM, locations)
  }else
  {
    print("Not implemented yet")
  }
  
  return(res)
}