
#' Class to represent a Renewal Equation Epidemic Model.
#'
#' @field name character. 
#' @field prms list. 
#' @field is.fitted logical. 
#' 
#' @return
#' @export
#'
#' @examples
#' 
setRefClass(
  Class = 'reem', 
  fields =  list(
    name      = "character",
    prms      = "list",
    is.fitted = "logical"
  ),
  
  # - - - - - - - - - - - - - - - - - -
  # - - - `reem` CLASS METHODS 
  # - - - - - - - - - - - - - - - - - -
  
  methods = 
    list(
      
      print_prms = function(){
        cat('\n--- Parameters for REEM `',name,'`\n')
        for(i in seq_along(prms)){
          cat(names(prms)[i],' = ', prms[[i]],'\n')
        }
        cat(' --------------------------\n')
      },
      
      
      simulate = function(deterministic){
        res = reem_simulate(prms, deterministic)  
        return(res)  
      }
      
    ))




