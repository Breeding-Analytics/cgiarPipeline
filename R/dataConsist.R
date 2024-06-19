#------------------------------------------------------------------------------
# Name for the module (to identify in the modifications table): dcRaw
# 
# Functions for data consistency
# 1. Check names
# 2. Check consistency
#    2.1. Only informative (list of warnings): comprehensive
#    2.2. Modifications: only for yields but potentially extended to other traits
#         To be implemented later
# 3. Compute variables? (maybe in some other module, after loading data)
# 
# Output for check consistency:
# - A data frame with a list of warnings for inconsistencies found
# - A data frame, same dimention as fieldbook, with indices to spot the values 
#   with inconsistencies. This could be used for:
#   - visualization and
#   - modifications => Tag inconsistencies and then change values to NA.
# 
# Scope: Only inconsistencies as outliers are already detected in the QA module
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# 1. Check names
# 1.1. Identify names
#      - Convert CO numbers to short labels
#      - Output: text listing all traits not recognized
# 1.2. Renaming (matching old with new names).
#      Proposal: Use the menu Other functions / Trait Transformations
# 
# Processes 1.1 and 1.2 will modify the data.
# Options:
# - Add new traits (columns) with the new names
# - Add matching list in the metadata table
#------------------------------------------------------------------------------

checkNames <- function(crop,
                       phenoDTfile, # The data object structure produced from bioflow
                       ) {
  
  # Get internal copy of data and convert to short labels
  
  mydata <- phenoDTfile$data$pheno
  
  if (crop == 'potato')
    mydata <- suppressWarnings(st4gi::convert.co.pt(mydata, 'co.to.labels'))
  
  if (crop == 'sweetpotato')
    mydata <- suppressWarnings(st4gi::convert.co.sp(mydata, 'co.to.labels'))
  
  # Check names (this will change all names to lower letters)
  
  if (crop == 'potato')
    tmp <- st4gi::check.names.pt(mydata)
  
  if (crop == 'sweetpotato')
    tmp <- st4gi::check.names.sp(mydata)
  
  # Output: check.names.xx will produce a warning with a list of names not recognized
  # so the user can match those with standard labels
  
  return(tmp)

}

matchNames <- function(crop, phenoDTfile) {
  
  # This should give the option to match column names for traits
  # This should be reactive? produce changes in the checkNames function output?
  
}

#------------------------------------------------------------------------------
# 2. Check consistency
#    - This will spot potential problems
#    - Only inconsistencies will be reported as outliers' detection is in the QA module
#    Outputs in the Dashboard:
#    - List of inconsistencies in data.frame format.
#    - Matrix representation of the data.frame with colors to spot potential problems:
#      - Red cells: inconsistencies (impossible values)
#      This will need modification of the functions to point out specific cells
#------------------------------------------------------------------------------

checkData <- function(crop,               # So far only potato and sweetpotato
                      phenoDTfile,        # The data.frame output from checkNames
                      f = 0,              # To avoid running outliers' detection
                      print.text = FALSE  # To avoid printing output
                      ) { 

  # Get internal copy of data after running checkNames and matchNames
  
  mydata <- phenoDTfile
  
  # Run check consistency

  if (crop == 'potato')
    output <- suppressWarnings(st4gi::check.data.pt(mydata, f, print.text = print.text))
  
  if (crop == 'sweetpotato')
    output <- suppressWarnings(st4gi::check.data.sp(mydata, f, print.text = print.text))
  
  # Output has 2 components:
  # $Inconsist.List: a data.frame with a list of all inconsistencies.
  # $Inconsist.Matrix: a data.frame with the positions in the fieldbook data frame where inconsistencies occur.
  
  return(output)
  
}
