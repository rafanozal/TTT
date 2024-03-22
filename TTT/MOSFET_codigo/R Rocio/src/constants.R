# Files and folders
{
 
  # From which folder do you read the data and into what folder write the results
  DATA_FOLDER        = file.path(paste(getwd(),"/../data/", sep = ""))
  RESULT_FOLDER      = file.path(paste(getwd(),"/../out/", sep = ""))
  
  # Names of all the files we have
  LOG_NOISE_FILENAME = "all_conditions_lognoise.txt"
  NO_NOISE_FILENAME  = "all_conditions_noise.txt"
  NO_NULL_FILENAME   = "all_conditions_noise_no0.txt"
  LARGO_FILENAME     = "lognoise_largo.txt"
  BACKFIT_FILENAME   = "datos_backfit.txt"
  
  # Filepath to those files
  LOG_NOISE_FILEPATH = file.path(paste(DATA_FOLDER, LOG_NOISE_FILENAME, sep = ""))
  NO_NOISE_FILEPATH  = file.path(paste(DATA_FOLDER, NO_NOISE_FILENAME,  sep = ""))
  NO_NULL_FILEPATH   = file.path(paste(DATA_FOLDER, NO_NULL_FILENAME,   sep = ""))
  LARGO_FILEPATH     = file.path(paste(DATA_FOLDER, LARGO_FILENAME,     sep = ""))
  BACKFIT_FILEPATH   = file.path(paste(DATA_FOLDER, BACKFIT_FILENAME,   sep = ""))
  
  # Create the result folder
  dir.create(file.path(RESULT_FOLDER),            showWarnings = FALSE) # You will always get a warning for this because this folder is probably already created, thus suppresing the warning system here
  
}