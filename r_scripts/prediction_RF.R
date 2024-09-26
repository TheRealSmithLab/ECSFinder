# Load necessary libraries
library(randomForest)
library(caret)



# Define the file paths from the command line arguments
args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
output_file <- args[2]
model_dir <- args[3]
# Load the pretrained Random Forest model
model_rf <- readRDS(paste0(model_dir,"/model_rf.rds"))

# Read the input data and handle potential issues
test_data <- read.csv(input_file, header = TRUE, stringsAsFactors = FALSE)
if (nrow(test_data) == 0) {
  stop("Input file is empty or not properly read.")
}

# Ensure there are no missing or additional columns
expected_cols <- c("min_energy", "pseudo_energy", "log_min_evalue", "covarying_bp", "MPI", "average_MFE_sample",  "sd_sample","zscore")

if (!all(expected_cols %in% colnames(test_data))) {
  stop("Input file is missing one or more required columns.")
}

# Prepare the input data
x_test <- test_data[, expected_cols]
# Replace -Inf in log_min_evalue with a very small value
x_test$log_min_evalue[is.infinite(x_test$log_min_evalue)] <- -1e6



if ("train" %in% class(model_rf)) {
  # Predict using the caret train Random Forest model
  predicted_probs <- predict(model_rf, newdata = x_test, type = "prob")[, 2]  
} else if ("randomForest" %in% class(model_rf)) {
  # Predict using the randomForest model directly
  predicted_probs <- predict(model_rf, newdata = x_test, type = "prob")[, 2]
} else {
  stop("Unknown model type: ", class(model_rf))
}

# Set a prediction threshold
threshold <- 0.75
predicted_class <- ifelse(predicted_probs >= threshold, "TP", "FP")

# Save the predictions to the output file without quotes
predictions <- data.frame(name_file = test_data$name_file, Predicted_Probabilities = predicted_probs, Predicted_Class = predicted_class)
# Check if the file exists
if (!file.exists(output_file)) {
  # File does not exist, write the table with headers
  write.table(predictions, output_file, sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE)
  cat("File created and data written.\n")
} else {
  # File exists, append the new data without column headers
  write.table(predictions, output_file, sep = ",", row.names = FALSE, col.names = FALSE, quote = FALSE, append = TRUE)
  cat("Data appended to existing file.\n")
}


# Print the predictions
print(predictions)
