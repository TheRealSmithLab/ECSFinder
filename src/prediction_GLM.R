# Load necessary libraries
library(glmnet)
library(caret)

# Load the pretrained model
model_lasso <- readRDS("/home/vandalovejoy/R_scripts/ECS_project/model_lasso.rds")

# Define the file paths from the command line arguments
args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
output_file <- args[2]

# Read the input data and handle potential issues
test_data <- read.csv(input_file, header = TRUE, stringsAsFactors = FALSE)
if (nrow(test_data) == 0) {
  stop("Input file is empty or not properly read.")
}

# Ensure there are no missing or additional columns
expected_cols <- c('log_min_evalue', 'covarying_bp', 'min_energy', 'pseudo_energy', 'MPI')
if (!all(expected_cols %in% colnames(test_data))) {
  stop("Input file is missing one or more required columns.")
}

# Prepare the input data
x_test <- as.matrix(test_data[, expected_cols])

# Check the class of the model_lasso object
if ("train" %in% class(model_lasso)) {
  # Predict using the caret train model
  predicted_probs <- predict(model_lasso, newdata = x_test, type = "prob")[, 2] # assuming binary classification
} else if ("glmnet" %in% class(model_lasso)) {
  # Predict using the glmnet model directly
  predicted_probs <- predict(model_lasso, newx = x_test, s = model_lasso$lambda.min, type = "response")
  predicted_probs <- as.vector(predicted_probs)
} else {
  stop("Unknown model type: ", class(model_lasso))
}

threshold <- 0.38
predicted_class <- ifelse(predicted_probs >= threshold, "TP", "FP")

# Save the predictions to the output file without quotes
predictions <- data.frame(name_file = test_data$name_file, Predicted_Probabilities = predicted_probs, Predicted_Class = predicted_class)
write.table(predictions, output_file, sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE)

# Print the predictions
print(predictions)
