library(car)       # Provides the vif() function
library(dplyr)     # Used for data processing

# Define a function to calculate VIF
calculate_vif <- function(file_path) {
  # Read the CSV file
  df <- read.csv(file_path, stringsAsFactors = FALSE)
  
  # Select numeric columns as environmental variables
  numeric_cols <- df %>% select(where(is.numeric)) %>% names()
  
  if (length(numeric_cols) < 2) {
    stop("At least two numeric variables are required to calculate VIF.")
  }
  
  # Keep only numeric columns and remove rows with missing values
  df_clean <- df[numeric_cols] %>% na.omit()
  
  if (nrow(df_clean) == 0) {
    stop("All numeric variables contain missing values, so VIF cannot be calculated.")
  }
  
  # Build a linear regression model using the first variable as the response
  # and the remaining variables as predictors.
  # Here, VIF is calculated among all predictor variables.
  # A common approach is to regress each variable against all other variables.
  # car::vif() can be directly applied to a linear model object.
  # Although any variable can be used as the response variable, vif() requires
  # a model with an intercept.
  # To avoid the response variable affecting interpretation, we use all variables
  # to construct a multiple regression model and select the first variable as
  # the response variable.
  # This is acceptable because VIF depends only on the correlation structure
  # among predictor variables, not on the response variable.
  
  formula <- as.formula(paste(numeric_cols[1], "~", paste(numeric_cols[-1], collapse = " + ")))
  model <- lm(formula, data = df_clean)
  
  # Calculate VIF
  vif_values <- vif(model)
  
  # Output the results
  result <- data.frame(
    Variable = names(vif_values),
    VIF = as.numeric(vif_values),
    Note = ifelse(vif_values > 10, "Severe collinearity",
                  ifelse(vif_values > 5, "Moderate collinearity", "Normal"))
  )
  
  print(result)
}

# Call the function; modify the file path according to your actual file
calculate_vif("DC_network_nutrient.csv")