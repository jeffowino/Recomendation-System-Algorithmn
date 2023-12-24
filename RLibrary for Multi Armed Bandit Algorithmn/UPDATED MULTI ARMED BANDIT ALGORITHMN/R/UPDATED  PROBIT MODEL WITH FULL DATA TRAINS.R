pkg.env <- new.env()

pkg.env$MAX_SURPRISE <- 5
pkg.env$MIN_SURPRISE <- pkg.env$MAX_SURPRISE * -1
pkg.env$PRIOR_VARIANCE <- 1.0

##Utilities

prior_mean <- function(prior_prob, beta, num_features){
  bias_mean <- qnorm(prior_prob) * (beta ** 2 + num_features)
  return(bias_mean)
}



clip <- function(x, a, b) {
  a + (x-a > 0)*(x-a) - (x-b > 0)*(x-b)
}

##Start of the Actual algorithmn
BOPR <- function(x, ...) UseMethod("BOPR")

BOPR.default <- function(x, y,
                                beta=0.05,
                                prior_prob=0.5,
                                epsilon=0.05,
                                subset=NULL,
                                ...)
{
  #init bias prior
  num_feat <- ncol(x)
  weight_mat <- matrix(nrow = 2,ncol = num_feat,data = 0, dimnames=list(c("mean", "variance"), colnames(x)))
  bias_mean <- prior_mean(prior_prob, beta, num_feat)
  
  #added bias mean and variance
  weight_mat[,seq(1,num_feat)] <- c(bias_mean, pkg.env$PRIOR_VARIANCE)
  if(!is.null(subset) & is.numeric(subset)){
    x <- x[subset,]
    y <- y[subset]
  }
  #train
  y_lab <- ifelse(y == 1, 1, -1)

  #active_mean_variance
  for(i in 1:nrow(x)){
    total_mean <- sum(weight_mat[1,])
    total_variance <- sum(weight_mat[2,]) + beta ** 2
    t <- y_lab[i] * total_mean / sqrt(total_variance)
    t <- clip(t, pkg.env$MIN_SURPRISE, pkg.env$MAX_SURPRISE )
    v <- dnorm(t) / pnorm(t)
    w <- v * (v + t)
    col_nm <- colnames(x)
    if(length(col_nm) != ncol(weight_mat))
      stop("num of features is not match.")
    for(j in col_nm){
        if(x[i,j] == 0) next
        mean_delta <- y_lab[i] * weight_mat['variance',j] / sqrt(total_variance) * v
        variance_multiplier <- 1.0 - weight_mat['variance',j] / total_variance * w
        updated <- c(
          mean = as.numeric(weight_mat['mean',j] + mean_delta),
          variance = as.numeric(weight_mat['variance',j] * variance_multiplier))
        #apply dynamics
        prior <- c(mean=0, variance=pkg.env$PRIOR_VARIANCE)
        adjusted_variance <- updated['variance'] * prior['variance'] /
        ((1.0 - epsilon) * prior['variance'] +
           epsilon * updated['variance'])
        adjusted_mean <- adjusted_variance * (
          (1.0 - epsilon) * updated['mean'] / updated['variance'] +
            epsilon * prior['mean'] / prior['variance'])

        weight_mat[,j] <- c(mean=adjusted_mean, variance=adjusted_variance)
    }

  }
  mdl_mat <- list(beta_matrix=weight_mat, beta=beta, prior_prob=prior_prob,epsilon=epsilon)
  class(mdl_mat) <- 'BOPR'
  return(mdl_mat)
}



#' @rdname BOPR
#' @export
BOPR.formula <- function(formula, data, subset=NULL, na.action =na.pass, beta=0.05,
                                prior_prob=0.5,
                                epsilon=0.05, ...)
{
  Call <- match.call()

  indx <- match(c("formula", "data", "subset"),
                names(Call), nomatch = 0L)
  if (indx[1] == 0L) stop("a 'formula' argument is required")
  temp <- Call[c(1L, indx)]      # only keep the arguments we wanted
  temp$na.action <- na.action    # This one has a default
  temp[[1L]] <- quote(stats::model.frame) # change the function called
  m <- eval.parent(temp)

  Terms <- attr(m, "terms")

  y <- model.extract(m, "response")

  m <- m[,-1,drop = FALSE]

  if(any(sapply(m, is.numeric) == TRUE)){
    stop("All variables need to be factor!")
  }

  character_vars <- lapply(m, class) == "character"
  m[, character_vars] <- lapply(m[, character_vars], as.factor)

  mdl_mat <- model.matrix(as.formula(paste("~",paste0(colnames(m), collapse=' + '))),
              m, contrasts.arg = lapply(m,contrasts, contrasts=FALSE))[,-1]

  object <- BOPR.default(x=mdl_mat, y=y, beta=beta, prior_prob=prior_prob,
                      epsilon=epsilon,subset=subset, ...)
  object$formula <- formula
  return(object)
}

###data is kicking in
.make_training_data <- function(formula, data, subset=NULL, na.action =na.pass)
{
  Call <- match.call()

  indx <- match(c("formula", "data", "subset"),
                names(Call), nomatch = 0L)
  if (indx[1] == 0L) stop("a 'formula' argument is required")
  temp <- Call[c(1L, indx)]      # only keep the arguments we wanted
  temp$na.action <- na.action    # This one has a default
  temp[[1L]] <- quote(stats::model.frame) # change the function called
  m <- eval.parent(temp)

  Terms <- attr(m, "terms")

  y <- model.extract(m, "response")

  m <- m[,-1,drop = FALSE]

  if(any(sapply(m, is.numeric) == TRUE)){
    stop("All variables need to be factor!")
  }

  character_vars <- lapply(m, class) == "character"
  m[, character_vars] <- lapply(m[, character_vars], as.factor)

  mdl_mat <- model.matrix(as.formula(paste("~",paste0(colnames(m), collapse=' + '))),
                          m, contrasts.arg = lapply(m,contrasts, contrasts=FALSE))[,-1]
  return(list(y=y, x=mdl_mat))
}

predict.BOPR <- function(object, newdata = NULL, type = "response", na.action = na.pass, ...){
  if(class(object) != 'BOPR'){
    stop('object is not BOPR class!')
  }
  if(any(class(newdata) == 'data.frame')) {
    variables <- labels(terms(object$formula))
    newdata_ <- newdata[,variables]
    character_vars <- lapply(newdata_, class) == "character"
    newdata_[, character_vars] <- lapply(newdata_[, character_vars], as.factor)
    mdl_mat <- model.matrix(as.formula(paste("~",paste0(variables, collapse=' + '))),
              newdata_, contrasts.arg = lapply(newdata_,contrasts, contrasts=FALSE))[,-1]
    if(!setequal(colnames(object$beta_matrix), colnames(mdl_mat))){
      inter_params <- intersect(colnames(object$beta_matrix), colnames(mdl_mat))
      beta_mat     <- object$beta_matrix[,inter_params]
      newdata.mat  <- mdl_mat[,inter_params]
    }else{
      beta_mat <- object$beta_matrix
      newdata.mat <- mdl_mat
    }
  }else if(any(class(newdata) == 'matrix')){
    newdata.mat <- newdata[,colnames(object$beta_matrix)]
  }else{
    stop('newdata must be one of matrix or data.frame!')
  }

  return(pnorm((newdata.mat %*% beta_mat['mean',])/(newdata.mat %*% beta_mat['variance',] + object$beta ** 2)))
}

#######
online_learning <- function(object, newdata=NULL, allow.new=TRUE){
  if(!any(class(newdata) == 'data.frame')){
    stop('newdata must be data.frame!')
  }
  res <- .make_training_data(object$formula, data = newdata)
  return(online_learning.matrix(object, x=res$x, y=res$y, allow.new=allow.new))
}

######
online_learning.matrix <- function(object, x = NULL, y=NULL, allow.new=TRUE){
  if(class(object) != 'BOPR'){
    stop('object is not BOPR class!')
  }
  y_lab <- ifelse(y == 1, 1, -1)

  #if there is other column in x
  if(!setequal(colnames(object$beta_matrix), colnames(x))){
    added_col <- setdiff(colnames(x), colnames(object$beta_matrix))
    print(sprintf('will added %s features!',length(added_col)))
    for(nm in added_col){
      bias_mean <- prior_mean(object$prior_prob, object$beta, ncol(object$beta_matrix) + 1)
      object$beta_matrix <- cbind(object$beta_matrix, c(mean=bias_mean, variance=pkg.env$PRIOR_VARIANCE))
      colnames(object$beta_matrix)[ncol(object$beta_matrix)] <- nm
    }
  }


  for(i in 1:nrow(x)){
    total_mean <- sum(object$beta_matrix[1,])
    total_variance <- sum(object$beta_matrix[2,]) + object$beta ** 2
    t <- y_lab[i] * total_mean / sqrt(total_variance)
    t <- clip(t, pkg.env$MIN_SURPRISE, pkg.env$MAX_SURPRISE )
    v <- dnorm(t) / pnorm(t)
    w <- v * (v + t)
    col_nm <- colnames(x)
    #if(length(col_nm) != ncol(object$beta_matrix))
    #  stop("num of features is not match.")
    for(j in col_nm){
      if(x[i,j] == 0) next
      mean_delta <- y_lab[i] * object$beta_matrix['variance',j] / sqrt(total_variance) * v
      variance_multiplier <- 1.0 - object$beta_matrix['variance',j] / total_variance * w
      updated <- c(
        mean = as.numeric(object$beta_matrix['mean',j] + mean_delta),
        variance = as.numeric(object$beta_matrix['variance',j] * variance_multiplier))
      #apply dynamics
      prior <- c(mean=0, variance=1)
      adjusted_variance <- updated['variance'] * prior['variance'] /
        ((1.0 - object$epsilon) * prior['variance'] +
           object$epsilon * updated['variance'])
      adjusted_mean <- adjusted_variance * (
        (1.0 - object$epsilon) * updated['mean'] / updated['variance'] +
          object$epsilon * prior['mean'] / prior['variance'])

      object$beta_matrix[,j] <- c(mean=adjusted_mean, variance=adjusted_variance)
    }

  }
  return(object)
}





setwd("~/ZIGMUNDPROJECT")

#source("C:\\Users\\USER\\OneDrive\\Documents\\ZIGMUNDPROJECT\\BOPR-master\\R\\util.R")

list.files()
Our_Data_Set=read.csv("Balanced Data.csv")
# Convert all columns to factors
Our_Data_Set <- data.frame(lapply(Our_Data_Set, as.factor))
idx  <- sample(1:nrow(Our_Data_Set))
first_train_set  <- Our_Data_Set[idx[1:10000],]
second_train_set  <- Our_Data_Set[idx[10001:20000],]


test_set <- Our_Data_Set[idx[21001:21909],]
#colnames(second_train_set)

# Get all variable names from the data
all_vars <- names(first_train_set)

# Remove 'addToCart' from the list
predictor_vars <- setdiff(all_vars, "addToCart")

# Create the formula dynamically
formula_string <- paste("addToCart ~", paste(predictor_vars, collapse = "+"))
formula_obj <- as.formula(formula_string)

# Fit the model using BOPR
#bopr_mdl <- BOPR(formula_obj, data = first_train_set, epsilon = 0.07)
#bopr_mdl_up  <- online_learning(bopr_mdl, second_train_set)

# Number of epochs
n_epochs <- 10

# Initialize BOPR model with first training set
bopr_mdl <- BOPR(formula_obj, data = first_train_set, epsilon = 0.07)

# Loop through each epoch
for(epoch in 1:n_epochs) {
  
  # Update the model with the second training set
  bopr_mdl_up <- online_learning(bopr_mdl, second_train_set)
  
  # Optional: print or log the epoch number and any relevant metrics
  print(paste("Epoch:", epoch))
  
  # Optional: evaluate your model on a validation set
  # validation_metrics <- evaluate_model(bopr_mdl, validation_set)
  # print(validation_metrics)
}


pred  <- predict(bopr_mdl_up, test_set)
pred 

typeof(test_set)
colnames(test_set)
library(ggplot2)

test_set$pred_2  <- predict(bopr_mdl_up, test_set)[,1]
ggplot(test_set, aes(pred_2)) + geom_density(aes(fill=factor(addToCart)), alpha=0.6) + xlim(0,1) + ggtitle("ε : 0.03, with  online learning")

 



# Define the user features
user_features <- list(
  browser_chrome = 0,
  browser_edge = 0,
  browser_firefox = 0,
  browser_other = 0,
  browser_safari = 1,
  os_android = 0,
  os_ios = 0,
  os_linux = 0,
  os_mac = 1,
  os_other = 0,
  os_windows = 0,
  part_of_day_evening = 1,
  part_of_day_morning = 0,
  part_of_day_night = 0,
  part_of_day_noon = 0,
  city_size_big = 0,
  city_size_normal = 1,
  city_size_small = 0,
  city_size_unknown = 0,
  city_size_very_big = 0,
  city_size_very_small = 0,
  referrer_internal = 0,
  referrer_other = 1,
  referrer_search = 1,
  brand_apple = 0,
  brand_huawei = 0,
  brand_other = 1,
  brand_samsung = 0,
  brand_xiaomi = 0,
  form_mobile = 0,
  form_other = 0,
  form_tablet = 1,
  part_of_week_weekday = 1,
  part_of_week_weekend = 0
  
)

variants <- c("variantName_achievment", "variantName_affiliation", "variantName_default", "variantName_power")

######Very strong one indeed 
# Convert user_features to a data frame
user_features_df <- as.data.frame(t(replicate(length(variants), unlist(user_features), simplify = "data.frame")))

# Set the variants
for (i in seq_along(variants)) {
  # Reset all variants to 0
  user_features_df[i, variants] <- 0
  
  # Set the current variant to 1
  user_features_df[i, variants[i]] <- 1
}

# Convert columns to factors if they were factors in the training set
for (col in colnames(user_features_df)) {
  if (is.factor(first_train_set[[col]])) {
    user_features_df[[col]] <- factor(user_features_df[[col]], levels = levels(first_train_set[[col]]))
  }
}

# Predict the conversion rate using the trained BOPR model
predicted_ctrs <- predict.BOPR(bopr_mdl_up, newdata = user_features_df)

# The predicted_ctrs vector contains the predicted conversion rates for each variant
predicted_ctrs



# Let's assume these are the variances for each variant's conversion rate
# You would need to obtain these from your model or assumptions about the data
variances <- rep(0.01, length(variants))####Have not used psoterior  variance
variances 
# Perform Thompson Sampling
samples <- rnorm(length(variants), mean = predicted_ctrs, sd = sqrt(variances))
samples
#samples <- rnorm(1000, mean = predicted_ctrs, sd = sqrt(variances))

# The action with the highest sampled value is selected
selected_variant_index <- which.max(samples)
selected_variant <- variants[selected_variant_index]

# The selected variant contains the name of the variant selected by Thompson Sampling
selected_variant

# Create a data frame containing the samples and corresponding variant names
samples_df <- data.frame(samples = samples, variant = rep(variants, each = length(samples)))
samples_df
samples
#############################
n_samples <- 4
samples_matrix <- matrix(nrow = n_samples, ncol = length(variants))
for (i in 1:length(variants)) {
  samples_matrix[,i] <- rnorm(n_samples, mean = predicted_ctrs[i], sd = sqrt(variances[i]))
}

samples_df <- as.data.frame(samples_matrix)
colnames(samples_df) <- variants

samples_long <- tidyr::gather(samples_df, "variant", "value")

ggplot(samples_long, aes(x = value, fill = variant)) +
  geom_density(alpha = 0.9) +
  xlim(0, 1.5) +
  ggtitle("Distribution of means")



############################
#Precison,Accuracy,Sensitity

# To make a decision, let's consider a threshold of 0.5 for the predicted probabilities
pred_class <- ifelse(pred > 0.5, 1, 0)

confusion_matrix <- table(test_set$addToCart, pred_class)
rownames(confusion_matrix) <- c("Actual 0", "Actual 1")
colnames(confusion_matrix) <- c("Predicted 0", "Predicted 1")
print(confusion_matrix)

accuracy <- sum(diag(confusion_matrix)) / sum(confusion_matrix)
print(paste("Accuracy:", accuracy))


library(pROC)

roc_obj <- roc(test_set$addToCart, pred)
auc_obj <- auc(roc_obj)
plot(roc_obj, main = paste("ROC Curve (AUC =", round(auc_obj, 2), ")"))


##Note that the accuracy and other metrics may not always be the best way to evaluate a model, especially if your dataset is imbalanced. The AUC-ROC curve is often a good metric for such datasets since it evaluates the model's ability to distinguish between the classes independent of a threshold.
auc_value <- auc(roc_obj)

print(paste("AUC Value:", auc_value))


####RECAL
# Recall (Sensitivity) calculation:
TP <- confusion_matrix[2, 2]  # True Positives (assuming positive class is in second position)
FN <- confusion_matrix[1, 2]  # False Negatives (assuming positive class is in second position)

recall <- TP / (TP + FN)

print(paste("Recall Value:", recall))



# Assuming the provided code is correct and included above

# Extract values from the confusion matrix
TP <- confusion_matrix[2, 2]  # True Positives 
FN <- confusion_matrix[1, 2]  # False Negatives 
FP <- confusion_matrix[2, 1]  # False Positives (assuming positive class is in second position)
TN <- confusion_matrix[1, 1]  # True Negatives (assuming positive class is in second position)

# Compute MCC avoiding integer overflow
MCC <- (TP * TN - FP * FN) / 
  (sqrt(TP + FP) * sqrt(TP + FN) * sqrt(TN + FP) * sqrt(TN + FN))

# Print MCC
print(paste("Matthews Correlation Coefficient (MCC):", MCC))




