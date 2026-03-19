## Homework 10

############################################################### Question 14.1

data <- read.csv("breast-cancer-wisconsin.data.txt", 
                 header=FALSE, 
                 na.strings="?")

colnames(data) <- c("ID", "Clump_thickness", "Uniformity_size", 
                    "Uniformity_shape", "Marginal_adhesion",
                    "Epithelial_size", "Bare_nuclei", 
                    "Bland_chromatin", "Normal_nucleoli", 
                    "Mitoses", "Class")

#View(data)
#table(data$Class[is.na(data$Bare_nuclei)])
#table(data$Bare_nuclei)
#sum(is.na(data$Bare_nuclei))

hist(data$Bare_nuclei, 
     breaks=seq(0.5, 10.5, by=1),  
     main="Distribution of Bare Nuclei",
     xlab="Bare Nuclei", 
     col="purple")

bare_mean <- mean(data$Bare_nuclei, na.rm=TRUE)

get_mode <- function(x) {
  x <- x[!is.na(x)]
  freq_table <- table(x)
  mode_value <- as.numeric(names(which.max(freq_table)))
  
  return(mode_value)
}

bare_mode <- get_mode(data$Bare_nuclei)
data_mode <- data
data_mode$Bare_nuclei[is.na(data_mode$Bare_nuclei)] <- bare_mode

#View(data_mode)
#sum(is.na(data_mode$Bare_nuclei))

############################################################### Question 14.2

complete <- data[!is.na(data$Bare_nuclei), ]
missing <- data[is.na(data$Bare_nuclei), ]

model <- lm(Bare_nuclei ~ . - ID - Class, data = complete)
predicted <- predict(model, newdata = missing)

#View(predicted)

predicted <-round(predicted)
predicted <- as.integer(pmin(pmax(predicted, 1), 10))

#View(predicted)

data_regression <- data
data_regression$Bare_nuclei[is.na(data_regression$Bare_nuclei)] <- predicted

#View(data_regression)

############################################################### Question 14.3

model_perturbation <- lm(Bare_nuclei ~ . - ID - Class, data = complete)
predicted <- predict(model_perturbation, newdata = missing)
sigma <- sigma(model_perturbation)

set.seed(42)  
noise <- rnorm(nrow(missing), mean = 0, sd = sigma)
predicted_perturbed <- predicted + noise

predicted_perturbed <- round(predicted_perturbed)
predicted_perturbed <- pmin(pmax(predicted_perturbed, 1), 10)

data_perturbation <- data
data_perturbation$Bare_nuclei[is.na(data_perturbation$Bare_nuclei)] <- predicted_perturbed

#View(data_perturbation)

############################################################### Question 14.4.1

##### SVM COMPARISONS 

library(kernlab)
library(caret)

data_mode$Class <- factor(data_mode$Class, levels = c(2, 4), labels = c("Benign", "Malignant"))
data_regression$Class <- factor(data_regression$Class, levels = c(2, 4), labels = c("Benign", "Malignant"))
data_perturbation$Class <- factor(data_perturbation$Class, levels = c(2, 4), labels = c("Benign", "Malignant"))

data_mode <- data_mode[, -1]
data_regression <- data_regression[, -1]
data_perturbation <- data_perturbation[, -1]

#View(data_mode)
#View(data_regression)
#View(data_perturbation)

ctrl <- trainControl(
  method = "cv", 
  number = 10,
  savePredictions = "final"
)

set.seed(42)
svm_mode <- train(Class ~ ., data = data_mode, 
                  method = "svmLinear2", 
                  trControl = ctrl,
                  tuneGrid = data.frame(cost = 1),
                  preProcess = c("center", "scale"))

set.seed(42)
svm_regression <- train(Class ~ ., data = data_regression, 
                        method = "svmLinear2", 
                        trControl = ctrl,
                        tuneGrid = data.frame(cost = 1),
                        preProcess = c("center", "scale"))

set.seed(42)
svm_perturbation <- train(Class ~ ., data = data_perturbation, 
                          method = "svmLinear2", 
                          trControl = ctrl,
                          tuneGrid = data.frame(cost = 1),
                          preProcess = c("center", "scale"))

## RESULTS:

# Accuracy
cat("Mode imputation accuracy:", svm_mode$results$Accuracy, "\n")
cat("Regression imputation accuracy:", svm_regression$results$Accuracy, "\n")
cat("Regression with perturbation accuracy:", svm_perturbation$results$Accuracy, "\n")

# Mode Imputation
cat("Mode Imputation\n")
print(table(Actual = svm_mode$pred$obs, Predicted = svm_mode$pred$pred))

# Regression Imputation
cat("Regression Imputation\n")
print(table(Actual = svm_regression$pred$obs, Predicted = svm_regression$pred$pred))

# Regression with Perturbation
cat("Regression with Perturbation\n")
print(table(Actual = svm_perturbation$pred$obs, Predicted = svm_perturbation$pred$pred))

##### kNN COMPARISONS

#View(data_mode)
#View(data_regression)
#View(data_perturbation)

# kNN on mode imputed data
set.seed(42)
knn_mode <- train(Class ~ ., data = data_mode, 
                  method = "knn", 
                  trControl = ctrl,
                  tuneGrid = data.frame(k = c(3, 5, 7, 9)),
                  preProcess = c("center", "scale"))

# kNN on regression imputed data
set.seed(42)
knn_regression <- train(Class ~ ., data = data_regression, 
                        method = "knn", 
                        trControl = ctrl,
                        tuneGrid = data.frame(k = c(3, 5, 7, 9)),
                        preProcess = c("center", "scale"))

# kNN on regression with perturbation
set.seed(42)
knn_perturbation <- train(Class ~ ., data = data_perturbation, 
                          method = "knn", 
                          trControl = ctrl,
                          tuneGrid = data.frame(k = c(3, 5, 7, 9)),
                          preProcess = c("center", "scale"))

# Best K for each
cat("Best K - Mode:", knn_mode$bestTune$k, "\n")
cat("Best K - Regression:", knn_regression$bestTune$k, "\n")
cat("Best K - Perturbation:", knn_perturbation$bestTune$k, "\n")

# Confusion matrices
cat("Mode Imputation kNN\n")
print(table(Actual = knn_mode$pred$obs, Predicted = knn_mode$pred$pred))

cat("Regression Imputation kNN\n")
print(table(Actual = knn_regression$pred$obs, Predicted = knn_regression$pred$pred))

cat("Regression with Perturbation kNN\n")
print(table(Actual = knn_perturbation$pred$obs, Predicted = knn_perturbation$pred$pred))

############################################################### Question 14.4.2

data_deleted <- data[complete.cases(data), ]
View(data_deleted)

data_deleted$Class <- factor(data_deleted$Class, levels = c(2, 4), labels = c("Benign", "Malignant"))
data_deleted <- data_deleted[, -1]
View(data_deleted)

# SVM
set.seed(42)
svm_deleted <- train(Class ~ ., data = data_deleted, 
                     method = "svmLinear2", 
                     trControl = ctrl,
                     tuneGrid = data.frame(cost = 1),
                     preProcess = c("center", "scale"))

# KNN
set.seed(42)
knn_deleted <- train(Class ~ ., data = data_deleted, 
                     method = "knn", 
                     trControl = ctrl,
                     tuneGrid = data.frame(k = c(3, 5, 7, 9)),
                     preProcess = c("center", "scale"))

## RESULTS:
cat("Best K:", knn_deleted$bestTune$k, "\n")

cat("Missing Values Deleted SVM\n")
print(table(Actual = svm_deleted$pred$obs, Predicted = svm_deleted$pred$pred))

cat("Missing Values Deleted kNN\n")
print(table(Actual = knn_deleted$pred$obs, Predicted = knn_deleted$pred$pred))

############################################################### Question 14.4.3

data_binary <- data
data_binary$Bare_nuclei_missing <- ifelse(is.na(data_binary$Bare_nuclei), 1, 0)
View(data_binary)

data_binary$Bare_nuclei[is.na(data_binary$Bare_nuclei)] <- bare_mode
View(data_binary)

data_binary$Class <- factor(data_binary$Class, levels = c(2, 4), labels = c("Benign", "Malignant"))
data_binary <- data_binary[, -1]
View(data_binary)

# SVM
set.seed(42)
svm_binary <- train(Class ~ ., data = data_binary, 
                    method = "svmLinear2", 
                    trControl = ctrl,
                    tuneGrid = data.frame(cost = 1),
                    preProcess = c("center", "scale"))

# kNN
set.seed(42)
knn_binary <- train(Class ~ ., data = data_binary, 
                    method = "knn", 
                    trControl = ctrl,
                    tuneGrid = data.frame(k = c(3, 5, 7, 9)),
                    preProcess = c("center", "scale"))

### RESULTS:
cat("Best K:", knn_binary$bestTune$k, "\n")

cat("Binary Variable SVM\n")
print(table(Actual = svm_binary$pred$obs, Predicted = svm_binary$pred$pred))

cat("Binary Variable kNN\n")
print(table(Actual = knn_binary$pred$obs, Predicted = knn_binary$pred$pred))


##### Plotting all differences:

results <- data.frame(
  Method = rep(c("Mode", "Regression", "Regression with Perturbation", "Deletion", "Binary Indicator"), each = 2),
  Classifier = rep(c("SVM", "KNN"), 5),
  Accuracy = c(
    svm_mode$results$Accuracy,
    knn_mode$results$Accuracy[which.max(knn_mode$results$Accuracy)],
    
    svm_regression$results$Accuracy,
    knn_regression$results$Accuracy[which.max(knn_regression$results$Accuracy)],
    
    svm_perturbation$results$Accuracy,
    knn_perturbation$results$Accuracy[which.max(knn_perturbation$results$Accuracy)],
    
    svm_deleted$results$Accuracy,
    knn_deleted$results$Accuracy[which.max(knn_deleted$results$Accuracy)],
    
    svm_binary$results$Accuracy,
    knn_binary$results$Accuracy[which.max(knn_binary$results$Accuracy)]
  )
)

print(results)
print(results$Method)

library(ggplot2)

ggplot(results, aes(x = Method, y = Accuracy, fill = Classifier)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  geom_text(aes(label = round(Accuracy, 4)), 
            position = position_dodge(width = 0.7), 
            vjust = -0.5, size = 3) +
  scale_fill_manual(values = c("SVM" = "purple1", "KNN" = "plum3")) +
  scale_x_discrete(limits = c("Mode", "Regression", "Regression with Perturbation", "Deletion", "Binary Indicator")) +
  coord_cartesian(ylim = c(0.95, 0.98)) +
  theme_minimal() +
  labs(title = "Classification Accuracy by Imputation Method",
       x = NULL,
       y = "Accuracy")






