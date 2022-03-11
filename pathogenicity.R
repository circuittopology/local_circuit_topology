# Single tree

install.packages("ISLR")
library(ISLR)
data(package="ISLR")
install.packages("tree")
require(tree)
install.packages("readxl")
library("readxl")
data <- read_excel("aa_E_nc_CT_CO_essential.xlsx")
pathogenic = data$pathogenic
data$aa1 <- as.factor(data$aa1)
data$aa2 <- as.factor(data$aa2)
data$essential <- as.factor(data$essential)
t = tree(pathogenic~., data=data)
plot(t)
text(t)


# Random Forests

install.packages("randomForest")
library(randomForest)
set.seed(17)
data_set_size <- floor(4*nrow(data)/5)
# Generate a random sample of "data_set_size" indexes
indexes <- sample(1:nrow(data), size = data_set_size)
# Assign the data to the correct sets
training <- data[indexes,]
validation1 <- data[-indexes,]
t <- as.factor(data$pathogenic)
rf_classifier = randomForest(factor(pathogenic)~., data=training,ntree=500,mtry=2,importance=TRUE)
prediction_for_table <- predict(rf_classifier,validation1[,-1])
table(t(validation1[,1]),predicted=prediction_for_table)


prediction_for_roc_curve <- predict(rf_classifier,validation1[,-1],type="prob")
pretty_colours <- c("#F8766D","#00BA38")
classes <- levels(validation1$pathogenic)

  # Define which observations belong to class[i]
  true_values <- ifelse(validation1[,1]==1,1,0)
  # Assess the performance of classifier for class[i]
  pred <- prediction(prediction_for_roc_curve[,2],true_values)
  perf <- performance(pred, "tpr", "fpr")
  
  plot(perf,main="ROC Curve",col=pretty_colours[1]) 
  # Calculate the AUC and print it to screen
  auc.perf <- performance(pred, measure = "auc")
  print(auc.perf@y.values)

  
  