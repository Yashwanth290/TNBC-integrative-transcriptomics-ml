############################################
# 🧬 1️⃣ SETUP
############################################
rm(list=ls())

dir.create("MP_AI_Project", showWarnings = FALSE)
setwd("MP_AI_Project")

dirs <- c("data_raw","data","features","models","results","plots")
for(d in dirs) dir.create(d, showWarnings = FALSE)

############################################
# 📦 2️⃣ PACKAGES
############################################
packages <- c("readr","dplyr","stringr","httr","jsonlite",
              "Peptides","caret","xgboost","pROC","ggplot2")

for(p in packages){
  if(!require(p, character.only=TRUE)){
    install.packages(p, dependencies=TRUE)
    library(p, character.only=TRUE)
  }
}

############################################
# 🧬 3️⃣ EXTRACT IDS
############################################
extract_ids <- function(file){
  lines <- readLines(file)
  headers <- lines[grepl("^>", lines)]
  ids <- stringr::str_extract(headers, "\\|([A-Z0-9]+)\\|")
  ids <- gsub("\\|","",ids)
  data.frame(Uniprot_ID = ids)
}

pos <- extract_ids("data_raw/pos.csv")
neg <- extract_ids("data_raw/neg.csv")

pos$label <- 1
neg$label <- 0

df <- unique(rbind(pos, neg))

############################################
# 🧬 4️⃣ LOAD UNIPROT
############################################
uniprot <- read.delim("data/uniprot.tsv", stringsAsFactors=FALSE)
uniprot <- uniprot[, c("Entry","Sequence","Length")]
colnames(uniprot) <- c("Uniprot_ID","sequence","Length")

df <- merge(df, uniprot, by="Uniprot_ID", all.x=TRUE)

############################################
# 🧬 5️⃣ FEATURE ENGINEERING
############################################
feature_fun <- function(seq){
  if(is.na(seq) || nchar(seq) < 10) return(rep(NA, 27))
  
  basic <- c(
    mw(seq), boman(seq), charge(seq),
    hydrophobicity(seq), instaIndex(seq),
    nchar(seq), pI(seq)
  )
  
  aa <- aaComp(seq)
  
  c(basic, aa)
}

features <- t(sapply(df$sequence, feature_fun))
features <- as.data.frame(features)

colnames(features) <- c(
  "MW","Boman","Charge","Hydro","Instability","Length","pI",
  paste0("AA_", names(aaComp("ACDEFGHIKLMNPQRSTVWY")))
)

df <- cbind(df, features)

############################################
# 🔗 6️⃣ PPI FEATURES
############################################
get_ppi <- function(id){
  url <- paste0("https://string-db.org/api/json/network?identifiers=", id)
  res <- try(httr::GET(url), silent=TRUE)
  
  if(inherits(res,"try-error") || httr::status_code(res)!=200) return(c(0,0))
  
  data <- try(jsonlite::fromJSON(httr::content(res,"text")), silent=TRUE)
  if(inherits(data,"try-error") || nrow(data)==0) return(c(0,0))
  
  c(nrow(data), mean(data$score))
}

ppi <- t(sapply(df$Uniprot_ID, get_ppi))
ppi <- as.data.frame(ppi)
colnames(ppi) <- c("PPI_Count","PPI_Score")

df <- cbind(df, ppi)

df$PPI_Log <- log1p(df$PPI_Count)
df$PPI_Strength <- df$PPI_Count * df$PPI_Score

############################################
# 🧬 7️⃣ INTERACTIONS
############################################
df$Hydro_Length <- df$Hydro * df$Length
df$Charge_pI <- df$Charge * df$pI

############################################
# 🧹 8️⃣ CLEAN DATA
############################################
df <- df[!is.na(df$sequence), ]

num_cols <- names(df)[sapply(df, is.numeric)]
for(col in num_cols){
  df[[col]][is.na(df[[col]])] <- median(df[[col]], na.rm=TRUE)
}

############################################
# 🧠 9️⃣ FEATURE MATRIX
############################################
num_cols <- names(df)[sapply(df, is.numeric)]
num_cols <- setdiff(num_cols, "label")

X <- df[, num_cols]

X <- data.frame(lapply(X, function(x) as.numeric(as.character(x))))
X[is.na(X)] <- 0

X <- as.data.frame(scale(X))
y <- df$label

############################################
# 🧠 🔟 TRAIN TEST SPLIT
############################################
set.seed(42)
train_idx <- caret::createDataPartition(y, p=0.8, list=FALSE)

trainX <- X[train_idx, ]
testX  <- X[-train_idx, ]

trainY <- y[train_idx]
testY  <- y[-train_idx]

############################################
# 🧠 1️⃣1️⃣ FEATURE SELECTION (FIXED)
############################################
rf_model <- randomForest::randomForest(
  x = trainX,
  y = as.factor(trainY),
  ntree = 200,
  importance = TRUE
)

imp <- importance(rf_model)

top_features <- rownames(imp)[order(imp[,1], decreasing=TRUE)[1:15]]

trainX <- trainX[, top_features]
testX  <- testX[, top_features]

############################################
# ⚖️ 1️⃣2️⃣ BALANCE
############################################
train_bal <- caret::upSample(
  x = trainX,
  y = as.factor(trainY)
)

trainX <- train_bal[, -ncol(train_bal)]
trainY <- as.numeric(train_bal$Class) - 1

############################################
# 🧠 1️⃣3️⃣ XGBOOST (TUNED)
############################################
dtrain <- xgb.DMatrix(as.matrix(trainX), label=trainY)
dtest  <- xgb.DMatrix(as.matrix(testX), label=testY)

model <- xgb.train(
  params = list(
    objective="binary:logistic",
    eval_metric="auc",
    max_depth=5,
    eta=0.03,
    subsample=0.9,
    colsample_bytree=0.9
  ),
  data = dtrain,
  nrounds = 400,
  verbose = 0
)

pred_prob <- predict(model, dtest)
pred_class <- ifelse(pred_prob > 0.5, 1, 0)

############################################
# 📊 1️⃣4️⃣ EVALUATION
############################################
roc_obj <- pROC::roc(testY, pred_prob)
auc_val <- pROC::auc(roc_obj)

conf <- caret::confusionMatrix(
  as.factor(pred_class),
  as.factor(testY)
)

print(auc_val)
print(conf$overall["Accuracy"])

############################################
# 📊 1️⃣5️⃣ PLOTS
############################################
png("plots/ROC.png")
plot(roc_obj, col="blue", lwd=3,
     main=paste("AUC =", round(auc_val,3)))
abline(a=0,b=1,lty=2,col="red")
dev.off()

############################################
# 💾 1️⃣6️⃣ SAVE
############################################
write.csv(data.frame(testY, pred_prob, pred_class),
          "results/final_predictions.csv",
          row.names=FALSE)

############################################
# 🧠 MULTI-MODEL TRAINING
############################################
library(caret)
library(xgboost)
library(randomForest)

set.seed(42)

control <- trainControl(
  method = "repeatedcv",
  number = 5,
  repeats = 3,
  classProbs = TRUE,
  summaryFunction = twoClassSummary
)

# Convert labels for caret
trainY_factor <- factor(ifelse(trainY==1,"Class1","Class0"))
testY_factor  <- factor(ifelse(testY==1,"Class1","Class0"))

############################################
# 🌲 RANDOM FOREST
############################################
rf_model <- train(
  x = trainX,
  y = trainY_factor,
  method = "rf",
  metric = "ROC",
  trControl = control,
  ntree = 300
)

############################################
# ⚡ XGBOOST
############################################
xgb_model <- train(
  x = trainX,
  y = trainY_factor,
  method = "xgbTree",
  metric = "ROC",
  trControl = control,
  tuneLength = 5
)

############################################
# 📈 LOGISTIC REGRESSION
############################################
glm_model <- train(
  x = trainX,
  y = trainY_factor,
  method = "glm",
  family = "binomial",
  metric = "ROC",
  trControl = control
)


############################################
# 📊 PREDICTIONS
############################################
rf_pred  <- predict(rf_model, testX, type="prob")[,2]
xgb_pred <- predict(xgb_model, testX, type="prob")[,2]
glm_pred <- predict(glm_model, testX, type="prob")[,2]


############################################
# 📊 ROC COMPARISON
############################################
library(pROC)

roc_rf  <- roc(testY, rf_pred)
roc_xgb <- roc(testY, xgb_pred)
roc_glm <- roc(testY, glm_pred)

png("plots/ROC_Comparison.png", width=1200, height=900)

plot(roc_xgb, col="red", lwd=3, main="Model Comparison ROC")
plot(roc_rf, col="blue", lwd=3, add=TRUE)
plot(roc_glm, col="green", lwd=3, add=TRUE)

legend("bottomright",
       legend=c(
         paste("XGB AUC:", round(auc(roc_xgb),3)),
         paste("RF AUC:", round(auc(roc_rf),3)),
         paste("GLM AUC:", round(auc(roc_glm),3))
       ),
       col=c("red","blue","green"), lwd=3)

dev.off()

cat("✅ ROC comparison saved\n")


############################################
# 🧠 BEST MODEL
############################################
auc_values <- c(
  XGB = auc(roc_xgb),
  RF  = auc(roc_rf),
  GLM = auc(roc_glm)
)

print(auc_values)

best_model_name <- names(which.max(auc_values))
cat("🏆 Best Model:", best_model_name, "\n")


############################################
# 🔍 SHAP INTERPRETATION
############################################
library(SHAPforxgboost)

if(best_model_name == "XGB"){
  
  final_model <- xgb_model$finalModel
  
  shap <- shap.values(
    xgb_model = final_model,
    X_train = as.matrix(trainX)
  )
  
  shap_long <- shap.prep(
    shap_contrib = shap$shap_score,
    X_train = trainX
  )
  
  png("plots/SHAP.png", width=1200, height=800)
  shap.plot.summary(shap_long)
  dev.off()
  
  cat("✅ SHAP generated for XGBoost\n")
}

############################################
# 📊 FINAL METRICS
############################################
best_pred <- switch(best_model_name,
                    "XGB" = xgb_pred,
                    "RF"  = rf_pred,
                    "GLM" = glm_pred)

best_class <- ifelse(best_pred > 0.5,1,0)

conf <- confusionMatrix(
  as.factor(best_class),
  as.factor(testY)
)

print(conf$overall)
print(conf$byClass)


############################################
# 📊 CONFUSION MATRIX (PLOT)
############################################
library(ggplot2)
library(caret)

conf <- confusionMatrix(as.factor(pred_class), as.factor(testY))

cm <- as.data.frame(conf$table)

png("plots/confusion_matrix.png", width=800, height=600)

ggplot(cm, aes(Prediction, Reference, fill=Freq)) +
  geom_tile() +
  geom_text(aes(label=Freq), color="white", size=6) +
  scale_fill_gradient(low="blue", high="red") +
  theme_minimal() +
  labs(title="Confusion Matrix")

dev.off()


############################################
# 📊 ROC CURVE
############################################
library(pROC)

roc_obj <- roc(testY, pred_prob)
auc_val <- auc(roc_obj)

png("plots/ROC.png")

plot(roc_obj, col="blue", lwd=3,
     main=paste("ROC Curve (AUC =", round(auc_val,3),")"))
abline(a=0,b=1,lty=2,col="red")

dev.off()


############################################
# 📊 ROC WITH METRICS
############################################
acc <- conf$overall["Accuracy"]

png("plots/ROC_with_metrics.png", width=900, height=700)

plot(roc_obj, col="darkgreen", lwd=3,
     main=paste(
       "AUC =", round(auc_val,3),
       "| Accuracy =", round(acc,3)
     ))

abline(a=0,b=1,lty=2,col="red")

dev.off()


############################################
# 📊 PRECISION-RECALL CURVE
############################################
library(PRROC)

pr <- pr.curve(
  scores.class0 = pred_prob[testY==1],
  scores.class1 = pred_prob[testY==0],
  curve = TRUE
)

png("plots/PR_curve.png")

plot(pr, main="Precision-Recall Curve")

dev.off()


############################################
# 📊 HISTOGRAM
############################################
png("plots/histogram.png")

ggplot(data.frame(pred_prob), aes(x=pred_prob)) +
  geom_histogram(bins=30, fill="steelblue", color="black") +
  theme_minimal() +
  labs(title="Prediction Probability Distribution")

dev.off()


############################################
# 📊 BOXPLOT
############################################
png("plots/boxplot.png")

ggplot(data.frame(Class=factor(testY), Prob=pred_prob),
       aes(x=Class, y=Prob)) +
  geom_boxplot(fill="orange") +
  theme_minimal() +
  labs(title="Prediction vs True Class")

dev.off()


############################################
# 📊 FEATURE IMPORTANCE
############################################
library(xgboost)

importance_matrix <- xgb.importance(model = model)

png("plots/feature_importance.png", width=1000, height=800)

xgb.plot.importance(importance_matrix)

dev.off()


############################################
# 📊 BAR PLOT (TOP FEATURES)
############################################
imp_df <- as.data.frame(importance_matrix)

png("plots/feature_bar.png", width=900, height=700)

ggplot(head(imp_df,10),
       aes(x=reorder(Feature, Gain), y=Gain)) +
  geom_bar(stat="identity", fill="purple") +
  coord_flip() +
  theme_minimal() +
  labs(title="Top Features")

dev.off()


############################################
# 📊 SHAP
############################################
library(SHAPforxgboost)

X_train_df <- as.data.frame(trainX)

shap <- shap.values(
  xgb_model = model,
  X_train = as.matrix(X_train_df)
)

shap_long <- shap.prep(
  shap_contrib = shap$shap_score,
  X_train = X_train_df
)

png("plots/SHAP.png", width=1200, height=800)

shap.plot.summary(shap_long)

dev.off()
cat("✅ All models trained\n")
cat("🎯 PIPELINE COMPLETED SUCCESSFULLY\n")