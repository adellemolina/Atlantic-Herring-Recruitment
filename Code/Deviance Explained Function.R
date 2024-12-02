

# basically its (total deviance - cv residual deviance)/total
print(gbm1)

# Hard coded versioin
(((BRT.ldev_wlags$self.statistics$mean.null)-(BRT.ldev_wlags$cv.statistics$deviance.mean))/(BRT.ldev_wlags$self.statistics$mean.null))*100 
(((brt.model$self.statistics$mean.null)-(brt.model$cv.statistics$deviance.mean))/(brt.model$self.statistics$mean.null))*100 


# Function to calculate deviance explained for BRT models created using gbm.sterp
dev.exp <- function(brt.model){
  (((brt.model$self.statistics$mean.null)-(brt.model$cv.statistics$deviance.mean))/(brt.model$self.statistics$mean.null))*100 
}
diff = model


dev.exp2 <- function(brt.model){
  abs(mean((brt.model$train.error-brt.model$cv_error)/brt.model$train.error))*100
}

dev.exp3 <- function(brt.model){
  abs(mean((brt.model$valid.error-brt.model$cv_error)/brt.model$valid.error))*100
}
