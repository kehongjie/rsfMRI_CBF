args <- as.numeric(commandArgs(TRUE))
idx_reg <- args[1] ## region indicator 
st <- args[2] ## starting voxel indicator
ed <- args[3] ## end voxel indicator 
print(paste("Voxel from", st, "to", ed))

library(data.table)
library(RANN)
library(e1071)
t1 <- Sys.time()

########## identify subjects & regions ########## 
sub_id <- read.table("/data/brutus_data34/Hongjie/ACP_subj_id.txt")
sub_id <- as.vector(unlist(sub_id))
reg_id <- read.table("/data/brutus_data34/Hongjie/ACP_region_id.txt")
reg_id <- as.vector(unlist(reg_id))

print(paste("This is region", reg_id[idx_reg]))
out_path1 <- paste("/data/brutus_data34/Hongjie/ACP_fit_object/", reg_id[idx_reg], "/",sep="")
out_path2 <- paste("/data/brutus_data34/Hongjie/ACP_prediction_result/", reg_id[idx_reg], "/",sep="")
dir.create(out_path1)
dir.create(out_path2)

vox_name <- read.table(paste("/data/brutus_data34/Hongjie/ACP_rsTS_by_voxel/", 
                             reg_id[idx_reg], "/", "voxel_name.txt", sep=""))
vox_name <- as.vector(unlist(vox_name))
print(paste("Number of voxels is", length(vox_name)))


########## other info ########## 
k <- 6 ## number of neighbors 

## get (x,y,z)
folder <- "/data/brutus_data9/BMA/CortRois_rsTS_cbf_ACP/" ## path for all regions
rts <- fread(file=paste(folder, sub_id[1], "/Ctx2_", reg_id[idx_reg], "_rsTS.txt", sep="")) ## other regions
data_coor <- rts[,1:3] ## (x,y,z) coords

## get GM partial occupancy data
all_ocp <- read.table(paste("/data/brutus_data34/Hongjie/ACP_GM_ocp/", 
                            reg_id[idx_reg], ".txt", sep=""), header=F)

## get CBF data
all_cbf <- read.table(paste("/data/brutus_data34/Hongjie/ACP_rsTS_by_voxel/", 
                            reg_id[idx_reg], "/", "CBF_subj_voxel.txt", sep=""))


## trainig/testing split
set.seed(123)
pos_tr <- sample(1:length(sub_id), floor(0.7*length(sub_id))) ## 70% training and 30 % testing
all_gamma <- c(0.05,0.1,0.5,1,2,5,10) ## you might want to tune these a bit
all_cost <- c(0.01,0.1,0.5,1,5,10,20,30,50,100)


# i <- 1
for (i in st:ed) { ## i is index for voxel, 1:length(vox_name)
  if (file.exists(paste(out_path1, vox_name[i], ".RData", sep="")) == FALSE) {
    print(paste(i, "-th voxel", sep=""))
    
    ## find neighbors
    nn <- nn2(data_coor, data_coor[i,], k+1) 
    # print(nn$nn.idx) ## positions of neighbors
    pos_nb <- as.numeric(nn$nn.idx)
    
    ## feature data from fft
    data_x <- NULL
    for (j in 1:(k+1)) {
      temp <- read.table(paste("/data/brutus_data34/Hongjie/ACP_sp_feature/", 
                               reg_id[idx_reg], "/", 
                               vox_name[pos_nb[j]],".txt",sep=""))
      data_x <- cbind(data_x, as.matrix(temp))
    }
    colnames(data_x) <- paste("feature", 1:((k+1)*8), sep="")
    
    ## combine x and y 
    data_vox <- cbind(cbf=all_cbf[,i], data_x)
    
    ## combine with partial occupancy data 
    data_ocp <- cbind(data_vox, ocp=all_ocp[,i])
    
    ## split training and testing
    data_tr1 <- data_vox[pos_tr, ] ## training set, only power
    data_te1 <- data_vox[-pos_tr, ] ## testing set, only power
    data_tr2 <- data_ocp[pos_tr, ] ## training set, w/ occupancy
    data_te2 <- data_ocp[-pos_tr, ] ## testing set, w/ occupancy
    
    pre_res <- matrix(NA, nrow(data_vox), 4) ## prediction results
    colnames(pre_res) <- c("training", "true", "model1", "model2")
    pre_res[pos_tr,1] <- 1 ## training set indicator
    pre_res[-pos_tr,1] <- 0
    pre_res[,2] <- data_vox[,"cbf"] ## true cbf values
    
    ## model w/o partial occupancy
    fit_cv1 <- tune.svm(cbf~., data=as.data.frame(data_tr1), gamma=all_gamma, 
                        cost=all_cost)
    fit1 <- svm(cbf~., data=data_tr1, gamma=fit_cv1$best.parameters[1], 
                cost=fit_cv1$best.parameters[2])
    pre_tr1 <- predict(fit1, data_tr1[,-1])
    pre_te1 <- predict(fit1, data_te1[,-1])
    
    ## model with partial occupancy
    fit_cv2 <- tune.svm(cbf~., data=as.data.frame(data_tr2), gamma=all_gamma,
                        cost=all_cost)
    fit2 <- svm(cbf~., data=data_tr2, gamma=fit_cv2$best.parameters[1],
                cost=fit_cv2$best.parameters[2])
    pre_tr2 <- predict(fit2, data_tr2[,-1])
    pre_te2 <- predict(fit2, data_te2[,-1])
    
    pre_res[pos_tr,3] <- pre_tr1
    pre_res[-pos_tr,3] <- pre_te1
    pre_res[pos_tr,4] <- pre_tr2
    pre_res[-pos_tr,4] <- pre_te2
    
    write.csv(pre_res, quote = F, row.names = F,
              file=paste(out_path2, vox_name[i], ".csv", sep=""))
    save(fit1, fit2, 
         file=paste(out_path1, vox_name[i], ".RData", sep=""))
  }
}

t2 <- Sys.time()
print(t2-t1)
print("All finished")

if (i == length(vox_name)) {
  write.table(NULL, file=paste(out_path1, "all_voxels_finished.txt", sep=""))
}



