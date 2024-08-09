args <- as.numeric(commandArgs(TRUE))
idx_reg <- args[1] ## region index

library(data.table)

########## identify subjects & regions ########## 
sub_id <- read.table("/data/brutus_data34/Hongjie/ACP_subj_id.txt")
sub_id <- as.vector(unlist(sub_id))
reg_id <- read.table("/data/brutus_data34/Hongjie/ACP_region_id.txt")
reg_id <- as.vector(unlist(reg_id))

########## identify working region ########## 
print(paste("This is region", reg_id[idx_reg]))
out_path <- paste("/data/brutus_data34/Hongjie/ACP_sp_feature/", reg_id[idx_reg], "/",sep="")
# dir.create(out_path)

all_coor <- read.table(paste("/data/brutus_data34/Hongjie/ACP_rsTS_by_voxel/", 
                             reg_id[idx_reg], "/", "voxel_name.txt", sep=""))
all_coor <- as.vector(unlist(all_coor))
print(paste("Number of voxels is", length(all_coor)))

########## extract feature ########## 
freq_use <- read.table("/data/brutus_data32/Hongjie/ACP_R199_power/frequency_vec.txt")
freq_use <- freq_use[,1] ## frequency of fft results 

midpt <- c(0.0400641, 0.120192, 0.200321, 0.280449, 
           0.360577, 0.440705, 0.520833, 0.600962)
bin_size <- midpt[4]-midpt[3]
bb <- matrix(c(midpt-bin_size/2, midpt+bin_size/2), ncol=2) ## alternative

t1 <- Sys.time()
for (j in 1:length(all_coor)) {
  if (file.exists(paste(out_path, all_coor[j], ".txt", sep="")) == FALSE) {
    print(j)
    pow_res <- fread(paste("/data/brutus_data34/Hongjie/ACP_power/", reg_id[idx_reg],
                           "/", all_coor[j], ".txt", sep=""))
    pow_res <- as.data.frame(pow_res)
    
    ## binning
    sp_ft <- matrix(NA, nrow(pow_res), nrow(bb))
    for (i in 1:nrow(bb)) {
      pos <- which(freq_use>=bb[i,1] & freq_use<=bb[i,2])
      sp_ft[,i] <- apply(pow_res[,pos], 1, mean)
    }
    
    write.table(sp_ft, col.names=F, row.names=F, quote=F, 
                file=paste(out_path, all_coor[j], ".txt", sep=""))
  }
}
t2 <- Sys.time()
print(t2-t1)

write.table(NULL, file=paste(out_path, "all_voxels_finished.txt", sep=""))

