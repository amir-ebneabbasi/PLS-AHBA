# Code to produce permuted indices using coordinates (vasa and hungarian)
# Sourced from frantisekvasa

rotate.parcellation = function(coord.l, coord.r, nrot=1000, method ='hungarian') {

  library(matrixStats)
  library(clue)
  
  if (!all(dim(coord.l)[2]==3,dim(coord.r)[2]==3)) {
    if (all(dim(coord.l)[1]==3,dim(coord.r)[1]==3)) {
      print('transposing coordinates to be of dimension nROI x 3')
      coord.l = t(coord.l)
      coord.r = t(coord.r)
    }
  }
  
  nroi.l = dim(coord.l)[1]  
  nroi.r = dim(coord.r)[1]   
  nroi = nroi.l+nroi.r     
  
  perm.id = array(0,dim=c(nroi,nrot)); 
  r = 0; c = 0; 
  
  I1 = diag(3); I1[1,1] = -1;

    while (r < nrot) {
    
    A = matrix(rnorm(9, mean = 0, sd = 1), nrow = 3, ncol = 3)
    qrdec = qr(A)      
    TL = qr.Q(qrdec)   
    temp = qr.R(qrdec)  
    TL = TL%*%diag(sign(diag(temp)))
    if (det(TL)<0) {
      TL[,1] = -TL[,1]
    }
    TR = I1 %*% TL %*% I1;
    coord.l.rot = coord.l %*% TL;
    coord.r.rot = coord.r %*% TR; 
    
    dist.l = array(0,dim=c(nroi.l,nroi.l));
    dist.r = array(0,dim=c(nroi.r,nroi.r));
    
    for (i in 1:nroi.l) {
      for (j in 1:nroi.l) {
        dist.l[i,j] = sqrt( sum( (coord.l[i,]-coord.l.rot[j,])^2 ) )
      }
    }
    for (i in 1:nroi.r) { 
      for (j in 1:nroi.r) {
        dist.r[i,j] = sqrt( sum( (coord.r[i,]-coord.r.rot[j,])^2 ) )
      }
    }
    
    if (method == 'vasa') {
      
      temp.dist.l = dist.l
      rot.l = c(); ref.l = c();
      for (i in 1:nroi.l) {
        ref.ix = which( rowMins(temp.dist.l,na.rm=T) == max(rowMins(temp.dist.l,na.rm=T),na.rm=T) )   
        rot.ix = which( temp.dist.l[ref.ix,] == min(temp.dist.l[ref.ix,],na.rm=T) )
        
        ref.l = c(ref.l,ref.ix) 
        rot.l = c(rot.l,rot.ix)
        
        temp.dist.l[,rot.ix] = NA 
        temp.dist.l[ref.ix,] = 0
      }
      
      temp.dist.r = dist.r;
      rot.r = c(); ref.r = c();
      for (i in 1:nroi.r) {
        ref.ix = which( rowMins(temp.dist.r,na.rm=T) == max(rowMins(temp.dist.r,na.rm=T),na.rm=T) )   
        rot.ix = which( temp.dist.r[ref.ix,] == min(temp.dist.r[ref.ix,],na.rm=T) )            
        
        ref.r = c(ref.r,ref.ix) 
        rot.r = c(rot.r,rot.ix)
        
        temp.dist.r[,rot.ix] = NA 
        temp.dist.r[ref.ix,] = 0       
      }
      
    } else if (method == 'hungarian') {
      
      rot.l = as.vector(solve_LSAP(dist.l, maximum=FALSE))
      ref.l = c(1:nroi.l) 
      
      rot.r = as.vector(solve_LSAP(dist.r, maximum=FALSE))
      ref.r = c(1:nroi.r) 
      
    } else {
      
      stop(paste('\'',method,'\' is not an accepted permutation method; valid options are \'vasa\' or \'hungarian\'',sep=''))
      
    }
    
    ref.lr = c(ref.l,nroi.l+ref.r); rot.lr = c(rot.l,nroi.l+rot.r);
    b = sort(ref.lr,index.return=T); 
    ref.lr.sort = ref.lr[b$ix]; rot.lr.sort = rot.lr[b$ix];
    
    if (!all(sort(rot.lr.sort,decreasing=F)==c(1:nroi))) {
      browser("permutation error")
    }
    
    if (!all(rot.lr.sort==c(1:nroi))) {
      r = r+1
      perm.id[,r] = rot.lr.sort 
    } else {
      c = c+1
      print(paste('map to itself n. ',toString(c),sep=''))
    }
    
    if (r%%10==0) print(paste('permutation ',toString(r),' of ',toString(nrot),sep=''))
    
  }
  
  return(perm.id)

}

# How to run
# Load the coordinates from the CSV file
# data <- read.csv("path/to/HCP_sphere.csv", header = FALSE)

# Split the data into left and right hemisphere coordinates
# coord.l <- data[1:34, ]  # First 34 rows for the left hemisphere
# coord.r <- data[35:68, ] # Next 34 rows for the right hemisphere

# coord.l <- data[1:180, ]  # First 180 rows for the left hemisphere
# coord.r <- data[181:360, ] # Next 180 rows for the right hemisphere

# Ensure the data is in matrix form
# coord.l <- as.matrix(coord.l)
# coord.r <- as.matrix(coord.r)

# Call the rotate.parcellation function
# perm.id <- rotate.parcellation(coord.l, coord.r)

# save permutation indices
# write.csv(perm.id, "path/to/output/HCP_hungarian.csv", row.names = FALSE)
