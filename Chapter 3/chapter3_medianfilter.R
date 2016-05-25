##
## Copyright 2016 Simon Chapple
##
## Packt: "Mastering Parallelism with R"
## Chapter 3 - Advanced MPI Grid Parallelism Median Filter
##
library(Rmpi)

# Useful constants
Height<-200; Width<-200; # Size of image
Dim<-2; # Square size of grid
N<-1; NE<-2; E<-3; SE<-4;
S<-5; SW<-6; W<-7; NW<-8;

worker_makeSquareGrid <- function(dim,comm)
{
  print(paste0("Base grid comm=",comm," dim=",dim)) 
  grid <- 1000 + dim    # assign comm handle for this size grid
  dims <- c(dim,dim)    # dimensions are 2D, size: dim X dim
  periods <- c(FALSE,FALSE)  # no wraparound at outermost edges
  if (mpi.cart.create(commold=comm,dims,periods,commcart=grid))  
  {
    return(grid)
  }
  return(-1) # An MPI error occurred
}

worker_initSpatialGrid <- function(dim,comm=Wcomm) {
  Gcomm <- worker_makeSquareGrid(dim,comm)
  myRank <- mpi.comm.rank(Gcomm)
  myUniverseRank <- mpi.comm.rank(1) # Lookup rank in cluster
  print(paste("myRank:",myRank))
  myCoords <- mpi.cart.coords(Gcomm,myRank,2)
  print(paste("myCoords:",myCoords))
  # (y^,x>) co-ordinate system
  myY <- myCoords[1]; myX <- myCoords[2]; 
  coords <- vector(mode="list", length=8)
  neighbors <- rep(-1,8)
  if (myY+1 < dim) { 
    neighbors[N] <- mpi.cart.rank(Gcomm,c(myY+1,myX))
  }
  if (myX+1 < dim && myY+1 < dim) { 
    neighbors[NE] <- mpi.cart.rank(Gcomm,c(myY+1,myX+1))
  }
  if (myX+1 < dim) { 
    neighbors[E] <- mpi.cart.rank(Gcomm,c(myY,myX+1))
  }
  if (myX+1 < dim && myY-1 >= 0) { 
    neighbors[SE] <- mpi.cart.rank(Gcomm,c(myY-1,myX+1))
  }
  if (myY-1 >= 0) {
    neighbors[S] <- mpi.cart.rank(Gcomm,c(myY-1,myX))
  }
  if (myX-1 >= 0 && myY-1 >= 0) { 
    neighbors[SW] <- mpi.cart.rank(Gcomm,c(myY-1,myX-1))
  }
  if (myX-1 >= 0) {
    neighbors[W] <- mpi.cart.rank(Gcomm,c(myY,myX-1))
  }
  if (myX-1 >= 0 && myY+1 < dim) { 
    neighbors[NW] <- mpi.cart.rank(Gcomm,c(myY+1,myX-1))
  }
  # Store reference for neighbor comms
  assign("Neighbors", neighbors, envir=.GlobalEnv)
  # Store reference for grid communicator
  assign("Gcomm", Gcomm, envir=.GlobalEnv)
  return(list(myY,myX,myUniverseRank))
}

worker_boundaryExchange <- function(img,neighbors,comm) {
  # More efficient to set-up non-blocking receives then sends
  neighbors <- Neighbors; comm <- Gcomm;
  print(comm); print(neighbors)
  
  # Set-up non-blocking receives for incoming boundary data
  # Local image tile has one pixel shared border
  len <- ncol(img)-2 
  rbuf <- vector(mode="list", length=8) # 8 receive buffers
  req <- 0
  for (i in 1:8) {
    if (neighbors[i]>=0) {
      rbuf[[i]] <- integer(length=len)
      tag <- mpi.any.tag()
      print(rbuf[[i]])
      mpi.irecv(rbuf[[i]],1,neighbors[i],tag,comm=comm,request=req)
      req <- req + 1
    }
  }

  edge <- ncol(img)-1 # image is square: ncol=nrow
  sbuf <- vector(mode="list", length=8) # 8 send buffers
  # non-block send my tile data boundaries to my neighbours
  if (neighbors[N]>=0) { # north
  	sbuf[[N]] <- img[2,2:edge]
    mpi.isend(sbuf[[N]],1,neighbors[N],N,comm=comm,request=req)
    req <- req + 1
  }
  if (neighbors[NE]>=0) { # ne
    sbuf[NE] <- img[2,edge] # top-right inner cell
    mpi.isend(sbuf[[NE]],1,neighbors[NE],NE,comm=comm,request=req)
    req <- req + 1
  }
  if (neighbors[E]>=0) { # east
    sbuf[[E]] <- img[2:edge,edge] # rightmost inner col
    mpi.isend(sbuf[[E]],1,neighbors[E],E,comm=comm,request=req)
    req <- req + 1
  }
  if (neighbors[SE]>=0) { # se
    sbuf[[SE]] <- img[edge,edge] # bottom-right inner cell
    mpi.isend(sbuf[[SE]],1,neighbors[SE],SE,comm=comm,request=req)
    req <- req + 1
  }
  if (neighbors[S]>=0) { # south
    sbuf[[S]] <- img[edge,2:edge] # bottom inner row
    mpi.isend(sbuf[[S]],1,neighbors[S],S,comm=comm,request=req)
    req <- req + 1
  }
  if (neighbors[SW]>=0) { # sw
    sbuf[[SW]] <- img[edge,2] # bottom-left inner cell
    mpi.isend(sbuf[[SW]],1,neighbors[SW],SW,comm=comm,request=req)
    req <- req + 1
  }
  if (neighbors[W]>=0) { # west
    sbuf[[W]] <- img[2:edge,2] # leftmost inner col
    mpi.isend(sbuf[[W]],1,neighbors[W],W,comm=comm,request=req)
    req <- req + 1
  }
  if (neighbors[NW]>=0) { # nw
    sbuf[[NW]] <- img[2,2] # top-left inner cell
    mpi.isend(sbuf[[NW]],1,neighbors[NW],NW,comm=comm,request=req)
    req <- req + 1
  }
  
  mpi.waitall(req) # Wait for all boundary comms to complete

  # Unpack received boundary data from neighbors into my image tile
  n <- ncol(img)
  if (neighbors[N]>=0) { # north
    img[1,2:edge] <- rbuf[[N]] # top row
  }
  if (neighbors[NE]>=0) { # ne
    img[1,n] <- rbuf[[NE]][1] # top-right cell
  }
  if (neighbors[E]>=0) { # east
    img[2:edge,n] <- rbuf[[E]] # rightmost column
  }
  if (neighbors[SE]>=0) { # se
    img[n,n] <- rbuf[[SE]][1] # bottom-right cell
  }
  if (neighbors[S]>=0) { # south
    img[n,2:edge] <- rbuf[[S]] # bottom row
  }
  if (neighbors[SW]>=0) { # sw
    img[n,1] <- rbuf[[SW]][1] # bottom-left cell
  }
  if (neighbors[W]>=0) { # west
    img[2:edge,1] <- rbuf[[W]] # leftmost column
  }
  if (neighbors[NW]>=0) { # nw
    img[1,1] <- rbuf[[NW]][1] # top-left cell
  }
  return(img)
}

medianFilterPixel3 <- function(y,x,img) {
  v <- vector("integer",9) # bottom-left to top-right
  v[1]<-img[y-1,x-1]; v[2]<-img[y-1,x]; v[3]<-img[y-1,x+1];
  v[4]<-img[y,  x-1]; v[5]<-img[y,  x]; v[6]<-img[y,  x+1];
  v[7]<-img[y+1,x-1]; v[8]<-img[y+1,x]; v[9]<-img[y+1,x+1];
  s <- sort(v); # sort by pixel value (default ascending)
  return (s[5]) # return the middle value of the nine
}

worker_gridApplyMedianFilter <- function(niters) {
  # Receive tile from Master on Rmpi default comm
  tile <- mpi.recv.Robj(0,1,comm=1,status=1)
  
  # Create local image with extra pixel boundary
  theight <- nrow(tile); iheight <- theight+2;
  twidth <- ncol(tile); iwidth <- twidth+2;
  print(paste("Received tile:",theight,twidth))
  img <- matrix(0L,nrow=iheight,ncol=iwidth)
  
  # Initialize borders with out-of-bound pixel values
  # These values will be sorted to the ends of the set of 9
  # and so will not interfere with the real image values
  img[1,1:iwidth] <- rep(c(-1,256),times=iwidth/2)
  img[1:iheight,1] <- rep(c(-1,256),times=iheight/2)
  img[iheight,1:iwidth] <- rep(c(256,-1),times=iwidth/2)
  img[1:iheight,iwidth] <- rep(c(256,-1),times=iheight/2)
  
  # Set internal bounded area to the received tile
  img[2:(theight+1),2:(twidth+1)] <- tile
  
  # Apply multi-pass image operation
  for (i in 1:niters) {
  	print(paste("Iteration",i))
    img <- worker_boundaryExchange(img)
    for (y in 2:theight+1) {
      for (x in 2:twidth+1) {
        img[y,x] <- medianFilterPixel3(y,x,img)
      }
    }
  }

  # Send processed tile to Master on default comm
  tile <- img[2:(theight+1),2:(twidth+1)]
  mpi.send.Robj(tile,0,2,comm=1)
}

############################################################
# Master co-ordinates creation and operation of the grid,
# but does not itself participate in any tile computation.

# Launch the Rmpi based grid with (dimXdim) worker processes
dim <- Dim;
np <- dim * dim # number of MPI processes in grid
mpi.spawn.Rslaves(
  Rscript=system.file("workerdaemon.R", package="Rmpi"),
  nslaves=np)

# Send all Master defined globals/functions to Workers
mpi.bcast.Robj2slave(all=TRUE) 

# Map grid co-ords to cluster rank assignment of the Workers
map <- mpi.remote.exec(worker_initSpatialGrid(),dim,
                       simplify=FALSE,comm=1)
workerRanks <- matrix(-1,nrow=dim,ncol=dim)
for (p in 1:length(map)) {
  y <- map[[p]][[1]]
  x <- map[[p]][[2]]
  rank <- map[[p]][[3]]
  print(paste0("Map ",p,": (",y,",",x,") => ",rank))
  workerRanks[y+1,x+1] <- rank
}

# We create large B/W image array with values in the range 101-111
height <- Height; width <- Width;
image1 <- matrix(sample(101:111,height*width,replace=TRUE),
                 height,width)
# We add a bit of white saturation noise (pixel value=255)
image1[height/6,width/6] <- 255
image1[height/5,width/5] <- 255
image1[height/4,width/4] <- 255
image1[height/3,width/3] <- 255
image1[height/2.1,width/2.1] <- 255
image1[height/1.1,width/1.1] <- 255
image1[height/1.2,width/1.2] <- 255
image1[height/1.3,width/1.3] <- 255
image1[height/1.4,width/1.4] <- 255
image1[height/1.5,width/1.5] <- 255

# Tell the workers to process the image (3 passes of the medianFilter)
# The Workers first wait to receive their local tile from the Master,
# will do their multi-pass image processing, then finally send their
# processed tiles back to the Master.
mpi.bcast.cmd(worker_gridApplyMedianFilter(3))
Start <- proc.time()

# We split the image into non-overlapping square grid tiles
# and distribute one per Worker
twidth <- width/dim # tile width
theight <- height/dim # tile height
for (ty in 0:(dim-1)) { # bottom-left to top-right
  sy <- (ty * theight) +1
  for (tx in 0:(dim-1)) {
    sx <- (tx * twidth) +1
    tile <- image1[sy:(sy+theight-1),sx:(sx+twidth-1)]
    # Send tile to the appropriate Worker
    worker <- workerRanks[ty+1,tx+1]
    mpi.send.Robj(tile,worker,1,comm=1)
    print(paste0("Sent tile to ", worker,
          " y=",sy,"-",sy+theight-1," x=",sx,"-",sx+twidth-1))
  }
}

# Create processed output image, initially blank
image2 <- matrix(0L,nrow=height,ncol=width)

# Master receives output tiles in sequence and unpacks
# each into its correct place to form the output image
for (ty in 0:(dim-1)) { # bottom-left to top-right
  sy <- (ty * theight) +1
  for (tx in 0:(dim-1)) {
    sx <- (tx * twidth) +1
    # Receive tile from the appropriate Worker
    worker <- workerRanks[ty+1,tx+1]
    tile <- mpi.recv.Robj(worker,2,comm=1)
    print(paste0("Received tile from ", worker,
          " y=",sy,"-",sy+theight-1," x=",sx,"-",sx+twidth-1))
    image2[sy:(sy+theight-1),sx:(sx+twidth-1)] <- tile
  }
}

# Ta da!
Finish <- proc.time()
print(paste("Image size:",Height,"x",Width," processed with",np,"Workers in",Finish[3]-Start[3],"elapsed seconds"))
print(paste("Noisy image max pixel value",max(image1))) # Saturated image=255
print(paste("Clean image max pixel value",max(image2))) # MedianFiltered image=111
mpi.close.Rslaves()
