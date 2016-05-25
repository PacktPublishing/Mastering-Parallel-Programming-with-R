##
## Copyright 2015 Simon Chapple
##
## Packt: "Mastering Parallelism with R"
## Chapter 1 - Aristotle's Number Puzzle
##

# teval() is a simple function to measure the execution time of some other function
teval <- function(...) {
  gc();
  start <- Sys.time()
  pstart <- proc.time()
  result <- eval(...)
  pfinish <- proc.time()
  finish <- Sys.time()
  wallclock = list(Start=format(start,"%H:%M:%OS3"),Finish=format(finish,"%H:%M:%OS3"))
  return ( list(PID=Sys.getpid(), WallClock=wallclock, Duration=pfinish-pstart, Result=result) )
}

# generateTriples() generates all possible combinations of initial three tile placements for the board
generateTriples <- function()
{
  triples <- list()
  for (x in 1:19) {
    for (y in 1:19) {
      if (y == x) next
      for (z in 1:19) {
        if (z == x || z == y || x+y+z != 38) next
        mirror <- FALSE
        reversed <- sprintf("%02d%02d%02d",z,y,x)
        for (t in triples) {
          if (reversed == t) {
            mirror <- TRUE
            break
          }
        }
        if (!mirror) {
          triples[length(triples)+1] <- sprintf("%02d%02d%02d",x,y,z)
        }
      }
    }
  }
  return (triples)
}

# solver(triple) finds the first puzzle solution based on given starting triple (if it exists)
solver <- function(triple)
{
  ##
  ## Some "global" variables for use in local functions
  ##
  doPrint      <- FALSE    # turn on/off any display output
  doProgress   <- FALSE    # turn on/off progress messages
  progressRate <- 10000    # report progress every nth tile placement
  tilesPlaced  <- 0        # accumulator for number of tile placements
  
  all_lines <- list(
    c(1,2,3),        c(1,4,8),        c(1,5,10,15,19),
    c(2,5,9,13),     c(2,6,11,16),    c(3,7,12),
    c(3,6,10,14,17), c(4,5,6,7),      c(4,9,14,18),
    c(7,11,15,18),   c(8,9,10,11,12), c(8,13,17),
    c(12,16,19),     c(13,14,15,16),  c(17,18,19)
  )
  
  cell_lines <- list (
    list( c(1,2,3),        c(1,4,8),        c(1,5,10,15,19) ), #Cell 1
    list( c(1,2,3),        c(2,5,9,13),     c(2,6,11,16)    ), #Cell 2
    list( c(1,2,3),        c(3,7,12),       c(3,6,10,14,17) ), #Cell 3
    list( c(4,5,6,7),      c(1,4,8),        c(4,9,14,18)    ), #Cell 4
    list( c(4,5,6,7),      c(2,5,9,13),     c(1,5,10,15,19) ), #Cell 5
    list( c(4,5,6,7),      c(2,6,11,16),    c(3,6,10,14,17) ), #Cell 6
    list( c(4,5,6,7),      c(3,7,12),       c(7,11,15,18)   ), #Cell 7
    list( c(1,4,8),        c(8,9,10,11,12), c(8,13,17)      ), #Cell 8
    list( c(2,5,9,13),     c(8,9,10,11,12), c(4,9,14,18)    ), #Cell 9
    list( c(8,9,10,11,12), c(3,6,10,14,17), c(1,5,10,15,19) ), #Cell 10
    list( c(8,9,10,11,12), c(2,6,11,16),    c(7,11,15,18)   ), #Cell 11
    list( c(8,9,10,11,12), c(2,6,11,16),    c(7,11,15,18)   ), #Cell 12
    list( c(2,5,9,13),     c(8,13,17),      c(13,14,15,16)  ), #Cell 13
    list( c(3,6,10,14,17), c(4,9,14,18),    c(13,14,15,16)  ), #Cell 14
    list( c(1,5,10,15,19), c(7,11,15,18),   c(13,14,15,16)  ), #Cell 15
    list( c(2,6,11,16),    c(13,14,15,16),  c(12,16,19)     ), #Cell 16
    list( c(8,13,17),      c(3,6,10,14,17), c(17,18,19)     ), #Cell 17
    list( c(4,9,14,18),    c(7,11,15,18),   c(17,18,19)     ), #Cell 18
    list( c(12,16,19),     c(17,18,19),     c(1,5,10,15,19) )  #Cell 19
  )
  
  #
  # sequence class
  # Implements a double-ended head/tail accessible stack.
  #
  sequence <- function()
  { 
    sequence <- new.env()
    sequence$.vector <- vector()
    sequence$size <- function() return( length(.vector) )
    sequence$pushHead <- function(value) .vector <<- c(.vector, value)
    sequence$pushTail <- function(value) .vector <<- c(value, .vector)
    sequence$popHead <- function() {
      value <- .vector[length(.vector)]
      .vector <<- .vector[-length(.vector)]
      return(value)
    }
    sequence$popTail <- function() {
      value <- .vector[1]
      .vector <<- .vector[-1]
      return(value)
    }
    environment(sequence$size)     <- as.environment(sequence)
    environment(sequence$pushHead) <- as.environment(sequence)
    environment(sequence$popHead)  <- as.environment(sequence)
    environment(sequence$pushTail) <- as.environment(sequence)
    environment(sequence$popTail)  <- as.environment(sequence)
    class(sequence) <- "sequence"
    return(sequence)
  }
  
  evaluateBoard <- function(board,target)
  {
    for (line in all_lines)
    {
      total <- 0
      for (cell in line) {
        total <- total + board[cell]
      }
      if (total != target) return(FALSE)
    }
    if (doPrint) {
      print(paste0("SUCCESS === ",tilesPlaced," at ",format(Sys.time(), "%H:%M:%OS3")))
      print(board)
    }
    return(TRUE)
  }

  evaluateCell <- function(board,target,cellplaced)
  {
    for (lines in cell_lines[cellplaced])
    {
      for (line in lines)
      {
        total <- 0
        checkExact <- TRUE
        for (cell in line)
        { # cells without tiles yet placed in a line mean target test can't apply
          if (board[cell] == 0) checkExact <- FALSE 
          else total <- total + board[cell]
        }
        if ( (checkExact && (total != target)) || total > target) return(FALSE)
      }
    }
    return(TRUE)
  }
  
  #
  # placeTiles
  #
  # placeTiles executes a depth first search via recursion, trying to place tiles
  # at the next empty cell until, either we have placed all tiles on the board with
  # all lines adding up to 38, or we can't find any remaining tile to place correctly.
  #  
  placeTiles <- function(cells,board,tilesRemaining)
  {
    for (cell in cells)
    {
      if (board[cell] != 0) next # Skip this cell if a tile already placed here
      
      maxTries <- tilesRemaining$size()
      for (t in 1:maxTries)
      {
        board[cell] = tilesRemaining$popHead()
        
        if (doPrint && doProgress && (tilesPlaced %% progressRate == 0)) {
          print(paste0("=== ",tilesPlaced," at ",format(Sys.time(), "%H:%M:%OS3")))
          print(board)
          print(tilesRemaining$.Data)
        }
        tilesPlaced <<- tilesPlaced + 1
        
        cellok <- evaluateCell(board,38,cell)
        if (cellok) {
          retval <- placeTiles(cells,board,tilesRemaining)
          if (retval$Success) return(retval)
        }
        tilesRemaining$pushTail(board[cell])
        board[cell] = 0
      }
      
      # We have tried all the available tiles for this cell without success
      return( list(Success = FALSE, Board = board, TilesPlaced = tilesPlaced) )
    }
    
    success <- evaluateBoard(board,38)
    return( list(Success = success, Board = board, TilesPlaced = tilesPlaced) )
  }
  
  teval <- function(...) {
    gc();
    start <- Sys.time()
    pstart <- proc.time()
    result <- eval(...)
    pfinish <- proc.time()
    finish <- Sys.time()
    wallclock = list(Start=format(start,"%H:%M:%OS3"),Finish=format(finish,"%H:%M:%OS3"))
    return ( list(PID=Sys.getpid(), WallClock=wallclock, Duration=pfinish-pstart, Result=result) )
  }
  
  ##
  ## Main body of solver
  ##
  tile1 <- as.integer(substr(triple,1,2))
  tile2 <- as.integer(substr(triple,3,4))
  tile3 <- as.integer(substr(triple,5,6))
  board <- c(tile1,tile2,tile3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
  cells <- c(4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19)
  tiles <- sequence()
  for (t in 1:19) {
    if (t == tile1 || t == tile2 || t == tile3) next
    tiles$pushHead(t)
  }
  result <- teval(placeTiles(cells,board,tiles))
  return ( list(Triple = triple, Result = result$Result, PID=result$PID,
                WallClock = result$WallClock, Duration = result$Duration) )
}

##
# run_serial() run the solver on a single core (classic R but slow)
run_serial <- function()
{
  tri <- generateTriples()
  #print(format(Sys.time(), "%H:%M:%OS3 - calling lapply(solver)"))
  results <- teval( lapply(tri,solver) )
  #print(format(Sys.time(), "%H:%M:%OS3 - lapply done"))
  results
}

##
# run_mc_fork(detectCores()) run the solver on as many cores as your system has (only works on OSX/Linux)
run_mc_fork <- function(cores)
{
  tri <- generateTriples()
  #print(format(Sys.time(), "%H:%M:%OS3 - calling mclapply(solver)"))
  results <- teval( mclapply(tri, solver, mc.cores=cores) )
  #print(format(Sys.time(), "%H:%M:%OS3 - mclapply done"))
  results
}

##
# run_aws_emr() runs the solver on AWS Elastic Map Reduce with your AWS account credentials (and costs money)
global_results <- list(1)
run_aws_emr <- function() # This function runs in the cloud and costs money!
{
  # I have had to significantly upgrade the original segue package
  # to work with the latest version of AWS API and EMR configuration.
  # Please use the segueSC version included in this book code directory.
  library(segueSC)
  
  ##firstcells <- as.list(1:19)
  triples <- generateTriples()
  tasks <- list()
  for (t in triples) {
    tasks[length(tasks)+1] = t
  }
  
  ## To get the answer as quickly as possible,
  ## we will use as many Hadoop nodes as there are board starting positions
  ## plus one more for the Hadoop master node
  numNodes <- 16
 
  ## Insert your account's security credentials here: AWS Key and Secret Key
  ## You can obtain these from Your Account/Security Credentials on the AWS console
  ## Never share this information with any other person
  setCredentials("ABCDEFGHIJKLMNOPQRST","abcdefghijklmnopqrstuvwxyz0123456789/ABC")
  print(format(Sys.time(), "%H:%M:%OS3 - creating cluster"))
  
  ## When we create the cluster it will by default be created in AWS region us-east
  ## using the availability zone 1c. If like me you live in the UK then you may
  ## want to change it to use AWS Ireland, i.e. eu-west and pick either 1a,1b or 1c
  ## You can find out which availabilty zones are up and running by checking the tab
  ## Service Health on the main EC2 service AWS console. At the time of writing all
  ## Amazon regions supported the EMR service.
  myCluster <- createCluster(numInstances=numNodes) #, enableDebugging=TRUE) #,location="eu-west-1a")
  print(format(Sys.time(), "%H:%M:%OS3 - cluster done"))
  
  ## Parallel calculation (mclapply):
  print(format(Sys.time(), "%H:%M:%OS3 - calling emrlapply"))
  results <- emrlapply(myCluster, tasks, solver, taskTimeout=10)
  print(format(Sys.time(), "%H:%M:%OS3 - emrlapply done"))
  print(results)
  global_results <<- results
  
  ## We must remember to save our bank balance!
  print(format(Sys.time(), "%H:%M:%OS3 - terminating cluster"))
  stopCluster(myCluster)
  print(format(Sys.time(), "%H:%M:%OS3 - cluster terminated"))

  results
}


## Instrumenting code
profileFn <- function(fn)       ## Turn on tracing for “fn”
{
  assign("profile.counter",0,envir=globalenv())
  trace(fn,quote(assign("profile.counter",get("profile.counter",envir=globalenv()) + 1,
                 envir=globalenv())), print=FALSE)
}
profileFnStats <- function(fn)  ## Get collected stats
{
  count <- get("profile.counter",envir=globalenv())
  return( list(Function=fn,Count=count) )
}
unprofileFn <- function(fn)     ## Turn off tracing and tidy up
{
  remove(list="profile.counter",envir=globalenv())
  untrace(fn)
}


## Graph AWS runtime results
plotDurations <- function(results,title="AWS EMR Solver Execution Profile")
{
    heights <- vector()
    names <- vector()
    colors <- vector()
    elapsedMax <- 0.0
    elapsedMin <- Inf
    elapsedSum <- 0.0
    fastest <- Inf
    solution <- ""
    eusdiffMax <- 0
    for (res in results) {
        names[length(names)+1] <- res$Triple
        elapsed <- res$Duration[3]
        elapsedSum <- elapsedSum + elapsed
        usrsys <- res$Duration[1] + res$Duration[2]
        eusdiff <- (elapsed - usrsys) / elapsed
        if (eusdiff > eusdiffMax) eusdiffMax = eusdiff
        if (elapsed < elapsedMin) elapsedMin <- elapsed
        if (elapsed > elapsedMax) elapsedMax <- elapsed
        heights[length(heights)+1] <- elapsed
        color <- "white"
        if (res$Result$Success) {
            color <- "black"
            if (elapsed < fastest) {
                fastest=elapsed
                solution=res$Triple
            }
        }
        colors[length(colors)+1] <- color
    }
    
    elapsedMean <- elapsedSum / (length(heights)*1.0)
    minMins <- elapsedMin / 60.0
    minSecs <- elapsedMin %% 60
    maxMins <- elapsedMax / 60.0
    maxSecs <- elapsedMax %% 60
    avgMins <- elapsedMean / 60.0
    avgSecs <- elapsedMean %% 60
    
    subtitle <- paste0("Boards=", length(heights),
    " -- Min=", sprintf("%.0f",minMins), "m", sprintf("%.0f",minSecs), "s",
    " Max=", sprintf("%.0f",maxMins), "m", sprintf("%.0f",maxSecs), "s",
    " Avg=", sprintf("%.0f",avgMins), "m", sprintf("%.0f",avgSecs), "s",
    " -- Fastest solution ",solution," in ",sprintf("%.0f",fastest),"s")
    barplot(heights,names.arg=names,col=colors,horiz=TRUE,space=2.0,
    main=title,
    sub=subtitle,
    xlab="Elapsed Time (s)",ylab="First Three Tiles", xlim=c(0,360),
    cex.names=0.7,cex.axis=0.75)
    return( list(Title=title, Stats=subtitle,
    Min=elapsedMin, Max=elapsedMax, Avg=elapsedMean,
    Total=elapsedSum, EUSMaxDiffPC=eusdiffMax * 100.0) )
}

##
## Packt: "Mastering Parallelism with R"
## Chapter 1 - Aristotle's Number Puzzle
##
## Copyright 2015 Simon Chapple
##


