#' @include geomod3D.R
NULL

# This file contains the functions necessary to train a GP with a genetic
# algorithm. In the future these functions might be turned into a class.

# Real-valued chromossomes always have values between 0 and 1. They are
# scaled back to input range at the time of evaluation.

# Some functions have a temperature parameter, which is also contained in the
# (0, 1) interval and decreases with the iteration number. It controls the
# amount of change in the chromossomes.

#### real encoding - mutation ####

GeneticRealMutationRadial <- function(parent, temperature){
  # Treats the chromossome as coordinate in Euclidean space, shifting it
  # a small distance in a random direction.

  minval <- parent - temperature
  minval[minval < 0] <- 0
  maxval <- parent + temperature
  maxval[maxval > 1] <- 1

  runif(n = length(parent),
        min = minval,
        max = maxval)

  # parent <- parent + runif(n = length(parent),
  #                          min = minval,
  #                          max = maxval)

  # parent[parent < 0] <- 0
  # parent[parent > 1] <- 1

  # parent
}

GeneticRealMutationParallel <- function(parent, temperature){
  # Treats the chromossome as individual values. The temperature parameter
  # controls the probability of change for each value. If a value is changed,
  # it can take any value in its range.

  # new values
  child <- runif(length(parent))

  # controls which values will be changed
  change <- sample(c(0, 1), size = length(parent), replace = T,
                   prob = c(1 - temperature, temperature))

  # output
  parent * (1 - change) + child * change
}


#### real encoding - crossover ####

GeneticRealCrossover1p <- function(parent1, parent2){
  xp <- sample(length(parent1) - 1, 1)

  child1 <- parent1; child2 <- parent2

  child1[1:xp] <- parent2[1:xp]
  child2[1:xp] <- parent1[1:xp]

  matrix(c(child1, child2), 2, length(parent1), byrow = T)
}

GeneticRealCrossover2p <- function(parent1, parent2){
  if (length(parent1) < 3)
    return(GeneticRealCrossover1p(parent1, parent2))

  xp1 <- sample(length(parent1) - 2, 1)
  xp2 <- sample((xp1 + 1):(length(parent1) - 1), 1)

  child1 <- parent1; child2 <- parent2

  child1[xp1:xp2] <- parent2[xp1:xp2]
  child2[xp1:xp2] <- parent1[xp1:xp2]

  matrix(c(child1, child2), 2, length(parent1), byrow = T)
}

GeneticRealCrossoverAverage <- function(parent1, parent2){
  frac <- runif(2)

  child1 <- frac[1] * parent1 + (1 - frac[1]) * parent2
  child2 <- frac[2] * parent2 + (1 - frac[2]) * parent1

  matrix(c(child1, child2), 2, length(parent1), byrow = T)
}

#### selection and fitness sharing ####
FindMinCount <- function(D){
  # For a given distance matrix, finds the minimum distance that gives at
  # least a count of 1 neighbor for all points.

  if (sum(D) == 0) return(0)

  diag(D) <- NA

  dmin <- min(D, na.rm = T)
  dmax <- max(D, na.rm = T)

  for (i in 0:1000) {
    x <- (dmax - dmin) * i / 1000 + dmin
    count <- rowSums(D <= x, na.rm = T)
    if (min(count) == 1) break
  }

  return(x)
}


ShareFitness <- function(popmatrix, fitvec, distance){
  # Shares fitness based on the number of neighbors within the radius
  # that gives at least 1 neighbor for each individual.
  #
  # Returns a vector with positive values representing relative probability
  # for selection.

  D <- as.matrix(dist(popmatrix, distance))

  x <- FindMinCount(D) # may be zero
  dmin <- quantile(D, probs = 0.1) + 1e-6

  sh <- matrix(1, nrow(D), ncol(D))
  sh <- sh - D / (x + dmin) # dmin avoids division by zero
  sh[sh < 0] <- 0

  w <- rowSums(sh) # weights according to number of neighbors and distance

  amp <- max(fitvec) - min(fitvec)
  if (amp == 0) amp <- 10
  fitvec <- fitvec - min(fitvec) + 0.1 * amp

  fitvec / w
}

#### the genetic algorithm ####
GeneticTrainingReal <- function(fitness, minval, maxval, popsize = 50,
                                initpop = ceiling(popsize / 5),
                                mutation = list(GeneticRealMutationRadial,
                                                GeneticRealMutationRadial),
                                crossover = list(GeneticRealCrossover1p,
                                                 GeneticRealCrossover2p,
                                                 GeneticRealCrossoverAverage),
                                mutprob = 0.35,
                                distance = "euclidean",
                                maxiter = 1000,
                                tol = 0.1, stopping = ceiling(maxiter / 5),
                                start = NULL,
                                blocks = rep(1, length(minval)), cycle = 25,
                                verbose = F){
  # This genetic algorithm works by selection two individuals from the
  # population, performing crossover and mutation. The new individuals replace
  # the ones with the smallest fitness.
  #
  # To save computation time, the population can be initialized with less
  # individuals than the maximum value, and filled as the main loop progress.
  #
  # Individuals are selected using fitness sharing to preserve diversity.

  # Checking
  if (length(minval) != length(maxval))
    stop("Maximum and minimum vectors have different lengths")
  if (!is.null(start) && length(start) != length(maxval))
    stop("Wrong length for starting point")
  if (length(blocks) != length(maxval))
    stop("Length of blocks must match chromossome length")

  # Initialization
  L <- length(maxval) # chromossome length

  popmatrix <- matrix(NA, popsize, L)
  popmatrix[1:initpop, ] <- runif(initpop * L)
  if (!is.null(start)){
    start <- (start - minval) / (maxval - minval + 1e-9)
    popmatrix[1, ] <- start
  }

  fitvec <- rep(NA, popsize)
  if (verbose) cat("Initializing population ")
  for (i in seq(initpop)) {
    if (verbose) cat(".")
    fitvec[i] <- fitness(popmatrix[i, ] * (maxval - minval) + minval)
  }
  if (verbose) cat("\n\n")
  fitmatrix <- matrix(NA, popsize, maxiter + 1)
  fitmatrix[, 1] <- fitvec

  bestfitness <- fitvec[which.max(fitvec)]
  bestsol <- popmatrix[which.max(fitvec), ]

  # Main loop
  stagnation <- 0
  temperature <- 1
  block_count <- rep(sort(unique(blocks)), each = cycle, length.out = maxiter)
  if (verbose) cat("Iteration 0 - Best fitness =", max(fitvec, na.rm = T))
  for (i in seq(maxiter)) {
    # selecting block
    current_block <- blocks == block_count[i]

    # selection
    selectable <- !is.na(fitvec)
    id <- sample(sum(selectable), size = 2,
                 prob = ShareFitness(popmatrix[selectable, current_block],
                                     fitvec[selectable],
                                     distance))
    tmp1 <- popmatrix[selectable, ][id[1], ]
    tmp2 <- popmatrix[selectable, ][id[2], ]

    # mutation/crossover
    crossfun <- crossover[[sample(length(crossover), 1)]]
    child <- crossfun(tmp1[current_block], tmp2[current_block])

    mutfun <- mutation[[sample(length(mutation), 1)]]
    if (runif(1) < mutprob) child[1, ] <- mutfun(child[1, ], temperature)
    if (runif(1) < mutprob) child[2, ] <- mutfun(child[2, ], temperature)

    # update population
    # replace minimum fitness
    id_upd <- order(fitvec, decreasing = F, na.last = F)

    # replace in random position (except the best solution)
    # pos <- which.max(fitvec)
    # w <- rep(1, popsize); w[pos] <- 0
    # id_upd <- sample(popsize, 2, prob = w)


    tmp1[current_block] <- child[1, ]
    tmp2[current_block] <- child[2, ]

    popmatrix[id_upd[1], ] <- tmp1
    popmatrix[id_upd[2], ] <- tmp2
    fitvec[id_upd[1]] <- fitness(tmp1 * (maxval - minval) + minval)
    fitvec[id_upd[2]] <- fitness(tmp2 * (maxval - minval) + minval)

    # log
    fitmatrix[, i + 1] <- fitvec
    if (verbose) cat("\rIteration", i, "- Best fitness =",
                     bestfitness)

    # stopping criterion
    ev <- max(fitvec, na.rm = T) - bestfitness #max(fitmatrix[, i], na.rm = T)
    if (ev < tol){
      stagnation <- stagnation + 1
      # temperature <- min(1, temperature + 0.01)
      temperature <- max(0.05, min(temperature + 0.01, 1 - i / maxiter))
    }
    else{
      stagnation <- 0
      # temperature <- max(0.05, 1 - i / maxiter)
      temperature <- 0.05
      bestfitness <- max(fitvec, na.rm = T)
      bestsol <- popmatrix[which.max(fitvec), ]
    }
    if (stagnation >= stopping){
      if (verbose) cat("\nTerminating training at iteration", i)
      break
    }
  }
  if (verbose) cat("\n")

  # output
  for (i in seq(popsize)) {
    popmatrix[i, ] <- popmatrix[i, ] * (maxval - minval) + minval
  }
  bestsol <- bestsol * (maxval - minval) + minval
  # bestfitness <- fitvec[which.max(fitvec)]
  # bestsol <- popmatrix[which.max(fitvec), ]

  list(bestsol = bestsol,
       bestfitness = bestfitness,
       lastpop = popmatrix,
       evolution = apply(fitmatrix, 2, function(y){
         if (all(is.na(y))) return(NA) else return(max(y, na.rm = T))
       }))
}
