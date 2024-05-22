#' Wyatt Petryshen
#' Yale University
#' January 27th, 2024
#'
#' Main function for the Metacommunity Model in Ecology from Leibhold and Chase 2018
#'
#' Objective is to model all four major processes: 1) Patch Dynamics; 2) Species Sorting; 3) Mass-Effect; 4) Neutral
#'
#' @param initial.metacommunity A list of vectors containing local populations that comprise the metacommunity
#' @param probability.immigration Probability of species immigration or speciation
#' @param death.rate A positive integer value for the death rate per time step
#' @param sim.generations Number of generations to simulate
#' @param metacommunity.distribution A boolean argument to allow for dynamic Jm
#' @param type Method for new species generation; defualt "Point"
#' @param lcommunity.env.response A vector containing the probability distribution of the  metacommunity
#' @param print.dialog Default False. Prints when simulation is complete and elpased time
#' @param save.file.path Default Null. Saves simulation to disk as an RDS file
#' @param cluster False. If running on cluster in terminal supress the return of results
#' @returns A numeric matrix of species identity
#' @example
#' untb(c(1,1,2,2), 0.1, 1, 100)
#' As a note regarding performance; all for loops and if statements are extremely slow
#' Changing to purrr:walk will improve for loop performance
#'
#' This model is a derivative of the NTBB model
# dependent functions from Rcpp
library(Rcpp)
sourceCpp("/Users/wyattpetryshen/Documents/GitHub/UNTB.Haar/Main/RCPP/nextGenSampler.cpp")
sourceCpp("/Users/wyattpetryshen/Documents/GitHub/UNTB.Haar/Main/RCPP/rcppSample.cpp")

setClass("NTBB Class",
         slots = c(date = "POSIXct", initial.abundance = "numeric",
                   probability.immigration = "numeric",
                   death.rate = "numeric",
                   sim.generations = "numeric",
                   sampling.scheme = "character",
                   lcommunity.env.response = "character",
                   results = "array"))

setClass("NTBB Save Class",
          slots = c(date = "POSIXct", current.generation = "numeric",
                    step.abundance = "list",
                    step.immigration = "vector",
                    step.death.rate = "vector",
                    sampling.scheme = "character",
                    meta.data = "character"))
# I will want to add some additional Rcpp functions to speed this up
fDn_no <- function(m, J){
  return(2*J*m)
}

fDn_ol <- function(m, J){
  return((m*(J-1))/(1-m))
}

assingment_prob <- function(x, Pr.r, Pr.s){
  return(rcppSample(x, Pr.r, Pr.s))
}

freq_probabilites <- function(x, uQ){
  # x is occurrence data at a timestep
  # uQ is the set of x
  uQ.Size <- length(uQ)
  prob.X <- numeric(length = uQ.Size)
  x.Size <- length(x)
  for(i in 1:uQ.Size){
    prob.X[i] = length(x[x == uQ[i]])/x.Size
  }
  return(prob.X)
}

RnextgenSampler <- function(x, density.Dependence = F){
  if(density.Dependence == T){
    uQ <- unique(x)
    prob.X <- freq_probabilites(x, uQ)
    rS = sample(uQ,1,prob = 1/prob.X)
  }
  else{
    rS = sample(x,1)
  }
  return(rS)
}

assignment_function <- function(x, nex.G, max.G, metacommunity.distribution, type){
  # draw replacement probabilities
  if(x == 1) {
    return(RnextgenSampler(nex.G))
  } else {
    if(type == "Point"){return(x)}
    else if(type == "Allan&Savage"){return(max.G + 1)}
    else if(type == "Jm"){return(sample(metacommunity.distribution, 1))}
  }
}

save_step_function <- function(save_path, MC_save_class){
  saveRDS(MC_save_class, save_path)
}

MCME <- function(initial.metacommunity.list,
  community.distances.matrix,
  probability.immigration.vector,
  death.rate.vector,
  sim.generations,
  lcommunity.env.response,
  type,
  metacommunity.distribution = NULL,
  print.dialog = F,
  save.file.path = NULL,
  cluster = F){
  # begin function
  format(Sys.time(), "%S")
  start_time <- Sys.time()

  # input paramter checks
  # number of local communities within the metacommunity
  # defined as a user input by initial.metacommunity.list
  num_Jl <- length(initial.metacommunity.list)

  ## check inital abundance for each local community within metacommunity
  if(num_Jl <= 0){
    return(print("Intial metacommunity must be provided"))
  }

  ## check that the distance matrix between local communities matchs the number of communities
  if(length(community.distances.matrix) != num_Jl){
    return(print("Distance matrix between local communities must eqaul number of intial local communities within the metacommunity"))
  }

  ## check lcommunity.env.response
  if(!is.null(lcommunity.env.response) & length(lcommunity.env.response) != sim.generations){
    stop(print("lcommunity.env.response must be equal length to simulation generations"))
  }

  # initiate time 0 local communities
  # nrow = habitat patchs; ncol = time steps
  d.Vector <- death.rate.vector # death rate; fixed valued
  J.Vector <- lengths(initial.metacommunity.list) # initial number of J individuals within the community
  # Current set of for initial immigration probabilities will only take in constant value
  # Could have the immigration probability evolve through time
  m.Vector <- probability.immigration.vector # resistance factor; probability of species immigration

  # Within the dJm framework we are dynamically defining Jm and saving into list objects
  # Pythonic method; may be much slower need to profile
  # Each row is a new generation; need to pad and transpose the file results; can also add all zeros and ignore for plotting

  # define list to store results
  master.List = list()
  # save inital population in master list
  master.List[[1]] = unname(initial.abundance.vector)

  # Each local community gets its own array to store trait values and species occurence
  Jl.names <- paste0("Jl.", num_Jl)


  # save initial configuration as a RDS
  initalsave <- new("NTBB Save Class",
  date = Sys.time(),
  step.abundance = initial.metacommunity.list,
  step.immigration = m.Vector,
  step.death.rate = d.Vector,
  sim.generations = sim.generations,
  sampling.scheme = sample.parms,
  meta.data = NULL)

  saveRDS(initalsave, save_path)

  if(print.dialog == T){print("Starting simulation")}
  # Simulation #



}
