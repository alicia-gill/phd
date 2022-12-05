#' Simulate epidemics with birth rate varying over time.
#'
#' Simulates a birth-death epidemic with varying birth rate and constant death rates and 1 infected on day 0.
#'
#' @param birth_rate_fn function of the birth rate over time of the epidemic.
#' @param death_rate death rate of the epidemic.
#' @param stop_time number of days to run the epidemic simulation for.
#'
#' @return data frame with: Index, Parent Index, Infection Status, Infection Time, Removal Time
#' @export
#'
#' @examples
#' epidemic_bt(birth_rate_fn = linear, death_rate = 0.1, stop_time = 50)
epidemic_bt <- function(birth_rate_fn, death_rate, stop_time) {
  #index labels everyone by order of entry
  #parent is the number of parent
  #status is infection status: I is infected, R is removed
  #inf_time is infection time
  #rem_time is removal time
  data <- list("index"=1,"parent"=0,"status"="I","inf_time"=0,"rem_time"=NA)

  n_inf <- 1 #number infected
  n_total <- 1 #total = number infected + number recovered
  current_time <- 0 #current time, initially 0

  inf <- c(1) #vector of infected people

  bt <- birth_rate_fn(t = current_time, stop_time = stop_time)

  #run loop while some people are still infected
  while (n_inf > 0) {
    #simulate event time
    event_time <- rexp(1, rate = n_inf * (bt + death_rate))

    #if time is past stop time, then stop
    current_time <- current_time + event_time
    if (current_time > stop_time) {
      break
    }

    #sample individual to be infected or recover
    #note: if sample is of length 1, then it samples from 1:x
    #so need to do case where length is 1 separately
    if (n_inf == 1) {
      i <- inf
    } else {
      i <- sample(inf, size = 1)
    }

    bt <- birth_rate_fn(t = current_time, stop_time = stop_time)
    a <- bt / (bt + death_rate)
    u <- runif(1)
    if (u < a) {
      #infection event
      n_inf <- n_inf + 1 #number infected increased
      n_total <- n_total + 1 #total number increases
      inf[n_inf] <- n_total #add m to n.inf vector

      #add row for new infected individual
      data$index[n_total] <- n_total
      data$parent[n_total] <- i
      data$status[n_total] <- "I"
      data$inf_time[n_total] <- current_time
      data$rem_time[n_total] <- NA
    } else {
      #recovery event
      data$status[i] <- "R"
      data$rem_time[i] <- current_time
      n_inf <- n_inf - 1
      inf <- subset(inf, inf != i) #remove i from inf
    }
  }

  data <- as.data.frame(data)
  return(data)
}
