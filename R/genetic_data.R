genetic_data <- function(ptree, stop_time) {
  t <- seq(0, stop_time, by=1)
  l <- stop_time + 1 #length of epidemic

  #day is 0:n_days
  data <- list("day" = t, "n_anc" = 0, "n_coal" = NA)

  ltt <- ape::ltt.plot.coords(ptree) #gives number of lineages through time
  ltt <- list("time"=ltt[,1],"N"=ltt[,2]) #make a list of time and number of lineages
  ltt$time <- stop_time + ltt$time #forwards in time
  if (ltt$time[1]==0) {
    ltt$time <- ltt$time[-1]
    ltt$N <- ltt$N[-1]
  }

  #number of lineages through time
  for (i in 1:l) {
    #the first occasion where ltt$time is larger than t[i]
    j <- which.max(ltt$time >= t[i])
    #at that time, the number of lineages is ltt$N
    data$n_anc[i] <- ltt$N[j]
  }

  #coalescence times
  coal_times <- c()
  j <- 0
  for (i in 1:(length(ltt$time)-1)) {
    if (ltt$N[i] < ltt$N[i+1]) {
      j <- j + 1
      coal_times[j] <- ltt$time[i]
    }
  }

  #number of coalescences per day
  for (i in 2:l) {
    data$n_coal[i] <- sum(coal_times>=t[i-1] & coal_times<t[i])
  }

  data <- as.data.frame(data)
  return(data)
}
