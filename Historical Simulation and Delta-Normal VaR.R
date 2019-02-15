# Historical Simulation and Delta-Normal VaR
# by Changjie Ma
# Feb. 2019

library(quantmod)
library(PortfolioAnalytics)
library(xts)
library(ggplot2)

# Global Paramaters
weights.portfolio <- c(1/4,1/4,1/4,1/4) # equal weights for the 4 indexes
mu <- c(0,0,0,0) # expect returns are zero
q <- 0.01 # calculate 1%-VaR

# get data from Yahoo
adj.close <- 6  # 6th field is adjusted close
equity.tickers <- c("^GSPC","^DJI","^IXIC","^RUT")
prices <- getSymbols(equity.tickers[1],from = "2004-01-01", to = "2018-12-31", source="yahoo", auto.assign=FALSE,
                     return.class="xts")[,adj.close]
for (i.tmp in 2:length(equity.tickers)) {
  prices.tmp <- getSymbols(equity.tickers[i.tmp],from = "2004-01-01", to = "2018-12-31", source="yahoo",
                           auto.assign=FALSE, return.class="xts")[,adj.close]
  prices <- cbind(prices, prices.tmp)
}
equity.names <- c("S.P.500","Dow.30","Nasdaq","Russell.2000")
colnames(prices) <- equity.names
returns <- diff(prices)/prices # calculate simple returns
prices <- fortify(prices)

# Part 1 Historical simulation method.
# Get positions of each day
one_m = 1000000
four_m = 4 * one_m
prices$'p_sp500' <- floor(one_m / prices$`S.P.500`)
prices$'p_dow30' <- floor(one_m / prices$`Dow.30`)
prices$'p_nas' <- floor(one_m / prices$Nasdaq)
prices$'p_russ' <- floor(one_m / prices$`Russell.2000`)

# Calculate the returns of each day and VaR-at-1%
for(row in 2:nrow(prices)){
  last_price <- prices[row-1, 'S.P.500'] * prices[row-1, 'p_sp500'] + prices[row-1, 'Dow.30'] * prices[row-1, 'p_dow30'] + prices[row-1, 'Nasdaq'] * prices[row-1, 'p_nas'] + prices[row-1, 'Russell.2000'] * prices[row-1, 'p_russ']
  current_price <- prices[row, 'S.P.500'] * prices[row-1, 'p_sp500'] + prices[row, 'Dow.30'] * prices[row-1, 'p_dow30'] + prices[row, 'Nasdaq'] * prices[row-1, 'p_nas'] + prices[row, 'Russell.2000'] * prices[row-1, 'p_russ']
  prices[row, 'return'] <- (current_price - last_price) / last_price
}
for(row in 2:nrow(prices)){
  var1 <- quantile(na.omit(prices[max(1,row-999):row, 'return']), 0.01)
  prices[row, 'var1'] <- var1
}
sub_prices <- prices[prices$Index > '2007-12-31' & prices$Index <= '2009-12-31',]
#ggplot(data=sub_prices, aes(x=Index, y=-var1, group=1)) +
#  geom_line()


# Part 2 Delta-Normal method with the equally-weighted covariance matrix estimator.
VaR.alldays.q2 <- c()
start <- length(index(returns['/2007-12-31'])) # the number of row for the start day in the 2-year horizon
for (i.tmp in index(returns['2008/2009'])){
  returns.data <- returns[(start-1000):(start-1)]
  variance.tmp <- c()
  for (j in 1:4){
    for (k in 1:4){
      variance.tmp <- c(variance.tmp,mean(returns.data[,j]*returns.data[,k]))
    }
  }
  stdev.tmp <- sqrt(mean(variance.tmp))
  VaR.tmp <- qnorm(q)*stdev.tmp
  start <- start+1
  VaR.alldays.q2<-c(VaR.alldays.q2,VaR.tmp[[1]])
}
VaR.dataframe.q2 <- xts(VaR.alldays.q2,order.by = index(returns['2008/2009']))
colnames(VaR.dataframe.q2)<-'1%-VaR'

# Part 3 Delta-Normal method using the exponentially-weighted covariance estimator.
lambda <- 0.94
n <- 100 # 0.94^100=0.002054875 << 0 so we choose 100 days of data
weights.data <- c() # Generate the weight vector for 100 days
for (tau.tmp in 1:n){
  weight <- lambda^(tau.tmp-1) # weights assigned to days
  weights.data <- c(weight,weights.data)
}
VaR.alldays.q3 <- c()
start <- length(index(returns['/2007-12-31'])) # the number of row for the start day in the 2-year horizon
for (i.tmp in index(returns['2008/2009'])){
  returns.data <- returns[(start-n):(start-1)]
  weights.dataframe <- xts(data.frame(weights.data),order.by = index(returns.data))
  variance.tmp <- c()
  for (j in 1:4){
    for (k in 1:4){
      variance.tmp <- c(variance.tmp,sum(returns.data[,j]*returns.data[,k]*weights.dataframe)*(1-lambda))
    }
  }
  stdev.tmp <- sqrt(mean(variance.tmp))
  VaR.tmp <- qnorm(q)*stdev.tmp
  VaR.alldays.q3<-c(VaR.alldays.q3,VaR.tmp[[1]])
  start <- start+1
}
VaR.dataframe.q3 <- xts(VaR.alldays.q3,order.by = index(returns['2008/2009']))
colnames(VaR.dataframe.q3)<-'1%-VaR'

# Part 4 Weighted historical simulation.
for (row in 1:nrow(prices)){
  if(prices[row, 'Index'] > '2009-12-31'){
    break
  }
  if(prices[row, 'Index'] <= '2007-12-31'){
    next
  }
  sub_df <- prices[(row-999):row, ]
  for (row2 in 1:nrow(sub_df)){
    sub_df[row2, 'weight'] <- 0.995^(1000 - row2)*(1-0.995) / (1-0.995^(1000))
  }
  ordered_sub_df <- sub_df[order(sub_df$return),]
  threshold <- 0.01
  sum <- 0
  for (row3 in 1:nrow(ordered_sub_df)){
    sum <- sum + ordered_sub_df[row3, 'weight']
    if(sum > threshold){
      prices[row, 'w_var1'] <- ordered_sub_df[row3, 'return']
      break
    }
  }
}

sub_prices <- prices[prices$Index > '2007-12-31' & prices$Index <= '2009-12-31',]
#ggplot(data=sub_prices, aes(x=Index, y=-w_var1, group=2)) +
#  geom_line()



sub_prices[,'var2'] <- VaR.alldays.q2
sub_prices[,'var3'] <- VaR.alldays.q3

# Plotting all graphs
ggplot(sub_prices, aes(Index)) + 
  geom_line(aes(y = var1, colour = "VaR1")) + 
  geom_line(aes(y = var2, colour = "VaR2")) + 
  geom_line(aes(y = var3, colour = "VaR3")) + 
  geom_line(aes(y = w_var1, colour = "VaR4"))


# Part 5 Delta-Normal VaR of another portfolio

# (a)

S0 <- 2647.08
K <- 2650
T <- 0.25
t <-0
r <- 0.01
d <- 0.02
sigma <- 0.18

d1 <- (log(S0/K)+(r-d+sigma^2/2)*(T-t))/(sigma*sqrt(T-t))
d2 <- d1-sigma*sqrt(T-t)

C <- S0*exp(-d*(T-t))*pnorm(d1)-K*exp(-r*(T-t))*pnorm(d2)
P <- K*exp(-r*(T-t))*pnorm(-d2)-S0*exp(-d*(T-t))*pnorm(-d1)
value.portfolio <- -50*100*C-50*100*P

delta.call <- exp(-d*(T-t))*pnorm(d1)
delta.put <- -exp(-d*(T-t))*pnorm(-d1)
gamma.call <- exp(-d*(T-t))*dnorm(d1)/(S0*sigma*sqrt(T-t))
gamma.put <- exp(-d*(T-t))*dnorm(d1)/(S0*sigma*sqrt(T-t))
vega.call <- exp(-d*(T-t))*S0*dnorm(d1)*sqrt(T-t)
vega.put <- exp(-d*(T-t))*S0*dnorm(d1)*sqrt(T-t)
theta.call <- -exp(-d*(T-t))*S0*sigma*dnorm(d1)/(2*sqrt(T-t))-r*K*exp(-r*(T-t))*pnorm(d2)+d*S0*exp(-d*(T-t))*pnorm(d1)
theta.put <- -exp(-d*(T-t))*S0*sigma*dnorm(d1)/(2*sqrt(T-t))+r*K*exp(-r*(T-t))*pnorm(-d2)-d*S0*exp(-d*(T-t))*pnorm(-d1)

delta <- -50*100*delta.call-50*100*delta.put
gamma <- -50*100*gamma.call-50*100*gamma.put
vega <- -50*vega.call-50*vega.put
theta <- -50*theta.call-50*theta.put
theta.daily <- theta/250
print(value.portfolio)
print(delta)
print(gamma)
print(vega)
print(theta)
print(theta.daily)

# (b)
stdev.s0 <- 0.01*S0
VaR.1.portfolio <- qnorm(0.05)*abs(delta)*stdev.s0
print(VaR.1.portfolio)

# (c)
stdev.vol <- 0.03*sigma
corr <- -0.7
stdev.portfolio <- sqrt((delta*stdev.s0)^2+(vega*100*stdev.vol)^2+2*corr*delta*stdev.s0*vega*100*stdev.vol)
VaR.2.portfolio <- qnorm(0.05)*stdev.portfolio
print(VaR.2.portfolio)




