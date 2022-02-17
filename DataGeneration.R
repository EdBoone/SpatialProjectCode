library(geoR)

# Data Generation

# The spackage used is geoR, and its specific function grf
# ?grf
# First try using geoR examples

# Notation:
#- phi is the range parameter
#- sigma^2 is the partial sill
#- tau^2 is the nugget effect

#- default cov.model = "matern"
#- default method = "cholesky"
#- default kappa = 0.5


### ---------- Simulation 0
# Simulation of bivariate data where 
# the two variables have different correlation functions (R)
# that is:
# Z1 =  mu1 + sigmasq1 * R1 + e1
# Z2 =  mu2 + sigmasq2 * R2 + e2
library(geoR)

set.seed(357)
mu1 <- mu2 <- 0
sigmasq1 <- 1; sigmasq2 <- 3
phi1 <- 0.1; phi2 <- 0.3
tausq1 <- tausq2 <- 0 

var1 <- grf(100, grid = "irreg", xlims = c(0, 1), ylims = c(0, 1),
             cov.pars=c(sigmasq1,phi1), nug=tausq1, mean=mu1)
## simulation step by step 
#no <-100
#R1 <- grf(100, cov.pars=c(1, phi1)) 
#e1 <- rnorm(no, sd=sqrt(tausq1))
#var1 <- list(coords=gr, data = mu1 + sqrt(sigmasq1)*R1$data + e1)
#class(var1) <- "geodata"
#plot(var1$data)

var2 <- grf(100, grid = "irreg", xlims = c(0, 1), ylims = c(0, 1),
            cov.pars=c(sigmasq2,phi2), nug=tausq2, mean=mu2)

plot(var1$data)
plot(var2$data)

simdata0 = list(x = var1$coords[, 1], y = var1$coords[, 2], 
                z1 = var1$data, z2 = var2$data)

# grid for prediction
gr1 <- expand.grid((0:100)/100, (0:100)/100)
names(gr1) <- c('x','y')


################################## Simulation on the grid

### ---------- Simulation 1 (as simulation 0 but on the grid)
# Simulation of bivariate data where 
# the two variables have different correlation functions (R)
# that is:
# Z1 =  mu1 + sigmasq1 * R1 + e1
# Z2 =  mu2 + sigmasq2 * R2 + e2

set.seed(357)
gr <- expand.grid((0:30)/30, (0:30)/30)
n <-  nrow(gr) 

mu1 <- mu2 <- 0
sigmasq1 <- 1; sigmasq2 <- 3
phi1 <- 0.1; phi2 <- 0.3
tausq1 <- tausq2 <- 0 

var1 <- grf(grid=gr, cov.pars=c(sigmasq1,phi1), nug=tausq1, mean=mu1)

var2 <- grf(grid=gr, cov.pars=c(sigmasq2,phi2), nug=tausq2, mean=mu2)

par(mfrow=c(2,2), mar=c(2,2,0,0), mgp=c(1.5,.6,0))
image(var1)
image(var2)
plot(var1$data)
plot(var2$data)

cor(var1$data, var2$data) 
simdata1 = list(x = var1$coords[, 1], y = var1$coords[, 2], 
                z1 = var1$data, z2 = var2$data)
#head(simdata1)


### ---------- Simulation 2
# Simulation of bivariate data where 
# the two variables have the same correlation functions (R)
# that is:
# Z1 =  mu1 + sigmasq1 * R + e1
# Z2 =  mu2 + sigmasq2 * R + e2
set.seed(357)
mu1 <- mu2 <- 0
sigmasq1 <- 1; sigmasq2 <- 4
phi <- 0.25
# tausq1 <- tausq2 <- 0
tausq1 <- 0.2; tausq2 <- 0.5 

## simulation of a common component
var <- grf(grid=gr, cov.pars=c(1, phi))

# simulation variable 1 
var1 <- var
var1$data <- mu1 + sqrt(sigmasq1) * var$data

# or with nugget effect
var1a <- var
var1a$data <- sqrt(sigmasq1) * var$data
var1a$data <- mu1 + var1a$data + rnorm(n, mean=0, sd=sqrt(tausq1))


## simulation variable 2
var2 <- var
var2$data <- mu2 + sqrt(sigmasq2) * var$data

# or with nugget effect
var2a <- var
var2a$data <- sqrt(sigmasq2) * var$data
var2a$data <- mu2 + var2a$data + rnorm(n, 0, sqrt(tausq2))


# plot data
par(mfrow=c(2,2), mar=c(1.7,1.7,0,0), mgp=c(1.5,0.5,0))
image(var1)
image(var2)
plot(var1$data)
plot(var2$data)

image(var1a)
image(var2a)
plot(var1a$data)
plot(var2a$data)

cor(var1$data, var2$data) 
cor(var1a$data, var2a$data) 

simdata2 = list(x = var1a$coords[, 1], y = var1a$coords[, 2], 
                z1 = var1a$data, z2 = var2a$data)










