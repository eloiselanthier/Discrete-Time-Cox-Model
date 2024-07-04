# Load necessary library
library(survival)

# Define parameters
k <- 3
nbBetas <- 3
interval_size <- 10

# Load data
data1 <- read.csv("Data_site_1.csv")
data2 <- read.csv("Data_site_2.csv")
data3 <- read.csv("Data_site_3.csv")

# Combine original data for later comparison
dataOG <- rbind(data1, data2, data3)

# Adjust time and order data
adjust_and_order <- function(data) {
  data$time[data$time == 0] <- 1
  data <- data[order(data$time), ]
  return(data)
}

data1 <- adjust_and_order(data1)
data2 <- adjust_and_order(data2)
data3 <- adjust_and_order(data3)

# Get maximum time from each site
data_max <- max(c(max(data1$time), max(data2$time), max(data3$time)))

# Create intervals
# intervals <- 1:data_max
intervals <- seq(1, data_max, by = interval_size)

# Initialize parameters
alpha <- rep(0, length(intervals))
betas <- rep(0, nbBetas)
combined_vector <- c(alpha, betas)

# Function to create stacked data
create_stacked_data <- function(data, intervals) {
  stacked_data <- data.frame()
  
  for (i in 1:nrow(data)) {
    for (m in intervals) {
      if (m < data$time[i] || (data$status[i] == 1 && m == data$time[i])) {
        y_m <- ifelse(data$time[i] > m, 0, data$status[i])
        new_sample <- c(m, y_m, data$X1[i], data$X2[i], data$X3[i])
        stacked_data <- rbind(stacked_data, new_sample)
      }
    }
  }
  
  colnames(stacked_data) <- c("time", "status", "X1", "X2", "X3")
  return(stacked_data)
}

# Create stacked data for each site
stacked_data1 <- create_stacked_data(data1, intervals)
stacked_data2 <- create_stacked_data(data2, intervals)
stacked_data3 <- create_stacked_data(data3, intervals)

# Create dummy matrices for each site
create_dummy_matrix <- function(stacked_data, intervals) {
  num_rows <- nrow(stacked_data)
  row_length <- length(intervals)
  df <- data.frame(matrix(0, nrow = num_rows, ncol = row_length))
  for (i in 1:num_rows) {
    col_value <- floor(stacked_data$time[i] / interval_size) + 1
    df[i, col_value] <- 1
  }
  return(df)
}

df1 <- create_dummy_matrix(stacked_data1, intervals)
df2 <- create_dummy_matrix(stacked_data2, intervals)
df3 <- create_dummy_matrix(stacked_data3, intervals)

# Create em_xi matrices
em_xi_1 <- cbind(df1, stacked_data1[, 3:5])
em_xi_2 <- cbind(df2, stacked_data2[, 3:5])
em_xi_3 <- cbind(df3, stacked_data3[, 3:5])

# Define sigmoid and its derivative
sigmoid <- function(x) {
  1 / (1 + exp(-x))
}

dsigma <- function(x) {
  exp(-x) / ((1 + exp(-x))^2)
}

# Calculate inside values and sigmoid values
inside <- function(combined_vector, em_xi) {
  t(combined_vector %*% t(em_xi))
}

y1 <- stacked_data1[, 2]
y2 <- stacked_data2[, 2]
y3 <- stacked_data3[, 2]

# Gradient calculation
calculate_gradient <- function(y, sigmoid, dsigma, em_xi) {
  grad <- (y - sigmoid) / (log(10) * sigmoid * (sigmoid - 1)) * dsigma
  grad_total <- colSums(grad * em_xi)
  return(grad_total)
}

# Hessian calculation
calculate_hessian <- function(y, sigmoid, dsigma, em_xi, inside) {
  hess_num <- exp(-inside) * (y*dsigma + sigmoid^2*dsigma - 2*y*sigmoid*dsigma) * (exp(-inside)+1) + exp(-2*inside) * (-exp(inside)+1) * (y-sigmoid) * sigmoid * (sigmoid-1)
  hess_denum <- log(10) * (exp(-inside)+1)^3 * sigmoid^2 * (sigmoid-1)^2
  hess <- hess_num/hess_denum
  
  hessian_matrix <- matrix(0, ncol(em_xi), ncol(em_xi))
  for (k in 1:length(y)){
    for (i in 1:ncol(em_xi)){
      for (j in 1:ncol(em_xi)){
        hessian_matrix[i, j] <- hessian_matrix[i, j] + hess[k] * (em_xi[k,i] * em_xi[k,j])
      }
    }
  }
  return(hessian_matrix)
}

# Newton-Raphson algorithm

inside1 <- inside(combined_vector, em_xi_1)
inside2 <- inside(combined_vector, em_xi_2)
inside3 <- inside(combined_vector, em_xi_3)
  
sigmoid1 <- sigmoid(inside1)
sigmoid2 <- sigmoid(inside2)
sigmoid3 <- sigmoid(inside3)
  
dsigma1 <- dsigma(inside1)
dsigma2 <- dsigma(inside2)
dsigma3 <- dsigma(inside3)

# First derivative
grad_total_1 <- calculate_gradient(y1, sigmoid1, dsigma1, em_xi_1)
grad_total_2 <- calculate_gradient(y2, sigmoid2, dsigma2, em_xi_2)
grad_total_3 <- calculate_gradient(y3, sigmoid3, dsigma3, em_xi_3)
gradients <- grad_total_1 + grad_total_2 + grad_total_3
  
# Second derivative
hess_total_1 <- calculate_hessian(y1, sigmoid1, dsigma1, em_xi_1, inside1)
hess_total_2 <- calculate_hessian(y2, sigmoid2, dsigma2, em_xi_2, inside2)
hess_total_3 <- calculate_hessian(y3, sigmoid3, dsigma3, em_xi_3, inside3)
hess_total <- hess_total_1 + hess_total_2 + hess_total_3


# Modify parameters  

hessian_inv <- solve(hess_total)
combined_vector <- combined_vector - t(hessian_inv %*% gradients)




# Print final betas
print(tail(combined_vector, 3))

# Fit Cox model using survival package for comparison
res.cox.Rounded <- coxph(Surv(time, status) ~ X1 + X2 + X3, data = dataOG)
print(coef(res.cox.Rounded))
