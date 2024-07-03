# Load necessary library
library(survival)

# Define parameters
k <- 3
nbBetas <- 3

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
intervals <- 1:data_max

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
    df[i, stacked_data$time[i]] <- 1
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

# Adam optimization https://optimization.cbe.cornell.edu/index.php?title=Adam
t <- 0
m <- 0
v <- 0
alpha <- 0.01
beta1 <- 0.9
beta2 <- 0.99
epsilon <- 1e-8

while(t < 3000){
  inside1 <- inside(combined_vector, em_xi_1)
  inside2 <- inside(combined_vector, em_xi_2)
  inside3 <- inside(combined_vector, em_xi_3)
  
  sigmoid1 <- sigmoid(inside1)
  sigmoid2 <- sigmoid(inside2)
  sigmoid3 <- sigmoid(inside3)
  
  dsigma1 <- dsigma(inside1)
  dsigma2 <- dsigma(inside2)
  dsigma3 <- dsigma(inside3)
  
  grad_total_1 <- calculate_gradient(y1, sigmoid1, dsigma1, em_xi_1)
  grad_total_2 <- calculate_gradient(y2, sigmoid2, dsigma2, em_xi_2)
  grad_total_3 <- calculate_gradient(y3, sigmoid3, dsigma3, em_xi_3)
  gradients <- grad_total_1 + grad_total_2 + grad_total_3
  
  t <- t + 1
  m <- beta1 * m + (1 - beta1) * gradients
  v <- beta2 * v + (1 - beta2) * gradients^2
  
  m_hat <- m / (1 - beta1^t)
  v_hat <- v / (1 - beta2^t)
  
  combined_vector <- combined_vector - alpha * m_hat / (sqrt(v_hat) + epsilon)
  # print(combined_vector[98:100])
}

# Print final betas
print(combined_vector[98:100])

# Fit Cox model using survival package for comparison
res.cox.Rounded <- coxph(Surv(time, status) ~ X1 + X2 + X3, data = dataOG)
print(coef(res.cox.Rounded))
