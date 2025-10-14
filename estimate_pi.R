# 方法 1：Monte Carlo 积分法
estimate_pi_integral <- function(N) {
  # 生成 N 个 [0,1] 区间内的均匀随机数
  x <- runif(N, min = 0, max = 1)
  
  # 对每个随机数 x 计算函数 f(x) = 4 / (1 + x^2)
  fx <- 4 / (1 + x*x)
  
  # 取这些 f(x) 的平均值作为 π 的估计值
  pi_hat <- mean(fx)
  
  return(pi_hat)
}

# 方法 2：Buffon's Needle 法
estimate_pi_buffon <- function(N, L = 1, t = 1) {
  # 检查条件：针的长度 L 必须小于等于间距 t
  if (L > t) stop("L must be <= t")
  
  # 随机生成 N 个针与平行线的夹角 θ
  theta <- runif(N, min = 0, max = pi/2)
  
  # 随机生成针中心到最近线的距离 d，d 服从 [0, t/2] 的均匀分布
  d <- runif(N, min = 0, max = t/2)
  
  # 判断针是否跨线，若 d <= (L/2)*sin(θ)，表示针跨线
  hits <- d <= (L/2) * sin(theta)
  
  # 计算跨线的比例 p_hat
  p_hat <- mean(hits)
  
  # 根据公式 π ≈ (2L)/(t*p_hat)
  pi_hat <- (2 * L) / (t * p_hat)
  
  return(pi_hat)
}

# 方法 3：Random Chord 法
estimate_pi_chord <- function(N) {
  # 在 [0, 2π] 内随机生成 N 对角度
  theta1 <- runif(N, min = 0, max = 2*pi)
  theta2 <- runif(N, min = 0, max = 2*pi)
  
  # 计算两点间的角度差
  dtheta <- abs(theta1 - theta2)
  
  # 若角度差超过 π，取 2π - dtheta（因为弦长取较短那条）
  dtheta <- ifelse(dtheta > pi, 2*pi - dtheta, dtheta)
  
  # 计算弦长 L
  L <- 2 * sin(dtheta / 2)
  
  # 根据期望公式 π = 4 / mean(L)
  pi_hat <- 4 / mean(L)
  
  return(pi_hat)
}


# # 运行例子
# set.seed(2025)  # 固定随机数种子，注释不固定

# 运行方法 1：Monte Carlo 积分法
cat("Method 1 (Monte Carlo Integration) - π estimation: ", estimate_pi_integral(10000), "\n")

# 运行方法 2：Buffon's Needle 法
cat("Method 2 (Buffon's Needle) - π estimation: ", estimate_pi_buffon(100000), "\n")

# 运行方法 3：Random Chord 法
cat("Method 3 (Random Chord) - π estimation: ", estimate_pi_chord(100000), "\n")