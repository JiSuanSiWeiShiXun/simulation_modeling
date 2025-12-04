# 麦当劳排队系统模拟 - Task 2

## 任务概述

设计一个综合仿真模型来估计麦当劳的等待时间，考虑多个服务窗口和取餐窗口的配置。目标是分析不同系统配置下的性能，并提供运营改进建议。

---

## 问题拆解

### 1. 系统假设

#### 1.1 顾客到达模式

基于真实快餐店运营情况，我们采用**非齐次泊松过程(NHPP)**，到达率随时间变化：

| 时间段 | 描述 | 到达率 (λ) | 理由 |
|-------|------|-----------|------|
| **10:00-11:30** | 早餐/上午 | 25-30 人/小时 | 适度的早晨客流 |
| **11:30-13:30** | 午餐高峰 | 45-50 人/小时 | 用餐高峰期 |
| **13:30-17:00** | 下午非高峰 | 8-12 人/小时 | 两餐之间低流量 |
| **17:00-19:00** | 晚餐高峰 | 40-45 人/小时 | 晚间就餐高峰 |
| **19:00-20:00** | 晚间收尾 | 15-20 人/小时 | 客流逐渐减少 |

**全天平均到达率**: 约21人/小时（与Task 1基准一致）

**到达过程**:
\\\
到达 ~ 非齐次泊松过程，到达率为 λ(t)
使用减薄算法(thinning algorithm)生成到达间隔时间
\\\

#### 1.2 服务配置方案

分析以下多种配置：

**配置 A: 单服务器串联（基准 - Task 1）**
- 1个服务窗口 (G) - 点餐付款
- 1个取餐窗口 (G) - 取餐
- 服务率: μ = 33.33/小时, μ = 20/小时

**配置 B: 多服务窗口**
- **2个服务窗口** (G, G) - 并行点餐
- 1个取餐窗口 (G)
- 服务率: μ = μ = 33.33/小时, μ = 20/小时

**配置 C: 多取餐窗口**
- 1个服务窗口 (G)
- **2个取餐窗口** (G, G) - 并行取餐
- 服务率: μ = 33.33/小时, μ = μ = 20/小时

**配置 D: 全多服务器系统**
- **2个服务窗口** (G, G)
- **2个取餐窗口** (G, G)
- 所有服务率同上

#### 1.3 顾客行为假设

1. **排队规则**: FIFO（先进先出）
2. **选队策略**: 多服务器时，顾客加入最短队列
3. **放弃**: 不考虑放弃（所有顾客都进入系统）
4. **中途离开**: 不考虑中途离开（顾客进入队列后不会离开）
5. **服务优先级**: 无优先级（所有顾客平等对待）

#### 1.4 服务时间分布

所有服务时间遵循**指数分布**：

| 服务点 | 平均服务时间 | 服务率 | 分布 |
|-------|------------|--------|------|
| 服务窗口 (G) | 1.8 分钟 | 33.33 人/小时 | Exp(μ) |
| 取餐窗口 (G) | 3.0 分钟 | 20 人/小时 | Exp(μ) |

**注意**: 当存在多个服务器时，每个服务器独立运行，服务率相同。

---

### 2. 关键事件和变量

#### 2.1 系统事件

仿真采用**事件驱动**方式，处理以下离散事件：

| 事件类型 | 描述 | 状态变化 |
|---------|------|---------|
| **E1: 顾客到达** | 新顾客进入系统 | - 加入服务队列<br>- 若服务器空闲则开始服务<br>- 生成下一个到达 |
| **E2: 服务开始 (G)** | 顾客开始点餐 | - 服务器变为忙碌<br>- 安排服务完成时间 |
| **E3: 服务完成 (G)** | 点餐完成，付款完成 | - 顾客移动到取餐队列<br>- 下一位顾客开始服务 |
| **E4: 取餐开始 (G)** | 顾客开始取餐 | - 取餐服务器变为忙碌<br>- 安排取餐完成时间 |
| **E5: 取餐完成 (G)** | 顾客取得食物 | - 顾客离开系统<br>- 下一位顾客开始取餐 |

#### 2.2 系统状态变量

**连续时间变量**:
- \	\ - 当前仿真时间（小时）
- \T_sim\ - 总仿真时长（通常为10小时）

**服务器状态变量**:
- \server1_busy[i]\ - 布尔值，服务窗口i是否忙碌
- \server2_busy[j]\ - 布尔值，取餐窗口j是否忙碌

**队列变量**:
- \Q1\ - 服务窗口队列（顾客ID列表）
- \Q2\ - 取餐窗口队列（顾客ID列表）
- \
1\ - 服务队列中的顾客数
- \
2\ - 取餐队列中的顾客数

**顾客特定变量**:
对于每位顾客 k:
- \rrival_time[k]\ - 进入系统时间
- \service_start[k]\ - 开始点餐时间
- \service_end[k]\ - 完成点餐时间
- \pickup_start[k]\ - 开始取餐时间
- \pickup_end[k]\ - 完成取餐时间
- \W1[k]\ - 在服务窗口的等待时间
- \W2[k]\ - 在取餐窗口的等待时间
- \W_total[k]\ - 系统总时间

**性能指标** (汇总):
- \W1\ - 服务窗口平均等待时间
- \W2\ - 取餐窗口平均等待时间
- \W_sys\ - 平均系统总时间
- \ρ1\ - 服务窗口利用率
- \ρ2\ - 取餐窗口利用率
- \L1\ - 服务队列平均长度
- \L2\ - 取餐队列平均长度

---

### 3. 流程图

#### 3.1 顾客流程（单个顾客视角）

\\\

   到达        λ(t) ~ NHPP

       
       

 服务队列 Q1   等待时间 W1

       
       

 服务窗口 G   服务时间 ~ Exp(μ)

       
       

 取餐队列 Q2   等待时间 W2

       
       

 取餐窗口 G   服务时间 ~ Exp(μ)

       
       

   离开      


系统总时间 = W1 + S1 + W2 + S2
\\\

#### 3.2 事件驱动仿真流程

\\\

         初始化系统                      
  - 设置 t = 0                          
  - 生成第一个到达                       
  - 所有服务器空闲                       
  - 队列为空                            

              
              
         
          事件循环     
                    
                                 
                                 
             
     找到下一个事件             
     min(t_arrival,            
         t_service,            
         t_pickup)             
             
                                
                                
                  
       推进时间                 
        t = t_next             
                  
                                
              
      事件类型？                
              
                              
      
  到达   服务  取餐       
         完成  完成       
      
                             
      
                
                
         
          t > T_sim?   
           且          
          系统为空？   
         
             
        是      否
             
      
         收集      
         结果      
      
\\\

#### 3.3 多服务器队列选择逻辑

**多个服务窗口时**:
\\\
顾客到达
    
    

 有空闲服务器吗？   

  是            否
               
  
 分配到     加入队列 
 空闲         Q1    
 服务器    

\\\

**队列分配策略**:
1. **单队列（推荐）**: 一个公共队列，由第一个空闲服务器服务
2. **最短队列**: 顾客加入最短的队列（可能导致不平衡）
3. **轮询**: 循环分配顾客

---

### 4. 程序设计

#### 4.1 程序结构

\\\
simulation_program/

 config/
    system_config.R          # 服务器配置
    arrival_patterns.R       # NHPP λ(t) 定义

 core/
    event_driven_engine.R    # 主仿真循环
    arrival_generator.R      # NHPP 到达生成
    server_manager.R         # 多服务器队列逻辑

 scenarios/
    config_A_baseline.R      # 单服务器串联
    config_B_multi_service.R # 多服务窗口
    config_C_multi_pickup.R  # 多取餐窗口
    config_D_full_multi.R    # 全多服务器

 analysis/
    performance_metrics.R    # 计算KPI
    statistical_analysis.R   # 置信区间、假设检验
    visualization.R          # 图表和绘图

 main_simulation.R             # 入口点
\\\

#### 4.2 核心算法（事件驱动）

**伪代码**:
\\\
# 初始化
t <- 0
event_calendar <- PriorityQueue()
event_calendar(Event(type=\"ARRIVAL\", time=generate_next_arrival(0)))

# 主循环
while (!(t > T_sim && system_empty())) {
    # 获取下一个事件
    event <- event_calendar()
    t <- event\
    
    # 处理事件
    if (event\ == \"ARRIVAL\") {
        customer <- create_customer(t)
        
        # 检查服务窗口
        idle_server <- find_idle_server(service_windows)
        if (!is.null(idle_server)) {
            assign_to_server(customer, idle_server)
            service_time <- rexp(1, mu1)
            event_calendar\(Event(
                type=\"SERVICE_COMPLETE\",
                time=t + service_time,
                customer=customer,
                server=idle_server
            ))
        } else {
            Q1\(customer)
        }
        
        # 生成下一个到达（如果在关门时间前）
        if (t < closing_time) {
            next_arrival <- generate_next_arrival(t)
            event_calendar\(Event(type=\"ARRIVAL\", time=next_arrival))
        }
    }
    
    else if (event\ == \"SERVICE_COMPLETE\") {
        customer <- event\
        server <- event\
        customer\ <- t
        
        # 释放服务器
        server\ <- FALSE
        
        # 检查队列是否有顾客
        if (!Q1\()) {
            next_customer <- Q1\()
            assign_to_server(next_customer, server)
            # 安排下一个完成时间
        }
        
        # 顾客移动到取餐
        idle_pickup <- find_idle_server(pickup_windows)
        if (!is.null(idle_pickup)) {
            assign_to_server(customer, idle_pickup)
            pickup_time <- rexp(1, mu2)
            event_calendar\(Event(
                type=\"PICKUP_COMPLETE\",
                time=t + pickup_time,
                customer=customer,
                server=idle_pickup
            ))
        } else {
            Q2\(customer)
        }
    }
    
    else if (event\ == \"PICKUP_COMPLETE\") {
        # 取餐完成的类似逻辑
        # 顾客离开，更新指标
        record_customer_data(customer)
    }
}

# 收集和分析结果
calculate_performance_metrics()
\\\

---

### 5. 数据收集和分析

#### 5.1 需要收集的数据

**每位顾客的数据**:
\\\
customer_data <- data.frame(
    customer_id,
    arrival_time,
    service_queue_join_time,
    service_start_time,
    service_end_time,
    pickup_queue_join_time,
    pickup_start_time,
    pickup_end_time,
    departure_time,
    wait_time_service,
    wait_time_pickup,
    total_time_in_system
)
\\\

**系统级数据**:
\\\
system_metrics <- data.frame(
    configuration,
    avg_wait_service,
    avg_wait_pickup,
    avg_total_time,
    max_wait_service,
    max_wait_pickup,
    service_utilization,
    pickup_utilization,
    avg_queue_length_service,
    avg_queue_length_pickup,
    throughput
)
\\\

#### 5.2 性能指标

**等待时间指标**:
- 服务窗口平均等待时间
- 取餐窗口平均等待时间
- 系统总平均时间

**利用率指标**:
- 服务器利用率 = 总服务时间 / (总仿真时间  服务器数量)

**队列长度指标**:
- 时间加权平均队列长度（利特尔法则）
- 观察到的最大队列长度

**服务水平指标**:
- 在目标时间内（如5分钟）服务的顾客百分比
- 第95百分位等待时间

#### 5.3 统计分析

**置信区间**:
- 运行 N = 100 次重复
- 计算每个指标的95%置信区间

**比较检验**:
- 使用配对t检验比较配置
- 多配置比较使用方差分析(ANOVA)

**敏感性分析**:
- 到达率变化 20%
- 服务率变化 20%
- 分析对等待时间的影响

---

### 6. 分析问题和洞察

#### 6.1 需要回答的关键问题

1. **瓶颈识别**:
   - 哪个阶段（服务 vs 取餐）导致更长的等待时间？
   - 瓶颈如何随不同配置变化？

2. **配置比较**:
   - 哪种配置最小化平均等待时间？
   - 增加服务器的成本效益如何？（边际效益递减？）

3. **高峰时段性能**:
   - 系统在午餐高峰（11:30-13:30）期间表现如何？
   - 系统是否达到稳态还是持续恶化？

4. **服务器利用率**:
   - 服务器是否得到有效利用？
   - 是否有可以减少的空闲时间？

5. **队列动态**:
   - 高峰时段队列增长多长？
   - 高峰后队列何时消散？

#### 6.2 预期洞察

**配置 A（基准 - 单链）**:
- 取餐窗口是瓶颈（μ = 20/小时 < λ = 21/小时）
- 取餐等待时间长，尤其是高峰期
- 取餐服务器利用率高（>90%）

**配置 B（多服务窗口）**:
- 显著减少服务等待时间
- 但取餐仍是瓶颈
- 总体改进：总时间减少约10-20%

**配置 C（多取餐窗口）**:
- 直接解决瓶颈
- 预期改进：总时间减少约40-50%
- 系统更平衡

**配置 D（全多服务器）**:
- 整体性能最佳
- 成本最高（共4个服务器）
- 相比配置C边际效益递减

#### 6.3 改进建议

基于仿真结果，我们将建议：

1. **短期改进**:
   - 高峰时段增加取餐窗口
   - 基于时段实施动态人员配置

2. **流程改进**:
   - 非高峰时段批量备餐
   - 优化厨房到窗口的工作流

3. **技术解决方案**:
   - 移动点餐减少服务窗口负荷
   - 数字显示屏提高取餐效率

4. **容量规划**:
   - 给定需求下的最优服务器数量
   - 额外服务器的成本效益分析

---

## 实施路线图

### 阶段 1: 单服务器基准（配置A）
- 实现带有NHPP到达的事件驱动仿真
- 与Task 1结果验证
- 建立基准指标

### 阶段 2: 多服务器扩展（配置B、C、D）
- 实现多服务器队列逻辑
- 测试所有配置
- 比较性能

### 阶段 3: 分析和可视化
- 生成综合性能报告
- 创建比较图表
- 统计显著性检验

### 阶段 4: 建议
- 综合发现
- 提供可操作的洞察
- 记录改进策略

---

## 预期交付成果

1. **代码**:
   - \event_driven_simulation.R\ - 核心仿真引擎
   - \configuration_scenarios.R\ - 所有4种配置
   - \nalysis_report.R\ - 统计分析和可视化

2. **数据文件**:
   - \config_A_results.csv\ - 基准性能数据
   - \config_B_results.csv\ - 多服务窗口结果
   - \config_C_results.csv\ - 多取餐窗口结果
   - \config_D_results.csv\ - 全多服务器结果
   - \comparison_summary.csv\ - 并排比较

3. **可视化**:
   - 等待时间分布（直方图）
   - 配置比较（条形图）
   - 队列长度随时间变化（折线图）
   - 利用率热图
   - 成本效益分析

4. **报告**:
   - 执行摘要
   - 方法论描述
   - 结果和发现
   - 运营改进建议

---

## 成功标准

仿真模型将被认为成功，如果它：

1.  准确建模顾客到达模式（具有真实高峰的NHPP）
2.  正确实现多服务器排队逻辑
3.  提供统计显著的结果（N  100次重复）
4.  清晰识别系统瓶颈
5.  展示各配置间可测量的改进
6.  提供有数据支持的可操作建议

---

## 参考文献

- **排队论**: Kleinrock, L. (1975). Queueing Systems, Volume 1: Theory
- **离散事件仿真**: Law, A. M. (2015). Simulation Modeling and Analysis
- **NHPP**: Ross, S. M. (2014). Introduction to Probability Models
