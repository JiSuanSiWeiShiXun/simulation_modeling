# Translation mapping for queuing system files
# This script converts all Chinese comments and output text to English

$translations = @{
    # Common terms
    "麦当劳排队系统模拟" = "McDonald's Queuing System Simulation"
    "麦当劳双服务窗口排队系统模拟" = "McDonald's Two-Server Queuing System Simulation"
    "事件驱动模拟法" = "Event-Driven Simulation Method"
    "事件驱动法" = "Event-Driven Method"
    "函数库" = "Function Library"
    "主程序" = "Main Program"
    "可视化" = "Visualization"
    "快速诊断" = "Quick Diagnostic"
    "对比分析" = "Comparative Analysis"
    
    # System parameters
    "系统参数" = "System Parameters"
    "顾客到达率" = "Customer arrival rate"
    "到达率" = "Arrival rate"
    "服务窗口服务率" = "Service window rate"
    "服务窗口" = "Service window"
    "取餐窗口服务率" = "Pickup window rate"
    "取餐窗口" = "Pickup window"
    "营业时间" = "Business hours"
    "营业时长" = "Business duration"
    "人/小时" = "customers/hour"
    "每小时" = "per hour"
    "分钟" = "minutes"
    "小时" = "hours"
    
    # Simulation terms
    "模拟时长" = "Simulation duration"
    "模拟" = "Simulation"
    "仿真" = "Simulation"
    "顾客" = "Customer"
    "顾客编号" = "Customer ID"
    "到达时间" = "Arrival time"
    "到达" = "Arrival"
    "离开" = "Departure"
    "队列长度" = "Queue length"
    "排队等待" = "joins queue"
    "等待时间" = "Wait time"
    "服务时间" = "Service time"
    "系统时间" = "System time"
    "系统停留时间" = "System time"
    "系统总时间" = "Total system time"
    
    # Window operations
    "窗口1" = "Window 1"
    "窗口2" = "Window 2"
    "开始" = "Start"
    "结束" = "End"
    "完成" = "Complete"
    "空闲" = "idle"
    "忙碌" = "busy"
    
    # Statistics
    "统计" = "Statistics"
    "平均" = "Average"
    "均值" = "Mean"
    "中位数" = "Median"
    "标准差" = "Standard deviation"
    "方差" = "Variance"
    "置信区间" = "Confidence interval"
    "百分位数" = "Percentile"
    "分位数" = "Percentile"
    "最大" = "Maximum"
    "最小" = "Minimum"
    "总" = "Total"
    
    # Task descriptions
    "任务" = "Task"
    "模拟2小时运营" = "Simulate 2-hour operation"
    "估计10小时内顾客平均系统时间" = "Estimate average customer system time over 10 hours"
    "估计加班时间" = "Estimate overtime"
    "加班时间" = "Overtime"
    "不需要加班" = "No overtime"
    "需要加班" = "Overtime required"
    
    # Output messages
    "生成" = "Generating"
    "保存" = "Saved"
    "已保存至" = "saved to"
    "详细结果" = "Detailed results"
    "详细记录" = "Detailed records"
    "汇总结果" = "Summary results"
    "可视化图表" = "visualization charts"
    "图表" = "Charts"
    "完成" = "Complete"
    "运行" = "Running"
    "次独立模拟" = "independent simulations"
    
    # Chart titles
    "顾客到达与离开" = "Customer Arrivals and Departures"
    "累计顾客数" = "Cumulative Customers"
    "系统中顾客数量" = "Customers in System"
    "系统停留时间分布" = "System Time Distribution"
    "等待时间比较" = "Wait Time Comparison"
    "两个窗口等待时间比较" = "Wait Time Comparison Between Two Windows"
    "服务时间对比" = "Service Time Comparison"
    "时间成分堆叠" = "Time Components Stacked"
    "顾客服务历程甘特图" = "Customer Service Journey Gantt Chart"
    "累积分布函数" = "Cumulative Distribution Function"
    "分布" = "Distribution"
    "频数" = "Frequency"
    "时间" = "Time"
    
    # File names
    "文件" = "Files"
    "生成的文件" = "Generated files"
    
    # Misc
    "说明" = "Description"
    "注" = "Note"
    "方法" = "Method"
    "使用" = "Using"
    "理论" = "Theoretical"
    "实际" = "Actual"
    "估计值" = "Estimate"
    "样本数量" = "Sample size"
    "极值" = "Extreme values"
    "概率" = "Probability"
}

Write-Host "Translation mapping loaded with $($translations.Count) entries"
Write-Host "Use this mapping to translate R files from Chinese to English"
