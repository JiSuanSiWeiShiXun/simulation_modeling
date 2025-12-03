# ==========================================
# 生成报告用的可视化表格
# ==========================================

# 读取数据
data <- read.csv("event_driven_2hr_results.csv")

cat("\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("生成Task a报告表格\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

# ==========================================
# 1. 生成HTML格式的专业表格
# ==========================================

# 准备显示数据（转换为分钟，保留2位小数）
display_data <- data.frame(
  顾客编号 = data$customer_id,
  到达时间 = sprintf("%.2f", data$arrival_time * 60),
  窗口1开始 = sprintf("%.2f", data$server1_start * 60),
  窗口1结束 = sprintf("%.2f", data$server1_end * 60),
  窗口2开始 = sprintf("%.2f", data$server2_start * 60),
  窗口2结束 = sprintf("%.2f", data$server2_end * 60),
  窗口1等待 = sprintf("%.2f", data$server1_wait * 60),
  窗口1服务 = sprintf("%.2f", data$server1_service * 60),
  窗口2等待 = sprintf("%.2f", data$server2_wait * 60),
  窗口2服务 = sprintf("%.2f", data$server2_service * 60),
  系统总时间 = sprintf("%.2f", data$total_time_in_system * 60)
)

# 计算统计值
total_customers <- nrow(data)
avg_system_time <- mean(data$total_time_in_system * 60)
avg_wait_time <- mean(data$total_wait_time * 60)
avg_wait1 <- mean(data$server1_wait * 60)
avg_wait2 <- mean(data$server2_wait * 60)
last_departure <- max(data$server2_end) * 60

# 生成HTML - 分段构建
html_header <- sprintf('<!DOCTYPE html>
<html>
<head>
    <meta charset="UTF-8">
    <title>Task a - 2小时运营详细记录</title>
    <style>
        body {
            font-family: "Microsoft YaHei", "Segoe UI", Arial, sans-serif;
            margin: 20px;
            background-color: #f5f5f5;
        }
        .container {
            max-width: 1400px;
            margin: 0 auto;
            background-color: white;
            padding: 30px;
            box-shadow: 0 2px 8px rgba(0,0,0,0.1);
            border-radius: 8px;
        }
        h1 {
            color: #2c3e50;
            text-align: center;
            border-bottom: 3px solid #3498db;
            padding-bottom: 15px;
        }
        .summary {
            background-color: #ecf0f1;
            padding: 20px;
            margin: 20px 0;
            border-radius: 5px;
            border-left: 5px solid #3498db;
        }
        .summary h2 {
            color: #2c3e50;
            margin-top: 0;
        }
        .summary-grid {
            display: grid;
            grid-template-columns: repeat(3, 1fr);
            gap: 15px;
            margin-top: 15px;
        }
        .summary-item {
            background-color: white;
            padding: 12px;
            border-radius: 4px;
            text-align: center;
        }
        .summary-label {
            font-size: 12px;
            color: #7f8c8d;
            margin-bottom: 5px;
        }
        .summary-value {
            font-size: 20px;
            font-weight: bold;
            color: #2c3e50;
        }
        table {
            width: 100%%;
            border-collapse: collapse;
            margin: 20px 0;
            font-size: 13px;
        }
        thead {
            background: linear-gradient(135deg, #667eea 0%%, #764ba2 100%%);
            color: white;
        }
        th {
            padding: 12px 8px;
            text-align: center;
            font-weight: 600;
            border: 1px solid #ddd;
        }
        td {
            padding: 10px 8px;
            text-align: center;
            border: 1px solid #ddd;
        }
        tbody tr:nth-child(odd) {
            background-color: #f8f9fa;
        }
        tbody tr:hover {
            background-color: #e3f2fd;
            transition: background-color 0.3s;
        }
        .highlight-time {
            font-weight: bold;
            color: #e74c3c;
        }
        .highlight-wait {
            background-color: #fff3cd;
        }
        .highlight-service {
            background-color: #d1ecf1;
        }
        .footer {
            margin-top: 30px;
            text-align: center;
            color: #7f8c8d;
            font-size: 12px;
            border-top: 1px solid #ddd;
            padding-top: 15px;
        }
        .note {
            background-color: #fff3cd;
            border-left: 4px solid #ffc107;
            padding: 15px;
            margin: 20px 0;
            border-radius: 4px;
        }
        .note-title {
            font-weight: bold;
            color: #856404;
            margin-bottom: 8px;
        }
    </style>
</head>
<body>
    <div class="container">
        <h1>Task a) 麦当劳排队系统 - 2小时运营详细记录</h1>
        
        <div class="summary">
            <h2>📊 运营统计摘要</h2>
            <div class="summary-grid">
                <div class="summary-item">
                    <div class="summary-label">总顾客数</div>
                    <div class="summary-value">%d 人</div>
                </div>
                <div class="summary-item">
                    <div class="summary-label">平均系统时间</div>
                    <div class="summary-value">%.2f 分钟</div>
                </div>
                <div class="summary-item">
                    <div class="summary-label">平均等待时间</div>
                    <div class="summary-value">%.2f 分钟</div>
                </div>
                <div class="summary-item">
                    <div class="summary-label">窗口1平均等待</div>
                    <div class="summary-value">%.2f 分钟</div>
                </div>
                <div class="summary-item">
                    <div class="summary-label">窗口2平均等待</div>
                    <div class="summary-value">%.2f 分钟</div>
                </div>
                <div class="summary-item">
                    <div class="summary-label">最后离开时间</div>
                    <div class="summary-value">%.2f 分钟</div>
                </div>
            </div>
        </div>
        
        <div class="note">
            <div class="note-title">⚠️ 说明</div>
            所有时间单位均为分钟。黄色背景表示等待时间，蓝色背景表示服务时间。粉色行表示总等待超过10分钟。
        </div>
        
        <table>
            <thead>
                <tr>
                    <th rowspan="2">顾客<br>编号</th>
                    <th rowspan="2">到达<br>时间</th>
                    <th colspan="2">服务窗口</th>
                    <th colspan="2">取餐窗口</th>
                    <th colspan="4">时间成分（分钟）</th>
                    <th rowspan="2">系统<br>总时间</th>
                </tr>
                <tr>
                    <th>开始</th>
                    <th>结束</th>
                    <th>开始</th>
                    <th>结束</th>
                    <th>窗口1<br>等待</th>
                    <th>窗口1<br>服务</th>
                    <th>窗口2<br>等待</th>
                    <th>窗口2<br>服务</th>
                </tr>
            </thead>
            <tbody>
', total_customers, avg_system_time, avg_wait_time, avg_wait1, avg_wait2, last_departure)

# 写入文件头
cat(html_header, file = "task_a_report_table.html")

# 添加数据行
for (i in 1:nrow(display_data)) {
  row <- display_data[i, ]
  
  # 判断是否高亮
  is_long_wait <- data$total_wait_time[i] * 60 > 10
  row_style <- if(is_long_wait) ' style="background-color: #ffe6e6;"' else ''
  
  html_row <- sprintf('
                <tr%s>
                    <td><strong>%s</strong></td>
                    <td>%s</td>
                    <td>%s</td>
                    <td>%s</td>
                    <td>%s</td>
                    <td>%s</td>
                    <td class="highlight-wait">%s</td>
                    <td class="highlight-service">%s</td>
                    <td class="highlight-wait">%s</td>
                    <td class="highlight-service">%s</td>
                    <td class="highlight-time">%s</td>
                </tr>',
    row_style,
    row$顾客编号, row$到达时间, row$窗口1开始, row$窗口1结束,
    row$窗口2开始, row$窗口2结束, row$窗口1等待, row$窗口1服务,
    row$窗口2等待, row$窗口2服务, row$系统总时间)
  
  cat(html_row, file = "task_a_report_table.html", append = TRUE)
}

# 写入文件尾
html_footer <- sprintf('
            </tbody>
        </table>
        
        <div class="footer">
            <p>生成时间: %s</p>
            <p>数据来源: event_driven_2hr_results.csv | 模拟方法: 事件驱动模拟</p>
        </div>
    </div>
</body>
</html>
', format(Sys.time(), "%Y-%m-%d %H:%M:%S"))

cat(html_footer, file = "task_a_report_table.html", append = TRUE)
cat("✓ HTML格式表格已保存至: task_a_report_table.html\n")

# ==========================================
# 2. 生成LaTeX格式表格
# ==========================================

latex_data <- display_data[1:min(30, nrow(display_data)), ]

latex_content <- sprintf('\\documentclass[12pt, a4paper]{article}
\\usepackage[UTF8]{ctex}
\\usepackage{booktabs}
\\usepackage{longtable}
\\usepackage{xcolor}
\\usepackage{colortbl}
\\usepackage[margin=2cm]{geometry}

\\begin{document}

\\title{\\textbf{Task a) 麦当劳排队系统模拟报告}}
\\author{2小时运营详细记录}
\\date{%s}
\\maketitle

\\section{运营统计摘要}

\\begin{itemize}
    \\item 总顾客数: %d 人
    \\item 平均系统时间: %.2f 分钟
    \\item 平均等待时间: %.2f 分钟
    \\item 窗口1平均等待: %.2f 分钟
    \\item 窗口2平均等待: %.2f 分钟
\\end{itemize}

\\section{顾客详细记录（前30位）}

\\begin{table}[h]
\\centering
\\small
\\begin{tabular}{ccccccccccc}
\\toprule
\\textbf{顾客} & \\textbf{到达} & \\multicolumn{2}{c}{\\textbf{服务窗口}} & \\multicolumn{2}{c}{\\textbf{取餐窗口}} & \\multicolumn{4}{c}{\\textbf{时间成分（分钟）}} & \\textbf{系统} \\\\
\\textbf{编号} & \\textbf{时间} & 开始 & 结束 & 开始 & 结束 & 窗口1 & 窗口1 & 窗口2 & 窗口2 & \\textbf{总时间} \\\\
& & & & & & 等待 & 服务 & 等待 & 服务 & \\\\
\\midrule
', format(Sys.Date(), "%Y年%m月%d日"), total_customers, avg_system_time, avg_wait_time, avg_wait1, avg_wait2)

for (i in 1:nrow(latex_data)) {
  row <- latex_data[i, ]
  latex_content <- paste0(latex_content, sprintf('%s & %s & %s & %s & %s & %s & %s & %s & %s & %s & %s \\\\\n',
                          row$顾客编号, row$到达时间, row$窗口1开始, row$窗口1结束,
                          row$窗口2开始, row$窗口2结束, row$窗口1等待, row$窗口1服务,
                          row$窗口2等待, row$窗口2服务, row$系统总时间))
}

latex_content <- paste0(latex_content, sprintf('\\bottomrule
\\end{tabular}
\\caption{前30位顾客的详细服务记录（时间单位：分钟）}
\\end{table}

\\textbf{注：}完整的%d位顾客数据请参见 event\\_driven\\_2hr\\_results.csv 文件。

\\end{document}
', total_customers))

writeLines(latex_content, "task_a_report_table.tex")
cat("✓ LaTeX格式表格已保存至: task_a_report_table.tex\n")

# ==========================================
# 3. 生成纯文本格式表格
# ==========================================

cat("\n")
cat(paste(rep("=", 140), collapse = ""), "\n")
cat("Task a) 顾客详细记录表格（前30位）\n")
cat(paste(rep("=", 140), collapse = ""), "\n")
cat(sprintf("%-6s %-8s %-8s %-8s %-8s %-8s | %-8s %-8s %-8s %-8s | %-10s\n",
            "顾客", "到达", "窗口1", "窗口1", "窗口2", "窗口2",
            "窗口1", "窗口1", "窗口2", "窗口2", "系统"))
cat(sprintf("%-6s %-8s %-8s %-8s %-8s %-8s | %-8s %-8s %-8s %-8s | %-10s\n",
            "编号", "时间", "开始", "结束", "开始", "结束",
            "等待", "服务", "等待", "服务", "总时间"))
cat(paste(rep("-", 140), collapse = ""), "\n")

for (i in 1:min(30, nrow(display_data))) {
  row <- display_data[i, ]
  cat(sprintf("%-6s %-8s %-8s %-8s %-8s %-8s | %-8s %-8s %-8s %-8s | %-10s\n",
              row$顾客编号, row$到达时间, row$窗口1开始, row$窗口1结束,
              row$窗口2开始, row$窗口2结束, row$窗口1等待, row$窗口1服务,
              row$窗口2等待, row$窗口2服务, row$系统总时间))
}
cat(paste(rep("=", 140), collapse = ""), "\n")
cat(sprintf("\n注：完整数据共 %d 位顾客，上表显示前30位\n", nrow(data)))

# ==========================================
# 4. 生成可视化表格图片
# ==========================================

png("task_a_report_table.png", width = 1600, height = 1600, res = 100)
par(mar = c(1, 1, 3, 1))

plot_data <- display_data[1:min(30, nrow(display_data)), ]

plot.new()
plot.window(xlim = c(0, 1), ylim = c(0, 1))

text(0.5, 0.97, "Task a) McDonald's Queuing System - 2-Hour Operation Detailed Records", 
     cex = 1.5, font = 2, col = "#2c3e50")

summary_y <- 0.91
text(0.5, summary_y, sprintf("Total Customers: %d | Avg System Time: %.2f min | Avg Wait Time: %.2f min | Window 2 Avg Wait: %.2f min",
                             total_customers, avg_system_time, avg_wait_time, avg_wait2),
     cex = 0.9, col = "#34495e")

n_rows <- nrow(plot_data)
n_cols <- 11
cell_height <- 0.8 / (n_rows + 2)
col_widths <- c(0.06, 0.08, 0.08, 0.08, 0.08, 0.08, 0.08, 0.08, 0.08, 0.08, 0.10)
start_x <- 0.02
start_y <- 0.85

y_pos <- start_y
x_pos <- start_x

headers <- c("Customer\nID", "Arrival\nTime", "Window 1\nStart", "Window 1\nEnd", 
             "Window 2\nStart", "Window 2\nEnd", "Window 1\nWait", "Window 1\nService",
             "Window 2\nWait", "Window 2\nService", "Total\nSystem Time")

for (j in 1:n_cols) {
  rect(x_pos, y_pos - cell_height, x_pos + col_widths[j], y_pos, 
       col = "#667eea", border = "white", lwd = 2)
  text(x_pos + col_widths[j]/2, y_pos - cell_height/2, headers[j], 
       cex = 0.65, col = "white", font = 2)
  x_pos <- x_pos + col_widths[j]
}

for (i in 1:n_rows) {
  y_pos <- y_pos - cell_height
  x_pos <- start_x
  
  bg_col <- if (i %% 2 == 1) "#f8f9fa" else "white"
  
  if (data$total_wait_time[i] * 60 > 10) {
    bg_col <- "#ffe6e6"
  }
  
  row_data <- as.character(plot_data[i, ])
  
  for (j in 1:n_cols) {
    cell_bg <- bg_col
    if (j == 7 || j == 9) cell_bg <- "#fff3cd"
    if (j == 8 || j == 10) cell_bg <- "#d1ecf1"
    
    rect(x_pos, y_pos - cell_height, x_pos + col_widths[j], y_pos, 
         col = cell_bg, border = "#ddd")
    
    text_col <- if (j == 11) "#e74c3c" else "black"
    text(x_pos + col_widths[j]/2, y_pos - cell_height/2, row_data[j], 
         cex = 0.6, col = text_col, font = if (j == 11) 2 else 1)
    
    x_pos <- x_pos + col_widths[j]
  }
}

text(0.5, 0.02, 
     sprintf("Note: All times in minutes. Yellow=Wait time, Blue=Service time, Pink rows=Wait >10 min. Total %d customers, showing first 30.", total_customers),
     cex = 0.7, col = "#7f8c8d")

dev.off()
cat("✓ Table image saved to: task_a_report_table.png\n")

cat("\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("报告表格生成完成！\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

cat("生成的文件:\n")
cat("  1. task_a_report_table.html - 交互式HTML表格（推荐用于网页报告）\n")
cat("  2. task_a_report_table.tex  - LaTeX格式表格（用于学术论文）\n")
cat("  3. task_a_report_table.png  - 表格图片（用于PPT/Word）\n")
cat("  4. 控制台输出              - 纯文本格式（用于快速查看）\n\n")

cat("使用建议:\n")
cat("  • HTML格式: 用浏览器打开，可以交互查看，适合电子报告\n")
cat("  • LaTeX格式: 用于学术论文，需要编译（xelatex task_a_report_table.tex）\n")
cat("  • PNG图片: 可直接插入Word/PPT，无需额外处理\n\n")
