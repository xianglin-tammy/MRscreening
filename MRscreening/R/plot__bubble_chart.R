
#' @title 绘制气泡图
#' @description 根据 MR 分析结果绘制气泡图，展示 Beta 值与 P 值的关系，并通过颜色和大小标记显著性和方向。
#' @param results 数据框，包含 MR 分析结果
#' @param p_column 字符串，表示 P 值所在的列名，默认 "mr.pval"
#' @param x_column 字符串，表示 Beta 值所在的列名，默认 "mr.b"
#' @param title 图表标题
#' @param x_lab x 轴标签，默认 "Beta"
#' @param y_lab y 轴标签，默认 "-log10(P值)"
#' @param labtext 图例中显示的大小说明，默认 "Beta绝对值"
#' @param name 用于显示标签的列名
#' @param col_func_negative 颜色映射函数（负方向）
#' @param col_func_positive 颜色映射函数（正方向）
#' @param width 图宽（单位：cm）
#' @param height 图高（单位：cm）
#' @param dpi 分辨率
#' @param output_file 输出文件路径，默认为 NULL 不输出
#' @param output_format 输出格式："png"、"jpeg" 或 "pdf"
#' @return 无返回值，绘图并打印或保存
#' @export
#' @name plot_bubble_chart
plot_bubble_chart <- function(results,
                              p_column="mr.pval",
                              x_column="mr.b",
                              title = "title",
                              x_lab="Beta",
                              y_lab="-log10(P值)",
                              labtext="Beta绝对值",
                              name="exposure.trait",
                              col_func_negative= colorRampPalette(c("lightblue1", "deepskyblue4")),
                              col_func_positive= colorRampPalette(c("pink", "firebrick4")),
                              width = 40,  # 图像宽度
                              height = 20, # 图像高度
                              dpi = 300,  # 分辨率
                              output_file = NULL,  # 输出文件参数
                              output_format = "png") {

  # 计算p值的负对数转换
  results$neg_log_p <- -log10(results[[p_column]])

  # 根据beta的绝对值设定气泡大小
  results$size <- abs(results[[x_column]]) * 10  # 可根据需要调整倍数

  # 通过颜色渐变的范围来映射beta值
  # 对于负值（x < 0），使用蓝色渐变
  results$color <- ifelse(results[[x_column]] < 0,
                          col_func_negative(100)[findInterval(abs(results[[x_column]]),
                                                          seq(0, max(abs(results[[x_column]])), length.out = 100))],
                          # 对于正值（x > 0），使用红色渐变
                          col_func_positive(100)[findInterval(abs(results[[x_column]]),
                                                         seq(0, max(abs(results[[x_column]])), length.out = 100))])
  # 如果指定了输出文件，打开图形设备
  if (!is.null(output_file)) {
    dir_path <- dirname(output_file)
    if (!dir.exists(dir_path)) {
      stop("指定的目录不存在: ", dir_path)
    }

    if (output_format == "png") {
      png(filename = output_file, width = width * 100, height = height * 100, res = dpi)
    } else if (output_format == "jpeg") {
      jpeg(filename = output_file, width = width * 100, height = height * 100, res = dpi)
    } else if (output_format == "pdf") {
      pdf(filename = output_file, width = width, height = height)
    } else {
      stop("不支持的图形格式，请选择 'png', 'jpeg' 或 'pdf'.")
    }
  }
  library(ggplot2)


  # 筛选出beta绝对值大于0.5的数据
  results$label <- ifelse(abs(results[[x_column]]) > 0.1, results[[name]], NA)

  # 绘制气泡图
p<- ggplot(results, aes(x = .data[[x_column]], y = neg_log_p, size = size, color = color)) +
    geom_point(alpha = 0.9, shape = 21, fill = results$color, color = "black") +  # 设置气泡透明度
    scale_size_continuous(range = c(3, 15)) +  # 设置气泡大小范围
    scale_x_continuous(limits = c(-1, 1)) +
    scale_y_continuous(limits = c(-100, NA),expand = expansion(mult = c(-0.1, 1)))+  # 调整扩展参数
    scale_color_identity() +  # 使用颜色
    labs(
      title = title,
      x = x_lab,
      y = y_lab,
      size = labtext
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      panel.grid = element_blank(),  # 移除背景网格
      panel.border = element_rect(color = "black", fill=NA, linewidth = 1),  # 添加面板边框
      plot.margin = margin(10, 10, 10, 10),
      axis.line = element_line(color = "black", linewidth = 0.5),
      axis.ticks = element_line(color = "black", linewidth = 0.5)  # 显示刻度线
    ) +
  # 自定义图例的样式
  guides(
    color = guide_legend(override.aes = list(shape = 21, fill = "white", color = "black")),  # 设置图例为白色圆形带黑色边框
    size = guide_legend(override.aes = list(shape = 21, fill = "white", color = "black"))  # 设置图例为白色圆形带黑色边框
  )+
  # 添加exposure标签
  geom_text(aes(label = label), na.rm = TRUE, vjust = -1, size = 3.5, color = "black")+  # 仅显示beta绝对值大于0.5的exposure名
  # 添加 P 值显著性阈值的虚线（例如 p = 0.05）
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey30", linewidth = 0.8)
  print(p)

  # 如果打开了图形设备，则关闭它
  if (!is.null(output_file)) {
    dev.off()
  }
}


