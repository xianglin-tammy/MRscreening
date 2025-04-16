#' 绘制分组森林图
#'
#' 该函数用于绘制带有分组信息和多个变量的森林图，支持误差条、点估计、暴露信息、OR区间与p值列的展示，并可导出为图片。
#'
#' @param results 数据框或列表，包含OR、置信区间和分组等信息。
#' @param title 图表标题。
#' @param tsize 标题字体大小。
#' @param asize 注释文本字体大小。
#' @param zero 零线位置，通常为1。
#' @param cex 字体缩放因子（暂未使用）。
#' @param boxsize 点大小缩放因子。
#' @param output_file 输出文件路径（如需保存图像）。
#' @param output_format 输出格式，可选 "png", "jpeg", "pdf"。
#' @param width 图像宽度（单位：厘米）。
#' @param height 图像高度（单位：厘米）。
#' @param dpi 分辨率，默认300。
#' @param zero_color 零线颜色。
#' @param zero_linetype 零线线型，默认 "dashed"。
#' @param point_shape 点形状（如 15 为实心正方形）。
#' @param point_color 点颜色。
#' @param errorbar_height 误差条高度。
#' @param errorbar_color 误差条颜色。
#' @param background_colors 背景交替颜色，默认白色与灰色。
#' @param or_column OR值列名。
#' @param or_loci_column OR下限列名。
#' @param or_upci_column OR上限列名。
#' @param p2_column 暴露信息列名。
#' @param p3_column OR与CI组合列名。
#' @param p4_column P值列名。
#' @param breaks X轴刻度线位置（数值向量）。
#' @param labels X轴刻度标签（与 breaks 对应）。
#' @param x_title X轴标题。
#' @param line_spacing 行间距设置（暂未使用）。
#' @param widths 各列宽度占比。
#'
#' @return 无返回值，绘图或保存图像。
#' @export


plot_forest <- function(results,
                        title = title,
                        tsize=8,
                        asize=3,
                        zero = 1,
                        cex = 2,
                        boxsize = 0.2,
                        output_file = NULL,  # 输出文件参数
                        output_format = "png",
                        width = 40,  # 图像宽度
                        height = 20, # 图像高度
                        dpi = 300,  # 分辨率
                        zero_color = "red",  # 零线颜色
                        zero_linetype = "dashed",  # 零线形状
                        point_shape = 15,  # 点形状
                        point_color = "black",  # 点颜色
                        errorbar_height = 0.2,  # 误差条高度
                        errorbar_color = "black",  # 误差条颜色
                        background_colors = c("white", "lightgrey"),  # 背景交替颜色
                        or_column = "or",   # 指定OR列名称
                        or_loci_column = "or_loci",   # 指定OR下限列名称
                        or_upci_column = "or_upci",   # 指定OR上限列名称
                        p2_column = "exposure.trait",         # p2使用的列
                        p3_column = "OR_CI",         # p3使用的列
                        p4_column = "mr.pval",         # p4使用的列
                        breaks = c(seq(min(df$mean), max(df$mean), by = (max(df$mean) - min(df$mean))/3)),
                        labels = sprintf("%.2f", seq(min(df$mean), max(df$mean), by = (max(df$mean) - min(df$mean)) / 3)),  # 格式化标签保留两位小数
                        x_title ="OR (95%CI)",
                        line_spacing = 1.2,  # 新增参数，默认略微加宽
                        widths = c(.25, 0.10, .08, .05)
) {
  # 检查数据的列名，确保列名匹配
  if (is.list(results)) {
    # 检查列是否存在
    if (!(or_column %in% names(results) && or_loci_column %in% names(results) && or_upci_column %in% names(results))) {
      stop(paste("The list must contain columns:", or_column, ",", or_loci_column, ",", or_upci_column))
    }
    df <- data.frame(mean = results[[or_column]],
                     lower = results[[or_loci_column]],
                     upper = results[[or_upci_column]])
  } else if (is.data.frame(results)) {
    # 检查列是否存在
    if (!(or_column %in% names(results) && or_loci_column %in% names(results) && or_upci_column %in% names(results))) {
      stop(paste("The data frame must contain columns:", or_column, ",", or_loci_column, ",", or_upci_column))
    }
    df <- results
  } else {
    stop("The input 'results' must be either a list or a data.frame.")
  }

  # 为 df 添加有效的行名，确保行名的数量和行数一致
  rownames(df) <- results$exposure.trait

  # 添加索引列以便填充交替颜色
  df$Index <- 1:nrow(df)

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

  # 绘制森林图
  library(ggplot2)
  p <- ggplot(df, aes(x = mean, y = Index)) +
    # 添加背景交替格子
    geom_tile(aes(x = 1, fill = factor(Index %% 2)),
              width = Inf, height = 1, alpha = 0.2) +
    scale_fill_manual(values = background_colors) +
    # 绘制误差条
    geom_errorbarh(aes(xmin = lower, xmax = upper), height = errorbar_height, color = errorbar_color) +
    # 绘制点
    geom_point(shape = point_shape, size = boxsize * 4, color = point_color) +

    # 使用geom_segment替代geom_vline，限制零线的上下范围
    annotate("segment", x = zero, xend = zero,
             y = min(df$Index), yend = max(df$Index),
             linetype = zero_linetype, color = zero_color)+
    scale_x_continuous(
      expand = c(0, 0),
      breaks = breaks,  # 自定义刻度位置
      labels = labels  # 自定义刻度标签
    ) +
    # 美化主题
    theme_minimal() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          legend.position = "none",
          axis.title.x = element_text(size =tsize),
          axis.title.y = element_blank(),  # 去除y轴标题
          axis.text.y = element_blank(),   # 去除y轴刻度文本
          axis.ticks.y = element_blank(),
          plot.margin = margin(0, 0, 0, 0),
          axis.text.x = element_text(size =tsize),  # 显示x轴文本
          axis.ticks.x = element_line(linewidth = 0.5),
          axis.line.x = element_line(linewidth = 0.5, color = "black"),
          panel.spacing = unit(1.2, "lines"))+  # 去除y轴刻度线
    labs(x = x_title) +
    scale_y_discrete(expand = expansion(add = 1))  # 移除y轴标签

  # 绘制其他图表
  library(dplyr)
  library(patchwork)
  library(tidyverse)

  p2 <- results %>%
    select({{p2_column}}) %>%
    mutate(row_number = row_number(),  # 为每一行创建一个行号
           color_group = row_number %% 2) %>%  # 每行交替颜色
    ggplot(aes(x = 1, y = row_number, fill = factor(color_group))) +  # 使用行号作为y轴
    # 绘制交替背景色
    geom_tile(width = Inf, height = 1, alpha = 0.2) +
    scale_fill_manual(values = background_colors, guide = "none") + # 指定交替颜色
    # 添加文本标签
    geom_text(aes(label = !!sym(p2_column)), size = asize, color = "black") +
    theme(panel.grid.major = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          plot.margin = unit(rep(0, 4), "cm"),
          panel.background = element_blank(),
          axis.title.y = element_blank(),  # 去除y轴标题
          axis.text.y = element_blank(),   # 去除y轴刻度文本
          axis.ticks.y = element_blank(), # 去除y轴刻度线
          panel.spacing = unit(1.2, "lines"))+  # 可以改为 1, 1.2, 1.5 试试效果


    labs(x = "Exposure trait") +
    scale_x_discrete(position = "top")+
    scale_y_discrete(expand = expansion(add = 1))


  p3 <- results %>%
    select({{p3_column}}) %>%  # 选择指定的列
    mutate(row_number = row_number(),  # 为每一行创建一个行号
           color_group = row_number %% 2) %>%  # 每行交替颜色
    ggplot(aes(x = 1, y = row_number, fill = factor(color_group))) +  # 使用行号作为y轴
    # 绘制交替背景色
    geom_tile(width = Inf, height = 1, alpha = 0.2) +
    scale_fill_manual(values = background_colors, guide = "none") + # 指定交替颜色
    # 添加文本标签
    geom_text(aes(label = !!sym(p3_column)), size = asize, color = "black") +
    theme(panel.grid.major = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          plot.margin = unit(rep(0, 4), "cm"),
          panel.background = element_blank(),
          axis.title.y = element_blank(),  # 去除y轴标题
          axis.text.y = element_blank(),   # 去除y轴刻度文本
          axis.ticks.y = element_blank(),
          panel.spacing = unit(1.2, "lines"))+   # 去除y轴刻度线
    labs(x = "OR (95%CI)") +
    scale_x_discrete(position = "top") +
    scale_y_discrete(expand = expansion(add = 1))




  library(scales)

  p4 <- results %>%
    select({{p4_column}}) %>%
    mutate(row_number = row_number(),  # 为每一行创建一个行号
           color_group = row_number %% 2) %>%  # 每行交替颜色
    ggplot(aes(x = 1, y = row_number, fill = factor(color_group))) +  # 使用行号作为y轴
    # 绘制交替背景色
    geom_tile(width = Inf, height = 1, alpha = 0.2) +
    scale_fill_manual(values = background_colors, guide = "none") + # 指定交替颜色
    # 添加文本标签并格式化为三位小数
    geom_text(aes(label = number(!!sym(p4_column), accuracy = 0.001)), size = asize, color = "black") +
    theme(panel.grid.major = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          plot.margin = unit(rep(0, 4), "cm"),
          panel.background = element_blank(),
          axis.title.y = element_blank(),  # 去除y轴标题
          axis.text.y = element_blank(),   # 去除y轴刻度文本
          axis.ticks.y = element_blank(),
          panel.spacing = unit(1.2, "lines")) +  # 去除y轴刻度线
    labs(x = "P value") +
    scale_x_discrete(position = "top") +
    scale_y_discrete(expand = expansion(add = 1))


  final_plot <- p2+p+p3+p4+
    plot_layout(widths = widths,guides = "collect")+# 坐标轴和标题
    ggtitle(title)+  # 添加标题
    theme(plot.margin = margin(0, 0, 0, 0),  # 去除外边距
          panel.spacing = unit(0, "lines"),   # 去除各个子图之间的间隙
          plot.title = element_text(hjust = 0.2, vjust = 1, face = "bold"))
  # 打印绘图
  print(final_plot)

  # 如果打开了图形设备，则关闭它
  if (!is.null(output_file)) {
    dev.off()
  }
}
