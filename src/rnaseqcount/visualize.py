# visualize.py
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from .core import CountData
from .stats import pca, sample_correlation, group_samples
from .algorithms import differential_test

def plot_pca(countdata: CountData, method: str = "TPM", filter_threshold: float = 0,
             label_name: str = None, grouped_by: str = None, groups: list = None,
             figsize: tuple = (10, 8), show_labels: bool = True):
    """
    绘制PCA散点图，支持按分组着色（与差异分析分组逻辑一致）。
    """
    from .stats import group_samples
    
    # 获取 PCA 坐标
    pca_df, pca_var_ratio = pca(countdata, method=method, n_components=2,
                                 filter_threshold=filter_threshold, label_name=label_name)
    
    # 如果指定了分组，计算每个样本的颜色组别
    if grouped_by is not None or groups is not None:
        group_dict = group_samples(countdata, grouped_by=grouped_by, groups=groups)
        # 构建样本名 -> 组名的映射
        sample_to_group = {}
        for grp, samples in group_dict.items():
            for s in samples:
                sample_to_group[s] = grp
        # 按照 pca_df 的索引顺序（即样本顺序）生成颜色标签
        color_labels = [sample_to_group.get(sample, 'Unknown') for sample in pca_df.index]
        # 检查是否只有一种分组（若所有样本属于同一组，则降级为单一颜色）
        if len(set(color_labels)) > 1:
            unique_groups = sorted(set(color_labels))
            # 使用 tab10 调色板，若分组数超过10则循环
            colors = plt.cm.tab10(np.linspace(0, 1, len(unique_groups)))
            color_map = {grp: colors[i] for i, grp in enumerate(unique_groups)}
            point_colors = [color_map[g] for g in color_labels]
            
            fig, ax = plt.subplots(figsize=figsize)
            ax.scatter(pca_df['PC1'], pca_df['PC2'], c=point_colors, s=50, alpha=0.7)
            # 添加图例
            from matplotlib.patches import Patch
            legend_elements = [Patch(facecolor=color_map[g], label=g) for g in unique_groups]
            ax.legend(handles=legend_elements, title=grouped_by or 'groups', loc='best')
        else:
            # 单一颜色（所有样本同一组或分组失败）
            fig, ax = plt.subplots(figsize=figsize)
            ax.scatter(pca_df['PC1'], pca_df['PC2'], s=50, alpha=0.7, color='steelblue')
    else:
        # 未指定分组，使用单一蓝色
        fig, ax = plt.subplots(figsize=figsize)
        ax.scatter(pca_df['PC1'], pca_df['PC2'], s=50, alpha=0.7, color='steelblue')
    
    ax.set_xlabel(f'PC1 ({pca_var_ratio[0]:.2%} variance)', fontsize=12)
    ax.set_ylabel(f'PC2 ({pca_var_ratio[1]:.2%} variance)', fontsize=12)
    ax.set_title(f'PCA Plot ({method}, threshold={filter_threshold})', fontsize=14, fontweight='bold')
    
    if show_labels:
        for sample, row in pca_df.iterrows():
            ax.annotate(sample, (row['PC1'], row['PC2']),
                        xytext=(5, 5), textcoords='offset points', fontsize=8)
    
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    return fig

def plot_sample_correlation(countdata: CountData, method: str = "TPM", filter_threshold: float = 0, 
                            label_name: str = None, cmap: str = 'coolwarm'):
    """
    绘制样本相关性热图，自动调整图形尺寸和标签字体以确保样本名完整显示。
    
    参数
    ----------
    countdata : CountData
        数据对象
    method : str
        归一化方法 ("TPM"/"FPKM"/"CPM")
    filter_threshold : float
        低表达过滤阈值
    label_name : str, optional
        metadata 中的列名，虽然传入此参数（用于内部计算相关性矩阵时可能重命名样本），
        但热图的坐标轴标签始终显示原始样本名（countdata.sample_names），避免重复标签重叠。
    cmap : str
        颜色映射
    
    返回
    -------
    matplotlib.figure.Figure
    """
    corr_matrix = sample_correlation(countdata, method=method, filter_threshold=filter_threshold, 
                                      label_name=label_name)
    n_samples = corr_matrix.shape[0]
    
    # 获取原始样本名（强制用于显示）
    original_sample_names = countdata.sample_names
    # 确保长度匹配（过滤不影响样本，但安全起见）
    if len(original_sample_names) != n_samples:
        # 若因过滤导致样本数变化（实际不会），则截取或使用 corr_matrix 的索引
        original_sample_names = original_sample_names[:n_samples]
    
    # 动态调整图形尺寸：每10个样本增加2英寸宽度和高度，最小8x6
    width = max(8, n_samples * 0.6)
    height = max(6, n_samples * 0.5)
    fig, ax = plt.subplots(figsize=(width, height))
    
    im = ax.imshow(corr_matrix, cmap=cmap, vmin=-1, vmax=1, aspect='auto')
    plt.colorbar(im, ax=ax, label='Pearson Correlation')
    
    # 设置刻度位置
    ax.set_xticks(np.arange(n_samples))
    ax.set_yticks(np.arange(n_samples))
    
    # 使用原始样本名作为标签
    x_labels = original_sample_names
    y_labels = original_sample_names
    
    # 根据样本数量调整字体大小
    base_fontsize = 10
    if n_samples > 30:
        base_fontsize = 6
    elif n_samples > 20:
        base_fontsize = 8
    
    # x轴标签旋转45度，右对齐，避免重叠
    ax.set_xticklabels(x_labels, rotation=45, ha='right', fontsize=base_fontsize)
    ax.set_yticklabels(y_labels, fontsize=base_fontsize)
    
    # 显示相关系数数值（样本数 ≤20 时显示，否则隐藏以保持可读性）
    if n_samples <= 20:
        for i in range(n_samples):
            for j in range(n_samples):
                ax.text(j, i, f'{corr_matrix.iloc[i, j]:.2f}', 
                        ha="center", va="center", color="black", fontsize=8)
    
    ax.set_title(f'Sample Correlation Heatmap ({method}, threshold={filter_threshold})', 
                 fontsize=14, fontweight='bold')
    ax.set_xlabel('Samples', fontsize=12)
    ax.set_ylabel('Samples', fontsize=12)
    
    plt.tight_layout()
    return fig

def plot_differential_test(countdata: CountData, method: str = "TPM", filter_threshold: float = 0, 
                           grouped_by: str = None, groups: list = None, 
                           log2_transform: bool = True, fdr_correction: bool = True, 
                           test_method: str = "t_test", paired: bool = False, 
                           significance_threshold: float = 0.05, log2fc_threshold: float = 2):
    ## 绘制火山图(Volcano Plot)展示差异表达结果
    ## significance_threshold: 显著性阈值
    ## log2fc_threshold: log2 fold change阈值，用于定义上调/下调基因
    ## 返回: matplotlib Figure对象
    result_df = differential_test(countdata, method=method, filter_threshold=filter_threshold, 
                                   grouped_by=grouped_by, groups=groups, 
                                   log2_transform=log2_transform, fdr_correction=fdr_correction, 
                                   test_method=test_method, paired=paired)
    fig, ax = plt.subplots(figsize=(10, 8))

    if fdr_correction:
        p_col = "fdr"
    else:
        p_col = "p_values"

    log2fc = result_df["log2_fold_changes"]
    p_vals = result_df[p_col]

    significant = p_vals < significance_threshold
    up_regulated = significant & (log2fc > log2fc_threshold)
    down_regulated = significant & (log2fc < -log2fc_threshold)
    not_significant = ~significant | ((log2fc >= -log2fc_threshold) & (log2fc <= log2fc_threshold))

    ax.scatter(log2fc[not_significant], -np.log10(p_vals[not_significant]), 
               alpha=0.5, s=10, c="grey", label='Not significant')
    ax.scatter(log2fc[up_regulated], -np.log10(p_vals[up_regulated]), 
               alpha=0.7, s=15, c="red", label=f'Up-regulated (log2FC > {log2fc_threshold})')
    ax.scatter(log2fc[down_regulated], -np.log10(p_vals[down_regulated]), 
               alpha=0.7, s=15, c="blue", label=f'Down-regulated (log2FC < -{log2fc_threshold})')

    ax.axhline(y=-np.log10(significance_threshold), color='red', linestyle='--', alpha=0.5)
    ax.axvline(x=log2fc_threshold, color='red', linestyle='--', alpha=0.5)
    ax.axvline(x=-log2fc_threshold, color='red', linestyle='--', alpha=0.5)

    ax.set_xlabel('log2(Fold Change)', fontsize=12, fontweight='bold')
    ax.set_ylabel('-log10(p-value)', fontsize=12, fontweight='bold')

    group_dict = group_samples(countdata, grouped_by=grouped_by, groups=groups)
    groups_list = list(group_dict.keys())
    test_name = "Paired " + test_method if paired else test_method
    correction_text = " (FDR corrected)" if fdr_correction else ""
    title = f'Volcano Plot - {test_name}{correction_text}\nGroups: {groups_list[0]} vs {groups_list[1]}'
    ax.set_title(title, fontsize=14, fontweight='bold')
    ax.legend(loc='best', frameon=True, fancybox=True, shadow=True)

    plt.tight_layout()
    return fig