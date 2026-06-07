import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import argparse

def plot_volcano_from_deseq2(
    csv_path: str,
    output_path: str,
    log2fc_threshold: float = 1.0,
    pvalue_threshold: float = 0.05,
    group1_label: str = "medium supplemented with filter-sterilized Enterococcus faecium supernatant (equivalent to 1x10^8 CFU/mL)",
    group2_label: str = "standard complete medium"
):
    """
    读取 DESeq2 输出的 CSV/TSV 文件，绘制火山图。
    文件应包含三列：gene_id, log2FoldChange, pvalue（制表符或逗号分隔）
    """
    # 自动检测分隔符
    with open(csv_path, 'r', encoding='utf-8') as f:
        first_line = f.readline()
        if '\t' in first_line:
            sep = '\t'
        else:
            sep = ','
    
    # 读取数据
    df = pd.read_csv(csv_path, encoding='utf-8', sep=sep)
    
    # 提取列（假设列名包含这些关键词，不区分大小写）
    col_mapping = {}
    for col in df.columns:
        col_lower = col.lower()
        if 'log2fold' in col_lower or 'log2fc' in col_lower:
            col_mapping['log2fc'] = col
        elif 'pvalue' in col_lower or 'p_val' in col_lower or 'pvalue' == col_lower:
            col_mapping['pval'] = col
        elif 'gene' in col_lower or 'id' in col_lower:
            col_mapping['gene'] = col
    
    # 如果没有自动匹配，则假定前三列顺序为 gene, log2FoldChange, pvalue
    if 'gene' not in col_mapping:
        gene_col = df.columns[0]
    else:
        gene_col = col_mapping['gene']
    if 'log2fc' not in col_mapping:
        log2fc_col = df.columns[1]
    else:
        log2fc_col = col_mapping['log2fc']
    if 'pval' not in col_mapping:
        pval_col = df.columns[2]
    else:
        pval_col = col_mapping['pval']
    
    # 去掉缺失值
    df = df[[gene_col, log2fc_col, pval_col]].dropna()
    log2fc = df[log2fc_col].values
    pvals = df[pval_col].values
    
    # 计算 -log10(pvalue)
    neg_log10_p = -np.log10(pvals)
    
    # 标记显著性
    significant = pvals < pvalue_threshold
    up_regulated = significant & (log2fc > log2fc_threshold)
    down_regulated = significant & (log2fc < -log2fc_threshold)
    not_significant = ~significant | ((log2fc >= -log2fc_threshold) & (log2fc <= log2fc_threshold))
    
    # 绘图
    fig, ax = plt.subplots(figsize=(10, 8))
    
    ax.scatter(log2fc[not_significant], neg_log10_p[not_significant],
               alpha=0.5, s=10, c="grey", label='Not significant')
    ax.scatter(log2fc[up_regulated], neg_log10_p[up_regulated],
               alpha=0.7, s=15, c="red", label=f'Up-regulated (log2FC > {log2fc_threshold})')
    ax.scatter(log2fc[down_regulated], neg_log10_p[down_regulated],
               alpha=0.7, s=15, c="blue", label=f'Down-regulated (log2FC < -{log2fc_threshold})')
    
    # 添加阈值线
    ax.axhline(y=-np.log10(pvalue_threshold), color='red', linestyle='--', alpha=0.5)
    ax.axvline(x=log2fc_threshold, color='red', linestyle='--', alpha=0.5)
    ax.axvline(x=-log2fc_threshold, color='red', linestyle='--', alpha=0.5)
    
    ax.set_xlabel('log2(Fold Change)', fontsize=12, fontweight='bold')
    ax.set_ylabel('-log10(p-value)', fontsize=12, fontweight='bold')
    
    title = f'Volcano Plot\n{group1_label} vs {group2_label}\n(p-value threshold: {pvalue_threshold})'
    ax.set_title(title, fontsize=14, fontweight='bold')
    ax.legend(loc='best', frameon=True, fancybox=True, shadow=True)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"火山图已保存至: {output_path}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="从 DESeq2 结果绘制火山图")
    parser.add_argument("--csv", required=True, help="DESeq2 输出的 CSV/TSV 文件路径")
    parser.add_argument("--output", default="volcano.png", help="输出图片路径（默认 volcano.png）")
    parser.add_argument("--log2fc_threshold", type=float, default=1.0, help="log2FC 阈值（默认 1.0）")
    parser.add_argument("--p_threshold", type=float, default=0.05, help="p 值阈值（默认 0.05）")
    args = parser.parse_args()
    
    plot_volcano_from_deseq2(
        csv_path=args.csv,
        output_path=args.output,
        log2fc_threshold=args.log2fc_threshold,
        pvalue_threshold=args.p_threshold
    )