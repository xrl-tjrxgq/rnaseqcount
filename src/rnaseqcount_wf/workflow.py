# workflow.py
import os
import pandas as pd
import numpy as np
import tempfile
from pathlib import Path

# 导入rnaseqcount库的函数
from rnaseqcount.core import CountData
from rnaseqcount.io import writein, writeout
from rnaseqcount.stats import calculate_TPM, filter_low_expression
from rnaseqcount.algorithms import differential_test, pick_significance
from rnaseqcount.visualize import plot_pca, plot_sample_correlation, plot_differential_test

# 导入数据模拟函数
from .data_simulator import generate_test_data_fixed


def run_workflow(
    data_source: str = "real",
    # 真实数据参数
    input_dir: str = None,           # 包含 counts.csv, metadata.csv, length.csv 的目录
    # 模拟数据参数
    n_genes: int = 1000,
    n_samples_per_group: int = 5,
    n_groups: int = 2,
    n_batches: int = 2,
    de_ratio: float = 0.2,
    seed: int = 42,
    # 分析参数
    method: str = "TPM",
    filter_threshold: float = 1.0,
    grouped_by: str = "condition",
    test_method: str = "t_test",
    paired: bool = False,
    fdr_correction: bool = True,
    significance_threshold: float = 0.05,
    log2fc_threshold: float = 1.0,
    # 输出参数
    output_dir: str = "./rnaseq_results",
    prefix: str = "workflow"
):
    """
    一键运行 RNA-seq 分析流程（模拟数据或真实数据）

    真实数据模式要求 input_dir 中包含：
        - counts.csv   (第一列 'id'，其余列为样本)
        - metadata.csv (至少包含 'sample_id' 和 'condition' 两列)
        - length.csv   (可选，包含 'id' 和 'length' 两列)
    """
    # 创建输出目录
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    # ---------- 1. 准备数据 ----------
    if data_source == "simulated":
        print(">>> 使用模拟数据 <<<")
        counts, metadata, length_dict = generate_test_data_fixed(
            n_genes=n_genes,
            n_samples_per_group=n_samples_per_group,
            n_groups=n_groups,
            n_batches=n_batches,
            de_ratio=de_ratio,
            seed=seed
        )
        cdata = CountData(counts=counts, metadata=metadata, length=length_dict)
        
    elif data_source == "real":
        print(">>> 使用真实数据 <<<")
        if input_dir is None:
            raise ValueError("真实数据模式必须提供 input_dir 参数")
        input_path = Path(input_dir)
        counts_path = input_path / "counts.csv"
        metadata_path = input_path / "metadata.csv"
        length_path = input_path / "length.csv" if (input_path / "length.csv").exists() else None
        
        if not counts_path.exists():
            raise FileNotFoundError(f"找不到 counts.csv 文件：{counts_path}")
        if not metadata_path.exists():
            raise FileNotFoundError(f"找不到 metadata.csv 文件：{metadata_path}")
        
        cdata = writein(
            counts_path=str(counts_path),
            metadata_path=str(metadata_path),
            length_path=str(length_path) if length_path else None
        )
    else:
        raise ValueError("data_source 必须是 'simulated' 或 'real'")

    print(f"数据加载完成: {cdata.n_genes} 基因, {cdata.n_samples} 样本")

    # 如果用户请求 TPM/FPKM 但没有长度文件，给出提示并降级为 CPM
    if method in ["TPM", "FPKM"] and cdata.length is None:
        print(f"警告：未提供基因长度文件，无法使用 {method} 标准化，将使用 CPM 代替。")
        method = "CPM"

    # ---------- 2. 标准化与过滤 ----------
    print(f"\n>>> 执行 {method} 标准化，过滤阈值 = {filter_threshold}")
    filtered_expr = filter_low_expression(cdata, method=method, filter_threshold=filter_threshold)
    print(f"过滤后保留 {filtered_expr.shape[0]} 个基因")

    # ---------- 3. PCA 分析 ----------
    print("\n>>> PCA 分析")
    fig_pca = plot_pca(
        cdata, method=method, filter_threshold=filter_threshold,
        grouped_by=grouped_by, label_name=None, show_labels=True
    )
    pca_path = os.path.join(output_dir, f"{prefix}_pca.png")
    writeout(fig_pca, pca_path.replace(".png", ""))

    # ---------- 4. 样本相关性热图 ----------
    print("\n>>> 样本相关性")
    fig_corr = plot_sample_correlation(
        cdata, method=method, filter_threshold=filter_threshold,
        label_name=grouped_by
    )
    corr_path = os.path.join(output_dir, f"{prefix}_correlation")
    writeout(fig_corr, corr_path)

    # ---------- 5. 差异表达分析 ----------
    print(f"\n>>> 差异表达分析 (分组列: {grouped_by}, 检验: {test_method})")
    diff_res = differential_test(
        cdata,
        method=method,
        filter_threshold=filter_threshold,
        grouped_by=grouped_by,
        log2_transform=True,
        fdr_correction=fdr_correction,
        test_method=test_method,
        paired=paired
    )
    diff_path = os.path.join(output_dir, f"{prefix}_diff_results.csv")
    writeout(diff_res, diff_path.replace(".csv", ""))

    # 显著差异基因筛选
    sig_genes = pick_significance(
        cdata,
        method=method,
        filter_threshold=filter_threshold,
        grouped_by=grouped_by,
        fdr_correction=fdr_correction,
        significance_threshold=significance_threshold
    )
    sig_path = os.path.join(output_dir, f"{prefix}_significant_genes.csv")
    writeout(sig_genes, sig_path.replace(".csv", ""))
    print(f"显著差异基因数: {len(sig_genes)}")

    # ---------- 6. 火山图 ----------
    print("\n>>> 绘制火山图")
    fig_volcano = plot_differential_test(
        cdata,
        method=method,
        filter_threshold=filter_threshold,
        grouped_by=grouped_by,
        fdr_correction=fdr_correction,
        test_method=test_method,
        significance_threshold=significance_threshold,
        log2fc_threshold=log2fc_threshold
    )
    volcano_path = os.path.join(output_dir, f"{prefix}_volcano.png")
    writeout(fig_volcano, volcano_path.replace(".png", ""))

    print(f"\n所有结果已保存至: {output_dir}")
    return cdata, diff_res, sig_genes