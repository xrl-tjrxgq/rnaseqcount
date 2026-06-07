# data_simulator.py
import numpy as np
import pandas as pd

def generate_test_data_fixed(
    n_genes: int = 1000,
    n_samples_per_group: int = 5,
    n_groups: int = 2,
    n_batches: int = 2,
    de_ratio: float = 0.3,  # 差异基因比例
    seed: int = 42
) -> tuple:
    """
    修正版：支持多组之间都有差异
    """
    np.random.seed(seed)
    
    # 1. 生成基因名称和长度
    gene_names = [f"gene_{i+1}" for i in range(n_genes)]
    gene_lengths = np.random.lognormal(mean=8, sigma=0.8, size=n_genes)
    gene_lengths = np.round(gene_lengths).astype(int)
    gene_lengths = np.clip(gene_lengths, 500, 50000)
    length_dict = {gene: length for gene, length in zip(gene_names, gene_lengths)}
    
    # 2. 生成样本元数据
    group_names = [f"group_{i+1}" for i in range(n_groups)]
    batch_names = [f"batch_{i+1}" for i in range(n_batches)]
    
    sample_names = []
    group_assignments = []
    batch_assignments = []
    sex_assignments = []
    
    for g_idx, group in enumerate(group_names):
        for s_idx in range(n_samples_per_group):
            sample_name = f"{group}_sample_{s_idx+1}"
            sample_names.append(sample_name)
            group_assignments.append(group)
            batch = np.random.choice(batch_names)
            batch_assignments.append(batch)
            sex = np.random.choice(['M', 'F'])
            sex_assignments.append(sex)
    
    metadata = pd.DataFrame({
        'condition': group_assignments,
        'batch': batch_assignments,
        'sex': sex_assignments
    }, index=sample_names)
    
    # 3. 生成计数矩阵
    counts_matrix = np.zeros((n_genes, len(sample_names)))
    
    # 基础表达水平
    base_expr = np.random.lognormal(mean=5, sigma=1.5, size=n_genes)
    base_expr = base_expr / base_expr.mean() * 1000
    
    # 为每个基因随机选择一个"基准组"（作为参照）
    # 其他组相对于基准组可能有差异
    for g_idx in range(n_genes):
        # 随机决定这个基因是否差异表达
        is_de = np.random.random() < de_ratio
        
        for s_idx, sample in enumerate(sample_names):
            group = group_assignments[s_idx]
            batch = batch_assignments[s_idx]
            
            # 基础表达量
            expr_mean = base_expr[g_idx]
            
            # 如果这个基因是差异表达基因，且不是基准组
            if is_de:
                # 随机选择一个组作为"处理组"（可以是任意组）
                # 这里简化：让 group_2 和 group_3 相对于 group_1 有差异
                if group == group_names[1]:  # group_2
                    fold_change = np.random.choice([-1, 1]) * np.random.uniform(1, 3.3)
                    expr_mean = expr_mean * (2 ** fold_change)
                elif n_groups > 2 and group == group_names[2]:  # group_3
                    fold_change = np.random.choice([-1, 1]) * np.random.uniform(1, 2.5)
                    expr_mean = expr_mean * (2 ** fold_change)
            
            # 批次效应
            if batch == batch_names[1]:
                expr_mean = expr_mean * np.random.uniform(0.8, 1.2)
            
            # 生成计数
            dispersion = 0.2
            if expr_mean > 0:
                overdispersion = np.random.gamma(shape=max(0.1, expr_mean/dispersion), 
                                                scale=dispersion)
                count = np.random.poisson(overdispersion)
            else:
                count = 0
            
            counts_matrix[g_idx, s_idx] = max(0, int(count))
    
    # 转换为DataFrame
    counts = pd.DataFrame(
        counts_matrix,
        index=gene_names,
        columns=sample_names
    )
    counts = counts.astype(int)
    
    # 添加低表达基因
    low_expr_genes = np.random.choice(n_genes, size=int(n_genes*0.1), replace=False)
    for g_idx in low_expr_genes:
        zero_samples = np.random.choice(len(sample_names), 
                                       size=int(len(sample_names)*0.8), 
                                       replace=False)
        for s_idx in zero_samples:
            counts.iloc[g_idx, s_idx] = np.random.poisson(2)
    
    return counts, metadata, length_dict