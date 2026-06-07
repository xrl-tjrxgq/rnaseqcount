# stats.py
import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
from scipy import stats
from scipy.stats import false_discovery_control
import matplotlib.pyplot as plt
from .core import CountData


def calculate_CPM(countdata: CountData) -> pd.DataFrame:
    ## 计算CPM (Counts Per Million) 归一化
    ## CPM = (counts / 样本总reads数) * 1,000,000
    ## 返回: CPM矩阵，基因×样本
    sample_reads_sum=pd.Series(countdata.sample_reads_sum)
    calculate_CPM=countdata._counts/sample_reads_sum*1000000
    calculate_CPM.index=countdata._counts.index
    return calculate_CPM


def calculate_FPK(countdata: CountData) -> pd.DataFrame:
    ## 计算FPK (Fragments Per Kilobase) - FPKM的中间步骤
    ## FPK = (counts / 基因长度) * 1000
    ## 需要length字典支持
    length_series=pd.Series(countdata._length)
    common_genes=countdata._counts.index.intersection(length_series.index)
    counts_aligned=countdata._counts.loc[common_genes]
    length_aligned=length_series[common_genes]
    calculate_FPK=pd.DataFrame()
    for sample in counts_aligned.columns:
        calculate_FPK[sample]=counts_aligned[sample]/length_aligned*1000
    calculate_FPK.index=counts_aligned.index
    calculate_FPK.index.name=countdata._counts.index.name
    return calculate_FPK


def calculate_FPKM(countdata: CountData) -> pd.DataFrame:
    ## 计算FPKM (Fragments Per Kilobase per Million)
    ## FPKM = FPK / (样本总reads数) * 1,000,000
    ## 返回: FPKM矩阵，基因×样本
    fpk=calculate_FPK(countdata)
    sample_reads_sum=countdata.sample_reads_sum
    calculate_FPKM=pd.DataFrame()
    for sample in fpk.columns:
        calculate_FPKM[sample]=fpk[sample]/sample_reads_sum[sample]*1000000
    calculate_FPKM.index=fpk.index
    calculate_FPKM.index.name=countdata._counts.index.name
    return calculate_FPKM


def calculate_TPM(countdata: CountData) -> pd.DataFrame:
    ## 计算TPM (Transcripts Per Million)
    ## TPM = (FPK / 样本总FPK) * 1,000,000
    ## 需要length字典支持
    ## 返回: TPM矩阵，基因×样本
    if countdata._length==None:
        raise ValueError("计算TPM需要length信息")
    fpk=calculate_FPK(countdata)
    sample_sums=fpk.sum(axis=0)
    calculate_TPM=pd.DataFrame()
    for sample in fpk.columns:
        calculate_TPM[sample]=fpk[sample]/sample_sums[sample]*1000000
    calculate_TPM.index=fpk.index
    calculate_TPM.index.name=countdata._counts.index.name
    return calculate_TPM


def filter_low_expression(countdata: CountData, method: str = "TPM", filter_threshold: float = 0) -> pd.DataFrame:
    ## 基于归一化表达值过滤低表达基因
    ## method: 归一化方法，可选"CPM"/"TPM"/"FPKM"
    ## filter_threshold: 表达阈值，保留至少在一个样本中表达量>该阈值的基因
    ## 返回: 过滤后的表达矩阵
    ## 打印: 保留基因比例
    if method == "CPM":
        to_be_filtered = calculate_CPM(countdata)
    elif method == "TPM":
        to_be_filtered = calculate_TPM(countdata)
    elif method == "FPKM":
        to_be_filtered = calculate_FPKM(countdata)
    else:
        raise ValueError("归一化方法必须选择CPM/TPM/FPKM")
    filtered = to_be_filtered[(to_be_filtered > filter_threshold).any(axis=1)]
    kept_gene_ratio = filtered.shape[0] / to_be_filtered.shape[0]
    print(f"kept_gene_ratio={kept_gene_ratio}")
    return filtered


def pca(countdata: CountData, method: str = "TPM", n_components: int = 2, 
        filter_threshold: float = 0, label_name: str = None):
    ## 执行主成分分析(PCA)
    ## method: 归一化方法，用于过滤低表达基因
    ## n_components: PCA主成分数量
    ## filter_threshold: 低表达过滤阈值
    ## label_name: metadata中的列名，用于将样本标签替换为分组标签
    ## 返回: (pca_df, pca_var_ratio) - PCA结果DataFrame和各主成分解释方差比例
    if method == "TPM":
        expr = filter_low_expression(countdata, method="TPM", filter_threshold=filter_threshold)
    elif method == "FPKM":
        expr = filter_low_expression(countdata, method="FPKM", filter_threshold=filter_threshold)
    elif method == "CPM":
        expr = filter_low_expression(countdata, method="CPM", filter_threshold=filter_threshold)
    else:
        raise ValueError("归一化方法必须选择CPM/TPM/FPKM")
    expr = expr.T
    if countdata._metadata is not None and label_name is not None:
        expr.index = countdata._metadata[label_name]
    pca = PCA(n_components=2)
    pca_result = pca.fit_transform(expr)
    pca_df = pd.DataFrame(pca_result, columns=[f'PC{i+1}' for i in range(n_components)], index=expr.index)
    pca_var_ratio = pca.explained_variance_ratio_
    return pca_df, pca_var_ratio


def sample_correlation(countdata: CountData, method: str = "TPM", 
                       filter_threshold: float = 0, label_name: str = None):
    ## 计算样本间Pearson相关系数矩阵
    ## method: 归一化方法
    ## filter_threshold: 过滤阈值
    ## label_name: 将样本名替换为metadata中的分组标签
    ## 返回: 相关系数矩阵DataFrame
    if method == "TPM":
        expr = filter_low_expression(countdata, method="TPM", filter_threshold=filter_threshold)
    elif method == "FPKM":
        expr = filter_low_expression(countdata, method="FPKM", filter_threshold=filter_threshold)
    elif method == "CPM":
        expr = filter_low_expression(countdata, method="CPM", filter_threshold=filter_threshold)
    else:
        raise ValueError("归一化方法必须选择CPM/TPM/FPKM")
    corr_matrix = expr.corr(method='pearson')
    if countdata._metadata is not None and label_name is not None:
        corr_matrix.columns = countdata._metadata[label_name]
        corr_matrix.index = countdata._metadata[label_name]
    return corr_matrix


def group_samples(countdata: CountData, grouped_by: str = None, groups: list = None) -> dict:
    ## 根据元数据将样本分组
    ## grouped_by: metadata中的列名，按该列的唯一值分组
    ## groups: 手动指定的组名列表，从metadata中筛选包含这些组名的样本
    ## 必须且只能选择grouped_by或groups两种方法其一
    ## 返回: {组名: [样本列表]}
    group_dict = {}
    if grouped_by is not None and groups is None:
        group_names = countdata._metadata[grouped_by].unique()
        for group_name in group_names:
            samples = []
            mask = (countdata._metadata[grouped_by] == group_name)
            samples = countdata._metadata[mask].index.tolist()
            group_dict[group_name] = samples
    elif grouped_by is None and groups is not None:
        group_names = groups
        for group_name in group_names:
            samples = []
            mask = (countdata._metadata == group_name).any(axis=1)
            samples = countdata._metadata[mask].index.tolist()
            group_dict[group_name] = samples
    else:
        raise ValueError("必须指定grouped_by或groups两者其一")
    return group_dict
    