# algorithms.py
import pandas as pd
import numpy as np
from scipy import stats
from scipy.stats import false_discovery_control
from .core import CountData
from .stats import filter_low_expression, group_samples


def differential_test(countdata: CountData, method: str = "TPM", filter_threshold: float = 0, 
                      grouped_by: str = None, groups: list = None, 
                      log2_transform: bool = True, fdr_correction: bool = True, 
                      test_method: str = "t_test", paired: bool = False):
    ## 执行差异表达分析（两组比较）
    ## method: 归一化方法
    ## filter_threshold: 过滤低表达基因的阈值
    ## grouped_by: metadata列名，自动确定两组
    ## groups: 手动指定两组名称的列表，如["control","treatment"]
    ## log2_transform: 是否对表达值进行log2(x+1)变换
    ## fdr_correction: 是否进行FDR多重检验校正
    ## test_method: 统计检验方法，"t_test"或"wilcoxon"
    ## paired: 是否配对检验
    ## 返回: DataFrame包含统计量、log2FC、p值、FDR校正p值、显著性标记
    if method == "TPM":
        expr = filter_low_expression(countdata, method="TPM", filter_threshold=filter_threshold)
    elif method == "FPKM":
        expr = filter_low_expression(countdata, method="FPKM", filter_threshold=filter_threshold)
    elif method == "CPM":
        expr = filter_low_expression(countdata, method="CPM", filter_threshold=filter_threshold)
    else:
        raise ValueError("归一化方法必须选择CPM/TPM/FPKM")

    if log2_transform:
        expr = np.log2(expr + 1)

    group_dict = group_samples(countdata, grouped_by=grouped_by, groups=groups)
    if len(group_dict.keys()) != 2:
        raise Exception("t-test要求分组数为2")
    for key in group_dict.keys():
        if len(group_dict[key]) < 2:
            raise Exception("每组需要至少两个样本")

    group1_samples = group_dict[list(group_dict.keys())[0]]
    group2_samples = group_dict[list(group_dict.keys())[1]]
    group1_name = list(group_dict.keys())[0]
    group2_name = list(group_dict.keys())[1]

    statistics = []
    p_values = []
    log2_fold_changes = []

    for idx, row in expr.iterrows():
        group1_vals = row[group1_samples].values
        group2_vals = row[group2_samples].values
        if log2_transform:
            lfc = np.mean(group2_vals) - np.mean(group1_vals)
        else:
            lfc = np.log2((np.mean(group2_vals) + 1) / (np.mean(group1_vals) + 1))
        log2_fold_changes.append(lfc)
        if test_method == "t_test":
            if paired == True:
                stat, p_val = stats.ttest_rel(group1_vals, group2_vals, nan_policy='omit')
            else:
                stat, p_val = stats.ttest_ind(group1_vals, group2_vals, nan_policy='omit')
        if test_method == "wilcoxon":
            if paired == True:
                stat, p_val = stats.wilcoxon(group1_vals, group2_vals, nan_policy='omit')
            else:
                stat, p_val = stats.mannwhitneyu(group1_vals, group2_vals, alternative='two-sided', nan_policy='omit')
        statistics.append(stat)
        p_values.append(p_val)

    result_df = pd.DataFrame({f"{test_method}": statistics, "log2_fold_changes": log2_fold_changes, 
                              "p_values": p_values}, index=expr.index)
    if fdr_correction:
        result_df["fdr"] = false_discovery_control(result_df["p_values"], method='bh')
        result_df['significance'] = countdata.significance_symbol(result_df["fdr"])
    else:
        result_df['significance'] = countdata.significance_symbol(result_df["p_values"])
    return result_df


def pick_significance(countdata: CountData, method: str = "TPM", filter_threshold: float = 0, 
                      grouped_by: str = None, groups: list = None, 
                      log2_transform: bool = True, fdr_correction: bool = True, 
                      test_method: str = "t_test", paired: bool = False, 
                      significance_threshold: float = 0.05):
    ## 筛选显著差异表达基因
    ## 参数同differential_test，额外参数:
    ## significance_threshold: 显著性阈值（p值或FDR）
    ## log2fc_threshold: log2 fold change阈值（本函数中定义但未使用该阈值进行筛选，仅筛选p值）
    ## 返回: 仅包含显著差异表达基因的DataFrame
    result_df = differential_test(countdata, method=method, filter_threshold=filter_threshold, 
                                   grouped_by=grouped_by, groups=groups, 
                                   log2_transform=log2_transform, fdr_correction=fdr_correction, 
                                   test_method=test_method, paired=paired)

    if fdr_correction:
        p_col = "fdr"
    else:
        p_col = "p_values"
    significance_df = result_df[result_df[p_col] < significance_threshold]
    return significance_df