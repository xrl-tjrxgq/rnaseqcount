# __init__.py
"""
RNA-seq数据分析工具包

提供RNA-seq数据读取、标准化、差异表达分析和可视化功能。
"""

# 导入核心类
from .core import CountData

# 导入IO功能
from .io import writein, writeout

# 导入统计分析功能
from .stats import (
    calculate_CPM,
    calculate_FPK,
    calculate_FPKM,
    calculate_TPM,
    filter_low_expression,
    pca,
    sample_correlation,
    group_samples
)

# 导入差异表达算法
from .algorithms import (
    differential_test,
    pick_significance
)

# 导入可视化功能
from .visualize import (
    plot_pca,
    plot_sample_correlation,
    plot_differential_test
)

# 定义包的版本
__version__ = "1.0.0"

# 定义公开的API接口
__all__ = [
    # 核心类
    "CountData",
    
    # IO功能
    "writein",
    "writeout",
    
    # 统计分析
    "calculate_CPM",
    "calculate_FPK",
    "calculate_FPKM",
    "calculate_TPM",
    "filter_low_expression",
    "pca",
    "sample_correlation",
    "group_samples",
    
    # 差异表达算法
    "differential_test",
    "pick_significance",
    
    # 可视化
    "plot_pca",
    "plot_sample_correlation",
    "plot_differential_test",
]