# core.py
import pandas as pd
import numpy as np


class CountData:
    """
    RNA-seq 计数数据的核心容器类
    存储counts矩阵,metadata和基因长度,提供基本的数据操作接口
    
    counts:pd.DataFrame
        index:基因名/特征ID
        columns:样本名
    metadata:pd.DataFrame
        index:样本名
        columns:样本属性
    length:dict
        key:基因名/特征ID
        value:基因长度
    """
    
    def __init__(
        self,
        counts: pd.DataFrame=None,
        metadata: pd.DataFrame=None,
        length: dict=None
    ):
        ## 构造函数，初始化CountData对象
        ## counts: 原始表达计数矩阵，基因×样本
        ## metadata: 样本元数据，包含样本属性信息
        ## length: 基因长度字典，用于FPKM/TPM计算
        self._validate(counts,metadata,length)
        self._counts=counts.copy() if counts is not None else None
        self._metadata=metadata.copy() if metadata is not None else None
        self._length=length.copy() if length is not None else None

    def _validate(self,counts,metadata,length):
        ## 内部验证函数，检查输入数据的有效性
        ## counts: 检查缺失值、负值、非整数值，以及是否与length/metadata索引匹配
        ## metadata: 检查缺失值
        if counts is not None:
            if counts.isnull().any().any():
                raise ValueError("counts 不能包含缺失值")
            if (counts<0).any().any():
                raise ValueError("counts 不能包含负值")
            if (counts % 1 != 0).any().any():
                raise ValueError("counts 不能包含非整数值")
            if length is not None:
                if not all(gene in counts.index for gene in length.keys()):
                    raise ValueError("counts 索引 与 length 键值不匹配")
            if metadata is not None:
                if not all(sample in metadata.index for sample in counts.columns):
                    raise ValueError("metadata 索引 与 counts 列名不匹配")
        if metadata is not None:
            if metadata.isnull().any().any():
                raise ValueError("metadata 不能包含缺失值")

    @property
    def counts(self) -> pd.DataFrame:
        ## 属性：返回counts矩阵的副本
        if self._counts is None:
            return None
        return self._counts.copy()
    
    @property
    def metadata(self) -> pd.DataFrame:
        ## 属性：返回metadata的副本
        if self._metadata is None:
            return None
        return self._metadata.copy()
    
    @property
    def length(self) -> dict:
        ## 属性：返回基因长度字典的副本
        if self._length is None:
            return None
        return self._length.copy()
    
    @property
    def shape(self) -> tuple:
        ## 属性：返回counts矩阵的维度 (基因数, 样本数)
        return self._counts.shape
    
    @property
    def n_genes(self) -> int:
        ## 属性：返回基因数量
        return self._counts.shape[0]
    
    @property
    def n_samples(self) -> int:
        ## 属性：返回样本数量
        return self._counts.shape[1]
    
    @property
    def gene_names(self) -> list:
        ## 属性：返回所有基因名称列表
        return self._counts.index.tolist()
    
    @property
    def sample_names(self) -> list:
        ## 属性：返回所有样本名称列表
        return self._counts.columns.tolist()

    @property
    def sample_reads_sum(self) -> dict:
        ## 属性：计算每个样本的总reads数
        ## 返回: {样本名: 总reads数}
        sum_dict={}
        for sample in self.sample_names:
            sum_dict[sample]=self._counts[sample].sum()
        return sum_dict

    @staticmethod
    def significance_symbol(values:pd.Series) -> list:
        ## 静态方法：根据p值/FDR值返回显著性标记符号
        ## values: p值或FDR值序列
        ## 返回: 符号列表 - 'ns'(p>0.05), '*'(p≤0.05), '**'(p≤0.01), '***'(p≤0.001), '****'(p≤0.0001)
        significance=[]
        for value in values:
            if value>0.05:
                significance.append('ns')
            elif value>0.01:
                significance.append('*')
            elif value>0.001:
                significance.append('**')
            elif value>0.0001:
                significance.append('***')
            else:
                significance.append('****')
        return significance