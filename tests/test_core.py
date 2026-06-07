# tests/test_core.py
"""测试 CountData 核心类"""
import pytest
import pandas as pd
import numpy as np

from rnaseqcount import CountData


class TestCountData:
    """测试 CountData 类的核心功能"""
    
    def setup_method(self):
        """准备标准测试数据"""
        self.counts = pd.DataFrame({
            "S1": [100, 200, 300],
            "S2": [150, 250, 350]
        }, index=["Gene1", "Gene2", "Gene3"])
        
        self.metadata = pd.DataFrame({
            "condition": ["control", "treatment"]
        }, index=["S1", "S2"])
        
        self.length = {
            "Gene1": 2000,
            "Gene2": 1500,
            "Gene3": 3000
        }
    
    # ==================== 正常初始化测试 ====================
    
    def test_init_all_params(self):
        """测试使用所有参数正常初始化"""
        cdata = CountData(counts=self.counts, metadata=self.metadata, length=self.length)
        
        assert cdata.counts is not None
        assert cdata.metadata is not None
        assert cdata.length is not None
    
    def test_init_counts_only(self):
        """测试只传入 counts"""
        cdata = CountData(counts=self.counts)
        assert cdata.counts is not None
        assert cdata.metadata is None
        assert cdata.length is None
    
    def test_init_metadata_only(self):
        """测试只传入 metadata"""
        cdata = CountData(metadata=self.metadata)
        assert cdata.counts is None
        assert cdata.metadata is not None
        assert cdata.length is None
    
    def test_init_length_only(self):
        """测试只传入 length"""
        cdata = CountData(length=self.length)
        assert cdata.counts is None
        assert cdata.metadata is None
        assert cdata.length is not None
    
    def test_init_no_params(self):
        """测试不传入任何参数"""
        cdata = CountData()
        assert cdata.counts is None
        assert cdata.metadata is None
        assert cdata.length is None
    
    # ==================== _validate 错误测试 ====================
    
    def test_counts_with_nan(self):
        """测试 counts 包含缺失值时抛出 ValueError"""
        bad_counts = self.counts.copy()
        bad_counts.iloc[0, 0] = np.nan
        
        with pytest.raises(ValueError, match="counts 不能包含缺失值"):
            CountData(counts=bad_counts)
    
    def test_counts_with_negative(self):
        """测试 counts 包含负值时抛出 ValueError"""
        bad_counts = self.counts.copy()
        bad_counts.iloc[0, 0] = -1
        
        with pytest.raises(ValueError, match="counts 不能包含负值"):
            CountData(counts=bad_counts)
    
    def test_counts_with_non_integer(self):
        """测试 counts 包含非整数值时抛出 ValueError"""
        bad_counts = self.counts.copy().astype(float) 
        bad_counts.iloc[0, 0] = 100.5
        
        with pytest.raises(ValueError, match="counts 不能包含非整数值"):
            CountData(counts=bad_counts)
    
    def test_counts_length_mismatch(self):
        """测试 counts 索引与 length 键值不匹配时抛出 ValueError"""
        bad_length = {
            "Gene1": 2000,
            "GeneX": 1500,  # Gene2 变成了 GeneX
            "Gene3": 3000
        }
        
        with pytest.raises(ValueError, match="counts 索引 与 length 键值不匹配"):
            CountData(counts=self.counts, length=bad_length)
    
    def test_metadata_counts_mismatch(self):
        """测试 metadata 索引与 counts 列名不匹配时抛出 ValueError"""
        bad_metadata = pd.DataFrame({
            "condition": ["control", "treatment"]
        }, index=["S1", "SX"])  # S2 变成了 SX
        
        with pytest.raises(ValueError, match="metadata 索引 与 counts 列名不匹配"):
            CountData(counts=self.counts, metadata=bad_metadata)
    
    def test_metadata_with_nan(self):
        """测试 metadata 包含缺失值时抛出 ValueError"""
        bad_metadata = self.metadata.copy()
        bad_metadata.iloc[0, 0] = np.nan
        
        with pytest.raises(ValueError, match="metadata 不能包含缺失值"):
            CountData(metadata=bad_metadata)
    
    # ==================== 属性测试 ====================
    
    def test_counts_property_returns_copy(self):
        """测试 counts 属性返回副本，修改不影响原对象"""
        cdata = CountData(counts=self.counts)
        counts_copy = cdata.counts
        counts_copy.iloc[0, 0] = 99999
        
        # 再次获取应该还是原来的值
        assert cdata.counts.iloc[0, 0] == 100
    
    def test_metadata_property_returns_copy(self):
        """测试 metadata 属性返回副本"""
        cdata = CountData(counts=self.counts, metadata=self.metadata)
        meta_copy = cdata.metadata
        meta_copy.iloc[0, 0] = "changed"
        
        assert cdata.metadata.iloc[0, 0] == "control"
    
    def test_length_property_returns_copy(self):
        """测试 length 属性返回副本"""
        cdata = CountData(length=self.length)
        len_copy = cdata.length
        len_copy["Gene1"] = 99999
        
        assert cdata.length["Gene1"] == 2000
    
    def test_shape_property(self):
        """测试 shape 属性返回正确的维度"""
        cdata = CountData(counts=self.counts)
        assert cdata.shape == (3, 2)
    
    def test_n_genes_property(self):
        """测试 n_genes 属性"""
        cdata = CountData(counts=self.counts)
        assert cdata.n_genes == 3
    
    def test_n_samples_property(self):
        """测试 n_samples 属性"""
        cdata = CountData(counts=self.counts)
        assert cdata.n_samples == 2
    
    def test_gene_names_property(self):
        """测试 gene_names 属性"""
        cdata = CountData(counts=self.counts)
        assert cdata.gene_names == ["Gene1", "Gene2", "Gene3"]
    
    def test_sample_names_property(self):
        """测试 sample_names 属性"""
        cdata = CountData(counts=self.counts)
        assert cdata.sample_names == ["S1", "S2"]
    
    def test_sample_reads_sum(self):
        """测试 sample_reads_sum 属性计算正确"""
        cdata = CountData(counts=self.counts)
        reads_sum = cdata.sample_reads_sum
        
        assert reads_sum["S1"] == 600  # 100 + 200 + 300
        assert reads_sum["S2"] == 750  # 150 + 250 + 350
    
    # ==================== 静态方法测试 ====================
    
    def test_significance_symbol(self):
        """测试 significance_symbol 静态方法"""
        p_values = pd.Series([0.5, 0.05, 0.01, 0.001, 0.0001, 0.00001])
        symbols = CountData.significance_symbol(p_values)
        
        assert symbols == ['ns', '*', '**', '***', '****', '****']
    


if __name__ == "__main__":
    pytest.main([__file__, "-v"])