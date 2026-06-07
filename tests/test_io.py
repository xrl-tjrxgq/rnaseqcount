# tests/test_io.py
"""测试文件读写功能"""
import pytest
import pandas as pd
import tempfile
import os
import matplotlib
matplotlib.use('Agg')  
import matplotlib.pyplot as plt

from rnaseqcount import writein, writeout, CountData


class TestIO:
    """测试IO模块"""
    
    def setup_method(self):
        """准备测试数据"""
        # 创建临时目录
        self.temp_dir = tempfile.mkdtemp()
        
        # 创建标准测试CSV文件
        self.counts_path = os.path.join(self.temp_dir, "counts.csv")
        self.metadata_path = os.path.join(self.temp_dir, "metadata.csv")
        self.length_path = os.path.join(self.temp_dir, "lengths.csv")
        
        # counts文件：标准格式，第一列为 "id"
        pd.DataFrame({
            "id": ["Gene1", "Gene2", "Gene3"],
            "S1": [100, 200, 300],
            "S2": [150, 250, 350],
            "S3": [120, 220, 320]
        }).to_csv(self.counts_path, index=False)
        
        # metadata文件：标准格式，第一列为 "run"
        pd.DataFrame({
            "run": ["S1", "S2", "S3"],
            "condition": ["control", "treatment", "treatment"],
            "batch": [1, 1, 2]
        }).to_csv(self.metadata_path, index=False)
        
        # length文件：标准格式，第一列 "id"，第二列 "length"
        pd.DataFrame({
            "id": ["Gene1", "Gene2", "Gene3"],
            "length": [2000, 1500, 3000]
        }).to_csv(self.length_path, index=False)
    
    def teardown_method(self):
        """清理临时文件"""
        import shutil
        shutil.rmtree(self.temp_dir)
    
    # ==================== writein() 正常读取测试 ====================
    
    def test_read_all_files(self):
        """测试同时读取所有三个文件"""
        cdata = writein(
            counts_path=self.counts_path,
            metadata_path=self.metadata_path,
            length_path=self.length_path
        )
        
        # 验证返回类型
        assert isinstance(cdata, CountData)
        
        # 验证 counts 数据
        assert cdata.counts is not None
        assert cdata.counts.shape == (3, 3)  # 3个基因，3个样本
        assert list(cdata.counts.index) == ["Gene1", "Gene2", "Gene3"]
        assert list(cdata.counts.columns) == ["S1", "S2", "S3"]
        assert cdata.counts.loc["Gene1", "S1"] == 100
        
        # 验证 metadata
        assert cdata.metadata is not None
        assert cdata.metadata.shape == (3, 2)  # 3个样本，2个属性
        assert list(cdata.metadata.index) == ["S1", "S2", "S3"]
        assert cdata.metadata.loc["S2", "condition"] == "treatment"
        
        # 验证 length
        assert cdata.length is not None
        assert isinstance(cdata.length, dict)
        assert cdata.length["Gene1"] == 2000
        assert cdata.length["Gene2"] == 1500
        assert len(cdata.length) == 3
    
    def test_read_counts_only(self):
        """测试只读取 counts 文件"""
        cdata = writein(counts_path=self.counts_path)
        
        assert isinstance(cdata, CountData)
        assert cdata.counts is not None
        assert cdata.counts.shape == (3, 3)
        assert cdata.metadata is None
        assert cdata.length is None
    
    def test_read_metadata_only(self):
        """测试只读取 metadata 文件"""
        cdata = writein(metadata_path=self.metadata_path)
        
        assert isinstance(cdata, CountData)
        assert cdata.counts is None
        assert cdata.metadata is not None
        assert cdata.metadata.shape == (3, 2)
        assert cdata.length is None
    
    def test_read_length_only(self):
        """测试只读取 length 文件"""
        cdata = writein(length_path=self.length_path)
        
        assert isinstance(cdata, CountData)
        assert cdata.counts is None
        assert cdata.metadata is None
        assert cdata.length is not None
        assert isinstance(cdata.length, dict)
        assert cdata.length["Gene3"] == 3000
    
    def test_read_counts_and_metadata(self):
        """测试读取 counts + metadata，不读取 length"""
        cdata = writein(
            counts_path=self.counts_path,
            metadata_path=self.metadata_path
        )
        
        assert cdata.counts is not None
        assert cdata.metadata is not None
        assert cdata.length is None
    
    def test_read_counts_and_length(self):
        """测试读取 counts + length，不读取 metadata"""
        cdata = writein(
            counts_path=self.counts_path,
            length_path=self.length_path
        )
        
        assert cdata.counts is not None
        assert cdata.metadata is None
        assert cdata.length is not None
    
    # ==================== writein() 边界情况测试 ====================
    
    def test_counts_id_case_insensitive(self):
        """测试 counts 文件第一列列名大小写不敏感（ID, Id, id 都应通过）"""
        # 测试 "ID"
        path_upper = os.path.join(self.temp_dir, "counts_ID.csv")
        pd.DataFrame({
            "ID": ["Gene1", "Gene2"],
            "S1": [100, 200]
        }).to_csv(path_upper, index=False)
        cdata = writein(counts_path=path_upper)
        assert cdata.counts is not None
        assert list(cdata.counts.index) == ["Gene1", "Gene2"]
        
        # 测试 "Id"
        path_mixed = os.path.join(self.temp_dir, "counts_Id.csv")
        pd.DataFrame({
            "Id": ["Gene1", "Gene2"],
            "S1": [100, 200]
        }).to_csv(path_mixed, index=False)
        cdata = writein(counts_path=path_mixed)
        assert cdata.counts is not None
    
    def test_length_column_case_insensitive(self):
        """测试 length 文件列名大小写不敏感（内部会转为小写）"""
        path_mixed = os.path.join(self.temp_dir, "lengths_mixed.csv")
        pd.DataFrame({
            "ID": ["Gene1", "Gene2"],
            "LENGTH": [2000, 1500]
        }).to_csv(path_mixed, index=False)
        
        cdata = writein(length_path=path_mixed)
        assert cdata.length is not None
        assert cdata.length["Gene1"] == 2000
    
    def test_no_files_provided(self):
        """测试所有参数都为 None"""
        cdata = writein()
        
        assert isinstance(cdata, CountData)
        assert cdata.counts is None
        assert cdata.metadata is None
        assert cdata.length is None
    
    def test_empty_dataframe(self):
        """测试读取空数据文件（只有列名，没有数据行）"""
        empty_counts = os.path.join(self.temp_dir, "empty_counts.csv")
        pd.DataFrame(columns=["id", "S1", "S2"]).to_csv(empty_counts, index=False)
        
        cdata = writein(counts_path=empty_counts)
        assert cdata.counts is not None
        assert cdata.counts.shape == (0, 2)  # 0行，2列（S1, S2）
    
    # ==================== writein() 异常测试 ====================
    
    def test_counts_not_csv_extension(self):
        """测试 counts 文件非 .csv 扩展名时抛出异常"""
        txt_path = os.path.join(self.temp_dir, "counts.txt")
        pd.DataFrame({
            "id": ["Gene1"],
            "S1": [100]
        }).to_csv(txt_path, index=False)
        
        with pytest.raises(Exception, match="仅支持.csv格式"):
            writein(counts_path=txt_path)
    
    def test_metadata_not_csv_extension(self):
        """测试 metadata 文件非 .csv 扩展名时抛出异常"""
        txt_path = os.path.join(self.temp_dir, "metadata.txt")
        pd.DataFrame({
            "run": ["S1"],
            "condition": ["ctrl"]
        }).to_csv(txt_path, index=False)
        
        with pytest.raises(Exception, match="仅支持.csv格式"):
            writein(metadata_path=txt_path)
    
    def test_length_not_csv_extension(self):
        """测试 length 文件非 .csv 扩展名时抛出异常"""
        txt_path = os.path.join(self.temp_dir, "lengths.txt")
        pd.DataFrame({
            "id": ["Gene1"],
            "length": [2000]
        }).to_csv(txt_path, index=False)
        
        with pytest.raises(Exception, match="仅支持.csv格式"):
            writein(length_path=txt_path)
    
    def test_counts_wrong_first_column(self):
        """测试 counts 文件第一列不是 id 时抛出 ValueError"""
        bad_path = os.path.join(self.temp_dir, "bad_counts.csv")
        pd.DataFrame({
            "gene_id": ["Gene1", "Gene2"],  # 错误的列名
            "S1": [100, 200]
        }).to_csv(bad_path, index=False)
        
        with pytest.raises(ValueError, match="第一列列名不是ID"):
            writein(counts_path=bad_path)
    
    def test_metadata_wrong_first_column(self):
        """测试 metadata 文件第一列不是 run 时抛出 ValueError"""
        bad_path = os.path.join(self.temp_dir, "bad_metadata.csv")
        pd.DataFrame({
            "sample": ["S1", "S2"],  # 错误的列名
            "condition": ["ctrl", "treat"]
        }).to_csv(bad_path, index=False)
        
        with pytest.raises(ValueError, match="第一列列名不是run"):
            writein(metadata_path=bad_path)
    
    def test_length_wrong_first_column(self):
        """测试 length 文件第一列不是 id 时抛出 ValueError"""
        bad_path = os.path.join(self.temp_dir, "bad_length1.csv")
        pd.DataFrame({
            "gene": ["Gene1", "Gene2"],  # 错误的列名
            "length": [2000, 1500]
        }).to_csv(bad_path, index=False)
        
        with pytest.raises(ValueError, match="第一列列名不是ID"):
            writein(length_path=bad_path)
    
    def test_length_wrong_second_column(self):
        """测试 length 文件第二列不是 length 时抛出 ValueError"""
        bad_path = os.path.join(self.temp_dir, "bad_length2.csv")
        pd.DataFrame({
            "id": ["Gene1", "Gene2"],
            "len": [2000, 1500]  # 错误的列名
        }).to_csv(bad_path, index=False)
        
        with pytest.raises(ValueError, match="第二列列名不是length"):
            writein(length_path=bad_path)
    
    # ==================== writeout() 正常测试 ====================
    
    def test_write_dataframe(self):
        """测试保存 DataFrame 为 CSV"""
        df = pd.DataFrame({
            "A": [1, 2, 3], 
            "B": ["x", "y", "z"]
        })
        output_path = os.path.join(self.temp_dir, "test_output")
        writeout(df, output_path)
        
        output_file = output_path + ".csv"
        assert os.path.exists(output_file)
        
        # 验证文件内容
        read_back = pd.read_csv(output_file)
        pd.testing.assert_frame_equal(read_back, df)
    
    def test_write_dataframe_no_index(self):
        """测试 DataFrame 保存时不包含行索引"""
        df = pd.DataFrame({
            "col1": [10, 20],
            "col2": [30, 40]
        }, index=["row1", "row2"])  # 有自定义索引
        
        output_path = os.path.join(self.temp_dir, "test_no_index")
        writeout(df, output_path)
        
        # 读取回来的文件不应该包含行索引列
        read_back = pd.read_csv(output_path + ".csv")
        assert "row1" not in read_back.columns
        assert list(read_back.columns) == ["col1", "col2"]
        assert read_back["col1"].tolist() == [10, 20]
    
    def test_write_figure(self):
        """测试保存 matplotlib Figure 为 PNG"""
        fig, ax = plt.subplots(figsize=(8, 6))
        ax.plot([1, 2, 3], [1, 4, 9], label="test line")
        ax.set_xlabel("X")
        ax.set_ylabel("Y")
        ax.set_title("Test Plot")
        ax.legend()
        
        output_path = os.path.join(self.temp_dir, "test_plot")
        writeout(fig, output_path)
        
        output_file = output_path + ".png"
        assert os.path.exists(output_file)
        assert os.path.getsize(output_file) > 0  # 文件非空
        
        plt.close(fig)
    



if __name__ == "__main__":
    pytest.main([__file__, "-v"])