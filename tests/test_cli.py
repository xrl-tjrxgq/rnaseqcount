# tests/test_cli.py
"""测试 CLI 命令行接口"""
import pytest
import sys
import os
import tempfile
import pandas as pd
from unittest.mock import patch, MagicMock
from io import StringIO

# 导入被测模块
from rnaseqcount.cli import main


class TestCLI:
    """测试命令行接口功能"""
    
    def setup_method(self):
        """准备测试数据和临时目录"""
        self.temp_dir = tempfile.mkdtemp()
        
        # 创建测试CSV文件
        self.counts_path = os.path.join(self.temp_dir, "counts.csv")
        self.metadata_path = os.path.join(self.temp_dir, "metadata.csv")
        self.length_path = os.path.join(self.temp_dir, "lengths.csv")
        
        # counts文件
        pd.DataFrame({
            "id": ["Gene1", "Gene2", "Gene3"],
            "S1": [100, 200, 300],
            "S2": [150, 250, 350],
            "S3": [120, 220, 320],
            "S4": [180, 280, 380]
        }).to_csv(self.counts_path, index=False)
        
        # metadata文件（两组：control, treatment）
        pd.DataFrame({
            "run": ["S1", "S2", "S3", "S4"],
            "condition": ["control", "control", "treatment", "treatment"]
        }).to_csv(self.metadata_path, index=False)
        
        # length文件
        pd.DataFrame({
            "id": ["Gene1", "Gene2", "Gene3"],
            "length": [2000, 1500, 3000]
        }).to_csv(self.length_path, index=False)
        
        self.output_prefix = os.path.join(self.temp_dir, "test_output")
    
    def teardown_method(self):
        """清理临时文件"""
        import shutil
        shutil.rmtree(self.temp_dir)
    
    def _run_cli(self, args_list):
        """辅助方法：运行CLI并捕获输出"""
        with patch.object(sys, 'argv', ['rnaseqcount'] + args_list):
            with patch('sys.stdout', new=StringIO()) as fake_stdout:
                with patch('sys.stderr', new=StringIO()) as fake_stderr:
                    try:
                        main()
                        return 0, fake_stdout.getvalue(), fake_stderr.getvalue()
                    except SystemExit as e:
                        return e.code, fake_stdout.getvalue(), fake_stderr.getvalue()
    
    # ==================== 基本功能测试 ====================
    
    def test_help_message(self):
        """测试 --help 显示帮助信息"""
        exit_code, stdout, stderr = self._run_cli(['--help'])
        assert exit_code == 0
        assert 'RNA-seq数据分析工具' in stdout
    
    def test_basic_load_data(self):
        """测试基本数据加载（只传 --counts）"""
        exit_code, stdout, stderr = self._run_cli([
            '--counts', self.counts_path,
            '--output', self.output_prefix
        ])
        assert exit_code == 0
        assert '数据加载成功' in stdout
        assert '3个基因' in stdout
        assert '4个样本' in stdout
    
    def test_load_all_files(self):
        """测试加载所有三个文件"""
        exit_code, stdout, stderr = self._run_cli([
            '--counts', self.counts_path,
            '--metadata', self.metadata_path,
            '--length', self.length_path,
            '--output', self.output_prefix
        ])
        assert exit_code == 0
        assert '数据加载成功' in stdout
    
    # ==================== PCA分析测试 ====================
    
    @patch('rnaseqcount.cli.plot_pca')
    @patch('rnaseqcount.cli.writeout')
    def test_pca_analysis(self, mock_writeout, mock_plot_pca):
        """测试 --pca 执行PCA分析"""
        mock_fig = MagicMock()
        mock_plot_pca.return_value = mock_fig
        
        exit_code, stdout, stderr = self._run_cli([
            '--counts', self.counts_path,
            '--metadata', self.metadata_path,
            '--pca',
            '--method', 'TPM',
            '--filter', '0',
            '--output', self.output_prefix
        ])
        
        assert exit_code == 0
        assert '执行PCA分析' in stdout
        assert 'PCA图已保存' in stdout
        mock_plot_pca.assert_called_once()
        mock_writeout.assert_called_with(mock_fig, f"{self.output_prefix}_pca")
    
    # ==================== 相关性分析测试 ====================
    
    @patch('rnaseqcount.cli.plot_sample_correlation')
    @patch('rnaseqcount.cli.writeout')
    def test_correlation_analysis(self, mock_writeout, mock_plot_corr):
        """测试 --correlation 执行相关性分析"""
        mock_fig = MagicMock()
        mock_plot_corr.return_value = mock_fig
        
        exit_code, stdout, stderr = self._run_cli([
            '--counts', self.counts_path,
            '--metadata', self.metadata_path,
            '--correlation',
            '--method', 'CPM',
            '--output', self.output_prefix
        ])
        
        assert exit_code == 0
        assert '计算样本相关性' in stdout
        assert '相关性热图已保存' in stdout
        mock_plot_corr.assert_called_once()
        mock_writeout.assert_called_with(mock_fig, f"{self.output_prefix}_correlation")
    
    # ==================== 差异表达分析测试 ====================
    
    @patch('rnaseqcount.cli.plot_differential_test')
    @patch('rnaseqcount.cli.differential_test')
    @patch('rnaseqcount.cli.writeout')
    def test_diff_analysis_with_grouped_by(self, mock_writeout, mock_diff, mock_plot_diff):
        """测试 --diff 使用 --grouped_by 分组"""
        mock_fig = MagicMock()
        mock_result = MagicMock()
        mock_plot_diff.return_value = mock_fig
        mock_diff.return_value = mock_result
        
        exit_code, stdout, stderr = self._run_cli([
            '--counts', self.counts_path,
            '--metadata', self.metadata_path,
            '--diff',
            '--grouped_by', 'condition',
            '--test_method', 't_test',
            '--output', self.output_prefix
        ])
        
        assert exit_code == 0
        assert '执行差异表达分析' in stdout
        assert '差异表达分析完成' in stdout
        
        # 验证调用参数
        mock_plot_diff.assert_called_once()
        call_kwargs = mock_plot_diff.call_args[1]
        assert call_kwargs['grouped_by'] == 'condition'
        assert call_kwargs['test_method'] == 't_test'
        assert call_kwargs['paired'] == False
        
        mock_diff.assert_called_once()
        mock_writeout.assert_any_call(mock_fig, f"{self.output_prefix}_volcano")
        mock_writeout.assert_any_call(mock_result, f"{self.output_prefix}_diff_results")
    
    @patch('rnaseqcount.cli.plot_differential_test')
    @patch('rnaseqcount.cli.differential_test')
    @patch('rnaseqcount.cli.writeout')
    def test_diff_analysis_with_groups(self, mock_writeout, mock_diff, mock_plot_diff):
        """测试 --diff 使用 --groups 手动指定两组"""
        mock_fig = MagicMock()
        mock_result = MagicMock()
        mock_plot_diff.return_value = mock_fig
        mock_diff.return_value = mock_result
        
        exit_code, stdout, stderr = self._run_cli([
            '--counts', self.counts_path,
            '--metadata', self.metadata_path,
            '--diff',
            '--groups', 'control', 'treatment',
            '--test_method', 'wilcoxon',
            '--paired',
            '--output', self.output_prefix
        ])
        
        assert exit_code == 0
        assert '执行差异表达分析' in stdout
        
        call_kwargs = mock_plot_diff.call_args[1]
        assert call_kwargs['groups'] == ['control', 'treatment']
        assert call_kwargs['test_method'] == 'wilcoxon'
        assert call_kwargs['paired'] == True
    
    def test_diff_without_group_info(self):
        """测试 --diff 但不指定分组信息时退出并报错"""
        exit_code, stdout, stderr = self._run_cli([
            '--counts', self.counts_path,
            '--metadata', self.metadata_path,
            '--diff',
            '--output', self.output_prefix
        ])
        
        assert exit_code == 1
        assert '错误: 差异表达分析需要指定--grouped_by或--groups' in stderr
    
    # ==================== 参数组合测试 ====================
    
    @patch('rnaseqcount.cli.plot_pca')
    @patch('rnaseqcount.cli.plot_sample_correlation')
    @patch('rnaseqcount.cli.plot_differential_test')
    @patch('rnaseqcount.cli.differential_test')
    @patch('rnaseqcount.cli.writeout')
    def test_multiple_analyses(self, mock_writeout, mock_diff, mock_plot_diff, 
                                mock_plot_corr, mock_plot_pca):
        """测试同时执行多个分析（PCA + 相关性 + 差异表达）"""
        mock_fig = MagicMock()
        mock_result = MagicMock()
        mock_plot_pca.return_value = mock_fig
        mock_plot_corr.return_value = mock_fig
        mock_plot_diff.return_value = mock_fig
        mock_diff.return_value = mock_result
        
        exit_code, stdout, stderr = self._run_cli([
            '--counts', self.counts_path,
            '--metadata', self.metadata_path,
            '--length', self.length_path,
            '--pca',
            '--correlation',
            '--diff',
            '--grouped_by', 'condition',
            '--method', 'FPKM',
            '--filter', '10',
            '--output', self.output_prefix
        ])
        
        assert exit_code == 0
        assert '执行PCA分析' in stdout
        assert '计算样本相关性' in stdout
        assert '执行差异表达分析' in stdout
        
        # 验证各分析都被调用
        mock_plot_pca.assert_called_once()
        mock_plot_corr.assert_called_once()
        mock_plot_diff.assert_called_once()
        mock_diff.assert_called_once()
        
        # 验证参数传递正确
        pca_kwargs = mock_plot_pca.call_args[1]
        assert pca_kwargs['method'] == 'FPKM'
        assert pca_kwargs['filter_threshold'] == 10.0
    
    # ==================== 参数验证测试 ====================
    
    def test_invalid_method_choice(self):
        """测试 --method 传入无效值时 argparse 报错"""
        exit_code, stdout, stderr = self._run_cli([
            '--counts', self.counts_path,
            '--method', 'INVALID',
            '--output', self.output_prefix
        ])
        
        assert exit_code == 2
        assert 'invalid choice' in stderr.lower() or '无效选择' in stderr
    
    def test_invalid_test_method_choice(self):
        """测试 --test_method 传入无效值时 argparse 报错"""
        exit_code, stdout, stderr = self._run_cli([
            '--counts', self.counts_path,
            '--diff',
            '--grouped_by', 'condition',
            '--test_method', 'anova',
            '--output', self.output_prefix
        ])
        
        assert exit_code == 2
        assert 'invalid choice' in stderr.lower() or '无效选择' in stderr


if __name__ == "__main__":
    pytest.main([__file__, "-v"])