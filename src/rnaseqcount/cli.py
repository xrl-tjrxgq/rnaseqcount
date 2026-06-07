# cli.py
import argparse
import sys
from .core import CountData
from .io import writein, writeout
from .stats import filter_low_expression, pca, sample_correlation, group_samples
from .algorithms import differential_test, pick_significance
from .visualize import plot_pca, plot_sample_correlation, plot_differential_test


def main():
    parser = argparse.ArgumentParser(description='RNA-seq数据分析工具')
    parser.add_argument('--counts', type=str, help='counts矩阵文件路径(.csv)')
    parser.add_argument('--metadata', type=str, help='样本元数据文件路径(.csv)')
    parser.add_argument('--length', type=str, help='基因长度文件路径(.csv)')
    parser.add_argument('--filter', type=float, default=0, help='低表达过滤阈值')
    parser.add_argument('--method', type=str, default='TPM', choices=['TPM', 'FPKM', 'CPM'], 
                       help='归一化方法')
    parser.add_argument('--pca', action='store_true', help='执行PCA分析')
    parser.add_argument('--correlation', action='store_true', help='计算样本相关性')
    parser.add_argument('--diff', action='store_true', help='执行差异表达分析')
    parser.add_argument('--grouped_by', type=str, help='metadata中的分组列名')
    parser.add_argument('--groups', type=str, nargs=2, help='手动指定两组名称，如 control treatment')
    parser.add_argument('--test_method', type=str, default='t_test', choices=['t_test', 'wilcoxon'],
                       help='统计检验方法')
    parser.add_argument('--paired', action='store_true', help='是否进行配对检验')
    parser.add_argument('--output', type=str, default='output', help='输出文件前缀')
    
    args = parser.parse_args()
    
    # 读取数据
    print("正在读取数据...")
    cdata = writein(counts_path=args.counts, metadata_path=args.metadata, length_path=args.length)
    print(f"数据加载成功: {cdata.shape[0]}个基因, {cdata.shape[1]}个样本")
    
    # PCA分析
    if args.pca:
        print("执行PCA分析...")
        fig = plot_pca(cdata, method=args.method, filter_threshold=args.filter)
        writeout(fig, f"{args.output}_pca")
        print("PCA图已保存")
    
    # 相关性分析
    if args.correlation:
        print("计算样本相关性...")
        fig = plot_sample_correlation(cdata, method=args.method, filter_threshold=args.filter)
        writeout(fig, f"{args.output}_correlation")
        print("相关性热图已保存")
    
    # 差异表达分析
    if args.diff:
        if not args.grouped_by and not args.groups:
            print("错误: 差异表达分析需要指定--grouped_by或--groups", file=sys.stderr)
            sys.exit(1)
        
        groups = None
        if args.groups:
            groups = list(args.groups)
        
        print(f"执行差异表达分析 (方法: {args.test_method})...")
        fig = plot_differential_test(cdata, method=args.method, 
                                     filter_threshold=args.filter,
                                     grouped_by=args.grouped_by, groups=groups,
                                     test_method=args.test_method, paired=args.paired)
        writeout(fig, f"{args.output}_volcano")
        
        result = differential_test(cdata, method=args.method,
                                   filter_threshold=args.filter,
                                   grouped_by=args.grouped_by, groups=groups,
                                   test_method=args.test_method, paired=args.paired)
        writeout(result, f"{args.output}_diff_results")
        print("差异表达分析完成")


if __name__ == "__main__":
    main()