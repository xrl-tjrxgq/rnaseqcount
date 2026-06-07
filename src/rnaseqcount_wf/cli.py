import argparse
from .workflow import run_workflow

def main():
    parser = argparse.ArgumentParser(description="RNA-seq 一键分析工作流")
    parser.add_argument("--data-source", choices=["simulated", "real"], default="simulated",
                        help="数据来源")
    # 真实数据参数
    parser.add_argument("--counts", type=str, help="counts 文件路径（真实数据）")
    parser.add_argument("--metadata", type=str, help="metadata 文件路径（真实数据）")
    parser.add_argument("--length", type=str, help="基因长度文件路径（真实数据）")
    # 模拟数据参数
    parser.add_argument("--n-genes", type=int, default=1000)
    parser.add_argument("--n-samples-per-group", type=int, default=5)
    parser.add_argument("--n-groups", type=int, default=2)
    parser.add_argument("--de-ratio", type=float, default=0.2)
    parser.add_argument("--seed", type=int, default=42)
    # 分析参数
    parser.add_argument("--method", choices=["TPM","FPKM","CPM"], default="TPM")
    parser.add_argument("--filter-threshold", type=float, default=1.0)
    parser.add_argument("--grouped-by", type=str, default="condition")
    parser.add_argument("--test-method", choices=["t_test","wilcoxon"], default="t_test")
    parser.add_argument("--paired", action="store_true")
    parser.add_argument("--no-fdr", dest="fdr_correction", action="store_false")
    parser.add_argument("--sig-threshold", type=float, default=0.05)
    parser.add_argument("--log2fc-threshold", type=float, default=1.0)
    # 输出
    parser.add_argument("--output-dir", default="./rnaseq_results")
    parser.add_argument("--prefix", default="workflow")

    args = parser.parse_args()

    run_workflow(
        data_source=args.data_source,
        counts_path=args.counts,
        metadata_path=args.metadata,
        length_path=args.length,
        n_genes=args.n_genes,
        n_samples_per_group=args.n_samples_per_group,
        n_groups=args.n_groups,
        de_ratio=args.de_ratio,
        seed=args.seed,
        method=args.method,
        filter_threshold=args.filter_threshold,
        grouped_by=args.grouped_by,
        test_method=args.test_method,
        paired=args.paired,
        fdr_correction=args.fdr_correction,
        significance_threshold=args.sig_threshold,
        log2fc_threshold=args.log2fc_threshold,
        output_dir=args.output_dir,
        prefix=args.prefix
    )

if __name__ == "__main__":
    main()