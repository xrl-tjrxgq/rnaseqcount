# 项目结构
src/
├── rnaseqcount/          # 核心分析库
│   ├── core.py           # CountData 容器类
│   ├── io.py             # 数据读写
│   ├── stats.py          # 标准化、过滤、PCA、相关性
│   ├── algorithms.py     # 差异表达分析
│   ├── visualize.py      # 绘图功能
│   ├── cli.py            # 命令行入口
│   └── __init__.py
├── rnaseqcount_wf/       # 一键工作流
│   ├── workflow.py       # run_workflow 函数
│   ├── data_simulator.py # 模拟数据生成
│   ├── cli.py            # 工作流命令行
│   └── __init__.py
└── rnaseqcount_prepare/  # 数据预处理工具
    ├── prepare.py        # 原始数据解析、长度获取
    └── __init__.py

# 安装
git clone https://github.com/yourname/rnaseqcount.git
cd rnaseqcount
pip install -e .
# 尚未上传github

# 使用说明
1.模拟数据测试
from rnaseqcount_wf import run_workflow
# 模拟数据（2组，每组5个样本，1000个基因，20%差异表达）
cdata, diff_res, sig_genes = run_workflow(
    data_source="simulated",
    n_genes=1000,
    n_samples_per_group=5,
    de_ratio=0.2,
    output_dir="./sim_results"
)

2.实际数据处理
2.1.原始数据预处理
from rnaseqcount_prepare import prepare_from_user_files
# 需要提供raw_count_matrix与metadata，输出格式标准化、用于下一步分析的counts, metadata, length_dict文件
result = prepare_from_user_files(
    metadata_path="SraRunTable.csv",          # 包含样本名和分组列
    counts_path="raw_counts.txt.gz",          # 原始计数矩阵（基因×样本）
    output_dir="./prepared_data",
    length_source='mygene',                   # 自动获取基因长度（符号）
    species='hsapiens',
    sample_id_col="Run",                      # metadata 中的样本名列
    condition_col="treatment"                 # metadata 中的分组列，一般分为2组，若非2组则需手动调整以进行比较分析
)
2.2.一键工作流
from rnaseqcount_wf import run_workflow
# 输入rnaseqcount_prepare包预处理结果输出路径，输出样本相关性热图、PCA散点图、火山图
run_workflow(
    data_source ='real',
    input_dir="./rnaseq_demo/real_prepare",  #即rnaseqcount_prepare输出的预处理结果输出路径
    method='TPM',
    filter_threshold=1.0,
    test_method='t_test',
    fdr_correction=True,
    significance_threshold=0.05,
    log2fc_threshold=1.0,
    output_dir='./rnaseq_demo/real_output',
    prefix='real_demo'
)

3.单个功能命令行使用

# 欢迎提交 Issue 和 Pull Request。如有任何问题，请联系项目维护者。