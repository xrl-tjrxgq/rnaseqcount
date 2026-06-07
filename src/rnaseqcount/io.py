# io.py
import pandas as pd
import matplotlib.pyplot as plt
from .core import CountData


def writein(counts_path: str = None, metadata_path: str = None, length_path: str = None):
    ## 从CSV文件读取数据并创建CountData对象
    ## counts_path: counts矩阵文件路径，第一列需命名为"id"
    ## metadata_path: 元数据文件路径，第一列需命名为"sample_id"
    ## length_path: 基因长度文件路径，第一列"id"，第二列"length"
    ## 返回: CountData实例
    # 标准读取（自动检测逗号分隔）
    counts_df = metadata_df = length_dict = None

    # 读取 counts
    if counts_path is not None:
        if not counts_path.endswith(".csv"):
            raise Exception("仅支持.csv格式")
        counts_df = pd.read_csv(counts_path)
        if counts_df.columns[0].lower() != "id":
            raise ValueError("counts文件第一列列名不是 'id'")
        counts_df.set_index(counts_df.columns[0], inplace=True)

    # 读取 metadata（第一列名应为 "sample_id"）
    if metadata_path is not None:
        if not metadata_path.endswith(".csv"):
            raise Exception("仅支持.csv格式")
        metadata_df = pd.read_csv(metadata_path)
        if metadata_df.columns[0].lower() != "sample_id":
            raise ValueError("metadata文件第一列列名不是 'sample_id'")
        metadata_df.set_index(metadata_df.columns[0], inplace=True)

    # 读取长度文件（可选）
    if length_path is not None:
        if not length_path.endswith(".csv"):
            raise Exception("仅支持.csv格式")
        length_df = pd.read_csv(length_path)
        length_df.columns = length_df.columns.str.lower()
        if length_df.columns[0] != "id":
            raise ValueError("长度文件第一列列名不是 'id'")
        if length_df.columns[1] != "length":
            raise ValueError("长度文件第二列列名不是 'length'")
        length_df.set_index(length_df.columns[0], inplace=True)
        length_dict = length_df['length'].to_dict()

    return CountData(counts=counts_df, metadata=metadata_df, length=length_dict)


def writeout(file, file_name:str):
    ## 将DataFrame或matplotlib Figure保存为文件
    ## file: pd.DataFrame或plt.Figure对象
    ## file_name: 输出文件名（不含扩展名）
    ## DataFrame保存为CSV，Figure保存为PNG
    if isinstance(file, pd.DataFrame):
        csv_file=f"{file_name}.csv"
        file.to_csv(csv_file, index=False, encoding='utf-8-sig')
        print(f"DataFrame已保存为: {csv_file}")
    if isinstance(file, plt.Figure):
        img_file=f"{file_name}.png"
        file.savefig(img_file, dpi=300, bbox_inches='tight')
        print(f"图片已保存为: {img_file}")# io.py