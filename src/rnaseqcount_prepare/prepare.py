import pandas as pd
import numpy as np
import gzip
import requests
import time
from pathlib import Path

try:
    import mygene
    MYGENE_AVAILABLE = True
except ImportError:
    MYGENE_AVAILABLE = False

def fetch_ensembl_gene_lengths(gene_ids, species="hsapiens"):
    """通过 Ensembl REST API 获取基因长度（适用于 Ensembl ID）"""
    if not gene_ids:
        return pd.DataFrame(columns=['id', 'length'])
    base_url = "https://rest.ensembl.org"
    lengths = []
    batch_size = 200
    for i in range(0, len(gene_ids), batch_size):
        batch = gene_ids[i:i+batch_size]
        ext = "/lookup/id"
        headers = {"Content-Type": "application/json", "Accept": "application/json"}
        data = {"ids": batch}
        response = requests.post(base_url + ext, headers=headers, json=data)
        if response.status_code != 200:
            print(f"警告：批次 {i//batch_size+1} 获取基因长度失败")
            continue
        result = response.json()
        for gid, info in result.items():
            if 'seq_region_name' in info and 'start' in info and 'end' in info:
                length = info['end'] - info['start'] + 1
                lengths.append({'id': gid, 'length': length})
            else:
                lengths.append({'id': gid, 'length': None})
        time.sleep(0.5)
    length_df = pd.DataFrame(lengths)
    if length_df.empty:
        return length_df
    if length_df['length'].isnull().any():
        median_len = length_df['length'].median()
        length_df['length'].fillna(median_len, inplace=True)
    return length_df

def fetch_gene_lengths_mygene(gene_symbols, species="human"):
    """通过 mygene 库获取基因长度（适用于基因符号，如 A1BG）"""
    if not MYGENE_AVAILABLE:
        raise ImportError("请安装 mygene 库：pip install mygene")
    mg = mygene.MyGeneInfo()
    lengths = []
    batch_size = 1000
    for i in range(0, len(gene_symbols), batch_size):
        batch = gene_symbols[i:i+batch_size]
        results = mg.querymany(batch,
                               scopes='symbol',
                               fields='genomic_pos',
                               species=species,
                               returnall=True,
                               verbose=False)
        for hit in results['out']:
            symbol = hit['query']
            length = None
            if 'genomic_pos' in hit and hit['genomic_pos']:
                gp = hit['genomic_pos']
                if isinstance(gp, list) and len(gp) > 0:
                    gp = gp[0]
                if isinstance(gp, dict) and 'end' in gp and 'start' in gp:
                    length = gp['end'] - gp['start'] + 1
            lengths.append({'id': symbol, 'length': length})
        for miss in results.get('missing', []):
            lengths.append({'id': miss, 'length': None})
        time.sleep(0.5)
    length_df = pd.DataFrame(lengths)
    if length_df.empty:
        return length_df
    if length_df['length'].isnull().any():
        median_len = length_df['length'].median()
        length_df['length'].fillna(median_len, inplace=True)
        print(f"警告：部分基因未获取到长度，已用中位数 {median_len} 填充")
    return length_df

def _find_column_ignore_case(df, col_name):
    """在 DataFrame 列名中查找忽略大小写且去除空格的匹配，返回原始列名"""
    col_name_stripped = col_name.strip().lower()
    for col in df.columns:
        col_stripped = str(col).strip().lower()
        if col_stripped == col_name_stripped:
            return col
    print(f"错误：找不到列 '{col_name}'。可用列名：{list(df.columns)}")
    return None

def prepare_from_user_files(
    metadata_path: str,
    counts_path: str,
    output_dir: str = ".",
    length_source: str = None,
    species: str = "hsapiens",
    sample_id_col: str = "run",
    condition_col: str = "condition"
) -> dict:
    """
    从用户提供的 metadata 文件和原始计数矩阵文件生成 rnaseqcount_wf 所需三个 CSV。
    输出的 metadata.csv 只包含两列：sample_id 和 condition（列名固定）。
    """
    # 读取 metadata（处理 BOM 和去除列名前后的空格）
    metadata = pd.read_csv(metadata_path, encoding='utf-8-sig')
    metadata.columns = metadata.columns.str.strip()
    out_path = Path(output_dir)
    out_path.mkdir(parents=True, exist_ok=True)

    print(f"读取 metadata: {metadata_path}")

    # 查找实际列名（不区分大小写）
    actual_sample_col = _find_column_ignore_case(metadata, sample_id_col)
    if actual_sample_col is None:
        raise ValueError(f"metadata 中找不到类似 '{sample_id_col}' 的列。")
    actual_condition_col = _find_column_ignore_case(metadata, condition_col)
    if actual_condition_col is None:
        raise ValueError(f"metadata 中找不到类似 '{condition_col}' 的列。")

    # 只保留两列，并重命名为固定名称 'sample_id' 和 'condition'
    metadata_out = pd.DataFrame({
        'sample_id': metadata[actual_sample_col],
        'condition': metadata[actual_condition_col]
    })

    # 检查分组是否只有两组（警告）
    groups = metadata_out['condition'].unique()
    if len(groups) != 2:
        print(f"警告：分组列包含 {len(groups)} 个不同值: {groups}")
        print("差异分析要求恰好两组，后续可使用 groups 参数手动指定。")
    else:
        print(f"检测到两组：{groups[0]} vs {groups[1]}")

    # 确保 sample_id 唯一
    if metadata_out['sample_id'].duplicated().any():
        raise ValueError("metadata 中的 sample_id 列有重复值，请确保每个样本唯一。")

    # 读取原始计数矩阵
    print(f"读取原始计数矩阵: {counts_path}")
    if counts_path.endswith('.gz'):
        with gzip.open(counts_path, 'rt') as f:
            first_line = f.readline().strip()
        if first_line.startswith('"') or first_line.startswith("'") or ',' in first_line:
            counts = pd.read_csv(counts_path, compression='gzip', index_col=0)
        else:
            counts = pd.read_csv(counts_path, compression='gzip', sep='\t', index_col=0)
    elif counts_path.endswith('.csv'):
        counts = pd.read_csv(counts_path, index_col=0)
    else:
        counts = pd.read_csv(counts_path, sep='\t', index_col=0)

    if not counts.dtypes.apply(lambda x: pd.api.types.is_integer_dtype(x)).all():
        print("警告：计数矩阵中包含非整数数值，可能不是原始 read counts。")

    sample_names_in_counts = counts.columns.tolist()
    print(f"计数矩阵包含 {len(sample_names_in_counts)} 个样本，{counts.shape[0]} 个基因。")

    # 验证 metadata 与 counts 样本名匹配
    metadata_samples = metadata_out['sample_id'].tolist()
    missing_in_counts = [s for s in metadata_samples if s not in sample_names_in_counts]
    extra_in_counts = [s for s in sample_names_in_counts if s not in metadata_samples]
    if missing_in_counts:
        raise ValueError(f"以下 metadata 中的样本在计数矩阵中找不到: {missing_in_counts}")
    if extra_in_counts:
        print(f"注意：计数矩阵中有 {len(extra_in_counts)} 个样本不在 metadata 中，将被忽略。")
        counts = counts[metadata_samples]
    else:
        print("metadata 与计数矩阵样本名完全匹配。")

    # 输出 counts.csv
    counts_out = counts.reset_index()
    counts_out.rename(columns={counts_out.columns[0]: 'id'}, inplace=True)
    counts_path_out = out_path / "counts.csv"
    counts_out.to_csv(counts_path_out, index=False)
    print(f"生成 counts.csv: {counts_path_out}")

    # 输出 metadata.csv（只含两列）
    metadata_path_out = out_path / "metadata.csv"
    metadata_out.to_csv(metadata_path_out, index=False)
    print(f"生成 metadata.csv: {metadata_path_out}")

    # 处理长度文件
    length_path_out = None
    gene_ids = counts.index.tolist()

    if length_source == 'ensembl':
        print(f"通过 Ensembl API 获取 {len(gene_ids)} 个基因的长度...")
        length_df = fetch_ensembl_gene_lengths(gene_ids, species=species)
        if not length_df.empty:
            length_path_out = out_path / "length.csv"
            length_df.to_csv(length_path_out, index=False)
            print(f"生成 length.csv: {length_path_out}")
        else:
            print("警告：未能获取任何基因长度，请检查基因 ID 是否为 Ensembl 格式。")
    elif length_source == 'mygene':
        species_mygene = 'human' if species.lower().startswith('h') else 'mouse'
        print(f"通过 mygene 获取 {len(gene_ids)} 个基因的长度...")
        length_df = fetch_gene_lengths_mygene(gene_ids, species=species_mygene)
        if not length_df.empty:
            length_path_out = out_path / "length.csv"
            length_df.to_csv(length_path_out, index=False)
            print(f"生成 length.csv: {length_path_out}")
        else:
            print("警告：未能获取任何基因长度，请检查基因 ID 是否为基因符号。")
    elif length_source is not None and Path(length_source).exists():
        print(f"从本地文件复制长度数据: {length_source}")
        length_df = pd.read_csv(length_source)
        if 'id' in length_df.columns and 'length' in length_df.columns:
            length_path_out = out_path / "length.csv"
            length_df.to_csv(length_path_out, index=False)
            print(f"生成 length.csv: {length_path_out}")
        else:
            print("错误：本地长度文件必须包含 'id' 和 'length' 列")
    elif length_source is not None:
        print(f"警告：无法识别的 length_source '{length_source}'，跳过长度文件生成。")
    else:
        print("未指定 length_source，跳过长度文件生成。后续分析请使用 CPM。")

    print("\n===== 处理完成 =====")
    print(f"counts 文件: {counts_path_out}")
    print(f"metadata 文件: {metadata_path_out}")
    if length_path_out:
        print(f"length 文件: {length_path_out}")
    else:
        print("未生成长度文件。")
    print("分组情况:")
    print(metadata_out['condition'].value_counts().to_string())

    return {
        'counts': str(counts_path_out),
        'metadata': str(metadata_path_out),
        'length': str(length_path_out) if length_path_out else None
    }