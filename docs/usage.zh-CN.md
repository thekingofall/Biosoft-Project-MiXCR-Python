# 用法
[English](usage.md) · [简体中文](usage.zh-CN.md) · [首页](../README.zh-CN.md)

## 安装
需要 Python 3.10 及以上版本和 NumPy。在仓库根目录创建并激活虚拟环境，执行 `python -m pip install ./python`。使用 `mixcrpy --version` 和 `mixcrpy --help` 检查安装。`python -m mixcrpy` 与命令行入口等价，可明确使用当前 Python 环境。

## 单端与双端
```bash
mixcrpy analyze reads.fastq.gz output/single --chains TRB
mixcrpy analyze reads_R1.fastq.gz output/paired --read2 reads_R2.fastq.gz --chains BCR
```

链类型可选 ALL、TCR、BCR、TRA、TRB、TRG、TRD、IGH、IGK、IGL，默认 ALL。FASTQ 需包含非空名称、核苷酸序列，以及等长的 Phred+33 质量字符串。

双端文件应记录数相等、名称对应，支持 CASAVA 描述和 /1、/2 后缀。确认文件已经按顺序配对、但名称不同时，可使用 `--pair-by-order`。原始名称参与确定性的模糊碱基处理。`--no-merge` 保留独立目标；支持的双端布局为 Opposite。

## 分开比对与组装
```bash
mixcrpy align reads.fastq output/aligned
mixcrpy assemble output/aligned.alignments.json output/clones
```

默认按 CDR3 组装。比对 JSON 保存链类型和 V 参考特征，供后续步骤读取。`--v-feature VGeneWithP` 可选择另一种 V 基因特征；默认 VRegionWithP。

## 多特征与基因分组
```bash
mixcrpy analyze full_reads.fastq output/cdr123 --features CDR1 CDR2 CDR3
mixcrpy analyze full_reads.fastq output/vdj --features VDJRegion --separate-by V J C
```

输入序列需要覆盖所选特征。`--separate-by` 将相应基因命中纳入克隆分组。验证过的配置包括 CDR3、CDR1+CDR2+CDR3、CDR3+FR4 和 VDJRegion。

## 部分 reads 拼接
```bash
mixcrpy align partial.fastq output/raw --allow-partial-alignments
mixcrpy assemble-partial output/raw.alignments.json output/overlaps
mixcrpy assemble output/overlaps.alignments.json output/clones
mixcrpy assemble-contigs output/overlaps.alignments.json output/clones.clones.json output/contigs
```

双端输入在比对命令中添加 `--read2 partial_R2.fastq`。`--allow-no-cdr3-part-alignments` 可保留符合条件且缺少 CDR3 邻接锚点的部分比对。`assemble-partial --overlapped-only` 输出成功拼接的记录。拼接结果保留所有参与 read 的 ID；克隆组装与 contig 步骤应使用同一份拼接后 JSON。

## Contig 重建
```bash
mixcrpy assemble-contigs output/sample.alignments.json output/sample.clones.json output/contigs
```

使用同一次分析、同一参考配置的比对与克隆文件。当前要求 CDR3 克隆，采用默认 FullSeqAssembler 配置。输出包括 `.contigs.fasta`、`.contigs.tsv`、`.contigs.json` 和 `.report.json`。一个克隆可能包含多个分支或目标。

## 与原版对照
可选的开发验证脚本需要自行取得的、获准使用的原版 MiXCR 3.0.13 JAR，以及 Java 环境。仓库提供便携测试的固定预期结果，原版可执行程序通过官方渠道获取。

```bash
python scripts/compare_example.py --jar /path/to/mixcr-3.0.13.jar --fastq examples/trb.fastq --output results/original-check
```

Java 不在 PATH 中时，添加 `--java /path/to/java`。脚本校验 JAR 的 SHA-256，分别运行两套程序的比对、克隆组装与导出，并保存命令、日志和比较结果。输出目录应使用新路径。检查范围以脚本声明字段为准；历史完整测试见[验证报告](validation.zh-CN.md)。

## 重跑 Python 测试
```bash
python -m pip install pytest
python -m pytest python/tests -q
```

测试读取仓库内由原版程序生成的 JSON 预期结果，运行这些检查无需 Java。重新生成合成示例可执行 `python scripts/make_example.py`。

## 使用边界
输入格式错误会终止命令。覆盖已有结果需要 `--force`。每次分析使用独立前缀；写入过程中被中断的分析可能留下不完整文件，建议换新前缀重跑。

命令行程序将比对记录保存在内存中。处理大数据前，先用有代表性的小样本评估内存和时间。速度等价性与临床性能尚未验证。参考库和参数会影响结果，应使用相同配置与 3.0.13 对照，并参考[验证范围](validation.zh-CN.md)解释差异。

原版 `.vdjca`、`.clns`、`.clna` 文件需要原版工具处理，本项目使用自有 JSON 中间格式。最新版功能与条款见[独立说明](latest-version.zh-CN.md)。
