# Worm genome assembly and annotation

本项目按物种分开保存。`solidus/` 是 S. solidus；`scotti/` 是 S. cotti。
`scotti` 沿用服务器和输入文件的目录拼写，物种名称为 **Schistocephalus cotti**。

| 入口 | 内容 | 状态 |
|---|---|---|
| [solidus/](solidus/README.md) | S. solidus 核基因组、mt、结构与功能注释、发布记录 | 原有发布记录；mt 注释存在跨物种核查中发现的待复查边界 |
| [scotti/](scotti/README.md) | S. cotti / scotti_2 的组装及注释入口 | mt 组装完成，注释修订候选待应用与验证 |
| [comparative-genomics/](comparative-genomics/README.md) | 多物种蛋白 QC、BUSCO、OrthoFinder 和功能比较 | 按该目录说明执行 |

```text
Worm-genome-assembly/
├── README.md
├── solidus/
│   ├── README.md
│   └── workflow-scripts/       # 原顶层 workflow-scripts，内容原样保留
├── scotti/
│   ├── README.md
│   └── mitochondria/
│       ├── README.md
│       ├── history_and_caveats.md
│       ├── curation_candidates.tsv
│       └── scripts/
└── comparative-genomics/
```

## 选择正确的物种流程

| 项目 | S. solidus | S. cotti |
|---|---|---|
| 本地代码 | solidus/workflow-scripts/ | scotti/mitochondria/scripts/ |
| 基因组原始长 reads | ONT 流程 | PacBio HiFi，scotti_2 同一个体两次运行合并 |
| mt 组装 | Flye nano-raw + Medaka | Flye pacbio-hifi；seed42 随机700 reads |
| 自组装 mt 长度 | 15,436 bp | 15,437 bp |
| 服务器组装路径 | /work/cyu/assembly/worm/ | /work/cyu/assembly/worm/scotti/ |
| 核注释路径 | /work/cyu/annotation/ssol-annotation/ | /work/cyu/annotation/scotti-annotation/ |

这些是本地代码目录的整理，**服务器数据路径没有改变**。
脚本仍使用原来的服务器输入路径和环境。部分 solidus 脚本是历史运行记录，
不是可从任意目录一键执行的流水线；请先阅读对应说明。
不要把 solidus 的 ONT、Medaka 参数直接用于 cotti HiFi。

S. cotti 的 cox3、trnH 和 rrnL 坐标修订尚未写入最终 GFF；
OH/OL 预测不能直接当作已确认的 D-loop。详见 cotti 的证据记录。
