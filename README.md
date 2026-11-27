[English](README.en.md)

# Amino Acid Check

Amino Acid Check 用于酶建模前的侧链检查。输入 PDB 后，程序为 His、Asp、Glu、Asn 和 Gln 生成含邻近残基的局部结构。可在离线 3D 视图中比较候选、查看判断依据、测距并记录人工复核结果。

## 下载与运行

从 [GitHub Releases](https://github.com/CXY-one/amino-acid-check/releases/latest) 下载 Windows 10/11 64 位便携 ZIP。完整解压后双击 `Amino Acid Check.exe`，并将 `_internal` 保留在程序旁。便携版已包含运行环境，可选用本机 PyMOL 以棒状打开当前候选。

界面内置简体中文和英文。在「设置 → 界面语言」中选择后点击「应用」即可切换。点击「打开语言文件夹」可以查看或添加 JSON 翻译文件。

## 检查一个结构

1. 拖入 `.pdb` 或 `.ent` 文件，选择输出位置，再点击「运行判断」。可在「设置」中选择残基类型、目标链和邻域距离。
2. 左表可按位点、残基、状态、判断或 pKa 筛选。右侧可切换 His 的 `HID`、`HIE`、`HIP`，Asp/Glu 的质子化状态，以及 Asn/Gln 的酰胺取向。长文本可悬停查看。
3. 用 3D 视图检查局部环境。左键旋转，中键平移，滚轮缩放。点击「测距」后依次选择两个原子。
4. 点击「确认，下一个」保存当前选择。「导出判断表」将所有位点导出到 Excel，按链分表，包含位点、残基、状态和判断四列。可用「打开已有结果」继续复核。

自动判断使用局部结构规则，供建模时参考。pKa 列显示 PROPKA 对 His、Asp 和 Glu 的估算。横线表示缺少判断或估值。

## 命令行

便携包内含 `hischeck-cli.exe`。下例检查 A 链 His，导出三种候选和预览图。输出目录应为尚未创建的目录。

```powershell
.\hischeck-cli.exe .\input.pdb --chain A --residue HIS --all-states --render -o .\output
```

`--propka` 添加 pKa 估值，重复 `--residue` 可检查多种残基。`summary.csv` 汇总状态，位点目录中的 `evidence.json` 记录判断依据，`review.json` 记录已确认的人工作业。

参见[开发者文档](docs/DEVELOPMENT.md)了解源码运行、测试、翻译和打包。项目源码采用 [MIT 许可证](LICENSE)，随附组件见[第三方声明](packaging/windows/THIRD-PARTY.txt)。
