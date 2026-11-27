# 开发者文档 / Developer guide

[简体中文 README](../README.md) · [English README](../README.en.md)

## 本地运行 / Local setup

Windows x64 上使用 Python 3.11。下面的 PowerShell 命令在源码根目录执行。桌面依赖包含 Qt WebEngine 和可选的 PROPKA 功能。只运行基本 CLI 时，安装 `requirements.txt` 即可。

Use Python 3.11 on Windows x64. Run these PowerShell commands at the root of the source checkout. The desktop dependencies include Qt WebEngine and PROPKA. Install only `requirements.txt` for the basic CLI.

```powershell
py -3.11 -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements.txt -r requirements-gui.txt
.\.venv\Scripts\python.exe -m hischeck.gui
.\.venv\Scripts\python.exe -m hischeck --help
```

也可双击 `launch_gui.vbs`，从 `.venv` 无控制台启动桌面版。设置 `AMINO_ACID_CHECK_LANGUAGE=en` 可在 CLI 中使用英文。

On Windows, `launch_gui.vbs` starts the desktop app from `.venv` without a console. Set `AMINO_ACID_CHECK_LANGUAGE=en` to use English in the CLI.

## 代码入口 / Code map

| 文件 / File | 用途 / Purpose |
| --- | --- |
| [`hischeck/__main__.py`](../hischeck/__main__.py) | CLI 参数、整次分析及结果写入 / CLI arguments, analysis, and run output |
| [`hischeck/structure.py`](../hischeck/structure.py) | PDB 读取、邻域提取、局部结构写入 / PDB parsing, neighborhoods, and local structure export |
| [`hischeck/analysis.py`](../hischeck/analysis.py), [`hischeck/sidechain_analysis.py`](../hischeck/sidechain_analysis.py) | His 与其余四类侧链的判断 / Rules for His and the other four residue types |
| [`hischeck/propka_estimates.py`](../hischeck/propka_estimates.py) | 独立的 PROPKA pKa 估算 / Separate PROPKA pKa estimates |
| [`hischeck/gui.py`](../hischeck/gui.py), [`hischeck/review.py`](../hischeck/review.py) | 桌面流程、人工复核、Excel 导出 / Desktop workflow, review, and Excel export |
| [`hischeck/view3d.py`](../hischeck/view3d.py), [`hischeck/web/`](../hischeck/web/) | 3D 场景数据与离线 3Dmol.js 界面 / Scene data and the offline 3Dmol.js viewer |

默认在目标官能团 5 Å 内搜索证据，并在 6 Å 范围内提取完整邻近残基。`run.json` 记录输入信息和参数，`summary.csv` 和 `evidence.json` 记录程序结果，`review.json` 记录人工选定的候选。桌面版调用同一 CLI 生成结果。PROPKA 数据供复核参考，规则判断单独计算。

The default radii are 5 Å for evidence and 6 Å for exporting complete neighboring residues, both measured from the target functional group. `run.json` records the input and settings. `summary.csv` and `evidence.json` record automatic results, and `review.json` records human selections. The desktop app runs the same CLI. PROPKA estimates are calculated separately from rule-based decisions.

## 测试和翻译 / Tests and translations

```powershell
.\.venv\Scripts\python.exe -m unittest discover -s tests -v
node --check hischeck/web/viewer.js
```

测试使用临时 PDB 数据。第二条命令需要 Node.js。翻译位于 [`hischeck/locales`](../hischeck/locales/)，可复制 `en.json` 或 `zh-CN.json`，修改语言名称和文案，保留消息键和 `{name}` 等格式占位符。在界面重新展开语言列表后可选择新增语言。

The tests generate their own PDB fixtures. The second command requires Node.js. To add a translation, copy `en.json` or `zh-CN.json` in [`hischeck/locales`](../hischeck/locales/), translate the name and messages, and keep the message keys and placeholders such as `{name}`. Reopen the language picker to find the new pack.

## 构建 Windows 便携版 / Windows build

已验证的构建环境为 Windows x64 和 Python 3.11.9。安装上述运行依赖后，使用版本锁定的构建依赖和 PyInstaller 配置。

The verified build environment is Windows x64 with Python 3.11.9. Install the runtime dependencies above, then use the pinned build tools and PyInstaller spec.

```powershell
.\.venv\Scripts\python.exe -m pip install -r packaging/windows/requirements-build.txt
.\.venv\Scripts\python.exe -m PyInstaller --noconfirm --clean packaging/windows/HisCheck.spec
```

`dist/Amino Acid Check/` 包含桌面和 CLI 可执行文件及共享的 `_internal`。分发时保留整目录，附上 [`README.txt`](../packaging/windows/README.txt)、[第三方声明及其所列许可](../packaging/windows/THIRD-PARTY.txt)和项目许可证。把 `hischeck/locales` 复制到可执行文件旁可提供可编辑语言包。发布前应从新路径打开解压后的程序，检查 3D 视图、语言切换和 Excel 导出。

`dist/Amino Acid Check/` contains the desktop and CLI executables and their shared `_internal` directory. Distribute the entire directory with [`README.txt`](../packaging/windows/README.txt), the [third-party notices and listed licenses](../packaging/windows/THIRD-PARTY.txt), and the project license. Copy `hischeck/locales` beside the executables to provide editable language packs. Before release, test the extracted application in a new location, including its 3D viewer, language switching, and Excel export.
