Amino Acid Check 0.9.2  Windows 便携版

适用系统为 64 位 Windows 10 和 Windows 11。请先完整解压，再双击 Amino Acid Check.exe。软件已包含运行环境和离线三维查看器，_internal 文件夹需放在两个 exe 旁边。

将 PDB 或 ENT 拖入窗口，或点击「选择文件」，然后点击「运行判断」并确认提示。程序建议将新结果保存到「文档\HisCheck」，输出目录可修改。可在「设置」中选择要处理的 His、Asp、Glu、Asn、Gln。右侧可在同一位点切换候选结构，通过 3D 视图查看原子、测距并查阅判断依据。

左表显示「位点」「残基」「状态」「判断」「pKa」，可按任一列筛选，鼠标悬停可查看单元格全文。Asp/Glu 去质子化、Asn/Gln 不翻转显示为「推荐」，质子化或翻转显示为「待复核」。pKa 由 PROPKA 估算，可在设置中关闭。选好状态后点击「确认，下一个」，左表显示「已核对」并转到下一个可见的待核对位点。确认按钮右侧的箭头会在下方悬浮显示「取消核对状态」。「打开已有结果」可继续核对。点击「导出判断表」保存 Excel，按链分表并导出全部位点的前四列，已核对的记录采用人工判断，其余记录采用程序判断。

选中位点后，点击「用 PyMOL 打开」可在本机 PyMOL 中以棒状显示当前候选。软件会自动搜索已安装的 PyMOL，也可在「设置」中粘贴程序路径或安装目录，点击「应用」后会记住该路径。PyMOL 是可选的独立软件。

内置简体中文和英文。在「设置 → 界面语言」选择后点击「应用」，界面立即切换。点击同一行的「打开语言文件夹」可访问 exe 旁的 locales。以 zh-CN.json 或 en.json 为模板添加 UTF-8 JSON 文件，保留消息键和格式占位符。重新展开语言下拉框即可选择新语言。

可在 PowerShell 中用同目录的命令行工具批量判断。
  .\hischeck-cli.exe .\input.pdb --chain A --residue HIS --residue ASP --residue GLU --residue ASN --residue GLN --render --all-states -o .\output
  .\hischeck-cli.exe --help

判断详情与候选状态记录在每个位点的 evidence.json，整次自动判断见 summary.csv。人工确认状态和对应的 PDB 文件名记录在各位点的 review.json。

Quick start

Extract the entire ZIP on 64-bit Windows 10 or 11, then open Amino Acid Check.exe. Keep _internal beside the executables. Open Settings, select English under Interface Language, and click Apply. To add a translation, click Open Language Folder and use en.json or zh-CN.json as a template.
