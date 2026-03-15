# simple_isac（中文版说明）

## 项目简介
本项目用于仿真 RIS 辅助的 ISAC（通信感知一体化）系统，重点分析不同基线策略下的通信-感知权衡。当前主口径为：固定雷达门限 7 dB，在满足感知约束前提下最大化通信总速率。

核心输出指标：
- 通信总速率（Sum-rate）
- 雷达信噪比（Radar SNR）

## 代码结构
- main_SR_all.m
  - 统一入口，一次运行输出 Pt/N/Gamma 三张图。
  - 支持模式：all / noris_only / noris_vs_fixed / noris_vs_fixed_nsp_vs_mrt。
- main_SR_P.m, main_SR_N.m, main_SR_gammat.m
  - 三张主图脚本（分别对应发射功率、RIS 单元数、门限扫描）。
- plot_spatial_heatmap_nsp_vs_mrt.m
  - 绘制 No-RIS / RIS+MRT / RIS+NSP 的空间热力图。
- run_simple_isac.m
  - 主实验脚本。
  - 完成多策略对比、全局相位扫描并输出图像。
- run_noris_only.m
  - 快速 No-RIS 流程（固定基站功率 1W）。
  - 仅优化 No-RIS 波束，终端输出阶段进度与迭代进度。
- run_noris_full.m
  - No-RIS 完整流程（优化 + 热力图 + 方向图 + eta 图 + 数据保存）。
- plot_eta_tradeoff.m
  - 独立绘制 eta 权衡图（固定 N=64）。
- utils/
  - compute_phi.m：RIS 相位计算（固定目标/随机）
  - design_w.m：通信与雷达波束设计（ZF + 雷达波束，支持 MRT/NSP，默认 NSP）
  - generate_baseline.m：基线信道与几何参数生成
  - generate_ris_parts.m：给定 N 时生成 RIS 相关信道部分
  - scan_phi.m：RIS 相位扫描优化
  - select_eta.m：满足雷达门限约束的 eta 选择
  - optimize_w_for_fixed_phi.m：固定 phi 下，按 SNR 约束最大化通信速率
  - run_joint_snr_optimization.m：联合优化（W 与 phi）
  - get_W_phi_SNR.m：联合优化核心迭代求解器
  - sum_rate.m：通信总速率计算
  - radar_snr.m：含 RIS 的雷达 SNR 计算
  - radar_snr_noris.m：无 RIS 的雷达 SNR 计算

## 当前模型说明（重要）
统一雷达门限：
- gamma_r = 7 dB（线性值 10^(7/10)）

项目中提供了两个可选开关，用于控制是否计入双向干扰：

1) includeRadarInterferenceInRate
- 位置：run_simple_isac.m、plot_eta_tradeoff.m
- 含义：是否将“雷达信号对通信用户的干扰”计入 Sum-rate 分母。
- false：不计入（便于与历史结果对齐）
- true：计入（更完整的并发 ISAC 模型）

2) includeCommInterferenceInRadar
- 位置：run_simple_isac.m、plot_eta_tradeoff.m
- 含义：是否将“通信信号对雷达回波的干扰”计入 Radar SNR 分母。
- false：不计入（便于与历史结果对齐）
- true：计入（更完整的并发 ISAC 模型）

建议：
- 需要与旧结果对比时：两个开关都设为 false。
- 需要物理更完整模型时：两个开关都设为 true。

优化口径说明：
- 当前主流程（run_simple_isac）中，fixed/random/scan/joint 策略均按“先满足雷达 SNR 门限，再优化通信总速率”的逻辑运行。
- eta 搜索与启发式 ZF+雷达波束设计仍保留在工程中，主要用于基线、初始化和对比分析。

三线模式说明（noris_vs_fixed_nsp_vs_mrt）：
- No-RIS：不使用 RIS，仅基站波束。
- RIS+MRT：固定 RIS 相位下，雷达波束采用 MRT。
- RIS+NSP：固定 RIS 相位下，雷达波束采用 NSP（通信零空间投影）。
- 当某点 eta=1 仍无法满足门限时，会提示不可行告警，曲线中该点按不可行处理。

## 已完成的关键更新
- 雷达门限统一固定为 7 dB（gamma_r = 10^(7/10)）。
- 主流程 fixed/random/scan 波束设计已切换为固定 phi 下的 SNR 约束速率优化。
- 雷达波束支持基于复合感知信道的 MRT 与 NSP，两者可切换（默认 NSP）。
- 归一化写法修复为 norm(vec) + 1e-12（避免方向偏移）。
- eta 搜索范围扩展到 [0, 1]，并加入不可行告警。
- scan_phi 相位网格去重（移除与 0 等价的 2pi 点）。
- sum_rate 支持可选雷达干扰项。
- radar_snr / radar_snr_noris 支持可选通信干扰项。
- generate_baseline.m 中移除函数内固定随机种子，随机性由主脚本统一控制。
- No-RIS 基线在主脚本中移到循环外只计算一次，避免重复计算与语义误导。

## 运行方法
### 方式零：统一入口（推荐）
运行 main_SR_all.m。

示例：
- main_SR_all()
- main_SR_all(7, 0:2:10, 'noris_only')
- main_SR_all(7, 0:2:10, 'noris_vs_fixed')
- main_SR_all(7, 0:2:10, 'noris_vs_fixed_nsp_vs_mrt')

说明：
- 第1个参数：固定门限（dB），用于 Pt/N 图。
- 第2个参数：门限扫描向量（dB），用于 Gamma 图。
- 第3个参数：运行模式。

在 noris_vs_fixed_nsp_vs_mrt 模式下会额外自动生成热力图。

### 方式一：完整主流程
运行 run_simple_isac.m。

将输出：
- N 对比图（速率、雷达 SNR）
- 全局相位扫描图
- eta 权衡图

### 方式二：仅绘制 eta 权衡图
运行 plot_eta_tradeoff.m。

### 方式三：快速 No-RIS 验证
运行 run_noris_only.m。

### 方式四：No-RIS 完整流程
运行 run_noris_full.m。

## 输出与归档规则
当前脚本已改为按“运行时间戳”自动建子目录保存结果，便于每次修改后对比。

统一入口 main_SR_all 的目录示例：
- out/main_SR_all_<mode>_g<fixed>dB_<timestamp>/
  - Pt/
  - N/
  - Gamma/
  - Heatmap/
  - README.txt

其他脚本目录示例：

- out/run_simple_isac_YYYYMMDD_HHMMSS/
- out/plot_eta_tradeoff_YYYYMMDD_HHMMSS/
- out/run_noris_full_YYYYMMDD_HHMMSS/

额外归档（手动打包）可放在：
- archives/out_snapshot_YYYYMMDD_HHMMSS.zip

仓库管理说明：
- archives/ 已加入 .gitignore（默认不纳入版本控制）
- MATLAB 临时文件 *.asv、*.m~ 已加入 .gitignore
- 建议将“代码变更”和“结果归档”分离管理：代码用 git，结果用 out/ 与 archives/

## 对比实验建议
可做两组实验对照：

A. 历史口径（理想化）
- includeRadarInterferenceInRate = false
- includeCommInterferenceInRadar = false

B. 完整并发口径（双向干扰）
- includeRadarInterferenceInRate = true
- includeCommInterferenceInRadar = true

然后对比两次运行生成的时间戳目录中的曲线差异。

## 注意事项
- 本项目默认随机种子在顶层脚本中设置（当前为 rng(1)），便于复现实验。
- 若做蒙特卡洛统计，请在外层循环中控制随机种子策略。
- 若更换系统参数（M/K/N/P/噪声功率等），建议保留输出目录时间戳机制，避免结果混淆。

## 快速检查清单
每次改完代码后建议按下面顺序执行：

1. 优先运行 main_SR_all(7,0:2:10,'noris_vs_fixed_nsp_vs_mrt')
2. 确认 out/ 下生成新的时间戳子目录
3. 对比新旧目录中的关键图（Pt/N/Gamma/Heatmap）
4. 如需留档，压缩 out/ 到 archives/out_snapshot_时间戳.zip
5. 提交代码前确认 git status 干净（结果文件不应进入提交）
