# simple_isac（中文版说明）

## 项目简介
本项目用于仿真 RIS 辅助的 ISAC（通信感知一体化）系统，重点分析在不同 RIS 相位策略、不同 RIS 单元数 N、不同功率分配系数 eta 下的通信性能与感知性能权衡关系。

核心输出指标：
- 通信总速率（Sum-rate）
- 雷达信噪比（Radar SNR）

## 代码结构
- run_simple_isac.m
  - 主实验脚本。
  - 完成多策略对比、全局相位扫描、eta 权衡分析并输出图像。
- plot_eta_tradeoff.m
  - 独立绘制 eta 权衡图（固定 N=64）。
- utils/
  - compute_phi.m：RIS 相位计算（固定目标/随机）
  - design_w.m：通信与雷达波束设计（ZF + MRT）
  - generate_baseline.m：基线信道与几何参数生成
  - generate_ris_parts.m：给定 N 时生成 RIS 相关信道部分
  - scan_phi.m：RIS 相位扫描优化
  - select_eta.m：满足雷达门限约束的 eta 选择
  - sum_rate.m：通信总速率计算
  - radar_snr.m：含 RIS 的雷达 SNR 计算
  - radar_snr_noris.m：无 RIS 的雷达 SNR 计算

## 当前模型说明（重要）
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

## 已完成的关键更新
- 雷达波束改为基于复合感知信道的 MRT（直达 + RIS 反射）。
- 归一化写法修复为 norm(vec) + 1e-12（避免方向偏移）。
- eta 搜索范围扩展到 [0, 1]，并加入不可行告警。
- scan_phi 相位网格去重（移除与 0 等价的 2pi 点）。
- sum_rate 支持可选雷达干扰项。
- radar_snr / radar_snr_noris 支持可选通信干扰项。
- generate_baseline.m 中移除函数内固定随机种子，随机性由主脚本统一控制。
- No-RIS 基线在主脚本中移到循环外只计算一次，避免重复计算与语义误导。

## 运行方法
### 方式一：完整主流程
运行 run_simple_isac.m。

将输出：
- N 对比图（速率、雷达 SNR）
- 全局相位扫描图
- eta 权衡图

### 方式二：仅绘制 eta 权衡图
运行 plot_eta_tradeoff.m。

## 输出与归档规则
当前脚本已改为按“运行时间戳”自动建子目录保存结果，便于每次修改后对比：

- out/run_simple_isac_YYYYMMDD_HHMMSS/
- out/plot_eta_tradeoff_YYYYMMDD_HHMMSS/

额外归档（手动打包）可放在：
- archives/out_snapshot_YYYYMMDD_HHMMSS.zip

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
