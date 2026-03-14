# Figure Caption Templates for simple_isac

## Usage
- Replace placeholders in braces, e.g., {N0}, {theta_t}, {delta_gain_db}.
- Keep notation consistent with the paper: RIS, No RIS, Sum-rate, Radar SNR, eta.
- If needed, keep one concise caption and move detailed interpretation to the main text.

## 1) Radar Spatial Beampattern

### 1.1 Short caption (for figure list)
Figure X. Radar spatial beampattern comparison versus angle under different strategies (No RIS, RIS-Fixed, RIS-Random, RIS-Scan-SNR).

### 1.2 Standard caption (recommended under figure)
Figure X shows the radar spatial beampattern over scanning angles from -90 deg to 90 deg for No RIS and RIS-assisted designs. The RIS-assisted fixed-phase strategy forms a sharper mainlobe around the target direction ({theta_t} deg), while maintaining lower sidelobe levels in non-target directions, especially near user directions ({theta_u_list} deg). This demonstrates that RIS improves directional energy focusing and suppresses undesired spatial leakage.

### 1.3 Chinese version (thesis style)
图X给出了角度范围 -90° 到 90° 下的雷达空间方向图，对比了 No RIS、RIS-Fixed、RIS-Random 和 RIS-Scan-SNR 等策略。可以看到，RIS-Fixed 在目标方向（{theta_t}°）形成更尖锐的主瓣，同时在非目标方向（尤其通信用户方向 {theta_u_list}° 附近）旁瓣更低，表明 RIS 能够提升目标方向能量聚焦并抑制空间泄露。

### 1.4 Quantitative sentence templates (optional)
- At the target direction ({theta_t} deg), the RIS-Fixed gain is approximately {target_gain_db} dB, which is {delta_target_vs_noris_db} dB higher than No RIS.
- Around user directions ({theta_u_list} deg), the average sidelobe level of RIS-Fixed is reduced by about {delta_leak_db} dB compared with No RIS.

## 2) Pareto Frontier (Radar SNR vs Sum-rate)

### 2.1 Short caption (for figure list)
Figure Y. Pareto frontier of ISAC tradeoff between Radar SNR and Sum-rate for RIS and No RIS schemes.

### 2.2 Standard caption (recommended under figure)
Figure Y presents the Pareto frontier in the Radar SNR-Sum-rate plane by sweeping eta from 0 to 1. Compared with the No RIS baseline, the RIS-assisted frontier is shifted outward toward the upper-right region, indicating an enlarged achievable performance region. This confirms that RIS provides simultaneous gains for sensing and communication under joint power allocation.

### 2.3 Chinese version (thesis style)
图Y展示了通过扫描 eta in [0,1] 得到的 Radar SNR 与 Sum-rate 帕累托前沿。与 No RIS 基线相比，RIS 曲线整体向右上方外扩，说明系统可达性能区域被显著扩大，即在相同通信速率下可获得更高雷达 SNR，或在相同雷达 SNR 下可获得更高通信速率。

### 2.4 Quantitative sentence templates (optional)
- For a fixed Sum-rate = {sr_ref} bps/Hz, RIS achieves about {delta_snr_db} dB higher Radar SNR than No RIS.
- For a fixed Radar SNR = {snr_ref_db} dB, RIS provides about {delta_sr} bps/Hz higher Sum-rate than No RIS.
- The area dominated by RIS over No RIS in the SNR-rate plane indicates a clear expansion of the ISAC performance boundary.

## 3) Cross-reference template in main text
As illustrated in Fig. X and Fig. Y, the RIS design not only improves point-wise sensing quality (higher target-direction beampattern gain) but also enhances the global communication-sensing tradeoff (outward-shifted Pareto frontier), validating the effectiveness of RIS-enabled joint ISAC optimization.

## 4) One-paragraph Chinese template (ready to paste)
如图X所示，在空间方向图中，RIS-Fixed 在目标方向形成更强主瓣并有效抑制非目标方向旁瓣，说明其具备更优的能量聚焦能力与更低的空间泄露。进一步地，如图Y所示，通过扫描 eta 得到的帕累托前沿表明 RIS 曲线整体相对 No RIS 向右上方外扩，验证了 RIS 在通信与感知联合优化中的双重增益：在不牺牲通信速率的条件下提升雷达探测质量，或在维持雷达性能的同时提高系统通信效率。
