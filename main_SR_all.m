function main_SR_all(varargin)
% main_SR_all - 统一入口：生成三张四方案对比图
% 用法:
%   main_SR_all()                          % 固定门限图=7 dB, 门限扫描=0:2:10 dB
%   main_SR_all(7, 0:1:12)                % 自定义固定门限与扫描范围
%   main_SR_all(5)                        % 固定门限=5 dB，扫描用默认范围
%   main_SR_all(7, 0:2:10, 'noris_only')  % 快速仅跑No-RIS
%   main_SR_all(7, 0:2:10, 'noris_vs_fixed') % 仅比较No-RIS与RIS固定相位
%   main_SR_all(7, 0:2:10, 'noris_vs_fixed_nsp_vs_mrt') % 比较No-RIS/RIS+MRT/RIS+NSP
%
% 输入:
%   gammatFixed_dB: 用于 main_SR_P / main_SR_N 的固定门限(dB)
%   gammatSweep_dB: 用于 main_SR_gammat 的门限扫描向量(dB)
%   runMode: 'all'(默认) / 'noris_only' / 'noris_vs_fixed' / 'noris_vs_fixed_nsp_vs_mrt'

if nargin >= 1 && ~isempty(varargin{1})
    gammatFixed_dB = varargin{1};
else
    gammatFixed_dB = 7;
end

if nargin >= 2 && ~isempty(varargin{2})
    gammatSweep_dB = varargin{2};
else
    gammatSweep_dB = 0:2:10;
end

if nargin >= 3 && ~isempty(varargin{3})
    runMode = char(varargin{3});
else
    runMode = 'all';
end

% 单次运行统一输出目录
fpath = which(mfilename);
if isempty(fpath)
    baseDir = pwd;
else
    baseDir = fileparts(fpath);
end
outRoot = fullfile(baseDir,'out');
if ~exist(outRoot,'dir'), mkdir(outRoot); end
runStamp = datestr(now,'yyyymmdd_HHMMSS');
gammatTag = strrep(num2str(gammatFixed_dB,'%g'),'.','p');
runDirName = sprintf('main_SR_all_%s_g%sdB_%s', runMode, gammatTag, runStamp);
runDir = fullfile(outRoot,runDirName);
if ~exist(runDir,'dir'), mkdir(runDir); end

outDirPt = fullfile(runDir,'Pt');
outDirN = fullfile(runDir,'N');
outDirGamma = fullfile(runDir,'Gamma');
outDirHeatmap = fullfile(runDir,'Heatmap');
if ~exist(outDirPt,'dir'), mkdir(outDirPt); end
if ~exist(outDirN,'dir'), mkdir(outDirN); end
if ~exist(outDirGamma,'dir'), mkdir(outDirGamma); end
if ~exist(outDirHeatmap,'dir'), mkdir(outDirHeatmap); end

fprintf('main_SR_all config: fixed=%.3f dB, sweep=[%s], mode=%s\n', ...
    gammatFixed_dB, num2str(gammatSweep_dB), runMode);
fprintf('Unified output directory: %s\n', runDir);

main_SR_P(gammatFixed_dB,runMode,outDirPt);
main_SR_N(gammatFixed_dB,runMode,outDirN);
main_SR_gammat(gammatSweep_dB,runMode,outDirGamma);
if strcmpi(runMode,'noris_vs_fixed_nsp_vs_mrt')
    plot_spatial_heatmap_nsp_vs_mrt(outDirHeatmap,gammatFixed_dB);
end

% 写入运行总览，提升目录可读性
readmePath = fullfile(runDir,'README.txt');
fid = fopen(readmePath,'w');
if fid > 0
    fprintf(fid,'main_SR_all unified run\n');
    fprintf(fid,'time: %s\n', datestr(now,'yyyy-mm-dd HH:MM:SS'));
    fprintf(fid,'mode: %s\n', runMode);
    fprintf(fid,'fixed threshold (dB): %.6f\n', gammatFixed_dB);
    fprintf(fid,'threshold sweep (dB): [%s]\n\n', num2str(gammatSweep_dB));
    fprintf(fid,'Subfolders:\n');
    fprintf(fid,'- Pt/: Sum-rate vs P_t outputs\n');
    fprintf(fid,'- N/: Sum-rate vs N outputs\n');
    fprintf(fid,'- Gamma/: Sum-rate vs Gamma_t outputs\n');
    fprintf(fid,'- Heatmap/: Spatial heatmaps\n');
    if strcmpi(runMode,'noris_vs_fixed')
        fprintf(fid,'\nMode note: this run only computes and plots two curves: RIS-Fixed and No-RIS.\n');
    elseif strcmpi(runMode,'noris_vs_fixed_nsp_vs_mrt')
        fprintf(fid,'\nMode note: this run computes and plots No-RIS, RIS+MRT, RIS+NSP and generates spatial heatmaps.\n');
    end
    fclose(fid);
end
end
