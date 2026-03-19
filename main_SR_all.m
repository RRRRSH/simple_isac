function main_SR_all(gammatFixed_dB, gammatSweep_dB)
% main_SR_all - 统一入口：生成四大主方案对比图 (Sum-rate 随 P_t, N, \Gamma_t 的变化)
% 用法:
%   main_SR_all()                          % 固定门限图=7 dB, 门限扫描=0:2:10 dB
%   main_SR_all(7, 0:1:12)                % 自定义固定门限与扫描范围

if nargin < 1 || isempty(gammatFixed_dB)
    gammatFixed_dB = 7;
end
if nargin < 2 || isempty(gammatSweep_dB)
    gammatSweep_dB = 0:2:10;
end

% 单次运行统一输出目录
fpath = which(mfilename);
if isempty(fpath), baseDir = pwd; else baseDir = fileparts(fpath); end
outRoot = fullfile(baseDir,'out');
if ~exist(outRoot,'dir'), mkdir(outRoot); end
runStamp = datestr(now,'yyyymmdd_HHMMSS');
gammatTag = strrep(num2str(gammatFixed_dB,'%g'),'.','p');
runDir = fullfile(outRoot,sprintf('main_SR_all_g%sdB_%s', gammatTag, runStamp));
if ~exist(runDir,'dir'), mkdir(runDir); end

outDirPt = fullfile(runDir,'Pt');
outDirN = fullfile(runDir,'N');
outDirGamma = fullfile(runDir,'Gamma');
if ~exist(outDirPt,'dir'), mkdir(outDirPt); end
if ~exist(outDirN,'dir'), mkdir(outDirN); end
if ~exist(outDirGamma,'dir'), mkdir(outDirGamma); end

fprintf('main_SR_all config: fixed=%.3f dB, sweep=[%s]\n', gammatFixed_dB, num2str(gammatSweep_dB));
fprintf('Unified output directory: %s\n', runDir);

main_SR_P(gammatFixed_dB, outDirPt);
main_SR_N(gammatFixed_dB, outDirN);
main_SR_gammat(gammatSweep_dB, outDirGamma);

% 写入运行总览
readmePath = fullfile(runDir,'README.txt');
fid = fopen(readmePath,'w');
if fid > 0
    fprintf(fid,'main_SR_all Unified Run (4 Baselines Mode)\n');
    fprintf(fid,'Time: %s\n', datestr(now,'yyyy-mm-dd HH:MM:SS'));
    fprintf(fid,'Fixed threshold (dB): %.6f\n', gammatFixed_dB);
    fprintf(fid,'Threshold sweep (dB): [%s]\n\n', num2str(gammatSweep_dB));
    fprintf(fid,'Subfolders:\n');
    fprintf(fid,'- Pt/: Sum-rate vs P_t\n');
    fprintf(fid,'- N/: Sum-rate vs N\n');
    fprintf(fid,'- Gamma/: Sum-rate vs Gamma_t\n');
    fclose(fid);
end
end
