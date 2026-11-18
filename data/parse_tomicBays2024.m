% parse_tomicBays2024.m
%
% OUTPUT (this script's only product):
%   tomicBays2024.mat   — structs ds (similarity), dp (perception), dm (memory)
%
% Canonical packaging of Tomic & Bays orientation CSVs for all model drivers.
% Perception/memory: divide target/stimulus/response by 2 before mod(., pi)
% (those columns were stored as 2 * orientation). Memory nontargets: fold only.
% Similarity `ds` is reused from a local parse_1 product (same fold as similarity).
%
% Prerequisites (local only; not on GitHub):
%   Run parse_1.m once so tomicBays.mat exists (provides ds + stimulus grid).
%   This script does not overwrite tomicBays.mat.
%
% CSV convention (Tomic; confirmed with collaborator):
%   Orientations sampled on [-pi/2, pi/2).
%   Perception/memory target and response stored as 2 * orientation on [-pi, pi).
%   Memory nontarget_* left on [-pi/2, pi/2).

clear;

outFile = 'tomicBays2024.mat';   % <-- sole output of this parse

pi = 3.1415926535;
decimalPoints = 4;

if ~isfile('tomicBays.mat')
   error('parse_tomicBays2024:missingParse1', ...
      ['Run parse_1.m first (local only). This script reuses ds/stimuli ', ...
      'from tomicBays.mat, then writes %s.'], outFile);
end
load tomicBays ds;

step = diff(ds.stimuli(1:2));
binsE = [(0 - step/2):step:(pi + step/2) pi];

%% perception: undo *2 on stimulus/response, then fold
T = readtable('oriPercData_clean.csv');
dp.participant = T.subject;
dp.response = mod(round(T.response, decimalPoints) / 2, pi);
dp.stimulus = mod(round(T.target, decimalPoints) / 2, pi);
[dp.nTrials, ~] = size(dp.participant);
dp.nParticipants = length(unique(dp.participant));
dp.sIdx = zeros(dp.nTrials, 1);
for i = 1:dp.nTrials
   dp.sIdx(i) = find(histcounts(dp.stimulus(i), binsE));
   choices = dp.response(i) + [0 -pi pi];
   [~, chooseIdx] = min(abs(choices - dp.stimulus(i)));
   dp.response(i) = choices(chooseIdx);
end
dp.sIdx(dp.sIdx == 73) = 1;
dp.stimuli = ds.stimuli;
dp.nStimuli = ds.nStimuli;

%% memory: undo *2 on target/response; nontargets fold only
T = readtable('oriWMData_clean.csv');
dm.response = mod(round(T.response, decimalPoints) / 2, pi);
dm.target = mod(round(T.target, decimalPoints) / 2, pi);
dm.participant = T.subject;
[dm.nTrials, ~] = size(dm.participant);

nt = [T.nontarget_1, T.nontarget_2, T.nontarget_3, T.nontarget_4, T.nontarget_5];
dm.nontarget = mod(round(nt, decimalPoints), pi);
dm.nNontargets = 5 - sum(isnan(dm.nontarget), 2);
dm.setSize = T.setsize;

dm.tIdx = zeros(dm.nTrials, 1);
for i = 1:dm.nTrials
   dm.tIdx(i) = find(histcounts(dm.target(i), binsE));
   choices = dm.response(i) + [0 -pi pi];
   [~, chooseIdx] = min(abs(choices - dm.target(i)));
   dm.response(i) = choices(chooseIdx);
end
dm.tIdx(dm.tIdx == 73) = 1;

dm.nIdx = nan(dm.nTrials, 5);
for i = 1:dm.nTrials
   for j = 1:5
      if ~isnan(dm.nontarget(i, j))
         dm.nIdx(i, j) = find(histcounts(dm.nontarget(i, j), binsE));
      end
   end
end
dm.nIdx(dm.nIdx == 73) = 1;

dm.stimuli = ds.stimuli;
dm.nStimuli = ds.nStimuli;

save(outFile, 'ds', 'dp', 'dm');
fprintf('Wrote %s (ds from parse_1; dp+dm with /2-before-mod for target/stimulus/response)\n', ...
   outFile);
