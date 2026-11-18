function [chains, info] = appendChainStruct(a, b)
%APPENDCHAINSTRUCT Concatenate coda structs along the chain (column) dim.
%
% Trinity coda fields are nSamples x nChains. If two waves differ in
% nSamples (e.g. Select seed 1500 vs weekend wave 2000), both are truncated
% to the shared minimum so concatenation can proceed.

info = struct('nSamplesUsed', [], 'truncatedA', false, 'truncatedB', false, ...
   'skippedFields', {{}}, 'onlyA', {{}}, 'onlyB', {{}});

if isempty(a)
   chains = b;
   if ~isempty(b)
      fn = fieldnames(b);
      info.nSamplesUsed = size(b.(fn{1}), 1);
   end
   return
end
if isempty(b)
   chains = a;
   fn = fieldnames(a);
   info.nSamplesUsed = size(a.(fn{1}), 1);
   return
end

fnA = fieldnames(a);
fnB = fieldnames(b);
common = intersect(fnA, fnB, 'stable');
info.onlyA = setdiff(fnA, fnB, 'stable');
info.onlyB = setdiff(fnB, fnA, 'stable');

nKeep = inf;
for i = 1:numel(common)
   nKeep = min(nKeep, size(a.(common{i}), 1));
   nKeep = min(nKeep, size(b.(common{i}), 1));
end
if ~isfinite(nKeep) || nKeep < 1
   error('appendChainStruct:nSamples', 'No common fields with samples to append');
end
info.nSamplesUsed = nKeep;

chains = struct();
for i = 1:numel(common)
   name = common{i};
   Aa = a.(name);
   Bb = b.(name);
   if size(Aa, 1) > nKeep
      info.truncatedA = true;
      Aa = Aa(1:nKeep, :);
   end
   if size(Bb, 1) > nKeep
      info.truncatedB = true;
      Bb = Bb(1:nKeep, :);
   end
   if size(Aa, 1) ~= size(Bb, 1)
      info.skippedFields{end+1} = name; %#ok<AGROW>
      continue
   end
   chains.(name) = [Aa Bb];
end

% Keep only-A / only-B fields truncated to nKeep (chain metadata may differ)
for i = 1:numel(info.onlyA)
   name = info.onlyA{i};
   Aa = a.(name);
   if size(Aa, 1) >= nKeep
      chains.(name) = Aa(1:nKeep, :);
   else
      chains.(name) = Aa;
   end
end
for i = 1:numel(info.onlyB)
   name = info.onlyB{i};
   Bb = b.(name);
   if size(Bb, 1) >= nKeep
      chains.(name) = Bb(1:nKeep, :);
   else
      chains.(name) = Bb;
   end
end
end
