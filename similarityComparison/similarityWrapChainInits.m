function init = similarityWrapChainInits(stim, xAinit, xBinit, xCinit, xDinit, nuInit, withPred)
%SIMILARITYWRAPCHAININITS Chain inits for wrap-copy similarity models.
%
%   mu is length nStimuli (wrap copies use wrapOff in JAGS). One nu[t] per trial,
%   matching perceptualReproduction_jags.txt.

init = struct(...
   'xA', xAinit, 'xB', xBinit, 'xC', xCinit, 'xD', xDinit, ...
   'nu', nuInit, 'mu', stim(:));
if nargin >= 7 && withPred
   init.xAp = xAinit;
   init.xBp = xBinit;
   init.xCp = xCinit;
   init.xDp = xDinit;
end
end
