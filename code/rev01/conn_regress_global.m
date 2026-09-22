function yHat = conn_regress_global(y, sscIdx)
% CONN_REGRESS_GLOBAL  Regress the short-separation-channel signal out of y.
%
% SYNTAX
%   yHat = conn_regress_global(y)
%   yHat = conn_regress_global(y, sscIdx)
%
% INPUT
%   y       Matrix of signals (nTimePoints x nChannels), all channels of the
%           montage, in Homer3 export order.
%   sscIdx  Optional vector of short-separation-channel column indices.
%           Default: the short channels of the 24-channel prefrontal Brite
%           montage, obtained from get_channels_from_template so that the
%           montage definition lives in exactly one place.
%
% OUTPUT
%   yHat    Filtered matrix with the SSC regressor removed, same size as y.
%
% NOTE ON REVISION 01 (Neurophotonics NPH-260108-1)
%   The submitted version hard-coded sscIdx = [3, 14]. Those are the short
%   channel indices of the MOTOR montage defined in get_channels_from_template
%   and were carried over from an earlier study. In the PREFRONTAL montage
%   used here the short channels are 3 (S03D01, right) and 24 (S10D08, left);
%   index 14 (S05D05) is a long channel over right superior frontal gyrus.
%   The submitted analysis therefore regressed one short channel together with
%   one long channel and never regressed the left short channel. The indices
%   are now derived from the montage template rather than hard-coded, and all
%   connectivity results were recomputed.
% ______________________________________________________________________________

if nargin < 2 || isempty(sscIdx)
    % Short channels are the montage channels that are not long channels
    [chIdxL, chIdxR] = get_channels_from_template('prefrontal');
    sscIdx = setdiff(1:size(y, 2), [chIdxL, chIdxR]);
end

% Guard against a montage/template mismatch
if isempty(sscIdx) || any(sscIdx > size(y, 2))
    error('conn_regress_global:badSSC', ...
          'Short-channel indices [%s] are not valid for %d channels.', ...
          num2str(sscIdx), size(y, 2));
end

%% Short-separation-channel regression
% Regressor: mean time course across the short-separation channels
yMean = mean(y(:, sscIdx), 2);
g = [ones(size(yMean)) yMean];

% Regularization parameter
alphaReg = 0.001;
lambda = alphaReg * max(diag(g' * g));

% Tikhonov regularized estimator
betaCoeff = pinv(g' * g + lambda * eye(size(g' * g))) * g' * y;

% Filtered signal
yHat = y - g * betaCoeff;

end

% EOF
