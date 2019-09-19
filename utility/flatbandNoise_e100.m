function [noiseVec] = flatbandNoise_e100(rms,wvec,phsvec,tvec)
%FLATBANDNOISE_E100 Creates noise with flat bandwidth at input time
%   rms:    RMS of noise
%   wvec:   Vector of randomly sampled frequencies
%   pasvec: Vector of randomly sampled phases for each term
%   tvec:   Vector of times you want the noise at

% Get number of components used
N = length(wvec);

% Make sure tvec is a column and wvec, phasvec are rows
if isrow(tvec)
    tvec = tvec.';
end

if ~isrow(phsvec)
    phsvec = phsvec.';
end

if ~isrow(wvec)
    wvec = wvec.';
end

% Compute the arguments as w*t+phi as a big matrix with components going
% across (columns) and time going down (rows). Then, take the cos of all of
% this, and sum across the rows to get field as a function of time.
% Finally, multiply by rms and a normalization factor
noiseVec = sqrt(2/N)*rms*(sum(cos(tvec*wvec + phsvec),2));

end

