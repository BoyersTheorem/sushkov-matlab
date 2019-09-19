function [meanCalc, stdErr, firstData, secndData, lastPt] = normIndivid_e200(experimentData, nBrt, nDrk, nSig)
%NORMINDIVID Normalize data by individually normalizing the averages 
%   separately and then averaging them together afterwards. Assumes
%   normalizing use bright and dark
%
%   Takes in as input an experimentData structure (which is just a renamed
%   gExperiment) to extract data from, and a list of integers for which
%   measurement is bright, dark, and the signal
%
%   Output is the mean first, followed by the stdErr since that's usually
%   what you'll care about. Also gives the full data in case you care,
%   preseparated, and lastPt
%
%   Version e200 8/5/2019: Now, you input row vectors into the indices
%   nBrt, nDrk, nSig for cases where within the same average you repeat a
%   signal. Often, this is done to decrease noise in bright and dark
%   measurements FOr calculations, in averages the results within an
%   average to get an overall 'bright', 'dark', 'signal', and then does the
%   usual procedure of normalizing, and averaging the 'averages'

% Get the number of repititions of each measurement
NMbrt = length(nBrt);
NMdrk = length(nDrk);
NMsig = length(nSig);

% First, determine the actual number of averages performed. Usually, you'll
% also do this in the analysis script, but having it here saves you from
% having to pass another value

% Assume that NumAv is the maxAv. Then, loop through the averages and check
% if the next average is all 0's. If so, then the current average is the last.
% If you never look ahead and see zeros, the number of averages is the max average.
MaxAv = experimentData.MetaData.Average; % max possible number of avgs
NumAv = MaxAv;
for j=1:MaxAv-1
    if experimentData.Data.AVE(1,j+1).X(1,1).xmean(1) == 0 % Look for an average being identically 0
        NumAv = j;
        break;
    end
end

% Determine the number of points in the sweep by taking length of the sweep
Nswe = length(experimentData.Data.Sweep(1).X);

% Determine the last point where averages were actually taken. Works
% the same as the MaxAv
lastPt = Nswe; % Assume the last point is the end of the sweep. Then, check if you find a zero earlier.
for idx = 1:Nswe-1
    % Look for the next point being 0
    if experimentData.Data.AVE(1,NumAv).X(1,1).xmean(idx+1) == 0
        lastPt = idx;   % If so, the current point is the last, and break out of the loop
        break;
    end
end

% Initialize output. Rows are the averages, columns are points in the sweep
fullData = zeros(NumAv, Nswe);

% Loop through the averages
for idxAv = 1:NumAv
    % Get each set of data by averaging over repeated measurements
    currBr = zeros(1,Nswe);
    currDk = zeros(1,Nswe);
    currSg = zeros(1,Nswe);
    
    for idxMes = 1:NMbrt
        currBr = currBr + experimentData.Data.AVE(idxAv).X(nBrt(idxMes)).xmean(1:Nswe)/NMbrt;
    end
    
    for idxMes = 1:NMdrk
        currDk = currDk + experimentData.Data.AVE(idxAv).X(nDrk(idxMes)).xmean(1:Nswe)/NMdrk;
    end
    
    for idxMes = 1:NMsig
        currSg = currSg + experimentData.Data.AVE(idxAv).X(nSig(idxMes)).xmean(1:Nswe)/NMsig;
    end
    
    % Normalize and insert
    fullData(idxAv,:) = (currSg - currDk)./(currBr - currDk);
end

% Separate the data into parts with the different averages, which can be
% used for calculation separately and then recombined
firstData = fullData(:, 1:lastPt);
if lastPt < Nswe
    secndData = fullData(1:end-1, lastPt+1:end);
else
    secndData = [];
end

% Calculate the mean. The extra input tells it to return a row, meaning
% average down the columns.
meanCalc = [mean(firstData,1), mean(secndData,1)];

% Calculate the standard error (note that the second [if exists] is 1
% shorter than the first). The second input is the weighting, 0 for
% default. The third column tells you to return a row, meaning taking the
% std. down the columns.
stdErr   = [std(firstData,0,1)/sqrt(NumAv), std(secndData,0,1)/sqrt(NumAv-1)];

end

