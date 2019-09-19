function [dataOut, stdvOut] = binAvgSeqntl_e200(dataIn, edgeDistribution, Nbins)
% The idea is the following. If I want to bin non-uniformly, I have to give
% some kind of distribution from which to draws the bins from. The way this
% works, is I sample linearly from the distribution, which will give me
% increasing points. These points them become the edges of my binning.
% Thus, a flatter region means that I have more closely spaced edges so
% that there will be fewer bins there. Qualitatively, the number of points
% per bin in a region is proportional to the derivative of the distribution
% in that region. By convention, the edgeDistribution has domain and range
% [0,1]. In practice I would like to do this without requiring a 'time'
% input for each 'signal' input, so I'll just assume all the input data is
% temporally order already, and then I can bin based on their indices,
% which runs from 1 to N.
%
% I then take the input independent variable and determine which bin each
% data point falls in using the built-in discretize. Once I have an 
% assignment of data point to bin, I can use that to average together the
% points in a given bin.

% Get the length of the data 
Npts    = length(dataIn);

% Sample the distribution to get the edges
Nedges   = Nbins + 1;                  % For fence posting, you have one more edge than bin
sampVec  = 0:1/(Nedges - 1):1;         % Vector to use to sample the distribution
edgesNrm = edgeDistribution(sampVec);  % Sample the distribution to get the edges
edgesAct = (Npts - 1)*edgesNrm + 1;    % Get the actual edges by multiplying by length and adding offset

% Use the edges to assign a bin to each data point
idxvec   = 1:Npts;                           % Create a vector of the indices
binAssgn = discretize(idxvec, edgesAct);    % Use builtin to get the edges

% Loop through the bins and determine the average and stddev of each.
dataOut  = zeros(1, Nbins);                  % Initialize the outputs
stdvOut  = zeros(1, Nbins);

for idx = 1:Nbins
    binPoints    = find(binAssgn == idx);    % For each bin point, find the indices where it occurs
    dataOut(idx) = mean(dataIn(binPoints));  % Extract those points, take mean, insert into bin
    stdvOut(idx) =  std(dataIn(binPoints));  % Do the same for the standard deviation
end

end

