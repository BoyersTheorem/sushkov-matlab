function [vecOut, errOut] = binAverage_e100(vec, Npbn, errType, varargin)
%BINAVERAGE_E100 Bins points together to reduce a vec in length by a factor
%of Npbn
%  
% vec      :   The vector you want to bin-average
% Npbn     :   Number of points per bin
% errType  :   'stddev' or 'stderr': How should the error be reported
% varargin :   Optionally, provide 'external err' to use for the extra
%              points that weren't binned.

% Check if the vector is divisable by the # of points per bin.
Npts = length(vec);
q    = floor(Npts/Npbn);
r    = mod(Npts, Npbn);

% Do you want the error to be std dev or std error?
switch errType
    case 'stddev'
        errFactor = 1;
    case 'stderr'
        errFactor = 1/sqrt(Npbn);
end

% Initialize according to whether the bins divide the input
if r == 0
    % Bins divide the points, so can include all points in the binning
    % initialize the outputs which will be combined
    rsltBeg = [];
    rsltMid = zeros(1,q);
    rsltEnd = [];
    
    % initialize the recast vector I'll actually use
    dataUse = vec;
    
    % Also initialize the error
    serrBeg = [];
    serrMid = zeros(1,q);
    serrEnd = [];
else
    % If no, exclude points at the begginning and end to have something
    % that is divisable.
    % Determine how many points to put at the beginning. Divide extra
    % points by two, put half at beginning, and the extra at the beggging 
    qExtra  = floor(r/2);
    rExtra  = mod(r, 2);
    rsltBeg = vec(1:(qExtra + rExtra));
    rsltMid = zeros(1,q);
    rsltEnd = vec((end - qExtra+1):end);
    
    dataUse = vec((qExtra + rExtra + 1):(end - qExtra));
    
    % For the error, if an 'external error' is provided with the variable
    % input, use that for the beginning and end. Else, just set them to
    % zero
    if nargin > 3
        errIn   = varargin{1};
        serrBeg = errIn(1:(qExtra + rExtra));
        serrEnd = errIn((end - qExtra+1):end);
    else
        serrBeg = zeros(1, qExtra + rExtra);
        serrEnd = zeros(1, qExtra);
    end
    
end

% Loop through the old vector
for idx = 1:q
    % Average together the points
    rsltMid(idx) = mean(dataUse( ((idx-1)*Npbn + 1):(idx*Npbn) ));
    
    % Calculate the standard deviation
    serrMid(idx) = std( dataUse( ((idx-1)*Npbn + 1):(idx*Npbn) ))*errFactor;
end

% Put the extra parts (if any) together to create the output
vecOut = [rsltBeg, rsltMid, rsltEnd];
errOut = [serrBeg, serrMid, serrEnd];

end

