%% Intro
% Script for analyzing echo data, based on Eric's version of echo_20-40-20, where measurements
% are bright, echo (20-40-20), dark, echo (20-40-60). 
%
% This version divides the differential echo data by the intrinsic noise of
% the NV (spin bath, resonances, etc.) so that the only decoherence is from
% the added noise. 
%
% Version 130 adds the threshold method

%% Initialize, load data, user input here

fileName = 'NV16_T2';
load(fileName);

% Hardware calibration values
whichNV    = 12;    % Which NV center to use, for determining the delta from the DC voltage.
toff       = .04;  % timing offset in us. What is the starting duration from sequence?
Tpi        = 32; % Time of 'stated' pi pulse in ns. This is the distance between peaks, not the 'actual' time used in a sequence, which is offset to account for pulse error
drVoltPi   = .1812; % Vpp Amplitude of AWG to achieve same amplitude on 'I' input as the 'Q' input to SRS

% Fitting Options
pinMean    = 1;
pinAmpl    = 1;
pinPower   = 1;
echoPower  = 1;
divideIntr = 0;     % 1 to divide the data by some 'intrinsic' decay with these param, 0 to fit data as normal

% Intrinsic fit parameters:
aIntr      = 0;
bIntr      = 1;
cIntr      = 70;   % Time in us
pIntr      = 1;

% Threshold level
threshold = exp(-1);

%% Extract parameters and data from files
StartT          = experimentData.MetaData.StartDateTime;     % scan start time
StartTime       = datestr(StartT,'yyyy-mm-dd HH:MM:SS');     % scan start time string
StartDay        = datestr(StartT,'yyyy-mm-dd');              % scan start day
Samp            = experimentData.MetaData.Samples;           % samples per point
StartCts        = experimentData.MetaData.InitialCounts/1e3; % starting counts
SeqName         = otherParam.SEQ.name;                       % Seq name
SeqName         = strrep(SeqName, '_','\_');                 % Seq modified for figure
AnalysisName    = mfilename;                                 % This script
AnalysisName    = strrep(AnalysisName, '_', '\_');
fileNameStr     = strrep(fileName, '_', '\_');

SRS1powe        = experimentData.Parameters.MWPower;         %
SRS1freq        = experimentData.Parameters.MWFreq;

driveVolt       = otherParam.AWG.ampl;   %Vpp                   % Amplitude (in Vpp) of the LZ drive
drivePeriod     = otherParam.AWG.period;                    % Period of the LZ drive in seconds
drVoltNum       = str2double(driveVolt);  %Vpp
drPerdNum       = str2double(drivePeriod);

dcAmpStr        = otherParam.AWG.ampl2;
dcOffStr        = otherParam.AWG.offset2;
dcAmp           = str2double(dcAmpStr);
dcOff           = str2double(dcOffStr);

% Determine actual number of averages.
MaxAv = experimentData.MetaData.Average; % max possible number of avgs
for j=1:MaxAv
    if experimentData.Data.AVE(1,j).X(1,1).xmean(1) == 0
        break;
    end
end
NumAv = j-1;

% Extract Data
x       = DataMat(1,:)*1e6; 
bright  = DataMat(2,:);
ramsey1 = DataMat(3,:);
dark    = DataMat(4,:);
ramsey2 = DataMat(5,:);

xsamp = min(x):(max(x)-min(x))/1000:max(x);

% Get average contrast
avgContrast = 100*(mean(bright) - mean(dark))/mean(bright);
contrStr = num2str(avgContrast,3);


%% Plot raw
% Custom Colors
purp = [ .8,  .2, .8]; % Richer purple
gref = [ .2,  .9, .2]; % 'Fuller' green - more saturated
cyaf = [ .2,  .8, .8]; % 'Fuller' cyan  - more saturated
oraf = [.95, .65, .1]; % fuller orange
redb = [.8,   .2, .2]; % 'Burnt" red

% Convert numbers to strings for title
SampStr      = num2str(Samp, '%5.0f');
NumAvStr     = num2str(NumAv,'%3.0f');
StartCtsStr  = num2str(round(StartCts, 3, 'significant'));

SRS1poweStr  = num2str(SRS1powe,'%4.2f');
SRS1freqStr  = num2str(SRS1freq*1e-9,'%7.5f');

drVoltstr    = num2str(drVoltNum, '%3.3f');
drPerdstr    = num2str(drPerdNum, '%3.3f');

% Create Raw Title String (rts)
rts = {
       ['Differential Echo.     Data file: ', fileNameStr];
       ['Seq: ', SeqName, ',  ','Analysis: ', AnalysisName]
       ['Start: ', StartTime, ',  ', 'Avgs: ', NumAvStr, ',  ', 'Samp/pt: ', SampStr,...
            ',  ', 'Start kcts: ', StartCtsStr];
       ['SRS1 Freq: ', SRS1freqStr, ' GHz,  ', 'SRS1 Power: ', SRS1poweStr];
       
       };

figure('pos', [100, 100, 1200, 600])
subplot(1,2,1)
hold on
grid on
xlabel('Added Time (\mus)')
ylabel('Raw counts')
title(rts)

plot(x, bright,    '.-r',                'MarkerSize', 13, 'Linewidth', 1.2)
plot(x, dark,      '.-b',                'MarkerSize', 13, 'Linewidth', 1.2)
plot(x, ramsey1,   '.-',  'Color', oraf, 'MarkerSize', 13, 'Linewidth', 1.2) 
plot(x, ramsey2,   '.-',  'Color', gref, 'MarkerSize', 13, 'Linewidth', 1.2)

legend('Bright', 'Dark', 'Ramsey 1', 'Ramsey 2', 'location', 'SouthEast')
hold off

%% Analysis  - fit to exponential
% Normalize
tTot = 2*(x + toff); % The total experiment length is the added time plus offset, times two since adds before and after the pi pulse
[echo1N, err1] = normIndivid(experimentData, 1, 3, 2);
[echo2N, err2] = normIndivid(experimentData, 1, 3, 4);

% Subtract data to get differential signal
echoDf = (echo1N - echo2N);
errDf  = sqrt(err1.^2 + err2.^2);

% Create intrinsic noise fit and divide it out
if divideIntr
    fIntrDcy = @(t) (aIntr + bIntr.*exp(-(t./cIntr).^pIntr));
    echoAdded = echoDf./fIntrDcy(tTot);
else
    echoAdded = echoDf;
end

% Prep the data
[tp, y1p] = prepareCurveData(tTot, echoAdded);

% Set bounds and start points
%%% First Ramsey (starts dark) %%%
if pinMean
    aLow = 0;
    aUpp = 0;
else
    aLow = -.3;
    aUpp =  .3;
end

if pinAmpl
    bLow = 1;
    bUpp = 1;
else
    bLow = 0;
    bUpp = 2;
end

if pinPower
    pLow = echoPower;
    pSta = echoPower;
    pUpp = echoPower;
else
    pLow = 0;
    pSta = echoPower;
    pUpp = 10;
end

%%%% First Echo (starts bright) %%%
lowBound = [  aLow,   bLow,     0, pLow];
startPt  = [    .5,     .5,    10, pSta];
upBound  = [  aUpp,   bUpp,  1000, pUpp];

% Set fittype and options

ft1 = fittype( 'a + b.*exp(-(x/c).^d)', 'independent', {'x'}, 'dependent', {'y'},...
               'coefficients', {'a', 'b', 'c', 'd'});
opts1            = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts1.Display    = 'Off';
opts1.StartPoint = startPt;
opts1.Lower      = lowBound;
opts1.Upper      = upBound;


% Perform fitting
[fit1, gof1] = fit(tp, y1p, ft1, opts1);

% Extract Parameters
a1 = fit1.a;
b1 = fit1.b;
c1 = fit1.c;
d1 = fit1.d;

% Get confidence intervals/errors
int1   = confint(fit1, 0.66);
int1_a = int1(:,1);
int1_b = int1(:,2);
int1_c = int1(:,3);
int1_d = int1(:,4);

a1err  = abs(int1_a(2) - int1_a(1))/2;
b1err  = abs(int1_b(2) - int1_b(1))/2;
c1err  = abs(int1_c(2) - int1_c(1))/2;
d1err  = abs(int1_d(2) - int1_d(1))/2;

T2_1 = c1;
T2_err1 = c1err;


%% Analysis - threshold method
% Calculated normalized differential signal above

% Subtract the threshold level and take the signum to find above and below.
aboveBelow = sign(echoAdded - threshold);

% Look for places where it exactly hit (0's) rare, but could happen
zeroIdxs   = find(aboveBelow == 0);
zeroTime   = tTot(zeroIdxs);

% Now look for places where a crossing occurs. First, subtract previous
% element from each to look for a change in sign
signChng   = diff(aboveBelow);

% Find places where the magnitude is 2, corresponding to a place where one
% is +-1 and the previous is -+1, a sign change. This finds the second part
% of the change, so then need to subtract by one to get the two parts of
% the sign change. Note that diff ignores the first point in the vector, so
% that the point x(j) - x(j-1) will get put in j-1, so that the output of
% diff corresponds to the first point.
nChanges   = length(find(abs(signChng) == 2)); % Find how many there are
chngIdxs   = zeros(2, nChanges);               % Initialize
chngIdxs(1,:)  = find(abs(signChng) == 2);     % Get the first point in each pair
chngIdxs(2,:)  = chngIdxs(1,:) + 1;            % Add 1 to get the second point

% Determine the times of crossings for each pairs
chngTimes  = (tTot(chngIdxs(1,:)) + tTot(chngIdxs(2,:)))/2;

% Average together the times 
crossMean  = mean([zeroTime, chngTimes]); % Put zeros together and take mean
crossStdv  =  std([zeroTime, chngTimes]);

% Get stepsize
dt         = tTot(2) - tTot(1);

%% Plot analysis
% Convert parameters to strings
a1str    = num2str(a1,     '%.3g');
b1str    = num2str(b1,     '%.3g');
c1str    = num2str(c1,     '%.3g');
d1str    = num2str(d1,     '%.3g');

a1Errstr = num2str(a1err,  '%.2g');
b1Errstr = num2str(b1err,  '%.2g');
c1Errstr = num2str(c1err,  '%.2g');
d1Errstr = num2str(d1err,  '%.2g');

apn      = num2str(pinMean);
bpn      = num2str(pinAmpl);
ppn      = num2str(pinPower); 
dvInt    = num2str(divideIntr);

T2_1str = num2str(T2_1, '%3.3f');
T2e1str = num2str(T2_err1, '%.2g');

% Threshold values
TthrshStr    = num2str(crossMean, '%3.3f');
TthrshErrStr = num2str(crossStdv, '%3.3f');
dtstr        = num2str(dt       , '%3.3f');

% Create fit plot
xsamp = 0:max(tTot)/(1000-1):max(tTot);

fitSim1 = feval(fit1, xsamp);

% Plot
ats = {
       ['Pins: a, b, p:  ', apn, ', ', bpn, ', ', ppn, ', ', '  Divide Intrinsic: ', dvInt];
       ['Fit model: a + b.*exp(-(x/c)^p)   p = ' d1str, ' \pm ', d1Errstr];
       ['a = ', a1str, ' \pm ', a1Errstr,'  b = ', b1str, ' \pm ', b1Errstr,...
        '  T_2 = ', T2_1str, ' \pm ', T2e1str];
       ['Threshold Method: ', TthrshStr, ' \pm ', TthrshErrStr, '   \Deltat = ', dtstr];
       };

subplot(1,2,2)
hold on
grid on
xlabel('Time during sequence (\mus)')
ylabel('Population')
title(ats)

errorbar(tTot, y1p, errDf, 'o', 'Color', cyaf, 'CapSize', 0, 'MarkerFaceColor', cyaf, 'MarkerEdgeColor', 'k', 'MarkerSize', 6, 'Linewidth', 1.2)
plot(xsamp, fitSim1, '-k',                                   'Linewidth', 1.5)

plot([0, max(tTot)], [threshold, threshold], '--',  'Color', [.4, .4, .4], 'Linewidth', 1);
plot(tTot(chngIdxs(1,:)), y1p(chngIdxs(1,:)), 'rx', 'MarkerSize', 16)
plot(tTot(chngIdxs(2,:)), y1p(chngIdxs(2,:)), 'rx', 'MarkerSize', 16)
plot(tTot(zeroIdxs), echoDf(zeroIdxs), 'rx', 'MarkerSize', 24)

axis([0, max(tTot), -.1, 1.1])
legend('Differential Echo', 'Fit', 'Threshold', 'Cross points', 'location', 'NorthEast')







