%% Intro
% This analyzes the duration needed to repolarize an NV spin after a pulse
% (typically pi/2) has been applied to it. For the sequence
% Laser_tau_repolarize_e100
%
% Version e200: 8/5/2019: This version is when there is a consistent
% (correct) dark measurement during the whole sweep.

%% User input

fileName = 'Laser_duration_8_6_13_3.mat';
fdt      = load(fileName);

tAvgStrt   = 3;            % Average the signal starting at this time
pinMean    = 0;            % Fix the fitting to go to 1

pErrTrgtdv = 1e-4;         % Targeted polarization (expressed as a difference from 1)
pErrTrgtbd = 1e-3;         % Targeted polarization (expressed as a difference from 1)

tOffset    = .0;           % Offset time to be added to duration to account fo polarization during readout
Nignore    = 0;           % Number of points at beginning to ignore in fitting due to pulse sequence issues with very short laser pulses.

Nsamp      = 2000;         % Number of samples for applications where you resample time

%% Header
% Some operators I might use in the analysis
Sx = [0, 1;...
      1, 0];
Sy = [0, -1i;...
      1i, 0];
Sz = [1,  0;...
      0, -1];
Ux = @(phi) expm(-1i*phi*Sx);
Uy = @(phi) expm(-1i*phi*Sy);
Uz = @(phi) expm(-1i*phi*Sz);

% Custom Colors
purp = [ .8,  .2, .8]; % Richer purple
gref = [ .2,  .9, .2]; % 'Fuller' green - more saturated
cyaf = [ .2,  .8, .8]; % 'Fuller' cyan  - more saturated
oraf = [.95, .65, .1]; % fuller orange
redb = [.8,   .2, .2]; % 'Burnt" red

% Brights
pur2 = [ .9,  .3, .9]; % Brighter purple
cya2 = [ .3,  .9, .9]; %
ora2 = [.98, .75, .2];

% Darks
cyad = [.02,  .5, .5];
orad = [.7,   .4, .015];
purd = [.4,  .02, .4];
gred = [.1,   .7  .1];

%% Extract parameters and data
StartT          = fdt.experimentData.MetaData.StartDateTime;     % scan start time
StartTime       = datestr(StartT,'yyyy-mm-dd HH:MM:SS');     % scan start time string
StartDay        = datestr(StartT,'yyyy-mm-dd');              % scan start day
Samp            = fdt.experimentData.MetaData.Samples;           % samples per point
StartCts        = fdt.experimentData.MetaData.InitialCounts/1e3; % starting counts
SeqName         = fdt.otherParam.SEQ.name;                       % Seq name
SeqName         = strrep(SeqName, '_','\_');                 % Seq modified for figure
AnalysisName    = mfilename;                                 % This script
AnalysisName    = strrep(AnalysisName, '_', '\_');
fileNameStr     = strrep(fileName, '_', '\_');

SRS1powe        = fdt.experimentData.Parameters.MWPower;         %
SRS1freq        = fdt.experimentData.Parameters.MWFreq;

% Determine actual number of averages.
MaxAv = fdt.experimentData.MetaData.Average; % max possible number of avgs
for j=1:MaxAv
    if fdt.experimentData.Data.AVE(1,j).X(1,1).xmean(1) == 0
        break;
    end
end
NumAv = j-1;

% Extract raw data
dt       = fdt.DataMat(1,:)*1e6 + tOffset; % Time in us. 
bright1  = fdt.DataMat(2,:);
bright2  = fdt.DataMat(3,:);
darkVar  = fdt.DataMat(4,:);
polrRaw  = fdt.DataMat(5,:);
dark1    = fdt.DataMat(6,:);
dark2    = fdt.DataMat(7,:);

%% Analysis - division based normalization
% Division based normalization
brtAvg   = (bright1 + bright2)/2;
drkAvg   = (dark1   +   dark2)/2;
polrNdv  =  polrRaw./brtAvg;

% Get mean after user defined starting point
avgStIdx = find((dt >= tAvgStrt), 1);
meanNdv  = mean(polrNdv(avgStIdx:end));
stdeNdv  = std(polrNdv(avgStIdx:end))/sqrt(length(polrNdv(avgStIdx:end)));

% Fit to exponential
if pinMean
    aLw  = 1;
    aUp  = 1;
else
    aLw  =  .5;
    aUp  = 1.5;
end

ftTypedv = fittype( 'a - b.*exp(-x/c)', 'independent', {'x'}, 'dependent', {'y'},...
               'coefficients', {'a', 'b', 'c'});
optsdv            = fitoptions( 'Method', 'NonlinearLeastSquares' );
optsdv.Display    = 'Off';
optsdv.StartPoint = [   .5,  .5,    10];
optsdv.Lower      = [  aLw,   0,     0];
optsdv.Upper      = [  aUp,   1, 10000];

[tFitdv,   yFitdv]  = prepareCurveData(dt(Nignore+1:end), polrNdv(Nignore+1:end));
[fitResdv,  gofdv]  = fit(tFitdv, yFitdv, ftTypedv, optsdv);

tvec   = dt(Nignore+1):(dt(end) - dt(Nignore+1))/(Nsamp -1):dt(end);

% Get results
aFitdv   = fitResdv.a;
bFitdv   = fitResdv.b;
cFitdv   = fitResdv.c;
fitCrvdv = feval(fitResdv, tvec);

% Get uncertainty
intResdv    = confint(fitResdv, 0.66);
intRes_adv  = intResdv(:,1);
aErrdv      = abs(intRes_adv(2) - intRes_adv(1))/2;

intRes_bdv  = intResdv(:,2);
bErrdv      = abs(intRes_bdv(2) - intRes_bdv(1))/2;

intRes_cdv  = intResdv(:,3);
cErrdv      = abs(intRes_cdv(2) - intRes_cdv(1))/2;

% Time for given polarization using fit
tTrgtdv     = cFitdv*log(bFitdv/pErrTrgtdv);

%% Analysis -  Bright-dark based normalization
% Normalize invididually:
nBrt = [1, 2];
nDrk = [5, 6];
nVar = 3;
nRep = 4;
[dVarNbd, dVarErr] = normIndivid_e200(fdt.experimentData, nBrt, nDrk, nVar);
[polrNbd, polrErr] = normIndivid_e200(fdt.experimentData, nBrt, nDrk, nRep);

% Get mean after user defined starting point
meanNbd  = mean(polrNbd(avgStIdx:end));
stdeNbd  = std(polrNbd(avgStIdx:end))/sqrt(length(polrNbd(avgStIdx:end)));

% Fit to exponential
if pinMean
    aLw  = 1;
    aUp  = 1;
else
    aLw  =  .5;
    aUp  = 1.5;
end

ftTypebd = fittype( 'a - b.*exp(-x/c)', 'independent', {'x'}, 'dependent', {'y'},...
               'coefficients', {'a', 'b', 'c'});
optsbd            = fitoptions( 'Method', 'NonlinearLeastSquares' );
optsbd.Display    = 'Off';
optsbd.StartPoint = [   .5,   1,    10];
optsbd.Lower      = [  aLw,  -2,     0];
optsbd.Upper      = [  aUp,   2, 10000];

[tFitbd,   yFitbd]  = prepareCurveData(dt(Nignore+1:end), polrNbd(Nignore+1:end));
[fitResbd,  gofbd]  = fit(tFitbd, yFitbd, ftTypebd, optsbd);

tvec   = dt(Nignore+1):(dt(end) - dt(Nignore+1))/(Nsamp -1):dt(end);

% Get results
aFitbd   = fitResbd.a;
bFitbd   = fitResbd.b;
cFitbd   = fitResbd.c;
fitCrvbd = feval(fitResbd, tvec);

% Get uncertainty
intResbd    = confint(fitResbd, 0.66);
intRes_abd  = intResbd(:,1);
aErrbd      = abs(intRes_abd(2) - intRes_abd(1))/2;

intRes_bbd  = intResbd(:,2);
bErrbd      = abs(intRes_bbd(2) - intRes_bbd(1))/2;

intRes_cbd  = intResbd(:,3);
cErrbd      = abs(intRes_cbd(2) - intRes_cbd(1))/2;

% Time for given polarization using fit
tTrgtbd     = cFitbd*log(bFitbd/pErrTrgtbd);


%% Plotting
% String conversions
SampStr      = num2str(Samp, '%5.0f');
NumAvStr     = num2str(NumAv,'%3.0f');
StartCtsStr  = num2str(round(StartCts, 3, 'significant'));

SRS1poweStr  = num2str(SRS1powe,'%4.2f');
SRS1freqStr  = num2str(SRS1freq*1e-9,'%7.5f');

tAvgSrStr   = num2str(tAvgStrt,   '%.2f');

advStr    = num2str(aFitdv,  '%.4g');
bdvStr    = num2str(bFitdv,  '%.4g');
cdvStr    = num2str(cFitdv,  '%.4g');
advErrStr = num2str(aErrdv,  '%.2g');
bdvErrStr = num2str(bErrdv,  '%.2g');
cdvErrStr = num2str(cErrdv,  '%.2g');
meanNdvStr    = num2str(meanNdv,    '%.3f');
stdeNdvStr    = num2str(stdeNdv,    '%.3f');
pErrTrgtdvStr = num2str(pErrTrgtdv, '%.3f');
tTrgtdvStr    = num2str(tTrgtdv,    '%.3f');

abdStr    = num2str(aFitbd,  '%.4g');
bbdStr    = num2str(bFitbd,  '%.4g');
cbdStr    = num2str(cFitbd,  '%.4g');
abdErrStr = num2str(aErrbd,  '%.2g');
bbdErrStr = num2str(bErrbd,  '%.2g');
cbdErrStr = num2str(cErrbd,  '%.2g');
meanNbdStr    = num2str(meanNbd,    '%.3f');
stdeNbdStr    = num2str(stdeNbd,    '%.3f');
pErrTrgtbdStr = num2str(pErrTrgtbd, '%.3f');
tTrgtbdStr    = num2str(tTrgtbd,    '%.3f');

rts   = {
       ['Laser repolarize sweep.   Data file: ', fileNameStr];
       ['Seq: ', SeqName, ',  ','Analysis: ', AnalysisName]
       ['Start: ', StartTime, ',  ', 'Avgs: ', NumAvStr, ',  ', 'Samp/pt: ', SampStr,...
            ',  ', 'Start kcts: ', StartCtsStr];
       ['SRS1 Freq: ', SRS1freqStr, ' GHz,  ', 'SRS1 Power: ', SRS1poweStr]
         };

dvts  = {
         ['Model: a - b*exp(-x/c)   a = ' advStr, ' \pm ', advErrStr];
         ['b = ', bdvStr, ' \pm ', bdvErrStr, '  c = ', cdvStr, ' \pm ', cdvErrStr];
         ['Target error: ', pErrTrgtdvStr, '  required time: ', tTrgtdvStr];
         ['Averaging start: ', tAvgSrStr, ' Mean = ', meanNdvStr, ' \pm ', stdeNdvStr];
         };     

bdts  = {
         ['Model: a - b*exp(-x/c)   a = ' abdStr, ' \pm ', abdErrStr];
         ['b = ', bbdStr, ' \pm ', bbdErrStr, '  c = ', cbdStr, ' \pm ', cbdErrStr];
         ['Target error: ', pErrTrgtbdStr, '  required time: ', tTrgtbdStr];
         ['Averaging start: ', tAvgSrStr, ' Mean = ', meanNbdStr, ' \pm ', stdeNbdStr];
         };     
     
     
figure('pos', [50, 50, 1400, 750], 'color', [1,1,1])
subplot(2,2,1)
hold on 
grid on

xlabel('Duration (\mus)');
ylabel('Raw Counts');

plot(dt, brtAvg,    '.-r',                'MarkerSize', 13, 'Linewidth', 1.2)
plot(dt, drkAvg,    '.-b',                'MarkerSize', 13, 'Linewidth', 1.2)
plot(dt, darkVar,   '.-k', 'Color', purp, 'MarkerSize', 13, 'Linewidth', 1.2)
plot(dt, polrRaw,   '.-',  'Color', gref, 'MarkerSize', 13, 'Linewidth', 1.2) 

title(rts);
legend('Bright avg.', 'Dark', 'Repolarize', 'location', 'SouthEast');


subplot(2,2, 2)
hold on 
grid on

xlabel('Duration (\mus)');
ylabel('Polarize/Bright');

plot(dt,   polrNdv,           '.-',  'Color', gref, 'MarkerSize', 13, 'Linewidth', 1.2) 
plot(tvec, fitCrvdv,          '-',  'Color',  'k', 'MarkerSize', 13, 'Linewidth', 1.2)

title(dvts);
legend('Repolarize', 'Fit', 'location', 'SouthEast')


subplot(2,2, [3,4])
hold on 
grid on

xlabel('Duration (\mus)');
ylabel('Polarize/Bright');

errorbar(dt,   polrNbd,  polrErr, 'o',  'Color', 'k', 'MarkerSize', 6, 'Linewidth', 1.2, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', gref, 'CapSize', 0) 
plot(    tvec, fitCrvbd,          '-',  'Color', 'r', 'MarkerSize', 13, 'Linewidth', 1.2)

title(bdts);
legend('Repolarize', 'Fit', 'location', 'SouthEast')












