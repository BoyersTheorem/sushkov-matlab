%% Initialize, load data, user input here

clear
fileName = 'newSetup_Qtrigger_delay_1';
load(fileName);

Tpi        = 40;   % Time of pi pulse in ns
Ppi        = 1.97; % Power in AMPR to achieve Tpi above. Should be for on resosnace as that Ppi depends only on drive strength

qVoltage   = 265.2/1000;  % Volts

%% Extract parameters and data from files
StartT          = experimentData.MetaData.StartDateTime;     % scan start time
StartTime       = datestr(StartT,'yyyy-mm-dd HH:MM:SS');     % scan start time string
StartDay        = datestr(StartT,'yyyy-mm-dd');              % scan start day
Samp            = experimentData.MetaData.Samples;           % samples per point
StartCts        = experimentData.MetaData.InitialCounts/1e3; % starting counts
SeqName         = otherParam.SEQ.name;                       % Seq name
SeqName         = strrep(SeqName, '_','\_');                 % Seq modified for figure
AnalysisName    = 'MW_delays_e100';                 % This script
AnalysisName    = strrep(AnalysisName, '_', '\_');
fileNameStr     = strrep(fileName, '_', '\_');

SRS1powe         = experimentData.Parameters.MWPower;         %
SRS1freq         = experimentData.Parameters.MWFreq;

piStrength = 1/(4*Tpi)*1e3; % Strength of pi pulse in MHz units

% Determine actual number of averages.
MaxAv = experimentData.MetaData.Average; % max number of avgs
for j=1:MaxAv
    if experimentData.Data.AVE(1,j).X(1,1).xmean(1) == 0
        break;
    end
end
NumAv = j-1;

% Extract Data
x         = DataMat(1,:)*1e9; % Time in ns
bright    = DataMat(2,:);
dark      = DataMat(3,:);
const     = DataMat(4,:);
sweep     = DataMat(5,:);

%% Analysis
% Which APD reads are the bright and dark
nBrt = 1;
nDrk = 2;

% Normalize data
[shrtDat, shrtStd] = normIndivid(experimentData, nBrt, nDrk, 3);
[longDat, longStd] = normIndivid(experimentData, nBrt, nDrk, 4);

%% Fitting 
% Model: stays constant, then linear rise, then stays constant
% Create model
modelDelay = @(x, xd) heaviside(xd - x(1)).*heaviside(x(2) - xd).*(1 - 0)./(x(2) - x(1)).*(xd - x(1)) + ...
                      heaviside(xd - x(2)).*1;

% Perform fitting to data - start and bounds
x0 = [ -20,  +40];
lb = [ -60,    0];
ub = [ +60,  100];

% Do a fitting with non-linear curve fit model. TBH, don't know why I used
% this instead of 'fit' built-in func
xresult = lsqcurvefit(modelDelay, x0, x, longDat, lb, ub);

% Sample the solution
Nsamp = 10000;
xsmp = min(x):(max(x)-min(x))/(Nsamp-1):max(x);
yfit  = modelDelay(xresult, xsmp);
%% Plot
% Custom Colors
purp = [ .8,  .2, .8]; % Richer purple
gref = [ .2,  .9, .2]; % 'Fuller' green - more saturated
cyaf = [ .2,  .8, .8]; % 'Fuller' cyan  - more saturated
oraf = [.95, .65, .1]; % fuller orange

% Convert numbers to strings for title
SampStr      = num2str(Samp, '%5.0f');
NumAvStr     = num2str(NumAv,'%3.0f');
StartCtsStr  = num2str(round(StartCts, 3, 'significant'));

SRS1poweStr  = num2str(SRS1powe,'%4.2f');
SRS1freqStr  = num2str(SRS1freq*1e-9,'%7.5f');
SRS1_40nsPi  = num2str(Ppi, '%3.2f');

x1str        = num2str(xresult(1), '%.3g');
x2str        = num2str(xresult(2), '%.3g');

% Create Raw Title String (rts)
rts = {
       ['Microwave Delay,  ', fileNameStr];
       ['Seq: ', SeqName, ',  ','Analysis: ', AnalysisName]
       ['Start: ', StartTime, ',  ', 'Avgs: ', NumAvStr, ',  ', 'Samp/pt: ', SampStr,...
            ',  ', 'Start kcts: ', StartCtsStr];
       ['SRS1 Freq: ', SRS1freqStr, ' GHz,  ', 'SRS1 Power: ', SRS1poweStr]
       
       };

figure('pos', [100, 100, 1200, 600])
subplot(1,2,1)
hold on
grid on
xlabel('Time Delay')
ylabel('Raw counts')
title(rts)

plot(x, bright, '.-r',                'MarkerSize', 18, 'Linewidth', 1.2)
plot(x, dark,   '.-b',                'MarkerSize', 18, 'Linewidth', 1.2)
plot(x, const,  '.-',  'Color', purp, 'MarkerSize', 18, 'Linewidth', 1.2)
plot(x, sweep,  '.-',  'Color', gref, 'MarkerSize', 18, 'Linewidth', 1.2)

legend('Bright', 'Dark', 'Short', 'Long', 'location', 'SouthEast')
hold off   

ats = {
       ['x1 = ', x1str, '   x2 = ', x2str];
       
      };


subplot(1,2,2)
hold on
grid on
xlabel('Time Delay')
ylabel('Normalized')
title(ats)

plot(   x, shrtDat, '.-',  'Color', purp, 'MarkerSize', 18, 'Linewidth', 1.2)
plot(   x, longDat, '.-',  'Color', gref, 'MarkerSize', 18, 'Linewidth', 1.2)
plot(xsmp,    yfit, '-k',                 'MarkerSize', 18, 'Linewidth', 2)

%axis([min(x),max(x),.8,1.05])

legend('Short', 'Long', 'Fit', 'location', 'SouthEast')
hold off   



