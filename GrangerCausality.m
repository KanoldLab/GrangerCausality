function Output = GrangerCausality(DFF, xy, MaxCells)
%Original formulation and coding by Alireza Sheikhattar & Behtash Babadi
%Implemented in Francis et al (2018). Neuron. 

%%Input
%DFF: fluorescence matrix of frame X trial X cell
%xy: spatial coordinates of cells
%MaxCells: max number of cells in Granger analysis

%% GC Parameters
% Parameter estimation settings
LL = 1000;
LR = 10;
LLcv = 100;
LRcv = LR;

% Cross-history covariance settings
%Window length
WH = 5;
%Samples per cell in the model
Mhc = 5;

% Self-history covariance settings
%Samples per cell in the model
Mhcs = 5;

% check for 15 diff gamma for each cell
gammacv = [0:0.2:3];

% Statistical Test: J-statistics based on FDR
%bigger value test higher stat val.
Md = Mhc;
%sig. level;
alpha = 0.01;

%determines which segemetns of analysis window are used to determine inhibitory vs excitatory links:
%limit length of window to avoid phase ambiguoity
Mlow = 2;

%% Gather data
DFFtemp=[];
if iscell(DFF)
    for ff = 1:length(DFF)
        DFFtemp = cat(2,DFFtemp,DFF{ff});
    end
    DFF=DFFtemp;
end

% Locate and remove Nan cells from recordings
Rnan = isnan(DFF);
cnan = mean(mean(Rnan,1),2);
ccorr(1,:) = squeeze(cnan == 1);
cpass = find(sum(ccorr,1) == 0);

% Locate and remove Nan trials from recordings
Rnan = isnan(DFF(:,:,cpass));
inan = sum(sum(Rnan,1),3);
rcorr = find(inan == 0);
DFF = DFF(:,rcorr,:);
xy = xy(cpass,:);
Ncells = min(size(DFF,3),MaxCells);

% Select a Subset of cells with highest response variability
varR = zeros(size(DFF,3),1);
if MaxCells < size(DFF,3)
    RR = size(DFF,2);
    for c = 1:size(DFF,3)
        varR(c) = mean(var(DFF(:,:,c)));
    end
    varS = (varR*RR)/sum(RR);
    [~,indm] = sort(varS,'descend');
    cellids = sort(indm(1:MaxCells)) ;
else
    cellids = [1:Ncells]';
end

%% GC Analysis: Multi-trial Multi-Cell Mod
[Devtotal GCJL  GCJH  wftotal  sig2ftotal  sig2rtotal  robstotal  rhatftotal  rhatrtotal  Bftotal  Brtotal ...
    cellids  gammaOpttotal XY Whc Whcv Whcs gammaOpt Jcost] = computeGrangerCausality(DFF, xy, Ncells, cellids, gammacv, ...
    alpha, Md, WH, Mhc, Mhcs, LL, LR,LLcv, LRcv, Mlow);

%% Output
%Paramaters
Output.Params.MaxCells = MaxCells;
Output.Params.Whc = Whc;
Output.Params.Whcs = Whcs;
Output.Params.LL = LL;
Output.Params.LR = LR;
Output.Params.WH = WH;
Output.Params.Mhc = Mhc;
Output.Params.Mhcs = Mhcs;
Output.Params.gammacv = gammacv;
Output.Params.Md = Md;
Output.Params.alpha = alpha;
Output.Params.Mlow = Mlow;

%Cross-validation
Output.CrossVal.gammaOpt = gammaOpt;
Output.CrossVal.Jcost = Jcost;
Output.CrossVal.Whcv = Whcv;
Output.CrossVal.LLcv = LLcv;
Output.CrossVal.LRcv = LRcv;

%GC analysis output
Output.Devtotal = Devtotal;
Output.GCJL = GCJL;
Output.GCJH = GCJH;
Output.wftotal = wftotal;
Output.sig2ftotal = sig2ftotal;
Output.sig2rtotal = sig2rtotal;
Output.robstotal = robstotal;
Output.rhatftotal = rhatftotal;
Output.rhatrtotal = rhatrtotal;
Output.Bftotal = Bftotal;
Output.Brtotal = Brtotal;
Output.cellids = cellids;
Output.gammaOpttotal = gammaOpttotal;