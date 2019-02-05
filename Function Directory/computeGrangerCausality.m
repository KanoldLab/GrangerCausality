function [Devtotal GCJL  GCJH  wftotal  sig2ftotal  sig2rtotal  robstotal  rhatftotal  rhatrtotal  Bftotal  Brtotal ...
    cellids  gammaOpttotal XY Whc Whcv Whcs gammaOpt Jcost] = computeGrangerCausality(DFF, xy, Ncells, cellids, gammacv, ...
    alpha, Md, WH, Mhc, Mhcs, LL, LR,LLcv, LRcv, Mlow)

Devtotal = zeros(Ncells,Ncells,1);
DLLtotal = Devtotal;
wftotal = cell(Ncells,1);
robstotal = cell(Ncells,1);
sig2ftotal = zeros(Ncells,1);
sig2rtotal = zeros(Ncells,Ncells,1);
Bftotal = zeros(Ncells,1);
Brtotal = zeros(Ncells,Ncells,1);
rhatftotal = cell(Ncells,1);
rhatrtotal = cell(Ncells,Ncells,1);
gammaOpttotal = zeros(Ncells,1);

%Selected cells
XY = xy(cellids,:);
DFF = DFF(:,:,cellids);
[Nframes,Ntrials,Ncells] = size(DFF);

% Cross-Hist Kernel: make sure that the length latency of interaction + latency of Ca maging
Whc = [1 WH*ones(1,Mhc-1)];
Lhc = sum(Whc);
Whcv = {Whc};

% Self-Hist Kernel
Whcs = [1 WH*ones(1,Mhcs-1)];
Whscv = {Whcs};
Lhcs = sum(Whcs);
Mf = (Ncells-1)*Mhc+Mhcs;
Mr = (Ncells-2)*Mhc+Mhcs;

% Cross-Validation
[gammaOpt,Jcost] = CrossValidModTrial2(DFF,gammacv,{Whc},{Whcs},LLcv,LRcv);
gammaOpttotal(:,1) = gammaOpt;

% Form Cross-History Covariate Matrices
Lm = max([Lhc,Lhcs]);
Np = Nframes - Lm;
Xcz = zeros(Np,Mhc,Ncells,Ntrials);
for r = 1:Ntrials
    %size [Nframes x Ncells]
    Xcz(:,:,:,r) = FormHistMatrixv2(squeeze(DFF(:,r,:)),Whc,Lm);
end

% Form Observation arrays for all cells
Robs = DFF(Lm+1:end,:,:);

% Zero-mean the observation vectors. This will remove the need for intercept (mu) estimation
for cc = 1:Ncells
    Robs(:,:,cc) = Robs(:,:,cc) - ones(Np,1)*mean(Robs(:,:,cc));
end
for ct = 1:Ncells
    N_indto = num2str(ct);
    cellC = 1:Ncells; cellC(ct) = [];
    gamma = gammaOpt(ct);
    
    % Form Self-History Covriate matrix
    Xcsz = zeros(Np,Mhc,Ntrials);
    for r = 1:Ntrials
        % size [Nframes x Ncells]
        Xcsz(:,:,r) = FormHistMatrixv2(squeeze(DFF(:,r,ct)),Whcs,Lm);
    end
    
    % Form Effective response vector
    robsz = Robs(:,:,ct);
    reff = robsz(:); Neff = numel(reff);
    
    % Full Model
    % Form standardized Full Model Covariate Matrix
    [Xneff,Dnf] = FormFullDesignv3(Xcz,Xcsz,ct);
    
    % Filtering/Estimation Setting
    [wnf,sig2f,Devf,Bf] = StaticEstim(Xneff,reff,gamma,LL,LR);
    
    % Save Tuning AR Coeffs full model for each cell
    wftotal{ct,1} = wnf;
    sig2ftotal(ct,1) = sig2f;
    robstotal{ct,1} = robsz;
    
    % Reconstructed observation vector full model
    rhatf = zeros(Np,Ntrials);
    for r = 1:Ntrials
        Xnf = Xneff((r-1)*Np+1:r*Np,:);
        rhatf(:,r) = Xnf*wnf;
    end
    rhatftotal{ct,1} = rhatf;
    Bftotal(ct,1) = Bf;
    
    for cf = cellC
        N_indfrom = num2str(cf) ;
        disp(['Estimating G-Causality from cell ' N_indfrom ' to cell ' N_indto ' ...'])
        
        % Reduced Model
        Xnefr = Xneff;
        cc = find(cellC == cf);
        Xnefr(:,(cc-1)*Mhc+Mhcs+1:cc*Mhc+Mhcs) = [];
        
        % Estimation Setting
        [wnr,sig2r,Devr,Br] = StaticEstim(Xnefr,reff,gamma,LL,LR);
        
        % Save estimated params/arrays
        sig2rtotal(ct,cf,1) = sig2r;
        
        % Reconstructed observation vector reduced model
        rhatr = zeros(Np,Ntrials);
        for r = 1:Ntrials
            Xnr = Xnefr((r-1)*Np+1:r*Np,:);
            rhatr(:,r) = Xnr*wnr;
        end
        rhatrtotal{ct,cf,1} = rhatr;
        Brtotal(ct,cf,1) = Br;
        
        % Granger Causality metric: Compute Difference of Deviance Statistic
        Devd = Devf - Devr;
        DB = Bf - Br;
        DLL = Devd - DB;
        Devtotal(ct,cf,1) = Devd;
        DLLtotal(ct,cf,1) = DLL;
        
    end
end

% Statistical Test: J-statistics based on FDR
DkSmth = Devtotal-Md;
[GCJL,GCJH] = FDRcontrolBHv5(Devtotal,DkSmth,alpha,Md,wftotal,Whc,Mlow);
