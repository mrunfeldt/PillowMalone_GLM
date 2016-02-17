
% "testStim" - spectrogram of test stimulus
% "testResp" - spike times of actual neuron in response to test stimulus.
    % * Each cell is trial
% Last modified 2016_02_16 by MJRunfeldt
% % % % % Run script % % % % 

clear all
% Set some directories we need (CHANGE THESE!)
GLMdir = '/Users/mel/Desktop/GLM/devel_GLMspiketools/'; % home directory for GLM code
%WorkingDir = pwd;  % current directory

DataParent = '/Users/mel/Desktop/GLM/dataFor_JPlillow/Nrn_1_crossClass/' ; % directory where data is

% % % Specify which cell to load
dmrfile = [DataParent,'dmr-50flo-40000fhi-4SM-150TM-40db-96kHz-96DF-30min.spr']; % Stimulus file.
nevfile = [DataParent,'nev_2009-07-08_13-13_ch16_curly_dmrlong.mat'];  % Spikes 

% % Add paths
% addpath(DataDir);  % Add data directory
cd(GLMdir); %setpaths; cd(WorkingDir); % this will add paths for GLM code

% % Set some global vars we will need
global RefreshRate;  % Stimulus refresh rate (Stim frames per second)
RefreshRate = 1000; % Stimuli are binned at 1ms
DTsim = .1; % bin size for computing log-likelihood & simulations (must evenly divide 1)

% % Load data % % % 
% Execute the wrapping function...
% "faxis" is the frequency axis for dmr: it comes from the "_param.mat"
% file that is loaded by "load_ripple" in "glmwrappor"
[nev_vec,dmr,t_vec,faxis_orig] = glmwrappor(dmrfile,nevfile); 

% MJR: Truncate Frequency axis to match sampling rate of future test stimulus % % % 
axLim=find(faxis_orig>1e3 & faxis_orig < 24e3) ; % truncate FREQUENCY AXIS
faxis = faxis_orig(axLim);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

dmr = dmr(axLim,:); 
nev_vec = nev_vec';dmr = dmr'-mean(dmr(:)); % normalize DMR

% % Rename key items we need
stim = dmr; clear dmr % dmr = stimulus. (time x freq)
tshift = 5; % number of bins by which to shift spikes back in time (due to latency of neuron).
spcount = [nev_vec(tshift+1:end);zeros(tshift,1)]; % spike count (binned at 1ms)
[slen,nkx] = size(stim); 

% % Compute the STA (i.e., STRF)
nkt = 64; % number of time lags in STA
sta = simpleSTC(stim,spcount,nkt);  % compute STA

% fAxis_ax = [1:10*round(length(faxis)/100):length(faxis)] ; % for freq axis label
% figure;imagesc(1:nkx,1:nkt,sta); set(gca,'xtick',fAxis_ax,'xticklabel',num2cell(round(faxis(fAxis_ax))))
% ylabel('Time before spike (ms)'); xlabel('Frequency [space] (Hz)'); title('STA = STRF')

% % % % Let's restrict stimulus frequency region based on STA
% MJR - let's not for now
% % iix = 81:125; % indices (frequency bins of spectrogram) to keep
% % stimShort = stim(:,iix); % reduced stimulus
% % nkx = size(stimShort,2);
% % sta = simpleSTC(stimShort,spcount,nkt);
% % figure; imagesc(1:nkx,1:nkt,sta);

% % Fit GLM filter in basis spanned by first few spatial singular vectors of STA
[u,s,v] = svd(sta); % singular value decomposition of STA

nbasis = 4; % number of singular vectors to keep
stabasis = v(:,1:nbasis);
stimproj = stim*stabasis; 
% MJR: stimulus (long time x freq) is reporojected into new frequency
% sub-space

% % Compute projected STA (just for fun)
staproj = simpleSTC(stimproj,spcount,nkt);
% figure;subplot(2,1,1);plot(diag(s),'o-'); title('singular Values - diag(s)')
% subplot(2,1,2);imagesc(staproj); title(['Reprojected stimulus - ',num2str(nbasis),'D'])

% % Now let's set up params for GLM in this reduced spatial stimulus
tsp = find(spcount);  % binarize spike counts to 0 / 1

%  Initialize params for fitting --------------
gg0 = makeFittingStruct_GLM(staproj,DTsim); % input is STA - stimulus filter (reduced in dimensionality)
gg0.tsp = tsp;
gg0.tspi = 1;
% remove spike-history filter
gg0.ih = [];
gg0.iht = [];
gg0.ihbas = [];
gg0.ihbasprs = [];

% % Now set basis for stimulus filter
kbasprs.neye = 0; % number of "identity" basis vectors
kbasprs.ncos = 16; % number of "raised cosine" basis vectors
kbasprs.kpeaks = 3*[1 (kbasprs.ncos)]; % location of first and last cosine basis vectors (relative to identity ones).
kbasprs.b = 1e4;  % set degree of nonlinear scaling of basis (higher -> more linear); must be positive.
[ktbasisorth,ktbasis] = makeBasis_StimKernel(kbasprs,nkt);

% Examine basis and STA
% figure;subplot(211); ttk = -nkt+1:0;
% plot(ttk,ktbasis); title('temporal basis for stimulus filter');
% subplot(212); plot(ttk,staproj, 'k',ttk, ktbasisorth*(ktbasisorth'*staproj), 'r--');
% title('STA (black) and best fit in basis (red)');xlabel('time before spike (ms)');

gg0.ktbas = ktbasisorth; % set basis for stim filter
gg0.kt = (ktbasisorth'*staproj); % set initial stim filter in basis
gg0.k = ktbasisorth*gg0.kt; % set stimulus filter itself
gg0.kbasprs = kbasprs; % set basis params for stim filter

% % Do ML fitting of filter, using full GLM with projected stimulus

%[logli0,rr0,tt0] = neglogli_GLM(gg0,stimproj); % Compute logli of initial params (if desired)
[logli0] = neglogli_GLM(gg0,stimproj); % we don't need rr0 or tt0 - MJR
logli0

% Do ML estimation of model params
opts = {'display', 'iter', 'maxiter', 100};
[gg1, negloglival1] = MLfit_GLM(gg0,stimproj,opts); % find ML estimate (requires optimization toolbox)

kfit1 = gg1.k*stabasis';  % estimated filter

% % % Now add spike-history filter
% % Compute spike-triggered spike history
binarycounts = min(spcount,1);  % truncate spike counts at above 1 to 1.
nsp = sum(binarycounts); % number of binary spikes.
tbins = 100; % number of time bins of spike history to use
yysta = simpleSTC(binarycounts(1:end-1)-mean(binarycounts),binarycounts(2:end),tbins); % create vectors of spike history
yysta = yysta-yysta(1);

% % Set up basis for post-spike filter
ihbasprs.ncols = 5;  % Number of basis vectors for post-spike kernel
ihbasprs.hpeaks = [.1 12];  % Peak location for first and last vectors
ihbasprs.b = .25;  % How nonlinear to make spacings
ihbasprs.absref = []; % absolute refractory period 
[iht,ihbasorth,ihbasis] = makeBasis_PostSpike(ihbasprs,1);
iht(1) = []; ihbasorth(1,:) = []; ihbasis(1,:) = []; % remove 0-bin for comparison
% 
% % % % Now inspect whether the basis can capture shape of spike-history sta
% % yystaIntrp = interp1((1:tbins)',flipud(yysta),iht,'linear',0);
% % figure; subplot(211); plot(iht,ihbasis);
% % subplot(212); plot(iht, iht*0, 'k', 1:tbins,flipud(yysta),iht, ihbasorth*(ihbasorth'*yystaIntrp));
% 
% % % % Make basis and insert them into fitting struct
[iht,ihbasorth,ihbasis] = makeBasis_PostSpike(ihbasprs,DTsim);
gg1.iht = iht;
gg1.ihbas = ihbasorth;
gg1.ihbasprs = ihbasprs;
gg1.ih = zeros(size(ihbasorth,2),1);


% % Do ML estimation of model params
[ggM, negloglivalSp] = MLfit_GLM(gg1,stimproj,opts); % find ML estimate


%% % % Now use Model to predict response to novel stimulus (Monkey call) % %

% % LOAD Call (test stimulus) spectrogram and spikes % % 
% spectrogram is processed to have spatial and temporal resolution as DMR
load([DataParent,'testStim'])

% % Uncomment the code below if you want to verify that the frequecy axis for training
% stimulus matches that of the test stimulus% % %
% dummy = load([DataParent,'faxis']); 
% figure;hold on; plot(faxis,'k-x','linewidth',2.5);plot(dummy.faxis,'r-x')


% % % Plot test stimulus spectrogram % % 
figure;pcolor([1:size(testStim,1)].*(1/RefreshRate),faxis,testStim');colormap jet; shading flat;colorbar;axis tight; 
xx=xlabel('Time (sec)');yy=ylabel('Frequency (Hz)');set(xx,'fontsize',18');set(yy,'fontsize',18')
tt=title('Test Stimulus = SQM Call');set(tt,'fontsize',18)

% Reproject Test stimulus into same dimensional space as training stimulus
newTest = testStim*stabasis;

% % LOAD Spikes (test response) 
load([DataParent,'testResp']) % "nrn" = sorted; "mua" is multiunit
nTrials = length(testResp);

spPred = cell(1,nTrials) ; 
for k = 1:nTrials % simulate trials
    %[spPred{k},vMpred] = simGLM(ggM,newTest); % PREDICT!
    [spPred{k},vMpred] = simGLMsingle_forcedReset(ggM,newTest);
     
% % % % Plot Per Trial Prediction % % % % % % % % % % % % % % % % % % % % % % % % % % 
%     figure;plot([1:length(vMpred)].*((1/RefreshRate)*DTsim),vMpred,'r','linewidth',1.5);
%     xx=xlabel('Time (sec)'); yy=ylabel('Vm');axis tight; tt=title('Prediction: press any key for next trial');
%     set(xx,'fontsize',18);set(yy,'fontsize',18);set(tt,'fontsize',18);pause;close
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
end

% % % Compare prediction to data % % 
time = [1:length(vMpred)*DTsim].*(1/RefreshRate) ; 
hData = zeros(nTrials,length(time)); hPred = hData ; % initialize

psthR = zeros(1,nTrials);
for a = 1:nTrials
    pred = spPred{a}.*(1/RefreshRate) ;  dat = testResp{a} ;
    hPred(a,:) = histc(pred,time).*RefreshRate; hData(a,:) = histc(dat,time).*RefreshRate ; % Bin
    dum=corrcoef(hPred(a,:),hData(a,:)); psthR(a) = dum(2);
end
figure;hist(psthR,20);leg=legend(['Mean R = ',num2str(mean(psthR))]); set(leg,'fontsize',18)
tt=title('PSTH Correlation Per Trial'); set(tt,'fontsize',18); 

fH = figure; 
subplot(2,1,1);hold on; plot(time,mean(hPred),'m','linewidth',3);
plot(time,mean(hPred)+std(hPred),'m','linewidth',0.5);plot(time,mean(hPred)-std(hPred),'m','linewidth',0.5)
tt=title('Prediction, Mean +/- Std over trials'); set(tt,'fontsize',18); axis tight
yy=ylabel('Instantaneous FR (Hz)'); set(yy,'fontsize',18); xx=xlabel('Time (sec)'); set(xx,'fontsize',18)
subplot(2,1,2);hold on;  plot(time,mean(hData),'b','linewidth',3);
plot(time,mean(hData)+std(hData),'b','linewidth',0.5);plot(time,mean(hData)-std(hData),'b','linewidth',0.5)
tt=title('Data, Mean +/- Std over trials'); set(tt,'fontsize',18); axis tight
yy=ylabel('Instantaneous FR (Hz)'); set(yy,'fontsize',18); xx=xlabel('Time (sec)'); set(xx,'fontsize',18)
set(fH, 'Position', [0 0 11e2 7e2]);movegui(fH,'center');
