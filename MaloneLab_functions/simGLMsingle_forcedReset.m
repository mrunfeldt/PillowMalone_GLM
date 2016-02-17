
% Modified by MJRunfeldt; 2015_12_21
% When more than

% ======================================================
function [tsp,Vmem,Ispk] = simGLMsingle_forcedReset(glmprs,Stim)
% Sub-function for simGLM.m
% Simulates the GLM point process model for a single (uncoupled) neuron
% using time-rescaling

% --------------- Check Inputs ----------------------------------
global RefreshRate;

dt = glmprs.dt;
nbinsPerEval = 100;  % Default number of bins to update for each spike

% % % When FR in window length(pastWindow) > FRmax, the Vm is "reset to
% what it would be if there was no spike-history filter % % % % 
pastWindow = 25 / dt ; % (ms) => samples
FRmax = 200; % in Hz
nSpikeMax = round(FRmax * (pastWindow*dt*1e-3) ); 

% Check that dt evenly divides 1
if (mod(1,dt) ~= 0), dt = 1/round(1/dt);
    fprintf(1, 'glmrunmod: reset dtsim = %.3f\n', dt);
end
if ~isfield(glmprs, 'nlfun'), nlfun = @exp;
else nlfun = glmprs.nlfun;
end

[slen,swid] = size(Stim); % length of stimulus
rlen = round(slen/dt);  % length of binned response

% -------------  Compute filtered resp to signal ----------------
k = glmprs.k;  % stim filter
if swid == size(k,2) % convolve if k is the same width as stimulus
    Vstm = sameconv(Stim,k); 
else 
    error('Mismatch between stimulus and filter size -- glmrunmod');
end
Vstm = [0; Vstm; 0]+glmprs.dc;

% -------------- Compute interpolated h current ----------------------
if isfield(glmprs,'ihbas');  % Check if h-kernel parametrized by a basis
    if isempty(glmprs.ih), ih = [];
    else ih = glmprs.ihbas*glmprs.ih;
    end
else ih = glmprs.ih;
end
if ~isempty(ih) % Interpolate ih current
    ihthi = [dt:dt:max(glmprs.iht)]';  % time points for sampling
    ihhi = interp1(glmprs.iht, ih, ihthi, 'linear', 0);
    hlen = length(ihhi);
else % No ih current
    hlen = 1; ihhi = 0; ihthi = dt;
    if ~isempty(glmprs.iht)
        fprintf('Warning: no post-spike kernel supplied but iht nonempty!!\n');
    end
end
    
% -------------  Static nonlinearity & spiking -------------------
nsp = 0;
tsp = zeros(round(slen/5),1);  % allocate space for spike times
Vmem = interp1([0:slen+1]',Vstm,[.5+dt:dt:slen+.5]', 'linear');
VmOG = Vmem; % MJR; Vm before spike-history filter

if (nargout > 2), Ispk = Vmem*0;
end

jbin = 1; % current time bin
nsp = 0; % number of spikes
tspnext = exprnd(1);  % time of next spike (in rescaled time); Random number from exponential distribution
rprev = 0;  % Integrated rescaled time up to current point
while jbin <= rlen
    iinxt = jbin:min(jbin+nbinsPerEval-1,rlen);
    rrnxt = nlfun(Vmem(iinxt))*dt/RefreshRate; % Cond Intensity
    rrcum = cumsum(rrnxt)+rprev; % integrated cond intensity
    
%     figure;subplot(3,1,1); plot(Vmem(iinxt),'linewidth',3);title('Vmem(iinxt)');
%     subplot(3,1,2);plot(rrnxt,'linewidth',3);title('rrnxt, "Cond Intensity"');
%     subplot(3,1,3);plot(rrcum,'linewidth',3);title('rrcum, "Integrated Cond Intensity"');pause; close 
 
    
    if (tspnext >= rrcum(end)) % No spike in this window
        jbin = iinxt(end)+1;
        rprev = rrcum(end);
    else   % Spike!
        ispk = iinxt(min(find(rrcum>=tspnext))); % time bin where spike occurred
        nsp = nsp+1;
        tsp(nsp,1) = ispk*dt; % spike time
        mxi = min(rlen, ispk+hlen); % max time affected by post-spike kernel
        iiPostSpk = ispk+1:mxi; % time bins affected by post-spike kernel
        
% % THIS SECTION modified by MJR % % %        
        pint = max(1,jbin-pastWindow) *dt;
        nSpikes = sum(tsp>pint & tsp<=ispk*dt); % number of spikes in window leadin
        if ~isempty(iiPostSpk) && nSpikes > nSpikeMax
            Vmem(iiPostSpk,:) = VmOG(iiPostSpk,:) ; % RESET 
        elseif ~isempty(iiPostSpk) && nSpikes <= nSpikeMax
            Vmem(iiPostSpk) = Vmem(iiPostSpk)+ihhi(1:mxi-ispk);
        end
        
        tspnext = exprnd(1);  % draw next spike time
        rprev = 0; % reset integrated intensity
        jbin = ispk+1;  % Move to next bin
        % --  Update # of samples per iter ---
        muISI = jbin/nsp;
        nbinsPerEval = max(20, round(1.5*muISI)); 
    end
end
tsp = tsp(1:nsp); % prune extra zeros
