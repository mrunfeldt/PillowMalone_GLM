function [nev_vec,dmr,t_vec,faxis] = glmwrappor(dmrfile,nevfile)

% Load event and trigger data...
load(nevfile);
fs = trigstruct.samprate;
if isfield(trigstruct,'trigsinsamples'); % MJR mod: 2015_11_19
    trigger = trigstruct.trigsinsamples ;
elseif isfield(trigstruct,'trigsinms');
    trigger = trigstruct.trigsinms*1e-3*trigstruct.samprate;
else disp ('!!!! No trigger in nev file !!!!!'); return
end


%"faxis" and more are loaded from "_param.mat" file that must be saved in
%same foder as the dmrfile and have a corresponding name % MJRunfeldt
%10_09_2015
% Load dmr and creating time vector...
[dmr,t_vec,faxis] = load_ripple(dmrfile,'float',trigger,fs); %



% Create event vector...
bin_edges = [t_vec (t_vec(end)+mean(diff(t_vec)))];
nev_vec = histc(nevstruct.timesinms,bin_edges);
nev_vec = nev_vec(1:end-1);