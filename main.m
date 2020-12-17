clear
close
%% preprocessing

load('potenziali.mat')
locs = 'layout.xyz';
fs = 250;
[n_conds, n_subs, n_chs, n_samples] = size(ave_tot);
t_vect = linspace(-0.5, 1, n_samples);
interp_data = ave_tot;
for i = 1: n_conds
    for j = 1: n_subs
        figure
        plot(t_vect, squeeze(ave_tot(i, j, :, :))')
        title(['cond:', num2str(i), '<---->', 'subj:', num2str(j)])
        pause
        close
        answer = inputdlg('bad cjs');
        badchs = str2num(answer{1,1}); 
        if ~isempty(badchs)
            ch = 1:28;
            ch(badchs) = [];
            outData = itab_spline_interpolation(squeeze(ave_tot(i, j, :, :)), ch, badchs, locs);
            interp_data(i, j, badchs, :) = outData;
            figure
            plot(t_vect, squeeze(interp_data(i, j, :, :))')
            title(['cond:', num2str(i), '<---->', 'subj:', num2str(j)])
        end
        pause
        close
        BADC{i, j} = badchs;  
    end
end

%% clustering
clear
load('potenziali_prepro.mat')
load('globaltemplates.mat')
[n_conds, n_subs, n_chs, n_samples] = size(interp_data);
t_vect = linspace(-0.5, 1, n_samples);
time = t_vect>-0.1 & t_vect<0.5;
gavg = NaN(n_conds, n_chs, sum(time));
for i = 1: n_conds
    figure
    for j = 1: n_subs
        subplot(9, 2, j)
        plot(t_vect(time), squeeze(interp_data(i, j, :, time))')
        title(num2str(j))
    end
end
for i = 1: n_conds
    gavg(i, :, :) = squeeze(mean(interp_data(i, :, :, time), 2)) ;
end


%% single subject backfitting
str.MS.segment_min_len = 1;
str.MS.showpeaks = 0;
str.MS.use_gfp_peaks = 0;
settings.plottopos = 0;
t_vect = linspace(-0.1, 0.5, 150);

DUR = [];
OCC = [];
GEV = [];
COV = [];
nch = 28;
for c = 1: 3
    fig=figure; 
%      maps = squeeze( condWiseProt(c,:,:));
%     maps = meanCindWise;
      maps = orderedMaps;
for i = 1:18
    eeg = squeeze(interp_data(c, i, :, :))';
    eeg1 = eeg;
%     for k = 1: 5
%         time1 = t_vect>intervalli(k, 1) & t_vect<intervalli(k, 2)
%         eeg1 = eeg(time1,:)
        str.MS.tf = size(eeg1,1);
        [ output ] = itab_computemsparameters(eeg1, maps, nch, str.MS, 0);
        DUR(c, i, :) = output.mean_dur;
        OCC(c, i, :) = output.freq;
         GEV(c, i, :) = output.gev;
          COV(c, i, :) = output.cov;
%         subplot(18,1, i)
%         %     itab_plotSegments(eeg', maps', EEG.chanlocs, output.states_sequence  , 250,  settings, t_vect(time)*1000)
%         plot_gfp_area(eeg',maps', EEG.chanlocs, output.states_sequence, 250,  settings, t_vect(time)*1000)
%         vline([0 300])
%     end
end
end



%% non parametric comparison
clear
load('potenziali_prepro.mat')
load('globaltemplates.mat')
[n_conds, n_subs, n_chs, n_samples] = size(interp_data);
t_vect = linspace(-0.5, 1, n_samples);
time =t_vect>=-0.5 & t_vect<=1
gavg = NaN(n_conds, n_chs, sum(time));
for i = 1: n_conds
    gavg(i, :, :) = squeeze(mean(interp_data(i, :, :, time), 2)) ;
end

EEG.chanlocs = readlocs('Layout_pre_stimulus_ok.xyz');
for i = 1: 4
   figure
   topoplot(orderedMaps(i,:), EEG.chanlocs, 'style', 'map',  'electrodes', 'off') 
end

str.pnts = 150;
str.srate = 250;
[ str ] = itab_createmsstr (str);
load('potenziali_prepro.mat')
load('globaltemplates.mat')
t_vect = linspace(-0.5, 1, size(interp_data, 4));
time = t_vect>-0.1 & t_vect<0.5;
data = interp_data(:, :, :, time);
t_vect_new = t_vect(time);
[ncond, nsub, nchan, nsample]  = size(data);
datar = reshape(data,ncond*nsub, nchan, nsample);
maps = orderedMaps;
time1 = t_vect_new>0 & t_vect_new<0.5;

for iter = 1: 1000
    ind = randperm(ncond*nsub);
    c1 = datar(ind(1:18), :, time1);
    c2 = datar(ind(19:36), :, time1);
    c3 = datar(ind(37:54), :, time1);
    gc1 = squeeze(mean(c1, 1))';
    gc2 = squeeze(mean(c2, 1))';
    gc3 = squeeze(mean(c3, 1))';

    eeg1 = gc1;
    eeg2 = gc2;
    eeg3 = gc3;
    str.MS.tf = size(eeg1,1);
    nch = size(eeg1,2);
    [ output1 ] = itab_computemsparameters(eeg1, maps, nch, str.MS, 0);
     [ output2 ] = itab_computemsparameters(eeg2, maps, nch, str.MS, 0);
      [ output3 ] = itab_computemsparameters(eeg3, maps, nch, str.MS, 0);
    dur(iter,:) = [output1.mean_dur'-output2.mean_dur',output1.mean_dur'-output3.mean_dur',output2.mean_dur'-output3.mean_dur'];
    occ(iter,:) = [output1.freq'-output2.freq',output1.freq'-output3.freq',output2.freq'-output3.freq'];
    gev(iter,:) = [output1.gev-output2.gev,output1.gev-output3.gev,output2.gev-output3.gev];
    cov(iter,:) = [output1.cov'-output2.cov',output1.cov'-output3.cov',output2.cov'-output3.cov'];
end




%% real data
% time2 = t_vect_new>0 & t_vect_new<0.4;
time2=t_vect>-0.5 & t_vect<=1
for i = 1:3
   
%     eeg = squeeze(gavg(i,:,time2))';
    eeg = squeeze(interp_data(i,16,:,time2))';
%     maps = squeeze( condWiseProt(i,:,:));
%     maps = meanCindWise;
    maps = orderedMaps;
    nch = 28;
    interp_correlation = 1;
    str.MS.segment_min_len = 5;
    str.MS.tf = size(eeg, 1);
    confstruct = str;
    [ output ] = itab_computemsparameters(eeg, maps, nch, str.MS, 1);
    dur_real(i,:) = output.mean_dur';
    occ_real(i,:) = output.freq';
    gev_real(i,:) = output.gev;
    cov_real(i,:) = output.cov';
    
    settings = struct;
    settings.plot_time = [];
    settings.plottopos = 0;
    settings.plotsegnos = 'first';
    figure
    itab_plotSegments(eeg', maps', EEG.chanlocs, output.states_sequence, 250,  settings, t_vect(time2)*1000)
   
    set(gca, 'FontSize', 30)
end

%% single subjet plot area
for i = 1:3
   
%     eeg = squeeze(gavg(i,:,time2))';
    eeg = squeeze(interp_data(i,16,:,time2))';
%     maps = squeeze( condWiseProt(i,:,:));
%     maps = meanCindWise;
    maps = orderedMaps;
    nch = 28;
    interp_correlation = 1;
    str.MS.segment_min_len = 5;
    str.MS.tf = size(eeg, 1);
    confstruct = str;
    [ output ] = itab_computemsparameters(eeg, maps, nch, str.MS, 1);
    dur_real(i,:) = output.mean_dur';
    occ_real(i,:) = output.freq';
    gev_real(i,:) = output.gev;
    cov_real(i,:) = output.cov';
    
    settings = struct;
    settings.plot_time = [];
    settings.plottopos = 0;
    settings.plotsegnos = 'first';
    figure
    itab_plotSegments(eeg', maps', EEG.chanlocs, output.states_sequence, 250,  settings, t_vect(time2)*1000)
   
    set(gca, 'FontSize', 30)
end
for j = 1: 18
    figure
    for i = 1:3
        
        %     eeg = squeeze(gavg(i,:,time2))';
        eeg = squeeze(interp_data(i,j,:,time2))';
        %     maps = squeeze( condWiseProt(i,:,:));
        %     maps = meanCindWise;
        maps = orderedMaps;
        nch = 28;
        interp_correlation = 1;
        str.MS.segment_min_len = 5;
        str.MS.tf = size(eeg, 1);
        confstruct = str;
        [ output ] = itab_computemsparameters(eeg, maps, nch, str.MS, 1);
        dur_real(i,:) = output.mean_dur';
        occ_real(i,:) = output.freq';
        gev_real(i,:) = output.gev;
        cov_real(i,:) = output.cov';
        
        settings = struct;
        settings.plot_time = [];
        settings.plottopos = 0;
        settings.plotsegnos = 'first';
        
        subplot(3, 1, i)
        itab_plotSegments(eeg', maps', EEG.chanlocs, output.states_sequence, 250,  settings, t_vect(time2)*1000)
        
        set(gca, 'FontSize', 30)
    end
end
%%
% p values
% duration
dur_real(isnan(dur_real))=0;
dur(isnan(dur))=0;
dA = dur(:, 1:4:end);
dB = dur(:, 2:4:end);
dC = dur(:, 3:4:end);
dD = dur(:, 4:4:end);

dAGvsIPS = dur_real(1,:)-dur_real(2,:);
dAGvsSHAM = dur_real(1,:)-dur_real(3,:);
dIPSvsSHAM = dur_real(2,:)-dur_real(3,:);

pd = find(dA(:,1)>dIPSvsSHAM(1))
length(pd)/1000

% duration AG vs IPS : A 

% frequency
occ_real(isnan(occ_real))=0;
occ(isnan(occ))=0;
fA = occ(:, 1:4:end);
fB = occ(:, 2:4:end);
fC = occ(:, 3:4:end);
fD = occ(:, 4:4:end);

fAGvsIPS = occ_real(1,:)-occ_real(2,:);
fAGvsSHAM = occ_real(1,:)-occ_real(3,:);
fIPSvsSHAM = occ_real(2,:)-occ_real(3,:);

pf = find(fC(:,2)>abs(fAGvsSHAM(3)))
length(pf)/1000

% variance
gev_real(isnan(gev_real))=0;
gev(isnan(gev))=0;
vA = gev(:, 1:4:end);
vB = gev(:, 2:4:end);
vC = gev(:, 3:4:end);
vD = gev(:, 4:4:end);

vAGvsIPS = gev_real(1,:)-gev_real(2,:);
vAGvsSHAM = gev_real(1,:)-gev_real(3,:);
vIPSvsSHAM = gev_real(2,:)-gev_real(3,:);

pv = find(vC(:,1)>vAGvsSHAM(3))
length(pv)/1000

% coverage
cov_real(isnan(cov_real))=0;
cov(isnan(cov))=0;
cA = cov(:, 1:4:end);
cB = cov(:, 2:4:end);
cC = cov(:, 3:4:end);
cD = cov(:, 4:4:end);

cAGvsIPS = cov_real(1,:)-cov_real(2,:);
cAGvsSHAM = cov_real(1,:)-cov_real(3,:);
cIPSvsSHAM = cov_real(2,:)-cov_real(3,:);

pc = find(cC(:,2)>cAGvsSHAM(3))
length(pc)/1000


