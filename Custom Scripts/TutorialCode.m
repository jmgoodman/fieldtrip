function TutorialCode

% from: https://www.fieldtriptoolbox.org/tutorial/spikefield/

%% cleanup
clear
clc
close all

%% adjust path
mfp = mfilename('fullpath');

[path1up,~,~] = fileparts(mfp);
[path2up,~,~] = fileparts(path1up);

addpath(path2up)

ft_defaults % this adds all the requisite directories

addpath(fullfile(path2up,'Data')) % and here's the custom directory that you made to house data

%%
% addpath(fullfile(path2up,'fileio'))
% addpath(fullfile(path2up,'utilities'))
% addpath(fullfile(path2up,'Data'))
%
% vv = version;
%
% versionformat    = 'R(\d){4}(a|b)';
% versioncondensed = regexp(vv,versionformat,'match');
% versioncondensed = versioncondensed{1}(2:(end-1));
%
% dd = dir(fullfile(path2up,'compat'));
% versionnames
%
% addpath(fullfile(path2up,'compat','matlablt2017b'));
% addpath(fullfile(path2up,'compat','matlablt2016b'));

%% load spiking data
filename         = 'p029_sort_final_01.nex';
spike            = ft_read_spike(filename);

cfg              = [];
cfg.spikechannel = {'sig002a_wf','sig003a_wf'};
spike            = ft_spike_select(cfg,spike);

%% load LFP data
% get the cfg.trl
cfg          = [];
cfg.dataset  = filename;
cfg.trialfun = 'trialfun_stimon_samples';
cfg          = ft_definetrial(cfg);

% read in the data in trials
cfg.channel   = {'AD01', 'AD02', 'AD03', 'AD04'}; % these channels contain the LFP
cfg.padding   = 10; % length to which we pad for filtering
cfg.dftfreq   = [60-1*(1/10):(1/10):60+1*(1/10) ]; % filter out 60 hz line noise
cfg.dftfilter = 'yes';
data_lfp      = ft_preprocessing(cfg); % read in the LFP

%% relationship between LFP and spike timestamps
% data_lfp.sampleinfo has firsttimestamp and such for each trial

% this is just "pseudocode"
% sample = double(ts - FirstTimeStamp) / double(TimeStampPerSample) + 1;

%% organize spikes into the same trial structure as the LFP data

cfg           = [];
cfg.dataset   = filename;
cfg.trialfun  = 'trialfun_stimon_samples';
cfg           = ft_definetrial(cfg);
trl           = cfg.trl;

cfg           = [];
cfg.hdr       = data_lfp.hdr; % contains information for conversion of samples to timestamps
cfg.trlunit   = 'samples';
cfg.trl       = trl; % now in samples
spikeTrials   = ft_spike_maketrials(cfg,spike);

%% organize into trial structure using a second method
cfg          = [];
cfg.dataset  = filename;
cfg.trialfun = 'trialfun_stimon'; % this was defined in the spike tutorial
cfg          = ft_definetrial(cfg);
cfg.timestampspersecond = 40000;
spikeTrials2 = ft_spike_maketrials(cfg,spike);

%%
data_all = ft_appendspike([],data_lfp, spike);

%% excise the portions of the LFP near the spike times - otherwise, the LFP contains the spikes and doing spike-field coherence starts to get weird...
cfg              = [];
cfg.method       = 'nan'; % replace the removed segment with nans
cfg.timwin       = [-0.002 0.002]; % remove 4 ms around every spike
cfg.spikechannel = spike.label{1};
cfg.channel      = data_lfp.label(2);
data_nan         = ft_spiketriggeredinterpolation(cfg, data_all);

cfg.method       = 'linear'; % remove the replaced segment with interpolation
data_i           = ft_spiketriggeredinterpolation(cfg, data_all);

%% plot the result
figure,
plot(data_i.time{1},data_i.trial{1}(2,:),'g-'), hold on, plot(data_i.time{1}, data_i.trial{1}(5,:),'r')
hold on
plot(data_nan.time{1},data_nan.trial{1}(2,:),'go')
hold on
plot(data_all.time{1},data_all.trial{1}(2,:),'k-')
xlabel('time (s)')

%% compute the spike-triggered average LFP (and plot it!)
% note: this INCLUDES the portions of data near every spike (so... what we did 2 steps ago was all for naught? am confused, tutorial...)
cfg              = [];
cfg.timwin       = [-0.25 0.25]; % take 400 ms
cfg.spikechannel = spike.label{1}; % first unit
cfg.channel      = data_lfp.label(1:4); % first four chans
cfg.latency      = [0.3 10];
staPost          = ft_spiketriggeredaverage(cfg, data_all);

% plot the sta
figure
plot(staPost.time, staPost.avg(:,:)')
legend(data_lfp.label)
xlabel('time (s)')
xlim(cfg.timwin)

%% STA for a strictly pre-stimulus window
% hmmm... this part of the tutorial is producing the same output as the previous part... but it SHOULD be giving a different result.
% dammit! now I have to debug something...

cfg              = [];
cfg.timwin       = [-0.25 0.25]; % take 400 ms
cfg.spikechannel = spike.label{1}; % first unit
cfg.channel      = data_lfp.label(1:4); % first four chans
cfg.latency      = [-2.75 0];
staPre           = ft_spiketriggeredaverage(cfg, data_all);

figure
plot(staPre.time, staPre.avg(:,:)')
legend(data_lfp.label)
xlabel('time (s)')
xlim(cfg.timwin)

%% subfunctions
    function trl = trialfun_stimon_samples(cfg)
        hdr   = ft_read_header(cfg.dataset);
        event = ft_read_event(cfg.dataset);
        correctresponse  = 10041;
        begintrial       = 10044;
        endtrial         = 10045;
        stimon           = 10030;
        distractorChange = 12000;
        targetChange     = 12001;
        attCnds          = 20001:20004; % att in/out by target change first/second
        E          = struct2cell(event);
        samples    = cell2mat(E(1,:)); % now in vector form
        value      = cell2mat(E(2,:));
        begmark    = find(value==begintrial); % loop through the trial beginnings
        endmark    = find(value==endtrial); % loop through the trial beginnings
        trl        = []; % initialize the cfg.trl
        for k=1:length(begmark)
            vals = value(begmark(k):endmark(k));
            if any(ismember(vals,attCnds)) && ~isempty(find(vals==correctresponse))
                % create the trl matrix in sample units
                samp = samples(begmark(k):endmark(k)); % in timestamp units
                beginSamp      = samp(find(vals==stimon));
                sampDistractor = samp(find(vals==distractorChange));
                sampTarget     = samp(find(vals==targetChange));
                endSamp        = min([sampTarget(:);sampDistractor(:)]); % limit until first change
                offset         = -round(hdr.Fs*2.75);
                trl            = [trl; [beginSamp+offset endSamp offset]];
            end
        end
        
    end

    function trl = trialfun_stimon(cfg)
        
        hdr   = ft_read_header(cfg.dataset);
        event = ft_read_event(cfg.dataset);
        correctresponse  = 10041;
        begintrial       = 10044;
        endtrial         = 10045;
        stimon           = 10030;
        distractorChange = 12000;
        targetChange     = 12001;
        attCnds          = 20001:20004;       % att in/out by target change first/second
        E          = struct2cell(event);
        samples    = cell2mat(E(1,:));        % now in vector form
        value      = cell2mat(E(2,:));
        timestamps = cell2mat(E(3,:));        % now in vector form
        begmark    = find(value==begintrial); % loop through the trial beginnings
        endmark    = find(value==endtrial);   % loop through the trial beginnings
        trl = [];
        for k=1:length(begmark)
            vals = value(begmark(k):endmark(k));
            if any(ismember(vals,attCnds)) && ~isempty(find(vals==correctresponse))
                ts = timestamps(begmark(k):endmark(k)); % in timestamp units
                beginTs      = ts(find(vals==stimon));
                tsDistractor = ts(find(vals==distractorChange));
                tsTarget     = ts(find(vals==targetChange));
                endTs        = min([tsTarget(:);tsDistractor(:)]);    % limit until first change
                offset       = - hdr.Fs*hdr.TimeStampPerSample*2.75;  % 40000 timestamps per second x 2.75 sec
                trl          = [trl; [beginTs+offset endTs offset]];
            end
        end
    end

%%
end