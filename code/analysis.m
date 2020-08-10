%% Extract and process raw data
clear all; close all;
p = loadProjectDataset('GCaMP_LearningGrating2AFC');
numSess = height(p.sessionInfo);

%Areas to extract activity
regions = {'Left VTA','Right VTA','Left DMS','Right DMS','Left NAc','Right NAc'};

%Get data from each mouse
t_sample_epoch = linspace(-0.5,0.8,100);
t_sample_wholetrial = linspace(-0.2,4,400);
warp_sizes = [50,100,20,100]; %number of elements for each epoch: pre-stim, stim-choice, choice-outcome, post-outcome
D = table;
for sess = 1:numSess
    fprintf('%d/%d %s\n',sess,numSess,p.sessionInfo.expRef{sess});
    
    %Load behav data
    b = getBehavData(p.sessionInfo.expRef{sess},p.sessionInfo.dataProfile{sess});
    
    %Get photometry data
    [F,t]=photometryAlign(p.sessionInfo.expRef{sess}, b.rewardTime(~isnan(b.rewardTime)));
%     
%     %Load photometry data
%     e=strsplit(p.sessionInfo.expRef{sess},'_');
%     thisDate = e{1};
%     thisSess = e{2};
%     thisSubj = e{3};
%     filepath = sprintf("\\\\QNAP-AL001.dpag.ox.ac.uk\\Data\\%s\\%s\\photoM\\%s_%s_%s_F.mat",thisSubj,thisDate,thisSubj,thisDate,thisSess);
%     load(filepath);
%     Water = photoMdata.TTL_1; % water TTL signal
%     TimeStamps=photoMdata.Time_s_;
%     
%     % Align dF/F times to block times
%     [~,reward_echo_photoM] = risetime(Water,TimeStamps);
%     switch(p.sessionInfo.dataProfile{sess})
%         case 'Grating2AFC_choiceWorld'
%             reward_echo_behav = b.rewardTime(~isnan(b.rewardTime));
%             stimTime = b.stimulusOnsetTime;
%         case 'Grating2AFC_choiceWorld_noTimeline'
%             reward_echo_behav = b.rewardTimeEstimate(~isnan(b.rewardTimeEstimate));
%             stimTime = b.stimulusOnsetTimeEstimate;
%             outcomeTime = nanmean([b.rewardTimeEstimate, b.punishSoundOnsetTimeEstimate],2);
%     end
%     
%     [~,photoM2block] = makeCorrection(reward_echo_behav, reward_echo_photoM, false); % find the correction factor
%     deltaFoverF= photoMdata.AnalogIn_2_dF_F0;
%     traceTime = applyCorrection(TimeStamps, photoM2block);  % photoM data is on the Timeline now
%     
%     b.stimTime = stimTime;
    b.outcomeTime = nanmean([b.rewardTime, b.punishSoundOnsetTime],2);
    
    b.dff_stim = nan(length(b.choice),length(t_sample_epoch),length(regions));
    b.dff_choice = nan(length(b.choice),length(t_sample_epoch),length(regions));
    b.dff_outcome = nan(length(b.choice),length(t_sample_epoch),length(regions));
    b.dff_timewarped = nan(length(b.choice),sum(warp_sizes),length(regions));
    b.dff_wholetrial = nan(length(b.choice),length(t_sample_wholetrial),length(regions));
    
    %for each possible area
    for a = 1:length(regions)
        
        %check whether this area is measured for any of the channels
        if strcmp(p.sessionInfo.photometryChannel2{sess}, regions{a})
            dff = zscore(F(:,1));
        elseif strcmp(p.sessionInfo.photometryChannel4{sess}, regions{a})
        	dff = zscore(F(:,2));
        else
            dff = nan(size(t));
        end
        
        b.dff_stim(:,:,a) = interp1(t,dff,b.stimulusOnsetTime + t_sample_epoch);
        b.dff_choice(:,:,a) = interp1(t,dff,b.choiceStartTime + t_sample_epoch);
        b.dff_outcome(:,:,a) = interp1(t,dff,b.outcomeTime + t_sample_epoch);
        b.dff_wholetrial(:,:,a) = interp1(t,dff, b.stimulusOnsetTime + t_sample_wholetrial );

        %get time-warped activity for each trial
        warp_samples = nan(length(b.choice), sum(warp_sizes));
        for tr = 1:length(b.choice)
            epoch1 = linspace(b.stimulusOnsetTime(tr)-0.5, b.stimulusOnsetTime(tr), warp_sizes(1));
            epoch2 = linspace(b.stimulusOnsetTime(tr), b.choiceStartTime(tr), warp_sizes(2));
            epoch3 = linspace(b.choiceStartTime(tr), b.outcomeTime(tr), warp_sizes(3));
            epoch4 = linspace(b.outcomeTime(tr), b.outcomeTime(tr)+1, warp_sizes(4));
            warp_samples(tr,:) = [epoch1, epoch2, epoch3, epoch4];
        end
        
        b.dff_timewarped(:,:,a) = interp1(t,dff,warp_samples);
    end
    
    %baseline subtract
    pre_stim_baseline = mean(b.dff_stim(:,t_sample_epoch<0,:),2);
    b.dff_stim = b.dff_stim - pre_stim_baseline;
    b.dff_outcome = b.dff_outcome - pre_stim_baseline;
    b.dff_outcome = b.dff_outcome - pre_stim_baseline;
    b.dff_wholetrial = b.dff_wholetrial - pre_stim_baseline;
    b.dff_timewarped = b.dff_timewarped - mean(b.dff_timewarped(:,1:warp_sizes(1),:),2);
%  
    %calculate RT
    b.RT = b.choiceStartTime - b.stimulusOnsetTime;
    b.expRef = repmat(p.sessionInfo.expRef(sess), length(b.choice), 1);
    
    %add wheel position data
    block = dat.loadBlock(p.sessionInfo.expRef{sess});
    b.wheel_stim =  block.inputSensorGain*interp1(block.inputSensorPositionTimes,block.inputSensorPositions,b.stimulusOnsetTime + t_sample_epoch);
    b.wheel_choice =  block.inputSensorGain*interp1(block.inputSensorPositionTimes,block.inputSensorPositions,b.choiceStartTime + t_sample_epoch);
    pre_stim_baseline = mean(b.wheel_stim(:,t_sample_epoch<0),2);
    b.wheel_stim = b.wheel_stim - pre_stim_baseline; %pre-stim baseline
    b.wheel_choice = b.wheel_choice - pre_stim_baseline; %pre-stim baseline

    %remove trials where the DA signal wasn't measured
    badIdx = isnan(nanmean(b.dff_stim(:,1,:),3)) | isnan(nanmean(b.dff_outcome(:,1,:),3));
    b(badIdx,:)=[];
    b.trialNumber = (1:length(b.choice))';
    
    %add to big table
    D = vertcat(D, b(:,{'expRef','trialNumber','contrastLeft','contrastRight',...
                        'choice','stimulusOnsetTime','choiceStartTime','outcomeTime',...
                        'RT','feedback','rewardVolume',...
                        'dff_stim','dff_choice','dff_outcome','dff_timewarped','dff_wholetrial',...
                        'wheel_stim','wheel_choice'}) );
end

save('../data/compiled.mat','p','D','t_sample_epoch','t_sample_wholetrial','regions','warp_sizes','-v7.3');
%% **: Load data and compute binned DA values
clear all; close all;
load('../data/compiled.mat');
numSess = height(p.sessionInfo);

stimWindow = [0.4 0.6];
% stimWindow = [0.3 0.4];
outcomeWindow = [0.35 0.4];

%Remove trials with very large RTs, so warping isn't so extreme.
D = D(D.RT<5,:);

%Extract activity over windows
D.dff_stim_binned = mean(D.dff_stim(:, stimWindow(1) <= t_sample_epoch & t_sample_epoch <= stimWindow(2),:),2);
D.dff_outcome_binned = mean(D.dff_outcome(:, outcomeWindow(1) <= t_sample_epoch & t_sample_epoch <= outcomeWindow(2),:),2);
D.dff_stim_binned(D.RT < stimWindow(2),:,:)=NaN; %set binned value to NaN if it contained the choice time

%Create contrast variable (useful later)
D.cDiff = D.contrastRight - D.contrastLeft;
num_stim_cond = length(unique(D.cDiff));

%add mouse information
D = innerjoin(D,p.sessionInfo);

%re-sort by mouse, session and trial number
D = sortrows(D,{'mouseName','sessionNum','trialNumber'});

%Compute zscored RT and psychometric weights
D.RTz = nan(size(D.RT));
p.sessionInfo.psych_B = nan(height(p.sessionInfo),1);
p.sessionInfo.psych_SL = nan(height(p.sessionInfo),1);
p.sessionInfo.psych_SR = nan(height(p.sessionInfo),1);

%Per session quantities
for sess = 1:numSess
    %zscore RT
    idx = strcmp(D.expRef,p.sessionInfo.expRef{sess});
    D.RTz(idx) = (D.RT(idx) - median(D.RT(idx)))/(1.486 * mad(D.RT(idx),1)); %modified z score
    
    %Compute psychometric weights for each session
    b = glmfit([D.contrastLeft(idx) D.contrastRight(idx)],D.choice(idx)=='Right choice', 'binomial','Constant','on');
    p.sessionInfo.psych_B(sess) = b(1);
    p.sessionInfo.psych_SL(sess) = b(2);
    p.sessionInfo.psych_SR(sess) = b(3);
end

%get the VTA activity (makes things easier later)
D.VTA_binned = nanmean( D.dff_stim_binned(:,:,contains(regions,'VTA')), 3);
%% Plot psychometric curves per session
g = gramm('x',D.cDiff,'y',double(D.choice=='Right choice'),'lightness',D.mouseName);
g.facet_grid(D.mouseName,D.sessionNum);
g.stat_summary('geom',{'point','line','errorbar'},'type',@grammCallbackBinomialConfidenceInterval);
g.set_layout_options('redraw',true,'redraw_gap',0.001);
g.axe_property('ylim',[0 1]);
g.geom_vline('xintercept',0,'style','k:');
g.geom_hline('yintercept',0.5,'style','k:');
g.set_names('x','Contrast','y','pR','column','','row','');
g.set_title('Psychometric curves');
% g.stat_glm('distribution','binomial','fullrange',true,'geom',{'line','black'});

figure('units','normalized','outerposition',[0 0 1 1],'color','w')
g.draw();

g.export('file_name','PsychometricCurves','export_path','../figures/behavioural/','file_type','pdf');
%% Plot wheel traces

%Calculate wheel velocity
D.wheelV_stim = diff(D.wheel_stim,[],2);
D.wheelV_choice = diff(D.wheel_choice,[],2);

%flip sign for right choices so velocity is always positive towards the
%chosen side
D.wheelV_stim(D.choice=='Right choice',:) = -D.wheelV_stim(D.choice=='Right choice',:);
D.wheelV_choice(D.choice=='Right choice',:) = -D.wheelV_choice(D.choice=='Right choice',:);

E = D(D.feedback=='Rewarded',:);

[G,G_labels]=groupsummary(E.wheelV_stim,{E.mouseName E.sessionNum E.choice},'mean');
g(1,1) = gramm('x',t_sample_epoch(1:end-1), 'y', G, 'color', G_labels{3});
g(1,1).facet_wrap(G_labels{1},'scale','independent');
g(1,1).stat_summary('type','ci','setylim',true);
g(1,1).geom_vline('xintercept',0,'style','k:');
g(1,1).set_title('Stimulus aligned');

[G,G_labels]=groupsummary(E.wheelV_choice,{E.mouseName E.sessionNum E.choice},'mean');
g(1,2) = gramm('x',t_sample_epoch(1:end-1), 'y', G, 'color', G_labels{3});
g(1,2).facet_wrap(G_labels{1},'scale','independent');
g(1,2).stat_summary('type','ci','setylim',true);
g(1,2).geom_vline('xintercept',0,'style','k:');
g(1,2).set_title('Choice aligned');

g.set_names('x','Time','y','Wheel Speed (??/sec)','column','','row','','color','');
g.set_layout_options('redraw',true,'redraw_gap',0.01);

figure('units','normalized','outerposition',[0 0 1 1],'color','w')
g.draw();

g.export('file_name','WheelVelocity','export_path','../figures/behavioural/','file_type','pdf');
%% Plot psych curves conditioned on history

%Conditioned on previous choice
clear g;
g = gramm('x',D.cDiff,'y',double(D.choice=='Right choice'),'color',circshift(D.choice,1));
g.facet_grid(D.mouseName,D.sessionNum);
g.stat_summary('geom',{'point','line','errorbar'},'type',@grammCallbackBinomialConfidenceInterval);
g.axe_property('ylim',[0 1]);
g.geom_vline('xintercept',0,'style','k:');
g.geom_hline('yintercept',0.5,'style','k:');
g.set_names('x','Contrast','y','p(Right)','column','','row','','color','');
g.set_title('Psychometric curve conditioned on previous choice');
% g.stat_glm('distribution','binomial','fullrange',true);
g.set_layout_options('redraw',true,'redraw_gap',0.005,'legend_width',0.09);
figure('units','normalized','outerposition',[0 0 1 1],'color','w')
g.draw();
g.export('file_name','PsychHist_PrevChoice','export_path','../figures/behavioural/','file_type','pdf');

%Conditioned on previous feedback for Right choice
clear g;
g = gramm('x',D.cDiff,'y',double(D.choice=='Right choice'),'color',circshift(D.feedback,1),'subset',circshift(D.choice,1)=='Right choice');
g.facet_grid(D.mouseName,D.sessionNum);
g.stat_summary('geom',{'point','line','errorbar'},'type',@grammCallbackBinomialConfidenceInterval);
g.axe_property('ylim',[0 1]);
g.geom_vline('xintercept',0,'style','k:');
g.geom_hline('yintercept',0.5,'style','k:');
g.set_names('x','Contrast','y','p(Right)','column','','row','','color','');
g.set_title('Psychometric curve conditioned on previous outcome for Right choice');
g.set_layout_options('redraw',true,'redraw_gap',0.005,'legend_width',0.09);
figure('units','normalized','outerposition',[0 0 1 1],'color','w')
g.draw();
g.export('file_name','PsychHist_PrevOutcomeRightChoice','export_path','../figures/behavioural/','file_type','pdf');

%% Plot RT medians across contrasts per session
g = gramm('x',D.cDiff,'y',D.RT,'color',D.choice,'subset',D.feedback=='Rewarded');
g.facet_grid(D.mouseName,D.sessionNum,'scale','free_y');
% g.stat_summary('geom',{'point','line','errorbar'},'type','quartile','setylim',true);
% custom_statfun = @(y) ([median(y), median(y)-mad(y,1), median(y)+mad(y,1)]);
g.stat_summary('geom',{'point','line','errorbar'},'setylim',true,'type',@(y) [median(y); median(y)-mad(y,1); median(y)+mad(y,1)]);
g.set_layout_options('redraw',true,'redraw_gap',0.005);
g.geom_vline('xintercept',0,'style','k:');
g.set_names('x','Contrast','y','RT','column','','row','');
% g.axe_property('ylim',[0 10]);
g.set_title('Reaction time (median +- mad) across all correct trials');

figure('units','normalized','outerposition',[0 0 1 1],'color','w')
g.draw();

g.export('file_name','ReactionTimeByContrast','export_path','../figures/behavioural','file_type','pdf');
%% Plot timewarped DA traces

%For each region separately
for a = 1:length(regions)
    
    %get list of sessions which have recording in this region
    idx = strcmp(p.sessionInfo.photometryChannel2, regions{a}) | strcmp(p.sessionInfo.photometryChannel4, regions{a});
    eRefs = p.sessionInfo.expRef(idx);
    
    E = D(ismember(D.expRef,eRefs),:);
    
    %time warped trace
    clear g;
    g(1,1) = gramm('x',1:sum(warp_sizes),'y',E.dff_timewarped(:,:,a),'color',E.cDiff,'subset',E.feedback=='Rewarded');
    g(1,1).facet_grid(E.sessionNum,E.mouseName,'scale','independent');
    g(1,1).stat_summary('setylim','true');
    g(1,1).set_color_options('map',0.9*RedWhiteBlue(floor(length(unique(E.cDiff))/2)));
    g(1,1).geom_vline('xintercept',cumsum(warp_sizes(1:3)),'style','k:');
    g(1,1).set_title('Rewarded trials separated by contrast');
    g(1,1).set_names('x','','y','zDFF','color','contrast','column','','row','');

    g.set_layout_options('redraw',true,'redraw_gap',0);
    g.axe_property('xtick','','xcolor','none');
    g.set_title(regions{a});
    
    figure('units','normalized','outerposition',[0.0042 0.0074 0.5380 0.9139],'color','w');
    g.draw();

    %match Y axes in columns and set limit to range of values
    numSess = groupsummary(E.sessionNum,{E.mouseName},'max');
    for c = 1:size(g.facet_axes_handles,2)
        ylim = [median(cellfun(@min,get(g.facet_axes_handles(1:numSess(c),c),'ylim'))),...
            median(cellfun(@max,get(g.facet_axes_handles(1:numSess(c),c),'ylim')))];
        
        set(g.facet_axes_handles(:,c),'ylim',ylim);
    end
    
%     g.export('file_name',sprintf('%s',strrep(regions{a},' ','')),'export_path','../figures/fluorescence_timewarped/','file_type','pdf');
    
end
%% Plot DA traces and tuning curves

%For each region separately
for a = 1:length(regions)
    
    %get list of sessions which have recording in this region
    idx = strcmp(p.sessionInfo.photometryChannel2, regions{a}) | strcmp(p.sessionInfo.photometryChannel4, regions{a});
    eRefs = p.sessionInfo.expRef(idx);
    
    E = D(ismember(D.expRef,eRefs),:);
  
    %traces
    clear g;
    g(1,1) = gramm('x',t_sample_epoch, 'y', E.dff_stim(:,:,a), 'color', E.cDiff, 'subset', E.feedback=='Rewarded');
    g(1,1).facet_grid(E.mouseName,E.sessionNum,'scale','free_y');
    g(1,1).set_names('x','','y','ZScore DA','color','C','column','','row','');
    g(1,1).stat_summary('setylim','true');
    g(1,1).set_title('DA after stimulus (rewarded trials)');
    g(1,1).set_color_options('map',0.9*RedWhiteBlue(floor(length(unique(E.cDiff))/2)));
    
    g(1,2) = gramm('x',t_sample_epoch, 'y', E.dff_outcome(:,:,a), 'color', E.cDiff, 'subset', E.feedback=='Rewarded');
    g(1,2).facet_grid(E.mouseName,E.sessionNum,'scale','free_y');
    g(1,2).set_names('x','','y','ZScore DA','color','C','column','','row','');
    g(1,2).stat_summary('setylim','true');
    g(1,2).set_title('DA after outcome (rewarded trials)');
    g(1,2).set_color_options('map',0.9*RedWhiteBlue(floor(length(unique(E.cDiff))/2)));
    
    g.geom_vline('xintercept',0,'style','k:');
    g.geom_hline('yintercept',0,'style','k:');
    g.set_layout_options('redraw',true,'redraw_gap',0.01);
    g.set_title(regions{a});
    
    figure('units','normalized','outerposition',[0 0 1 1],'color','w')
    g.draw();
    
    set(g(1,1).facet_axes_handles(:,2:end),'ytick','','ycolor','none');
    set(g(1,1).facet_axes_handles(1:end-1,:),'xtick','','xcolor','none');
    set(g(1,2).facet_axes_handles(:,2:end),'ytick','','ycolor','none');
    set(g(1,2).facet_axes_handles(1:end-1,:),'xtick','','xcolor','none');
%     g.export('file_name',sprintf('%s',strrep(regions{a},' ','')),'export_path','../figures/fluorescence/','file_type','pdf');
    
    %tuning curves
    clear g;
    g(1,1) = gramm('x',E.cDiff,'y',E.dff_stim_binned(:,:,a),'subset',E.feedback=='Rewarded');
    g(1,1).facet_grid(E.mouseName,E.sessionNum,'scale','free_y');
    g(1,1).geom_vline('xintercept',0,'style','k:');
    g(1,1).stat_summary('geom',{'point','line','errorbar'},'setylim','true');
    g(1,1).set_names('x','Contrast','y',sprintf('%0.2f-%0.2f @ stim',stimWindow(1),stimWindow(2)),'column','','row','','color','');
    g(1,1).set_title(sprintf('Tuning curves of binned [%0.2f-%0.2f] post-stim dF/F for rewarded trials',stimWindow(1),stimWindow(2)));
    
    g(1,2) = gramm('x',E.cDiff,'y',E.dff_outcome_binned(:,:,a),'color',E.feedback);
    g(1,2).facet_grid(E.mouseName,E.sessionNum,'scale','free_y');
    g(1,2).geom_vline('xintercept',0,'style','k:');
    g(1,2).stat_summary('geom',{'point','line','errorbar'},'setylim','true');
    g(1,2).set_names('x','Contrast','y',sprintf('%0.2f-%0.2f @ outcome',outcomeWindow(1),outcomeWindow(2)),'column','','row','','color','');
    g(1,2).set_title(sprintf('Tuning curves of binned [%0.2f-%0.2f] post-outcome dF/F',outcomeWindow(1),outcomeWindow(2)));
    g.set_title(regions{a});

    g.set_layout_options('redraw',true,'redraw_gap',0.01);
    
    figure('units','normalized','outerposition',[0 0 1 1],'color','w');
    g.draw();
    set(g(1,1).facet_axes_handles(:,2:end),'ytick','','ycolor','none');
    set(g(1,1).facet_axes_handles(1:end-1,:),'xtick','','xcolor','none');
    set(g(1,2).facet_axes_handles(:,2:end),'ytick','','ycolor','none');
    set(g(1,2).facet_axes_handles(1:end-1,:),'xtick','','xcolor','none');
%     g.export('file_name',sprintf('%s',strrep(regions{a},' ','')),'export_path','../figures/tuning_curves','file_type','pdf');

end
%% Plot DA, Psych and zRT plots for each mouse

[mice,~,miceID] = unique(D.mouseName);
for m = 1:max(miceID)
    
    E =  D(strcmp(D.mouseName,mice{m}),:);
    
    clear g;
    %Psych curve
    g(1,1) = gramm('x',E.cDiff,'y',double(E.choice=='Right choice'));
    g(1,1).facet_grid([],E.sessionNum,'column_labels',true);
    g(1,1).stat_summary('geom',{'point','line','errorbar'},'type',@grammCallbackBinomialConfidenceInterval);
    g(1,1).axe_property('ylim',[0 1]);
    g(1,1).geom_vline('xintercept',0,'style','k:');
    g(1,1).geom_hline('yintercept',0.5,'style','k:');
    g(1,1).set_names('x','Contrast','y','p(Right) + 95CI','Column','');
    g(1,1).set_layout_options('legend',true);
    
    %RT distribution
    g(2,1) = gramm('x',E.cDiff,'y',E.RT,'color',E.choice,'subset',E.feedback=='Rewarded');
    g(2,1).facet_grid([],E.sessionNum,'column_labels',false);
    g(2,1).stat_summary('geom',{'point','line','errorbar'},'type',@(y) [median(y); median(y)-mad(y,1); median(y)+mad(y,1)],'setylim',true);
    g(2,1).set_names('x','Contrast','y','RT median+MAD','color','');
    
    
    %which regions this mouse has
    thisMouseRegions = unique([E.photometryChannel2;E.photometryChannel4]);
    thisMouseRegions(cellfun(@isempty, thisMouseRegions)) = [];
    
    for a = 1:length(thisMouseRegions)
        x = 2+a;
        g(x,1) = gramm('x',t_sample_epoch, 'y', E.dff_stim(:,:,strcmp(regions,thisMouseRegions{a})), 'color', E.cDiff, 'subset', E.feedback=='Rewarded');
        g(x,1).facet_grid([],E.sessionNum,'column_labels',false);
        g(x,1).set_names('x','Time','y',sprintf('%s @ stim',thisMouseRegions{a}),'color','C');
        g(x,1).stat_summary('setylim',true);
        g(x,1).set_color_options('map',0.9*RedWhiteBlue(floor(length(unique(E.cDiff))/2)));
        g(x,1).geom_vline('xintercept',0,'style','k:');
        g(x,1).geom_hline('yintercept',0,'style','k:');
    end
    
    for a = 1:length(thisMouseRegions)
        x = 2+length(thisMouseRegions)+a;
        g(x,1) = gramm('x',t_sample_epoch, 'y', E.dff_outcome(:,:,strcmp(regions,thisMouseRegions{a})), 'color', E.cDiff, 'subset', E.feedback=='Rewarded');
        g(x,1).facet_grid([],E.sessionNum,'column_labels',false);
        g(x,1).set_names('x','Time','y',sprintf('%s @ reward',thisMouseRegions{a}),'color','C');
        g(x,1).stat_summary('setylim',true);
        g(x,1).set_color_options('map',0.9*RedWhiteBlue(floor(length(unique(E.cDiff))/2)));
        g(x,1).geom_hline('yintercept',0,'style','k:');
        g(x,1).geom_vline('xintercept',0,'style','k:');
    end
    
    g.set_title( mice{m} );
    g.set_layout_options('redraw',true,'redraw_gap',0.01,'legend_width',0.1);

    figure('units','normalized','outerposition',[0 0 1 1],'color','w')
    g.draw();

    g.export('file_name',mice{m},'export_path','../figures/','file_type','png');

end

%% Plot DA decoding

for a = 1:length(regions)
    
    %get list of sessions which have recording in this region
    idx = strcmp(p.sessionInfo.photometryChannel2, regions{a}) | strcmp(p.sessionInfo.photometryChannel4, regions{a});
    
    g(1,1) = gramm('x',t_sample_epoch,'y',p.sessionInfo.choiceDecoder(idx,:,a),'lightness',p.sessionInfo.sessionNum(idx));
%     g(1,1).facet_grid(p.sessionInfo.mouseName(idx),p.sessionInfo.sessionNum(idx));
    g(1,1).facet_grid(p.sessionInfo.mouseName(idx),[]);
    g(1,1).geom_line();
    g(1,1).stat_summary();
    g(1,1).set_names('x','Move onset','y','Choice dec','color','','row','','column','');
    g(1,1).set_title('Choice decoder (controlling for stimulus');
     
    g(1,2) = gramm('x',t_sample_epoch,'y',p.sessionInfo.rewardDecoder(idx,:,a),'lightness',p.sessionInfo.sessionNum(idx));
%     g(1,2).facet_grid(p.sessionInfo.mouseName(idx),p.sessionInfo.sessionNum(idx));
    g(1,2).facet_grid(p.sessionInfo.mouseName(idx),[]);
    g(1,2).geom_line();
     g(1,2).stat_summary();
    g(1,2).set_names('x','Outcome onset','y','Reward dec','color','','row','','column','');
    g(1,2).set_title('Reward decoder (controlling for stimulus');
     
    g.set_layout_options('redraw',true,'redraw_gap',0.01);
    g.set_title(regions{a});
        
    g.geom_vline('xintercept',0,'style','k:');
    g.geom_hline('yintercept',0.5,'style','k:');
    
    
    figure('units','normalized','outerposition',[0 0 1 1],'color','w');
    g.draw();

    g.export('file_name',sprintf('%s_Decoding',strrep(regions{a},' ','')),'export_path','../figures/','file_type','pdf');

end
%% Plot Binned DA over trials for each mouse, concetenated sessions
numTrialsSmoothWindow = 200;

[mice,~,miceID] = unique(D.mouseName);
D.trialNumberConcatenated = nan(size(D.trialNumber));
for m = 1:length(mice)
    idx = strcmp(D.mouseName,mice{m});
    D.trialNumberConcatenated(idx) = (1:sum(idx))';
end

clear g;
g(1,1) = gramm('x',D.trialNumberConcatenated,'y',D.DAstim_binned,'color',D.cDiff,'subset',D.feedback=='Rewarded');
g(1,1).facet_grid(D.mouseName,[],'scale','independent');
g(1,1).stat_smooth('method','moving','lambda',numTrialsSmoothWindow);
g(1,1).set_names('x','Trial number','y',sprintf('%0.2f-%0.2f @ stim',stimWindow(1),stimWindow(2)));
g(1,1).axe_property('ylim',[-0.5 2.5]);
g(1,1).set_title('DA after stimulus');


g(1,2) = gramm('x',D.trialNumberConcatenated,'y',D.DAoutcome_binned,'color',D.cDiff,'subset',D.feedback=='Rewarded');
g(1,2).facet_grid(D.mouseName,[],'scale','independent');
g(1,2).stat_smooth('method','moving','lambda',numTrialsSmoothWindow);
g(1,2).set_names('x','Trial number','y',sprintf('%0.2f-%0.2f @ outcome',outcomeWindow(1),outcomeWindow(2)));
g(1,2).axe_property('ylim',[-0.5 3.5]);
g(1,2).set_title('DA after Reward');

g.set_color_options('map',0.9*RedWhiteBlue(floor(num_stim_cond/2)));
g.set_layout_options('redraw',true,'redraw_gap',0.01);

figure('units','normalized','outerposition',[0 0 1 1],'color','w')
g.draw();

%add session boundaries
for m = 1:max(miceID)
    tn = D.trialNumber(miceID==m);
    
    arrayfun(@(x) xline(g(1).facet_axes_handles(m), x, 'k:'), find(diff(tn)<0));
    arrayfun(@(x) xline(g(2).facet_axes_handles(m), x, 'k:'), find(diff(tn)<0));
end

g.export('file_name',sprintf('DFF_correctTrials_%dTrialSmoothing',numTrialsSmoothWindow),'export_path','../figures/','file_type','pdf');
%% Plot trace of contralateral stimulus decoding

g = gramm('x',t_sample_epoch, 'y', p.auROC);
g.facet_grid(p.mouseName,p.sessionNum,'scale','free_y');
g.geom_line();
g.geom_vline('xintercept',0,'style','k:');
g.geom_hline('yintercept',0.5,'style','k:');
g.set_names('x','Stim onset','y','Stim deooding','column','','row','');
g.set_layout_options('redraw',true,'redraw_gap',0.01);
g.set_title('Contralateral stimulus decoding');

figure('units','normalized','outerposition',[0 0 1 1],'color','w')
g.draw();
g.export('file_name','DecodingContraStimulus_perSession','export_path','../figures/','file_type','pdf');
%% HighC: Block trials and correlate performance and binned DA

numTrialsWindow = 30;
C_condition = 1;

%Get only highest contrast trials
E = D(abs(D.cDiff) == C_condition,:);

%re-label trials in sequence and then create blocks
[mice,~,miceID] = unique(E.mouseName);
E.trialNumberConcatenated = nan(size(E.trialNumber));
E.block = nan(size(E.trialNumber));
E.blockGroup = nan(size(E.trialNumber));
for m = 1:length(mice)
    idx = strcmp(E.mouseName,mice{m});
    E.trialNumberConcatenated(idx) = (1:sum(idx))';
    
    numBlocks = floor(sum(idx)/numTrialsWindow);
    E.block(idx) = discretize(E.trialNumberConcatenated(idx),numBlocks);
    E.blockGroup(idx) = discretize(E.block(idx),3);
end
E.blockGroup = categorical(E.blockGroup,1:3,{'1/3','2/3','3/3'});

%Performance
[G_PERF,G1]=groupsummary(E.feedback=='Rewarded',{E.mouseName,E.blockGroup,E.block}, 'mean');

%Dopamine at the time of stimulus
[G_DASTIM]=groupsummary(E.DAstim_binned,{E.mouseName,E.blockGroup,E.blockGroup,E.block}, 'mean');

%Dopamine at the time of outcome, for rewarded trials only
rewardedIdx = E.feedback=='Rewarded';
[G_DAOUT,G2]=groupsummary(E.DAoutcome_binned(rewardedIdx),{E.mouseName(rewardedIdx),E.blockGroup(rewardedIdx),E.blockGroup(rewardedIdx),E.block(rewardedIdx)}, 'mean');

clear g;
g(1,1) = gramm('x',G_PERF,'y',G_DASTIM,'color',G1{2});
g(1,1).set_names('x','Performance','y','DA @ Stim','color','Training stage');
g(1,1).facet_grid([],G1{1});
g(1,1).stat_glm('disp_fit',true);
g(1,1).geom_point();
g(1,1).set_title(sprintf('Dopamine after stimulus (C=%0.2f) [%0.2f-%0.2f], correlated with performance in blocks of %d trials',C_condition,stimWindow(1),stimWindow(2),numTrialsWindow));

g(2,1) = gramm('x',G_PERF,'y',G_DAOUT,'color',G2{2});
g(2,1).set_names('x','Performance','y','DA @ Reward','color','Training stage');
g(2,1).facet_grid([],G2{1});
g(2,1).stat_glm('disp_fit',true);
g(2,1).geom_point();
g(2,1).set_title(sprintf('Dopamine after reward (C=%0.2f) [%0.2f-%0.2f], correlated with performance in blocks of %d trials',C_condition,outcomeWindow(1),outcomeWindow(2),numTrialsWindow));
g.set_layout_options('redraw',true,'redraw_gap',0.01);

figure('units','normalized','outerposition',[0 0 1 1],'color','w')
g.draw();
g.export('file_name','CorrelatePerformanceDopamine','export_path','../figures/','file_type','pdf');
%% Blocked stimulus decoding

numTrialsWindow = 150;

decL = cell(5,1); 
decR = cell(5,1); 
mouse = cell(5,1);
block = cell(5,1);
perfL = cell(5,1);
perfR = cell(5,1);
rtL = cell(5,1);
rtR = cell(5,1);
for m = 1:5 %for each mouse
    Dm = D(D.mouseID==m,:);
    
    %go through each session and remove trials which aren't at the
    %highest contrast
    [maxC,sessID] = groupsummary(abs(Dm.cDiff),Dm.sessionID,'max');
    for sess = 1:max(sessID)
        Dm(Dm.sessionID==sess & abs(Dm.cDiff)~=maxC(sess),:)=[];
    end
    
    %Split into blocks
    numTrials = height(Dm);
    numBlocks = floor(numTrials/numTrialsWindow);
    blockID = discretize( 1:numTrials, numBlocks )';
    
    for b = 1:numBlocks
        idx = blockID==b;
        
        activity = Dm.DAstim_binned(idx);
        choice = Dm.choice(idx);
        
        %Decode whether stimulus was present on the Left, NOT
        %controlling for choice (since we're at high contrast
        %conditions so there aren't enough error trials
        [decL{m}(b,1),p] = binaryDecoder(activity, Dm.contrastLeft(idx)>0, {ones(size(choice))},'numShuffles',2000);
        if ~(p < 0.025 || p > 0.975)
            decL{m}(b,1) = 0.5;
        end
        
        [decR{m}(b,1),p] = binaryDecoder(activity, Dm.contrastRight(idx)>0, {ones(size(choice))},'numShuffles',2000);
        if ~(p < 0.025 || p > 0.975)
            decR{m}(b,1) = 0.5;
        end
        
        %performance L and R
        fb =  Dm.feedback(idx);
        perfL{m}(b,1) = nanmean(fb(choice=='Left choice')=='Rewarded');
        perfR{m}(b,1) = nanmean(fb(choice=='Right choice')=='Rewarded');
        
        
        %RT L and R
        RT = Dm.RT(idx);
        rtL{m}(b,1) = nanmedian(RT(choice=='Left choice'));
        rtR{m}(b,1) = nanmedian(RT(choice=='Right choice'));
        
    end
    
    mouse{m} = repmat(allExpRefs(m,1),numBlocks,1);
    block{m} = (1:numBlocks)';
end

%concatenate
mouseCat = cat(1,mouse{:});
blockCat = cat(1,block{:});
decLCat = cat(1,decL{:});
decRCat = cat(1,decR{:});
perfLCat = cat(1,perfL{:});
perfRCat = cat(1,perfR{:});
rtLCat = cat(1,rtL{:});
rtRCat = cat(1,rtR{:});

clear g;
side = repmat({'Left stim','Right stim'},length(decLCat),1);
g(1,1) = gramm('x',repmat(blockCat,2,1),'y',[decLCat;decRCat],'color',side(:));
g(1,1).facet_grid([],repmat(mouseCat,2,1));
g(1,1).geom_hline('yintercept',0.5,'style','k--');
g(1,1).geom_line();
g(1,1).axe_property('ylim',[0 1]);
g(1,1).set_names('x','Block','y','Stimulus Left decoder');

ch = repmat({'Left choice','Right choice'},length(perfLCat),1);
g(2,1) = gramm('x',repmat(blockCat,2,1),'y',[perfLCat;perfRCat],'color',ch(:));
g(2,1).facet_grid([],repmat(mouseCat,2,1));
g(2,1).geom_line();
g(2,1).axe_property('ylim',[0 1]);
g(2,1).set_names('x','Block','y','Performance');

g(3,1) = gramm('x',repmat(blockCat,2,1),'y',[rtLCat;rtRCat],'color',ch(:));
g(3,1).facet_grid([],repmat(mouseCat,2,1));
g(3,1).geom_line();
g(3,1).set_names('x','Block','y','Median RT');


figure; g.draw();
%% Correlate DA and RT asymmetries 

%For each session, compute the slope of VTA vs C, and RT vs C and compare
%those across sessions
p.sessionInfo.RTvC = NaN(height(p.sessionInfo),2);
p.sessionInfo.VTAvC = NaN(height(p.sessionInfo),2);
for sess = 1:height(p.sessionInfo)
    E = D( strcmp(D.expRef, p.sessionInfo.expRef{sess}), :);

    idxL = E.cDiff<=0;
    idxR = E.cDiff>=0;
    
    if length(unique(E.cDiff(idxR)))>1 & ~all(isnan(E.VTA_binned))
        %Calculate slope of VTA DA vs Contrast
        bL = glmfit(abs(E.cDiff(idxL)),E.VTA_binned(idxL),'normal','constant','on');
        bR = glmfit(abs(E.cDiff(idxR)),E.VTA_binned(idxR),'normal','constant','on');
        p.sessionInfo.VTAvC(sess,:) = [bL(2) bR(2)];
        
%         yfit = glmval(bL, abs(E.cDiff(idxL)), 'identity');
%         plot( abs(E.cDiff(idxL)), E.VTA_binned(idxL), 'ko', abs(E.cDiff(idxL)), yfit, 'r-');
        
        %Calculate slope of RT vs Contrast
        bL = glmfit(abs(E.cDiff(idxL)),E.RT(idxL),'normal','constant','on');
        bR = glmfit(abs(E.cDiff(idxR)),E.RT(idxR),'normal','constant','on');
        p.sessionInfo.RTvC(sess,:) = [bL(2) bR(2)];
        
%         g = gramm('x',E.cDiff,'y',E.choice=='Right choice');
%         g.stat_summary('geom','point');
%         figure; g.draw();
%         yfit = glmval([p.sessionInfo.psych_B(sess) p.sessionInfo.psych_SL(sess) p.sessionInfo.psych_SR(sess)]', [E.contrastLeft E.contrastRight], 'logit');
%         plot(g.facet_axes_handles, E.cDiff, yfit, 'ko');
% %         
        
    end
end


clear g;
g(1,1) = gramm('x',-p.sessionInfo.psych_SL,'y',p.sessionInfo.VTAvC(:,1),'color',p.sessionInfo.mouseName);
g(1,1).geom_point();
g(1,1).set_names('x','SL parameter','y','VTA slope on CL');
g(1,2) = gramm('x',p.sessionInfo.psych_SR,'y',p.sessionInfo.VTAvC(:,2),'color',p.sessionInfo.mouseName);
g(1,2).geom_point();
g(1,2).set_names('x','SR parameter','y','VTA slope on CR');
g(1,3) = gramm('x',-p.sessionInfo.psych_SL - p.sessionInfo.psych_SR,'y',p.sessionInfo.VTAvC(:,1)-p.sessionInfo.VTAvC(:,2),'color',p.sessionInfo.mouseName);
g(1,3).geom_point();
g(1,3).set_names('x','SL-SR parameter','y','VTA_L - VTA_R');
g(1,4) = gramm('x',p.sessionInfo.psych_B,'y',p.sessionInfo.VTAvC(:,1)-p.sessionInfo.VTAvC(:,2),'color',p.sessionInfo.mouseName);
g(1,4).geom_point();
g(1,4).set_names('x','Bias parameter','y','VTA_L - VTA_R');
figure; g.draw();
g(1,1).update('color',[]);
g(1,1).stat_glm('disp_fit',true);
g(1,2).update('color',[]);
g(1,2).stat_glm('disp_fit',true);
g(1,3).update('color',[]);
g(1,3).stat_glm('disp_fit',true);
g(1,4).update('color',[]);
g(1,4).stat_glm('disp_fit',true);
g.draw();


clear g;
g(1,1) = gramm('x', -p.sessionInfo.RTvC(:,1),...
          'y', p.sessionInfo.VTAvC(:,1),...
          'color', p.sessionInfo.mouseName);
g(1,2) = gramm('x', -p.sessionInfo.RTvC(:,2),...
          'y', p.sessionInfo.VTAvC(:,2),...
          'color', p.sessionInfo.mouseName);
g(1,3) = gramm('x', p.sessionInfo.RTvC(:,1)-p.sessionInfo.RTvC(:,2),...
          'y', p.sessionInfo.VTAvC(:,2)-p.sessionInfo.VTAvC(:,1),...
          'color', p.sessionInfo.mouseName);
g(1,1).geom_point();
g(1,2).geom_point();
g(1,3).geom_point();
g(1,1).set_names('x','RT_L vs C_L slope','y','VTA_L vs C_L slope');
g(1,2).set_names('x','RT_R vs C_R slope','y','VTA_R vs C_R slope');
g(1,3).set_names('x','RT L - R slopes','y','VTA R - L slopes');
g.geom_vline('xintercept',0,'style','k:');
g.geom_hline('yintercept',0,'style','k:');
g.set_layout_options('redraw',true,'redraw_gap',0.01);
figure; g.draw();
g(1,1).update('color',[]);
g(1,1).stat_glm('disp_fit',true);
g(1,2).update('color',[]);
g(1,2).stat_glm('disp_fit',true);
g(1,3).update('color',[]);
g(1,3).stat_glm('disp_fit',true);
g.draw();


clear g;
g(1,1) = gramm('x',abs(D.cDiff),'y',D.VTA_binned,'color', D.cDiff>0, 'subset', D.cDiff~=0);
g(1,1).stat_summary();
g(1,1).set_names('x','abs(C)','y','VTA');
g(1,1).facet_grid(D.mouseName,D.sessionNum);

g(1,2) = gramm('x',abs(D.cDiff),'y',D.RT,'color', D.cDiff>0, 'subset', D.cDiff~=0);
g(1,2).stat_summary();
g(1,2).set_names('x','abs(C)','y','RT');
g(1,2).facet_grid(D.mouseName,D.sessionNum);

figure; g.draw();

%For each session compute the L-R difference in DA and RT for easy contrast
%correct trials only
[mice,~,miceID] = unique(D.mouseName);
M = {};
RT = [];
DA = [];
PERF = {};
for m = 1:max(miceID)
    E =  D(strcmp(D.mouseName,mice{m}) & abs(D.cDiff)>=0.5,:);
    E.RT( E.RT < stimWindow(2) ) = NaN;
    E.sessionStage = discretize(E.sessionNum,3);
    
    if ~all(isnan(E.VTA_binned))
        DAR=groupsummary(E.VTA_binned(E.cDiff>0),E.sessionNum(E.cDiff>0),'median');
        DAL=groupsummary(E.VTA_binned(E.cDiff<0),E.sessionNum(E.cDiff<0),'median');
%         da = (DAR - DAL);
        da = (DAR - DAL)./(mean([DAR, DAL],2) + 0.05);
        DA = [DA; da];
        
        RTR=groupsummary(E.RT(E.cDiff>0),E.sessionNum(E.cDiff>0),'median');
        RTL=groupsummary(E.RT(E.cDiff<0),E.sessionNum(E.cDiff<0),'median');
        RT = [RT; RTL-RTR];
        
        M = [M; repmat(mice(m),length(DAR),1)];
    end
%     pCorrect = groupsummary(E.feedback=='Rewarded',E.sessionNum,'mean');
%     PERF = [PERF; categorical(pCorrect>=median(pCorrect),[0 1],{'Low performance','High performance'})];
end


g = gramm('x',RT,'y',DA,'color',M);
% g.facet_grid([],M');
g.set_layout_options('redraw',true,'redraw_gap',0.01);
% g.stat_glm('disp_fit',true);c
% g.stat_ellipse('type','95percentile');
g.geom_point();
g.geom_vline('xintercept',0,'style','k:');
g.geom_hline('yintercept',0,'style','k:');
g.set_names('x','\Delta RT (Right choice faster)','y','\Delta DA (Right stim larger)');
figure; g.draw();
g.update('color',[]);
g.stat_glm('disp_fit',true);
g.draw();

%% Plot timewarped DA traces in DMS, focusing on sessions with and without behavioural history effects
% 
% %Mouse name, sessWithHist, sessWithoutHist effects
% sess = {'ALK074',[4 6],[7 9 10 11];
%         'ALK083',[1 2 3 7],[5 9 10 12];
%         'ALK084',[1 2 3 5],[6 7 8];
%         'MMM003',[2 4 5 6],[3 9 10];
%         'MMM009',[1 2 3 5],[4 6 7 8];
%         'MMM010',[1 2],[3 4]};
% %go through and assign sessions
% D.histEffectLabel = nan(height(D),1);
% for s = 1:size(sess,1)
%     idx = strcmp(D.mouseName,sess{s,1}) & ismember(D.sessionNum, sess{s,2} );
%     D.histEffectLabel(idx) = 1;
%     
%     idx = strcmp(D.mouseName,sess{s,1}) & ismember(D.sessionNum, sess{s,3} );
%     D.histEffectLabel(idx) = 0;
% end
% D.histEffectLabel = categorical(D.histEffectLabel,[0 1],{'NoHistEffect','HistEffect'});
%         

%For each region separately
for a = 3:4
    
    %get list of sessions which have recording in this region
    idx = strcmp(p.sessionInfo.photometryChannel2, regions{a}) | strcmp(p.sessionInfo.photometryChannel4, regions{a});
    eRefs = p.sessionInfo.expRef(idx);
    
    E = D(ismember(D.expRef,eRefs),:);
    
    E.earlyLate = categorical(E.sessionNum>4,[0 1],{'Sess 1-4','Sess 5+'});
    
    
    %remove undefined sessions
%     E(isundefined( E.histEffectLabel),:) = [];
    
    %time warped trace
    clear g;
    g(1,1) = gramm('x',1:sum(warp_sizes),'y',E.dff_timewarped(:,:,a),'color',E.cDiff,'subset',E.feedback=='Rewarded');
    g(1,1).facet_grid(E.earlyLate,E.mouseName,'scale','independent');
    g(1,1).stat_summary('setylim','true');
    g(1,1).set_color_options('map',0.9*RedWhiteBlue(floor(length(unique(E.cDiff))/2)));
    g(1,1).geom_vline('xintercept',cumsum(warp_sizes(1:3)),'style','k:');
    g(1,1).set_title('Rewarded trials separated by contrast');
    g(1,1).set_names('x','','y','zDFF','color','contrast','column','','row','');
    g(1,1).axe_property('xtick','','xcolor','none');

    g(2,1) = gramm('x',E.cDiff,'y',double(E.choice=='Right choice'),'color',circshift(E.choice,1),'subset',circshift(E.feedback,1)=='Rewarded');
    g(2,1).facet_grid(E.earlyLate,E.mouseName,'scale','independent');
    g(2,1).stat_summary('geom',{'point','line','errorbar'},'type',@grammCallbackBinomialConfidenceInterval);
    g(2,1).axe_property('ylim',[0 1]);
    g(2,1).geom_vline('xintercept',0,'style','k:');
    g(2,1).geom_hline('yintercept',0.5,'style','k:');
    g(2,1).set_names('x','C','y','pR','column','','row','','color','');

    g.set_layout_options('redraw',true,'redraw_gap',0);
    g.set_title(regions{a});
    
    figure('units','normalized','outerposition',[0.0042 0.0074 0.5380 0.9139],'color','w');
    g.draw();
end
%% Fit and plot kernels to stimulus/choice/reward

p = loadProjectDataset('GCaMP_LearningGrating2AFC');


%Remove trials with RT>3
D = D(D.RT<3 & D.RT>0,:);

% Create spline basis set
numSplines = 9;
% splineBasis = create_bspline_basis([-0.2 4],numSplines);
% splineBasisMat = full(eval_basis( t_sample_wholetrial ,splineBasis));
splineBasis = create_bspline_basis([0 2],numSplines);
splineBasisMat = full(eval_basis( 0:mean(diff(t_sample_wholetrial)):2 ,splineBasis));
numSplineTimepoints = size(splineBasisMat,1);

%Create boxcar windows for LeftStimulus events and convolve with splines
D.stimLeftEvent = zeros(height(D), length(t_sample_wholetrial));
D.stimRightEvent = zeros(height(D), length(t_sample_wholetrial));
D.choiceLeftEvent = zeros(height(D), length(t_sample_wholetrial));
D.choiceRightEvent = zeros(height(D), length(t_sample_wholetrial));
D.rewardEvent = zeros(height(D), length(t_sample_wholetrial));

D.stimLeftConvSpline = nan( height(D) , length(t_sample_wholetrial), numSplines );
D.stimRightConvSpline = nan( height(D) , length(t_sample_wholetrial), numSplines );
D.choiceLeftConvSpline = nan( height(D) , length(t_sample_wholetrial), numSplines );
D.choiceRightConvSpline = nan( height(D) , length(t_sample_wholetrial), numSplines );
D.rewardConvSpline = nan( height(D) , length(t_sample_wholetrial), numSplines );

idx = find( t_sample_wholetrial > 0,1,'first'); %stimulus time
D.stimLeftEvent(D.contrastLeft>0,idx)=1;
D.stimRightEvent(D.contrastRight>0,idx)=1;
[~,idx]=max(t_sample_wholetrial > D.RT,[],2); %choice time
D.choiceLeftEvent(sub2ind( size(D.choiceLeftEvent),find(D.choice=='Left choice'),idx(D.choice=='Left choice')))=1;
D.choiceRightEvent(sub2ind( size(D.choiceRightEvent),find(D.choice=='Right choice'),idx(D.choice=='Right choice')))=1;
[~,idx]=max(t_sample_wholetrial > (D.outcomeTime - D.stimulusOnsetTime),[],2); %reward time
D.rewardEvent(sub2ind( size(D.rewardEvent),find(D.feedback=='Rewarded'),idx(D.feedback=='Rewarded')))=1;

%now convolve with splines
for sp = 1:numSplines
    w = conv2(D.stimLeftEvent,splineBasisMat(:,sp)');
    D.stimLeftConvSpline(:, :, sp) = w(:,1:length(t_sample_wholetrial));
    
    w = conv2(D.stimRightEvent,splineBasisMat(:,sp)');
    D.stimRightConvSpline(:, :, sp) = w(:,1:length(t_sample_wholetrial));
    
    w = conv2(D.choiceLeftEvent,splineBasisMat(:,sp)');
    D.choiceLeftConvSpline(:, :, sp) = w(:,1:length(t_sample_wholetrial));
    
    w = conv2(D.choiceRightEvent,splineBasisMat(:,sp)');
    D.choiceRightConvSpline(:, :, sp) = w(:,1:length(t_sample_wholetrial));
    
    w = conv2(D.rewardEvent,splineBasisMat(:,sp)');
    D.rewardConvSpline(:, :, sp) = w(:,1:length(t_sample_wholetrial));
end

kernels = table;

kernelLabel = {'contrastLeft Kernel','contrastRight Kernel','choice left Kernel','choice right Kernel','Reward Kernel'};
%for each region
for a = 1:length(regions)
        
    %Get activity in one region
    Y = D.dff_wholetrial(:,:,a);
    
    %Separate for each session
    sess = unique(D.expRef(~isnan(Y(:,1))));
    for s = 1:length(sess)
        idx = strcmp(D.expRef,sess{s});
        Ysess = Y(idx,:);
        SL = reshape(D.stimLeftConvSpline(idx,:,:), numel(Ysess), numSplines);
        SR = reshape(D.stimRightConvSpline(idx,:,:), numel(Ysess), numSplines);
        ChL = reshape(D.choiceLeftConvSpline(idx,:,:), numel(Ysess), numSplines);
        ChR = reshape(D.choiceRightConvSpline(idx,:,:), numel(Ysess), numSplines);
        R = reshape(D.rewardConvSpline(idx,:,:), numel(Ysess), numSplines);
        X = [SL, SR, ChL, ChR, R];
        coeff = glmfit(X,Ysess(:),'normal','constant','off');

        for evts = 1:length(kernelLabel)
            K = splineBasisMat*coeff( (1:numSplines) + (evts-1)*numSplines );
            
            row = table;
            row.mouseName = unique(D.mouseName(idx));
            row.sessionNumber = unique(D.sessionNum(idx));
            row.kernelType = kernelLabel(evts);
            row.kernelValue = K';
            row.region = regions(a);
            kernels = [kernels; row];
        end
    
    end
end

g = gramm('x', 0:mean(diff(t_sample_wholetrial)):2, 'y', kernels.kernelValue);
g.facet_grid(kernels.region, kernels.kernelType,'scale','free_y');
g.geom_hline('yintercept',0,'style','k:');
g.geom_line();
g.set_names('x','Time from event (sec)','y','Kernel value','row','','column','');
figure; g.draw();
g.update();
g.set_color_options('chroma',0,'lightness',30);
g.stat_summary('setylim',true);
g.draw();


%% Create plots for presentation

%Plot simulations with different levels of sensitivity
figure;
CL = [linspace(1,0,50)'; zeros(50,1)];
CR = [zeros(50,1); linspace(0,1,50)'];

subplot(1,3,1);
yfit = glmval([0 -10 3]', [CL CR], 'logit');
plot(CR-CL,yfit);
xline(0); yline(0.5);

subplot(1,3,2);
yfit = glmval([0 -3 3]', [CL CR], 'logit');
plot(CR-CL,yfit);
xline(0); yline(0.5);

subplot(1,3,3);
yfit = glmval([0 -3 10]', [CL CR], 'logit');
plot(CR-CL,yfit);
xline(0); yline(0.5);




%Last several sessions of ALK074: psych and RT plots
E=D(strcmp(D.mouseName,'ALK074') & D.sessionNum>15,:);
clear g;
g(1,1) = gramm('x',E.cDiff,'y',double(E.choice=='Right choice'));
% g(1,1) = gramm('x',E.cDiff,'y',double(E.feedback=='Rewarded'));
g(1,1).stat_summary('geom',{'point','line','errorbar'},'type',@grammCallbackBinomialConfidenceInterval);
g(1,1).axe_property('ylim',[0 1]);
g(1,1).geom_vline('xintercept',0,'style','k:');
g(1,1).geom_hline('yintercept',0.5,'style','k:');
g(1,1).set_names('x','Contrast','y','p(Right) and 95% Binomial CI','column','','row','');
g(1,1).set_title('Psychometric curve');

g(2,1) = gramm('x',E.cDiff,'y',E.RT,'subset',E.feedback=='Rewarded');
% g(2,1) = gramm('x',E.cDiff,'y',E.RT,'color',E.choice);
g(2,1).stat_summary('geom',{'point','line','errorbar'},'setylim',true,'type',@(y) [median(y); median(y)-mad(y,1); median(y)+mad(y,1)]);
g(2,1).geom_vline('xintercept',0,'style','k:');
g(2,1).set_names('x','Contrast','y','RT median and MAD','column','','row','');
g(2,1).set_title('Reaction time (rewarded trials)');

g.axe_property('xtick',unique(E.cDiff));

g.set_layout_options('redraw',true,'redraw_gap',0.001);

figure('color','w')
g.draw();


% Plot progress of trainng across all mice (performance over days);
clear g;
%Perf (high C) over days
g(1,1) = gramm('x',D.sessionNum,'y',D.feedback=='Rewarded','group',D.mouseName,'subset',abs(D.cDiff)>=0.5);
g(1,1).stat_summary('geom',{'line','point'});
g(1,1).geom_hline('yintercept',0.5,'style','k:');
g(1,1).set_names('x','Days','y','Performance on easy contrasts','column','','row','');
g(1,1).axe_property('ylim',[0 1]);

%RT over days
g(2,1) = gramm('x',D.sessionNum,'y',D.RT,'group',D.mouseName,'subset',abs(D.cDiff)>=0.5 & D.feedback=='Rewarded');
% g(2,1) = gramm('x',D.sessionNum,'y',D.RT,'lightness',D.mouseName);
g(2,1).stat_summary('geom',{'line','point'},'type','quartile');
g(2,1).set_names('x','Days','y','median RT','column','','row','');


%L-R performance
idxR = D.cDiff>=0.5;
idxL = D.cDiff<=0.5;
[perfR,perflabel] = groupsummary(D.feedback(idxR)=='Rewarded',{D.sessionNum(idxR) D.mouseName(idxR)},'mean');
[perfL] = groupsummary(D.feedback(idxL)=='Rewarded',{D.sessionNum(idxL) D.mouseName(idxL)},'mean');
perf = abs(perfL - perfR);
g(1,2) = gramm('x',perflabel{1},'y',perf,'group',perflabel{2});
g(1,2).stat_summary('geom',{'line','point'});
g(1,2).set_names('x','Days','y','|Acc_L - Acc_R|');
g(1,2).geom_hline('yintercept',0,'style','k:');

%L-R RT
idxR = D.cDiff>=0.5 & D.feedback=='Rewarded';
idxL = D.cDiff<=0.5 & D.feedback=='Rewarded';
[rtR,rtlabel] = groupsummary(D.RT(idxR),{D.sessionNum(idxR) D.mouseName(idxR)},'median','IncludeEmptyGroups',true);
[rtL] = groupsummary(D.RT(idxL),{D.sessionNum(idxL) D.mouseName(idxL)},'median','IncludeEmptyGroups',true);
rt = abs(rtL - rtR);
g(2,2) = gramm('x',rtlabel{1},'y',rt,'group',rtlabel{2});
g(2,2).stat_summary('geom',{'line','point'});
g(2,2).set_names('x','Days','y','|RT_L - RT_R|');
g(2,2).geom_hline('yintercept',0,'style','k:');
% g.axe_property('xlim',[0 12]);
figure; g.draw();

%add in average over mice
idx = abs(D.cDiff)>=0.5;
[perf1,perf1label] = groupsummary(D.feedback(idx)=='Rewarded',{D.sessionNum(idx) D.mouseName(idx)},'mean');
plot( g(1,1).facet_axes_handles, groupsummary(perf1,perf1label{1},'median'), 'k-','linewidth',4);

[rt1,rt1label] = groupsummary(D.RT(idx),{D.sessionNum(idx) D.mouseName(idx)},'median');
plot( g(2,1).facet_axes_handles, groupsummary(rt1,rt1label{1},'median'), 'k-','linewidth',4);


% plot( g(1,2).facet_axes_handles, groupsummary(perf,perflabel{1},'mean'), 'k-','linewidth',4);
% plot( g(2,2).facet_axes_handles, groupsummary(rt,rtlabel{1},'mean'), 'k-','linewidth',4);
plot( g(1,2).facet_axes_handles, groupsummary(perf,perflabel{1},'median'), 'k-','linewidth',4);
plot( g(2,2).facet_axes_handles, groupsummary(rt,rtlabel{1},'median'), 'k-','linewidth',4);



% Psych SL and SR
g(1,1) = gramm('x',p.sessionInfo.sessionNum,'y',-p.sessionInfo.psych_SL,'lightness',p.sessionInfo.mouseName); 
g(1,1).geom_line(); 
g(1,1).geom_point();
g(1,1).set_names('x','days','y','SL parameter');
g(1,2) = gramm('x',p.sessionInfo.sessionNum,'y',p.sessionInfo.psych_SR,'lightness',p.sessionInfo.mouseName); 
g(1,2).geom_line(); 
g(1,2).geom_point();
g(1,2).set_names('x','days','y','SR parameter');
g(1,3) = gramm('x',p.sessionInfo.sessionNum,'y',-p.sessionInfo.psych_SL - p.sessionInfo.psych_SR,'lightness',p.sessionInfo.mouseName); 
g(1,3).geom_line(); 
g(1,3).geom_point(); 
g(1,3).set_names('x','days','y','SL-SR parameter');
g.set_layout_options('redraw',true,'redraw_gap',0.01); 
g.axe_property('xlim',[0 12]);
figure; g.draw()

% Plot progress of DA_L - DA_R over days

% Define contra and ipsi VTA responses
E = D( contains(D.photometryChannel2,'VTA'),:);
E.VTA_binned = nanmean( E.dff_stim_binned(:,:,contains(regions,'VTA')), 3);
E.cDiff_contraIpsi = E.cDiff; %positive = ipsi, negative = contra
E.cDiff_contraIpsi( strcmp(E.photometryChannel2,'Left VTA') ) = -E.cDiff_contraIpsi( strcmp(E.photometryChannel2,'Left VTA') ) ;
clear g;
g(1,1) = gramm('x',E.sessionNum,'y',E.VTA_binned,'lightness',E.mouseName,'subset',abs(E.cDiff)>=0.5);
g(1,1).stat_summary('geom',{'line','point','errorbar'});
g(1,1).geom_hline('yintercept',0,'style','k:');
g(1,1).set_names('x','Days','y','dF/F (z-score) post-stim','column','','row','');

idxR = E.cDiff>=0.5 & E.feedback=='Rewarded';
idxL = E.cDiff<=0.5 & E.feedback=='Rewarded';
[daR,dalabel] = groupsummary(E.VTA_binned(idxR),{E.sessionNum(idxR) E.mouseName(idxR)},'mean');
[daL] = groupsummary(E.VTA_binned(idxL),{E.sessionNum(idxL) E.mouseName(idxL)},'mean');
da = (daL - daR)./(mean([daL daR],2)+0.05); %normalise by mean response each day
da = abs(da);
g(2,1) = gramm('x',dalabel{1},'y',da,'lightness',dalabel{2});
g(2,1).stat_summary('geom',{'line','point'});
g(2,1).set_names('x','Days','y','|DA_L - DA_R|');
g(2,1).geom_hline('yintercept',0,'style','k:');
g.axe_property('xlim',[0 12]);
figure; g.draw();

idx = abs(E.cDiff)>=0.5;
[da1,da1label] = groupsummary(E.VTA_binned(idx),{E.sessionNum(idx) E.mouseName(idx)},'mean');
plot( g(1,1).facet_axes_handles, 1:13, groupsummary(da1,da1label{1},'mean'), 'k-','linewidth',4);
plot( g(2,1).facet_axes_handles, 1:13, groupsummary(da,dalabel{1},'mean'), 'k-','linewidth',4);


% tuning curves
clear g;
g = gramm('x',E.cDiff,'y',E.VTA_binned,'subset',E.feedback=='Rewarded','lightness',E.mouseName);
g.facet_grid(E.mouseName, E.sessionNum,'scale','free_y');
g.stat_summary('geom',{'line','point','errorbar'},'setylim',true);
g.set_names('x','Contrast','y','dF/F (z-score) post-stim','row','','column','');
g.set_layout_options('redraw',true,'redraw_gap',0.01);
figure; g.draw();


% tuning curves with slope per side
clear g;
E = E(E.cDiff~=0,:);
g = gramm('x',abs(E.cDiff),'y',E.VTA_binned,'subset',E.feedback=='Rewarded','color',sign(E.cDiff));
g.facet_grid(E.mouseName, E.sessionNum,'scale','free_y');
% g.stat_summary('geom',{'line','point','errorbar'},'setylim',true);
g.stat_summary('geom',{'point','errorbar'},'setylim',true);
g.stat_glm('disp_fit',false);
g.set_names('x','ABS(Contrast)','y','dF/F (z-score) post-stim','row','','column','','color','Stimulus side');
g.set_layout_options('redraw',true,'redraw_gap',0.01);
figure; g.draw();

%% Fit history model 

d = struct;
d.contrastLeft = D.contrastLeft;
d.contrastRight = D.contrastRight;
d.choice = double(D.choice=='Right choice');
[mice,~,d.subjectID]=unique(D.mouseName);
[~,~,d.sessionID] = unique([d.subjectID D.sessionNum],'rows');

%currentWinEffect
win = zeros(size(D.choice));
win(D.feedback=='Rewarded' & D.choice=='Right choice') = +1;
win(D.feedback=='Rewarded' & D.choice=='Left choice') = -1;

%currentLoseEffect
lose = zeros(size(D.choice));
lose(D.feedback=='Unrewarded' & D.choice=='Right choice') = +1;
lose(D.feedback=='Unrewarded' & D.choice=='Left choice') = -1;

d.prevWin = circshift(win,-1); % -1=L, +1=R previous trial rewarded
d.prevLose = circshift(lose,-1); % -1=L, +1=R previous trial rewarded

%Remove first trials
d = structfun(@(f) f( D.trialNumber>1, :) , d,  'UniformOutput', 0);

%fit
fit = stan_fitModel('Hierarchical_Logistic_OutcomeHistory',d,'\\QNAP-AL001.dpag.ox.ac.uk\PZatka-Haas\GCaMPLearning\data\fit_model.mat');
% fit = load("\\QNAP-AL001.dpag.ox.ac.uk\PZatka-Haas\GCaMPLearning\data\fit_model.mat");

