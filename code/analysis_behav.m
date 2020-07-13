
%% Load IBL dataset
D = loadIBLdataset();
D.cDiff = D.contrastRight - D.contrastLeft;

%% Plot performance, RT, and behav asymmetries over days

clear g;
%Performance at high contrast
g(1,1) = gramm('x',D.sessionNum,'y',D.feedback=='Rewarded','group',D.mouseName,'subset',abs(D.cDiff)>=0.5);
g(1,1).stat_summary('geom',{'line','point'});
g(1,1).geom_hline('yintercept',0.5,'style','k:');
g(1,1).set_names('x','Days','y','Performance on easy contrasts','column','','row','');
g(1,1).axe_property('ylim',[0 1]);

%RT at high contrast for correct trials
g(2,1) = gramm('x',D.sessionNum,'y',D.RT,'group',D.mouseName,'subset',abs(D.cDiff)>=0.5 & D.feedback=='Rewarded');
g(2,1).stat_summary('geom',{'line','point'},'type','quartile');
g(2,1).set_names('x','Days','y','median RT','column','','row','');

%Performance asymmetry
idxR = D.cDiff>=0.5;
idxL = D.cDiff<=0.5;
[perfR,perflabel] = groupsummary(D.feedback(idxR)=='Rewarded',{D.sessionNum(idxR) D.mouseName(idxR)},'mean','IncludeEmptyGroups',true);
[perfL] = groupsummary(D.feedback(idxL)=='Rewarded',{D.sessionNum(idxL) D.mouseName(idxL)},'mean','IncludeEmptyGroups',true);
perf = abs(perfL - perfR)./mean([perfL perfR],2);
g(1,2) = gramm('x',perflabel{1},'y',perf,'group',perflabel{2});
g(1,2).stat_summary('geom',{'line','point'});
g(1,2).set_names('x','Days','y','|Acc_L - Acc_R|/avg(Perf)');
g(1,2).geom_hline('yintercept',0,'style','k:');

% RT asymmetry
idxR = D.cDiff>=0.5 & D.feedback=='Rewarded';
idxL = D.cDiff<=0.5 & D.feedback=='Rewarded';
[rtR,rtlabel] = groupsummary(D.RT(idxR),{D.sessionNum(idxR) D.mouseName(idxR)},'median','IncludeEmptyGroups',true);
[rtL] = groupsummary(D.RT(idxL),{D.sessionNum(idxL) D.mouseName(idxL)},'median','IncludeEmptyGroups',true);
rt = abs(rtL - rtR);
g(2,2) = gramm('x',rtlabel{1},'y',rt,'group',rtlabel{2});
g(2,2).stat_summary('geom',{'line','point'});
g(2,2).set_names('x','Days','y','|RT_L - RT_R|');
g(2,2).geom_hline('yintercept',0,'style','k:');
figure; g.draw();

%add in average over mice
idx = abs(D.cDiff)>=0.5;
[perf1,perf1label] = groupsummary(D.feedback(idx)=='Rewarded',{D.sessionNum(idx) D.mouseName(idx)},'mean');
g(1,1).update('x',  perf1label{1}, 'y', perf1);
g(1,1).stat_summary('type','95percentile');
g(1,1).set_color_options('map',[0 0 0]);
g(1,1).draw();

[rt1,rt1label] = groupsummary(D.RT(idx),{D.sessionNum(idx) D.mouseName(idx)},'median');
g(2,1).update('x',  rt1label{1}, 'y', rt1);
g(2,1).stat_summary('type','95percentile');
g(2,1).set_color_options('map',[0 0 0]);
g(2,1).draw();


g(1,2).update('x', perflabel{1}, 'y', perf);
g(1,2).stat_summary('type','95percentile');
g(1,2).set_color_options('map',[0 0 0]);
g(1,2).draw();

g(2,2).update('x', rtlabel{1}, 'y', rt);
g(2,2).stat_summary('type','95percentile');
g(2,2).set_color_options('map',[0 0 0]);
g(2,2).draw();

% 
% plot( g(1,2).facet_axes_handles, groupsummary(perf,perflabel{1},'median'), 'k-','linewidth',4);
% plot( g(2,2).facet_axes_handles, groupsummary(rt,rtlabel{1},'median'), 'k-','linewidth',4);
