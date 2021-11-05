% Script for running the supplementary analysis of the AVSRT-Autism Study.

%   Dependencies:
%       RaceModel: https://github.com/mickcrosse/RaceModel
%       PERMUTOOLS: https://github.com/mickcrosse/PERMUTOOLS
%       RSE-box: https://github.com/tomotto/RSE-box
%       M3 Toolbox: https://github.com/canlab/MediationToolbox
%       MES: https://github.com/hhentschke/measures-of-effect-size-toolbox

%   References:
%      [1] Crosse MJ, Foxe JJ, Molholm S (2019) Developmental Recovery of 
%          Impaired Multisensory Processing in Autism and the Cost of 
%          Switching Sensory Modality. bioRxiv 565333.
%      [2] Crosse MJ, Foxe JJ, Tarrit K, Freedman EG, Molholm S (submitted) 
%          Resolution of Impaired Multisensory Processing in Autism and the 
%          Cost of Switching Sensory Modality.

% Copyright (c) 2021, Mick Crosse <mickcrosse@gmail.com>

%% Fig. S1 | Methods

close all
clc

% Set parameters
cond = 5;
task = 1;
clear pl

rmvPos = zeros(nSubjs,1);
rmvNeg = zeros(nSubjs,1);
for i = 1:nSubjs
    idx = squeeze(CDFs(i,cond,task,:)>0);
    rmvPos(i) = sum(CDFs(i,cond,task,idx));
    idx = squeeze(CDFs(i,cond,task,:)<0);
    rmvNeg(i) = sum(CDFs(i,cond,task,idx));
end

% Define parameters
xaxis = prob/100;

nBins = 20;

xvals = linspace(perLim(exSub,1),perLim(exSub,2),21);

% Generate figure
axes('parent',figure,'fontname',fnt,'fontsize',fntSz)
set(gcf,'position',[20,50,550,650],'paperpositionmode','auto')

subplot(3,1,1)
hold all
x1 = 100; x2 = 700;
y1 = 0; y2 = 0.2;
for i = 1:2
    plot([perLim(exSub,i),perLim(exSub,i)],[y1,y2],'--k')
end
for i = 2:20
    plot([xvals(i),xvals(i)],[y1,y2],'-','color',[0.75,0.75,0.75])
end
for i = 1:3
    pl(i) = histogram(spcloneRTs(spcloneConds==i),nBins,'facecolor',cols1(i,:),...
        'facealpha',1,'Normalization','probability');
end
set(gca,'layer','top','tickdir','out')
xlim([x1,x2])
ylim([y1,y2])
% text(x1-150,y2+8,'C','fontname',fnt,'fontsize',fntSz+8)
xlabel('Reaction time (ms)','fontname',fnt,'fontsize',fntSz)
ylabel('Probability','fontname',fnt,'fontsize',fntSz)
l = legend(pl,'AV','A','V');
p = get(l,'position');
set(l,'box','off','fontname',fnt,'fontsize',fntSz-2,'position',[p(1)-0.22,p(2)+0.025,p(3),p(4)])
axis square

subplot(3,1,2)
hold all
x1 = 0; x2 = 1;
y1 = 0; y2 = 1;
% for i = 1:20
%     plot([i/20,i/20],[y1+0.005,y2],'-','color',[0.75,0.75,0.75])
% end
for i = 1:3
    pl(i) = plot(xaxis,squeeze(CDFs(exSub,i,1,:)),'-','linewidth',...
        lnWth,'color',cols3(i,:),'markersize',3,'markerfacecolor',cols1(i,:));
end
pl(4) = plot(xaxis,squeeze(CDFs(exSub,6,1,:)),':','linewidth',...
    lnWth,'color',[0,0,0],'markersize',3,'markerfacecolor',cols1(i,:));
set(gca,'tickdir','out')
xlim([x1,x2])
ylim([y1,y2])
% text(x1-0.44,y2+0.1,'D','fontname',fnt,'fontsize',fntSz+8)
xlabel('Quantile','fontname',fnt,'fontsize',fntSz)
ylabel({'Cumulative';'probability'},'fontname',fnt,'fontsize',fntSz)
l = legend(pl,'AV','A','V','Race'); p = get(l,'position');
% set(l,'box','off','fontname',fnt,'fontsize',fntSz-2,'position',[p(1)-0.23,p(2)-0.1,p(3),p(4)])
set(l,'box','off','fontname',fnt,'fontsize',fntSz-2,'position',[p(1)-0.23,p(2)-0.09,p(3),p(4)])
axis square

h = subplot(3,1,3);
hold all
x1 = 0; x2 = 0.6;
y1 = 0; y2 = 1;
pl(2) = plot(xaxis(1:8),squeeze(mean(CDFs(exSub,8,1,1:8),3)),'-','color',[0.6,0.6,0.6],'linewidth',lnWth);
pl(1) = plot(xaxis(1:1:end),squeeze(mean(CDFs(exSub,4,6:8,1:1:end),3)),'x','linewidth',1.5,'markersize',8,'color',cols1(2,:),'markerfacecolor',cols1(2,:));
pl(3) = plot(xaxis,squeeze(mean(CDFs(exSub,6,6:8,:),3)),'--k','linewidth',lnWth);
set(gca,'tickdir','out')
xlim([x1,x2])
ylim([y1,y2])
% text(x1-0.44,y2+0.1,'E','fontname',fnt,'fontsize',fntSz+8)
xlabel('Quantile','fontname',fnt,'fontsize',fntSz)
ylabel({'Cumulative';'probability'},'fontname',fnt,'fontsize',fntSz)
l = legend(pl,'Sim.','Miller','Raab'); p = get(l,'position');
posnew = p; posnew(1) = posnew(1)-0.23; posnew(2) = posnew(2)-0.1;
set(l,'box','off','fontname',fnt,'fontsize',fntSz-2,'position',posnew)
axis square

set(gcf,'color','w')


%% Fig. S2 | Percentage subjects with multisensory facilitation

close all

perRMV = zeros(2,4,7);

for i = 1:2
    for j = 1:4
        idx = grp_vec==i & dev_vec==j;
        perRMV(i,j,:) = sum(rmv(idx,2:8)>0)/sum(idx)*100;
    end
end

figure,
subplot(2,1,1)
bar(squeeze(perRMV(1,:,:))','edgecolor','none')
ylabel({'AV > Race';'(% subjects)'})
title('NT')
text(-0.5,100,'a','fontsize',12,'fontweight','bold')

subplot(2,1,2)
bar(squeeze(perRMV(2,:,:))','edgecolor','none')
l = legend(age_groups);
set(l,'box','off')
xlabel('Quantile')
ylabel({'AV > Race';'(% subjects)'})
title('ASD')
text(-0.5,100,'b','fontsize',12,'fontweight','bold')

set(gcf,'color','w')

%% Fig. S3 | Group comparisons in Gain

% close all;
clc;

% Set parameters
nperm = 1e4;
lnWth = 2;
fig_label = true;
age_match = true;

% Compute RMV moving mean and run perm test
ageRng = 6:26;
window = 7;
x1 = ageRng(1); x2 = ageRng(end);

% Preallocate memory
g = zeros(1,numel(ageRng));
ci = zeros(2,numel(ageRng));
bf10 = zeros(1,numel(ageRng));

% Loop through age bins
for i = 1:numel(ageRng)
    
    switch age_match
        case true
            if i>13
                [idx1,idx2] = movmeanmatch(i,grp_vec,sex_vec,age_vec,[],ageRng,window,false);
            else
                [idx1,idx2] = movmeanmatch(i,grp_vec,sex_vec,age_vec,piq_vec,ageRng,window,false);
            end
        case false
            idx1 = find(grp_vec==1 & age_vec>ageRng(i)-window/2 & age_vec<ageRng(i)+window/2);
            idx2 = find(grp_vec==2 & age_vec>ageRng(i)-window/2 & age_vec<ageRng(i)+window/2);
    end
    
    gainNT = gain(idx1);
    gainASD = gain(idx2);
            
    % Effect size
    d = mes(gainNT,gainASD,'hedgesg','isDep',0,'nBoot',nperm);
    g(i) = d.hedgesg;
    ci(:,i) = d.hedgesgCi;

    % Bayes factor analysis
    bf10(i) = 1/bf.ttest2(gainNT,gainASD);
    
end

% Generate figure
axes('parent',figure,'fontname',fnt,'fontsize',fntSz)
set(gcf,'position',[20,50,500,450],'paperpositionmode','auto')

y1 = 0; y2 = 1.2;
plt = subplot(2,8,2:7);
hold on
plot(ageRng,g,'color',cols1(5,:),'linewidth',lnWth);
plot(ageRng,ci,':','color',cols1(5,:),'linewidth',lnWth);
title({'Moving Gain';'Group Comparisons'},'fontname',fnt,'fontsize',fntSz,'fontweight','normal')
ylabel({'Effect size';'(Hedge''s \it{g})'},'fontname',fnt,'fontsize',fntSz)
set(gca,'layer','bottom','box','off','tickdir','out',...
    'xtick',6:2:26,'ytick',0:0.2:1.2)
xlim([x1,x2])
ylim([y1,y2])
grid on
if fig_label
    text(x1-(x2-x1)*0.15,y2,'a','fontname',fnt,'fontsize',fntSz+2,...
        'fontweight','bold','verticalalign','bottom')
end

y1 = 0; y2 = 5;
plt = subplot(2,8,10:15); 
hold on
plot(ageRng,bf10,'color',cols1(6,:),'linewidth',lnWth);
xlabel('Age (yrs)','fontname',fnt,'fontsize',fntSz)
ylabel({'Bayes factor';'(BF_0_1)'},'fontname',fnt,'fontsize',fntSz)
set(gca,'layer','bottom','box','off','tickdir','out',...
    'xtick',6:2:26,'ytick',0:y2)
xlim([x1,x2])
ylim([y1,y2])
grid on
if fig_label
    text(x1-(x2-x1)*0.15,y2,'b','fontname',fnt,'fontsize',fntSz+2,...
        'fontweight','bold','verticalalign','bottom')
end

set(gcf,'color','w')

%% Fig. S4a | Modelling multisensory behavior Bayes factors

close all;
clc;

% Generate figure
axes('parent',figure,'fontname',fnt,'fontsize',fntSz)
set(gcf,'position',[20,50,650,600],'paperpositionmode','auto')

% Panel C & D
models = {'Model 1A','Model 1V','Model 2A','Model 2V'};
labels = {'a','b'};
sp = [1,2,5,6;3,4,7,8]+8;
x1 = 0; x2 = 6;
y1 = 0; y2 = 10;
ytick = 2:2:10;
ctr = 1;
dev1 = 1;
dev2 = 2;
for n = 1:2
    
    for i = 1:numel(models)
        
        rVals = zeros(5,2); 
        pVals = zeros(5,2);
        bf01 = zeros(5,2);
        
        for j = 1:nRatios
            
            % Deine predicted and empirical benefits
            benMod = simCPs(:,i,j,:)-max(CDFs(:,2:3,1,:),[],2);
            benMod = squeeze(trapz(prob,benMod,4));
            
            for m = 1:2
                x = benMod(grp_vec==m & dev_vec==dev1);
                y = benEmp(grp_vec==m & dev_vec==dev1);
                [rVals(j,m),stats] = permtest_corr(x,y,'type',type,'nperm',nperm);
                bf01(j,m) = 1/bf.corr(x,y);
                pVals(j,m) = stats.p;
            end
            
        end
        
        h = subplot(4,4,sp(n,i));
        hold all
        plot([0,6],[0,0],':k')
        pl = plot(1:5,bf01,'-','linewidth',2.5);
        for j = 1:nRatios
            if pVals(j,1)<0.05
                text(j,rVals(j,1)+0.02,'*','fontname',fnt,'fontsize',fntSz+5,'horizontalalign','center')
            end
            if pVals(j,2)<0.05
                text(j,rVals(j,2)+0.02,'*','fontname',fnt,'fontsize',fntSz+5,'horizontalalign','center')
            end
        end
        set(gca,'tickdir','out','xtick',1:5,'xticklabel',{'0','0.25','0.50','0.75','1.00'},...
            'fontname',fnt,'fontsize',fntSz)
        ylim([y1,y2])
        xlim([x1,x2])
        xtickangle(45)
        title(models{i},'fontname',fnt,'fontsize',fntSz)
        if n==2
            set(gca,'ytick',ytick)
            p = get(h,'position');
            set(h,'position',[p(1)+0.04,p(2),p(3),p(4)])
        else
            set(gca,'ytick',ytick)
            p = get(h,'position');
            set(h,'position',[p(1)-0.02,p(2),p(3),p(4)])
        end
        if i==1
            if fig_label
                text(x1-(x2-x1)*0.5,y2,labels{n},'fontname',fnt,'fontsize',fntSz+2,...
                    'fontweight','bold','verticalalign','bottom')
            end
            text(7,y2+0.25,[age_groups{dev1},' yrs'],'horizontalalign','center',...
                'fontweight','normal','fontname',fnt,'fontsize',fntSz)
        end
        if i<3
            set(gca,'xticklabel',[])
            p = get(h,'position');
            set(h,'position',[p(1),p(2)+0.02,p(3),p(4)])
        else
            p = get(h,'position');
            set(h,'position',[p(1),p(2)+0.02,p(3),p(4)])
        end
        if i==2 || i==4
            set(gca,'yticklabel',[])
            p = get(h,'position');
            set(h,'position',[p(1)-0.02,p(2),p(3),p(4)])
        elseif i==3
            p = get(h,'position');
            set(h,'position',[p(1),p(2),p(3),p(4)])
            xlabel('Race\leftarrow{\itp}\rightarrowBias','fontname',fnt,'fontsize',fntSz)
            ylabel('Correlation ({\itr})','fontname',fnt,'fontsize',fntSz)
        end
        if ctr==1
            l = legend(pl,diag_groups);
            p = get(l,'position');
            set(l,'box','off','fontname',fnt,'fontsize',fntSz-2,...
                'position',[p(1)+0.02,p(2)+0.03,p(3),p(4)])
        elseif ctr==4
            dev1 = dev2;
        end
        
        ctr = ctr+1;
        
    end
end


set(gcf,'color','w')

%% Fig. S4 | Modelling multisensory behavior in teenagers and adults

close all;
clc;

% Generate figure
axes('parent',figure,'fontname',fnt,'fontsize',fntSz)
set(gcf,'position',[20,50,650,600],'paperpositionmode','auto')

% Panel C & D
models = {'Model 1A','Model 1V','Model 2A','Model 2V'};
labels = {'a','b'};
sp = [1,2,5,6;3,4,7,8]+8;
x1 = 0; x2 = 6;
y1 = -0.2; y2 = 0.8;
ytick = -0.2:0.2:0.8;
ctr = 1;
dev1 = 3;
dev2 = 4;
for n = 1:2
    
    for i = 1:numel(models)
        
        rVals = zeros(5,2); pVals = zeros(5,2);
        
        for j = 1:nRatios
            
            % Deine predicted and empirical benefits
            benMod = simCPs(:,i,j,:)-max(CDFs(:,2:3,1,:),[],2);
            benMod = squeeze(trapz(prob,benMod,4));
            
            for m = 1:2
                x = benMod(grp_vec==m & dev_vec==dev1);
                y = benEmp(grp_vec==m & dev_vec==dev1);
                [rVals(j,m),stats] = permtest_corr(x,y,'type',type,'nperm',nperm);
                pVals(j,m) = stats.p;
            end
            
        end
        
        h = subplot(4,4,sp(n,i));
        hold all
        plot([0,6],[0,0],':k')
        pl = plot(1:5,rVals,'-','linewidth',2.5);
        for j = 1:nRatios
            if pVals(j,1)<0.05
                text(j,rVals(j,1)+0.02,'*','fontname',fnt,'fontsize',fntSz+5,'horizontalalign','center')
            end
            if pVals(j,2)<0.05
                text(j,rVals(j,2)+0.02,'*','fontname',fnt,'fontsize',fntSz+5,'horizontalalign','center')
            end
        end
        set(gca,'tickdir','out','xtick',1:5,'xticklabel',{'0','0.25','0.50','0.75','1.00'},...
            'fontname',fnt,'fontsize',fntSz)
        ylim([y1,y2])
        xlim([x1,x2])
        xtickangle(45)
        title(models{i},'fontname',fnt,'fontsize',fntSz)
        if n==2
            set(gca,'ytick',ytick)
            p = get(h,'position');
            set(h,'position',[p(1)+0.04,p(2),p(3),p(4)])
        else
            set(gca,'ytick',ytick)
            p = get(h,'position');
            set(h,'position',[p(1)-0.02,p(2),p(3),p(4)])
        end
        if i==1
            if fig_label
                text(x1-(x2-x1)*0.5,y2,labels{n},'fontname',fnt,'fontsize',fntSz+2,...
                    'fontweight','bold','verticalalign','bottom')
            end
            text(7,y2+0.25,[age_groups{dev1},' yrs'],'horizontalalign','center',...
                'fontweight','normal','fontname',fnt,'fontsize',fntSz)
        end
        if i<3
            set(gca,'xticklabel',[])
            p = get(h,'position');
            set(h,'position',[p(1),p(2)+0.02,p(3),p(4)])
        else
            p = get(h,'position');
            set(h,'position',[p(1),p(2)+0.02,p(3),p(4)])
        end
        if i==2 || i==4
            set(gca,'yticklabel',[])
            p = get(h,'position');
            set(h,'position',[p(1)-0.02,p(2),p(3),p(4)])
        elseif i==3
            p = get(h,'position');
            set(h,'position',[p(1),p(2),p(3),p(4)])
            xlabel('Race\leftarrow{\itp}\rightarrowBias','fontname',fnt,'fontsize',fntSz)
            ylabel('Correlation ({\itr})','fontname',fnt,'fontsize',fntSz)
        end
        if ctr==1
            l = legend(pl,diag_groups);
            p = get(l,'position');
            set(l,'box','off','fontname',fnt,'fontsize',fntSz-2,...
                'position',[p(1)+0.02,p(2)+0.03,p(3),p(4)])
        elseif ctr==4
            dev1 = dev2;
        end
        
        ctr = ctr+1;
        
    end
end


set(gcf,'color','w')

%% Fig. S5 | Group comparisons in MSE

close all;
clc;

% Set parameters
nperm = 1e4;
lnWth = 2;
fig_label = true;
age_match = true;

% Compute RMV moving mean and run perm test
ageRng = 6:26;
window = 7;
x1 = ageRng(1); x2 = ageRng(end);

% Generate figure
axes('parent',figure,'fontname',fnt,'fontsize',fntSz)
set(gcf,'position',[20,50,800,425],'paperpositionmode','auto')

% Preallocate memory
g = zeros(1,numel(ageRng));
ci = zeros(2,numel(ageRng));
bf10 = zeros(1,numel(ageRng));

% Loop through age bins
for i = 1:numel(ageRng)
    
    switch age_match
        case true
            if i>13
                [idx1,idx2] = movmeanmatch(i,grp_vec,sex_vec,age_vec,[],ageRng,window,false);
            else
                [idx1,idx2] = movmeanmatch(i,grp_vec,sex_vec,age_vec,piq_vec,ageRng,window,false);
            end
        case false
            idx1 = find(grp_vec==1 & age_vec>ageRng(i)-window/2 & age_vec<ageRng(i)+window/2);
            idx2 = find(grp_vec==2 & age_vec>ageRng(i)-window/2 & age_vec<ageRng(i)+window/2);
    end
    
    costNT = cost(idx1,2);
    costASD = cost(idx2,2);
            
    % Effect size
    d = mes(costNT,costASD,'hedgesg','isDep',0,'nBoot',nperm);
    g(i) = d.hedgesg;
    ci(:,i) = d.hedgesgCi;
    
    % Bayes factor analysis
    bf10(i) = 1/bf.ttest2(costNT,costASD);
    
end

y1 = 0; y2 = 1.2;
plt = subplot(2,2,1);
hold on
plot(ageRng,g,'color',cols1(5,:),'linewidth',lnWth);
plot(ageRng,ci,':','color',cols1(5,:),'linewidth',lnWth);
title({'Moving MSE (V\rightarrowA)';'Group Comparisons'},...
    'fontname',fnt,'fontsize',fntSz,'fontweight','normal')
ylabel({'Effect size';'(Hedge''s \it{g})'},'fontname',fnt,'fontsize',fntSz)
set(gca,'layer','top','box','off','tickdir','out',...
    'xtick',6:2:26,'ytick',0:0.2:2)
xlim([x1,x2])
ylim([y1,y2])
grid on
if fig_label
    text(x1-(x2-x1)*0.15,y2,'a','fontname',fnt,'fontsize',fntSz+2,...
        'fontweight','bold','verticalalign','bottom')
end

y1 = 0; y2 = 5;
plt = subplot(2,2,3); 
hold on
plot(ageRng,bf10,'color',cols1(6,:),'linewidth',lnWth);
xlabel('Age (yrs)','fontname',fnt,'fontsize',fntSz)
ylabel({'Bayes factor';'(BF_0_1)'},'fontname',fnt,'fontsize',fntSz)
set(gca,'layer','top','box','off','tickdir','out',...
    'xtick',6:2:26,'ytick',0:5)
xlim([x1,x2])
ylim([y1,y2])
grid on
if fig_label
    text(x1-(x2-x1)*0.15,y2,'b','fontname',fnt,'fontsize',fntSz+2,...
        'fontweight','bold','verticalalign','bottom')
end

% Preallocate memory
g = zeros(1,numel(ageRng));
ci = zeros(2,numel(ageRng));
bf10 = zeros(1,numel(ageRng));

% Loop through age bins
for i = 1:numel(ageRng)
    
    switch age_match
        case true
            if i>13
                [idx1,idx2] = movmeanmatch(i,grp_vec,sex_vec,age_vec,[],ageRng,window,false);
            else
                [idx1,idx2] = movmeanmatch(i,grp_vec,sex_vec,age_vec,piq_vec,ageRng,window,false);
            end
        case false
            idx1 = find(grp_vec==1 & age_vec>ageRng(i)-window/2 & age_vec<ageRng(i)+window/2);
            idx2 = find(grp_vec==2 & age_vec>ageRng(i)-window/2 & age_vec<ageRng(i)+window/2);
    end
    
    costNT = cost(idx1,3);
    costASD = cost(idx2,3);
            
    % Effect size
    d = mes(costNT,costASD,'hedgesg','isDep',0,'nBoot',nperm);
    g(i) = d.hedgesg;
    ci(:,i) = d.hedgesgCi;

    % Bayes factor analysis
    bf10(i) = 1/bf.ttest2(costNT,costASD);
    
end

y1 = 0; y2 = 1.2;
plt = subplot(2,2,2);
hold on
plot(ageRng,g,'color',cols1(5,:),'linewidth',lnWth);
plot(ageRng,ci,':','color',cols1(5,:),'linewidth',lnWth);
title({'Moving MSE (A\rightarrowV)';'Group Comparisons'},...
    'fontname',fnt,'fontsize',fntSz,'fontweight','normal')
set(gca,'layer','top','box','off','tickdir','out',...
    'xtick',6:2:26,'ytick',0:0.2:2)
xlim([x1,x2])
ylim([y1,y2])
grid on

y1 = 0; y2 = 5;
plt = subplot(2,2,4); 
hold on
plot(ageRng,bf10,'color',cols1(6,:),'linewidth',lnWth);
xlabel('Age (yrs)','fontname',fnt,'fontsize',fntSz)
set(gca,'layer','top','box','off','tickdir','out',...
    'xtick',6:2:26,'ytick',0:5)
xlim([x1,x2])
ylim([y1,y2])
grid on

set(gcf,'color','w')


%% Fig S6 | Contribution of MSE to RMV

close all;
clc;

xaxis = prob;
nperm = 1e4;

axes('Parent',figure,'FontName',fnt);
set(gcf,'Position',[20,50,600,600],'PaperPositionMode','auto')

% Panel A
x1 = -0.15; x2 = 0.25;
y1 = -0.25; y2 = 0.15;

cond = 7;
task = 3;

% Define race model violation
gain = squeeze(trapz(xaxis,CDFs(:,cond,task,:),4));
swCost = squeeze(CDFs(:,1:3,2,:)-CDFs(:,1:3,3,:));
scArea = squeeze(trapz(xaxis,swCost,3));

ctr = 1;
for i = 1:2
    
    idx = find(grp_vec==i);
    for j = 1:3
        
        sc = scArea(idx,j);
        rmv = gain(idx);
        age = age_vec(idx);
        hr = hitRate(idx);
        
        [r,p] = partialcorr(sc,rmv,age);
        disp([num2str(r^2),',',num2str(p)])
        
        h = subplot(4,3,ctr);
        
        hold all
        scatter(sc,rmv,mkrSz,'filled','markerfacecolor',cols1(i,:))
        H = lsline;
        set(H,'Color','k','LineWidth',lnWth,'LineStyle','-')
        set(gca,'box','off','tickdir','out')
        if p < 0.05
            text(x2,y1,['{\itR}^2 = ',num2str(round(r^2*1e2)/1e2),'*'],...
                'horizontalalign','right','verticalalign','bottom')
        else
            text(x2,y1,['{\itR}^2 = ',num2str(round(r^2*1e2)/1e2)],...
                'horizontalalign','right','verticalalign','bottom')
        end
        if i==1 && j==1
            text(x1-0.2,y2,'a','FontName',fnt,'FontSize',fntSz+2,...
                'fontweight','bold','verticalalignment','bottom')
        end
        if i==1 && j==2
            title({'Multisensory Gain on Switch Trials';switch_conds{j}},...
                'FontName',fnt,'FontSize',fntSz)
        elseif i==1
            title(switch_conds{j},'FontName',fnt,'FontSize',fntSz)
        end
        if i==2
            p = get(h,'position');
            set(h,'position',[p(1),p(2)+0.035,p(3),p(4)])
        else
            set(gca,'xticklabel',[])
        end
        if i==2 && j==1
            xlabel('MSE (AUC)','FontName',fnt,'FontSize',fntSz)
            ylabel('Gain (AUC)','FontName',fnt,'FontSize',fntSz)
        end
        if intersect(j,2:4)
            set(gca,'yticklabel',[])
        end
        if j==1
            text(x1-30,-5,diag_groups{i},'FontName',fnt,'FontSize',fntSz,...
                'horizontalalignment','center','fontweight','normal')
        end
        p = get(h,'position');
        set(h,'position',[p(1),p(2),p(3),p(4)])
        xlim([x1,x2])
        ylim([y1,y2])
        axis square
        ctr = ctr+1;
        
    end
        
end

% Panel B

cond = 7;
task = 2;

gain = squeeze(trapz(xaxis,CDFs(:,cond,task,:),4));

for i = 1:2
    
    %         idx = find(grp_vec==i & dev_vec==4);
    idx = find(grp_vec==i);
    for j = 1:3
        
        sc = scArea(idx,j);
        rmv = gain(idx);
        age = age_vec(idx);
        hr = hitRate(idx);
        
        [r,p] = partialcorr(sc,rmv,age);
        disp([num2str(r^2),',',num2str(p)])
        
        h = subplot(4,3,ctr);
        
        hold all
        scatter(sc,rmv,mkrSz,'filled','markerfacecolor',cols1(i,:))
        H = lsline;
        set(H,'Color','k','LineWidth',lnWth,'LineStyle','-')
        set(gca,'box','off','tickdir','out')
        if p < 0.05
            text(x2,y1,['{\itR}^2 = ',num2str(round(r^2*1e2)/1e2),'*'],...
                'horizontalalign','right','verticalalign','bottom')
        else
            text(x2,y1,['{\itR}^2 = ',num2str(round(r^2*1e2)/1e2)],...
                'horizontalalign','right','verticalalign','bottom')
        end
        if i==1 && j==1
            text(x1-0.2,y2,'b','FontName',fnt,'FontSize',fntSz+2,...
                'fontweight','bold','verticalalignment','bottom')
        end
        if i==1 && j==2
            title({'Multisensory Gain on Repeat Trials';switch_conds{j}},...
                'FontName',fnt,'FontSize',fntSz)
        elseif i==1
            title(switch_conds{j},'FontName',fnt,'FontSize',fntSz)
        end
        if i==2
            p = get(h,'position');
            set(h,'position',[p(1),p(2)+0.035,p(3),p(4)])
        else
            set(gca,'xticklabel',[])
        end
        if i==2 && j==1
            xlabel('MSE (AUC)','FontName',fnt,'FontSize',fntSz)
            ylabel('Gain (AUC)','FontName',fnt,'FontSize',fntSz)
        end
        if intersect(j,2:4)
            set(gca,'yticklabel',[])
        end
        if j==1
            text(x1-30,-5,diag_groups{i},'FontName',fnt,...
                'FontSize',fntSz,'horizontalalignment','center','fontweight','normal')
        end
        p = get(h,'position');
        set(h,'position',[p(1),p(2)-0.05,p(3),p(4)])
        xlim([x1,x2])
        ylim([y1,y2])
        axis square
        ctr = ctr+1;
        
    end
        
end

set(gcf,'color','w')


%% Fig S7a | Group differences & recovery

close all;
clc;

% Set parameters
age_match = true;
nperm = 1e4;
lnWth = 2;
fig_label = true;
cond = 7;
task = 2;
xaxis = prob;
fdr_method = 'bh'; % 'bh' or 'bky' or 'none'

gainRe = squeeze(trapz(xaxis,CDFs(:,cond,task,:),4));

% Compute RMV moving mean and run perm test
ageRng = 6:26;
window = 7;
x1 = ageRng(1); x2 = ageRng(end);
switch test
    case 'ver'
        y1 = -0.06; y2 = 0.03;
    case 'hor'
        y1 = -60; y2 = 20;
end

switch test
    case 'ver'
        units = ' (AUC)';
        xlab1 = 'Quantile';
        xlab2 = 'Quantile';
    case 'hor'
        units = ' (ms)';
        xlab1 = 'RT (ms)';
        xlab2 = 'Percentile';
end


% Preallocate memory
rmvAvg = zeros(numel(ageRng),2);
errNT = zeros(numel(ageRng),2);
errASD = zeros(numel(ageRng),2);
p = zeros(1,numel(ageRng));
h = zeros(1,numel(ageRng));
g = zeros(1,numel(ageRng));
ci = zeros(2,numel(ageRng));
bf10 = zeros(1,numel(ageRng));

% Loop through age bins
for i = 1:numel(ageRng)

    switch age_match
        case true
            if i>13
                [idx1,idx2] = movmeanmatch(i,grp_vec,sex_vec,age_vec,[],ageRng,window,false);
            else
                [idx1,idx2] = movmeanmatch(i,grp_vec,sex_vec,age_vec,piq_vec,ageRng,window,false);
            end
        case false
            idx1 = find(grp_vec==1 & age_vec>ageRng(i)-window/2 & age_vec<ageRng(i)+window/2);
            idx2 = find(grp_vec==2 & age_vec>ageRng(i)-window/2 & age_vec<ageRng(i)+window/2);
    end
    
    gainNT = gainRe(idx1);
    gainASD = gainRe(idx2);
    
    % Compute RMV and CIs
    rmvAvg(i,1) = squeeze(mean(gainNT));
    rmvci = bootci(nperm,@mean,gainNT);
    errNT(i,1) = rmvAvg(i,1)-rmvci(1);
    errNT(i,2) = rmvci(2)-rmvAvg(i,1);
    rmvAvg(i,2) = squeeze(mean(gainASD));
    rmvci = bootci(nperm,@mean,gainASD);
    errASD(i,1) = rmvAvg(i,2)-rmvci(1);
    errASD(i,2) = rmvci(2)-rmvAvg(i,2);
    
    % Compare variance before running permutation tests
    [~,~,orig] = permtest_f2(gainNT,gainASD,nperm);
    if mean(orig.h)>0.5
        disp(['Unequal variance: ',num2str(ageRng(i)),' years'])
    end
    
    % Run permutation tests at each quantile
    [~,~,orig] = permtest_t2(gainNT,gainASD,'nperm',nperm);
    bf10(i) = 1/bf.ttest2(gainNT,gainASD);
    d = mes(gainNT,gainASD,'hedgesg','isDep',0,'nBoot',nperm);
    g(i) = d.hedgesg;
    ci(:,i) = d.hedgesgCi;
    
    p(i) = orig.p;
    h(i) = orig.h;
    
end

% Apply FDR correction
switch fdr_method
    case 'bky'
        hx = fdr_bky(p,alpha);
    case 'bh'
        hx = fdr_bh(p,alpha,'pdep');
    case 'none'
        hx = h;
end

% Plot data

% Generate figure
axes('parent',figure,'fontname',fnt,'fontsize',fntSz)
set(gcf,'position',[20,50,500,450],'paperpositionmode','auto')

plt = subplot(2,8,2:7);
hold all
plot(ageRng,rmvAvg(:,1)+1e3,'color',cols1(1,:),'linewidth',lnWth)
plot(ageRng,rmvAvg(:,2)+1e3,'color',cols1(2,:),'linewidth',lnWth)
xaxis = linspace(ageRng(1),ageRng(end),length(ageRng));
for k = 2
    if k == 1
        hCorr = h;
        shade = 0.925;
    else
        hCorr = hx;
        shade = 0.85;
    end
    try
        sigX = xaxis(hCorr==1);
        for l = 1:length(sigX)
            X = [sigX(l)-0.5,sigX(l)-0.5,sigX(l)+0.5,sigX(l)+0.5];
            Y = [y1,y2,y2,y1];
            C = ones(1,3)*shade;
            fill(X,Y,C,'edgecolor','none')
        end
    catch
    end
end
plot([x1,x2],[0,0],':k')
[l,p] = boundedline(ageRng,rmvAvg(:,1),errNT,ageRng,rmvAvg(:,2),errASD,...
    'alpha','cmap',cols1(1:2,:));
outlinebounds(l,p);
plot(ageRng,rmvAvg(:,1),'color',cols1(1,:),'linewidth',lnWth)
plot(ageRng,rmvAvg(:,2),'color',cols1(2,:),'linewidth',lnWth)
ylabel({'Moving mean';['gain {\itk}=7',units]},'fontname',fnt,'fontsize',fntSz)
l = legend(diag_groups);
set(l,'box','off','location','southeast','fontname',fnt,'fontsize',fntSz-2)
switch test
    case 'ver'
        set(gca,'layer','top','box','off','tickdir','out',...
            'xtick',6:2:26,'ytick',-0.06:0.02:0.04)
    case 'hor'
        set(gca,'layer','top','box','off','tickdir','out',...
            'xtick',6:2:26,'ytick',-60:20:20)
end
xlim([x1,x2])
ylim([y1,y2])
if fig_label
    text(x1-(x2-x1)*0.15,y2,'a','fontname',fnt,'fontsize',fntSz+2,...
        'fontweight','bold','verticalalign','bottom')
end

y1 = 0; y2 = 1.2;
yy1 = 0; yy2 = 6;
plt = subplot(2,8,10:15); 
hold on
[ax,h1,h2] = plotyy(ageRng,g,ageRng,bf10,'plot');
set(h1,'color',cols1(5,:),'linewidth',lnWth)
set(h2,'color',cols1(6,:),'linewidth',lnWth)
l = legend('H''s \it{g}','BF_0_1');
set(l,'box','off','location','northwest','fontname',fnt,'fontsize',fntSz-2)
pos = get(l,'pos');
set(l,'position',[pos(1)-0.025,pos(2)+0.02,pos(3),pos(4)+0.02])
% hold(ax(1));
plot(ax(1),ageRng,ci,':','color',cols1(5,:),'linewidth',lnWth);
xlabel('Age (yrs)','fontname',fnt,'fontsize',fntSz)
ylabel(ax(1),{'Effect size';'(Hedge''s \it{g})'},'fontname',fnt,'fontsize',fntSz)
ylabel(ax(2),{'Bayes factor';'(BF_0_1)'},'fontname',fnt,'fontsize',fntSz)
set(gca,'layer','bottom','box','off','tickdir','out',...
    'xtick',6:2:26,'ytick',0:0.2:1.2)
set(ax(2),'ytick',0:6)
set(ax,'tickdir','out')
xlim(ax(1),[x1,x2])
xlim(ax(2),[x1,x2])
ylim(ax(1),[y1,y2])
ylim(ax(2),[yy1,yy2])
grid on
if fig_label
    text(x1-(x2-x1)*0.15,y2,'b','fontname',fnt,'fontsize',fntSz+2,...
        'fontweight','bold','verticalalign','bottom')
end
set(ax(1),'ycolor','k')
set(ax(2),'ycolor','k')

set(gcf,'color','w')

%% Fig S7b | Group differences & recovery

close all;
clc;

% Set parameters
age_match = true;
nperm = 1e4;
lnWth = 2;
fig_label = true;
cond = 7;
task = 2;
xaxis = prob;
fdr_method = 'bh'; % 'bh' or 'bky' or 'none'

gainRe = squeeze(trapz(xaxis,CDFs(:,cond,task,:),4));

% Compute RMV moving mean and run perm test
ageRng = 6:26;
window = 7;
x1 = ageRng(1); x2 = ageRng(end);
switch test
    case 'ver'
        y1 = -0.06; y2 = 0.03;
    case 'hor'
        y1 = -60; y2 = 20;
end

switch test
    case 'ver'
        units = ' (AUC)';
        xlab1 = 'Quantile';
        xlab2 = 'Quantile';
    case 'hor'
        units = ' (ms)';
        xlab1 = 'RT (ms)';
        xlab2 = 'Percentile';
end


% Preallocate memory
rmvAvg = zeros(numel(ageRng),2);
errNT = zeros(numel(ageRng),2);
errASD = zeros(numel(ageRng),2);
p = zeros(1,numel(ageRng));
h = zeros(1,numel(ageRng));
g = zeros(1,numel(ageRng));
ci = zeros(2,numel(ageRng));
bf10 = zeros(1,numel(ageRng));

% Loop through age bins
for i = 1:numel(ageRng)

    switch age_match
        case true
            if i>13
                [idx1,idx2] = movmeanmatch(i,grp_vec,sex_vec,age_vec,[],ageRng,window,false);
            else
                [idx1,idx2] = movmeanmatch(i,grp_vec,sex_vec,age_vec,piq_vec,ageRng,window,false);
            end
        case false
            idx1 = find(grp_vec==1 & age_vec>ageRng(i)-window/2 & age_vec<ageRng(i)+window/2);
            idx2 = find(grp_vec==2 & age_vec>ageRng(i)-window/2 & age_vec<ageRng(i)+window/2);
    end
    
    gainNT = gainRe(idx1);
    gainASD = gainRe(idx2);
    
    % Compute RMV and CIs
    rmvAvg(i,1) = squeeze(mean(gainNT));
    rmvci = bootci(nperm,@mean,gainNT);
    errNT(i,1) = rmvAvg(i,1)-rmvci(1);
    errNT(i,2) = rmvci(2)-rmvAvg(i,1);
    rmvAvg(i,2) = squeeze(mean(gainASD));
    rmvci = bootci(nperm,@mean,gainASD);
    errASD(i,1) = rmvAvg(i,2)-rmvci(1);
    errASD(i,2) = rmvci(2)-rmvAvg(i,2);
    
    % Compare variance before running permutation tests
    [~,~,orig] = permtest_f2(gainNT,gainASD,nperm);
    if mean(orig.h)>0.5
        disp(['Unequal variance: ',num2str(ageRng(i)),' years'])
    end
    
    % Run permutation tests at each quantile
    [~,~,orig] = permtest_t2(gainNT,gainASD,'nperm',nperm);
    bf10(i) = 1/bf.ttest2(gainNT,gainASD);
    d = mes(gainNT,gainASD,'hedgesg','isDep',0,'nBoot',nperm);
    g(i) = d.hedgesg;
    ci(:,i) = d.hedgesgCi;
    
    p(i) = orig.p;
    h(i) = orig.h;
    
end

% Apply FDR correction
switch fdr_method
    case 'bky'
        hx = fdr_bky(p,alpha);
    case 'bh'
        hx = fdr_bh(p,alpha,'pdep');
    case 'none'
        hx = h;
end

% Plot data

% Generate figure
axes('parent',figure,'fontname',fnt,'fontsize',fntSz)
set(gcf,'position',[20,50,500,650],'paperpositionmode','auto')

plt = subplot(2,8,2:7);
hold all
plot(ageRng,rmvAvg(:,1)+1e3,'color',cols1(1,:),'linewidth',lnWth)
plot(ageRng,rmvAvg(:,2)+1e3,'color',cols1(2,:),'linewidth',lnWth)
xaxis = linspace(ageRng(1),ageRng(end),length(ageRng));
for k = 2
    if k == 1
        hCorr = h;
        shade = 0.925;
    else
        hCorr = hx;
        shade = 0.85;
    end
    try
        sigX = xaxis(hCorr==1);
        for l = 1:length(sigX)
            X = [sigX(l)-0.5,sigX(l)-0.5,sigX(l)+0.5,sigX(l)+0.5];
            Y = [y1,y2,y2,y1];
            C = ones(1,3)*shade;
            fill(X,Y,C,'edgecolor','none')
        end
    catch
    end
end
plot([x1,x2],[0,0],':k')
[l,p] = boundedline(ageRng,rmvAvg(:,1),errNT,ageRng,rmvAvg(:,2),errASD,...
    'alpha','cmap',cols1(1:2,:));
outlinebounds(l,p);
plot(ageRng,rmvAvg(:,1),'color',cols1(1,:),'linewidth',lnWth)
plot(ageRng,rmvAvg(:,2),'color',cols1(2,:),'linewidth',lnWth)
ylabel({'Moving mean';['gain {\itk}=7',units]},'fontname',fnt,'fontsize',fntSz)
l = legend(diag_groups);
set(l,'box','off','location','southeast','fontname',fnt,'fontsize',fntSz-2)
switch test
    case 'ver'
        set(gca,'layer','top','box','off','tickdir','out',...
            'xtick',6:2:26,'ytick',-0.06:0.02:0.04)
    case 'hor'
        set(gca,'layer','top','box','off','tickdir','out',...
            'xtick',6:2:26,'ytick',-60:20:20)
end
xlim([x1,x2])
ylim([y1,y2])
if fig_label
    text(x1-(x2-x1)*0.15,y2,'a','fontname',fnt,'fontsize',fntSz+2,...
        'fontweight','bold','verticalalign','bottom')
end

y1 = 0; y2 = 1;
plt = subplot(2,8,10:15); 
hold on
plot(ageRng,g,'color',cols1(5,:),'linewidth',lnWth);
plot(ageRng,bf10,'color',cols1(6,:),'linewidth',lnWth);
l = legend('H''s \it{g}');
l = legend('BF_0_1');
set(l,'box','off','location','northwest','fontname',fnt,'fontsize',fntSz-2)
pos = get(l,'pos');
set(l,'position',[pos(1)-0.025,pos(2)+0.02,pos(3),pos(4)+0.02])
% hold(ax(1));
plot(ax(1),ageRng,ci,':','color',cols1(5,:),'linewidth',lnWth);
xlabel('Age (yrs)','fontname',fnt,'fontsize',fntSz)
ylabel(ax(1),{'Effect size';'(Hedge''s \it{g})'},'fontname',fnt,'fontsize',fntSz)
ylabel(ax(2),{'Bayes factor';'(BF_0_1)'},'fontname',fnt,'fontsize',fntSz)
set(gca,'layer','bottom','box','off','tickdir','out',...
    'xtick',6:2:26)
set(ax,'tickdir','out')
xlim(ax(1),[x1,x2])
xlim(ax(2),[x1,x2])
ylim(ax(1),[y1,y2])
ylim(ax(2),[yy1,yy2])
grid on
if fig_label
    text(x1-(x2-x1)*0.15,y2,'b','fontname',fnt,'fontsize',fntSz+2,...
        'fontweight','bold','verticalalign','bottom')
end

set(gcf,'color','w')

%% Table S1. Race model test with Miller's bound

clc;

nperm = 1e4;

quant = 2:8;
rmv = squeeze(CDFs(:,1,1,:)-CDFs(:,8,1,:));

dvals  = zeros(2,4,7);
cis  = zeros(2,4,2,7);
pvals  = zeros(2,4,7);

for i = 1:2
    for j = 1:4
        idx = grp_vec==i & dev_vec==j;
        [tstat,corx] = st_tmaxperm(rmv(idx,quant),[],'nperm',nperm,'tail','right');
        d = mes(squeeze(CDFs(idx,1,1,quant)),squeeze(CDFs(idx,8,1,quant)),'hedgesg','isDep',1,'nBoot',nperm);
        dvals(i,j,:) = d.hedgesg;
        cis(i,j,:,:) = d.hedgesgCi;
        pvals(i,j,:) = corx.p;
    end
end

for k = 1:length(quant)
    for i = 1:2
        for j = 1:4
            if pvals(i,j,k)<0.05
                fprintf('%1.2f[%1.1f,%1.1f]*\t',dvals(i,j,k),cis(i,j,1,k),cis(i,j,2,k))
            else
                fprintf('%1.2f[%1.1f,%1.1f]\t',dvals(i,j,k),cis(i,j,1,k),cis(i,j,2,k))
            end
%             fprintf('p = %1.3f, d = %1.2f\t',pvals(i,j,k),dvals(i,j,k))
        end
    end
    fprintf('\n')
end

%% Table S2. Race model test using repeat trials only

clc;

nperm = 1e4;

quant = 2:8;
rmv = squeeze(CDFs(:,1,2,:)-CDFs(:,6,2,:));

dvals  = zeros(2,4,7);
cis  = zeros(2,4,2,7);
pvals  = zeros(2,4,7);

for i = 1:2
    for j = 1:4
        idx = grp_vec==i & dev_vec==j;
        [tstat,corx] = permtest_t(rmv(idx,quant),[],'nperm',nperm,'tail','right');
        d = mes(squeeze(CDFs(idx,1,2,quant)),squeeze(CDFs(idx,6,2,quant)),'hedgesg','isDep',1,'nBoot',nperm);
        dvals(i,j,:) = d.hedgesg;
        cis(i,j,:,:) = d.hedgesgCi;
        pvals(i,j,:) = corx.p;
    end
end

for k = 1:length(quant)
    for i = 1:2
        for j = 1:4
            if pvals(i,j,k)<0.05
                fprintf('%1.2f[%1.1f,%1.1f]*\t',dvals(i,j,k),cis(i,j,1,k),cis(i,j,2,k))
            else
                fprintf('%1.2f[%1.1f,%1.1f]\t',dvals(i,j,k),cis(i,j,1,k),cis(i,j,2,k))
            end
%             fprintf('p = %1.3f, d = %1.2f\t',pvals(i,j,k),dvals(i,j,k))
        end
    end
    fprintf('\n')
end

%% Table S3. Race model test with Miller's bound using repeat trials only

clc;

nperm = 1e4;

quant = 2:8;
rmv = squeeze(CDFs(:,1,2,:)-CDFs(:,8,2,:));

dvals  = zeros(2,4,7);
cis  = zeros(2,4,2,7);
pvals  = zeros(2,4,7);

for i = 1:2
    for j = 1:4
        idx = grp_vec==i & dev_vec==j;
        [tstat,corx] = st_tmaxperm(rmv(idx,quant),[],'nperm',nperm,'tail','right');
        d = mes(squeeze(CDFs(idx,1,2,quant)),squeeze(CDFs(idx,8,2,quant)),'hedgesg','isDep',1,'nBoot',nperm);
        dvals(i,j,:) = d.hedgesg;
        cis(i,j,:,:) = d.hedgesgCi;
        pvals(i,j,:) = corx.p;
    end
end

for k = 1:length(quant)
    for i = 1:2
        for j = 1:4
            if pvals(i,j,k)<0.05
                fprintf('%1.2f[%1.1f,%1.1f]*\t',dvals(i,j,k),cis(i,j,1,k),cis(i,j,2,k))
            else
                fprintf('%1.2f[%1.1f,%1.1f]\t',dvals(i,j,k),cis(i,j,1,k),cis(i,j,2,k))
            end
%             fprintf('p = %1.3f, d = %1.2f\t',pvals(i,j,k),dvals(i,j,k))
        end
    end
    fprintf('\n')
end

%% Fig 2b: Group mean RTs for each condition

close all
clc

quant = 3;

rmvMeth = 'pos'; % 'quant' or 'all' or 'pos'

if strcmp(rmvMeth,'quant')
    rmvArea = squeeze(CDFs(:,5,1,quant));
elseif strcmp(rmvMeth,'all')
    rmvArea = squeeze(trapz(CDFs(:,5,1,:),4));
elseif strcmp(rmvMeth,'pos')
    rmvArea = zeros(nSubjs,1);
    for i = 1:nSubjs
        % Only consider values above zero
        rmvIdx = find(squeeze(CDFs(i,5,1,:)>0));
        % Ignore RMV starting after 6th quantile
        if ~isempty(rmvIdx)
            if rmvIdx(1)>6
                rmvIdx = [];
            end
        end
        % Ignore RMV after break in positive values
        rmvGap = find(diff(rmvIdx)>1,1);
        if ~isempty(rmvGap)
            rmvIdx = rmvIdx(1:rmvGap);
        end
        % Compute area using trapezoidal method
        if ~isempty(rmvIdx)
            if rmvIdx(1)>1
                x1 = 1-abs(CDFs(i,5,1,rmvIdx(1)-1))/(CDFs(i,5,1,rmvIdx(1))+abs(CDFs(i,5,1,rmvIdx(1)-1)));
            end
            if rmvIdx(end)<length(prob)
                x2 = 1-abs(CDFs(i,5,1,rmvIdx(end)+1))/(CDFs(i,5,1,rmvIdx(end))+abs(CDFs(i,5,1,rmvIdx(end)+1)));
            end
            rmvArea(i) = trapz([1-x1;(1:length(rmvIdx))';length(rmvIdx)+x2],[0;squeeze(CDFs(i,5,1,rmvIdx));0]);
        end
        clear rmvIdx rmvGap
    end
end

% Preallocate memory
idx = cell(2,4);
RTrmv = zeros(2,4,3);
RTnormv = zeros(2,4,3);
err = zeros(2,2,4,3);
ci = zeros(2,2,4,3);
nperm = 1e4;

% Compute group RT central tendancy and 95% CIs
ctMeth = 'mean'; % 'mean', 'med', 'mode'
for i = 1:2
    for j = 1:4
        idx{i,j} = find(grp_vec==i & dev_vec==j & rmvArea>0);
        if strcmp(ctMeth,'mean')
            RTrmv(i,j,:) = squeeze(mean(RTct(idx{i,j},:,1),1));
            %             ci(:,i,j,:) = bootci(nperm,@mean,squeeze(RTct(idx{i,j},:,1)));
        elseif strcmp(ctMeth,'med')
            RTrmv(i,j,:) = squeeze(median(RTct(idx{i,j},:,1),1));
            %             ci(:,i,j,:) = bootci(nperm,@median,squeeze(RTct(idx{i,j},:,1)));
        elseif strcmp(ctMeth,'mode')
            RTrmv(i,j,:) = squeeze(mode(round(RTct(idx{i,j},:,1)),1));
            %             ci(:,i,j,:) = bootci(nperm,@mode,squeeze(round(RTct(idx{i,j},:,1))));
        end
        idx{i,j} = find(grp_vec==i & dev_vec==j & rmvArea==0);
        if strcmp(ctMeth,'mean')
            RTnormv(i,j,:) = squeeze(mean(RTct(idx{i,j},:,1),1));
            %             ci(:,i,j,:) = bootci(nperm,@mean,squeeze(RTct(idx{i,j},:,1)));
        elseif strcmp(ctMeth,'med')
            RTnormv(i,j,:) = squeeze(median(RTct(idx{i,j},:,1),1));
            %             ci(:,i,j,:) = bootci(nperm,@median,squeeze(RTct(idx{i,j},:,1)));
        elseif strcmp(ctMeth,'mode')
            RTnormv(i,j,:) = squeeze(mode(round(RTct(idx{i,j},:,1)),1));
            %             ci(:,i,j,:) = bootci(nperm,@mode,squeeze(round(RTct(idx{i,j},:,1))));
        end
    end
end

% Compute error
% err(1,:,:,:) = RTavg-squeeze(ci(1,:,:,:));
% err(2,:,:,:) = squeeze(ci(2,:,:,:))-RTavg;

% Generate figure
axes('parent',figure,'fontname',fnt,'fontsize',fntSz)
set(gcf,'position',[20,150,500,400],'paperpositionmode','auto')
for i = 1:2
    subplot(2,2,i)
    hold on
    for j = 1:3
        plot(1:4,squeeze(RTrmv(i,:,j)),'-','linewidth',lnWth-0.5,'color',cols1(j,:))
    end
    for j = 1:3
        plot(1:4,squeeze(RTnormv(i,:,j)),'--','linewidth',lnWth-0.5,'color',cols1(j,:))
    end
    set(gca,'box','off','xtick',1:4,'xticklabel',age_groups,'tickdir','out','fontname',fnt,'fontsize',fntSz)
    title(diag_groups{i},'fontname',fnt,'fontsize',fntSz)
    %     xlabel('Age (yrs)','fontname',fnt,'fontsize',fntSz)
    if i==1
        ylabel('RT (ms)','fontname',fnt,'fontsize',fntSz)
        l = legend(conditions);
        %         pos = get(l,'position');
        %         posnew = pos; posnew(1) = posnew(1)+0.03; posnew(2) = posnew(2)-0.08;
        set(l,'box','off','location','northeast','fontname',fnt,'fontsize',fntSz-2,'position',posnew)
    end
    xlim([0,5])
    %     ylim([200,600])
    xtickangle(45)
    axis square
end
for i = 1:2
    subplot(2,2,i+2)
    hold on
    plot(1:4,squeeze(RTnormv(i,:,:)-RTrmv(i,:,:)),'linewidth',lnWth-0.5)
    set(gca,'box','off','xtick',1:4,'xticklabel',age_groups,'tickdir','out','fontname',fnt,'fontsize',fntSz)
    %     title(diag_groups{i},'fontname',fnt,'fontsize',fntSz)
    xlabel('Age (yrs)','fontname',fnt,'fontsize',fntSz)
    if i==1
        ylabel('RT (ms)','fontname',fnt,'fontsize',fntSz)
        l = legend(conditions);
        %         pos = get(l,'position');
        %         posnew = pos; posnew(1) = posnew(1)+0.03; posnew(2) = posnew(2)-0.08;
        set(l,'box','off','location','northeast','fontname',fnt,'fontsize',fntSz-2,'position',posnew)
    end
    xlim([0,5])
    ylim([-50,150])
    xtickangle(45)
    axis square
end

set(gcf,'color','w')

clear idx RTavg err

%% Fig 6: Sex Differences

% close all;
clc;

% Define parameters
age_match = true;
label = {'A','B'};
alpha = 0.05;
nperm = 1e4;
quant = 2:length(prob)-1;
cond = 7;
task = 6:8;
xaxis = prob/100;

% Define race model violation 
gain = squeeze(trapz(xaxis,mean(CDFs(:,cond,task,:),3),4));

% Initialize counter
ctr = 1;

% Generate figure
axes('Parent',figure,'FontName',fnt);
set(gcf,'Position',[20,50,700,550],'PaperPositionMode','auto')
for i = 1:2
    
    x1 = 0; x2 = 1;
    y1 = -0.14; y2 = 0.14;
    
    % Panels A & B
    for j = 1:4
        
        switch age_match
            case true
                if j==4
                    [idx1,idx2] = sexmatch(i,j,grp_vec,dev_vec,age_vec,sex_vec,[],false);
                else
                    [idx1,idx2] = sexmatch(i,j,grp_vec,dev_vec,age_vec,sex_vec,piq_vec,false);
                end
            case false
                idx1 = find(sex_vec==1 & grp_vec==i & dev_vec==j);
                idx2 = find(sex_vec==2 & grp_vec==i & dev_vec==j);
        end
        
%         if age_match == true
%             
%             % Get ages of all subjects in each group
%             sz1 = sum(sex_vec==1 & grp_vec==i & dev_vec==j);
%             sz2 = sum(sex_vec==2 & grp_vec==i & dev_vec==j);
%             if sz1>sz2
%                 grp1 = 1; grp2 = 2;
%             else
%                 grp1 = 2; grp2 = 1;
%             end
%             idx1_all = find(sex_vec==grp1 & grp_vec==i & dev_vec==j);
%             idx2_all = find(sex_vec==grp2 & grp_vec==i & dev_vec==j);
%             age1 = age_vec(idx1_all);
%             age2 = age_vec(idx2_all);
%             
%             % Sex- and age-match controls and patients
%             idx = zeros(length(age2),1);
%             for k = 1:length(age2)
%                 [~,idx(k)] = min(abs(age1-age2(k)));
%                 age1(idx(k)) = NaN;
%             end
%             if sz1>sz2
%                 idx1 = idx1_all(idx);
%                 idx2 = idx2_all;
%             else
%                 idx1 = idx2_all;
%                 idx2 = idx1_all(idx);
%             end
%             
%         elseif age_match == false
%             
%             idx1 = find(sex_vec==1 & grp_vec==i & dev_vec==j);
%             idx2 = find(sex_vec==2 & grp_vec==i & dev_vec==j);
%             
%         end
        
        % Compute RMV and 95% CIs
        rmvM = squeeze(mean(CDFs(idx1,cond,task,:),3));
        ci = bootci(nperm,@mean,rmvM);
        errM = zeros(length(prob),2);
        errM(:,1) = mean(rmvM)-ci(1,:);
        errM(:,2) = ci(2,:)-mean(rmvM);
        
        rmvF = squeeze(mean(CDFs(idx2,cond,task,:),3));
        ci = bootci(nperm,@mean,rmvF);
        errF = zeros(length(prob),2);
        errF(:,1) = mean(rmvF)-ci(1,:);
        errF(:,2) = ci(2,:)-mean(rmvF);
        
        [~,corx,orig] = permtest_t2(rmvM(:,quant),rmvF(:,quant),'nperm',nperm);
        
        subplot(3,4,ctr)
        hold all
        xaxis = prob;
        plot(xaxis,mean(rmvF)+1e3,'color',cols1(5,:),'LineWidth',lnWth)
        plot(xaxis,mean(rmvM)+1e3,'color',cols1(4,:),'LineWidth',lnWth)
        for k = 2
            if k == 1
                hCorr = [0,orig.h,0];
                shade = 0.9;
            else
                hCorr = [0,corx.h,0];
                shade = 0.8;
            end
            try
                sigX = xaxis(hCorr==1);
                for l = 1:length(sigX)
                    X = [sigX(l)-0.025,sigX(l)-0.025,sigX(l)+0.025,sigX(l)+0.025];
                    Y = [y1,y2,y2,y1];
                    C = ones(1,3)*shade;
                    fill(X,Y,C,'edgecolor','none')
                end
            catch
            end
        end
        [hl,hp] = boundedline(xaxis,mean(rmvF),errF,xaxis,mean(rmvM),...
            errM,'alpha','cmap',cols1([5,4],:)); outlinebounds(hl,hp);
        plot(xaxis,mean(rmvF),'color',cols1(5,:),'LineWidth',lnWth)
        plot(xaxis,mean(rmvM),'color',cols1(4,:),'LineWidth',lnWth)
        set(gca,'layer','top','box','off','xaxislocation','origin')
        set(gca,'tickdir','both','xticklabel',[],'ytick',-0.1:0.1:0.1)
        xlim([x1,x2])
        ylim([y1,y2])
        axis square
        if i==1
            title([age_groups{j},' yrs'],'FontName',fnt,'FontSize',fntSz)
        end
        if j==1
            ylabel('Residuals','FontName',fnt,'FontSize',fntSz)
            text(x1-0.45,y2+0.03,label{i},'FontName',fnt,'FontSize',fntSz+8)
        elseif intersect(j,2:4)
            set(gca,'yticklabel',[])
        end
        if i==1 && j==1
            l = legend('Fem.','Male');
            set(l,'box','off','location','northeast','FontName',fnt,'FontSize',fntSz-2)
        elseif j==3
            text(x1-0.15,y2+0.05,diag_groups{i},'FontName',fnt,...
                'FontSize',fntSz,'Fontweight','bold','horizontalalignment','center')
        end
        
        % Plot RMV-area
        p = get(gca,'Position');
        h = axes('Parent',gcf,'Position',[p(1)+.1 p(2)+.008 p(3)-.1 p(4)-.14]);
        hold on
        [~,~,orig] = permtest_t2(gain(idx1),gain(idx2),'nperm',nperm);
        bar(h,1,mean(gain(idx2)),0.5,'facecolor',cols1(5,:));
        bar(h,2,mean(gain(idx1)),0.5,'facecolor',cols1(4,:));
        if orig.h==1
            text(1.5,-0.0003,'*','FontSize',fntSz+2,'horizontalalignment','center','FontName',fnt)
        end
        set(h,'box','on','XTick',[],'YTick',[]);
        yyaxis left
        ylabel('Area','FontName',fnt,'FontSize',fntSz-2,'verticalalignment','middle')
        xlim([0.5,2.5])
        ylim([-0.0007,0.0003])
        err1 = std(gain(idx2))/sqrt(length(idx2));
        err2 = std(gain(idx1))/sqrt(length(idx1));
        errorbar(1:2,[mean(gain(idx2)),mean(gain(idx1))],[err1,err2],[err1,err2],'linestyle','none','color','k')
        yyaxis right
        set(h,'ytick',-0.05:0.05:0,'yticklabel',{'-.05','0'},{'YColor'},{'k'});
        ylim([-0.0007,0.0003]) % need this a second time
        
        h = 0;
        
        % Run permutation test with tmax correction
        if h==0
            [tstat,~,orig,stats] = permtest_t2(gain(idx2),gain(idx1),'nperm',nperm,[],[],'equal');
            d = mes(gain(idx2),gain(idx1),'hedgesg','nBoot',nperm);
            fprintf('RMIm = %0.2f ± %0.2f, ',mean(gain(idx1)),std(gain(idx1)))
            fprintf('RMIf = %0.2f ± %0.2f, ',mean(gain(idx2)),std(gain(idx2)))
            fprintf('t(%d) = %0.2f, ',stats.df,tstat)
            fprintf('p = %0.3f, ',orig.p)
            fprintf('d = %0.2f, ',d.hedgesg)
            fprintf('95CI [%0.2f, %0.2f]\n',d.hedgesgCi)
            d = d.hedgesg;
        elseif h==1
            [tstat,~,orig,stats] = permtest_t2(gain(idx2),gain(idx1),'nperm',nperm,[],[],'equal');
            d = mes(gain(idx2),gain(idx1),'glassdelta','nBoot',nperm);
            fprintf('RMIm = %0.2f ± %0.2f, ',mean(gain(idx1)),std(gain(idx1)))
            fprintf('RMIf = %0.2f ± %0.2f, ',mean(gain(idx2)),std(gain(idx2)))
            fprintf('t(%d) = %0.2f, ',stats.df,tstat)
            fprintf('p = %0.3f, ',orig.p)
            fprintf('d = %0.2f, ',d.glassdelta)
            fprintf('95CI [%0.2f, %0.2f]\n',d.glassdeltaCi)
            d = d.glassdelta;
        end
        
        ctr = ctr+1;
        
    end
    
    % Panel C
    window = 7;
    ageRng = 6:30;
    x1 = ageRng(1); x2 = ageRng(end);
    y1 = -0.05; y2 = 0.05;
    rmvAvg = zeros(length(ageRng),2);
    err = zeros(length(ageRng),2,2);
    p = zeros(length(ageRng),1);
    h = zeros(length(ageRng),1);
    
    for j = 1:length(ageRng)
        
        if age_match == true
            
            % Get ages of all subjects in each group
            sz1 = sum(sex_vec==1 & grp_vec==i & age_vec>ageRng(j)-window/2 & age_vec<ageRng(j)+window/2);
            sz2 = sum(sex_vec==2 & grp_vec==i & age_vec>ageRng(j)-window/2 & age_vec<ageRng(j)+window/2);
            if sz1>sz2
                grp1 = 1; grp2 = 2;
            else
                grp1 = 2; grp2 = 1;
            end
            idx1_all = find(sex_vec==grp1 & grp_vec==i & age_vec>ageRng(j)-window/2 & age_vec<ageRng(j)+window/2);
            idx2_all = find(sex_vec==grp2 & grp_vec==i & age_vec>ageRng(j)-window/2 & age_vec<ageRng(j)+window/2);
            age1 = age_vec(idx1_all);
            age2 = age_vec(idx2_all);
            
            % Sex- and age-match controls and patients
            idx = zeros(length(age2),1);
            for k = 1:length(age2)
                [~,idx(k)] = min(abs(age1-age2(k)));
                age1(idx(k)) = NaN;
            end
            if sz1>sz2
                idx1 = idx1_all(idx);
                idx2 = idx2_all;
            else
                idx1 = idx2_all;
                idx2 = idx1_all(idx);
            end
            
        elseif age_match == false
            
            idx1 = find(sex_vec==1 & grp_vec==i & age_vec>ageRng(j)-window/2 & age_vec<ageRng(j)+window/2);
            idx2 = find(sex_vec==2 & grp_vec==i & age_vec>ageRng(j)-window/2 & age_vec<ageRng(j)+window/2);
            
        end
        
        rmvAvg(j,1) = squeeze(mean(gain(idx1)));
        rmvAvg(j,2) = squeeze(mean(gain(idx2)));
        
        if length(idx1)>1
            ci = bootci(nperm,@mean,gain(idx1));
            err(j,1,1) = rmvAvg(j,1)-ci(1);
            err(j,2,1) = ci(2)-rmvAvg(j,1); clear ci
            ci = bootci(nperm,@mean,gain(idx2));
            err(j,1,2) = rmvAvg(j,2)-ci(1);
            err(j,2,2) = ci(2)-rmvAvg(j,2); clear ci
        else
            disp('ah for feck sake')
            err(j,1,1) = 0;
            err(j,2,1) = 0;
            err(j,1,2) = 0;
            err(j,2,2) = 0;
        end
        
        % Run permutation test
        if length(idx1)>1
            [~,~,orig] = permtest_t2(gain(idx1),gain(idx2),'nperm',nperm);
            p(j) = orig.p;
            h(j) = orig.h;
        else
            p(j) = 1;
            h(j) = 0;
        end
        
    end
    
    hx = fdr_bky(p);
    
    if i == 1
        h1 = subplot(3,4,9:10);
    else
        h1 = subplot(3,4,11:12);
    end
    hold all
    plot(ageRng,rmvAvg(:,2)+1e3,'color',cols1(5,:),'LineWidth',lnWth)
    plot(ageRng,rmvAvg(:,1)+1e3,'color',cols1(4,:),'LineWidth',lnWth)
    xaxis = linspace(ageRng(1),ageRng(end),length(ageRng));
    for k = 2
        if k == 1
            hCorr = h;
            shade = 0.925;
        else
            hCorr = hx;
            shade = 0.85;
        end
        try
            sigX = xaxis(hCorr==1);
            for l = 1:length(sigX)
                X = [sigX(l)-0.5,sigX(l)-0.5,sigX(l)+0.5,sigX(l)+0.5];
                Y = [y1,y2,y2,y1];
                C = ones(1,3)*shade;
                fill(X,Y,C,'edgecolor','none')
            end
        catch
        end
    end
    plot([x1,x2],[0,0],'--','color',[0.6,0.6,0.6])
    [hl,hp] = boundedline(ageRng,rmvAvg(:,2),err(:,:,2),ageRng,rmvAvg(:,1),...
        err(:,:,1),'alpha','cmap',cols1(1:2,:));
    outlinebounds(hl,hp);
    plot(ageRng,rmvAvg(:,2),'color',cols1(5,:),'LineWidth',lnWth)
    plot(ageRng,rmvAvg(:,1),'color',cols1(4,:),'LineWidth',lnWth)
    set(gca,'layer','top','box','off','tickdir','out',...
        'ytick',-0.05:0.05:0.05,'yticklabel',{'-.05','0','.05'})
    title(diag_groups{i},'FontName',fnt,'FontSize',fntSz)
    xlabel('Age (yrs)','FontName',fnt,'FontSize',fntSz)
    if i==1
        ylabel('Violation','FontName',fnt,'FontSize',fntSz)
        text(x1-5,y2+0.005,'C','FontName',fnt,'FontSize',fntSz+8)
    elseif i==2
        l = legend(sex_groups);
        set(l,'box','off','location','northwest','FontName',fnt,'FontSize',fntSz-2)
        set(gca,'yticklabel',[])
    end
    xlim([x1,x2])
    ylim([y1,y2])
    p = get(h1,'position');
    set(h1,'position',[p(1),p(2)-0.03,p(3),p(4)])
    
end

% % Plot x-axis scale
% h = subplot(3,4,9);
% set(gca,'layer','top','box','off','ycolor','none','color','none')
% set(gca,'tickdir','both','xtick',0:0.5:1,'xticklabel',{'0','.5','1'})
% xlabel('Quantile')
% axis square
% pos = get(h,'position');
% posnew = pos; posnew(2) = posnew(2)+0.26; set(h,'position',posnew)

set(gcf,'color','w')

%% Plot MDS

close all;
clc;

% close all;
clc;

cond = 7;
task = 6:8;
% nVars = 11;
nVars = 29;
idx = grp_vec==1 | grp_vec==2;
data = zeros(sum(idx),nVars);
quant = 2:length(prob)-1;
xaxis = prob/100;

mseAV = squeeze(trapz(xaxis,mean(CDFs(:,1,2,:)-CDFs(:,1,3,:),2),4));
mseA = squeeze(trapz(xaxis,mean(CDFs(:,2,2,:)-CDFs(:,2,3,:),2),4));
mseV = squeeze(trapz(xaxis,mean(CDFs(:,3,2,:)-CDFs(:,3,3,:),2),4));
% rmv = squeeze(trapz(xaxis,mean(CDFs(:,cond,task,:),3),4));
rmv = squeeze(mean(CDFs(:,cond,task,quant),3));

data(:,1) = age_vec(idx);
data(:,2) = piq_vec(idx);
data(:,3) = viq_vec(idx);
data(:,4) = fsiq_vec(idx);
data(:,5) = Fscore(idx);
data(:,6) = mseAV(idx);
data(:,7) = mseA(idx);
data(:,8) = mseV(idx);
data(:,9) = racefit(idx,1);
data(:,10) = racefit(idx,2);
% data(:,11) = rmv(idx);
data(:,11:29) = rmv(idx,:);

nDims = 3;
x1 = cell(1,5);
idx = cell(1,5);

dist = 'euclidean'; % 'euclidean' 'seuclidean' 'correlation' 'spearman'
start = 'random'; % 'cmdscale' 'random'
crit = 'stress'; % 'stress' 'sstress' 'metricstress' 'metricsstress'

for i = 1:5
    
    if i==5
        idx{i} = 1:nSubjs;
    else
        idx{i} = dev_vec==i;
    end
    X1 = data(idx{i},:);
    
    % Compute Euclidean distances between each pair of subjects
    D1 = pdist(X1,dist);
    
    % Perform non-metric MDS in nFeats dimenstions to evaluate how many
    % dimensitons required to reduce Kruskals stress < 0.1
    stress = zeros(1,5);
    for j = 1:5
        [~,stress(j)] = mdscale(D1,j,'criterion',crit,'start',start);
    end
    figure, plot([0,6],[0.1,0.1],'--k')
    hold on, plot(stress,'-or','linewidth',3), xlim([0,6]),ylim([0,1])
    xlabel('Dimensions'); ylabel('Kruskalls Stress');
    
    %     % Perform non-metric MDS, incrementing dimensions until Kruskals stress < 0.1
    %     for j = 1:5
    %         [x1{i},stress] = mdscale(D1,j,'criterion',crit,'start',start);
    %         if stress < 0.1
    %             break
    %         end
    %     end
    
    % Perform non-metric MDS in 2 dimensions
    [x1{i},~] = mdscale(D1,nDims,'criterion',crit,'start',start);
    
    clear D1 stress1
    
end

%% Generate figure

close all;
nDims = 3;

axes('parent',figure,'fontname',fnt,'fontsize',fntSz)
set(gcf,'position',[20,50,1200,300],'paperpositionmode','auto')

for i = 1:5
    
    if i==5
        idx{i} = 1:nSubjs;
    else
        idx{i} = find(dev_vec==i);
    end
    
    subplot(1,5,i)
    if nDims==2
        scatter(x1{i}(:,1),x1{i}(:,2),'MarkerEdgeColor','none');
        for j = 1:length(x1{i})
            text(x1{i}(j,1),x1{i}(j,2),['S',num2str(j)],'color',cols1(grp_vec(idx{i}(j)),:),'horizontal','left',...
                'vertical','bottom','horizontalalign','center','verticalalign',...
                'middle','FontName',fnt,'FontSize',8,'Fontweight','bold')
        end
    elseif nDims==3
        scatter3(x1{i}(:,1),x1{i}(:,2),x1{i}(:,3),'MarkerEdgeColor','none');
        for j = 1:length(x1{i})
            text(x1{i}(j,1),x1{i}(j,2),x1{i}(j,3),['S',num2str(j)],'color',cols1(grp_vec(idx{i}(j)),:),'horizontal','left',...
                'vertical','bottom','horizontalalign','center','verticalalign',...
                'middle','FontName',fnt,'FontSize',8,'Fontweight','bold')
        end
        zlim([min(x1{i}(:,3)),max(x1{i}(:,3))])
        
    end
    axis square
    
    if i==3
        xlabel('MDS 1')
    elseif i==1
        ylabel('MDS 2')
    end
    
    %     xlim([-1,1])
    %     ylim([-0.4,0.4])
    xlim([min(x1{i}(:,1)),max(x1{i}(:,1))])
    ylim([min(x1{i}(:,2)),max(x1{i}(:,2))])
    
end

set(gcf,'color','w')

%% Mediation analysis

close all;
clc;

xaxis = prob/100;

% Define race model violation
rmv = squeeze(CDFs(:,1,1,:)-mean(CDFs(:,6,6:8,:),3));
gain = squeeze(trapz(xaxis,rmv,2));
mse = squeeze(trapz(xaxis,mean(CDFs(:,1:3,2,:)-CDFs(:,1:3,3,:),2),4));

% Define variables
subj = 1:nSubjs;
grp = grp_vec;
sex = sex_vec;
age = age_vec;
dev = dev_vec;
rmv = gain;

% Convert non-numeric variables to categorical
subj = categorical(subj','Ordinal',0);
grp = categorical(grp,[1,2],{'TD','ASD'},'Ordinal',0);
sex = categorical(sex,[1,2],{'Male','Fem'},'Ordinal',0);
dev = categorical(dev,[1,2,3,4],{'6-9','10-12','13-17','18-40'},'Ordinal',1);

% Create dataset
ds = dataset(subj,grp,sex,age,dev,rmv,mse);

lm1 = fitlm(ds,'mse ~ age');
lm2 = fitlm(ds,'rmv ~ age + mse');
lm3 = fitlm(ds,'rmv ~ age');

paths = [lm1.Coefficients{2,1},lm2.Coefficients{3,1},lm2.Coefficients{2,1},lm3.Coefficients{2,1},lm1.Coefficients{2,1}*lm2.Coefficients{3,1}]
SE = [lm1.Coefficients{2,2},lm2.Coefficients{3,2},lm2.Coefficients{2,2},lm3.Coefficients{2,2},0]
pVals = [lm1.Coefficients{2,4},lm2.Coefficients{3,4},lm2.Coefficients{2,4},lm3.Coefficients{2,4},0]

