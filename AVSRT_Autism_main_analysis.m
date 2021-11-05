% Script for running the main analysis of the AVSRT-Autism Study.

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

clear;
clc;

% Directory
direc = '.\Data\AVSRT-Autism\';
subj_info = xlsread([direc,'Subject Info\AVSRT_ID_Sex_Age.xlsx']);

% Subject groups
subj_group = {'TD Child','TD Adult','ASD Child','ASD Adult'};

% Set parameters
nSubjs = 411;
minAge = 6;
maxAge = 40;
minIQ = 80;
minRT = 100;
maxRT = 2e3;
minISI = 950;
maxISI = 3150;
minPer = 2.5;
maxPer = 97.5;
minTrials = 20;
startTrials = 3;
factor = 3;
outliers = 'cond';
testDate = false;

% Preallocate memory
faRate = zeros(nSubjs,1);
hitRate = zeros(nSubjs,1);
Fscore = zeros(nSubjs,1);
AHR = zeros(nSubjs,1);
VHR = zeros(nSubjs,1);
perLim = zeros(nSubjs,2);
numRTs = zeros(nSubjs,3);
lwrISI = zeros(nSubjs,1);
testYear = zeros(nSubjs,1);

grp_vec = zeros(nSubjs,1);
sex_vec = zeros(nSubjs,1);
age_vec = zeros(nSubjs,1);
dev_vec = zeros(nSubjs,1);
piq_vec = zeros(nSubjs,1);
viq_vec = zeros(nSubjs,1);
fsiq_vec = zeros(nSubjs,1);
ados_vec = zeros(nSubjs,1);

RT_cell = cell(nSubjs,1);
RTprev_cell = cell(nSubjs,1);
cond_cell = cell(nSubjs,1);
task_cell = cell(nSubjs,1);
type_cell = cell(nSubjs,1);
ISI_cell = cell(nSubjs,1);
subj_cell = cell(nSubjs,1);
grp_cell = cell(nSubjs,1);
sex_cell = cell(nSubjs,1);
age_cell = cell(nSubjs,1);
dev_cell = cell(nSubjs,1);
piq_cell = cell(nSubjs,1);
viq_cell = cell(nSubjs,1);
fsiq_cell = cell(nSubjs,1);
ados_cell = cell(nSubjs,1);

% Initialize variables
ctr = 1;
ids = [];
allRTs = [];
fastRTs = [];
slowRTs = [];

% Process data
fprintf('Processing data')

% Loop through groups
for i = 1:length(subj_group), fprintf('.')
    
    % Get subject IDs
    id = dir([direc,'Behaviour\',subj_group{i}]); id(1:2) = [];
    
    % Define diagnosis group (TD=1, ASD=2)
    switch subj_group{i}
        case subj_group(1:2)
            grp = 1;
        case subj_group(3:4)
            grp = 2;
    end
    
    % Loop through subjects
    for j = 1:size(id,1)
        
        % Load data
        load([direc,'Behaviour\',subj_group{i},'\',id(j).name,'\',...
            id(j).name,'_rt\',id(j).name,'_RTs.mat'],'RTs','R2s','ISIs',...
            'tasks','conds','nStims','nResps','nResps2');
        
        % Get testing date
        if testDate == true
            files = dir(['K:\Data\AVSRT\',subj_group{i},'\',id(j).name,'\*.bdf']);
            testYear(ctr) = str2double(files(1).date(8:11));
        end
        
        % Compute some parameters
        nR2s = nansum(R2s);
        nBlocks = length(nStims);
        nTrials = sum(nStims);
        
        % Compute percentage of fast and slow RTs
        allRTs = [allRTs;RTs(~isnan(RTs) & RTs<4e3)];
        fastRTs = [fastRTs;sum(RTs<minRT)/length(~isnan(RTs))*100];
        slowRTs = [slowRTs;sum(RTs>maxRT)/length(~isnan(RTs))*100];
        
        % Define previous RTs
        RTprev = [NaN;RTs(1:end-1)];
        RTprev(nStims(1:end-1)+1) = NaN;
        
        % Re-number conditions (AV=1, A=2, V=3)
        conds = conds-2;
        
        % Find number of hits
        noRT = isnan(RTs) | RTs>maxRT;
        misses = sum(noRT);
        hits = nTrials-misses;
        aHits = sum(conds==2 & noRT==0);
        vHits = sum(conds==3 & noRT==0);
        
        % Compute precision and recall
        falarms = sum(nResps)-hits+sum(RTs<minRT);
        precis = hits./(hits+falarms);
        recall = hits./(hits+misses);
        
        % Compute false alarm rate and hit rate
        faRate(ctr) = falarms/nTrials*100;
        hitRate(ctr) = hits/nTrials*100;
        Fscore(ctr) = 2*(precis.*recall)./(precis+recall);
        
        % Compute A and V hit rates as percentage of each other
        AHR(ctr) = aHits/vHits*100;
        VHR(ctr) = vHits/aHits*100;
        
        % Index first 3 trials of each block (training trials)
        startTrials = [];
        endTrials = [0;cumsum(nStims(1:end-1))];
        for k = 1:nBlocks
            startTrials = [startTrials;endTrials(k)+(1:startTrials)'];
        end
        
        % Index bad trials (misses, double-presses, fast/slow RTs, training)
        badRTs = unique([find(isnan(RTs));find(RTs<minRT|RTs>maxRT);...
            find(R2s==true);find(ISIs<minISI|ISIs>maxISI);startTrials]);
        
        % Find bad previous trials
        badprevRTs = badRTs(badRTs<length(RTs))+1;
        RTprev(badprevRTs) = NaN;
        
        % Get rid of bad trials
        RTs(badRTs) = []; %#ok<*SAGROW>
        RTprev(badRTs) = [];
        ISIs(badRTs) = [];
        tasks(badRTs) = [];
        conds(badRTs) = [];
        clear badRTs badprevRTs
        
        % Get raw data for example subject
        if strcmp(id(j).name,'12161')
            spcloneRTs = RTs;
            spcloneConds = conds;
        end
        
        % Use middle 95th percentile of RTs
        switch outliers
            case 'cond'
                perCond = zeros(3,2);
                idxConds = cell(1,3);
                for k = 1:3
                    idxConds{k} = find(conds==k);
                    perCond(k,:) = prctile(RTs(conds==k),[minPer,maxPer]);
                    badRTs{k} = idxConds{k}(RTs(conds==k)<perCond(k,1) | RTs(conds==k)>perCond(k,2));
                end
                perLim(ctr,:) = [min(perCond(:,1)),max(perCond(:,2))];
                badRTs = cell2mat(badRTs(:));
            case 'all'
                perLim(ctr,:) = prctile(RTs,[minPer,maxPer]);
                badRTs = find(RTs<perLim(ctr,1) | RTs>perLim(ctr,2));
        end
        
        % Find bad previous trials
        badprevRTs = badRTs(badRTs<length(RTs))+1;
        RTprev(badprevRTs) = NaN;
        
        % Get rid of bad trials
        RTs(badRTs) = [];
        RTprev(badRTs) = [];
        ISIs(badRTs) = [];
        tasks(badRTs) = [];
        conds(badRTs) = [];
        
        % Compute # of RTs per condition
        for k = 1:3
            numRTs(ctr,k) = sum(conds==k);
        end
        
        % Compute average ISI
        lwrISI(ctr) = nanmin(ISIs);
        
        % Get subject sex and age
        sex = subj_info((subj_info(:,1)==str2double(id(j).name)),3);
        age = subj_info((subj_info(:,1)==str2double(id(j).name)),4);
        piq = max(subj_info((subj_info(:,1)==str2double(id(j).name)),5));
        viq = subj_info((subj_info(:,1)==str2double(id(j).name)),6);
        fsiq = subj_info((subj_info(:,1)==str2double(id(j).name)),7);
        ados = subj_info((subj_info(:,1)==str2double(id(j).name)),8);
        
        % Define developmental group (6-9=1, 10-12=2, 13-17=3, 18-40=4)
        if age>=minAge && age<10
            dev = 1;
        elseif age>=10 && age<13
            dev = 2;
        elseif age>=13 && age<18
            dev = 3;
        elseif age>=18 && age<maxAge
            dev = 4;
        else
            dev = NaN;
        end
        
        % Define task type (swtich=1, repeat=2, neither=NaN)
        type = tasks;
        type(tasks==2|tasks==3|tasks==6|tasks==8) = 1;
        type(tasks==1|tasks==5|tasks==9) = 2;
        type(tasks==4|tasks==7) = NaN;
        
        % Store data in vectors
        sex_vec(ctr) = sex;
        grp_vec(ctr) = grp;
        dev_vec(ctr) = dev;
        age_vec(ctr) = age;
        piq_vec(ctr) = piq;
        viq_vec(ctr) = viq;
        fsiq_vec(ctr) = fsiq;
        ados_vec(ctr) = ados;
        
        % Store data in cell arrays
        RT_cell{ctr} = RTs;
        RTprev_cell{ctr} = RTprev;
        ISI_cell{ctr} = ISIs;
        cond_cell{ctr} = conds;
        task_cell{ctr} = tasks;
        type_cell{ctr} = type;
        subj_cell{ctr} = ctr*ones(size(RTs));
        grp_cell{ctr} = grp*ones(size(RTs));
        sex_cell{ctr} = sex*ones(size(RTs));
        age_cell{ctr} = age*ones(size(RTs));
        dev_cell{ctr} = dev*ones(size(RTs));
        piq_cell{ctr} = piq*ones(size(RTs));
        viq_cell{ctr} = viq*ones(size(RTs));
        fsiq_cell{ctr} = fsiq*ones(size(RTs));
        ados_cell{ctr} = ados*ones(size(RTs));
        
        % Increment counter
        ctr = ctr+1;
        
        clear RTs RTprev ISIs R2s badRTs
        clear conds tasks type trials
        clear nStims nResps nResps2
        clear sex age dev piq viq fsiq ados
        
    end
    
    ids = [ids;id];
    clear id group
    
end

% Update number of subjects
nSubjs = ctr-1;

% Store variables for plotting
tmp1 = faRate;
tmp2 = Fscore;
midRTs = cell2mat(RT_cell);

% Index bad subjects
fprintf('\n\nBad subjects:\n')

% Identify subjects too young/old
badIDs = isnan(age_vec) | age_vec<minAge | age_vec>maxAge;
for i = find(badIDs)'
    fprintf('%s - Age: %.1f yrs\n',ids(i).name,age_vec(i))
end
idx = badIDs;

% Identify subjects with low PIQ
tmp = piq_vec; tmp(idx) = NaN;
badIDs = tmp<minIQ | (isnan(piq_vec) & grp_vec==2 & dev_vec~=4);
for i = find(badIDs)'
    fprintf('%s - PIQ: %d\n',ids(i).name,piq_vec(i))
end
idx = idx | badIDs; clear tmp

% Identify subjects with high false alarm rate (excessive button presses)
tmp = faRate; tmp(idx) = NaN;
thr = nanmean(tmp)+factor*nanstd(tmp);
badIDs = tmp>thr;
for i = find(badIDs)'
    fprintf('%s - FA: %.0f%%\n',ids(i).name,faRate(i))
end
idx = idx | badIDs; clear tmp thr

% Identify subjects with low hit rate (not attending)
tmp = Fscore; tmp(idx) = NaN;
thr = nanmean(tmp)-factor*nanstd(tmp);
badIDs = tmp<thr;
for i = find(badIDs)'
    fprintf('%s - F1: %.2f\n',ids(i).name,Fscore(i))
end
idx = idx | badIDs; clear tmp thr

% Identify subjects with low V hit rate (eyes closed)
tmp = VHR; tmp(idx) = NaN;
tmp(tmp>1000) = NaN;
thr = nanmean(tmp)-factor*nanstd(tmp);
badIDs = tmp<thr;
for i = find(badIDs)'
    fprintf('%s - V-HR: %.0f%%\n',ids(i).name,VHR(i))
end
idx = idx | badIDs; clear tmp thr

% Identify subjects with low A hit rate (not listening)
tmp = AHR; tmp(idx) = NaN;
tmp(tmp>1000) = NaN;
thr = nanmean(tmp)-factor*nanstd(tmp);
badIDs = tmp<thr;
for i = find(badIDs)'
    fprintf('%s - A-HR: %.0f%%\n',ids(i).name,AHR(i))
end
idx = idx | badIDs; clear tmp thr

% Identify subjects with < 20 RTs (min required number)
tmp = numRTs; tmp(idx,:) = NaN;
badIDs = min(tmp,[],2)<minTrials;
for i = find(badIDs)'
    fprintf('%s - #RTs: %d, %d, %d\n',ids(i).name,numRTs(i,:)<minTrials)
end
idx = idx | badIDs; clear tmp thr

% Identify subjects with min ISIs > 2000 ms (should be 1000-3000 ms)
tmp = lwrISI; tmp(idx,:) = NaN;
badIDs = tmp>2e3;
for i = find(badIDs)'
    fprintf('%s - Min ISI: %.0f ms\n',ids(i).name,lwrISI(i))
end
idx = idx | badIDs; clear tmp thr badIDs

% Compute total bad subjects
nBad = sum(idx)-length(ctr:nSubjs);
fprintf('Subjects rejected: %d/%d (%.1f%%)\n',nBad,nSubjs,nBad/nSubjs*100)

% Get rid of bad subjects
ids(idx) = [];
grp_vec(idx) = [];
sex_vec(idx) = [];
age_vec(idx) = [];
dev_vec(idx) = [];
piq_vec(idx) = [];
viq_vec(idx) = [];
fsiq_vec(idx) = [];
ados_vec(idx) = [];
faRate(idx) = [];
hitRate(idx) = [];
Fscore(idx) = [];
AHR(idx) = [];
VHR(idx) = [];
numRTs(idx,:) = [];
fastRTs(idx) = [];
slowRTs(idx) = [];
perLim(idx,:) = [];
testYear(idx) = [];
lwrISI(idx) = [];

subj_cell(idx) = [];
grp_cell(idx) = [];
sex_cell(idx) = [];
age_cell(idx) = [];
dev_cell(idx) = [];
piq_cell(idx) = [];
viq_cell(idx) = [];
fsiq_cell(idx) = [];
ados_cell(idx) = [];
RT_cell(idx) = [];
RTprev_cell(idx) = [];
ISI_cell(idx) = [];
task_cell(idx) = [];
type_cell(idx) = [];
cond_cell(idx) = [];

% Update number of subjects
nSubjs = length(ids);
exSub = find(strcmp({ids.name},'12161'));

%% Plot results of outlier correction procedure

close all;
clc;

% Print percentage of fast/slow outliers
fprintf('Percentage of fast outliers: %1.1f ± %1.1f\n',mean(fastRTs),std(fastRTs));
fprintf('Percentage of slow outliers: %1.1f ± %1.1f\n',mean(slowRTs),std(slowRTs));

% Plot all and middle 95th percentile RTs
figure
subplot(2,1,1), histogram(allRTs,150), title('Before'), xlim([0,maxRT])
subplot(2,1,2), histogram(midRTs,100), title('After'), xlim([0,maxRT])

% Plot response rate of subjects within 95th percentile
figure
grp = [zeros(size(tmp1')),ones(size(faRate'))];
subplot(2,3,[1,4]), boxplot([tmp1;faRate]',grp)
title('False Alarm Rate'), ylabel('%'), set(gca,'xticklabel',{'Before','After'})
subplot(2,3,2:3), histogram(tmp1,60), title('Before'), xlim([0,250])
subplot(2,3,5:6), histogram(faRate,20), title('After'), xlim([0,250])

% Plot hit rate of subjects within 95th percentile
figure
grp = [zeros(size(tmp2')),ones(size(Fscore'))];
subplot(2,3,[1,4]), boxplot([tmp2;Fscore]',grp)
title('F_1 score'), ylabel('%'), set(gca,'xticklabel',{'Before','After'})
subplot(2,3,2:3), histogram(tmp2,40), title('Before'), xlim([0,1])
subplot(2,3,5:6), histogram(Fscore,20), title('After'), xlim([0,1])

%% Compute participant demographics

clc;

% Total number of subjects
fprintf('Number of subjects: %d (%d F)\n',length(ids),sum(sex_vec==2));
fprintf('Number of NTs: %d (%d F)\n',sum(grp_vec==1),sum(grp_vec==1 & sex_vec==2));
fprintf('Number of ASDs: %d (%d F)\n',sum(grp_vec==2),sum(grp_vec==2 & sex_vec==2));

% Age ranges
fprintf('\nNT age range: %1.1f-%1.1f yrs\n',min(age_vec(grp_vec==1)),max(age_vec(grp_vec==1)));
fprintf('ASD age range: %1.1f-%1.1f yrs\n',min(age_vec(grp_vec==2)),max(age_vec(grp_vec==2)));

% Number of subjects
fprintf('\nNumber of subjects:\n')
for i = 1:2
    for j = 1:4
        n = sum(grp_vec==i & dev_vec==j);
        fprintf('%d\t',n);
    end
    fprintf('\n');
end

% Number of females
fprintf('\nNumber of females:\n')
for i = 1:2
    for j = 1:4
        idx = grp_vec==i & dev_vec==j;
        fprintf('%d\t',length(find(sex_vec(idx)==2)));
    end
    fprintf('\n');
end

% Number of IQs
fprintf('\nNumber of IQs:\n')
for i = 1:2
    for j = 1:4
        idx = grp_vec==i & dev_vec==j;
        fprintf('%d\t',sum(~isnan(piq_vec(idx))));
    end
    fprintf('\n');
end

% Number of IQs
fprintf('\nPercentage of IQs:\n')
for i = 1:2
    idx = grp_vec==i;
    fprintf('%.1f\t',sum(~isnan(piq_vec(idx)))/sum(idx)*100);
    fprintf('\n');
end

% Mean age
fprintf('\nMean age (SD):\n')
for i = 1:2
    for j = 1:4
        idx = grp_vec==i & dev_vec==j;
        fprintf('%1.1f (%1.1f)\t',nanmean(age_vec(idx)),nanstd(age_vec(idx)));
    end
    fprintf('\n');
end

% Mean F1 score
fprintf('\nMean F1 score (SD):\n')
for i = 1:2
    for j = 1:4
        idx = find(grp_vec==i & dev_vec==j);
        fprintf('%1.2f (%1.2f)\t',nanmean(Fscore(idx)),nanstd(Fscore(idx)));
    end
    fprintf('\n');
end

% Mean False alarm rate
fprintf('\nMean False alarm rate (SD):\n')
for i = 1:2
    for j = 1:4
        idx = find(grp_vec==i & dev_vec==j);
        fprintf('%1.2f (%1.2f)\t',nanmean(faRate(idx)/100),nanstd(faRate(idx)/100));
    end
    fprintf('\n');
end

% Mean misses
fprintf('\nMean misses (SD):\n')
for i = 1:2
    for j = 1:4
        idx = find(grp_vec==i & dev_vec==j);
        fprintf('%1.2f (%1.2f)\t',nanmean(1-hitRate(idx)/100),nanstd(1-hitRate(idx)/100));
    end
    fprintf('\n');
end

% Mean PIQ
fprintf('\nMean PIQ (SD):\n')
for i = 1:2
    for j = 1:4
        idx = grp_vec==i & dev_vec==j;
        fprintf('%1.1f (%1.1f)\t',nanmean(piq_vec(idx)),nanstd(piq_vec(idx)));
    end
    fprintf('\n');
end

% Mean VIQ
fprintf('\nMean VIQ (SD):\n')
for i = 1:2
    for j = 1:4
        idx = grp_vec==i & dev_vec==j;
        fprintf('%1.1f (%1.1f)\t',nanmean(viq_vec(idx)),nanstd(viq_vec(idx)));
    end
    fprintf('\n');
end

% Mean FSIQ
fprintf('\nMean FSIQ (SD):\n')
for i = 1:2
    for j = 1:4
        idx = grp_vec==i & dev_vec==j;
        fprintf('%1.1f (%1.1f)\t',nanmean(fsiq_vec(idx)),nanstd(fsiq_vec(idx)));
    end
    fprintf('\n');
end

% Mean ADOS score
fprintf('\nMean ADOS score (SD):\n')
for i = 1:2
    for j = 1:4
        idx = grp_vec==i & dev_vec==j;
        fprintf('%1.1f (%1.1f), n = %d\t',nanmean(ados_vec(idx)),nanstd(ados_vec(idx)),sum(~isnan(ados_vec(idx))));
    end
    fprintf('\n');
end

%% Test for group differences in age, performance and IQ

close all;
clc;

% Set params
nperm = 1e4;
age_match = false;

% Test for age difference
fprintf('Age:\n')
for i = 1:4
    
    idx1 = grp_vec==1 & dev_vec==i & ~isnan(age_vec);
    idx2 = grp_vec==2 & dev_vec==i & ~isnan(age_vec);
    
    [tstat,~,orig,stats] = permtest_t2(age_vec(idx1),age_vec(idx2),'nperm',nperm);
    d = mes(age_vec(idx1),age_vec(idx2),'hedgesg','isDep',0,'nBoot',nperm);
    fprintf('t(%d) = %0.2f, ',stats.df,tstat)
    fprintf('p = %0.3f, ',orig.p)
    fprintf('g = %0.2f, ',d.hedgesg)
    fprintf('95CI [%0.2f, %0.2f]\n',d.hedgesgCi)
    
end

% Test for F1 score difference
fprintf('\nF1 score:\n')
for i = 1:4
    
    switch age_match
        case true
            if i==4
                [idx1,idx2] = agematch(i,grp_vec,dev_vec,sex_vec,age_vec,[],false);
            else
                [idx1,idx2] = agematch(i,grp_vec,dev_vec,sex_vec,age_vec,piq_vec,false);
            end
        case false
            idx1 = find(grp_vec==1 & dev_vec==i & ~isnan(Fscore));
            idx2 = find(grp_vec==2 & dev_vec==i & ~isnan(Fscore));
    end
    
    [tstat,~,orig,stats] = permtest_t2(Fscore(idx1),Fscore(idx2),'nperm',nperm);
    d = mes(Fscore(idx1),Fscore(idx2),'hedgesg','isDep',0,'nBoot',nperm);
    fprintf('t(%d) = %0.2f, ',stats.df,tstat)
    fprintf('p = %0.3f, ',orig.p)
    fprintf('g = %0.2f, ',d.hedgesg)
    fprintf('95CI [%0.2f, %0.2f]\n',d.hedgesgCi)
    
end

% Test for PIQ difference
fprintf('\nPIQ:\n')
for i = 1:4
    
    switch age_match
        case true
            [idx1,idx2] = agematch(i,grp_vec,dev_vec,sex_vec,age_vec,piq_vec,false);
        case false
            idx1 = find(grp_vec==1 & dev_vec==i & ~isnan(piq_vec));
            idx2 = find(grp_vec==2 & dev_vec==i & ~isnan(piq_vec));
    end
    
    [tstat,~,orig,stats] = permtest_t2(piq_vec(idx1),piq_vec(idx2),'nperm',nperm);
    d = mes(piq_vec(idx1),piq_vec(idx2),'hedgesg','isDep',0,'nBoot',nperm);
    fprintf('t(%d) = %0.2f, ',stats.df,tstat)
    fprintf('p = %0.3f, ',orig.p)
    fprintf('g = %0.2f, ',d.hedgesg)
    fprintf('95CI [%0.2f, %0.2f]\n',d.hedgesgCi)
    
end

% Test for VIQ difference
fprintf('\nVIQ:\n')
for i = 1:4
    
    switch age_match
        case true
            [idx1,idx2] = agematch(i,grp_vec,dev_vec,sex_vec,age_vec,piq_vec,false);
        case false
            idx1 = find(grp_vec==1 & dev_vec==i & ~isnan(viq_vec));
            idx2 = find(grp_vec==2 & dev_vec==i & ~isnan(viq_vec));
    end
    
    [tstat,~,orig,stats] = permtest_t2(viq_vec(idx1),viq_vec(idx2),'nperm',nperm);
    d = mes(viq_vec(idx1),viq_vec(idx2),'hedgesg','isDep',0,'nBoot',nperm);
    fprintf('t(%d) = %0.2f, ',stats.df,tstat)
    fprintf('p = %0.3f, ',orig.p)
    fprintf('g = %0.2f, ',d.hedgesg)
    fprintf('95CI [%0.2f, %0.2f]\n',d.hedgesgCi)
    
end

% Test for FSIQ difference
fprintf('\nFSIQ:\n')
for i = 1:4
    
    switch age_match
        case true
            [idx1,idx2] = agematch(i,grp_vec,dev_vec,sex_vec,age_vec,piq_vec,true);
        case false
            idx1 = find(grp_vec==1 & dev_vec==i & ~isnan(fsiq_vec));
            idx2 = find(grp_vec==2 & dev_vec==i & ~isnan(fsiq_vec));
    end
    
    [tstat,~,orig,stats] = permtest_t2(fsiq_vec(idx1),fsiq_vec(idx2),'nperm',nperm);
    d = mes(fsiq_vec(idx1),fsiq_vec(idx2),'hedgesg','isDep',0,'nBoot',nperm);
    fprintf('t(%d) = %0.2f, ',stats.df,tstat)
    fprintf('p = %0.3f, ',orig.p)
    fprintf('g = %0.2f, ',d.hedgesg)
    fprintf('95CI [%0.2f, %0.2f]\n',d.hedgesgCi)
    
end

%% Compute CDFs and RMV

nConds = 4;
nConds2 = 11;
nTasks = 8;
test = 'ver';
dep = 0;
nSamps = 1e4;

% Define CDF percentiles
prob = 0:0.05:1;

% Pre-allocate memory
RTper = zeros(nSubjs,numel(prob));
RTs = cell(nSubjs,nConds,nTasks);

% Separate RTs by condition and trial type
for i = 1:nSubjs
    RTper(i,:) = linspace(perLim(i,1),perLim(i,2),numel(prob));
    for j = 1:3
        RTs{i,j,1} = RT_cell{i}(cond_cell{i}==j); % 1. mix
        RTs{i,j,2} = RT_cell{i}((task_cell{i}==1 | task_cell{i}==5 | task_cell{i}==9) & cond_cell{i}==j); % 2. repeat
        RTs{i,j,3} = RT_cell{i}((task_cell{i}==2 | task_cell{i}==3 | task_cell{i}==6 | task_cell{i}==8) & cond_cell{i}==j); % 3. switch
        RTs{i,j,4} = RT_cell{i}((task_cell{i}==2 | task_cell{i}==6 | task_cell{i}==8) & cond_cell{i}==j); % 4. switch A->AV
        RTs{i,j,5} = RT_cell{i}((task_cell{i}==3 | task_cell{i}==6 | task_cell{i}==8) & cond_cell{i}==j); % 5. switch V->AV
        RTs{i,j,6} = RT_cell{i}((task_cell{i}==1 | task_cell{i}==4 | task_cell{i}==7) & cond_cell{i}==j); % 6. AV prior
        RTs{i,j,7} = RT_cell{i}((task_cell{i}==2 | task_cell{i}==5 | task_cell{i}==8) & cond_cell{i}==j); % 7. A prior
        RTs{i,j,8} = RT_cell{i}((task_cell{i}==3 | task_cell{i}==6 | task_cell{i}==9) & cond_cell{i}==j); % 8. V prior
        RTs{i,j,9} = RT_cell{i}((task_cell{i}==1 | task_cell{i}==5 | task_cell{i}==9) & cond_cell{i}==j & ISI_cell{i}<1500); % 9. repeat fast ISI
        RTs{i,j,10} = RT_cell{i}((task_cell{i}==2 | task_cell{i}==3 | task_cell{i}==6 | task_cell{i}==8) & cond_cell{i}==j & ISI_cell{i}<1500); % 10. switch fast ISI
        RTs{i,j,11} = RT_cell{i}((task_cell{i}==1 | task_cell{i}==5 | task_cell{i}==9) & cond_cell{i}==j & ISI_cell{i}>2500); % 11. repeat slow ISI
        RTs{i,j,12} = RT_cell{i}((task_cell{i}==2 | task_cell{i}==3 | task_cell{i}==6 | task_cell{i}==8) & cond_cell{i}==j & ISI_cell{i}>2500); % 12. switch slow ISI
        RTs{i,j,13} = RT_cell{i}(cond_cell{i}==j & ISI_cell{i}>2500); % 13. mix slow ISI
    end
end

% Simulate race distribution by randomly sampling from A & V RTs
% simultaneously and taking minimum value
for i = 1:nSubjs
    
    RTs{i,4,1} = min([RTs{i,2,1}(randi(length(RTs{i,2,1}),nSamps,1)),RTs{i,3,1}(randi(length(RTs{i,3,1}),nSamps,1))],[],2);
    RTs{i,4,2} = min([RTs{i,2,2}(randi(length(RTs{i,2,2}),nSamps,1)),RTs{i,3,2}(randi(length(RTs{i,3,2}),nSamps,1))],[],2);
    RTs{i,4,3} = min([RTs{i,2,3}(randi(length(RTs{i,2,3}),nSamps,1)),RTs{i,3,3}(randi(length(RTs{i,3,3}),nSamps,1))],[],2);
    RTs{i,4,4} = min([RTs{i,2,4}(randi(length(RTs{i,2,4}),nSamps,1)),RTs{i,3,4}(randi(length(RTs{i,3,4}),nSamps,1))],[],2);
    RTs{i,4,5} = min([RTs{i,2,5}(randi(length(RTs{i,2,5}),nSamps,1)),RTs{i,3,5}(randi(length(RTs{i,3,5}),nSamps,1))],[],2);
    RTs{i,4,6} = min([RTs{i,2,6}(randi(length(RTs{i,2,6}),nSamps,1)),RTs{i,3,6}(randi(length(RTs{i,3,6}),nSamps,1))],[],2);
    RTs{i,4,7} = min([RTs{i,2,7}(randi(length(RTs{i,2,7}),nSamps,1)),RTs{i,3,7}(randi(length(RTs{i,3,7}),nSamps,1))],[],2);
    RTs{i,4,8} = min([RTs{i,2,8}(randi(length(RTs{i,2,8}),nSamps,1)),RTs{i,3,8}(randi(length(RTs{i,3,8}),nSamps,1))],[],2);
    %     RTs{i,4,9} = min([RTs{i,2,9}(randi(length(RTs{i,2,9}),nSamps,1)),RTs{i,3,9}(randi(length(RTs{i,3,9}),nSamps,1))],[],2);
    %     RTs{i,4,10} = min([RTs{i,2,10}(randi(length(RTs{i,2,10}),nSamps,1)),RTs{i,3,10}(randi(length(RTs{i,3,10}),nSamps,1))],[],2);
    %     RTs{i,4,11} = min([RTs{i,2,11}(randi(length(RTs{i,2,11}),nSamps,1)),RTs{i,3,11}(randi(length(RTs{i,3,11}),nSamps,1))],[],2);
    %     RTs{i,4,12} = min([RTs{i,2,12}(randi(length(RTs{i,2,12}),nSamps,1)),RTs{i,3,12}(randi(length(RTs{i,3,12}),nSamps,1))],[],2);
    %     RTs{i,4,13} = min([RTs{i,2,13}(randi(length(RTs{i,2,13}),nSamps,1)),RTs{i,3,13}(randi(length(RTs{i,3,13}),nSamps,1))],[],2);
    
end

% Compute CDFs for each subject, condition, trial type
CDFs = zeros(nSubjs,nConds2,nTasks,numel(prob));
for i = 1:nSubjs
    for j = 1:nTasks
        [CDFs(i,6,j,:),CDFs(i,2,j,:),CDFs(i,3,j,:),CDFs(i,1,j,:)] = ...
            racemodel(RTs{i,2,j},RTs{i,3,j},RTs{i,1,j},...
            prob,'lim',perLim(i,:),'dep',dep,'test',test);
    end
end

% Compute number of RTs for each subject, condition, trial type
nRTs = zeros(nSubjs,nConds,nTasks,numel(prob));
for i = 1:nSubjs
    for j = 1:nConds
        for k = 1:nTasks
            for n = 2:numel(prob)
                nRTs(i,j,k,n) = sum(RTs{i,j,k}<=RTper(i,n));
            end
        end
    end
end

% Compute RMV based on simulation
CDFs(:,5,:,:) = CDFs(:,1,:,:) - CDFs(:,4,:,:);

switch test
    
    case 'ver'
        
        % Compute RMV based on Raab's model
        CDFs(:,7,:,:) = CDFs(:,1,:,:) - CDFs(:,6,:,:);
        
        % Compute Miller's bound
        CDFs(:,8,:,:) = CDFs(:,2,:,:) + CDFs(:,3,:,:); CDFs(CDFs>1) = 1;
        
        % Compute RMV based on Miller's bound
        CDFs(:,9,:,:) = CDFs(:,1,:,:) - CDFs(:,8,:,:);
        
        % Compute empirical benefit based on Grice's bound
        CDFs(:,10,:,:) = CDFs(:,1,:,:) - max(CDFs(:,2:3,:,:),[],2);
        
        % Compute predicted benefit (Raab's model vs. Grice's bound)
        CDFs(:,11,:,:) = CDFs(:,6,:,:) - max(CDFs(:,2:3,:,:),[],2);
        
    case 'hor'
        
        % Compute RMV based on Raab's model
        CDFs(:,7,:,:) = CDFs(:,6,:,:) - CDFs(:,1,:,:);
        
        % Compute Miller's bound
        CDFs(:,8,:,:) = CDFs(:,2,:,:) + CDFs(:,3,:,:);
        
        % Compute RMV based on Miller's bound
        CDFs(:,9,:,:) = CDFs(:,1,:,:) - CDFs(:,8,:,:);
        
        % Compute empirical benefit based on Grice's bound
        CDFs(:,10,:,:) = CDFs(:,1,:,:) - max(CDFs(:,2:3,:,:),[],2);
        
        % Compute predicted benefit (Raab's model vs. Grice's bound)
        CDFs(:,11,:,:) = CDFs(:,6,:,:) - max(CDFs(:,2:3,:,:),[],2);
        
end

%% Compute multisnesory benefits and gain

meth = 'all'; % 'all' 'pos' 'quant'
quant = 3;

% Define race model
race = mean(CDFs(:,6,6:8,:),3);
grice = CDFs(:,2:3,1,:);

% Compute multisensory benefits
switch test
    case 'ver'
        benPred = race - max(grice,[],2);
        benPred = squeeze(trapz(prob,benPred,4));
        benEmp = CDFs(:,1,1,:) - max(grice,[],2);
        benEmp = squeeze(trapz(prob,benEmp,4));
    case 'hor'
        benPred = min(grice,[],2) - race;
        benPred = squeeze(trapz(prob,benPred,4));
        benEmp = min(grice,[],2) - CDFs(:,1,1,:);
        benEmp = squeeze(trapz(prob,benEmp,4));
end

% Compute violation of channel-dependent race model
switch test
    case 'ver'
        rmv = squeeze(CDFs(:,1,1,:) - race);
    case 'hor'
        rmv = squeeze(race - CDFs(:,1,1,:));
end

% Compute modality switch costs
switch test
    case 'ver'
        msc = squeeze(CDFs(:,1:3,2,:) - CDFs(:,1:3,3,:));
    case 'hor'
        msc = squeeze(CDFs(:,1:3,3,:) - CDFs(:,1:3,2,:));
end

% Compute multisensory gain
switch meth
    case 'quant'
        gain = rmv(:,quant);
        cost = msc(:,quant);
    case 'all'
        gain = trapz(prob,rmv,2);
        cost = trapz(prob,msc,3);
end

%% Define plotting parameters

shades = [0,0.25,0.5,0.75];
colors = [0.0000    0.4470    0.7410
    0.8500    0.3250    0.0980
    0.9290    0.6940    0.1250
    0.4940    0.1840    0.5560
    0.4660    0.6740    0.1880
    0.3010    0.7450    0.9330
    0.6350    0.0780    0.1840];

cols1 = colors;
cols2 = linspecer(3);
cols3 = [0,0,0;0.3,0.3,0.3;0.6,0.6,0.6];

age_groups = {'6–9','10–12','13–17','18–40'};
diag_groups = {'NT','ASD'};
sex_groups = {'Female','Male'};
conditions = {'AV','A','V'};
switch_conds = {'V/A\rightarrowAV','V\rightarrowA','A\rightarrowV'};
labels = {'a','b'};

alpha = 0.05;
fnt = 'Helvetica';
fntSz = 10;
mkrSz = 10;
set(0,'DefaultAxesTitleFontWeight','normal')

%% Compare Race Model Approaches

clc;

% Set params
test = 1;
nperm = 1e4;
quant = 2:numel(prob)-1;

if test == 1
    % Channel-independent race model vs. channel-dependent race model
    [tstat,stats,orig] = permtest_t(squeeze(mean(CDFs(:,6,6:8,quant),3)),squeeze(mean(CDFs(:,6,1,quant),3)));
    d = mes(squeeze(mean(CDFs(:,6,6:8,quant),3)),squeeze(mean(CDFs(:,6,1,quant),3)),'hedgesg','isDep',1,'nBoot',nperm);
    d.hedgesg
elseif test == 2
    % Simulated race distribution vs. Raab's model
    [tstat,corx] = permtest_t(squeeze(mean(CDFs(:,6,6:8,quant),3)),squeeze(mean(CDFs(:,4,6:8,quant),3)));
    d = mes(squeeze(mean(CDFs(:,6,6:8,quant),3)),squeeze(mean(CDFs(:,4,6:8,quant),3)),'hedgesg','isDep',1,'nBoot',nperm);
    d.hedgesg
elseif test == 2
    % Simulated race distribution vs. Miller's bound
    [tstat,corx] = permtest_t(squeeze(mean(CDFs(:,8,1,quant),3)),squeeze(mean(CDFs(:,4,6:8,quant),3)));
    d = mes(squeeze(mean(CDFs(:,8,1,quant),3)),squeeze(mean(CDFs(:,4,6:8,quant),3)),'hedgesg','isDep',1,'nBoot',nperm);
    d.hedgesg
end

%% RESULTS SECTION

% 1. Reaction times and multisensory benefits

%% 1.1 Single trial reaction times (LME)
% A linear mixed-effects analysis was used to examine the effect of group,
% age and stimulus condition on RTs. Subject, ISI and preceding modality
% were included as random factors, with slope adjustments for condition.

close all hidden;
clc;

% Convert cell arrays to vectors
RT = cell2mat(RT_cell);
RTprev = cell2mat(RTprev_cell);
ISI = cell2mat(ISI_cell);
type = cell2mat(type_cell);
task = cell2mat(task_cell);
cond = cell2mat(cond_cell);
grp = cell2mat(grp_cell);
subj = cell2mat(subj_cell);
sex = cell2mat(sex_cell);
age = cell2mat(age_cell);
dev = cell2mat(dev_cell);
piq = cell2mat(piq_cell);
viq = cell2mat(viq_cell);
fsiq = cell2mat(fsiq_cell);
ados = cell2mat(ados_cell);

% Get previous modality
prev = task;
prev(prev==1|prev==4|prev==7) = 1;
prev(prev==2|prev==5|prev==8) = 2;
prev(prev==3|prev==6|prev==9) = 3;
prev(isnan(RTprev)) = NaN;

% Round RTs
ISI = round(ISI);

% Convert categorical variables to nominal or ordinal data
subj = categorical(subj,'Ordinal',0);
grp = categorical(grp,[1,2],{'NT','ASD'},'Ordinal',0);
dev = categorical(dev,[1,2,3,4],{'6-9','10-12','13-17','18-40'},'Ordinal',1);
cond = categorical(cond,[1,2,3],{'AV','A','V'},'Ordinal',0);
prev = categorical(prev,[1,2,3],{'AV','A','V'},'Ordinal',0);
type = categorical(type,[1,2],{'Switch','Repeat'},'Ordinal',0);
task = categorical(task,'Ordinal',0);

% Create dataset containing variables
ds = dataset(subj,grp,sex,age,dev,cond,prev,task,type,RT,RTprev,ISI,...
    piq,viq,fsiq,ados);

% Fit linear mixed-effects model by ML criterion
disp('Fitting linear mixed-effects model...')
lm = fitlme(ds,'RT ~ grp + age * cond + (1+cond|subj) + (1+cond|prev) + (1|ISI)');

% Run ANOVA using Satterthwaite approximation
disp('Performing hypothesis tests on fixed effect terms...')
stats = anova(lm);

% Display LME and ANOVA stats
disp('Linear mixed-effects analysis:')
disp(lm.Rsquared)
disp(lm)
disp('ANOVA on fixed effect terms:')
disp(stats)

% Compute residuals stats
r = residuals(lm,'ResidualType','Raw');
pr = residuals(lm,'ResidualType','Pearson');
st = residuals(lm,'ResidualType','Standardized');

figure,
subplot(2,2,1), histfit(r), xlabel('Residuals'), ylabel('Frequency')
title('Histogram of residual with fitted normal density')
subplot(2,2,2), qqplot(r)
% subplot(2,2,3), plotResiduals(lm,'fitted')
subplot(2,2,4), plotResiduals(lm,'lagged')

figure, boxplot([r,pr,st])

% Check normality of residuals
mc_normalitytest(r);

% % Plot intercept and slope CIs and relationship
% figure, plotIntervals(ds,'RT ~ cond','subj')
% [~,~,RE] = randomEffects(lme);
% figure, scatter(RE.Estimate(1:2:end),RE.Estimate(2:2:end))
% [r,p] = corr(RE.Estimate(1:2:end),RE.Estimate(2:2:end))

%% 1.2 Multisensory benefits (LM)
% A linear regression analysis was used to examine the effects of diagnosis
% and age on multisensory benefit.

close all;
clc;

% Set params
ben_type = 'emp'; % 'pred' or 'emp'
stand_vars = false;

% Define variables
subj = 1:nSubjs;
grp = grp_vec;
sex = sex_vec;
age = age_vec;
dev = dev_vec;
switch ben_type
    case 'pred'
        rse = benPred;
    case 'emp'
        rse = benEmp;
end

% Standardize continuous variabiles
if stand_vars
    rse = zscore(rse);
    age = zscore(age);
end

% Convert non-numeric variables to categorical
subj = categorical(subj','Ordinal',0);
grp = categorical(grp,[1,2],{'NT','ASD'},'Ordinal',0);
sex = categorical(sex,[1,2],{'Male','Fem'},'Ordinal',0);
dev = categorical(dev,[1,2,3,4],{'6-9','10-12','13-17','18-40'},'Ordinal',1);

% Create dataset
ds = table(subj,grp,sex,age,dev,rse);

% Fit linear model to dataset
lm = fitlm(ds,'rse ~ grp + age');

% Perform hypothesis tests on fixed effect terms
disp('Performing hypothesis tests on fixed effect terms...')
stats = anova(lm);

% Display LME and ANOVA stats
disp(lm)
disp(stats)

% Define residuals
r = lm.Residuals.Raw;
pr = lm.Residuals.Pearson;
st = lm.Residuals.Standardized;

% Bayes factor analysis
models = {'rse ~ grp + age','rse ~ grp','rse ~ age'};
bf01 = 1./bf.anova(ds,models);
fprintf('\nBayes factor analysis:')
fprintf('\nBF01 (full): %d',bf01(1))
fprintf('\nBF01 (grp): %d',bf01(1)/bf01(3))
fprintf('\nBF01 (age): %d\n',bf01(1)/bf01(2))

% % Bayes factor analysis
% models = {'gain ~ grp*age','gain ~ grp + age:grp','gain ~ age + grp:age'};
% bf01 = 1./bf.anova(ds,models);
% fprintf('\nBayes factor analysis:')
% fprintf('\nBF01 (full): %d',bf01(1))
% fprintf('\nBF01 (grp): %d',bf01(1)/bf01(3))
% fprintf('\nBF01 (age): %d\n\n',bf01(1)/bf01(2))

% Check normality of residuals
disp('Normality test on residuals:')
mc_normalitytest(r);

figure,
subplot(2,2,1), histfit(r), xlabel('Residuals'), ylabel('Frequency')
title('Histogram of residuals with fitted normal density')
subplot(2,2,2), qqplot(r)
subplot(2,2,3), plotResiduals(lm,'fitted')
subplot(2,2,4), plotResiduals(lm,'lagged')

figure, boxplot([r,pr,st])

%% 1.3 RT variability (LM)
% A linear regression analysis was used to examine the effects of diagnosis
% age and condition on RT variability.

close all;
clc;

% Set params
ct_meth = 'mad1';

% Compute RT variability for each subject, condition, trial type
RTvar = zeros(nSubjs,3);
for i = 1:nSubjs
    for j = 1:3
        if strcmp(ct_meth,'std')
            RTvar(i,j) = std(RTs{i,j,1});
        elseif strcmp(ct_meth,'iqr')
            RTvar(i,j) = iqr(RTs{i,j,1});
        elseif strcmp(ct_meth,'mad0')
            RTvar(i,j) = mad(RTs{i,j,1},0);
        elseif strcmp(ct_meth,'mad1')
            RTvar(i,j) = mad(RTs{i,j,1},1);
        elseif strcmp(ct_meth,'range')
            RTvar(i,j) = range(RTs{i,j,1});
        end
    end
end
RTvar = reshape(RTvar,[nSubjs*3,1]);

% Define variables
subj = [1:nSubjs,1:nSubjs,1:nSubjs];
grp = [grp_vec;grp_vec;grp_vec];
sex = [sex_vec;sex_vec;sex_vec];
age = [age_vec;age_vec;age_vec];
dev = [dev_vec;dev_vec;dev_vec];
cond = [1*ones(nSubjs,1);2*ones(nSubjs,1);3*ones(nSubjs,1)];

% Convert non-numeric variables to categorical
subj = categorical(subj','Ordinal',0);
grp = categorical(grp,[1,2],{'NT','ASD'},'Ordinal',0);
sex = categorical(sex,[1,2],{'Male','Fem'},'Ordinal',0);
dev = categorical(dev,[1,2,3,4],{'6-9','10-12','13-17','18-40'},'Ordinal',1);
cond = categorical(cond,[1,2,3],{'AV','A','V'},'Ordinal',0);

% Create dataset
ds = dataset(subj,grp,sex,age,dev,cond,RTvar);

% Fit linear model to dataset
lm = fitlm(ds,'RTvar ~ grp + age + cond');

% Perform hypothesis tests on fixed effect terms
stats = anova(lm);

% Display LME and ANOVA stats
disp(lm)
fprintf('\nANOVA on fixed effect terms:\n')
disp(stats)

% Define residuals
r = lm.Residuals.Raw;
pr = lm.Residuals.Pearson;
st = lm.Residuals.Standardized;

% Check normality of residuals
disp('Normality test on residuals:')
mc_normalitytest(r);

figure,
subplot(2,2,1), histfit(r), xlabel('Residuals'), ylabel('Frequency')
title('Histogram of residuals with fitted normal density')
subplot(2,2,2), qqplot(r)
subplot(2,2,3), plotResiduals(lm,'fitted')
subplot(2,2,4), plotResiduals(lm,'lagged')

figure, boxplot([r,pr,st])

%% Fig. 1 | RT variability and multisensory benefits

close all;
clc;

% Set params
ct_meth = 'med'; % 'mean' 'med' 'mode'
rVal = 'r'; % 'r' 'rho' 'R2'
nperm = 1e4;
age_match = false;
test_norm = false;
fig_label = true;

% Compute RT central tendancy for each subject, condition, trial type
RTct = zeros(nSubjs,3,4);
for i = 1:nSubjs
    for j = 1:3
        for k = 1:nTasks
            if strcmp(ct_meth,'mean')
                RTct(i,j,k) = mean(RTs{i,j,k});
            elseif strcmp(ct_meth,'med')
                RTct(i,j,k) = median(RTs{i,j,k});
            elseif strcmp(ct_meth,'mode')
                RTct(i,j,k) = mode(round(RTs{i,j,k}));
            end
        end
    end
end

% Preallocate memory
idx = cell(2,4);
RTavg = zeros(2,4,3);
err = zeros(2,2,4,3);
ci = zeros(2,2,4,3);

% Compute group RT central tendancy and 95% CIs
for i = 1:2
    for j = 1:4
        idx{i,j} = find(grp_vec==i & dev_vec==j);
        if test_norm
            for k = 1:3
                mc_normalitytest(RTct(idx{i,j},k,1))
            end
        end
        if strcmp(ct_meth,'mean')
            RTavg(i,j,:) = squeeze(mean(RTct(idx{i,j},:,1)));
            ci(:,i,j,:) = bootci(nperm,@mean,squeeze(RTct(idx{i,j},:,1)));
        elseif strcmp(ct_meth,'med')
            RTavg(i,j,:) = squeeze(median(RTct(idx{i,j},:,1)));
            ci(:,i,j,:) = bootci(nperm,@median,squeeze(RTct(idx{i,j},:,1)));
        elseif strcmp(ct_meth,'mode')
            RTavg(i,j,:) = squeeze(mode(round(RTct(idx{i,j},:,1))));
            ci(:,i,j,:) = bootci(nperm,@mode,squeeze(round(RTct(idx{i,j},:,1))));
        end
    end
end

% Compute error
err(1,:,:,:) = RTavg-squeeze(ci(1,:,:,:));
err(2,:,:,:) = squeeze(ci(2,:,:,:))-RTavg;

% Generate figure
axes('parent',figure,'fontname',fnt,'fontsize',fntSz)
set(gcf,'pos',[20,150,800,700])

% Panel A
x1 = 0; x2 = 5;
y1 = 200; y2 = 600;
lnWth = 1.5;
for i = 1:2
    plt = subplot(3,4,i);
    hold on
    for j = 1:3
        plot(1:4,squeeze(RTavg(i,:,j))+1e3,'linewidth',lnWth,'color',cols3(j,:))
    end
    for j = 1:3
        errorbar(1:4,RTavg(i,:,j),squeeze(err(1,i,:,j)),squeeze(err(2,i,:,j)),'color',cols3(j,:),'linewidth',lnWth)
    end
    set(gca,'box','off','tickdir','out','fontname',fnt,'fontsize',fntSz,...
        'xtick',1:4,'xticklabel',age_groups,'ytick',200:100:600)
    title(diag_groups{i},'fontname',fnt,'fontsize',fntSz)
    if i==1
        if fig_label
            text(x1-(x2-x1)*0.5,y2,'a','fontname',fnt,'fontsize',fntSz+2,...
                'fontweight','bold','verticalalign','bottom')
        end
        xlabel('Age (yrs)','fontname',fnt,'fontsize',fntSz)
        ylabel('Median RT (ms)','fontname',fnt,'fontsize',fntSz)
        l = legend(conditions);
        pos = get(l,'pos');
        pos(1) = pos(1)+0.01; pos(2) = pos(2)+0.03;
        set(l,'box','off','fontname',fnt,'fontsize',fntSz-2,'pos',pos)
        pos = get(plt,'pos');
        set(plt,'pos',[pos(1),pos(2)+0.03,pos(3),pos(4)]);
    else
        set(gca,'yticklabel',[])
        pos = get(plt,'pos');
        set(plt,'pos',[pos(1)-0.02,pos(2)+0.03,pos(3),pos(4)]);
    end
    xlim([x1,x2])
    ylim([y1,y2])
    xtickangle(45)
    axis square
end

% Panel B (left)
switch test
    case 'ver'
        x1 = 0; x2 = 1;
        xlab1 = 'Quantile';
    case 'hor'
        x1 = 185; x2 = 525;
        xlab1 = 'RT (ms)';
end
y1 = 0; y2 = 1;

plt = subplot(3,5,4);
hold all
switch test
    case 'ver'
        fill([prob,fliplr(prob)],[squeeze(mean(CDFs(exSub,6,6:8,:),3))',fliplr(squeeze(max(CDFs(exSub,2:3,1,:),[],2))')],[1,0.8,0.7],'EdgeColor','None');
        for i = 1:3
            plot(prob,squeeze(CDFs(exSub,i,1,:)),'linewidth',lnWth,'color',cols3(i,:))
        end
        plot(prob,squeeze(mean(CDFs(exSub,6,6:8,:),3)),':k','linewidth',lnWth)
        set(gca,'tickdir','out','fontname',fnt,'fontsize',fntSz,...
            'xtick',0:0.25:1,'xticklabel',{'0','','0.5','','1'},...
            'ytick',0:0.25:1,'yticklabel',{'0','','0.5','','1'})
    case 'hor'
        fill([squeeze(mean(CDFs(exSub,6,6:8,:),3))',fliplr(squeeze(min(CDFs(exSub,2:3,1,:),[],2))')],[prob,fliplr(prob)],[1,0.8,0.7],'EdgeColor','None');
        for i = 1:3
            plot(squeeze(CDFs(exSub,i,1,:)),prob,'linewidth',lnWth,'color',cols3(i,:))
        end
        plot(squeeze(mean(CDFs(exSub,6,6:8,:),3)),prob,':k','linewidth',lnWth)
        set(gca,'tickdir','out','fontname',fnt,'fontsize',fntSz,...
            'xtick',200:100:500,...
            'ytick',0:0.25:1,'yticklabel',{'0','','0.5','','1'})
end
xlim([x1,x2]), ylim([y1,y2])
title({'Predicted';'Benefit'},'fontname',fnt,'fontsize',fntSz)
if fig_label
    text(x1-(x2-x1)*0.5,y2,'b','fontname',fnt,'fontsize',fntSz+2,...
        'fontweight','bold','verticalalign','bottom')
end
xlabel(xlab1,'fontname',fnt,'fontsize',fntSz)
ylabel({'Cumulative';'probability'},'fontname',fnt,'fontsize',fntSz)
pos = get(plt,'pos');
set(plt,'pos',[pos(1),pos(2)+0.03,pos(3),pos(4)]);
axis square

% Panel B (right)
plt = subplot(3,5,5);
hold all
switch test
    case 'ver'
        fill([prob,fliplr(prob)],[squeeze(CDFs(exSub,1,1,:))',fliplr(squeeze(max(CDFs(exSub,2:3,1,:),[],2))')],[1,0.8,0.7],'EdgeColor','None');
        for i = 1:3
            pl(i) = plot(prob,squeeze(CDFs(exSub,i,1,:)),'linewidth',lnWth,'color',cols3(i,:));
        end
        pl(4) = plot(prob,squeeze(mean(CDFs(exSub,6,6:8,:),3)),':k','linewidth',lnWth);
        set(gca,'tickdir','out','fontname',fnt,'fontsize',fntSz,...
            'xtick',0:0.25:1,'xticklabel',{'0','','0.5','','1'},...
            'ytick',0:0.25:1,'yticklabel',[])
    case 'hor'
        fill([squeeze(CDFs(exSub,1,1,:))',fliplr(squeeze(min(CDFs(exSub,2:3,1,:),[],2))')],[prob,fliplr(prob)],[1,0.8,0.7],'EdgeColor','None');
        for i = 1:3
            pl(i) = plot(squeeze(CDFs(exSub,i,1,:)),prob,'linewidth',lnWth,'color',cols3(i,:));
        end
        pl(4) = plot(squeeze(mean(CDFs(exSub,6,6:8,:),3)),prob,':k','linewidth',lnWth);
        set(gca,'tickdir','out','fontname',fnt,'fontsize',fntSz,...
            'xtick',200:100:500,...
            'ytick',0:0.25:1,'yticklabel',[])
end
xlim([x1,x2]), ylim([y1,y2])
title({'Empirical';'Benefit'},'fontname',fnt,'fontsize',fntSz)
l = legend(pl,[conditions,'Race']);
pos = get(l,'pos');
pos(1) = pos(1)-0.08; pos(2) = pos(2)-0.17;
set(l,'box','off','fontname',fnt,'fontsize',fntSz-2,'pos',pos,...
    'orientation','horizontal')
pos = get(plt,'pos');
set(plt,'pos',[pos(1),pos(2)+0.03,pos(3),pos(4)]);
axis square

% Panel C
switch test
    case 'ver'
        x1 = -0.035; x2 = 0.15;
        y1 = -0.035; y2 = 0.15;
        units = ' (AUC)';
    case 'hor'
        x1 = -55; x2 = 215;
        y1 = -55; y2 = 100;
        units = ' (ms)';
end
ctr = 5;
for i = 1:4
    
    switch age_match
        case true
            [idx1,idx2] = agematch(i,grp_vec,dev_vec,age_vec,piq_vec,false);
        case false
            idx1 = find(grp_vec==1 & dev_vec==i);
            idx2 = find(grp_vec==2 & dev_vec==i);
    end
    
    subplot(3,4,i+4)
    hold all
    scatter(benPred(idx1),benEmp(idx1),20,cols1(1,:),'filled')
    scatter(benPred(idx2),benEmp(idx2),20,cols1(2,:),'filled')
    plot([x1,x2],[0,0],':k')
    plot([0,0],[y1,y2],':k')
    plot([x1,x2],[x1,x2],'--k')
    xlim([x1,x2]), ylim([y1,y2])
    switch test
        case 'ver'
            set(gca,'tickdir','out','fontname',fnt,'fontsize',fntSz,...
                'xtick',-0.05:0.05:0.15,'ytick',-0.05:0.05:0.15)
        case 'hor'
            set(gca,'tickdir','out','fontname',fnt,'fontsize',fntSz,...
                'xtick',0:100:200,'ytick',-50:50:100)
    end
    title([age_groups{i},' yrs'],'fontname',fnt,'fontsize',fntSz)
    if i==1
        xlabel({'Predicted';['benefit',units]},'fontname',fnt,'fontsize',fntSz)
        ylabel({'Empirical';['benefit',units]},'fontname',fnt,'fontsize',fntSz)
        if fig_label
            text(x1-(x2-x1)*0.5,y2,'c','fontname',fnt,'fontsize',fntSz+2,...
                'fontweight','bold','verticalalign','bottom')
        end
        for j = 1:2
            text(x1+x2*0.05,y2-(j*y2*0.13),diag_groups{j},'fontname',fnt,...
                'fontsize',fntSz,'color',cols1(j,:),'FontWeight','Bold')
        end
    else
        set(gca,'yticklabel',[])
    end
    axis square
    ctr = ctr+1;
    
end

ct_meth = 'mad1';

% Compute RT variability for each subject, condition, trial type
RTvar = zeros(nSubjs,3,1);
for i = 1:nSubjs
    for j = 1:3
        for k = 1
            if strcmp(ct_meth,'std')
                RTvar(i,j,k) = std(RTs{i,j,k});
            elseif strcmp(ct_meth,'iqr')
                RTvar(i,j,k) = iqr(RTs{i,j,k});
            elseif strcmp(ct_meth,'mad0')
                RTvar(i,j,k) = mad(RTs{i,j,k},0);
            elseif strcmp(ct_meth,'mad1')
                RTvar(i,j,k) = mad(RTs{i,j,k},1);
            elseif strcmp(ct_meth,'range')
                RTvar(i,j,k) = range(RTs{i,j,k});
            end
        end
    end
end

% Preallocate memory
idx = cell(2,4);
RTavg = zeros(2,4,3);
err = zeros(2,2,4,3);
ci = zeros(2,2,4,3);

ct_meth = 'mean';

% Compute group RT central tendancy and 95% CIs
for i = 1:2
    for j = 1:4
        idx{i,j} = find(grp_vec==i & dev_vec==j);
        if test_norm
            for k = 1:3
                mc_normalitytest(RTvar(idx{i,j},k,1))
            end
        end
        if strcmp(ct_meth,'mean')
            RTavg(i,j,:) = squeeze(mean(RTvar(idx{i,j},:,1)));
            ci(:,i,j,:) = bootci(nperm,@mean,squeeze(RTvar(idx{i,j},:,1)));
        elseif strcmp(ct_meth,'median')
            RTavg(i,j,:) = squeeze(median(RTvar(idx{i,j},:,1)));
            ci(:,i,j,:) = bootci(nperm,@median,squeeze(RTvar(idx{i,j},:,1)));
        elseif strcmp(ct_meth,'mode')
            RTavg(i,j,:) = squeeze(mode(round(RTvar(idx{i,j},:,1))));
            ci(:,i,j,:) = bootci(nperm,@mode,squeeze(round(RTvar(idx{i,j},:,1))));
        end
    end
end

% Compute error
err(1,:,:,:) = RTavg-squeeze(ci(1,:,:,:));
err(2,:,:,:) = squeeze(ci(2,:,:,:))-RTavg;

% Panel D
x1 = 0; x2 = 5;
y1 = 20; y2 = 130;
lnWth = 1.5;
for i = 1:2
    plt = subplot(3,4,i+8);
    hold on
    for j = 1:3
        plot(1:4,squeeze(RTavg(i,:,j))+1e3,'linewidth',lnWth,'color',cols3(j,:))
    end
    for j = 1:3
        errorbar(1:4,RTavg(i,:,j),squeeze(err(1,i,:,j)),squeeze(err(2,i,:,j)),'color',cols3(j,:),'linewidth',lnWth)
    end
    set(gca,'box','off','tickdir','out','fontname',fnt,'fontsize',fntSz,...
        'xtick',1:4,'xticklabel',age_groups)
    title(diag_groups{i},'fontname',fnt,'fontsize',fntSz)
    if i==1
        if fig_label
            text(x1-(x2-x1)*0.5,y2,'d','fontname',fnt,'fontsize',fntSz+2,...
            'fontweight','bold','verticalalign','bottom')
        end
        l = legend(conditions);
        pos = get(l,'pos');
        pos(1) = pos(1)+0.01; pos(2) = pos(2)-0.01;
        set(l,'box','off','fontname',fnt,'fontsize',fntSz-2,'pos',pos)
        ylabel('MAD of RTs (ms)','fontname',fnt,'fontsize',fntSz)
        xlabel('Age (yrs)','fontname',fnt,'fontsize',fntSz)
        pos = get(plt,'pos');
        set(plt,'pos',[pos(1),pos(2)-0.01,pos(3),pos(4)]);
    else
        set(gca,'yticklabel',[])
        pos = get(plt,'pos');
        set(plt,'pos',[pos(1)-0.02,pos(2)-0.01,pos(3),pos(4)]);
    end
    xlim([x1,x2])
    ylim([y1,y2])
    xtickangle(45)
    axis square
end


% Panel E
switch test
    case 'ver'
        x1 = 0; x2 = 300;
        y1 = -0.035; y2 = 0.23;
        units = ' (AUC)';
    case 'hor'
        x1 = 0; x2 = 300;
        y1 = -55; y2 = 300;
        units = ' (ms)';
end
benefits = {'Predicted','Empirical'};
for i = 1:2
    
    switch age_match
        case true
            [idx1,idx2] = agematch(i,grp_vec,dev_vec,sex_vec,age_vec,piq_vec,false);
        case false
            idx1 = find(grp_vec==1);
            idx2 = find(grp_vec==2);
    end
    
    xNT = squeeze(mean(RTvar(idx1,2:3,1),2));
    xASD = squeeze(mean(RTvar(idx2,2:3,1),2));
    if i==1
        yNT = benPred(idx1);
        yASD = benPred(idx2);
    else
        yNT = benEmp(idx1);
        yASD = benEmp(idx2);
    end
    
    plt = subplot(3,4,i+10);
    hold all
    scatter(xNT,yNT,15,cols1(1,:),'filled')
    scatter(xASD,yASD,15,cols1(2,:),'filled')
    xlim([x1,x2])
    ylim([y1,y2])
    switch test
        case 'ver'
            set(gca,'tickdir','out','fontname',fnt,'fontsize',fntSz,...
                'ytick',0:0.05:0.2,'xtick',0:100:300)
        case 'hor'
            set(gca,'tickdir','out','fontname',fnt,'fontsize',fntSz,...
                'xtick',0:100:300,'ytick',0:100:300)
    end
    title(benefits{i},'fontname',fnt,'fontsize',fntSz)
    if i==1
        text(x2+50,y2+y2*0.2,'Variability Rule','fontname',fnt,'fontsize',fntSz,...
            'fontweight','normal','verticalalign','bottom','horizontalalign','center')
        ylabel(['Benefit',units],'fontname',fnt,'fontsize',fntSz)
        xlabel({'MAD of unisensory';'RTs (ms)'},'fontname',fnt,'fontsize',fntSz)
        if fig_label
            text(x1-(x2-x1)*0.4,y2,'e','fontname',fnt,'fontsize',fntSz+2,...
                'fontweight','bold','verticalalign','bottom')
        end
        pos = get(plt,'pos');
        set(plt,'pos',[pos(1)+0.02,pos(2)-0.01,pos(3),pos(4)]);
        for j = 1:2
            text(x1+x2*0.05,y2-(j-1)*(y2-y1)*0.1,diag_groups{j},'fontname',fnt,...
                'fontsize',fntSz,'color',cols1(j,:),'FontWeight','Bold',...
                'verticalalign','top')
        end
    else
        set(gca,'yticklabel',[])
        pos = get(plt,'pos');
        set(plt,'pos',[pos(1),pos(2)-0.01,pos(3),pos(4)]);
    end
    
    switch rVal
        case 'r'
            [r,stats] = permtest_corr(xNT,yNT,'nperm',nperm,...
                'type','pearson','tail','right');
            if stats.p<0.05, mkr = '*'; else, mkr = ''; end
            text(x2,y2,['{\itr} = ',num2str(round(r*1e2)/1e2),mkr],...
                'horizontalalign','right','verticalalign','top','color',colors(1,:))
            [r,stats] = permtest_corr(xASD,yASD,'nperm',nperm,...
                'type','pearson','tail','right');
            if stats.p<0.05, mkr = '*'; else, mkr = ''; end
            text(x2,y2-(y2-y1)*0.1,['{\itr} = ',num2str(round(r*1e2)/1e2),mkr],...
                'horizontalalign','right','verticalalign','top','color',colors(2,:))
        case 'rho'
            [r,stats] = permtest_corr(xNT,yNT,'nperm',nperm,...
                'type','spearman','tail','right');
            if stats.p<0.05, mkr = '*'; else, mkr = ''; end
            text(x2,y2,['{\itr}_s = ',num2str(round(r*1e2)/1e2),mkr],...
                'horizontalalign','right','verticalalign','top','color',colors(1,:))
            [r,stats] = permtest_corr(xASD,yASD,'nperm',nperm,...
                'type','spearman','tail','right');
            if stats.p<0.05, mkr = '*'; else, mkr = ''; end
            text(x2,y2-(y2-y1)*0.1,['{\itr}_s = ',num2str(round(r*1e2)/1e2),mkr],...
                'horizontalalign','right','verticalalign','top','color',colors(2,:))
        case 'R2'
            [r,stats] = permtest_corr(xNT,yNT,'nperm',nperm,...
                'type','pearson','tail','right');
            r, r^2, stats.p
            if stats.p<0.05, mkr = '*'; else, mkr = ''; end
            text(x2,y2,['{\itR}^2 = ',num2str(round(r^2*1e2)/1e2),mkr],...
                'horizontalalign','right','verticalalign','top','color',colors(1,:))
            [r,stats] = permtest_corr(xASD,yASD,'nperm',nperm,...
                'type','pearson','tail','right');
            r, r^2, stats.p
            if stats.p<0.05, mkr = '*'; else, mkr = ''; end
            text(x2,y2-(y2-y1)*0.1,['{\itR}^2 = ',num2str(round(r^2*1e2)/1e2),mkr],...
                'horizontalalign','right','verticalalign','top','color',colors(2,:))
    end
        
    axis square

end

set(gcf,'color','w')

%% RESULTS SECTION

% 2. Testing the race model

%% 2.1 Race model Test

clc;

nperm = 1e4;

quant = 2:20;

dvals  = zeros(2,4,numel(quant));
cis  = zeros(2,4,2,numel(quant));
pvals  = zeros(2,4,numel(quant));

for i = 1:2
    for j = 1:4
        idx = grp_vec==i & dev_vec==j;
        [tstat,stats] = permtest_t(rmv(idx,quant),[],'nperm',nperm,'tail','right');
        d = mes(squeeze(CDFs(idx,1,1,quant)),squeeze(mean(CDFs(idx,6,6:8,quant),3)),'hedgesg','isDep',1,'nBoot',nperm);
        dvals(i,j,:) = d.hedgesg;
        cis(i,j,:,:) = d.hedgesgCi;
        pvals(i,j,:) = stats.p;
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


%% Fig. 2: Race Model Violation

close all;
clc;

% Set parameters
lnWth = 2;
nperm = 1e4;
quant = 2:10;
dev_grps = false;
fig_label = true;

% Generate figure
axes('parent',figure,'fontname',fnt,'fontsize',fntSz)
set(gcf,'pos',[20,50,850,600])

% Panels A & B
x1 = 0; x2 = 1;
switch test
    case 'ver'
        y1 = -0.15; y2 = 0.15;
        units = '';
        xaxis = 'Quantile';
    case 'hor'
        y1 = -20; y2 = 20;
        units = ' (ms)';
        xaxis = 'Percentile';
end
ctr = 1;
for i = 1:2
    
    for j = 1:4
        
        % Index groups
        idx = grp_vec==i & dev_vec==j;
        
        % Compute RMV and CIs
        rmvi = rmv(idx,:);
        ci = bootci(nperm,@mean,rmvi);
        err = zeros(length(prob),2);
        err(:,1) = mean(rmvi)-ci(1,:);
        err(:,2) = ci(2,:)-mean(rmvi);
        
        % Run permutation test with tmax correction
        [~,stats,orig] = permtest_t(rmvi(:,quant),[],'nperm',nperm,'tail','right');
        
        plt = subplot(3,5,ctr);
        hold all
        
        for k = 2
            if k == 1
                hCorr = [0,orig.h,0];
                shade = 0.925;
            else
                hCorr = [0,stats.h,0];
                shade = 0.85;
            end
            try
                sigX = prob(hCorr==1);
                for l = 1:length(sigX)
                    X = [sigX(l)-0.025,sigX(l)-0.025,sigX(l)+0.025,sigX(l)+0.025];
                    Y = [y1,y2,y2,y1];
                    C = ones(1,3)*shade;
                    fill(X,Y,C,'edgecolor','none')
                end
            catch
            end
        end
        [l,p] = boundedline(prob,mean(rmvi,1),err,'alpha','cmap',cols1(i,:));
        outlinebounds(l,p);
        plot(prob,mean(rmvi,1),'color',cols1(i,:),'linewidth',lnWth)
        switch test
            case 'ver'
                set(gca,'layer','top','box','off','xaxislocation','origin',...
                    'fontname',fnt,'fontsize',fntSz,...
                    'tickdir','both','xtick',0:0.25:1,'xticklabel',[],...
                    'ytick',-0.15:0.05:0.15,'yticklabel',{'','-0.1','','0','','0.1',''})
            case 'hor'
                set(gca,'layer','top','box','off','xaxislocation','origin',...
                    'fontname',fnt,'fontsize',fntSz,...
                    'tickdir','both','xtick',0:0.25:1,'xticklabel',[],...
                    'ytick',-20:10:20)
        end
        xlim([x1,x2])
        ylim([y1,y2])
        axis square
        if i==1 && j==3
            title({diag_groups{i};[age_groups{j},' yrs']},'fontname',fnt,'fontsize',fntSz)
        elseif i==2 && j==3
            title(diag_groups{i},'fontname',fnt,'fontsize',fntSz)
        elseif i==1
            title([age_groups{j},' yrs'],'fontname',fnt,'fontsize',fntSz)
        end
        if j==1
            ylabel({'\Delta Cum.';['probability',units]},'fontname',fnt,'fontsize',fntSz)
            if fig_label
                text(x1-(x2-x1)*0.5,y2,labels{i},'fontname',fnt,...
                    'fontsize',fntSz+2,'fontweight','bold',...
                    'verticalalign','bottom')
            end
        else
            set(gca,'yticklabel',[])
        end
        
        if i==2
            pos = get(plt,'pos');
            set(plt,'pos',[pos(1),pos(2)+0.025,pos(3),pos(4)])
        end
        
        % Plot x-axis
        if ctr==6
            pos = get(gca,'pos');
            axes('parent',gcf,'pos',[pos(1),pos(2)-0.075,pos(3),pos(4)]);
            hold on
            set(gca,'layer','top','box','off','ycolor','none','color','none')
            set(gca,'tickdir','both','xtick',0:0.25:1,'xticklabel',{'0','','0.5','','1'})
            xlabel(xaxis,'fontname',fnt,'fontsize',fntSz)
            axis square
        end
        
        ctr = ctr+1;
        
    end
    
    % Plot all ages
    plt = subplot(3,5,ctr);
    hold all
    for j = 1:4
        idx = grp_vec==i & dev_vec==j;
        plot(prob,mean(rmv(idx,:)),'color',[shades(j),shades(j),0],'linewidth',lnWth)
    end
    set(gca,'layer','top','box','off','tickdir','both','xaxislocation','origin',...
        'fontname',fnt,'fontsize',fntSz,...
        'xtick',0:0.25:1,'xticklabel',[],'ytick',-0.15:0.05:0.15,'yticklabel',[])
    if i==1
        title('All Ages','fontname',fnt,'fontsize',fntSz)
    else
        pos = get(plt,'pos');
        set(plt,'pos',[pos(1),pos(2)+0.025,pos(3),pos(4)])
        l = legend(age_groups);
        pos = get(l,'pos');
        set(l,'box','off','orientation','vertical',...
            'fontname',fnt,'fontsize',fntSz-2,...
            'pos',[pos(1)+0.02,pos(2)+0.1075,pos(3),pos(4)])
    end
    xlim([x1,x2])
    ylim([y1,y2])
    axis square
    
    ctr = ctr+1;
    
end

% Compute MSE and R values
mseVals = zeros(nSubjs,nSubjs);
for i = 1:nSubjs
    for j = 1:nSubjs
        mseVals(i,j) = sqrt(mean((rmv(i,2:end-1)-rmv(j,2:end-1)).^2));
    end
end
rVals = permtest_corr(rmv(:,2:end-1)',[],'type','rankit');

% Compute averages across different age groups
if dev_grps==true
    
    dim = length(age_groups);
    
    mseMat = zeros(dim,dim);
    rMat = NaN(dim,dim);
    for i = 1:4
        for j = 1:4
            tmp1 = [];
            tmp2 = [];
            for k = find(grp_vec==1 & dev_vec==i)
                idx = grp_vec==2 & dev_vec==j;
                tmp1 = [tmp1;mean(mseVals(k,idx))];
                tmp2 = [tmp2;mean(rVals(k,idx))];
            end
            mseMat(i,j) = mean(tmp1);
            rMat(i,j) = mean(tmp2);
        end
    end
    
elseif dev_grps==false
    
    inc = 2.5;
    ran = 6:inc:18.5;
    dim = numel(ran);
    
    x = 0;
    mseMat = NaN(dim,dim);
    rMat = NaN(dim,dim);
    for i = ran
        x = x+1;
        y = 0;
        for j = ran
            y = y+1;
            tmp1 = [];
            tmp2 = [];
            for k = find(grp_vec==2 & age_vec>=i & age_vec<i+inc)
                idx = grp_vec==1 & age_vec>=j & age_vec<j+inc;
                tmp1 = [tmp1;mean(mseVals(k,idx))];
                tmp2 = [tmp2;mean(rVals(k,idx))];
            end
            mseMat(x,y) = mean(tmp1);
            rMat(x,y) = mean(tmp2);
        end
    end
    
end

[~,imse] = min(mseMat,[],2);
[~,ir] = max(rMat,[],2);

% Panel C (right)
plt = subplot(3,5,14:15);
hold all
imagesc(rMat)
plot(0.5:dim+0.5,0.5:dim+0.5,'--k')
plot(ir,1:dim,'-r','linewidth',2)
if dev_grps==true
    set(gca,'box','on','layer','top','ydir','normal','fontname',fnt,'fontsize',fntSz,...
        'xtick',1:5,'xticklabel',age_groups,...
        'ytick',1:5,'yticklabel',age_groups)
elseif dev_grps==false
    set(gca,'box','on','layer','top','ydir','normal','fontname',fnt,...
        'fontsize',fntSz,'xtick',1:6,'ytick',1:6,...
        'xticklabel',{'6–8.5','11–13.5','15–18.5','18.5–21'},...
        'yticklabel',{'6–8.5','8.5–11','11–13.5','13.5–15','15–18.5','18.5–21'})
%         'yticklabel',{'6–9','9–12','12–15','15–18','18–21'})
%         'xticklabel',{'9–12','15–18','21–24','27–30','33–36'},...
%         'yticklabel',{'9–12','15–18','21–24','27–30','33–36'})
end
xtickangle(45)
title({'Correlation';'Coefficient (RIN)'},'fontname',fnt,'fontsize',fntSz)
xlim([0.5,dim+0.5])
ylim([0.5,dim+0.5])
axis square
c = colorbar;
title(c,'{\itr}_{RIN}','fontname',fnt,'fontsize',fntSz)
colormap bone
pos = get(plt,'pos');
set(plt,'pos',[pos(1),pos(2)+0.01,pos(3),pos(4)])

% Panel C (left)
plt = subplot(3,5,12:13);
hold all
imagesc(mseMat)
plot(0.5:dim+0.5,0.5:dim+0.5,'--k')
plot(imse,1:dim,'-r','linewidth',2)
if dev_grps==true
    set(gca,'box','on','layer','top','ydir','normal','fontname',fnt,'fontsize',fntSz,...
        'xtick',1:5,'xticklabel',age_groups,...
        'ytick',1:5,'yticklabel',age_groups)
elseif dev_grps==false
    set(gca,'box','on','layer','top','ydir','normal','fontname',fnt,...
        'fontsize',fntSz,'xtick',1:6,'ytick',1:6,...
        'xticklabel',{'6–8.5','8.5–11','11–13.5','13.5–15','15–18.5','18.5–21'},...
        'yticklabel',{'6–8.5','8.5–11','11–13.5','13.5–15','15–18.5','18.5–21'})
%         'xticklabel',{'9–12','15–18','21–24','27–30','33–36'},...
%         'yticklabel',{'9–12','15–18','21–24','27–30','33–36'})
end
xtickangle(45)
ylabel('ASD age (yrs)','fontname',fnt,'fontsize',fntSz)
xlabel('NT age (yrs)','fontname',fnt,'fontsize',fntSz)
title({'Root Mean';'Squared Error'},'fontname',fnt,'fontsize',fntSz)
if fig_label
    text(0-dim*0.45,dim+0.5,'c','fontname',fnt,'fontsize',fntSz+2,...
        'fontweight','bold','verticalalign','bottom')
end
xlim([0.5,dim+0.5])
ylim([0.5,dim+0.5])
axis square
c = colorbar;
title(c,'RMSE','fontname',fnt,'fontsize',fntSz)
colormap bone
pos = get(plt,'pos');
set(plt,'pos',[pos(1)+0.02,pos(2)+0.01,pos(3),pos(4)])


set(gcf,'color','w')

rLwr = tril(rMat); rLwr = rLwr(:);
rUpr = triu(rMat); rUpr = rUpr(:);

[tstat,~,orig,stats] = permtest_t(rLwr,rUpr,'nperm',nperm);
d = mes(rLwr,rUpr,'hedgesg','isDep',1,'nBoot',nperm);
fprintf('Rlwr = %0.2f ± %0.2f, ',mean(rLwr),std(rLwr))
fprintf('Rupr = %0.2f ± %0.2f, ',mean(rUpr),std(rUpr))
fprintf('t(%d) = %0.2f, ',stats.df,tstat)
fprintf('p = %0.3f, ',orig.p)
fprintf('g = %0.2f, ',d.hedgesg)
fprintf('95CI [%0.2f, %0.2f], ',d.hedgesgCi)
bf01 = 1/bf.ttest(rLwr,rUpr);
fprintf('BF = %0.2f\n',bf01)

mseLwr = tril(mseMat); mseLwr = mseLwr(:);
mseUpr = triu(mseMat); mseUpr = mseUpr(:);

[tstat,~,orig,stats] = permtest_t(mseLwr,mseUpr,'nperm',nperm);
d = mes(mseLwr,mseUpr,'hedgesg','isDep',1,'nBoot',nperm);
fprintf('MSElwr = %0.2f ± %0.2f, ',mean(mseLwr),std(mseLwr))
fprintf('MSEupr = %0.2f ± %0.2f, ',mean(mseUpr),std(mseUpr))
fprintf('t(%d) = %0.2f, ',stats.df,tstat)
fprintf('p = %0.3f, ',orig.p)
fprintf('g = %0.2f, ',d.hedgesg)
fprintf('95CI [%0.2f, %0.2f], ',d.hedgesgCi)
bf01 = 1/bf.ttest(mseLwr,mseUpr);
fprintf('BF = %0.2f\n',bf01)

%% RESULTS SECTION

% 3. Delayed multisensory development in autism

%% 3.1 Linear regression analysis

close all;
clc;

% Define variables
subj = 1:nSubjs;
grp = grp_vec;
sex = sex_vec;
age = age_vec;
dev = dev_vec;

% Convert non-numeric variables to categorical
subj = categorical(subj','Ordinal',0);
grp = categorical(grp,[1,2],{'TD','ASD'},'Ordinal',0);
sex = categorical(sex,[1,2],{'Male','Fem'},'Ordinal',0);
dev = categorical(dev,[1,2,3,4],{'6-9','10-12','13-17','18-40'},'Ordinal',1);

% Create dataset
ds = table(subj,grp,sex,age,dev,gain);

% Fit linear model to dataset
disp('Fitting linear model...')
lm = fitlm(ds,'gain ~ grp * age');

% Perform hypothesis tests on fixed effect terms
disp('Performing hypothesis tests on fixed effect terms...')
stats = anova(lm);

% Display LME and ANOVA stats
disp(lm)
disp(stats)

% Define residuals
r = lm.Residuals.Raw;
pr = lm.Residuals.Pearson;
st = lm.Residuals.Standardized;

% Bayes factor analysis
models = {'gain ~ grp*age','gain ~ age + grp:age','gain ~ grp + age:grp','gain ~ grp + age'};
bf01 = 1./bf.anova(ds,models);
fprintf('\nBayes factor analysis:')
fprintf('\nBF01 (full): %d',bf01(1))
fprintf('\nBF01 (grp): %d',bf01(1)/bf01(2))
fprintf('\nBF01 (age): %d',bf01(1)/bf01(3))
fprintf('\nBF01 (grp:age): %d\n\n',bf01(1)/bf01(4))

% Check normality of residuals
disp('Normality test on residuals:')
mc_normalitytest(r);

figure,
subplot(2,2,1), histfit(r), xlabel('Residuals'), ylabel('Frequency')
title('Histogram of residuals with fitted normal density')
subplot(2,2,2), qqplot(r)
subplot(2,2,3), plotResiduals(lm,'fitted')
subplot(2,2,4), plotResiduals(lm,'lagged')

figure, boxplot([r,pr,st])

%% 3.2 Regression analysis

clc;

% Set params
nperm = 1e4;
fdr_method = 'bh'; % 'bh' or 'bky' or 'none'

% Initialize variables
ctr = 1;
p = zeros(1,4);
h = zeros(1,4);

for i = 1:2
    idx = grp_vec==i & dev_vec~=4;
    [r,stats] = permtest_corr(gain(idx),age_vec(idx),'nperm',nperm);
    bf01 = 1/bf.corr(gain(idx),age_vec(idx));
    fprintf('r = %1.3f, R2 = %1.3f, p = %1.4f, BF01 = %1.13f\n',r,r^2,stats.p,bf01)
    p(ctr) = stats.p;
    idx = grp_vec==i & dev_vec==4;
    [r,stats] = permtest_corr(gain(idx),age_vec(idx),'nperm',nperm);
    bf01 = 1/bf.corr(gain(idx),age_vec(idx));
    fprintf('r = %1.3f, R2 = %1.3f, p = %1.5f, BF01 = %1.13f\n',r,r^2,stats.p,bf01)
    p(ctr+1) = stats.p;
    ctr = ctr+2;
end

% Apply FDR correction on p-values
switch fdr_method
    case 'bky'
        hx = fdr_bky(p,alpha);
    case 'bh'
        hx = fdr_bh(p,alpha,'pdep');
    case 'none'
        hx = h;
end
disp(['h: ',num2str(hx)])

%% Fig. 3: Group differences & recovery

close all;
clc;

% Set parameters
nperm = 1e4;
lnWth = 2;
age_match = true;
fig_label = true;
fdr_method = 'bh'; % 'bh' or 'bky' or 'none'
rVal = 'R2'; % 'r' or 'rho' or 'R2'

switch test
    case 'ver'
        quant = 2:length(prob)-1;
        units = ' (AUC)';
        xlab1 = 'Quantile';
        xlab2 = 'Quantile';
    case 'hor'
        quant = 1:length(prob)-1;
        units = ' (ms)';
        xlab1 = 'RT (ms)';
        xlab2 = 'Percentile';
end

% Preallocate memory
pVals = zeros(1,4);
gainPos = zeros(nSubjs,1);
gainNeg = zeros(nSubjs,1);

% Define negative and positive RMV
for i = 1:nSubjs
    gainPos(i) = getauc(prob,rmv(i,:),'pos');
    gainNeg(i) = getauc(prob,rmv(i,:),'neg');
end

% Generate figure
axes('parent',figure,'fontname',fnt,'fontsize',fntSz)
set(gcf,'pos',[20,20,850,800])

% Panel A
plt = subplot(4,5,1);
hold all
y1 = 0; y2 = 1;
switch test
    case 'ver'
        x1 = 0; x2 = 1;
        fill([prob,fliplr(prob)],...
            [squeeze(CDFs(exSub,1,1,:))',fliplr(squeeze(mean(CDFs(exSub,6,6:8,:),3))')],...
            [1,0.8,0.7],'EdgeColor','None');
        for i = 1:3
            pl(i) = plot(prob,squeeze(CDFs(exSub,i,1,:)),'-',...
                'linewidth',1.5,'color',cols3(i,:));
        end
        pl(4) = plot(prob,squeeze(mean(CDFs(exSub,6,6:8,:),3)),':k','linewidth',1.5);
        set(gca,'tickdir','out','xtick',0:0.25:1,'xticklabel',{'0','','0.5','','1'},...
            'ytick',0:0.25:1,'yticklabel',{'0','','0.5','','1'})
    case 'hor'
        x1 = 185; x2 = 525;
        fill([squeeze(CDFs(exSub,1,1,:))',fliplr(squeeze(mean(CDFs(exSub,6,6:8,:),3))')],...
            [prob,fliplr(prob)],[1,0.8,0.7],'EdgeColor','None');
        for i = 1:3
            pl(i) = plot(squeeze(CDFs(exSub,i,1,:)),prob,'-','linewidth',1.5,'color',cols3(i,:));
        end
        pl(4) = plot(squeeze(mean(CDFs(exSub,6,6:8,:),3)),prob,':k','linewidth',1.5);
        set(gca,'tickdir','out','xtick',200:100:500,...
            'ytick',0:0.25:1,'yticklabel',{'0','','0.5','','1'})
end
xlim([x1,x2])
ylim([y1,y2])
title('Multisensory Gain','fontname',fnt,'fontsize',fntSz,'fontweight','normal')
xlabel(xlab1,'fontname',fnt,'fontsize',fntSz)
ylabel({'Cumulative';'probability'},'fontname',fnt,'fontsize',fntSz)
if fig_label
    text(x1-(x2-x1)*0.5,y2,'a','fontname',fnt,'fontsize',fntSz+2,...
        'fontweight','bold','verticalalign','bottom')
end
l = legend(pl,'AV','A','V','Race');
pos = get(l,'pos');
set(l,'box','off','fontname',fnt,'fontsize',fntSz-2,...
    'pos',[pos(1)-0.02,pos(2)-0.035,pos(3),pos(4)])
pos = get(plt,'pos');
set(plt,'pos',[pos(1)-0.05,pos(2)+0.01,pos(3),pos(4)])
axis square

% Panel B
plt = subplot(4,5,6);
hold all
idx = 1:nSubjs;
switch test
    case 'ver'
        x1 = 0; x2 = 0.15;
        y1 = 0; y2 = 0.1;
    case 'hor'
        x1 = 0; x2 = 180;
        y1 = 0; y2 = 35;
end
scatter(abs(gainNeg(idx)),gainPos(idx),10,'filled',...
    'markerfacecolor',cols1(4,:))
set(gca,'tickdir','out')
xlim([x1,x2])
ylim([y1,y2])
xlabel('AUC < 0','fontname',fnt,'fontsize',fntSz)
ylabel('AUC > 0','fontname',fnt,'fontsize',fntSz)
if fig_label
    text(x1-(x2-x1)*0.5,y2,'b','fontname',fnt,'fontsize',fntSz+2,...
        'verticalalign','bottom','fontweight','bold')
end
[r,stats] = permtest_corr(abs(gainNeg(idx)),gainPos(idx),'type','spearman');
if stats.p<0.05, mkr = '*'; else mkr = ''; end
text(x2,y2*0.9,['{\itr}_s = ',num2str(round(r*1e2)/1e2),mkr],...
    'fontname',fnt,'fontsize',fntSz,'horizontalalign','right')
pos = get(plt,'pos');
set(plt,'pos',[pos(1)-0.05,pos(2),pos(3),pos(4)])
axis square

% Initialize variables
ctr = 2;
h = zeros(1,4);
p = zeros(1,4);

% Loop though developmental groups
for i = 1:4
    
    switch age_match
        case true
            if i==4
                [idx1,idx2] = agematch(i,grp_vec,dev_vec,sex_vec,age_vec,[],false);
            else
                [idx1,idx2] = agematch(i,grp_vec,dev_vec,sex_vec,age_vec,piq_vec,false);
            end
        case false
            idx1 = find(grp_vec==1 & dev_vec==i);
            idx2 = find(grp_vec==2 & dev_vec==i);
    end
    
    femaleNTs = sum(sex_vec(idx1)==2);
    disp(['NTs: ',num2str(femaleNTs)])
    femaleASDs = sum(sex_vec(idx2)==2);
    disp(['ASDs: ',num2str(femaleASDs)])
    
    % Compute RMV
    rmvNT = rmv(idx1,:);
    gainNT = gain(idx1);
    rmvASD = rmv(idx2,:);
    gainASD = gain(idx2);
    
    % Compute CIs
    ci = bootci(nperm,@mean,rmvNT);
    errNT = zeros(length(prob),2);
    errNT(:,1) = mean(rmvNT)-ci(1,:);
    errNT(:,2) = ci(2,:)-mean(rmvNT);
    clear ci
    ci = bootci(nperm,@mean,rmvASD);
    errASD = zeros(length(prob),2);
    errASD(:,1) = mean(rmvASD)-ci(1,:);
    errASD(:,2) = ci(2,:)-mean(rmvASD);
    clear ci
    
    % Compare variance before running permutation tests
    [~,stats] = permtest_f2(rmvNT(:,quant),rmvASD(:,quant));
    
    % Run permutation tests at each quantile
    if mean(stats.h)<=0.5
        [~,stats,orig] = permtest_t2(rmvNT(:,quant),rmvASD(:,quant),...
            'nperm',nperm);
    elseif mean(stats.h)>0.5
        [~,stats,orig] = permtest_t2(rmvNT(:,quant),rmvASD(:,quant),...
            'nperm',nperm,'varx','unequal');
        disp('Unequal variance')
    end
    
    % Panel C
    switch test
        case 'ver'
            y1 = -0.15; y2 = 0.15;
        case 'hor'
            y1 = -20; y2 = 20;
    end
    x1 = 0; x2 = 1;
    plt = subplot(4,5,ctr);
    hold all
    plot(prob,mean(rmvNT)+1e3,'color',cols1(1,:),'linewidth',lnWth)
    plot(prob,mean(rmvASD)+1e3,'color',cols1(2,:),'linewidth',lnWth)
    for k = 2
        if k == 1
            hCorr = [orig.h,0];
            switch test
                case 'ver'
                    hCorr = [0,hCorr];
            end
            shade = 0.925;
        else
            hCorr = [stats.h,0];
            switch test
                case 'ver'
                    hCorr = [0,hCorr];
            end
            shade = 0.85;
        end
        try
            sigX = prob(hCorr==1);
            for l = 1:length(sigX)
                X = [sigX(l)-0.025,sigX(l)-0.025,sigX(l)+0.025,sigX(l)+0.025];
                Y = [y1,y2,y2,y1];
                C = ones(1,3)*shade;
                fill(X,Y,C,'edgecolor','none')
            end
        catch
        end
    end
    [hl,hp] = boundedline(prob,mean(rmvNT),errNT,prob,mean(rmvASD),errASD,...
        'alpha','cmap',cols1(1:2,:)); outlinebounds(hl,hp);
    plot(prob,mean(rmvNT),'color',cols1(1,:),'linewidth',lnWth)
    plot(prob,mean(rmvASD),'color',cols1(2,:),'linewidth',lnWth)
    switch test
        case 'ver'
            set(gca,'layer','top','box','off','xaxislocation','origin',...
                'tickdir','both','xtick',0:0.25:1,'xticklabel',[],...
                'ytick',-0.15:0.05:0.15,'yticklabel',{'','-0.1','','0','','0.1',''})
        case 'hor'
            set(gca,'layer','top','box','off','xaxislocation','origin',...
                'tickdir','both','xtick',0:0.25:1,'xticklabel',[],...
                'ytick',-20:10:20,'yticklabel',{'','-10','0','10',''})
    end
    xlim([x1,x2])
    ylim([y1,y2])
    axis square
    title([age_groups{i},' yrs'],'fontname',fnt,'fontsize',fntSz,'fontweight','normal')
    if i==1
        if fig_label
            text(x1-(x2-x1)*0.5,y2,'c','fontname',fnt,'fontsize',fntSz+2,...
                'fontweight','bold','verticalalign','bottom')
        end
        switch test
            case 'ver'
                ylabel({'\Delta Cum.';'probability'},'fontname',fnt,'fontsize',fntSz)
            case 'hor'
                ylabel({'\Delta RT (ms)'},'fontname',fnt,'fontsize',fntSz)
        end
    elseif intersect(i,2:4)
        set(gca,'yticklabel',[])
    end
    if i==1
        l = legend(diag_groups);
        pos = get(l,'pos');
        set(l,'box','off','fontname',fnt,'fontsize',fntSz-2,...
            'pos',[pos(1)+0.05,pos(2)+0.03,pos(3),pos(4)])
    elseif i==4
        pos = get(gca,'pos');
        axes('parent',gcf,'pos',[pos(1)+0.03,pos(2)-0.01,pos(3),pos(4)]);
        hold on
        set(gca,'layer','top','box','off','ycolor','none','color','none')
        set(gca,'tickdir','both','xtick',0:0.25:1,'xticklabel',{'0','','0.5','','1'})
        xlabel(xlab2,'fontname',fnt,'fontsize',fntSz)
        axis square
    end
    pos = get(plt,'pos');
    set(plt,'pos',[pos(1)+0.03,pos(2)+0.01,pos(3),pos(4)])
    
    % Compare variance before running permutation test
    [~,~,orig] = permtest_f2(gainNT,gainASD,nperm);
    
    % Run permutation test
    if orig.h==0
        [tstat,~,orig,stats] = permtest_t2(gainNT,gainASD,'nperm',nperm);
        d = mes(gainNT,gainASD,'hedgesg','isDep',0,'nBoot',nperm);
        fprintf('AgeNT: %1.1f (%1.1f), ',mean(age_vec(idx1)),std(age_vec(idx1)));
        fprintf('AgeASD: %1.1f (%1.1f), ',mean(age_vec(idx2)),std(age_vec(idx2)));
        fprintf('RMVnt = %0.2f ± %0.2f, ',mean(gainNT),std(gainNT))
        fprintf('RMVasd = %0.2f ± %0.2f, ',mean(gainASD),std(gainASD))
        fprintf('t(%d) = %0.2f, ',stats.df,tstat)
        fprintf('p = %0.3f, ',orig.p)
        fprintf('g = %0.2f, ',d.hedgesg)
        fprintf('95CI [%0.2f, %0.2f], ',d.hedgesgCi)
        bf01 = 1/bf.ttest2(gainNT,gainASD,'scale',1);
        fprintf('BF = %0.2f\n',bf01)
        d = d.hedgesg;
    elseif orig.h==1
        [tstat,~,orig,stats] = permtest_t2(gainNT,gainASD,'nperm',nperm,'varx','unequal');
        d = mes(gainNT,gainASD,'glassdelta','isDep',0,'nBoot',nperm);
        n = numel(idx1);
        d.glassdelta = d.glassdelta*gamma((n-1)/2)/(sqrt((n-1)/2)*gamma((n-2)/2));
        fprintf('*AgeNT: %1.1f (%1.1f)',mean(age_vec(idx1)),std(age_vec(idx1)));
        fprintf('AgeASD: %1.1f (%1.1f)',mean(age_vec(idx2)),std(age_vec(idx2)));
        fprintf('RMVnt = %0.2f ± %0.2f, ',mean(gainNT),std(gainNT))
        fprintf('RMVasd = %0.2f ± %0.2f, ',mean(gainASD),std(gainASD))
        fprintf('t(%0.2f) = %0.2f, ',stats.df,tstat)
        fprintf('p = %0.3f, ',orig.p)
        fprintf('delta = %0.2f, ',d.glassdelta)
        fprintf('95CI [%0.2f, %0.2f]\n',d.glassdeltaCi)
        bf01 = 1/bf.ttest2(gainNT,gainASD,'scale',1);
        fprintf('BF = %0.2f\n',bf01)
        d = d.glassdelta;
    end
    
    % Panel D
    switch test
        case 'ver'
            y1 = -0.13; y2 = 0.13;
        case 'hor'
            y1 = -180; y2 = 50;
    end
    x1 = 0.5; x2 = 2.5;
    plt = subplot(4,5,ctr+5);
    hold on
    plot([x1,x2],[0,0],':k')
    data = {gainNT,gainASD};
    boxdotplot([1,2],data,cols1(1:2,:),'iqr',nperm,[0.5,0.5,0.5],0.5,0.5,15,0.2,2);
    if orig.p<0.05
        plot([1.2,1.8],[y2,y2]-0.03,'k')
        plot([1.2,1.2],[y2-y2*0.1,y2]-0.03,'k')
        plot([1.8,1.8],[y2-y2*0.1,y2]-0.03,'k')
    end
    if orig.p<0.001
        text(1.5,y2-0.06,'***','horizontalalign','center',...
            'verticalalign','bottom','fontname',fnt,'fontsize',fntSz+8)
    elseif orig.p<0.01
        text(1.5,y2-0.06,'**','horizontalalign','center',...
            'verticalalign','bottom','fontname',fnt,'fontsize',fntSz+8)
    elseif orig.p<0.05
        text(1.5,y2-0.06,'*','horizontalalign','center',...
            'verticalalign','bottom','fontname',fnt,'fontsize',fntSz+8)
    end
    switch test
        case 'ver'
            set(gca,'box','off','tickdir','out','xtick',1:2,...
                'xticklabels',diag_groups,'ytick',-0.1:0.1:0.1)
        case 'hor'
            set(gca,'box','off','tickdir','out','xtick',1:2,...
                'xticklabels',diag_groups,'ytick',-200:50:50)
    end
    ylim([y1,y2])
    xlim([x1,x2])
    axis square
    if i==1
        if fig_label
            text(x1-(x2-x1)*0.5,y2,'d','fontname',fnt,'fontsize',fntSz+2,...
                'fontweight','bold','verticalalign','bottom')
        end
        ylabel({'Multisensory';['gain',units]},'fontname',fnt,'fontsize',fntSz,...
            'horizontalalign','center')
    end
    if intersect(i,2:4)
        set(gca,'yticklabel',[])
    end
    pos = get(plt,'pos');
    set(plt,'pos',[pos(1)+0.03,pos(2),pos(3),pos(4)])
    
    h(i) = orig.h;
    p(i) = orig.p;
    
    ctr = ctr+1;
    
end

% Apply FDR correction on p-values
switch fdr_method
    case 'bky'
        hx = fdr_bky(p,alpha);
    case 'bh'
        hx = fdr_bh(p,alpha,'pdep');
    case 'none'
        hx = h;
end
disp(['h: ',num2str(hx)])

% Panel E
switch test
    case 'ver'
        y1 = -0.13; y2 = 0.1;
    case 'hor'
        y1 = -200; y2 = 50;
end
x1 = 5; x2 = 40;
plt = subplot(4,5,11:12); hold all
for i = 1:2
    idx = find(grp_vec==i);
    scatter(age_vec(idx),gain(idx),20,'filled',...
        'markerfacecolor',cols1(i,:))
    switch test
        case 'ver'
            set(gca,'box','off','tickdir','out','ytick',-0.1:0.1:0.1)
        case 'hor'
            set(gca,'box','off','tickdir','out','ytick',-200:50:50)
    end
    switch rVal
        case 'r'
            [r,stats] = permtest_corr(age_vec(idx),gain(idx),...
                'nperm',nperm,'type','pearson');
            if stats.p<0.05, mkr = '*'; else mkr = ''; end
            text(x2,y1+((3-i)*(y2-y1)*0.1),['{\itr} = ',num2str(round(r*1e2)/1e2),mkr],...
                'horizontalalign','right','color',colors(i,:))
        case 'rho'
            [r,stats] = permtest_corr(age_vec(idx),gain(idx),...
                'nperm',nperm,'type','spearman');
            if stats.p<0.05, mkr = '*'; else mkr = ''; end
            text(x2,y1+((3-i)*(y2-y1)*0.11),['{\itr}_s = ',num2str(round(r*1e2)/1e2),mkr],...
                'horizontalalign','right','color',colors(i,:))
        case 'R2'
            [r,stats] = permtest_corr(age_vec(idx),gain(idx),...
                'nperm',nperm,'type','pearson');
            if stats.p<0.05, mkr = '*'; else mkr = ''; end
            text(x2,y1+((3-i)*(y2-y1)*0.1),['{\itR}^2 = ',num2str(round(r^2*1e2)/1e2),mkr],...
                'horizontalalign','right','color',colors(i,:))
    end
    if i==2
        xlabel('Age (yrs)','fontname',fnt,'fontsize',fntSz)
        ylabel({'Multisensory';['gain',units]},'fontname',fnt,'fontsize',fntSz)
        if fig_label
            text(x1-(x2-x1)*0.25,y2,'e','fontname',fnt,'fontsize',fntSz+2,...
                'fontweight','bold','verticalalign','bottom')
        end
        for j = 1:2
            text(x1+x2*0.02,y2-(j*y2*0.2),diag_groups{j},'fontname',fnt,...
                'fontsize',fntSz,'color',cols1(j,:),'FontWeight','Bold')
        end
    end
    xlim([x1,x2])
    ylim([y1,y2])
    pos = get(plt,'pos');
    set(plt,'pos',[pos(1)-0.025,pos(2)-0.01,pos(3),pos(4)+0.01])
end

% Panel F

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

% Preallocate memory
rmvAvg = zeros(numel(ageRng),2);
errNT = zeros(numel(ageRng),2);
errASD = zeros(numel(ageRng),2);
p = zeros(1,numel(ageRng));
h = zeros(1,numel(ageRng));
bf01 = zeros(1,numel(ageRng));
hedgesg = zeros(1,numel(ageRng));

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
    
    disp(['N: ',num2str(numel(idx1)),', ',num2str(numel(idx2))])
    
    % Compute RMV and CIs
    rmvAvg(i,1) = squeeze(mean(gainNT));
    ci = bootci(nperm,@mean,gainNT);
    errNT(i,1) = rmvAvg(i,1)-ci(1);
    errNT(i,2) = ci(2)-rmvAvg(i,1); clear ci
    rmvAvg(i,2) = squeeze(mean(gainASD));
    ci = bootci(nperm,@mean,gainASD);
    errASD(i,1) = rmvAvg(i,2)-ci(1);
    errASD(i,2) = ci(2)-rmvAvg(i,2); clear ci
    
    % Compare variance before running permutation tests
    [~,~,orig] = permtest_f2(gainNT,gainASD,nperm);
    if mean(orig.h)>0.5
        disp(['Unequal variance: ',num2str(ageRng(i)),' years'])
    end
    
    % Run permutation tests at each quantile
    [~,~,orig] = permtest_t2(gainNT,gainASD,'nperm',nperm);
    d = mes(gainNT,gainASD,'hedgesg','isDep',0,'nBoot',nperm);
    hedgesg(i) = d.hedgesg;
    bf01(i) = 1/bf.ttest2(gainNT,gainASD);
        
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
plt = subplot(4,5,13:15);
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
[hl,hp] = boundedline(ageRng,rmvAvg(:,1),errNT,ageRng,rmvAvg(:,2),errASD,...
    'alpha','cmap',cols1(1:2,:)); outlinebounds(hl,hp);
plot(ageRng,rmvAvg(:,1),'color',cols1(1,:),'linewidth',lnWth)
plot(ageRng,rmvAvg(:,2),'color',cols1(2,:),'linewidth',lnWth)
xlabel('Age (yrs)','fontname',fnt,'fontsize',fntSz)
ylabel({'Moving mean';['gain {\itk}=7',units]},'fontname',fnt,'fontsize',fntSz)
if fig_label
    text(x1-(x2-x1)*0.15,y2,'f','fontname',fnt,'fontsize',fntSz+2,...
        'fontweight','bold','verticalalign','bottom')
end
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
pos = get(plt,'pos');
set(plt,'pos',[pos(1)+0.025,pos(2)-0.02,pos(3),pos(4)+0.02])

set(gcf,'color','w')

%% RESULTS SECTION

% 4. Predicting multisensory benefits

%% 4.1 ANCOVA Analysis

close all hidden;
clc;

nperm = 1e4;

for i = 1:2
    idx = grp_vec==i;
    dev = categorical(dev_vec(idx),'Ordinal',0);
    ben = benPred(idx); 
    emp_ben = benEmp(idx);
    [~,atab,ctab,stats] = aoctool(ben,emp_ben,dev);
    bf01_1 = 1/bf.bfFromF(cell2mat(atab(2,5)),cell2mat(atab(2,2)),cell2mat(atab(5,2)),sum(idx));
    bf01_2 = 1/bf.bfFromF(cell2mat(atab(3,5)),cell2mat(atab(2,2)),cell2mat(atab(5,2)),sum(idx));
    bf01_3 = 1/bf.bfFromF(cell2mat(atab(4,5)),cell2mat(atab(2,2)),cell2mat(atab(5,2)),sum(idx));
    fprintf('%s:\t\tF(%d,%d) = %1.3f, p = %1.5f, r = %1.3f, R2 = %1.3f, BF01 = %1.19f\n',atab{2,1},cell2mat(atab(2,2)),cell2mat(atab(5,2)),cell2mat(atab(2,5)),cell2mat(atab(2,6)),sqrt(cell2mat(atab(2,3))/(cell2mat(atab(2,3))+cell2mat(atab(5,3)))),cell2mat(atab(2,3))/(cell2mat(atab(2,3))+cell2mat(atab(5,3))),bf01_1)
    fprintf('%s:\t\tF(%d,%d) = %1.3f, p = %1.5f, r = %1.3f, R2 = %1.3f, BF01 = %1.19f\n',atab{3,1},cell2mat(atab(3,2)),cell2mat(atab(5,2)),cell2mat(atab(3,5)),cell2mat(atab(3,6)),sqrt(cell2mat(atab(3,3))/(cell2mat(atab(3,3))+cell2mat(atab(5,3)))),cell2mat(atab(3,3))/(cell2mat(atab(3,3))+cell2mat(atab(5,3))),bf01_2)
    fprintf('%s:\tF(%d,%d) = %1.3f, p = %1.5f, r = %1.3f, R2 = %1.3f, BF01 = %1.19f\n\n',atab{4,1},cell2mat(atab(4,2)),cell2mat(atab(5,2)),cell2mat(atab(4,5)),cell2mat(atab(4,6)),sqrt(cell2mat(atab(4,3))/(cell2mat(atab(4,3))+cell2mat(atab(5,3)))),cell2mat(atab(4,3))/(cell2mat(atab(4,3))+cell2mat(atab(5,3))),bf01_3)
    figure, disp(multcompare(stats))
    for j = 1:4
        idx = grp_vec==i & dev_vec==j;
        [r,stats] = permtest_corr(benPred(idx),benEmp(idx),'nperm',nperm);
        bf01 = 1/bf.corr(benPred(idx),benEmp(idx));
        fprintf('r = %1.3f, R2 = %1.3f, p = %1.5f, BF01 = %1.10f\n',r,r^2,stats.p,bf01)
    end
    fprintf('\n')
end

%% Fig. 4: Modelling multisensory behaviour

close all;
clc;

% Set parameters
bias = [2,2,2;3,3,3;2,2,3;3,2,3];
pBias = 0:0.25:1;
nModels = size(bias,1);
nRatios = numel(pBias);
type = 'pearson';
nperm = 1e4;
fig_label = true;

% Compute bias models
race = mean(CDFs(:,6,6:8,:),3);
simCPs = zeros(nSubjs,nModels,nRatios,numel(prob));
for j = 1:nModels
    model = (CDFs(:,bias(j,1),6,:) + CDFs(:,bias(j,2),7,:) + ...
        CDFs(:,bias(j,3),8,:))/3;
    for k = 1:nRatios
        simCPs(:,j,k,:) = (1-pBias(k))*race + pBias(k)*model;
    end
end

% Generate figure
axes('parent',figure,'fontname',fnt,'fontsize',fntSz)
set(gcf,'pos',[20,50,650,600])

for i = 1:2
    
    % Panel A
    switch test
        case 'ver'
            x1 = 0; x2 = 0.15;
            y1 = -0.04; y2 = 0.15;
        case 'hor'
            x1 = 0; x2 = 100;
            y1 = 0; y2 = 100;
    end
    plt = subplot(3,3,i);
    hold all
    for j = 1:4
        x = benPred(grp_vec==i & dev_vec==j);
        y = benEmp(grp_vec==i & dev_vec==j);
        s = scatter(x,y,15,[shades(j),shades(j),0],'filled');
    end
    for j = 1:4
        x = benPred(grp_vec==i & dev_vec==j);
        y = benEmp(grp_vec==i & dev_vec==j);
        pfit = polyfit(x,y,1);
        xi = linspace(x1,x2) ;
        yi = pfit(1)*xi+pfit(2);
        pl(j) = plot(xi,yi,'color',[shades(j),shades(j),0],'linewidth',2);
    end
    xlim([x1,x2]), ylim([y1,y2])
    set(gca,'layer','top','tickdir','out','xtick',0:0.05:0.15,...
        'ytick',0:0.05:0.15,'fontname',fnt,'fontsize',fntSz)
    title([diag_groups{i}],'fontname',fnt,'fontsize',fntSz)
    pos = get(plt,'pos');
    if i==1
        xlabel({'Predicted';'benefit (AUC)'},'fontname',fnt,'fontsize',fntSz)
        ylabel({'Empirical';'benefit (AUC)'},'fontname',fnt,'fontsize',fntSz)
        if fig_label
            text(x1-(x2-x1)*0.5,y2,'a','fontname',fnt,'fontsize',fntSz+2,...
                'fontweight','bold','verticalalign','bottom')
        end
        set(plt,'pos',[pos(1)-0.02,pos(2),pos(3),pos(4)])
        l = legend(pl,age_groups);
        pos = get(l,'pos');
        set(l,'box','off','fontname',fnt,'fontsize',fntSz-2,...
            'orientation','horizontal','pos',[pos(1)+0.11,pos(2)-0.26,pos(3),pos(4)])
    else
        set(gca,'yticklabel',[])
        set(plt,'pos',[pos(1)-0.05,pos(2),pos(3),pos(4)])
    end
    axis square
    
    % Panel B
    plt = subplot(3,6,i+4);
    x1 = 0.4; x2 = 4.6;
    y1 = -0.1; y2 = 1;
    hold on
    for j = 1:4
        x = benPred(grp_vec==i & dev_vec==j);
        y = benEmp(grp_vec==i & dev_vec==j);
        [r,stats] = permtest_corr(x,y,'type',type,'nperm',nperm);
        disp(['r = ',num2str(r^2)]); disp(['p = ',num2str(stats.p)])
        bar(j,r,'facecolor',[shades(j),shades(j),0],'edgecolor','none')
        if stats.p<0.05
            text(j,r+0.01,'*','fontname',fnt,'fontsize',fntSz+5,...
                'horizontalalign','center')
        end
    end
    pos = get(plt,'pos');
    set(plt,'pos',[pos(1)+0.02,pos(2),pos(3),pos(4)])
    set(gca,'layer','top','tickdir','out','ytick',0:0.2:1,...
        'xtick',1:4,'xticklabel',age_groups,'fontname',fnt,'fontsize',fntSz)
    xtickangle(45)
    title([diag_groups{i}],'fontname',fnt,'fontsize',fntSz)
    if i==1
        ylabel('Correlation ({\itr})','fontname',fnt,'fontsize',fntSz)
        xlabel('Age (yrs)','fontname',fnt,'fontsize',fntSz)
        if fig_label
            text(x1-(x2-x1)*0.7,y2,'b','fontname',fnt,'fontsize',fntSz+2,...
                'fontweight','bold','verticalalign','bottom')
        end
    else
        set(gca,'yticklabel',[])
    end
    xlim([x1,x2])
    ylim([y1,y2])
    
end

% Panel C & D
models = {'Model 1A','Model 1V','Model 2A','Model 2V'};
labels = {'c','d'};
sp = [1,2,5,6;3,4,7,8]+8;
x1 = 0; x2 = 6;
y1 = -0.2; y2 = 0.6;
ytick = -0.2:0.2:0.6;
ctr = 1;
dev1 = 1;
dev2 = 2;
for n = 1:2
    
    for i = 1:nModels
        
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
        
        plt = subplot(4,4,sp(n,i));
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
            pos = get(plt,'pos');
            set(plt,'pos',[pos(1)+0.04,pos(2),pos(3),pos(4)])
        else
            set(gca,'ytick',ytick)
            pos = get(plt,'pos');
            set(plt,'pos',[pos(1)-0.02,pos(2),pos(3),pos(4)])
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
            pos = get(plt,'pos');
            set(plt,'pos',[pos(1),pos(2)+0.02,pos(3),pos(4)])
        else
            pos = get(plt,'pos');
            set(plt,'pos',[pos(1),pos(2)+0.02,pos(3),pos(4)])
        end
        if i==2 || i==4
            set(gca,'yticklabel',[])
            pos = get(plt,'pos');
            set(plt,'pos',[pos(1)-0.02,pos(2),pos(3),pos(4)])
        elseif i==3
            pos = get(plt,'pos');
            set(plt,'pos',[pos(1),pos(2),pos(3),pos(4)])
            xlabel('Race\leftarrow{\itp}\rightarrowBias','fontname',fnt,'fontsize',fntSz)
            ylabel('Correlation ({\itr})','fontname',fnt,'fontsize',fntSz)
        end
        if ctr==1
            l = legend(pl,diag_groups);
            pos = get(l,'pos');
            set(l,'box','off','fontname',fnt,'fontsize',fntSz-2,...
                'pos',[pos(1)+0.02,pos(2)+0.03,pos(3),pos(4)])
        elseif ctr==4
            dev1 = dev2;
        end
        
        ctr = ctr+1;
        
    end
end

set(gcf,'color','w')

%% Fig. 5: Sesnory dominance

close all;
clc;

% Set parameters
bias = [2,2,3;3,2,3;2,2,2;3,3,3];
pBias = 0:0.25:1;
nModels = size(bias,1);
nRatios = numel(pBias);
fig_label = true;

% Compute bias model
race = mean(CDFs(:,6,6:8,:),3);
simCPs = zeros(nSubjs,nModels,nRatios,numel(prob));
for j = 1:nModels
    model = (CDFs(:,bias(j,1),6,:) + CDFs(:,bias(j,2),7,:) + ...
        CDFs(:,bias(j,3),8,:))/3;
    for k = 1:nRatios
        simCPs(:,j,k,:) = (1-pBias(k))*race + pBias(k)*model;
    end
end

% Generate figure
axes('parent',figure,'fontname',fnt,'fontsize',fntSz)
set(gcf,'pos',[20,150,400,450],'paperposmode','auto')

% Panel C & D
y1 = -0.2; y2 = 0.6;
x1 = 0; x2 = 5;
for n = 1:2
    
    rVals = zeros(4,2); pVals = zeros(4,2);
    
    for i = 3:nModels
        
        benMod = simCPs(:,i,5,:)-max(CDFs(:,2:3,1,:),[],2);
        benMod = squeeze(trapz(prob,benMod,4));
        
        for j = 1:4
            x = benMod(grp_vec==n & dev_vec==j);
            y = benEmp(grp_vec==n & dev_vec==j);
            [rVals(j,i-2),stats] = permtest_corr(x,y); pVals(j,i-2) = stats.p;
        end
        
    end
    
    plt = subplot(2,2,n);
    hold all
    plot([0,5],[0,0],':k'); pl = [];
    pl(1) = plot(1:4,rVals(:,1),'-k','linewidth',2);
    pl(2) = plot(1:4,rVals(:,2),'--k','linewidth',2);
    plot(1:4,rVals(:,1),'o','markerfacecolor','k','markeredgecolor','k','markersize',5);
    plot(1:4,rVals(:,2),'o','markerfacecolor','k','markeredgecolor','k','markersize',5);
    set(gca,'tickdir','out','xtick',1:4,'xticklabel',age_groups,...
        'fontname',fnt,'fontsize',fntSz)
    ylim([y1,y2])
    xlim([x1,x2])
    xtickangle(45)
    title(diag_groups{n},'fontname',fnt,'fontsize',fntSz)
    axis square
    if n==1
        ylabel('Correlation ({\itr})','fontname',fnt,'fontsize',fntSz)
    else
        set(gca,'yticklabel',[])
        l = legend(pl,'A-bias','V-bias');
        pos = get(l,'pos');
        set(l,'box','off','pos',[pos(1)+0.045,pos(2)+0.04,pos(3),pos(4)],...
            'fontname',fnt,'fontsize',fntSz-2)
    end
    if n==1
        if fig_label
            text(x1-(x2-x1)*0.35,y2,'a','fontname',fnt,'fontsize',fntSz+2,...
                'fontweight','bold','verticalalign','bottom')
        end
    end
    
end

y1 = 0; y2 = 2;
x1 = 0; x2 = 5;
for n = 1:2
    
    rVals = zeros(4,2); pVals = zeros(4,2);
    
    for i = 1:2
        
        for j = 1:4
            idx = grp_vec==n & dev_vec==j;
            if i==1
                tmp = squeeze(CDFs(idx,1,2,:)-CDFs(idx,1,4,:));
                tmp1 = squeeze(CDFs(idx,1,2,:)-CDFs(idx,1,3,:));
                tmp2 = mean(trapz(prob,tmp1,2));
                rVals(j,i) = mean(trapz(prob,tmp,2))/tmp2;
            else
                tmp = squeeze(CDFs(idx,1,2,:)-CDFs(idx,1,5,:));
                tmp1 = squeeze(CDFs(idx,1,2,:)-CDFs(idx,1,3,:));
                tmp2 = mean(trapz(prob,tmp1,2));
                rVals(j,i) = mean(trapz(prob,tmp,2))/tmp2;
            end
        end
        
    end
    
    plt = subplot(2,2,n+2);
    hold all
    pl(1) = plot(1:4,rVals(:,1),'-k','linewidth',2);
    pl(2) = plot(1:4,rVals(:,2),'--k','linewidth',2);
    plot(1:4,rVals(:,1),'o','markerfacecolor','k','markeredgecolor','k','markersize',5);
    plot(1:4,rVals(:,2),'o','markerfacecolor','k','markeredgecolor','k','markersize',5);
    set(gca,'tickdir','out','xtick',1:4,'xticklabel',age_groups,...
        'fontname',fnt,'fontsize',fntSz)
    ylim([y1,y2])
    xlim([x1,x2])
    xtickangle(45)
    pos = get(plt,'pos');
    set(gca,'pos',[pos(1),pos(2)+0.05,pos(3),pos(4)])
    axis square
    if n==1
        xlabel('Age (yrs)','fontname',fnt,'fontsize',fntSz)
        ylabel('MSE (%)','fontname',fnt,'fontsize',fntSz)
        ylabel('Normalized MSE','fontname',fnt,'fontsize',fntSz)
    else
        set(gca,'yticklabel',[])
        l = legend(pl,'A\rightarrowAV','V\rightarrowAV');
        pos = get(l,'pos');
        set(l,'box','off','pos',[pos(1)+0.045,pos(2)+0.05,pos(3),pos(4)],...
            'fontname',fnt,'fontsize',fntSz-2)
    end
    if n==1
        if fig_label
            text(x1-(x2-x1)*0.35,y2,'b','fontname',fnt,'fontsize',fntSz+2,...
                'fontweight','bold','verticalalign','bottom')
        end
    end
    
end

set(gcf,'color','w')

%% RESULTS SECTION

% 6. Modality switch effects

%% 6.1 Linear regression analysis

close all;
clc;

% Define variables
subj = [1:nSubjs,1:nSubjs,1:nSubjs]';
grp = [grp_vec;grp_vec;grp_vec];
sex = [sex_vec;sex_vec;sex_vec];
age = [age_vec;age_vec;age_vec];
dev = [dev_vec;dev_vec;dev_vec];
cond = [ones(nSubjs,1);2*ones(nSubjs,1);3*ones(nSubjs,1)];
sc = cost(:);

% Convert non-numeric variables to categorical
subj = categorical(subj,'Ordinal',0);
grp = categorical(grp,[1,2],{'TD','ASD'},'Ordinal',0);
sex = categorical(sex,[1,2],{'Male','Fem'},'Ordinal',0);
dev = categorical(dev,[1,2,3,4],{'6-9','10-12','13-17','18-40'},'Ordinal',1);
cond = categorical(cond,[1,2,3],{'AV','A','V'},'Ordinal',0);

% Create dataset
ds = table(subj,grp,sex,age,dev,cond,sc);

% Fit linear model to dataset
disp('Fitting linear model...')
lm = fitlm(ds,'sc ~ grp + age + cond');

% Perform hypothesis tests on fixed effect terms
disp('Performing hypothesis tests on fixed effect terms...')
stats = anova(lm);

% Display LME and ANOVA stats
disp(lm)
disp(stats)

% Define residuals
r = lm.Residuals.Raw;
pr = lm.Residuals.Pearson;
st = lm.Residuals.Standardized;

% Bayes factor analysis
models = {'sc ~ grp + age + cond','sc ~ age + cond','sc ~ grp + cond','sc ~ grp + age'};
bf01 = 1./bf.anova(ds,models);
fprintf('\nBayes factor analysis:')
fprintf('\nBF01 (full): %d',bf01(1))
fprintf('\nBF01 (grp): %d',bf01(1)/bf01(2))
fprintf('\nBF01 (age): %d\',bf01(1)/bf01(3))
fprintf('\nBF01 (cond): %d\n\n',bf01(1)/bf01(4))
bf01 = 1/bf.bfFromF(F,df1,df2,n);
fprintf('\nBF01 (cond_A): %d\n\n',bf01(1)/bf01(4))
bf01 = 1/bf.bfFromF(F,df1,df2,n);
fprintf('\nBF01 (cond_V): %d\n\n',bf01(1)/bf01(4))

figure,
subplot(2,2,1), histfit(r), xlabel('Residuals'), ylabel('Frequency')
title('Histogram of residuals with fitted normal density')
subplot(2,2,2), qqplot(r)
subplot(2,2,3), plotResiduals(lm,'fitted')
subplot(2,2,4), plotResiduals(lm,'lagged')

figure, boxplot([r,pr,st])

mc_normalitytest(r);

%% 6.2 Permutation test à la Williams et al. (2013)

close all;
clc;

ageLwr = 8;
ageUpr = 15;
ct_meth = 'mean';
nperm = 1e4;
nTasks = 13;
tail = 'left';
age_match = true;
isi_type = 'long';

% Compute RT central tendancy for each subject, condition, trial type
RTct = zeros(nSubjs,3,nTasks);
for i = 1:nSubjs
    for j = 1:3
        for k = 1:nTasks
            if strcmp(ct_meth,'mean')
                RTct(i,j,k) = mean(RTs{i,j,k});
            elseif strcmp(ct_meth,'med')
                RTct(i,j,k) = median(RTs{i,j,k});
            elseif strcmp(ct_meth,'mode')
                RTct(i,j,k) = mode(round(RTs{i,j,k}));
            end
        end
    end
end

ageRng = 12;
window = 8;
i = 1;
switch age_match
    case true
        [idx1,idx2] = movmeanmatch(i,grp_vec,sex_vec,age_vec,piq_vec,ageRng,window,false);
    case false
        idx1 = find(grp_vec==1 & age_vec>ageRng(i)-window/2 & age_vec<ageRng(i)+window/2);
        idx2 = find(grp_vec==2 & age_vec>ageRng(i)-window/2 & age_vec<ageRng(i)+window/2);
end

switch isi_type
    case 'all'
        costNT = squeeze(RTct(idx1,1:3,3)-RTct(idx1,1:3,2)); % all ISIs
        costASD = squeeze(RTct(idx2,1:3,3)-RTct(idx2,1:3,2)); % all ISIs
    case 'long'
        costNT = squeeze(RTct(idx1,1:3,12)-RTct(idx1,1:3,11)); % long ISIs
        costASD = squeeze(RTct(idx2,1:3,12)-RTct(idx2,1:3,11)); % long ISIs
    case 'short'
        costNT = squeeze(RTct(idx1,1:3,10)-RTct(idx1,1:3,9)); % short ISIs
        costASD = squeeze(RTct(idx2,1:3,10)-RTct(idx2,1:3,9)); % short ISIs
end

[~,stats] = permtest_f2(costNT,costASD,'nperm',nperm,'tail',tail);
stats.h = 0;

% Run permutation test with tmax correction
if stats.h==0
    [tstat,stats,orig,param] = permtest_t2(costNT,costASD,'nperm',nperm,...
        'varx','equal','tail',tail);
    d = mes(costNT,costASD,'hedgesg','isDep',0,'nBoot',nperm);
    for j = 1:3
        bf01 = 1/bf.ttest2(costNT(:,j),costASD(:,j),'scale',1,'tail',tail);
        fprintf('MSEnt = %0.3f ± %0.3f, ',mean(costNT(:,j)),std(costNT(:,j)))
        fprintf('MSEasd = %0.3f ± %0.3f, ',mean(costASD(:,j)),std(costASD(:,j)))
        fprintf('t(%d) = %0.2f, ',param.df(j),tstat(j))
        fprintf('p = %0.3f, ',orig.p(j))
        fprintf('g = %0.2f, ',d.hedgesg(j))
        fprintf('95CI [%0.2f, %0.2f], ',d.hedgesgCi(:,j))
        fprintf('BF01 = %0.2f\n',bf01)
    end
elseif stats.h==1
    [tstat,stats,orig,param] = permtest_t2(costNT,costASD,'nperm',nperm,...
        'varx','unequal','tail',tail);
    d = mes(costNT,costASD,'glassdelta','isDep',0,'nBoot',nperm);
    d.glassdelta = d.glassdelta*gamma((length(costNT(:,j))-1)/2)/(sqrt((length(costNT(:,j))-1)/2)*gamma((length(costNT(:,j))-2)/2));
    for j = 1:3
        bf01 = 1/bf.ttest2(costNT(:,j),costASD(:,j),'scale',1,'tail',tail);
        fprintf('MSEnt = %0.3f ± %0.3f, ',mean(costNT(:,j)),std(costNT(:,j)))
        fprintf('MSEasd = %0.3f ± %0.3f, ',mean(costASD(:,j)),std(costASD(:,j)))
        fprintf('t(%0.2f) = %0.2f, ',param.df(j),tstat(j))
        fprintf('p = %0.3f, ',orig.p(j))
        fprintf('g = %0.2f, ',d.glassdelta(j))
        fprintf('95CI [%0.2f, %0.2f], ',d.glassdeltaCi(:,j))
        fprintf('BF01 = %0.2f\n',bf01)
    end
end

%%
sc1 = squeeze(RTct(idx1,3,10)-RTct(idx1,3,9)); % short ISIs
sc2 = squeeze(RTct(idx1,3,12)-RTct(idx1,3,11)); % long ISIs
sc3 = squeeze(RTct(idx2,3,10)-RTct(idx2,3,9)); % short ISIs
sc4 = squeeze(RTct(idx2,3,12)-RTct(idx2,3,11)); % long ISIs
sc = [sc1;sc2;sc3;sc4];

% Define variables
subs = (1:nSubjs)';
n = numel(idx1);
subj = [subs(idx1);subs(idx1);subs(idx2);subs(idx2)];
grp = [grp_vec(idx1);grp_vec(idx1);grp_vec(idx2);grp_vec(idx2)];
sex = [sex_vec(idx1);sex_vec(idx1);sex_vec(idx2);sex_vec(idx2)];
age = [age_vec(idx1);age_vec(idx1);age_vec(idx2);age_vec(idx2)];
dev = [dev_vec(idx1);dev_vec(idx1);dev_vec(idx2);dev_vec(idx2)];
isi = [ones(n,1);2*ones(n,1);ones(n,1);2*ones(n,1)];

% Convert non-numeric variables to categorical
subj = categorical(subj,'Ordinal',0);
grp = categorical(grp,[1,2],{'TD','ASD'},'Ordinal',0);
sex = categorical(sex,[1,2],{'Male','Fem'},'Ordinal',0);
dev = categorical(dev,[1,2,3,4],{'6-9','10-12','13-17','18-40'},'Ordinal',1);
isi = categorical(isi,[1,2],{'Short','Long'},'Ordinal',1);

% Create dataset
ds = dataset(subj,grp,sex,age,dev,isi,sc);

% Fit linear model to dataset
disp('Fitting linear model...')
% lm = fitlm(ds,'sc ~ age + grp + isi');
lm = fitlm(ds,'sc ~ isi');

% Perform hypothesis tests on fixed effect terms
disp('Performing hypothesis tests on fixed effect terms...')
stats = anova(lm);

% Display LME and ANOVA stats
disp(lm)
disp(stats)

% Define residuals
r = lm.Residuals.Raw;
pr = lm.Residuals.Pearson;
st = lm.Residuals.Standardized;

figure,
subplot(2,2,1), histfit(r), xlabel('Residuals'), ylabel('Frequency')
title('Histogram of residuals with fitted normal density')
subplot(2,2,2), qqplot(r)
subplot(2,2,3), plotResiduals(lm,'fitted')
subplot(2,2,4), plotResiduals(lm,'lagged')

figure, boxplot([r,pr,st])

mc_normalitytest(r);

%% Long vs short ISIs

clc;

tail = 'both';

sc1 = squeeze(RTct(idx1,3,3)-RTct(idx1,3,2)); % all ISIs
sc2 = squeeze(RTct(idx1,3,12)-RTct(idx1,3,11)); % long ISIs
sc3 = squeeze(RTct(idx2,3,3)-RTct(idx2,3,2)); % all ISIs
sc4 = squeeze(RTct(idx2,3,12)-RTct(idx2,3,11)); % long ISIs

scAll = sc1;
scLong = sc2;

j=1;
[tstat,stats,orig,param] = permtest_t(scAll,scLong,'nperm',nperm,...
    'varx','equal','tail',tail);
d = mes(scAll,scLong,'hedgesg','isDep',1,'nBoot',nperm);
bf01 = 1/bf.ttest(scAll(:,j),scLong(:,j),'scale',1,'tail',tail);
fprintf('MSEall = %0.3f ± %0.3f, ',mean(scAll(:,j)),std(scAll(:,j)))
fprintf('MSElong = %0.3f ± %0.3f, ',mean(scLong(:,j)),std(scLong(:,j)))
fprintf('t(%d) = %0.2f, ',param.df(j),tstat(j))
fprintf('p = %0.3f, ',orig.p(j))
fprintf('g = %0.2f, ',d.hedgesg(j))
fprintf('95CI [%0.2f, %0.2f], ',d.hedgesgCi(:,j))
fprintf('BF01 = %0.2f\n',bf01)

scAll = sc3;
scLong = sc4;

j=1;
[tstat,stats,orig,param] = permtest_t(scAll,scLong,'nperm',nperm,...
    'varx','equal','tail',tail);
d = mes(scAll,scLong,'hedgesg','isDep',1,'nBoot',nperm);
bf01 = 1/bf.ttest(scAll(:,j),scLong(:,j),'scale',1,'tail',tail);
fprintf('MSEall = %0.3f ± %0.3f, ',mean(scAll(:,j)),std(scAll(:,j)))
fprintf('MSElong = %0.3f ± %0.3f, ',mean(scLong(:,j)),std(scLong(:,j)))
fprintf('t(%d) = %0.2f, ',param.df(j),tstat(j))
fprintf('p = %0.3f, ',orig.p(j))
fprintf('g = %0.2f, ',d.hedgesg(j))
fprintf('95CI [%0.2f, %0.2f], ',d.hedgesgCi(:,j))
fprintf('BF01 = %0.2f\n',bf01)


%% Fig. 6: Modality switch effects

close all;
clc;

% Set parameters
nperm = 1e4;
ct_meth = 'mean';
lnWth = 2;
age_match = true;
fig_label = true;
fdr_method = 'bh'; % 'bh' or 'bky' or 'none'

% Generate figure
axes('parent',figure,'fontname',fnt,'fontsize',fntSz)
set(gcf,'pos',[20,50,800,400])

% Panel A
x1 = 0.5; x2 = 3.5;
y1 = -0.01; y2 = 0.105;
for i = 1:4
    
    switch age_match
        case true
            if i==4
                [idx1,idx2] = agematch(i,grp_vec,dev_vec,sex_vec,age_vec,[],false);
            else
                [idx1,idx2] = agematch(i,grp_vec,dev_vec,sex_vec,age_vec,piq_vec,false);
            end
        case false
            idx1 = find(grp_vec==1 & dev_vec==i);
            idx2 = find(grp_vec==2 & dev_vec==i);
    end
    
    costNT = cost(idx1,:);
    costASD = cost(idx2,:);
    
    [~,~,stats] = permtest_f2(costNT,costASD);
    
    % Run permutation test with tmax correction
    if stats.h==0
        [tstat,stats,~,param] = permtest_t2(costNT,costASD,'nperm',nperm);
        d = mes(costNT,costASD,'hedgesg','isDep',0,'nBoot',nperm);
        for j = 1:3
            bf01 = 1/bf.ttest2(costNT(:,j),costASD(:,j));
            fprintf('SCnt = %0.3f ± %0.3f, ',mean(costNT(:,j)),std(costNT(:,j)))
            fprintf('SCasd = %0.3f ± %0.3f, ',mean(costASD(:,j)),std(costASD(:,j)))
            fprintf('t(%d) = %0.2f, ',param.df(j),tstat(j))
            fprintf('p = %0.3f, ',stats.p(j))
            fprintf('g = %0.2f, ',d.hedgesg(j))
            fprintf('95CI [%0.2f, %0.2f], ',d.hedgesgCi(:,j))
            fprintf('BF01 = %0.2f\n',bf01)
        end
    elseif stats.h==1
        [tstat,stats,~,param] = permtest_t2(costNT,costASD,'nperm',nperm,'varx','unequal');
        d = mes(costNT,costASD,'glassdelta','isDep',0,'nBoot',nperm);
        %             d.glassdelta = d.glassdelta*gamma((length(areaNT(:,j))-1)/2)/(sqrt((length(areaNT(:,j))-1)/2)*gamma((length(areaNT(:,j))-2)/2));
        for j = 1:3
            bf01 = 1/bf.ttest2(costNT(:,j),costASD(:,j));
            fprintf('SCnt = %0.3f ± %0.3f, ',mean(costNT(:,j)),std(costNT(:,j)))
            fprintf('SCasd = %0.3f ± %0.3f, ',mean(costASD),std(costASD))
            fprintf('t(%0.2f) = %0.2f, ',param.df(j),tstat(j))
            fprintf('p = %0.3f, ',stats.p(j))
            fprintf('g = %0.2f, ',d.glassdelta(j))
            fprintf('95CI [%0.2f, %0.2f], ',d.glassdeltaCi(:,j))
            fprintf('BF01 = %0.2f\n',bf01)
        end
    end
    
    % Compute group SC central tendancy and 95% CIs
    if strcmp(ct_meth,'mean')
        scNTavg = mean(costNT);
        scASDavg = mean(costASD);
        ciNT = bootci(nperm,@mean,costNT);
        ciASD = bootci(nperm,@mean,costASD);
    elseif strcmp(ct_meth,'med')
        scNTavg = median(costNT);
        scASDavg = median(costASD);
        ciNT = bootci(nperm,@median,costNT);
        ciASD = bootci(nperm,@median,costASD);
    elseif strcmp(ct_meth,'mode')
        scNTavg = mode(costNT);
        scASDavg = mode(costASD);
        ciNT = bootci(nperm,@mode,costNT);
        ciASD = bootci(nperm,@mode,costASD);
    end
    errNT = zeros(2,3);
    errASD = zeros(2,3);
    errNT(1,:) = scNTavg-ciNT(1,:);
    errNT(2,:) = ciNT(2,:)-scNTavg;
    errASD(1,:) = scASDavg-ciASD(1,:);
    errASD(2,:) = ciASD(2,:)-scASDavg;
    
    subplot(2,4,i)
    hold all
    plot(1:3,scNTavg+1e3,'color',cols1(1,:),'linewidth',lnWth)
    plot(1:3,scASDavg+1e3,'color',cols1(2,:),'linewidth',lnWth)
    plot([x1,x2],[0,0],':k','linewidth',1)
    errorbar((1:3)-0.1,scNTavg,errNT(1,:),errNT(2,:),'color',cols1(1,:),...
        'linewidth',lnWth,'linestyle','none','marker','.','MarkerSize',10)
    errorbar((1:3)+0.1,scASDavg,errASD(1,:),errASD(2,:),'color',cols1(2,:),...
        'linewidth',lnWth,'linestyle','none','marker','.','MarkerSize',10)
    %     errorbar((1:3),scNTavg,errNT(1,:),errNT(2,:),'color',cols1(1,:),...
    %         'linewidth',lnWth,'linestyle','-')
    %     errorbar((1:3),scASDavg,errASD(1,:),errASD(2,:),'color',cols1(2,:),...
    %         'linewidth',lnWth,'linestyle','-')
    set(gca,'layer','top','box','off','tickdir','out','xtick',1:3,...
        'xticklabel',switch_conds,'ytick',0:0.05:0.15)
    title([age_groups{i},' yrs'],'fontname',fnt,'fontsize',fntSz)
    for j = 1:3
        if stats.p(j) < 0.001
            text(j,scNTavg(j)+errNT(1,j)+0.01,'***','fontname',fnt,...
                'fontsize',fntSz+8,'horizontalalign','center')
        elseif stats.p(j) < 0.01
            text(j,scNTavg(j)+errNT(1,j)+0.01,'**','fontname',fnt,...
                'fontsize',fntSz+8,'horizontalalign','center')
        elseif stats.p(j) < 0.05
            text(j,scNTavg(j)+errNT(1,j)+0.01,'*','fontname',fnt,...
                'fontsize',fntSz+8,'horizontalalign','center')
        end
    end
    if i==1
        ylabel('MSE (AUC)','fontname',fnt,'fontsize',fntSz)
        if fig_label
            text(x1-(x2-x1)*0.5,y2,'a','fontname',fnt,'fontsize',fntSz+2,...
                'fontweight','bold','verticalalign','bottom')
        end
        l = legend(diag_groups);
        pos = get(l,'pos');
        set(l,'box','off','fontname',fnt,'fontsize',fntSz-2,...
            'pos',[pos(1)-0.03,pos(2)+0.04,pos(3),pos(4)])
    else
        set(gca,'yticklabel',[])
    end
    xlim([x1,x2]), ylim([y1,y2])
    xtickangle(45)
    axis square
    
end

% Panel B

% Compute moving mean MSE and run perm test
ageRng = 6:26;
window = 7;
numAge = numel(ageRng);
x1 = ageRng(1); x2 = ageRng(end);
y1 = 0; y2 = 0.1;
ctr = 3;
for n = 2:3
    
    % Preallocate memory
    scAvg = zeros(numAge,2);
    errNT = zeros(numAge,2);
    errASD = zeros(numAge,2);
    p = zeros(1,numel(ageRng));
    h = zeros(1,numel(ageRng));
    dataNT = nan(81,numel(ageRng));
    dataASD = nan(81,numel(ageRng));
    
    % Loop through age bins
    for i = 1:length(ageRng)
        
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
        
        costNT = cost(idx1,n);
        costASD = cost(idx2,n);
        
        % Compute mean MSC and CI
        scAvg(i,1) = mean(costNT);
        ci = bootci(nperm,@mean,costNT);
        errNT(i,1) = scAvg(i,1)-ci(1);
        errNT(i,2) = ci(2)-scAvg(i,1); clear ci
        
        scAvg(i,2) = mean(costASD);
        ci = bootci(nperm,@mean,costASD);
        errASD(i,1) = scAvg(i,2)-ci(1);
        errASD(i,2) = ci(2)-scAvg(i,2); clear ci
        
        % Run permutation test
        [~,stats] = permtest_t2(costNT,costASD,nperm);
        p(i) = stats.p;
        h(i) = stats.h;
        
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
    subplot(2,4,n+ctr:n+ctr+1)
    hold all
    plot(ageRng,scAvg(:,1)+1e3,'color',cols1(1,:),'linewidth',lnWth)
    plot(ageRng,scAvg(:,2)+1e3,'color',cols1(2,:),'linewidth',lnWth)
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
    [hl,hp] = boundedline(ageRng,scAvg(:,1),errNT,ageRng,scAvg(:,2),errASD,...
        'alpha','cmap',cols1(1:2,:)); outlinebounds(hl,hp);
    plot(ageRng,scAvg(:,1),'color',cols1(1,:),'linewidth',lnWth)
    plot(ageRng,scAvg(:,2),'color',cols1(2,:),'linewidth',lnWth)
    set(gca,'layer','top','box','off','tickdir','out',...
        'xtick',6:2:26,'ytick',0:0.02:0.1)
    if n==2
        xlabel('Age (yrs)','fontname',fnt,'fontsize',fntSz)
        title('Visual\rightarrowAuditory','fontname',fnt,'fontsize',fntSz)
        ylabel({'Moving mean';'MSE {\itk}=7 (AUC)'},'fontname',fnt,'fontsize',fntSz)
        l = legend(diag_groups);
        set(l,'box','off','location','northwest','fontname',fnt,'fontsize',fntSz-2)
        if fig_label
            text(x1-(x2-x1)*0.2,y2,'b','fontname',fnt,'fontsize',fntSz+2,...
                'fontweight','bold','verticalalign','bottom')
        end
    else
        title('Auditory\rightarrowVisual','fontname',fnt,'fontsize',fntSz)
        set(gca,'yticklabel',[])
    end
    xlim([x1,x2])
    ylim([y1,y2])
    
    ctr = ctr+1;
    
end

set(gcf,'color','w')

%% RESULTS SECTION

% 7. Increased channel dependency in ASD

%% 7.1 Fit Otto's model

% Set parameters
qnts = 2:length(prob);
nQnts = length(prob);
laterfit = zeros(nSubjs,2,2);
laterCDFs = zeros(nSubjs,2,nQnts);
racefit = zeros(nSubjs,2);
race = zeros(nSubjs,nQnts);
x = zeros(nSubjs,nQnts);
logL = zeros(nSubjs,1);

for i = 1:nSubjs
    
    % Get subjects data
    data = RTs(i,1:3,1);
    
    % Create quantiles
    x(i,:) = linspace(perLim(i,1)/1e3,perLim(i,2)/1e3,nQnts);
    
    % Fitting the LATER model to the single signal conditions
    for j = 1:2
        laterfit(i,j,:) = fitLater(data{j+1}/1e3);
        laterCDFs(i,j,qnts) = laterCDF(x(i,qnts),squeeze(laterfit(i,j,1)),squeeze(laterfit(i,j,2)));
    end
    
    % Fitting the RACE model to the redundant signals condition
    [racefit(i,:),logL(i)] = fitRace(data{1}/1e3,squeeze(laterfit(i,:,1)),squeeze(laterfit(i,:,2)));
    race(i,qnts) = raceCDF(x(i,qnts),squeeze(laterfit(i,:,1)),squeeze(laterfit(i,:,2)),racefit(i,1),racefit(i,2));
    
end

rho = racefit(:,1);
eta = racefit(:,2);

%% 7.2 Linear regression analysis

close all;
clc;

% Define variables
subj = 1:nSubjs;
grp = grp_vec;
sex = sex_vec;
age = age_vec;
dev = dev_vec;

% Convert non-numeric variables to categorical
subj = categorical(subj','Ordinal',0);
grp = categorical(grp,[1,2],{'TD','ASD'},'Ordinal',0);
sex = categorical(sex,[1,2],{'Male','Fem'},'Ordinal',0);
dev = categorical(dev,[1,2,3,4],{'6-9','10-12','13-17','18-40'},'Ordinal',1);

% Create dataset
ds = dataset(subj,grp,sex,age,dev,rho,eta);

% Fit linear model to dataset
fprintf('Fitting linear model...\n')
lm = fitlm(ds,'rho ~ grp + age');
disp(lm)

% Perform hypothesis tests on fixed effect terms
fprintf('Performing hypothesis tests on fixed effect terms...\n\n')
stats = anova(lm);
disp(stats)

% Define residuals
r = lm.Residuals.Raw;
pr = lm.Residuals.Pearson;
st = lm.Residuals.Standardized;

fprintf('Performing goodness-of-fit tests on residuals...\n\n')
mc_normalitytest(r);

figure,
subplot(2,2,1), histfit(r), xlabel('Residuals'), ylabel('Frequency')
title('Histogram of residuals with fitted normal density')
subplot(2,2,2), qqplot(r)
subplot(2,2,3), plotResiduals(lm,'fitted')
subplot(2,2,4), plotResiduals(lm,'lagged')

figure, boxplot([r,pr,st])

%% 7.3 Permutation test

% clc
nperm = 1e4;

ctr = 1;
for i = 1:4
    
    paramNT = rho_cell{1,ctr};
    paramASD = rho_cell{1,ctr+1};
    
    [~,~,orig] = permtest_f2(paramNT,paramASD,'nperm',nperm);
    bf01 = 1/bf.ttest2(paramNT,paramASD);
    
    % Run permutation test with tmax correction
    if orig.h==0
        [tstat,~,orig,stats] = permtest_t2(paramNT,paramASD,'nperm',nperm,'varx','equal');
        d = mes(paramNT,paramASD,'hedgesg','isDep',0,'nBoot',nperm);
        fprintf('SCnt = %0.3f ± %0.3f, ',mean(paramNT),std(paramNT))
        fprintf('SCasd = %0.3f ± %0.3f, ',mean(paramASD),std(paramASD))
        fprintf('t(%d) = %0.2f, ',stats.df,tstat)
        fprintf('p = %0.3f, ',orig.p)
        fprintf('g = %0.2f, ',d.hedgesg)
        fprintf('95CI [%0.2f, %0.2f], ',d.hedgesgCi)
        fprintf('BF01 = %0.2f\n',bf01)
    elseif orig.h==1
        disp('Yurt!')
        [tstat,~,orig,stats] = permtest_t2(paramNT,paramASD,'nperm',nperm,'varx','unequal');
        d = mes(paramNT,paramASD,'glassdelta','isDep',0,'nBoot',nperm);
        n = numel(paramNT);
        d.glassdelta = d.glassdelta*gamma((n-1)/2)/(sqrt((n-1)/2)*gamma((n-2)/2));
        fprintf('SCnt = %0.3f ± %0.3f, ',mean(paramNT),std(paramNT))
        fprintf('SCasd = %0.3f ± %0.3f, ',mean(paramASD),std(paramASD))
        fprintf('t(%0.2f) = %0.2f, ',stats.df,tstat)
        fprintf('p = %0.3f, ',orig.p)
        fprintf('g = %0.2f, ',d.glassdelta)
        fprintf('95CI [%0.2f, %0.2f], ',d.glassdeltaCi)
        fprintf('BF01 = %0.2f\n',bf01)
    end
    
    ctr = ctr+2;
    
end


%% Fig. 7: LATER Race Model

close all;
clc;

% Define parameters
nperm = 1e4;
xaxis = prob;
idxs = cell(2,4);
rho_cell = cell(1,8);
eta_cell = cell(1,8);
ctr = 0;
rVal = 'r';
age_match = true;
fig_label = true;

for i = 1:4
    
    switch age_match
        case true
            if i==4
                [idx1,idx2] = agematch(i,grp_vec,dev_vec,sex_vec,age_vec,[],false);
            else
                [idx1,idx2] = agematch(i,grp_vec,dev_vec,sex_vec,age_vec,piq_vec,false);
            end
        case false
            idx1 = find(grp_vec==1 & dev_vec==i);
            idx2 = find(grp_vec==2 & dev_vec==i);
    end
    
    rho_cell{i+ctr} = rho(idx1);
    rho_cell{i+ctr+1} = rho(idx2);
    eta_cell{i+ctr} = eta(idx1);
    eta_cell{i+ctr+1} = eta(idx2);
    
    idxs{1,i} = idx1;
    idxs{2,i} = idx2;
    
    ctr = ctr+1;
    
end

factor = 1;
quant = 2:length(prob);
quant2 = 1:factor:length(prob)*factor-(factor-1);

idx = grp_vec==1 & dev_vec==4;

col = [cols1(1:2,:);cols1(1:2,:);cols1(1:2,:);cols1(1:2,:)];

a2afreq = squeeze(nRTs(:,2,7,quant)./sum(nRTs(:,2,7:8,quant),3))*100;
a2vfreq = squeeze(nRTs(:,3,7,quant)./sum(nRTs(:,3,7:8,quant),3))*100;
v2vfreq = squeeze(nRTs(:,3,8,quant)./sum(nRTs(:,3,7:8,quant),3))*100;
v2afreq = squeeze(nRTs(:,2,8,quant)./sum(nRTs(:,2,7:8,quant),3))*100;

% Generate figure
axes('parent',figure,'fontname',fnt,'fontsize',fntSz)
set(gcf,'pos',[20,50,850,450],'paperposmode','auto')

% Panel A
plt = subplot(2,4,1);
x1 = 48; x2 = 67;
y1 = 20; y2 = 52;
hold on
plotpolyci(nanmean(a2afreq(idx,:),1),nanmean(a2vfreq(idx,:),1),2)
for i = 1:length(quant)
    scatter(nanmean(a2afreq(idx,i),1),nanmean(a2vfreq(idx,i),1),mkrSz+20,[1,1,1]*i/length(quant),'fill','markeredgecolor','k')
end
set(gca,'box','off','tickdir','out')
xlim([x1,x2]);
ylim([y1,y2]);
title('Chan. Dependency','fontname',fnt,'fontsize',fntSz)
xlabel('A\rightarrowA frequency (%)','fontname',fnt,'fontsize',fntSz)
ylabel('A\rightarrowV frequency (%)','fontname',fnt,'fontsize',fntSz)
if fig_label
    text(x1-6,y2,'a','fontname',fnt,'fontsize',fntSz+2,...
        'fontweight','bold','verticalalign','bottom')
end
axis square
pos = get(plt,'pos');
set(plt,'pos',[pos(1)-0.05,pos(2),pos(3),pos(4)]);
pos = get(gca,'pos');
axes('parent',gcf,'pos',[pos(1)+0.09,pos(2)+0.11,pos(3)-0.08,pos(4)-0.025]);
hold on
set(gca,'layer','top','box','off','xcolor','none','ycolor','none','color','none')
c = colorbar('location','north');
colormap gray
xlabel(c,'Quantile','fontname',fnt,'fontsize',fntSz-2)
axis square

% Panel B
idx = exSub;
x1 = 0; x2 = 1;
y1 = 0; y2 = 1;
plt = subplot(2,4,5);
hold on
plot(xaxis,squeeze(mean(laterCDFs(idx,1,quant2),1)),'color',cols3(2,:),'linewidth',1.5)
plot(xaxis,squeeze(mean(laterCDFs(idx,2,quant2),1)),'color',cols3(3,:),'linewidth',1.5)
plot(xaxis,mean(race(idx,quant2),1),'r','linewidth',1.5)
for i = 1:3
    scatter(xaxis,squeeze(mean(CDFs(idx,i,1,:),1)),15,cols3(i,:),'filled')
end
set(gca,'box','off','tickdir','out')
xlim([x1,x2]);
ylim([y1,y2]);
title('Otto''s Model','fontname',fnt,'fontsize',fntSz)
xlabel('Quantile','fontname',fnt,'fontsize',fntSz)
ylabel({'Cumulative';'probability'},'fontname',fnt,'fontsize',fntSz)
if fig_label
    text(x1-0.4,y2,'b','fontname',fnt,'fontsize',fntSz+2,...
        'fontweight','bold','verticalalign','bottom')
end
text(0.6,0.65,['{\rho} = ',num2str(round(rho(exSub)*100)/100)],'fontname',fnt,'fontsize',fntSz)
text(0.6,0.55,['{\eta} = ',num2str(round(eta(exSub)*100)/100)],'fontname',fnt,'fontsize',fntSz)
l = legend('A_f_i_t','V_f_i_t','Model');
pos = get(l,'pos');
set(l,'box','off','location','southeast','fontname',fnt,...
    'fontsize',fntSz-2,'pos',[pos(1)-0.04,pos(2)-0.15,pos(3),pos(4)])
axis square
pos = get(plt,'pos');
set(plt,'pos',[pos(1)-0.05,pos(2),pos(3),pos(4)]);

% Panel C
x1 = 0; x2 = 12;
y1 = -0.4; y2 = 0.8;
plt = subplot(2,4,2:3);
hold all
scatter([1e3,2e3],[1e3,2e3],[],cols1(1,:),'filled');
scatter([1e3,2e3],[1e3,2e3],[],cols1(2,:),'filled');
plot([0,12],[0,0],':k')
boxdotplot([1,2,4,5,7,8,10,11],eta_cell,col,'iqr',nperm,[0.5,0.5,0.5],0.5,0.5,15,0.2,2);
xlim([x1,x2]);
ylim([y1,y2]);
set(gca,'box','off','tickdir','out','xtick',[1.5,4.5,7.5,10.5],...
    'xticklabel',age_groups,'ytick',-0.4:0.4:0.8)
title('Best-fitting Model Parameters','fontname',fnt,'fontsize',fntSz)
ylabel('Noise (\eta)','fontname',fnt,'fontsize',fntSz)
if fig_label
    text(x1-2,y2,'c','fontname',fnt,'fontsize',fntSz+2,...
        'fontweight','bold','verticalalign','bottom')
end
pos = get(plt,'pos');
set(plt,'pos',[pos(1),pos(2)+0.03,pos(3),pos(4)-0.06]);
for i = 1:2
    text(0.25,y2-0.1*(y2-y1)*i,diag_groups{i},'color',cols1(i,:),'fontweight','bold');
end

% Panel D
x1 = 0; x2 = 12;
y1 = -1; y2 = 1;
plt = subplot(2,4,6:7);
hold all
plot([0,12],[0,0],':k')
boxdotplot([1,2,4,5,7,8,10,11],rho_cell,col,'iqr',nperm,[0.5,0.5,0.5],0.5,0.5,15,0.2,2);
text(1.5,1.1,'*','fontsize',fntSz+8,'horizontalalign','center')
text(4.5,1.1,'**','fontsize',fntSz+8,'horizontalalign','center')
text(7.5,1.1,'**','fontsize',fntSz+8,'horizontalalign','center')
xlim([x1,x2]);
ylim([y1,y2]);
set(gca,'box','off','tickdir','out','xtick',[1.5,4.5,7.5,10.5],'xticklabel',age_groups)
ylabel('Correlation (\rho)','fontname',fnt,'fontsize',fntSz)
if fig_label
    text(x1-2,y2,'d','fontname',fnt,'fontsize',fntSz+2,...
        'fontweight','bold','verticalalign','bottom')
end
xlabel('Age (yrs)','fontname',fnt,'fontsize',fntSz)
pos = get(plt,'pos');
set(plt,'pos',[pos(1),pos(2)+0.03,pos(3),pos(4)-0.06]);

ct_meth = 'mad1';

% Compute RT variability
RTvar = zeros(nSubjs,3,4);
for i = 1:nSubjs
    for j = 1:3
        for k = 1:nTasks
            if strcmp(ct_meth,'std')
                RTvar(i,j,k) = std(RTs{i,j,k});
            elseif strcmp(ct_meth,'iqr')
                RTvar(i,j,k) = iqr(RTs{i,j,k});
            elseif strcmp(ct_meth,'mad0')
                RTvar(i,j,k) = mad(RTs{i,j,k},0);
            elseif strcmp(ct_meth,'mad1')
                RTvar(i,j,k) = mad(RTs{i,j,k},1);
            elseif strcmp(ct_meth,'range')
                RTvar(i,j,k) = range(RTs{i,j,k});
            end
        end
    end
end

% Preallocate memory
idx = cell(2,4);
RTavg = zeros(2,4,3);
err = zeros(2,2,4,3);
ci = zeros(2,2,4,3);

ct_meth = 'median';

% Compute group RT central tendancy and 95% CIs
for i = 1:2
    for j = 1:4
        idx{i,j} = find(grp_vec==i & dev_vec==j);
%         for k = 1:3
%             mc_normalitytest(RTct(idx{i,j},k,1))
%         end
        if strcmp(ct_meth,'mean')
            RTavg(i,j,:) = squeeze(mean(RTvar(idx{i,j},:,1)));
            ci(:,i,j,:) = bootci(nperm,@mean,squeeze(RTvar(idx{i,j},:,1)));
        elseif strcmp(ct_meth,'median')
            RTavg(i,j,:) = squeeze(median(RTvar(idx{i,j},:,1)));
            ci(:,i,j,:) = bootci(nperm,@median,squeeze(RTvar(idx{i,j},:,1)));
        elseif strcmp(ct_meth,'mode')
            RTavg(i,j,:) = squeeze(mode(round(RTvar(idx{i,j},:,1))));
            ci(:,i,j,:) = bootci(nperm,@mode,squeeze(round(RTvar(idx{i,j},:,1))));
        end
    end
end

% Compute error
err(1,:,:,:) = RTavg-squeeze(ci(1,:,:,:));
err(2,:,:,:) = squeeze(ci(2,:,:,:))-RTavg;

% Panel E
x1 = -1; x2 = 1;
y1 = -0.11; y2 = 0.4;
lnWth = 1.5;
for i = 1:2
    idx = grp_vec==i;
    plt = subplot(2,4,i*4);
    hold on
    for j = 1:3
        
        scatter(rho(idx),cost(idx,j),10,cols3(j,:),'filled')
    
        switch rVal
            case 'r'
                [r,stats] = permtest_corr(rho(idx),cost(idx,j),...
                    'nperm',nperm,'type','pearson','tail','both')
                if stats.p<0.05, mkr = '*'; else mkr = ''; end
                text(x2,y2-((j-1)*(y2-y1)*0.1),['{\itr} = ',num2str(round(r*1e2)/1e2),mkr],...
                    'horizontalalign','right','verticalalign','top',...
                    'color',cols3(j,:))
            case 'rho'
                [r,stats] = permtest_corr(rho(idx),cost(idx,j),...
                    'nperm',nperm,'type','spearman','tail','both');
                if stats.p<0.05, mkr = '*'; else mkr = ''; end
                text(x2,y2-((j-1)*(y2-y1)*0.1),['{\itr}_s = ',num2str(round(r*1e2)/1e2),mkr],...
                    'horizontalalign','right','verticalalign','top',...
                    'color',cols3(j,:))
            case 'R2'
                [r,stats] = permtest_corr(rho(idx),cost(idx,j),...
                    'nperm',nperm,'type','pearson','tail','both');
                if stats.p<0.05, mkr = '*'; else mkr = ''; end
                text(x2,y2-((j-1)*(y2-y1)*0.1),['{\itR}^2 = ',num2str(round(r^2*1e2)/1e2),mkr],...
                    'horizontalalign','right','verticalalign','top',...
                    'color',cols3(j,:))
        end
    
    end
    
    set(gca,'box','off','tickdir','out','ytick',-0.1:0.1:0.4,'xtick',-1:0.5:1,...
        'fontname',fnt,'fontsize',fntSz)
    pos = get(plt,'pos');
    set(plt,'pos',[pos(1)+0.03,pos(2),pos(3),pos(4)]);
    if i==1
        if fig_label
            text(x1-(x2-x1)*0.4,y2,'e','fontname',fnt,'fontsize',fntSz+2,...
                'fontweight','bold','verticalalign','bottom')
        end
        for j = 1:3
            text(x1+(x2-x1)*0.025,y2-(j-1)*(y2-y1)*0.1,conditions{j},'color',cols3(j,:),...
                'fontweight','bold','horizontalalign','left',...
                'verticalalign','top');
        end
    else
        xlabel('Correlation (\rho)','fontname',fnt,'fontsize',fntSz)
        ylabel('MSE (AUC)','fontname',fnt,'fontsize',fntSz)
    end
    title(diag_groups{i},'fontname',fnt,'fontsize',fntSz)
    xlim([x1,x2])
    ylim([y1,y2])
    axis square
    
end


set(gcf,'color','w')

%% RESULTS SECITON

% 8. Linking MSE to RMV

%% 8.1 Partial correlaiton analysis: MSE vs Gain

clc;
format long

for i = 1:2
    
    idx = grp_vec==i;    
    gainSw = trapz(prob,squeeze(CDFs(idx,7,3,:)),2);
    gainRe = trapz(prob,squeeze(CDFs(idx,7,2,:)),2);

    [r,p] = partialcorr(gainSw,cost(idx,:),age_vec(idx))
    R2 = r.^2
    [r,p] = partialcorr(gainRe,cost(idx,:),age_vec(idx))
    R2 = r.^2

end

format short

%% 8.1 Partial correlaiton analysis: ADOS vs Gain

clc;
format long

for i = 2
    
    idx = grp_vec==i & ~isnan(ados_vec);
%     idx = grp_vec==i & ~isnan(ados_vec) & ados_vec>0;

%     [r,p] = partialcorr(gain(idx),ados_vec(idx),age_vec(idx))
%     R2 = r.^2
%     [r,stats] = permtest_corr(gain(idx),ados_vec(idx))
%     R2 = r.^2
    
    [r,p] = partialcorr(cost(idx,2:3),ados_vec(idx),age_vec(idx))
    R2 = r.^2
    [r,stats] = permtest_corr(cost(idx,2:3),ados_vec(idx))
    R2 = r.^2

end

format short


%% Fig. 8A,B: Linking MSE to RMV

close all;
clc;

xaxis = prob;
nperm = 1e4;
fig_label = true;
lnWth = 1.5;

axes('parent',figure,'fontname',fnt);
set(gcf,'pos',[20,50,1000,650])

% Panel A

x1 = 0; x2 = 1;
y1 = -0.2; y2 = 0.2;

ctr = 2;
for i = 1:2
    plt = subplot(5,4,ctr);
    hold on
    for j = 1:4
        idx = grp_vec==i & dev_vec==j;
        rmvi = squeeze(mean(CDFs(idx,7,2,:)));
        plot(xaxis,rmvi,'linewidth',lnWth,'color',[shades(j),shades(j),0])
    end
    set(gca,'layer','top','box','off','xaxislocation','origin')
    set(gca,'tickdir','both','xtick',0:0.25:1,'xticklabel',[],'ytick',-0.2:0.1:0.2)
    set(gca,'yticklabel',[])
    axis square
    xlim([x1,x2])
    ylim([y1,y2])
    pos = get(plt,'pos');
    if i==1
        title('Repeat Trials','fontname',fnt,'fontsize',fntSz)
        set(plt,'pos',[pos(1),pos(2),pos(3),pos(4)])
    else
        set(plt,'pos',[pos(1),pos(2)+0.035,pos(3),pos(4)])
    end
    ctr = ctr+4;
end

% Plot x-axis
pos = get(gca,'pos');
axes('parent',gcf,'pos',[pos(1),pos(2)-0.025,pos(3),pos(4)]);
hold on
set(gca,'layer','top','box','off','tickdir','both',...
    'xtick',0:0.25:1,'xticklabel',{'0','','0.5','','1'},'ycolor','none','color','none')
xlabel('Quantile','fontname',fnt,'fontsize',fntSz)
axis square

ctr = 1;
for i = 1:2
    plt = subplot(5,4,ctr);
    hold on
    for j = 1:4
        idx = grp_vec==i & dev_vec==j;
        rmvi = squeeze(mean(CDFs(idx,7,3,:)));
        plot(xaxis,rmvi,'linewidth',lnWth,'color',[shades(j),shades(j),0])
    end
    set(gca,'layer','top','box','off','xaxislocation','origin')
    set(gca,'tickdir','both','xtick',0:0.25:1,'xticklabel',[],'ytick',-0.2:0.1:0.2)
    set(gca,'yticklabel',{'','-0.1','0','0.1',''})
    axis square
    xlim([x1,x2])
    ylim([y1,y2])
    ylabel({diag_groups{i};'\DeltaCP'},'fontname',fnt,'fontsize',fntSz)
    pos = get(plt,'pos');
    if i==1
        title('Switch Trials','fontname',fnt,'fontsize',fntSz)
        if fig_label
            text(x1-(x2-x1)*0.6,y2,'a','fontname',fnt,'fontsize',fntSz+2,...
                'fontweight','bold','verticalalign','bottom')
        end
        set(plt,'pos',[pos(1)+0.09,pos(2),pos(3),pos(4)])
    else
        set(plt,'pos',[pos(1)+0.09,pos(2)+0.035,pos(3),pos(4)])
        l = legend(age_groups);
        pos = get(l,'pos');
        set(l,'box','off','fontsize',fntSz-2,...
            'pos',[pos(1)+0.02,pos(2)-0.095,pos(3),pos(4)])
    end
    ctr = ctr+4;
end

% Panel B
x1 = -0.2; x2 = 0.2;
y1 = -0.2; y2 = 0.2;
for i = 1:1
    plt = subplot(3,4,i+2);
    idx = grp_vec==i;
    gainRe = trapz(prob,squeeze(CDFs(idx,7,2,:)),2);
    gainSw = trapz(prob,squeeze(CDFs(idx,7,3,:)),2);
    [tstat,stats,~,param] = permtest_t(gainSw,gainRe,'nperm',nperm);
    d = mes(gainSw,gainRe,'hedgesg','isDep',1,'nBoot',nperm);
    bf01 = 1/bf.ttest(gainSw,gainRe);
    fprintf('t(%d) = %0.2f, ',param.df,tstat)
    fprintf('p = %0.4f, ',stats.p)
    fprintf('g = %0.2f, ',d.hedgesg)
    fprintf('95CI [%0.2f, %0.2f], ',d.hedgesgCi)
    fprintf('BF01 = %0.15f\n',bf01)
    fprintf('Percentage: %0.1f%%\n',sum(gainSw>gainRe)/numel(gainSw)*100)
    
    hold all
    scatter(gainRe,gainSw,10,cols1(i,:),'filled')
    plot([0,0],[y1,y2],':k')
    plot([x1,x2],[0,0],':k')
    plot([x1,x2],[y1,y2],'--k')
    xlim([x1,x2])
    ylim([y1,y2])
    axis square
    set(gca,'tickdir','out','xtick',-0.2:0.1:0.2,'ytick',-0.2:0.1:0.2,...
        'fontname',fnt,'fontsize',fntSz)
    title(diag_groups{i},'fontname',fnt,'fontsize',fntSz)
    
    if i==1
        xlabel('Repeat trials','fontname',fnt,'fontsize',fntSz)
        ylabel('Switch trials','fontname',fnt,'fontsize',fntSz)
        text(x2-0.05,y2+0.08,'Multisensory Gain','fontname',fnt,'fontsize',fntSz)
        if fig_label
            text(x1-(x2-x1)*0.35,y2,'b','fontname',fnt,'fontsize',fntSz+2,...
                'fontweight','bold','verticalalign','bottom')
        end
    else
        set(gca,'yticklabel',[])
    end
    pos = get(plt,'pos');
    set(plt,'pos',[pos(1),pos(2)-0.025,pos(3),pos(4)])
end

% % plt = subplot(3,4,7:9);
% plt = subplot(3,4,8:12);
% imshow('C:\Users\mcrosse\Pictures\Figures\Projects\AVSRT\mediation_model.png')

set(gcf,'color','w')

%% Fig. 8C,D: Mediation analysis

clc;

% Set parameters
nperm = 1e4;

for i = 2
    
    idx = grp_vec==i & dev_vec~=4;
%     idx = grp_vec==i;
    
    % Define race model violation
    X = zscore(age_vec(idx));
    Y = zscore(gain(idx));
    M = zscore(mean(cost(idx,2:3),2));
            
    [~,stats] = mediation(X,Y,M,'boot','bootsamples',nperm);
    disp(stats)
    
    [r,p] = partialcorr(gain(idx),cost(idx,:),age_vec(idx))

end

%% Fig. 8E,F: Correlation matrix

close all;
clc;

% Define variables
vars1 = {'Age','PIQ','VIQ','FSIQ','ADOS','F_1','MSE_A','MSE_V',...
    'Ben_P_r_e_d','Ben_E_m_p','Gain','\it\rho','\it\eta'};
vars2 = {{'Age';'\downarrow'},{'PIQ';'\downarrow'},{'VIQ';'\downarrow'},...
    {'FSIQ';'\downarrow'},{'ADOS';'\downarrow'},{'F_1';'\downarrow'},...
    {'MSE_A_V';'\downarrow'},{'MSE_A';'\downarrow'},{'MSE_V';'\downarrow'},...
    {'Gain';'\downarrow'},{'{\it\rho}';'\downarrow'},{'{\it\eta}';'\downarrow'}};
mseAV = squeeze(trapz(prob,mean(CDFs(:,1,2,:)-CDFs(:,1,3,:),2),4));
mseA = squeeze(trapz(prob,mean(CDFs(:,2,2,:)-CDFs(:,2,3,:),2),4));
mseV = squeeze(trapz(prob,mean(CDFs(:,3,2,:)-CDFs(:,3,3,:),2),4));

% Define parameters
nVars = numel(vars1);
nperm = 1e4;
type = 'pearson';
correct_fwer = false;

x1 = 0.5; x2 = 13.5;
y1 = 0.5; y2 = 13.5;

for n = 1:numel(diag_groups)
    
    % Index group
    idx = grp_vec==n;
    
    % Preallocate memory
    data = zeros(sum(idx),nVars);
    
    % Construct data matrix
    data(:,1) = age_vec(idx);
    data(:,2) = piq_vec(idx);
    data(:,3) = viq_vec(idx);
    data(:,4) = fsiq_vec(idx);
    data(:,5) = ados_vec(idx);
    data(:,6) = Fscore(idx);
    data(:,7) = mseA(idx);
    data(:,8) = mseV(idx);
    data(:,9) = benPred(idx);
    data(:,10) = benEmp(idx);
    data(:,11) = gain(idx);
    data(:,12) = racefit(idx,1);
    data(:,13) = racefit(idx,2);
    
    % Compute ADOS correlations
    [r1,stats1,orig1] = permtest_corr(data,[],'nperm',nperm,'rows','complete','type',type);
    
    % Compute IQ correlations w/ ADOS scores set to zero
    data(isnan(data(:,5)),5) = 0;
    [r2,stats2,orig2] = permtest_corr(data,[],'nperm',nperm,'rows','complete','type',type);
    
    % Compute all other correlations
    [r,stats,orig] = permtest_corr(data,[],'nperm',nperm,'rows','all','type',type);
    
    if correct_fwer
        p = stats.p;
        p(isnan(r)) = stats2.p(isnan(r)); r(isnan(r)) = r2(isnan(r));
        p(:,5) = stats1.p(:,5); r(:,5) = r1(:,5);
        p(5,:) = stats1.p(5,:); r(5,:) = r1(5,:);
    else
        p = orig.p;
        p(isnan(r)) = orig2.p(isnan(r)); r(isnan(r)) = r2(isnan(r));
        p(:,5) = orig1.p(:,5); r(:,5) = r1(:,5);
        p(5,:) = orig1.p(5,:); r(5,:) = r1(5,:);
    end
        
    % Extract lower triangle
    triMat = tril(r);
    triMat2 = tril(r.^2);
    triMat2(triMat2==0) = 1;
    triMat2(isnan(triMat2)) = 1;
    
    % Create figure
    figure
    hold all
    imagesc(triMat2)
    plot([0,5.5],[5.5,5.5],':r','linewidth',2)
    plot([5.5,5.5],[5.5,13.5],':r','linewidth',2)
    plot([0,11.5],[11.5,11.5],':r','linewidth',2)
    plot([11.5,11.5],[11.5,13.5],':r','linewidth',2)
    plot([5.5,6],[5.5,5],':r','linewidth',2)
    plot([11.5,12],[11.5,11],':r','linewidth',2)
    
    for i = 1:size(data,2)
        for j = 1:size(data,2)
            if i==j && j~=nVars
                text(i-0.3,j,vars1{i},'horizontalalign','left','fontname',fnt,'fontsize',fntSz)
            end
            if triMat(i,j)==0
                if triMat2(j,i)>0.6
                    fntCol = 'k';
                else
                    fntCol = 'w';
                end
                if p(i,j)<0.05
                    if isnan(r(i,j))
                        text(i,j,'-','horizontalalign','center','fontname',fnt,'fontsize',fntSz-1,'color',fntCol)
                    else
                        text(i,j,[num2str(round(r(i,j)*100)/100),'*'],'horizontalalign','center','fontname',fnt,'fontsize',fntSz-1,'color',fntCol)
                    end
                else
                    text(i,j,num2str(round(r(i,j)*100)/100),'horizontalalign','center','fontname',fnt,'fontsize',fntSz-1,'color',fntCol)
                end
            end
        end
    end
    text(x2,8,'*{\itp} < 0.05','fontname',fnt,'fontsize',fntSz,...
        'horizontalalign','right')
    set(gca,'ydir','reverse','xtick',[3,8.5,12.5],'xticklabels',{'Participant measures','Empirical measures','Model params'},'ytick',2:size(data,2),'yticklabels',vars1(2:end),'fontname',fnt,'fontsize',fntSz)
    title(diag_groups{n},'fontname',fnt,'fontsize',fntSz)
    xlim([x1,x2])
    ylim([y1,y2])
    colormap bone
    pos = get(gca,'pos');
    axes('parent',gcf,'pos',[pos(1)+0.26,pos(2)+0.375,pos(3),pos(4)-0.4]);
    hold on
    set(gca,'layer','top','box','off','xcolor','none','ycolor','none','color','none')
    
    % Colorbar
    c = colorbar('location','east');
    title(c,'{\itR}^2','fontname',fnt,'fontsize',fntSz)
    axis square
    
    set(gcf,'color','w')
    
end

%% 1.2 Multisensory benefits (LM)
% A linear regression analysis was used to examine the effects of diagnosis
% and age on multisensory benefit.

close all;
clc;

% Define variables
subj = 1:nSubjs;
grp = grp_vec;
sex = sex_vec;
age = age_vec;
dev = dev_vec;
piq = piq_vec;
viq = viq_vec;
fsiq = fsiq_vec;
ados = ados_vec;
mseAV = cost(:,1);
mseA = cost(:,2);
mseV = cost(:,3);

% Convert non-numeric variables to categorical
subj = categorical(subj','Ordinal',0);
grp = categorical(grp,[1,2],{'NT','ASD'},'Ordinal',0);
sex = categorical(sex,[1,2],{'Male','Fem'},'Ordinal',0);
dev = categorical(dev,[1,2,3,4],{'6-9','10-12','13-17','18-40'},'Ordinal',1);

% Create dataset
ds = dataset(subj,grp,sex,age,dev,piq,viq,fsiq,ados,gain,mseAV,mseA,mseV,...
    benPred,benEmp,eta,rho);

% Fit linear model to dataset
% lm = fitlm(ds,'viq ~ gain + grp');
lm = fitlm(ds,'ados ~ mseA * age');

% Perform hypothesis tests on fixed effect terms
stats = anova(lm);

% Display LME and ANOVA stats
disp(lm)
fprintf('\nANOVA on fixed effect terms:\n')
disp(stats)

% Define residuals
r = lm.Residuals.Raw;
pr = lm.Residuals.Pearson;
st = lm.Residuals.Standardized;

% Check normality of residuals
disp('Normality test on residuals:')
mc_normalitytest(r);

figure,
subplot(2,2,1), histfit(r), xlabel('Residuals'), ylabel('Frequency')
title('Histogram of residuals with fitted normal density')
subplot(2,2,2), qqplot(r)
subplot(2,2,3), plotResiduals(lm,'fitted')
subplot(2,2,4), plotResiduals(lm,'lagged')

figure, boxplot([r,pr,st])

