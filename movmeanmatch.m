function [idx1,idx2] = movmeanmatch(i,grp,sex,age,piq,ageRng,window,plotflag)


% Get relevant controls and patients features
if isempty(piq)
    idx1 = find(grp==1 & age>ageRng(i)-window/2 & age<ageRng(i)+window/2);
    idx2 = find(grp==2 & age>ageRng(i)-window/2 & age<ageRng(i)+window/2);
    controls = [sex(idx1),age(idx1)];
    patients = [sex(idx2),age(idx2)];
else
    idx1 = find(grp==1 & age>ageRng(i)-window/2 & age<ageRng(i)+window/2 & ~isnan(piq));
    idx2 = find(grp==2 & age>ageRng(i)-window/2 & age<ageRng(i)+window/2 & ~isnan(piq));
    controls = [sex(idx1),age(idx1),piq(idx1)];
    patients = [sex(idx2),age(idx2),piq(idx2)];
end

% Perform KNN search for closest controls
k = sum(idx1);
[idx,d] = knnsearch(controls,patients,'k',k,'distance','seuclidean');

% Redistribute duplicates
for n = 1:1e3
    if length(unique(idx(:,1)))==length(idx(:,1))
        disp('Yurt!')
        break
    end
    for j = 1:size(idx,1)
       ties = find(idx(:,1)==idx(j,1));
       if length(ties)>1 && d(j,1)~=min(d(ties,1))
           idx(j,:) = circshift(idx(j,:),-1);
           d(j,:) = circshift(d(j,:),-1);
       end
    end
end

% Return index of controls
idx1 = idx1(idx(:,1));


if plotflag
    
    age1 = controls(idx(:,1),2);
    age2 = age(idx2);
    
    if ~isempty(piq)
        piq1 = controls(idx(:,1),3);
        piq2 = piq(idx2);
    else
        piq1 = zeros(size(age1));
        piq2 = ones(size(age2));
    end

    figure, hold all
    scatter(age1,piq1,50,'b','filled')
    scatter(age2,piq2,50,'r')
    for j = 1:size(idx,1)
        plot([age1(j),age2(j)],[piq1(j),piq2(j)],'k')
    end
    xlabel('Age (yrs)')
    ylabel('PIQ')
    title('KNN Output')
    
end



