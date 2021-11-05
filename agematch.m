function [idx1,idx2] = agematch(devi,grp,dev,sex,age,iq,plotknn)


% Get control/patient data
if isempty(iq)
    idx1 = find(grp==1 & dev==devi);
    idx2 = find(grp==2 & dev==devi);
    controls = [sex(idx1),age(idx1)];
    patients = [sex(idx2),age(idx2)];
else
    idx1 = find(grp==1 & dev==devi & ~isnan(iq));
    idx2 = find(grp==2 & dev==devi & ~isnan(iq));
    controls = [sex(idx1),age(idx1),iq(idx1)];
    patients = [sex(idx2),age(idx2),iq(idx2)];
end

% Perform KNN search for closest match
k1 = numel(idx1);
k2 = numel(idx2);
if k1>=k2
    [idx,d] = knnsearch(controls,patients,'k',k1,'distance','seuclidean');
elseif k1<k2
    [idx,d] = knnsearch(patients,controls,'k',k2,'distance','seuclidean');
end

% Redistribute duplicates
while numel(unique(idx(:,1)))<size(idx,1)
    for j = 1:size(idx,1)
       ties = find(idx(:,1)==idx(j,1));
       if numel(ties)>1 && d(j,1)~=min(d(ties,1))
           idx(j,:) = circshift(idx(j,:),-1);
           d(j,:) = circshift(d(j,:),-1);
       elseif numel(ties)>1
           ties(ties==j) = [];
           for k = 1:numel(ties)
               idx(ties(k),:) = circshift(idx(ties(k),:),-1);
               d(ties(k),:) = circshift(d(ties(k),:),-1);
           end
       end
    end
end

% Return index of controls
if k1>=k2
    idx1 = idx1(idx(:,1));
elseif k1<k2
    idx2 = idx2(idx(:,1));
end

if plotknn
    
    if k1>=k2
        age1 = controls(idx(:,1),1);
        age2 = age(idx2);
    elseif k1<k2
        age1 = age(idx1);
        age2 = patients(idx(:,1),1);
    end
    
    if ~isempty(iq)
        if k1>=k2
            piq1 = controls(idx(:,1),2);
            piq2 = iq(idx2);
        elseif k1<k2
            piq1 = iq(idx1);
            piq2 = patients(idx(:,1),2);
        end
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
    title('KNN Search')
    
end
