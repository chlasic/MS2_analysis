function [MS2sumIntData, perATSsync, minATSsync, ATS1ATS2int, ATS1ATS2intRand] = MS2sumInt(MS2norm)


%%% Summed ATS intensity over time
%%% A new variable 'MS2sumIntData' gets generated.
%%% The variable has 4 rows in a cell format | 1: nuc info | 2: summed ATS intensity over time|.
%%% In the 1st row in cell,
%%% | 1: Image ID | 2: Nucleus ID | 3: Location (um from distal end) | 4: Total # ATS/nuc.
%%% In the 2nd row,
%%% | Normalized summed ATS signal intensity in each point |
%%% In the 3rd row,
%%% | Total # of ATS/nuc seen at each time point |
%%% In the 4th,
%%% | Pearson's correlation coefficient(s) between ATS in the same nuc|.
%%% In the 5th,
%%% | Pearson's r between ATS in the same nuc only when both ATS are ON |.
%%% In the 6th,
%%% | Person's r with random data.
%%% In the 7th,
%%% | Person's r with random data (the data pooled for each nucleus).


%%% The variable 'perATSsync' contains 8 rows, each row records:
%%% 1st row: % both ATS are ON,
%%% 2nd:     % both ATS are OFF,
%%% 3rd:     % ATS are synchronized (ON or OFF),
%%% 4th:     % ATS are async,
%%% 5-8th: Repetition of 1-4th row but with probability calculated from original observation.






MS2sumIntData = cell(7,size(MS2norm,2));
perATSsync = zeros(4,999999);
minATSsync = zeros(4,999999);
ATS1ATS2int = cell(1,size(MS2norm,2));
ATS1ATS2intRand = cell(1,size(MS2norm,2));
loc2 = 1;


for i=1:size(MS2norm,2)  
    loc = 1;
    locATS12 = 1;
    for j=1:MS2norm{1,i}(1,end)
        indx = find(MS2norm{1,i}(1,:) == j);
        if ~isempty(indx)
            MS2sumIntData{1,i}(1,loc) = i;
            MS2sumIntData{1,i}(2,loc) = j;
            MS2sumIntData{1,i}(3,loc) = MS2norm{1,i}(2,indx(1));
            MS2sumIntData{1,i}(4,loc) = length(indx);

            MS2sumIntData{2,i}(:,loc) = sum(MS2norm{2,i}(:,indx),2);
            MS2sumIntData{3,i}(:,loc) = sum(MS2norm{5,i}(:,indx(:)),2);
            
            IS = length(indx);
            
            %%% Correlation betw. ATS in a nuc (using all time points)
            if IS == 2
                CorrT = zeros(1,1);
                CorrT(1) = corr(MS2norm{2,i}(:,indx(1)), MS2norm{2,i}(:,indx(2)));
                MS2sumIntData{4,i}(1,loc) = CorrT;
            elseif IS == 3
                CorrT = zeros(3,1);
                CorrT(1) = corr(MS2norm{2,i}(:,indx(1)), MS2norm{2,i}(:,indx(2)));
                CorrT(2) = corr(MS2norm{2,i}(:,indx(3)), MS2norm{2,i}(:,indx(2)));
                CorrT(3) = corr(MS2norm{2,i}(:,indx(1)), MS2norm{2,i}(:,indx(3)));
                MS2sumIntData{4,i}(1:3,loc) = CorrT; 
            elseif IS == 4
                CorrT = zeros(6,1);
                CorrT(1) = corr(MS2norm{2,i}(:,indx(1)), MS2norm{2,i}(:,indx(2)));
                CorrT(2) = corr(MS2norm{2,i}(:,indx(1)), MS2norm{2,i}(:,indx(3)));
                CorrT(3) = corr(MS2norm{2,i}(:,indx(1)), MS2norm{2,i}(:,indx(4)));
                CorrT(4) = corr(MS2norm{2,i}(:,indx(2)), MS2norm{2,i}(:,indx(3)));
                CorrT(5) = corr(MS2norm{2,i}(:,indx(2)), MS2norm{2,i}(:,indx(4)));
                CorrT(6) = corr(MS2norm{2,i}(:,indx(4)), MS2norm{2,i}(:,indx(3)));
                MS2sumIntData{4,i}(1:6,loc) = CorrT; 
            end
            
            %%% Correlation of ATSs in a nuc ONLY WHEN TRX FIRING.
            if IS == 2
                CorrT = zeros(1,1);
                aVal = 1; bVal = 2;
                tempPos = MS2norm{5,i}(:,indx(aVal)) .* MS2norm{5,i}(:,indx(bVal));
                MS2temp1 = MS2norm{2,i}(tempPos>0,indx(aVal));
                MS2temp2 = MS2norm{2,i}(tempPos>0,indx(bVal));
                if ~isempty(MS2temp1) && ~isempty(MS2temp2)
                    CorrT(1) = corr(MS2temp1, MS2temp2);
                    MS2temp2rand = MS2temp2(randperm(length(MS2temp2)));
                    CorrTrand = corr(MS2temp1, MS2temp2rand);
                else
                    CorrT = 0;
                    CorrTrand = 0;
                end
                MS2sumIntData{5,i}(1,loc) = CorrT;
                MS2sumIntData{7,i}(1,loc) = CorrTrand;
            elseif IS == 3
                CorrT = zeros(3,1);
                MS2temp1=[];
                MS2temp2=[];
                for k=1:3
                    if k == 1
                        aVal = 1; bVal = 2;
                    elseif k == 2
                        aVal = 3; bVal = 2;
                    elseif k == 3
                        aVal = 1; bVal = 3;
                    end
                    tempPos = MS2norm{5,i}(:,indx(aVal)) .* MS2norm{5,i}(:,indx(bVal));
                    MS2temp11 = MS2norm{2,i}(tempPos>0,indx(aVal));
                    MS2temp22 = MS2norm{2,i}(tempPos>0,indx(bVal));
                    if ~isempty(MS2temp11) && ~isempty(MS2temp22)
                        CorrT(k) = corr(MS2temp11, MS2temp22);
                        MS2temp1=[MS2temp1;MS2temp11];
                        MS2temp2=[MS2temp2;MS2temp22];
                    end   
                end
                if ~isempty(MS2temp1) && ~isempty(MS2temp2)
                    CorrTnew = corr(MS2temp1, MS2temp2);
                    MS2temp2rand = MS2temp2(randperm(length(MS2temp2)));
                    CorrTrand = corr(MS2temp1, MS2temp2rand);
                else
                    CorrTnew = 0;
                    CorrTrand = 0;
                end
                MS2sumIntData{5,i}(1:3,loc) = CorrTnew; 
                MS2sumIntData{7,i}(1:3,loc) = CorrTrand;
            elseif IS == 4
                CorrT = zeros(6,1);
                MS2temp1=[];
                MS2temp2=[];
                for k=1:6
                    if k == 1
                        aVal = 1; bVal = 2;
                    elseif k == 2
                        aVal = 1; bVal = 3;
                    elseif k == 3
                        aVal = 1; bVal = 4;
                    elseif k == 4
                        aVal = 2; bVal = 3;
                    elseif k == 5
                        aVal = 2; bVal = 4;
                    elseif k == 6
                        aVal = 3; bVal = 4;
                    end
                    tempPos = MS2norm{5,i}(:,indx(aVal)) .* MS2norm{5,i}(:,indx(bVal));
                    MS2temp11 = MS2norm{2,i}(tempPos>0,indx(aVal));
                    MS2temp22 = MS2norm{2,i}(tempPos>0,indx(bVal));
                    if ~isempty(MS2temp11) && ~isempty(MS2temp22)
                        CorrT(k) = corr(MS2temp11, MS2temp22);
                        MS2temp1=[MS2temp1;MS2temp11];
                        MS2temp2=[MS2temp2;MS2temp22];
                    end
                end
                if ~isempty(MS2temp1) && ~isempty(MS2temp2)
                    CorrTnew = corr(MS2temp1, MS2temp2);
                    MS2temp2rand = MS2temp2(randperm(length(MS2temp2)));
                    CorrTrand = corr(MS2temp1, MS2temp2rand);
                else
                    CorrTnew = 0;
                    CorrTrand = 0;
                end
                MS2sumIntData{5,i}(1:6,loc) = CorrTnew; 
                MS2sumIntData{7,i}(1:6,loc) = CorrTrand;
            end
            
            if IS > 1
                if isempty(MS2temp2)
                    MS2temp2rand = MS2temp2;
                end
                ATS1ATS2int{1,i}(locATS12:locATS12+length(MS2temp1)-1,:) = [MS2temp1 MS2temp2];
                ATS1ATS2intRand{1,i}(locATS12:locATS12+length(MS2temp1)-1,:) = [MS2temp1 MS2temp2rand];
                locATS12 = locATS12+length(MS2temp1);
            end
            
            %%% Correlation with randomly paired ATS (two ATS from different nuclei).
            %%% Data are recorded on 6th row of 'MS2sumIntData'.
            if IS > 1
                combi = nchoosek(1:IS, 2);
%                 combi = [(1:IS)' (1:IS)'];
                CorrT = zeros(size(combi,1),1);
                for k = 1:size(combi,1)
                    randPos = combi(k,1);
                    while randPos == combi(k,1)
                        randPos = randi(size(MS2norm{1,i},2),1);
                    end
                    testSet1 = MS2norm{5,i}(:,indx(combi(k,1)));
                    testSet2 = MS2norm{5,i}(:,randPos);
                    tempPos =  testSet1 .* testSet2;
                    MS2temp1 = MS2norm{2,i}(tempPos>0,indx(combi(k,1)));
                    MS2temp2 = MS2norm{2,i}(tempPos>0,randPos);
%                     MS2temp1 = MS2norm{2,i}(tempPos>-1,indx(combi(k,1)));
%                     MS2temp2 = MS2norm{2,i}(tempPos>-1,randPos);
                    if ~isempty(MS2temp1) && ~isempty(MS2temp2)
                        CorrT(k) = corr(MS2temp1, MS2temp2);
                    end
                end
                MS2sumIntData{6,i}(1:size(combi,1),loc) = CorrT; 
            end
            
            
            %%% --- % of synchronized ATS bursting or resting % probability calculation.
            %%% This part generates 'perATSsync' which has 8 rows.
            if IS > 1
                combi = nchoosek(1:IS, 2);
                for k = 1:size(combi,1)
                    %%% TrxONnum records # of ON-state ATS of the two ATS looked at now.
                    tempPos1 = MS2norm{5,i}(:,indx(combi(k,1)));
                    pPos1 = sum(tempPos1>0)/length(tempPos1)*100;
                    tempPos2 = MS2norm{5,i}(:,indx(combi(k,1)));
                    pPos2 = sum(tempPos2>0)/length(tempPos2)*100;
                    TrxSumNum = MS2norm{5,i}(:,indx(combi(k,1))) + MS2norm{5,i}(:,indx(combi(k,2)));
                    tot = length(TrxSumNum)/100;
                    PerTsum = [sum(TrxSumNum ==2)/tot; sum(TrxSumNum == 0)/tot; ...
                        sum(TrxSumNum ==0 | TrxSumNum ==2)/tot; sum(TrxSumNum ==1)/tot];
                    minTsum = [sum(TrxSumNum ==2); sum(TrxSumNum == 0); ...
                        sum(TrxSumNum ==0 | TrxSumNum ==2); sum(TrxSumNum ==1)];
                    perATSsync(1:4,loc2) = PerTsum;
                    perATSsync(5,loc2) = pPos1 * pPos2 / 100;
                    perATSsync(6,loc2) = (100-pPos1) * (100-pPos2) / 100;
                    perATSsync(7,loc2) = perATSsync(5,loc2)+perATSsync(6,loc2);
                    perATSsync(8,loc2) = 100 - perATSsync(7,loc2);
                    
                    
                    DurON = MS2norm{3,i}(:,indx(combi(k,1)));
                    DurON = DurON(2:2:end);
                    CouON = sum(DurON>0);
                    mpr = 5 / (CouON + 1);
                    
                    minATSsync(1:4,loc2) = minTsum * mpr;
                    minATSsync(5,loc2) = pPos1 * pPos2 / 100 * tot * mpr;
                    minATSsync(6,loc2) = (100-pPos1) * (100-pPos2) / 100 * tot * mpr;
                    minATSsync(7,loc2) = (perATSsync(5,loc2)+perATSsync(6,loc2)) * tot * mpr;
                    minATSsync(8,loc2) = (100 - perATSsync(7,loc2)) * tot * mpr;
                    
                    loc2 = loc2+1;
                end    
            end
            loc=loc+1;
        end
    end
end

perATSsync(:,loc2:end) = [];
minATSsync(:,loc2:end) = [];

