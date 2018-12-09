function MS2norm2 = MS2normalize(MS2norm)

for i=1:size(MS2norm,2)
    temp = MS2norm{2,i};
    temp(temp<-7000) = 0;
    MS2norm{2,i} = temp;
end


mVal = mean(mean(MS2norm{2,1}(MS2norm{2,1}>0)));

for i = 2:size(MS2norm,2)
    currM = mean(mean(MS2norm{2,i}(MS2norm{2,i}>0)));
    MS2norm{2,i} = MS2norm{2,i} * mVal / currM;
end
    

for i=1:size(MS2norm,2)
    for j=1:size(MS2norm{1,i},2)
        if j ==1
            MS2norm{1,i}(3,j) = 1;
        else
            if MS2norm{1,i}(2,j-1) == MS2norm{1,i}(2,j) && MS2norm{1,i}(1,j-1) == MS2norm{1,i}(1,j)
                MS2norm{1,i}(3,j) = MS2norm{1,i}(3,j-1)+1;
            else
                MS2norm{1,i}(3,j) = 1;
            end
        end
    end
end
    


%%% put 0s (rest) & 1s (burst) for the data in 2nd row of MS2norm.
TimeInterval = 5    ;    %%% Time interval in minutes

for i=1:size(MS2norm,2)
    BurstRest = zeros(size(MS2norm{2,i},1), size(MS2norm{2,i},2));
    durationTP = floor(MS2norm{3,i}/TimeInterval);
    for j=1:size(durationTP,2)
        allSize = sum(durationTP(:,j));
        zeroOne = zeros(allSize,1);
        loc = 1;
        if durationTP(1,j) ~= 0
            preLoc = loc;
            loc = loc+durationTP(1,j);
            zeroOne(preLoc:loc-1,1) = 0;
        end
        
        for k=2:size(durationTP,1)
            if durationTP(k,j) == 0
                break
            end
            preLoc = loc;
            loc = loc+durationTP(k,j);
            
            if mod(k,2) == 0
                zeroOne(preLoc:loc-1,1) = 1;
            else
                zeroOne(preLoc:loc-1,1) = 0;
            end
        end
        
        if allSize < size(BurstRest,1)
            BurstRest(1:allSize,j) = zeroOne;
        elseif allSize > size(BurstRest,1)
            BurstRest(:,j) = zeroOne(1:size(BurstRest,1),1);
        else
            BurstRest(:,j) = zeroOne;
        end
    end
    
    MS2norm{5,i} = BurstRest;
end
     
MS2norm2 = MS2norm;