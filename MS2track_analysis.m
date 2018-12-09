%%% Analyze manually tracked MS2 data (sygl-1 ATS over time, preset time gap is 5 min).
%%%
%%% First, manually copy&paste MS2 measurements (MS excel) into 'MS2trk'.
%%% The cell 'MS2trk' contains two rows: 1st = (nuc ID & um distance from distal tip).
%%%                                      2nd = (normalized sig intensity every time points).
%%%                                      3rd = (ATS bursting period including imcomplete events).
%%%                                      4th = (ATS bursting period only for complete events).
%%%                                      5th = 0: rest, 1: burst -- for the 2nd column (MS2 intensity recording).
%%% The MS2trk 3rd column comtains ATS bursting duration (ON or OFF) at the beginning or end that do not show starting/ending time points.
%%% MS2trk serves as a backup variable for MS2norm. Use 'MS2norm' for all analyses.
%%% Raw data (raw signal intensity, bg intensity, etc. are stored at '121417 MS2 manual tracking.xlsx'.


%%% Quick visualization of MS2 signal over time
imgn =    1   ;
atsn =    2   ;

figure('pos',[200 500 500 300]), plot(1:size(MS2trk{2,imgn}(:,atsn),1), MS2trk{2,imgn}(:,atsn));




%% Normalize ATS intensity throughout gonads (mandatory step)
%%% This includes assigning ATS # in MS2norm.
MS2norm = MS2normalize(MS2trk);



%% Apply smooth filter to raw MS2 data (optional)
MS2smooth = cell(size(MS2norm,1),size(MS2norm,2));

for i = 1:size(MS2norm,2)
    for j = 1:size(MS2norm{2,i},2)
        MS2smooth{2,i}(:,j) = smooth(1:size(MS2norm{2,i}(:,j),1),MS2norm{2,i}(:,j),'lowess');
    end
end

MS2smooth(1,:) = MS2norm(1,:);




%% display results
for i = 1:size(MS2trk,2)
    for j = 1:size(MS2trk{2,i},2)

        imgn =    i   ;
        atsn =    j   ;

        figure('pos',[200 300 1500 500])
        hold on
        plot(((1:size(MS2norm{2,imgn}(:,atsn),1))-1)*5, MS2norm{2,imgn}(:,atsn), 'k');
        plot(-100:300, zeros(length(-100:300),1), 'k:');
        
        for k=1:13
            plot([25*k 25*k], [-5000 5000], 'b:');
        end

        axis([ 0  (size(MS2norm{2,imgn}(:,atsn),1)-1)*5   min(MS2norm{2,imgn}(:,atsn))-100  max(MS2norm{2,imgn}(:,atsn))+100  ])
        xticks(0:5:(size(MS2trk{2,imgn}(:,atsn),1)-1)*5)
        xtickangle(45)
        grid on
        distN = num2str(MS2trk{1,i}(2,j));
        title(strcat('distDE =  ',distN, ' um'));
        fprintf('\n\t\t%d/%d ATS in %d/%d gonad.\n\n', j,size(MS2trk{2,i},2) ,i,size(MS2trk,2) ); 
        box on
        
        pause
        close all
    end
end






%% plot mean ATS intensity within gonad (confirming no spacial pattern)
%%% Visualize ATS intensity curve in each gonad
MS2norm = MS2norm_N2;

for i = 1:size(MS2norm,2)
    plot(mean(MS2norm{2,i}));
    pause
        close all
end



%%% scatter plot: distDE (distance from distal end) & ATS intensity
MS2norm = MS2norm_N2;


distATS = MS2norm{1,1}(2,:);
intATS = mean(MS2norm{2,1});

for i = 2:size(MS2norm,2)
    distATS = [distATS MS2norm{1,i}(2,:)];
    intATS = [intATS mean(MS2norm{2,i})];
end

% gets mean of transcriptional activity
gap = 1;
mLine = zeros(1,20);
seL = zeros(1,20);
for i=1:20
    cpool = intATS(distATS > (i-1)*3 & distATS < i*3);
    if isempty(cpool)
        cpool = 0;
    end
    mLine(i) = mean(cpool);
    seL(i) = std(cpool)/sqrt(length(cpool));
end


normz = 0;
figure('pos',[300 200 400 500])
hold off
plot(distATS, intATS+normz, 'k.', 'markersize', 15);     % L4: x2
axis([-1 61  0 4000 ])
hold on
plot(1:3:60, mLine+normz, 'c', 'linewidth', 2);

corn = corr(distATS', intATS')




%% ***** Barplot & Scatterplot: distDE & Trx firing duration *****

%%%%% Put all duration data together in one variable 'TrxONduration'.
MS2norm = MS2norm_N2;

TrxONduration = cell(1,size(MS2norm,2));
TrxOFFduration = cell(1,size(MS2norm,2));

%%% 3: All data points were included (even incomplete ATS burst cycles)
%%% 4: only with full measurement (flanking ATS burst peaks excluded)
UsingSet =  4  ; 

for i = 1:size(MS2norm,2)
    ONcount = 1;
    OFFcount = 1;
    tempStoreON = zeros(4,9999999);
    tempStoreOFF = zeros(4,9999999);
    for j = 1:size(MS2norm{UsingSet,i},2)
        tempON = MS2norm{UsingSet,i}(2:2:end,j);
        tempON(tempON == 0) = [];
        
        if ~isempty(tempON)
            tempStoreON(4,ONcount:ONcount+length(tempON)-1) = tempON;
            for k = 1:3
                tempStoreON(k,ONcount:ONcount+length(tempON)-1) = MS2norm{1,i}(k,j);
            end
            ONcount = ONcount + length(tempON);
        end
        
        tempOFF = MS2norm{UsingSet,i}(1:2:end,j);
        tempOFF(tempOFF == 0) = [];
        
        if ~isempty(tempOFF)
            tempStoreOFF(4,OFFcount:OFFcount+length(tempOFF)-1) = tempOFF;
            for k = 1:3
                tempStoreOFF(k,OFFcount:OFFcount+length(tempOFF)-1) = MS2norm{1,i}(k,j);
            end
            OFFcount = OFFcount + length(tempOFF);
        end
    end
    
    tempStoreON(:,ONcount:end) = [];
    tempStoreOFF(:,OFFcount:end) = [];
        
    TrxONduration{i} = tempStoreON;
    TrxOFFduration{i} = tempStoreOFF;
end
TrxONsum = cell2mat(TrxONduration);
TrxOFFsum = cell2mat(TrxOFFduration);


%%% boxplot
crs=5;
xcr = 0:crs:40;
bins = length(xcr);

ONDurATSpCR = zeros(2, bins);
OFFDurATSpCR = zeros(2, bins);

for j=1:bins   
    xr = [crs*(j-1) crs*j];
    CurrInxON = TrxONsum(2,:) >= xr(1) & TrxONsum(2,:) < xr(2);
    allATSon = sum(CurrInxON);
    
    ONDurATSpCR(1,j) = mean(TrxONsum(4,CurrInxON));
    ONDurATSpCR(2,j) = std(TrxONsum(4,CurrInxON))/sqrt(allATSon);

    CurrInxOFF = TrxOFFsum(2,:) >= xr(1) & TrxOFFsum(2,:) < xr(2);
    allATSoff = sum(CurrInxOFF);
    
    OFFDurATSpCR(1,j) = mean(TrxOFFsum(4,CurrInxOFF));
    OFFDurATSpCR(2,j) = std(TrxOFFsum(4,CurrInxOFF))/sqrt(allATSoff);  
end



%----------bargraph  (ON duration)
figure('pos',[ 800 200 420 500])

hold off
bar(xcr+2.5, ONDurATSpCR(1,:))
hold on
errorbar(xcr+2.5, ONDurATSpCR(1,:), ONDurATSpCR(2,:), 'k.', 'linewidth', 2);
axis([ -1 61 0 100 ])
box on
ylabel('ON duration of \itsygl-1\rm ATS (min)' , 'fontsize',15);
xlabel('\itum\rm from distal end', 'fontsize',15);


%----------bargraph   (OFF duration)
figure('pos',[ 800 200 420 500])

hold off
bar(xcr+2.5, OFFDurATSpCR(1,:))
hold on
errorbar(xcr+2.5, OFFDurATSpCR(1,:), OFFDurATSpCR(2,:), 'k.', 'linewidth', 2);
axis([ -1 61 0 100 ])
box on
ylabel('ON duration of \itsygl-1\rm ATS (min)' , 'fontsize',15);
xlabel('\itum\rm from distal end', 'fontsize',15);










%%%%% =====================================================
%%%%% ==========   Further analyze ATS ==================== 



%%% (1) Summed ATS intensity over time
%%% The function 'MS2sumInt' records summed ATS intensity over time.
%%% A new variable 'MS2sumIntData' gets generated.
%%% The variable has 3 rows in a cell format | 1: nuc info | 2: summed ATS intensity over time|.
%%% In the first cell,
%%% | 1: Image ID | 2: Nucleus ID | 3: Location (um from distal end) | 4: Total # ATS/nuc.
%%% In the second cell,
%%% | Normalized summed ATS signal intensity in each point |
%%% In the third cell,

[MS2sumIntData_N2, perATSsync_N2, minATSsync_N2, ATS1ATS2int_N2, ATS1ATS2intRand_N2] = MS2sumInt(MS2norm_N2);
[MS2sumIntData_q224, perATSsync_q224, minATSsync_q224, ATS1ATS2int_q224, ATS1ATS2intRand_q224] = MS2sumInt(MS2norm_q224);

%%% (1) Correlation between ATS -- two loci within a nucleus over time.

%%% The variable 'allATSsInt' contains 3 columns:
%%% |1: ATS1 intensity in a nuc| 2: ATS2 intensity in the same nuc | 3: difference (2-1)|
allATSsInt_N2 = cell2mat(ATS1ATS2int_N2');
allATSsInt_q224 = cell2mat(ATS1ATS2int_q224');
allATSsInt_N2(:,3) = allATSsInt_N2(:,2) - allATSsInt_N2(:,1);
allATSsInt_q224(:,3) = allATSsInt_q224(:,2) - allATSsInt_q224(:,1);
allATSsIntRand_N2 = cell2mat(ATS1ATS2intRand_N2');
allATSsIntRand_q224 = cell2mat(ATS1ATS2intRand_q224');
allATSsIntRand_N2(:,3) = allATSsIntRand_N2(:,2) - allATSsIntRand_N2(:,1);
allATSsIntRand_q224(:,3) = allATSsIntRand_q224(:,2) - allATSsIntRand_q224(:,1);



% ------- boxplot for sig. intensity difference of two ATS within a nuc & random pairing.
figure,bxplot(allATSsInt_N2(:,3), allATSsIntRand_N2(:,3));
[~,p]=ttest2(allATSsInt_q224(:,3), allATSsIntRand_q224(:,3))
axis([0.2 2.8 -6000 6000])


% -------- violin plot for Difference in sig. intensity betw. two loci & random pairing.
kVal1 = allATSsInt_N2(:,3);
kVal2 = allATSsIntRand_N2(:,3);
kVal1(:,2) = 1;
kVal2(:,2) = 2;
kVal = [kVal1; kVal2];
figure('pos', [300 200 400 500])
distributionPlot(kVal(:,1), 'groups', kVal(:,2));
hold on
mV = [mean(kVal1(:,1)) mean(kVal2(:,1))];
sdV = [std(kVal1(:,1)) std(kVal2(:,1))];
errorbar(1:2,mV,[0 0], 'g-', 'linewidth', 2 )
errorbar(1:2, mV, sdV,  'g-', 'linewidth', 2)
axis([0.2 2.8 -6000 6000])

figure,bxplot(allATSsInt_q224(:,3), allATSsIntRand_q224(:,3));
[~,p]=ttest2(allATSsInt_q224(:,3), allATSsIntRand_q224(:,3))
axis([0.2 2.8 -6000 6000])

% -------- scatterplot for sig. intensity of two ATS within a nuc.
kVal = allATSsInt_N2;
figure('pos', [300 200 400 500])
hold off
plot(kVal(:,1), kVal(:,2), 'k.', 'markersize', 15);
corr(kVal(:,1), kVal(:,2))








% --------- Analysis for  ATS ON/OFF sync and async (binary)
%%% The variable 'perATSsync'
%%% 1st row: % both ATS are ON,
%%% 2nd:     % both ATS are OFF,
%%% 3rd:     % ATS are synchronized (ON or OFF),
%%% 4th:     % ATS are async,
%%% 5-8th: Repetition of 1-4th row but with probability calculated from original observation.

col = 1;
[a,b] = ttest2(perATSsync_N2(col,:), perATSsync_N2(col+4,:))


meanATSsync = mean(perATSsync_N2,2);
stdATSsync = [std(perATSsync_N2(1,:)); std(perATSsync_N2(2,:)); std(perATSsync_N2(3,:)); std(perATSsync_N2(4,:)); ...
    std(perATSsync_N2(5,:)); std(perATSsync_N2(6,:)); std(perATSsync_N2(7,:)); std(perATSsync_N2(8,:))];
seATSsync = stdATSsync/sqrt(size(perATSsync_N2,2)-1);

% --------- bargraph for binary % two ATS sync/async
figure('pos', [300 200 400 500])
yval = [meanATSsync(1) meanATSsync(5); meanATSsync(2) meanATSsync(6); meanATSsync(3) meanATSsync(7); meanATSsync(4) meanATSsync(8)];
yval2 = [meanATSsync(1) meanATSsync(5) meanATSsync(2) meanATSsync(6) meanATSsync(3) meanATSsync(7) meanATSsync(4) meanATSsync(8)];
erval = [seATSsync(1) seATSsync(5) seATSsync(2) seATSsync(6) seATSsync(3) seATSsync(7) seATSsync(4) seATSsync(8)];
bar(yval)
hold on
offset = .15;
errorbar([ 1-offset 1+offset 2-offset 2+offset 3-offset 3+offset 4-offset 4+offset],yval2,erval,  'k.', 'linewidth', 2)
axis([0.5 4.5 0 100]);



%%% -------- 2nd method (plot minutes instead of %) ------------
col = 4;
[a,b] = ttest2(minATSsync_N2(col,:), minATSsync_N2(col+4,:))


meanmATSsync = mean(minATSsync_N2,2);
stdmATSsync = [std(minATSsync_N2(1,:)); std(minATSsync_N2(2,:)); std(minATSsync_N2(3,:)); std(minATSsync_N2(4,:)); ...
    std(minATSsync_N2(5,:)); std(minATSsync_N2(6,:)); std(minATSsync_N2(7,:)); std(minATSsync_N2(8,:))];
semATSsync = stdATSsync/sqrt(size(minATSsync_N2,2)-1);

% --------- bargraph for binary % two ATS sync/async
figure('pos', [300 200 400 500])
yvalm = [meanmATSsync(1) meanmATSsync(5); meanmATSsync(2) meanmATSsync(6); meanmATSsync(4) meanmATSsync(8)];
yval2m = [meanmATSsync(1) meanmATSsync(5) meanmATSsync(2) meanmATSsync(6) meanmATSsync(4) meanmATSsync(8)];
ervalm = [semATSsync(1) semATSsync(5) semATSsync(2) semATSsync(6) semATSsync(4) semATSsync(8)];

bar(yvalm)
hold on
offset = .15;
errorbar([ 1-offset 1+offset 2-offset 2+offset 3-offset 3+offset],yval2m,ervalm,  'k.', 'linewidth', 2)
axis([0.5 3.5 0 50]);






%%% (2) Difference in ave. ATS intensity between two consecutive trx bursting events (the same locus).  
%%% A matrix variable 'MSnorm_ATSrel' is generated.
%%% In the variable, each data point is a mean ATS intensity during one trx bursting event.
%%% Each column is for each ATS (loci) with one or more bursting.
%%% An array 'MS2norm_ATSONdiff' records mean ATS intensity difference between two neighboring bursting event.
%%% 'MS2norm_ATSONdiffRand' randomizes ATS mean intensity (multi-events) within a nucleus.


MS2norm = MS2norm_N2;

MS2norm_ATSrel = zeros(99,99999);
MS2norm_ATSdur = zeros(99,99999);
MS2norm_ATSrelLoc = zeros(1,99999);
loc = 1;
for i = 1:size(MS2norm,2)
    for j = 1:size(MS2norm{3,i},2)
        tempATSave = [];   % record ATS intensity
        tempATSave2 = [];  % record burst duration
        for k = 1:floor(size(MS2norm{3,i},1)/2)
            pos = k*2;
            if MS2norm{3,i}(pos,j) > 0
                ONrange1 = floor(sum(MS2norm{3,i}(1:pos-1,j))/5);
                if ONrange1 < 1
                    ONrange1 = 1;
                end
                ONrange2 = floor(sum(MS2norm{3,i}(1:pos,j))/5);
                if ONrange2 > size(MS2norm{2,i},1)
                    ONrange2 = size(MS2norm{2,i},1);
                end
                temp = mean(MS2norm{2,i}(ONrange1:ONrange2,j));
                temp2 = MS2norm{3,i}(pos,j);
                tempATSave = [tempATSave; temp];
                tempATSave2 = [tempATSave2; temp2];
            end
        end
        if size(tempATSave,1) > 1
            MS2norm_ATSrel(1:size(tempATSave,1),loc) = tempATSave;
            MS2norm_ATSdur(1:size(tempATSave,1),loc) = tempATSave2;
            MS2norm_ATSrelLoc(1,loc) = MS2norm{1,i}(2,j);
            loc = loc+1;
        end
    end
end

MS2norm_ATSrel(:,loc:end) = [];
MS2norm_ATSdur(:,loc:end) = [];
MS2norm_ATSrelLoc(:,loc:end) = [];

MS2norm_ATSrel_N2 = MS2norm_ATSrel;
MS2norm_ATSdur_N2 = MS2norm_ATSdur;

MS2norm_ATSrel_q224 = MS2norm_ATSrel;



% ------ correlation betw. two burst INTENSITIES
%%% 'MS2norm_ATStwoEv' contains 3 columns
%%% |1: Intensity of 1st event | 2: Intensity of 2nd event | 3: difference (2-1)|
MS2norm_ATSrel = MS2norm_ATSrel_N2;

MS2norm_ATStwoEv = zeros(9999999,2);
MS2norm_ATStwoEvRand = zeros(9999999,2);

pos = 1;
for i = 1:size(MS2norm_ATSrel,2)
    idx = find(MS2norm_ATSrel(:,i) == 0,1)-1;
    rLoc = randperm(idx);
    randT = MS2norm_ATSrel(rLoc,i);
    temp = [MS2norm_ATSrel(1:idx-1,i) MS2norm_ATSrel(2:idx,i)];
    tempR = [randT(1:idx-1) randT(2:idx)];
    MS2norm_ATStwoEv(pos:pos+idx-2,1:2) = temp;
    MS2norm_ATStwoEv(pos:pos+idx-2,4) = MS2norm_ATSrelLoc(i);  %%% record um distance from DE in 4th column.
    MS2norm_ATStwoEvRand(pos:pos+idx-2,:) = tempR;
    
    pos = pos+idx-1;
end
MS2norm_ATStwoEv(pos:end,:) = [];
MS2norm_ATStwoEvRand(pos:end,:) = [];

%--- 2nd method for randomization
randLoc = randperm(size(MS2norm_ATStwoEv,1));
MS2norm_ATStwoEvRand2 = [ MS2norm_ATStwoEv(:,1) MS2norm_ATStwoEv(randLoc,2) ];

% ------ get difference
MS2norm_ATStwoEv(:,3) = MS2norm_ATStwoEv(:,2) - MS2norm_ATStwoEv(:,1);
MS2norm_ATStwoEvRand(:,3) = MS2norm_ATStwoEvRand(:,2) - MS2norm_ATStwoEvRand(:,1);
MS2norm_ATStwoEvRand2(:,3) = MS2norm_ATStwoEvRand2(:,2) - MS2norm_ATStwoEvRand2(:,1);


% -------- scatterplot for ATS intensity of 1st event vs. 2nd event.
kVal1 = MS2norm_ATStwoEv(MS2norm_ATStwoEv(:,4) >= 0 & MS2norm_ATStwoEv(:,4) < 10,:);
kVal2 = MS2norm_ATStwoEv(MS2norm_ATStwoEv(:,4) >= 10 & MS2norm_ATStwoEv(:,4) < 20,:);
kVal3 = MS2norm_ATStwoEv(MS2norm_ATStwoEv(:,4) >= 20 & MS2norm_ATStwoEv(:,4) < 30,:);
kVal4 = MS2norm_ATStwoEv(MS2norm_ATStwoEv(:,4) >= 30 & MS2norm_ATStwoEv(:,4) < 65,:);
figure('pos', [300 200 400 500])
hold on
plot(kVal1(:,1), kVal1(:,2), 'k.', 'markersize', 15);
plot(kVal2(:,1), kVal2(:,2), 'r.', 'markersize', 15);
plot(kVal3(:,1), kVal3(:,2), 'g.', 'markersize', 15);
plot(kVal4(:,1), kVal4(:,2), 'b.', 'markersize', 15);
axis([0 5000 0 4100])


valNow = kVal4;
corr(valNow(:,1), valNow(:,2))




% -------- barplot for difference
figure('pos', [300 200 400 500])
histogram(kVal(:,3));
box on

% -------- violin plot
kVal = MS2norm_ATStwoEvRand2;
figure('pos', [300 200 400 500])
distributionPlot(kVal(:,3));
hold on
mV = mean(kVal(:,3));
sdV = std(kVal(:,3));
errorbar(1,mV,0, 'g-', 'linewidth', 2 )
errorbar(1, mV, sdV,  'g-', 'linewidth', 2)

box on


MS2norm_ATStwoEv_N2 = MS2norm_ATStwoEv;
MS2norm_ATStwoEvRand2_N2 = MS2norm_ATStwoEvRand2;




%%% (3) Correlation between ON/OFF duration of consecutive trx. events.
%%% This code generates 4 variables. e.g. below.
%%% MS2norm_ATSdurON records: ON duration then next ON duration
%%% MS2norm_ATSdurOFF records: OFF duration then next OFF duration ...
%%%
%%% Each variable contains 4 columns:
%%% | 1: um from distal end to ATS | 2: Duration of 1st event (min)| 3: Duration of 2nd event (min) | 4: Difference in duration|||
MS2norm = MS2norm_N2;

MS2norm_ATSdurON = zeros(9999999,4);
MS2norm_ATSdurOFF = zeros(9999999,4);
MS2norm_ATSdurONoff = zeros(9999999,4);
MS2norm_ATSdurOFFon = zeros(9999999,4);
locon = 1;
locoff = 1;
loconOFF = 1;
locoffON = 1;
for i = 1:size(MS2norm,2)
    for j = 1:size(MS2norm{3,i},2)
        tempATSave = [];
        for k = 1:floor(size(MS2norm{3,i},1)/2)
            dat = MS2norm{3,i}(:,j);
            posOFF = k*2-1;
            posON = k*2;
            %%% Two neighboring OFF states
            if posOFF+2 <= size(MS2norm{3,i},1)
                if dat(posOFF) ~= 0 && dat(posOFF+2) ~=0
                    MS2norm_ATSdurOFF(locoff,1) = MS2norm{1,i}(2,j);
                    MS2norm_ATSdurOFF(locoff,2:3) = [dat(posOFF) dat(posOFF+2)];
                    MS2norm_ATSdurOFF(locoff,4) = dat(posOFF) - dat(posOFF+2);
                    locoff = locoff+1;
                end
            end
            %%% Two neighboring ON states
            if posON+2 <= size(MS2norm{3,i},1)
                if dat(posON) ~= 0 && dat(posON+2) ~=0
                    MS2norm_ATSdurON(locon,1) = MS2norm{1,i}(2,j);
                    MS2norm_ATSdurON(locon,2:3) = [dat(posON) dat(posON+2)];
                    MS2norm_ATSdurON(locon,4) = dat(posON) - dat(posON+2);
                    locon = locon+1;
                end
            end
            %%% Consecutive OFF then ON states
            if posOFF+1 <= size(MS2norm{3,i},1)
                if dat(posOFF) ~= 0 && dat(posOFF+1) ~=0
                    MS2norm_ATSdurOFFon(locoffON,1) = MS2norm{1,i}(2,j);
                    MS2norm_ATSdurOFFon(locoffON,2:3) = [dat(posOFF) dat(posOFF+1)];
                    MS2norm_ATSdurOFFon(locoffON,4) = dat(posOFF) - dat(posOFF+1);
                    locoffON = locoffON+1;
                end
            end
            %%% Consecutive ON then OFF states
            if posON+1 <= size(MS2norm{3,i},1)
                if dat(posON) ~= 0 && dat(posON+1) ~=0
                    MS2norm_ATSdurONoff(loconOFF,1) = MS2norm{1,i}(2,j);
                    MS2norm_ATSdurONoff(loconOFF,2:3) = [dat(posON) dat(posON+1)];
                    MS2norm_ATSdurONoff(loconOFF,4) = dat(posON) - dat(posON+1);
                    loconOFF = loconOFF+1;
                end
            end
        end
    end
end

MS2norm_ATSdurON(locon:end,:) = [];
MS2norm_ATSdurOFF(locoff:end,:) = [];
MS2norm_ATSdurONoff(loconOFF:end,:) = [];
MS2norm_ATSdurOFFon(locoffON:end,:) = [];
      

% ------ scatter plot and Pearson's r
kVal = MS2norm_ATSdurONoff;
figure('pos', [300 200 400 500])
hold off
plot(kVal(:,2), kVal(:,3), 'k.', 'markersize', 15);  

corr(kVal(:,2), kVal(:,3))



% -------- scatterplot for ON duration of 1st event vs. 2nd event.
kVal = MS2norm_ATSdurOFFon;
kVal1 = kVal(kVal(:,1) >= 0 & kVal(:,1) < 10,:);
kVal2 = kVal(kVal(:,1) >= 10 & kVal(:,1) < 20,:);
kVal3 = kVal(kVal(:,1) >= 20 & kVal(:,1) < 30,:);
kVal4 = kVal(kVal(:,1) >= 30 & kVal(:,1) < 65,:);
figure('pos', [300 200 400 500])
hold on
plot(kVal1(:,2), kVal1(:,3), 'k.', 'markersize', 15);
plot(kVal2(:,2), kVal2(:,3), 'r.', 'markersize', 15);
plot(kVal3(:,2), kVal3(:,3), 'g.', 'markersize', 15);
plot(kVal4(:,2), kVal4(:,3), 'b.', 'markersize', 15);
box on
corr(kVal(:,2), kVal(:,3))



valNow = kVal4;

corr(valNow(:,2), valNow(:,3))


% ------- boxplot for duration difference
kVal = MS2norm_ATSdurON;
figure('pos', [300 200 400 500])
hist(kVal(:,4), 40);
box on

% ------- violin plot
kVal = MS2norm_ATSdurOFF;
figure('pos', [300 200 160 500])
distributionPlot(kVal(:,4));
hold on
mV = mean(kVal(:,4));
sdV = std(kVal(:,4));
errorbar(1,mV,0, 'g-', 'linewidth', 2 )
errorbar(1, mV, sdV,  'g-', 'linewidth', 2)
box on




%%% ======== (4) # ATS per nucleus in space and time
%%% The variable 'meanATStime' records:
%%% |1: um distance from Distal end | 2: average # ATS per nucleus | 3: ave. summed ATS intensity per nuc|

MS2sumIntData = MS2sumIntData_N2;

meanATStime = zeros(99999999,2);
loc = 1;

for i=1:size(MS2sumIntData,2)
    nucDist = MS2sumIntData{1,i}(3,:);
    meanATSnum = mean(MS2sumIntData{3,i},1);
    meanATSint = zeros(1,length(nucDist));
    for j=1:length(nucDist)
        meanATSint(j) = mean(MS2sumIntData{2,i}(MS2sumIntData{3,i}(:,j)>0,j));
    end
    siz = length(nucDist);
    meanATStime(loc:loc+siz-1,1:3) = [nucDist' meanATSnum' meanATSint'];
    loc = loc+siz;
end

meanATStime(loc:end,:) = [];

% calculate ave. # ATS and summed ATS intensity
crs=5;

xcr = 0:crs:60;
bins = length(xcr);

resultNumATSpCR = zeros(3,bins);
resultIntATSpCR = zeros(3,bins);

for j=1:bins     % 1 - 13th cell row
    xr = [crs*(j-1) crs*j];
    sSize = sum(meanATStime(meanATStime(:,1) >=xr(1) & meanATStime(:,1) < xr(2),2));
    if sSize ~= 0
        tempNmean = mean(meanATStime(meanATStime(:,1) >=xr(1) & meanATStime(:,1) < xr(2),2));
        tempNstd = std(meanATStime(meanATStime(:,1) >=xr(1) & meanATStime(:,1) < xr(2),2));
        tempNse = tempNstd/sqrt(sSize);

        tempImean = mean(meanATStime(meanATStime(:,1) >=xr(1) & meanATStime(:,1) < xr(2),3));
        tempIstd = std(meanATStime(meanATStime(:,1) >=xr(1) & meanATStime(:,1) < xr(2),3));
        tempIse = tempNstd/sqrt(sSize);

        resultNumATSpCR(:,j) = [tempNmean; tempNstd; tempNse];
        resultIntATSpCR(:,j) = [tempImean; tempIstd; tempIse];
    else
        resultNumATSpCR(:,j) = 0;
        resultIntATSpCR(:,j) = 0;
    end
end

% gets mean of transcriptional activity
gap = 1;
mLine = zeros(1,30);
seL = zeros(1,30);
for i=1:30
    cpool = meanATStime(meanATStime(:,1) >= (i-1)*2 & meanATStime(:,1) < i*2,3);
    if isempty(cpool)
        cpool = 0;
    end
    mLine(i) = mean(cpool);
    seL(i) = std(cpool)/sqrt(length(cpool));
end

%----------bargraph for #ATS/cell averaged over time
figure('pos', [300 200 350 500])
y_max = 2;
hold off
bar(xcr+crs/2, resultNumATSpCR(1,:))
hold on
errorbar(xcr+crs/2, resultNumATSpCR(1,:), resultNumATSpCR(3,:), 'k.', 'linewidth', 2);
axis([ -4 64 0 y_max ])
xticks(0:10:100)
box on
% ylabel('# \itsygl-1\rm ATS per cell' , 'fontsize',15);
% xlabel('\itum\rm from distal end', 'fontsize',15);

%-----------scatter plot for summed ATS intensity per cell averaged over time.
figure('pos', [300 200 400 500])
hold off
plot(meanATStime(:,1), meanATStime(:,3), 'k.', 'markersize', 15);     % L4: x2
y_max2 = 8000;
axis([-1 61  0 y_max2 ])
% xlabel('Radius of Nucleus (\mu\itm\rm)' , 'fontsize',15);
% ylabel('DAPI signal intensity (\ita.u.\rm)', 'fontsize',15); 
% text(1.2, 10*1e5, strcat('r = ',num2str(round(corrAll,2,'significant'))), 'color', 'r', 'fontsize', 20);
hold on
plot(1:2:60, mLine, 'c', 'linewidth', 2);
ylabel('Summed ATS intensities (a.u.)' , 'fontsize',15);
xlabel('\itum\rm from distal end', 'fontsize',15);
















