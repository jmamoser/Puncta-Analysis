%close all;
%clear all;
%%% filepaths %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
projectpath='\\spencerstorage.int.colorado.edu\LabShare\IXMicroImages-goodNames\Justin\'; %path where the processed data will be stored; create a folder called "Data" here
imagepath=  '\\spencerstorage.int.colorado.edu\LabShare\IXMicroImages-goodNames\Justin\';  %where the images are stored, generally
experimentpath='MANSI_ANALYSIS_20170104\';
% experimentpath='MA56-20160211-VitC-NAC-pilot3colorMCF10A_782\';
datadir=([projectpath,experimentpath,'Data_MA72\']);

%%% analysis settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
motheroption=0; %0:no gating 1:mothers 2:no mothers
daughteroption=1; %0:no gating 1:daughters 2:no daughters
framesaftermitosis = 20; % categorizes based on puncta existence between POM and these many frames after mitosis.
framenum = 120; % number of frames of the entire movie to look at.
quiescentanalysis=0;
if quiescentanalysis
    motheroption=2; daughteroption=2;
end
IFoption=0; %0: No IF data 1: IF data
framesperhr=5;
smoothingsize=5;
interval = 6:20; upperint = 60; lowerint = 100-upperint;
fprintf(['\nInterval: [' num2str(interval(1)) ':' num2str(interval(end)) ']'])
fprintf(['\nUpper percentile: ' num2str(upperint)]);
fprintf(['\nLower percentile: ' num2str(lowerint) '\n']);
%%% pool data from multiple wells %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rowmat=[2 3 5];
fprintf(['Well: ' num2str(rowmat) '\n'])
colmat=[2];
sitemat=[1];
tracesCDK2inc = [];
CDK2incPOM = [];
CDK2incPuncta = [];

tracesCDK2low = [];
CDK2lowPOM = [];
CDK2lowPuncta = [];

tracesCDK2emerg = [];
CDK2emergPOM = [];
CDK2emergPuncta = [];

for row=rowmat
    for col=colmat
        for site=sitemat
            shot=[num2str(row),'_',num2str(col),'_',num2str(site)];
            load([datadir,'Row_',num2str(row),'_Col_',num2str(col),'_allgoodcdk2cells.mat']);%load tracedata
            
            tracesCDK2inc=[tracesCDK2inc;cdk2inctraces];
            CDK2incPOM=[CDK2incPOM, cdk2incframeofmitosis];
            CDK2incPuncta=[CDK2incPuncta;cdk2incpuncta];
            
            tracesCDK2low=[tracesCDK2low;cdk2lowtraces];
            CDK2lowPOM=[CDK2lowPOM, cdk2lowframeofmitosis];
            CDK2lowPuncta=[CDK2lowPuncta;cdk2lowpuncta];
            
            tracesCDK2emerg=[tracesCDK2emerg;cdk2emergtraces];
            CDK2emergPOM=[CDK2emergPOM, cdk2emergframeofmitosis];
            CDK2emergPuncta=[CDK2emergPuncta;cdk2emergpuncta];
            
        end
    end
end

%% % align traces %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
frames=1:309; INC = 1; LOW = 0; EMERG = 2;

cdk2incframesmat=repmat(frames, size(tracesCDK2inc,1), 1);
cdk2incmitosismat=repmat(CDK2incPOM, length(frames), 1)';
cdk2incshiftedframes=cdk2incframesmat-cdk2incmitosismat;
cdk2incshiftedframes=cdk2incshiftedframes./framesperhr ;
tracesCDK2inc_24 = tracesCDK2inc(:,1:309);
punctaCDK2inc_24 = CDK2incPuncta(:,1:309); %JM
INCclass(1:size(punctaCDK2inc_24),1) = INC;

cdk2lowframesmat=repmat(frames, size(tracesCDK2low,1), 1);
cdk2lowmitosismat=repmat(CDK2lowPOM, length(frames), 1)';
cdk2lowshiftedframes=cdk2lowframesmat-cdk2lowmitosismat;
cdk2lowshiftedframes=cdk2lowshiftedframes./framesperhr ;
tracesCDK2low_24 = tracesCDK2low(:,1:309);
punctaCDK2low_24 = CDK2lowPuncta(:,1:309); %JM
LOWclass(1:size(punctaCDK2low_24),1) = LOW;

cdk2emergframesmat=repmat(frames, size(tracesCDK2emerg,1), 1);
cdk2emergmitosismat=repmat(CDK2emergPOM, length(frames), 1)';
cdk2emergshiftedframes=cdk2emergframesmat-cdk2emergmitosismat;
cdk2emergshiftedframes=cdk2emergshiftedframes./framesperhr ;
tracesCDK2emerg_24 = tracesCDK2emerg(:,1:309);
punctaCDK2emerg_24 = CDK2emergPuncta(:,1:309); %JM
EMERGclass(1:size(punctaCDK2emerg_24),1) = EMERG;

timestore       = [cdk2incshiftedframes; cdk2lowshiftedframes; cdk2emergshiftedframes];
CDK2tracestore  = [tracesCDK2inc_24; tracesCDK2low_24; tracesCDK2emerg_24];
POMstore        = [cdk2incmitosismat; cdk2lowmitosismat; cdk2emergmitosismat];
punctastore     = [punctaCDK2inc_24; punctaCDK2low_24; punctaCDK2emerg_24];
classstore      = [INCclass; LOWclass; EMERGclass];

numtraces = length(CDK2incPOM) + length(CDK2lowPOM) + length(CDK2emergPOM);

%% % Plotting based on Classification 2:


%%% for aligning time
numframes = 309;
lengthleftinc=ones(length(CDK2incPOM),1)*NaN; lengthrightinc=lengthleftinc;
lengthleftlow=ones(length(CDK2lowPOM),1)*NaN; lengthrightlow=lengthleftlow;
lengthleftemerg=ones(length(CDK2emergPOM),1)*NaN; lengthrightemerg=lengthleftemerg;
incpunctatraceright=ones(length(CDK2incPOM),numframes)*NaN;
incpunctatraceleft=ones(length(CDK2incPOM),numframes)*NaN;
lowpunctatraceright=ones(length(CDK2lowPOM),numframes)*NaN;
lowpunctatraceleft=ones(length(CDK2lowPOM),numframes)*NaN;
emergpunctatraceright=ones(length(CDK2emergPOM),numframes)*NaN;
emergpunctatraceleft=ones(length(CDK2emergPOM),numframes)*NaN;
tracedataincleft=ones(length(CDK2incPOM),numframes)*NaN;
tracedataincright=ones(length(CDK2incPOM),numframes)*NaN;
tracedatalowleft=ones(length(CDK2lowPOM),numframes)*NaN;
tracedatalowright=ones(length(CDK2lowPOM),numframes)*NaN;
tracedataemergleft=ones(length(CDK2emergPOM),numframes)*NaN;
tracedataemergright=ones(length(CDK2emergPOM),numframes)*NaN;

for i=1:length(CDK2incPOM)
    
    plotdata=tracesCDK2inc_24(i,:);
    realframes=find(~isnan(tracesCDK2inc_24(i,:)));
    plotdata(realframes)=smooth(tracesCDK2inc_24(i,realframes),smoothingsize);
    firstframe=find(~isnan(tracesCDK2inc_24(i,:)),1,'first');
    lastframe=find(~isnan(tracesCDK2inc_24(i,:)),1,'last');
    lengthleftinc(i)=CDK2incPOM(i)-firstframe+1;
    lengthrightinc(i)=lastframe-CDK2incPOM(i)+1;
    incpunctatraceleft(i,numframes-lengthleftinc(i)+1:numframes)=CDK2incPuncta(i,firstframe:CDK2incPOM(i));
    incpunctatraceright(i,1:lengthrightinc(i))=CDK2incPuncta(i,CDK2incPOM(i):lastframe);
    tracedataincleft(i,numframes-lengthleftinc(i)+1:numframes)=tracesCDK2inc_24(i,firstframe:CDK2incPOM(i));
    tracedataincright(i,1:lengthrightinc(i))=tracesCDK2inc_24(i,CDK2incPOM(i):lastframe);
    
end

lengthleftincmax=floor(prctile(lengthleftinc,90));
rightlengthincmax=floor(prctile(lengthrightinc,90));
alignedtime1=(-lengthleftincmax+1:rightlengthincmax);
incpunctatraceleft=incpunctatraceleft(:,numframes-lengthleftincmax+1:numframes);
incpunctatraceright=incpunctatraceright(:,1:rightlengthincmax);
incpuncta_datatotal = [incpunctatraceleft incpunctatraceright];
tracedataincleft=tracedataincleft(:,numframes-lengthleftincmax+1:numframes);
tracedataincright=tracedataincright(:,1:rightlengthincmax);
incdatatotal=[tracedataincleft tracedataincright];

%20170113 found bug; data from column 99 gets duplicated in column 100.
%Take out this column for incdatatotal, incpuncta_datatotal, also do so for
%the low and emerg.
incdatatotal_edit = incdatatotal(:,[1:lengthleftincmax lengthleftincmax+2:end]);
incpuncta_datatotal_edit = incpuncta_datatotal(:,[1:lengthleftincmax lengthleftincmax+2:end]);
alignedtime1_edit = alignedtime1(1:end-1);

for i=1:length(CDK2lowPOM)
    
    plotdata=tracesCDK2low_24(i,:);
    realframes=find(~isnan(tracesCDK2low_24(i,:)));
    plotdata(realframes)=smooth(tracesCDK2low_24(i,realframes),smoothingsize);
    firstframe=find(~isnan(tracesCDK2low_24(i,:)),1,'first');
    lastframe=find(~isnan(tracesCDK2low_24(i,:)),1,'last');
    lengthleftlow(i)=CDK2lowPOM(i)-firstframe+1;
    lengthrightlow(i)=lastframe-CDK2lowPOM(i)+1;
    lowpunctatraceleft(i,numframes-lengthleftlow(i)+1:numframes)=CDK2lowPuncta(i,firstframe:CDK2lowPOM(i));
    lowpunctatraceright(i,1:lengthrightlow(i))=CDK2lowPuncta(i,CDK2lowPOM(i):lastframe);
    tracedatalowleft(i,numframes-lengthleftlow(i)+1:numframes)=tracesCDK2low_24(i,firstframe:CDK2lowPOM(i));
    tracedatalowright(i,1:lengthrightlow(i))=tracesCDK2low_24(i,CDK2lowPOM(i):lastframe);
end

lengthleftlowmax=floor(prctile(lengthleftlow,90));
rightlengthlowmax=floor(prctile(lengthrightlow,90));
alignedtime2=(-lengthleftlowmax+1:rightlengthlowmax);
lowpunctatraceleft=lowpunctatraceleft(:,numframes-lengthleftlowmax+1:numframes);
lowpunctatraceright=lowpunctatraceright(:,1:rightlengthlowmax);
lowpuncta_datatotal = [lowpunctatraceleft lowpunctatraceright];
tracedatalowleft=tracedatalowleft(:,numframes-lengthleftlowmax+1:numframes);
tracedatalowright=tracedatalowright(:,1:rightlengthlowmax);
lowdatatotal=[tracedatalowleft tracedatalowright];

lowdatatotal_edit = lowdatatotal(:,[1:lengthleftlowmax lengthleftlowmax+2:end]);
lowpuncta_datatotal_edit = lowpuncta_datatotal(:,[1:lengthleftlowmax lengthleftlowmax+2:end]);
alignedtime2_edit = alignedtime2(1:end-1);


for i=1:length(CDK2emergPOM)
    
    plotdata=tracesCDK2emerg_24(i,:);
    realframes=find(~isnan(tracesCDK2emerg_24(i,:)));
    plotdata(realframes)=smooth(tracesCDK2emerg_24(i,realframes),smoothingsize);
    firstframe=find(~isnan(tracesCDK2emerg_24(i,:)),1,'first');
    lastframe=find(~isnan(tracesCDK2emerg_24(i,:)),1,'last');
    lengthleftemerg(i)=CDK2emergPOM(i)-firstframe+1;
    lengthrightemerg(i)=lastframe-CDK2emergPOM(i)+1;
    emergpunctatraceleft(i,numframes-lengthleftemerg(i)+1:numframes)=CDK2emergPuncta(i,firstframe:CDK2emergPOM(i));
    emergpunctatraceright(i,1:lengthrightemerg(i))=CDK2emergPuncta(i,CDK2emergPOM(i):lastframe);
    tracedataemergleft(i,numframes-lengthleftemerg(i)+1:numframes)=tracesCDK2emerg_24(i,firstframe:CDK2emergPOM(i));
    tracedataemergright(i,1:lengthrightemerg(i))=tracesCDK2emerg_24(i,CDK2emergPOM(i):lastframe);
end


lengthleftemergmax=floor(prctile(lengthleftemerg,90));
rightlengthemergmax=floor(prctile(lengthrightemerg,90));
alignedtime3=(-lengthleftemergmax+1:rightlengthemergmax);
emergpunctatraceleft=emergpunctatraceleft(:,numframes-lengthleftemergmax+1:numframes);
emergpunctatraceright=emergpunctatraceright(:,1:rightlengthemergmax);
emergpuncta_datatotal = [emergpunctatraceleft emergpunctatraceright];
tracedataemergleft=tracedataemergleft(:,numframes-lengthleftemergmax+1:numframes);
tracedataemergright=tracedataemergright(:,1:rightlengthemergmax);
emergdatatotal=[tracedataemergleft tracedataemergright];

emergdatatotal_edit = emergdatatotal(:,[1:lengthleftemergmax lengthleftemergmax+2:end]);
emergpuncta_datatotal_edit = emergpuncta_datatotal(:,[1:lengthleftemergmax lengthleftemergmax+2:end]);
alignedtime3_edit = alignedtime3(1:end-1);


alignedtime1 = alignedtime1./framesperhr; %inc time
alignedtime2 = alignedtime2./framesperhr; %low time
alignedtime3 = alignedtime3./framesperhr; %emerging time

alignedtime1_edit = alignedtime1_edit./framesperhr; %inc time
alignedtime2_edit = alignedtime2_edit./framesperhr; %low time
alignedtime3_edit = alignedtime3_edit./framesperhr; %emerging time



plotdata2_mean = nanmean(incpuncta_datatotal); %inc puncta mean
plotdata3_mean = nanmean(lowpuncta_datatotal); %low puncta mean
plotdata4_mean = nanmean(emergpuncta_datatotal); %emerging puncta mena
plotCdk2emergdata_median = nanmedian(emergdatatotal);
plotCdk2emergdata_mean = nanmean(emergdatatotal);
plotCdk2incdata_median = nanmedian(incdatatotal);
plotCdk2incdata_mean = nanmean(incdatatotal);
plotCdk2lowdata_median = nanmedian(lowdatatotal);

%% 20170113: Pull out Cell ID for each 
% load([datadir,'Row_' num2str(rowmat) '_Col_' num2str(colmat) '_' 'keyinfo_Cdk2traces.mat']); %loads genealogy, tracedata and jitters. We only care about genealogy
% %load([datadir,'20170113Aligned']);
% count = 1;
% for i = 1:5%length(incdatatotal)
%     pattern = tracesCDK2inc(i,~isnan(tracesCDK2inc(i,:)));
%     disp(pattern)
%     pattern1 = incdatatotal(i,~isnan(incdatatotal(i,:)));
%     disp(pattern1)
% %     if length(pattern) > 20
% %         pattern = pattern(1:20);
% %     end
%     for j = 1:length(tracesCDK2)
%         str = tracesCDK2(j,:);
%         x = strfind(str,pattern);
%         %count = count+1;
%         if ~isempty(x)
%             disp(x)
%             store(count,1) = x;
%             store(count,2) = i;
%             store(count,3) = j;
%             count = count+1;
%             continue
%         end      
%     end 
% end
% 
% for i = 1:length(lowdatatotal)
%     
% end
% for i = 1:length(emergdatatotal)
%     
% end


%Create matrix for sum, mean, median puncta in given interval after mitosis
%%
punctaMatrix = NaN(length(incdatatotal_edit)+length(lowdatatotal_edit)+length(emergdatatotal_edit),4);
inczero = find(alignedtime1_edit == 0);
for i = 1:size(incdatatotal,1)
    if sum(isnan(incpuncta_datatotal_edit(i,inczero+interval))) < 1 && nansum(incpuncta_datatotal_edit(i,:)) >= 1 %quick and dirty way of making sure that all the cells included in this analysis have the 53BP1 sensor
        %the above line makes sure there aren't any NaNs in the relevant
        %interval after the align time, and that the cell definitely
        %expresses the 53BP1mChy sensor
        %if ~isempty(incpuncta_datatotal(i,inczero+interval))
            
            puncta = incpuncta_datatotal_edit(i,inczero+interval);
            %fprintf('\nPuncta:\n')
            %disp(puncta)
            nonzeros = find(puncta > 0); %returns indices of actual puncta
            %disp(nonzeros)
            if ~isempty(nonzeros) && length(nonzeros) >= 2
                for k = 1:length(nonzeros)-1
                    if  nonzeros(k) > 1 && nonzeros(k) < length(puncta)%k >1 && k < length(nonzeros)
                        if puncta(nonzeros(k)+1) == 0 && puncta(nonzeros(k)-1) == 0
                            %if k > 1
                            %if puncta(nonzeros(k)-1) == 0
                            puncta(nonzeros(k))=0;
                            %end
                            %end
                        end
                    end
                end
            %else
                %puncta(nonzeros) = 0;
            end
            %fprintf('\nAmended:\n')
            %disp(puncta)
            punctaMatrix(i,1) = nansum(puncta);
            punctaMatrix(i,2) = nanmean(puncta);
            punctaMatrix(i,3) = nanmedian(puncta);
            punctaMatrix(i,4) = nanstd(puncta);
    else
        punctaMatrix(i,1:4) = NaN;
    end
end

x=i;
lowzero = find(alignedtime2_edit == 0);
for i = 1:size(lowdatatotal,1)
    if sum(isnan(lowpuncta_datatotal_edit(i,lowzero+interval))) < 1 && nansum(lowpuncta_datatotal_edit(i,:)) >= 1 %quick and dirty way of making sure that all the cells included in this analysis have the 53BP1 sensor
        %the above line makes sure there aren't any NaNs in the relevant
        %interval after the align time, and that the cell definitely
        %expresses the 53BP1mChy sensor
        %if ~isempty(incpuncta_datatotal(i,inczero+interval))
            
            puncta = lowpuncta_datatotal_edit(i,lowzero+interval);
            %disp(puncta)
            nonzeros = find(puncta > 0); %returns indices of actual puncta
            %disp(nonzeros)
            if ~isempty(nonzeros)
                for k = 1:length(nonzeros)-1
                    if  nonzeros(k) > 1 && nonzeros(k) < length(puncta)
                        if puncta(nonzeros(k)+1) == 0 && puncta(nonzeros(k)-1) == 0
                            %if k > 1
                            %if puncta(nonzeros(k)-1) == 0
                            puncta(nonzeros(k))=0;
                            %end
                            %end
                        end
                    end
                end
            end
            punctaMatrix(x+i,1) = nansum(puncta);
            punctaMatrix(x+i,2) = nanmean(puncta);
            punctaMatrix(x+i,3) = nanmedian(puncta);
            punctaMatrix(x+i,4) = nanstd(puncta);
    else
        punctaMatrix(x+i,1:4) = NaN;
    end
end

x=x+i;
emergzero = find(alignedtime3_edit == 0);
for i = 1:size(emergdatatotal,1)
    if sum(isnan(emergpuncta_datatotal_edit(i,emergzero+interval))) < 1 && nansum(emergpuncta_datatotal_edit(i,:)) >= 1 %quick and dirty way of making sure that all the cells included in this analysis have the 53BP1 sensor
        %the above line makes sure there aren't any NaNs in the relevant
        %interval after the align time, and that the cell definitely
        %expresses the 53BP1mChy sensor
        %if ~isempty(incpuncta_datatotal(i,inczero+interval))
            
            puncta = emergpuncta_datatotal_edit(i,emergzero+interval);
            %fprintf('\nPuncta:\n')
            %disp(puncta)
            nonzeros = find(puncta > 0); %returns indices of actual puncta
            %disp(nonzeros)
            if ~isempty(nonzeros)
                for k = 1:length(nonzeros)-1
                    if  nonzeros(k) > 1 && nonzeros(k) < length(puncta)
                        if puncta(nonzeros(k)+1) == 0 && puncta(nonzeros(k)-1) == 0
                            %if k > 1
                            %if puncta(nonzeros(k)-1) == 0
                            puncta(nonzeros(k))=0;
                            %end
                            %end
                        end
                    end
                end
            end
            %fprintf('\nAmended:\n')
            %disp(puncta)
            punctaMatrix(x+i,1) = nansum(puncta);
            punctaMatrix(x+i,2) = nanmean(puncta);
            punctaMatrix(x+i,3) = nanmedian(puncta);
            punctaMatrix(x+i,4) = nanstd(puncta);
    else
        punctaMatrix(x+i,1:4) = NaN;
    end
end

%now rank each column in descending order, much as we did before
[critSum,ixsum] = sort(punctaMatrix(:,1),'descend');
[critMean,ixmean] = sort(punctaMatrix(:,2),'descend');
[critMed,ixmed] = sort(punctaMatrix(:,3),'descend');
[critStd,ixstd] = sort(punctaMatrix(:,4),'descend');

validSum = find(~isnan(critSum));
validMean = find(~isnan(critMean));
validMed = find(~isnan(critMed));
validStd = find(~isnan(critStd));

critSumval = critSum(validSum); ixsum=ixsum(validSum);
critMeanval = critMean(validMean); ixmean=ixmean(validMean);
critMedval = critMed(validMed); ixmed=ixmed(validMed);
critStdval = critStd(validStd); ixstd=ixstd(validStd);

%now to plot the top 95th and bottom 5th intervals
topSum = prctile(critSumval,upperint);
bottomSum = prctile(critSumval,lowerint);
%topSum = 1;
%bottomSum = 0;
tops = find(critSumval >= topSum);
bottoms = find(critSumval <= bottomSum);
fprintf(['\n' num2str(length(tops)) '\n'])
fprintf([num2str(length(bottoms)) '\n'])
inc_hipuncta = 0; inc_lowpuncta = 0;
low_hipuncta = 0; low_lowpuncta = 0;
emerg_hipuncta = 0; emerg_lowpuncta = 0;
figure
suptitle(['Grouped according to sum of puncta in interval [' num2str(interval(1)) ':' num2str(interval(end)) ']'])
subplot(2,2,1)
title(['Top ' num2str(upperint) '% of cells with puncta, >= ' num2str(topSum)])
xlabel('Time after anaphase (hours)')
ylabel('CDK2 activity')
subplot(2,2,2)
title(['Top ' num2str(upperint) '% of cells with puncta, >= ' num2str(topSum)])
xlabel('Time after anaphase (hours)')
ylabel('Number of puncta')
subplot(2,2,3)
title(['Bottom ' num2str(lowerint) '% of cells with puncta, <= ' num2str(bottomSum)])
xlabel('Time after anaphase (hours)')
ylabel('CDK2 activity')
subplot(2,2,4)
title(['Bottom ' num2str(lowerint) '% of cells with puncta, <=' num2str(bottomSum)])
xlabel('Time after anaphase (hours)')
ylabel('Number of puncta')
for i = 1:(length(tops) + length(bottoms))
    %first things first, get the index for that cell sin the in the
    %original trace and puncta matrices
    %tops holds the index in the critSum matrix
    
    if i <= length(tops)
        index = ixsum(tops(i));
        plot = 0;
        %color = 'r';
    elseif i <= (length(bottoms) + length(tops))
        index = ixsum(bottoms(i-length(tops)));
        plot = 2;
        %color = 'b';
    end
    
    %index will ultimately be used to point towards one of three matrices:
    %incdata, lowdata, emergdata and their corresponding puncta matrices.
    %As such, if the index is less than or equal to the size of inc, then
    %use it to pull from inc. If bigger than inc but less than or equal to
    %the sum of inc and low, pull from low, etc.
    %However, this makes a problem, as index is for all data together, so
    %to deal with this, if index > = the number of inc traces, then we'll
    %need to subtract that number, and so on.
    if i <= (length(tops) + length(bottoms))
        if index <= size(incdatatotal,1)
            time = alignedtime1_edit;
            CDK2 = incdatatotal_edit(index,:);
            puncta = incpuncta_datatotal_edit(index,:);
            color = 'b';
            if i <= length(tops)
                inc_hipuncta = inc_hipuncta+1;
            end
            if i > length(tops)
                inc_lowpuncta = inc_lowpuncta+1;
            end
        elseif index > size(incdatatotal,1) && index <= (size(incdatatotal,1) + size(lowdatatotal,1))
            index = index - size(incdatatotal_edit,1);
            time = alignedtime2_edit;
            CDK2 = lowdatatotal_edit(index,:);
            puncta = lowpuncta_datatotal_edit(index,:);
            color = 'r';
            if i <= length(tops)
                low_hipuncta = low_hipuncta+1;
            end
            if i > length(tops)
                low_lowpuncta = low_lowpuncta+1;
            end
        elseif index > (size(incdatatotal_edit,1) + size(lowdatatotal_edit,1))
            index = index - (size(incdatatotal_edit,1) + size(lowdatatotal_edit,1));
            time = alignedtime3_edit;
            CDK2 = emergdatatotal_edit(index,:);
            puncta = emergpuncta_datatotal_edit(index,:);
            color = 'g';
            if i <= length(tops)
                emerg_hipuncta = emerg_hipuncta+1;
            end
            if i > length(tops)
                emerg_lowpuncta = emerg_lowpuncta+1;
            end
            
        end
        
        
        subplot(2,2,1+plot)
        line(time,CDK2,'color',color)
        xlim([-25 20])
        ylim([0.2 2])
        
        hold on;
        subplot(2,2,2+plot)
        line(time,puncta,'color',color)
        xlim([-25 20])
        ylim([0 10])
        hold on;
    end
end
fprintf('\nSUM')
fprintf(['\n Number of CDK2inc in top 95th percentile: ' num2str(inc_hipuncta)])
fprintf(['\n Number of CDK2low in top 95th percentile: ' num2str(low_hipuncta)])
fprintf(['\n Number of CDK2emerg in top 95th percentile: ' num2str(emerg_hipuncta)])

fprintf(['\n Number of CDK2inc in bottom 5th percentile: ' num2str(inc_lowpuncta)])
fprintf(['\n Number of CDK2low in bottom 5th percentile: ' num2str(low_lowpuncta)])
fprintf(['\n Number of CDK2emerg in bottom 5th percentile: ' num2str(emerg_lowpuncta)])

fprintf(['\n Percentage of top percentile that''s CDK2inc: ' num2str(inc_hipuncta/(inc_hipuncta + low_hipuncta + emerg_hipuncta)) '%'])
fprintf(['\n Percentage of bottom percentile that''s CDK2inc: ' num2str(inc_lowpuncta/(inc_lowpuncta + low_lowpuncta + emerg_lowpuncta)) '%'])

topMean = prctile(critMeanval,upperint);
bottomMean = prctile(critMeanval,lowerint);
%topMean = 1;
%bottomMean = 0;

tops = find(critMeanval >= topMean);
bottoms = find(critMeanval <= bottomMean);
inc_hipuncta = 0; inc_lowpuncta = 0;
low_hipuncta = 0; low_lowpuncta = 0;
emerg_hipuncta = 0; emerg_lowpuncta = 0;
figure
suptitle(['Grouped according to mean puncta in interval [' num2str(interval(1)) ':' num2str(interval(end)) ']'])
subplot(2,2,1)
title(['Top ' num2str(upperint) '% of cells with puncta, >= ' num2str(topMean)])
xlabel('Time after anaphase (hours)')
ylabel('CDK2 activity')
subplot(2,2,2)
title(['Top ' num2str(upperint) '% of cells with puncta, >= ' num2str(topMean)])
xlabel('Time after anaphase (hours)')
ylabel('Number of puncta')
subplot(2,2,3)
title(['Bottom ' num2str(lowerint) '% of cells with puncta, <= ' num2str(bottomMean)])
xlabel('Time after anaphase (hours)')
ylabel('CDK2 activity')
subplot(2,2,4)
title(['Bottom ' num2str(lowerint) '% of cells with puncta, <=' num2str(bottomMean)])
xlabel('Time after anaphase (hours)')
ylabel('Number of puncta')
for i = 1:(length(tops) + length(bottoms))
    %first things first, get the index for that cell sin the in the
    %original trace and puncta matrices
    %tops holds the index in the critSum matrix
    
    if i <= length(tops)
        index = ixmean(tops(i));
        plot = 0;
        %color = 'r';
    elseif i <= (length(bottoms) + length(tops))
        index = ixmean(bottoms(i-length(tops)));
        plot = 2;
        %color = 'b';
    end
    
    %index will ultimately be used to point towards one of three matrices:
    %incdata, lowdata, emergdata and their corresponding puncta matrices.
    %As such, if the index is less than or equal to the size of inc, then
    %use it to pull from inc. If bigger than inc but less than or equal to
    %the sum of inc and low, pull from low, etc.
    %However, this makes a problem, as index is for all data together, so
    %to deal with this, if index > = the number of inc traces, then we'll
    %need to subtract that number, and so on.
    if i <= (length(tops) + length(bottoms))
        if index <= size(incdatatotal_edit,1)
            time = alignedtime1_edit;
            CDK2 = incdatatotal_edit(index,:);
            puncta = incpuncta_datatotal_edit(index,:);
            color = 'b';
            if i <= length(tops)
                inc_hipuncta = inc_hipuncta+1;
            end
            if i > length(tops)
                inc_lowpuncta = inc_lowpuncta+1;
            end
        elseif index > size(incdatatotal_edit,1) && index <= (size(incdatatotal_edit,1) + size(lowdatatotal_edit,1))
            index = index - size(incdatatotal_edit,1);
            time = alignedtime2_edit;
            CDK2 = lowdatatotal_edit(index,:);
            puncta = lowpuncta_datatotal_edit(index,:);
            color = 'r';
            if i <= length(tops)
                low_hipuncta = low_hipuncta+1;
            end
            if i > length(tops)
                low_lowpuncta = low_lowpuncta+1;
            end
        elseif index > (size(incdatatotal_edit,1) + size(lowdatatotal_edit,1))
            index = index - (size(incdatatotal_edit,1) + size(lowdatatotal_edit,1));
            time = alignedtime3_edit;
            CDK2 = emergdatatotal_edit(index,:);
            puncta = emergpuncta_datatotal_edit(index,:);
            color = 'g';
            if i <= length(tops)
                emerg_hipuncta = emerg_hipuncta+1;
            end
            if i > length(tops)
                emerg_lowpuncta = emerg_lowpuncta+1;
            end
            
        end
        
        
        subplot(2,2,1+plot)
        line(time,CDK2,'color',color)
        xlim([-25 20])
        ylim([0.2 2])
        
        hold on;
        subplot(2,2,2+plot)
        line(time,puncta,'color',color)
        xlim([-25 20])
        ylim([0 10])
        hold on;
    end
end
fprintf('\nMEAN')
fprintf(['\n Number of CDK2inc in top 95th percentile: ' num2str(inc_hipuncta)])
fprintf(['\n Number of CDK2low in top 95th percentile: ' num2str(low_hipuncta)])
fprintf(['\n Number of CDK2emerg in top 95th percentile: ' num2str(emerg_hipuncta)])

fprintf(['\n Number of CDK2inc in bottom 5th percentile: ' num2str(inc_lowpuncta)])
fprintf(['\n Number of CDK2low in bottom 5th percentile: ' num2str(low_lowpuncta)])
fprintf(['\n Number of CDK2emerg in bottom 5th percentile: ' num2str(emerg_lowpuncta)])

fprintf(['\n Percentage of top percentile that''s CDK2inc: ' num2str(inc_hipuncta/(inc_hipuncta + low_hipuncta + emerg_hipuncta)) '%'])
fprintf(['\n Percentage of bottom percentile that''s CDK2inc: ' num2str(inc_lowpuncta/(inc_lowpuncta + low_lowpuncta + emerg_lowpuncta)) '%'])

%
topMed = prctile(critMedval,upperint);
bottomMed = prctile(critMedval,lowerint);
%topMed = 1;
%bottomMed = 0;

tops = find(critMedval >= topMed);
bottoms = find(critMedval <= bottomMed);
inc_hipuncta = 0; inc_lowpuncta = 0;
low_hipuncta = 0; low_lowpuncta = 0;
emerg_hipuncta = 0; emerg_lowpuncta = 0;
figure
suptitle(['Grouped according to median number of puncta in interval [' num2str(interval(1)) ':' num2str(interval(end)) ']'])
subplot(2,2,1)
title(['Top ' num2str(upperint) '% of cells with puncta, >= ' num2str(topMed)])
xlabel('Time after anaphase (hours)')
ylabel('CDK2 activity')
subplot(2,2,2)
title(['Top ' num2str(upperint) '% of cells with puncta, >= ' num2str(topMed)])
xlabel('Time after anaphase (hours)')
ylabel('Number of puncta')
subplot(2,2,3)
title(['Bottom ' num2str(lowerint) '% of cells with puncta, <= ' num2str(bottomMed)])
xlabel('Time after anaphase (hours)')
ylabel('CDK2 activity')
subplot(2,2,4)
title(['Bottom ' num2str(lowerint) '% of cells with puncta, <=' num2str(bottomMed)])
xlabel('Time after anaphase (hours)')
ylabel('Number of puncta')
for i = 1:(length(tops) + length(bottoms))
    %first things first, get the index for that cell sin the in the
    %original trace and puncta matrices
    %tops holds the index in the critSum matrix
    
    if i <= length(tops)
        index = ixmed(tops(i));
        plot = 0;
        %color = 'r';
    elseif i <= (length(bottoms) + length(tops))
        index = ixmed(bottoms(i-length(tops)));
        plot = 2;
        %color = 'b';
    end
    
    %index will ultimately be used to point towards one of three matrices:
    %incdata, lowdata, emergdata and their corresponding puncta matrices.
    %As such, if the index is less than or equal to the size of inc, then
    %use it to pull from inc. If bigger than inc but less than or equal to
    %the sum of inc and low, pull from low, etc.
    %However, this makes a problem, as index is for all data together, so
    %to deal with this, if index > = the number of inc traces, then we'll
    %need to subtract that number, and so on.
    if i <= (length(tops) + length(bottoms))
        if index <= size(incdatatotal_edit,1)
            time = alignedtime1_edit;
            CDK2 = incdatatotal_edit(index,:);
            puncta = incpuncta_datatotal_edit(index,:);
            color = 'b';
            if i <= length(tops)
                inc_hipuncta = inc_hipuncta+1;
            end
            if i > length(tops)
                inc_lowpuncta = inc_lowpuncta+1;
            end
        elseif index > size(incdatatotal_edit,1) && index <= (size(incdatatotal_edit,1) + size(lowdatatotal_edit,1))
            index = index - size(incdatatotal_edit,1);
            time = alignedtime2_edit;
            CDK2 = lowdatatotal_edit(index,:);
            puncta = lowpuncta_datatotal_edit(index,:);
            color = 'r';
            if i <= length(tops)
                low_hipuncta = low_hipuncta+1;
            end
            if i > length(tops)
                low_lowpuncta = low_lowpuncta+1;
            end
        elseif index > (size(incdatatotal_edit,1) + size(lowdatatotal_edit,1))
            index = index - (size(incdatatotal_edit,1) + size(lowdatatotal_edit,1));
            time = alignedtime3_edit;
            CDK2 = emergdatatotal_edit(index,:);
            puncta = emergpuncta_datatotal_edit(index,:);
            color = 'g';
            if i <= length(tops)
                emerg_hipuncta = emerg_hipuncta+1;
            end
            if i > length(tops)
                emerg_lowpuncta = emerg_lowpuncta+1;
            end
            
        end
        
        
        subplot(2,2,1+plot)
        line(time,CDK2,'color',color)
        xlim([-25 20])
        ylim([0.2 2])
        
        hold on;
        subplot(2,2,2+plot)
        line(time,puncta,'color',color)
        xlim([-25 20])
        ylim([0 10])
        hold on;
    end
end
fprintf('\nMEDIAN')
fprintf(['\n Number of CDK2inc in top 95th percentile: ' num2str(inc_hipuncta)])
fprintf(['\n Number of CDK2low in top 95th percentile: ' num2str(low_hipuncta)])
fprintf(['\n Number of CDK2emerg in top 95th percentile: ' num2str(emerg_hipuncta)])

fprintf(['\n Number of CDK2inc in bottom 5th percentile: ' num2str(inc_lowpuncta)])
fprintf(['\n Number of CDK2low in bottom 5th percentile: ' num2str(low_lowpuncta)])
fprintf(['\n Number of CDK2emerg in bottom 5th percentile: ' num2str(emerg_lowpuncta)])

fprintf(['\n Percentage of top percentile that''s CDK2inc: ' num2str(inc_hipuncta/(inc_hipuncta + low_hipuncta + emerg_hipuncta)) '%'])
fprintf(['\n Percentage of bottom percentile that''s CDK2inc: ' num2str(inc_lowpuncta/(inc_lowpuncta + low_lowpuncta + emerg_lowpuncta)) '%'])

topStd = prctile(critStdval,upperint);
bottomStd = prctile(critStdval,lowerint);
%topMed = 1;
%bottomMed = 0;

tops = find(critStdval >= topMed);
bottoms = find(critStdval <= bottomMed);
inc_hipuncta = 0; inc_lowpuncta = 0;
low_hipuncta = 0; low_lowpuncta = 0;
emerg_hipuncta = 0; emerg_lowpuncta = 0;
figure
suptitle(['Grouped according to standard deviation of the number of puncta in interval [' num2str(interval(1)) ':' num2str(interval(end)) ']'])
subplot(2,2,1)
title(['Top 95% of cells with puncta, >= ' num2str(topMed)])
xlabel('Time after anaphase (hours)')
ylabel('CDK2 activity')
subplot(2,2,2)
title(['Top 95% of cells with puncta, >= ' num2str(topMed)])
xlabel('Time after anaphase (hours)')
ylabel('Number of puncta')
subplot(2,2,3)
title(['Bottom 5% of cells with puncta, <= ' num2str(bottomMed)])
xlabel('Time after anaphase (hours)')
ylabel('CDK2 activity')
subplot(2,2,4)
title(['Bottom 5% of cells with puncta, <=' num2str(bottomMed)])
xlabel('Time after anaphase (hours)')
ylabel('Number of puncta')
for i = 1:(length(tops) + length(bottoms))
    %first things first, get the index for that cell sin the in the
    %original trace and puncta matrices
    %tops holds the index in the critSum matrix
    
    if i <= length(tops)
        index = ixstd(tops(i));
        plot = 0;
        %color = 'r';
    elseif i <= (length(bottoms) + length(tops))
        index = ixstd(bottoms(i-length(tops)));
        plot = 2;
        %color = 'b';
    end
    
    %index will ultimately be used to point towards one of three matrices:
    %incdata, lowdata, emergdata and their corresponding puncta matrices.
    %As such, if the index is less than or equal to the size of inc, then
    %use it to pull from inc. If bigger than inc but less than or equal to
    %the sum of inc and low, pull from low, etc.
    %However, this makes a problem, as index is for all data together, so
    %to deal with this, if index > = the number of inc traces, then we'll
    %need to subtract that number, and so on.
    if i <= (length(tops) + length(bottoms))
        if index <= size(incdatatotal,1)
            time = alignedtime1;
            CDK2 = incdatatotal(index,:);
            puncta = incpuncta_datatotal(index,:);
            color = 'b';
            if i <= length(tops)
                inc_hipuncta = inc_hipuncta+1;
            end
            if i > length(tops)
                inc_lowpuncta = inc_lowpuncta+1;
            end
        elseif index > size(incdatatotal,1) && index <= (size(incdatatotal,1) + size(lowdatatotal,1))
            index = index - size(incdatatotal,1);
            time = alignedtime2;
            CDK2 = lowdatatotal(index,:);
            puncta = lowpuncta_datatotal(index,:);
            color = 'r';
            if i <= length(tops)
                low_hipuncta = low_hipuncta+1;
            end
            if i > length(tops)
                low_lowpuncta = low_lowpuncta+1;
            end
        elseif index > (size(incdatatotal,1) + size(lowdatatotal,1))
            index = index - (size(incdatatotal,1) + size(lowdatatotal,1));
            time = alignedtime3;
            CDK2 = emergdatatotal(index,:);
            puncta = emergpuncta_datatotal(index,:);
            color = 'g';
            if i <= length(tops)
                emerg_hipuncta = emerg_hipuncta+1;
            end
            if i > length(tops)
                emerg_lowpuncta = emerg_lowpuncta+1;
            end
            
        end
        
        
        subplot(2,2,1+plot)
        line(time,CDK2,'color',color)
        xlim([-25 20])
        ylim([0.2 2])
        
        hold on;
        subplot(2,2,2+plot)
        line(time,puncta,'color',color)
        xlim([-25 20])
        ylim([0 10])
        hold on;
    end
end
fprintf('\nSTANDARD DEVIATION')
fprintf(['\n Number of CDK2inc in top 95th percentile: ' num2str(inc_hipuncta)])
fprintf(['\n Number of CDK2low in top 95th percentile: ' num2str(low_hipuncta)])
fprintf(['\n Number of CDK2emerg in top 95th percentile: ' num2str(emerg_hipuncta)])

fprintf(['\n Number of CDK2inc in bottom 5th percentile: ' num2str(inc_lowpuncta)])
fprintf(['\n Number of CDK2low in bottom 5th percentile: ' num2str(low_lowpuncta)])
fprintf(['\n Number of CDK2emerg in bottom 5th percentile: ' num2str(emerg_lowpuncta)])

fprintf(['\n Percentage of top percentile that''s CDK2inc: ' num2str(inc_hipuncta/(inc_hipuncta + low_hipuncta + emerg_hipuncta)) '%'])
fprintf(['\n Percentage of bottom percentile that''s CDK2inc: ' num2str(inc_lowpuncta/(inc_lowpuncta + low_lowpuncta + emerg_lowpuncta)) '%'])
fprintf('\n')