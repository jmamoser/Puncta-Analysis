close all;
clear all;
%%% filepaths %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
projectpath='\\spencerstorage.int.colorado.edu\LabShare\IXMicroImages-goodNames\Mansi\'; %path where the processed data will be stored; create a folder called "Data" here
imagepath=  '\\spencerstorage.int.colorado.edu\LabShare\IXMicroImages-goodNames\Mansi\';  %where the images are stored, generally
experimentpath='MA61-20160419-mChyBP1-drugs_1556\'; 
% experimentpath='MA56-20160211-VitC-NAC-pilot3colorMCF10A_782\';
datadir=([projectpath,experimentpath,'Data\']);

%%% analysis settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
motheroption=0; %0:no gating 1:mothers 2:no mothers
daughteroption=1; %0:no gating 1:daughters 2:no daughters
framesaftermitosis = 20; % categorizes based on puncta existence between POM and these many frames after mitosis.
framenum = 120; % number of frames of the entire movie to look at.
quiescentanalysis=0;
numoftraces = 400;
numoftraces2 = 100;
if quiescentanalysis
    motheroption=2; daughteroption=2;
end
IFoption=0; %0: No IF data 1: IF data
framesperhr=5;

%%% pool data from multiple wells %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rowmat=[1:2];
colmat=[1];
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

% %%% categorize traces %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% numtraces=size(tracestats,1);
% calldelay=25; %must be >= minlengthtrace
% callthresh=0.6;
% categorization=ones(numtraces,1)*NaN;
% valstore=ones(numtraces,1)*NaN;
% 
% for cc=1:numtraces
%    calltime=tracestats(cc,1)+calldelay-1; % #frames after mitosis
%    if(calltime < size(tracesCDK2,2)-calldelay)
%         categorization(cc)=tracesCDK2(cc,calltime)>callthresh;
%         valstore(cc)=min(tracesCDK2(cc,tracestats(cc,2)-2:tracestats(cc,2)));
%     end
% end


% Cdk2inc=find(categorization==1);
% Cdk2low=find(categorization==0);
%%% Categorize based on puncta within 5 hours of mitosis



%%% align traces %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
frames=1:120;
cdk2incframesmat=repmat(frames, size(tracesCDK2inc,1), 1);
cdk2incmitosismat=repmat(CDK2incPOM, length(frames), 1)';
cdk2incshiftedframes=cdk2incframesmat-cdk2incmitosismat;
cdk2incshiftedframes=cdk2incshiftedframes./framesperhr ;
tracesCDK2inc_24 = tracesCDK2inc(:,1:120);

cdk2lowframesmat=repmat(frames, size(tracesCDK2low,1), 1);
cdk2lowmitosismat=repmat(CDK2lowPOM, length(frames), 1)';
cdk2lowshiftedframes=cdk2lowframesmat-cdk2lowmitosismat;
cdk2lowshiftedframes=cdk2lowshiftedframes./framesperhr ;
tracesCDK2low_24 = tracesCDK2low(:,1:120);

cdk2emergframesmat=repmat(frames, size(tracesCDK2emerg,1), 1);
cdk2emergmitosismat=repmat(CDK2emergPOM, length(frames), 1)';
cdk2emergshiftedframes=cdk2emergframesmat-cdk2emergmitosismat;
cdk2emergshiftedframes=cdk2emergshiftedframes./framesperhr ;
tracesCDK2emerg_24 = tracesCDK2emerg(:,1:120);


% [alignedtime,aligneddata]=aligntraces_4(tracesCDK2,POI,tracestats,motherstats,daughteroption);
% alignedtime=alignedtime/framesperhr;

%justin addition
% drugframe = 1;
% 
% drugPOI = POI-drugframe; %number of frames between drug addition and mitosis. If the drug was added BEFORE mitosis, drugPOI > 0; if the drug was added after mitosis, drugPOI <0
% tzeroframe = find(alignedtime==0); %zero on the x axis in the aligned mitosis plots 
% aligneddrug = tzeroframe - drugPOI; %index of aligned time at which drug was added
% smoothingsize=5;
% ymin=0.2; ymax=3;
% figure (300), hold on;
% goodcdkcounter=0;
% goodtraceindex=[];
% inccounter = 0;
% lowcounter = 0;

%%% plot traces based on classification on with or without puncta %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numtraces = length(CDK2incPOM) + length(CDK2lowPOM) + length(CDK2emergPOM);

smoothingsize=5;
ymin=0.2; ymax=2.5;
figure(1), hold on;
figure(2), hold on;
traces_wo_puncta = 0;
traces_with_puncta = 0;
cdk2low_wo_puncta = 0;
cdk2low_with_puncta = 0;
cdk2emerg_wo_puncta = 0;
cdk2emerg_with_puncta = 0;
cdk2inc_with_puncta = 0;
cdk2inc_wo_puncta = 0;
total_traces_after_gating = 0;
indices = 1:numtraces;
%for pc=1:numtraces
withcount = 1; without_count = 1;
while withcount <= numoftraces || without_count <= numoftraces
     pc = randi(length(indices));
     if ((pc <= length(CDK2incPOM)) && (CDK2incPOM(pc)+ framesaftermitosis < framenum))
        total_traces_after_gating = total_traces_after_gating+1;
        
        if isempty(find(CDK2incPuncta(pc,CDK2incPOM(pc)+6 : CDK2incPOM(pc)+framesaftermitosis))) && withcount <= numoftraces% +6 stands for 6th frame after POM.
            
            figure(1)
                plotdata=tracesCDK2inc_24(pc,:);
                realframes=find(~isnan(tracesCDK2inc_24(pc,:)));
                plotdata(realframes)=smooth(tracesCDK2inc_24(pc,realframes),smoothingsize);
                tracecolor='b';
                cdk2inc_wo_puncta = cdk2inc_wo_puncta+1;
                traces_wo_puncta = traces_wo_puncta +1;
                withcount = withcount + 1;
                line(cdk2incshiftedframes(pc,:), plotdata, 'color', tracecolor)
        elseif without_count <= numoftraces
                figure(2)
                plotdata=tracesCDK2inc_24(pc,:);
                realframes=find(~isnan(tracesCDK2inc_24(pc,:)));
                plotdata(realframes)=smooth(tracesCDK2inc_24(pc,realframes),smoothingsize);
                tracecolor= 'b';
                cdk2inc_with_puncta = cdk2inc_with_puncta + 1;
                traces_with_puncta = traces_with_puncta +1;
                without_count = without_count + 1;
                line(cdk2incshiftedframes(pc,:), plotdata, 'color', tracecolor)
        end
        
        
    elseif ((pc > length(CDK2incPOM)) && (pc <= length(CDK2incPOM)+ length(CDK2lowPOM)) && (CDK2lowPOM(pc-length(CDK2incPOM))+ framesaftermitosis < framenum))
            ei = pc-length(CDK2incPOM);
            total_traces_after_gating = total_traces_after_gating+1;
        if isempty(find(CDK2lowPuncta(ei,CDK2lowPOM(ei)+6 : CDK2lowPOM(ei)+framesaftermitosis))) && withcount <= numoftraces% +6 stands for 6th frame after POM.if isempty(find(CDK2lowPuncta(ei,(CDK2lowPOM(ei)+6)):(CDK2lowPOM(ei)+framesaftermitosis)))% +6 stands for 6th frame after POM.
                figure(1)
                plotdata=tracesCDK2low_24(ei,:);
                realframes=find(~isnan(tracesCDK2low_24(ei,:)));
                plotdata(realframes)=smooth(tracesCDK2low_24(ei,realframes),smoothingsize);
                tracecolor='r';
                cdk2low_wo_puncta = cdk2low_wo_puncta+1;
                traces_wo_puncta = traces_wo_puncta +1;
                withcount = withcount + 1;
                line(cdk2lowshiftedframes(ei,:), plotdata, 'color', tracecolor)
        elseif without_count <= numoftraces
                figure(2)
                plotdata=tracesCDK2low_24(ei,:);
                realframes=find(~isnan(tracesCDK2low_24(ei,:)));
                plotdata(realframes)=smooth(tracesCDK2low_24(ei,realframes),smoothingsize);
                tracecolor= 'r';
                cdk2low_with_puncta = cdk2low_with_puncta + 1;
                traces_with_puncta = traces_with_puncta +1;
                without_count = without_count + 1;
                line(cdk2lowshiftedframes(ei,:), plotdata, 'color', tracecolor)
        end
        
        
   elseif ((pc > length(CDK2incPOM)+ length(CDK2lowPOM)) && (CDK2emergPOM(pc-length(CDK2incPOM)-length(CDK2lowPOM))+ framesaftermitosis < framenum))
            ai = pc-length(CDK2incPOM)-length(CDK2lowPOM);
            total_traces_after_gating = total_traces_after_gating+1;
        if isempty(find(CDK2emergPuncta(ai,CDK2emergPOM(ai)+6 : CDK2emergPOM(ai)+framesaftermitosis))) && withcount <= numoftraces% +6 stands for 6th frame after POM.
                figure(1)
                plotdata=tracesCDK2emerg_24(ai,:);
                realframes=find(~isnan(tracesCDK2emerg_24(ai,:)));
                plotdata(realframes)=smooth(tracesCDK2emerg_24(ai,realframes),smoothingsize);
                tracecolor='g';
                cdk2emerg_wo_puncta = cdk2emerg_wo_puncta+1;
                traces_wo_puncta = traces_wo_puncta +1;
                withcount = withcount + 1;
                line(cdk2emergshiftedframes(ai,:), plotdata, 'color', tracecolor)
         elseif without_count <= numoftraces
                figure(2)
                plotdata=tracesCDK2emerg_24(ai,:);
                realframes=find(~isnan(tracesCDK2emerg_24(ai,:)));
                plotdata(realframes)=smooth(tracesCDK2emerg_24(ai,realframes),smoothingsize);
                tracecolor= 'g';
                cdk2emerg_with_puncta = cdk2emerg_with_puncta + 1;
                traces_with_puncta = traces_with_puncta +1;
                without_count = without_count + 1;
                line(cdk2emergshiftedframes(ai,:), plotdata, 'color', tracecolor)
        end
       
       
      
    end
   
end
xlim([-20 20])
ylabel('CDK2 Activity');
xlabel('Time(hrs) Relative to Anaphase');


pct_Cdk2low = ((cdk2low_wo_puncta+cdk2low_with_puncta)/total_traces_after_gating).*100;
pct_cdk2low_wo_puncta = (cdk2low_wo_puncta/traces_wo_puncta)*100;
pct_cdk2low_with_puncta = (cdk2low_with_puncta/traces_with_puncta)*100;

pct_cdk2emerg_wo_puncta = (cdk2emerg_wo_puncta/traces_wo_puncta)*100;
pct_cdk2emerg_with_puncta = (cdk2emerg_with_puncta/traces_with_puncta)*100;

pct_cdk2low_with_puncta_wrt_cdk2low = (cdk2low_with_puncta/(cdk2low_wo_puncta+cdk2low_with_puncta))*100;
pct_cdk2inc_with_puncta_wrt_cdk2inc = (cdk2inc_with_puncta/(cdk2inc_with_puncta+cdk2inc_wo_puncta))*100;


set(gcf,'color','w','PaperPosition',[0 0 8 6]);
figure(1)
title([num2str(rowmat) '\_' num2str(colmat),'  ,%CDK2low of Total: ' num2str(pct_Cdk2low),' # Traces w/o puncta: ', num2str(traces_wo_puncta),' %CDK2Low of Traces w/o puncta: ',num2str(pct_cdk2low_wo_puncta)], 'fontsize', 8);

%title(  ['Row: ' num2str(rowmat) '  Col: ' num2str(colmat) '       Total Traces: ' num2str(numtraces)  '    %CDK2low of Total Traces:  ' num2str(pct_Cdk2low)] , 'fontsize', 16)
%saveas(gcf, [datadir, 'Row_' num2str(rowmat) '_Col_' num2str(colmat) '_CDK2alignedToMitosisredblue.fig']);
figure(2)
title([num2str(rowmat) '\_' num2str(colmat),'  ,%CDK2low of Total: ' num2str(pct_Cdk2low),' # Traces with puncta: ', num2str(traces_with_puncta),' ,%CDK2Low of Traces with puncta: ',num2str(pct_cdk2low_with_puncta)], 'fontsize', 8);


%%% Plotting based on Classification 2:


%%% for aligning time
numframes = 120;
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

indices1 = [1:size(tracesCDK2inc_24,1)];
count = 1;
while count <= numoftraces2 
        i = randi(length(indices1));
        plotdata=tracesCDK2inc_24(i,:);
        realframes=find(~isnan(tracesCDK2inc_24(i,:)));
        plotdata(realframes)=smooth(tracesCDK2inc_24(i,realframes),smoothingsize);
        figure(3)
        subplot(2,1,1)
        line(cdk2incshiftedframes(i,:),plotdata);
        title('CDK2 inc');
        firstframe=find(~isnan(tracesCDK2inc_24(i,:)),1,'first');
        lastframe=find(~isnan(tracesCDK2inc_24(i,:)),1,'last');
        lengthleftinc(i)=CDK2incPOM(i)-firstframe+1;
        lengthrightinc(i)=lastframe-CDK2incPOM(i)+1;
        incpunctatraceleft(i,numframes-lengthleftinc(i)+1:numframes)=CDK2incPuncta(i,firstframe:CDK2incPOM(i));
        incpunctatraceright(i,1:lengthrightinc(i))=CDK2incPuncta(i,CDK2incPOM(i):lastframe);
        tracedataincleft(i,numframes-lengthleftinc(i)+1:numframes)=tracesCDK2inc_24(i,firstframe:CDK2incPOM(i));
        tracedataincright(i,1:lengthrightinc(i))=tracesCDK2inc_24(i,CDK2incPOM(i):lastframe);
        count = count+1;
        
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

indices2 = [1:size(tracesCDK2low_24,1)];
count = 1;
while count <= numoftraces2 
        i = randi(length(indices2));
        plotdata=tracesCDK2low_24(i,:);
        realframes=find(~isnan(tracesCDK2low_24(i,:)));
        plotdata(realframes)=smooth(tracesCDK2low_24(i,realframes),smoothingsize);
        figure(4)
        subplot(2,1,1)
        line(cdk2lowshiftedframes(i,:),plotdata);
        title('CDK2 low');
        firstframe=find(~isnan(tracesCDK2low_24(i,:)),1,'first');
        lastframe=find(~isnan(tracesCDK2low_24(i,:)),1,'last');
        lengthleftlow(i)=CDK2lowPOM(i)-firstframe+1;
        lengthrightlow(i)=lastframe-CDK2lowPOM(i)+1;
        lowpunctatraceleft(i,numframes-lengthleftlow(i)+1:numframes)=CDK2lowPuncta(i,firstframe:CDK2lowPOM(i));
        lowpunctatraceright(i,1:lengthrightlow(i))=CDK2lowPuncta(i,CDK2lowPOM(i):lastframe);
        tracedatalowleft(i,numframes-lengthleftlow(i)+1:numframes)=tracesCDK2low_24(i,firstframe:CDK2lowPOM(i));
        tracedatalowright(i,1:lengthrightlow(i))=tracesCDK2low_24(i,CDK2lowPOM(i):lastframe);
        count = count+1;
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



%for i=1:length(CDK2emergPOM)
indices3 = [1:size(tracesCDK2emerg_24,1)];
count = 1;
while count <= numoftraces2 
        i = randi(length(indices3));
        plotdata=tracesCDK2emerg_24(i,:);
        realframes=find(~isnan(tracesCDK2emerg_24(i,:)));
        plotdata(realframes)=smooth(tracesCDK2emerg_24(i,realframes),smoothingsize);
        figure(5)
        subplot(2,1,1)
        line(cdk2emergshiftedframes(i,:),plotdata);
        title('CDK2 emerg')
        firstframe=find(~isnan(tracesCDK2emerg_24(i,:)),1,'first');
        lastframe=find(~isnan(tracesCDK2emerg_24(i,:)),1,'last');
        lengthleftemerg(i)=CDK2emergPOM(i)-firstframe+1;
        lengthrightemerg(i)=lastframe-CDK2emergPOM(i)+1;
        emergpunctatraceleft(i,numframes-lengthleftemerg(i)+1:numframes)=CDK2emergPuncta(i,firstframe:CDK2emergPOM(i));
        emergpunctatraceright(i,1:lengthrightemerg(i))=CDK2emergPuncta(i,CDK2emergPOM(i):lastframe);
        tracedataemergleft(i,numframes-lengthleftemerg(i)+1:numframes)=tracesCDK2emerg_24(i,firstframe:CDK2emergPOM(i));
        tracedataemergright(i,1:lengthrightemerg(i))=tracesCDK2emerg_24(i,CDK2emergPOM(i):lastframe);
        count = count+1;
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

alignedtime1 = alignedtime1./framesperhr;
alignedtime2 = alignedtime2./framesperhr;
alignedtime3 = alignedtime3./framesperhr;

plotdata2_mean = nanmean(incpuncta_datatotal);
plotdata3_mean = nanmean(lowpuncta_datatotal);
plotdata4_mean = nanmean(emergpuncta_datatotal);


for i= 1: size(incdatatotal,2)
    error1(i) = nanstd(incpuncta_datatotal(:,i))/sqrt(length(~isnan(incdatatotal(:,i))));
      
end

for i = 1: size(lowdatatotal,2)
    error2(i) = nanstd(lowpuncta_datatotal(:,i))/sqrt(length(~isnan(lowdatatotal(:,i))));
end

for i = 1: size(emergdatatotal,2)
    error3(i) = nanstd(emergpuncta_datatotal(:,i))/sqrt(length(~isnan(emergdatatotal(:,i))));
end
% 
figure(3)
subplot(2,1,2)
errorbar(alignedtime1,plotdata2_mean,error1);


figure(4)
subplot(2,1,2)
errorbar(alignedtime2,plotdata3_mean,error2);


figure(5)
subplot(2,1,2)
errorbar(alignedtime3,plotdata4_mean,error3);


figure(6)
errorbar(alignedtime1,plotdata2_mean,error1,'b');hold on;
errorbar(alignedtime2,plotdata3_mean,error2,'r');hold on;
errorbar(alignedtime3,plotdata4_mean,error3,'g');hold off;
%         subplot(2,1,2)

%%% record number of traces %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('%0.0f traces\n',numtraces);