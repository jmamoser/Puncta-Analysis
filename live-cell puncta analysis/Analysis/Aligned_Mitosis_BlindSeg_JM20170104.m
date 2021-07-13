close all;
%clear all;
%%% filepaths %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
projectpath='Y:\IXMicroImages-goodNames\Justin\'; %path where the processed data will be stored; create a folder called "Data" here
imagepath=  'Y:\IXMicroImages-goodNames\Justin\';  %where the images are stored, generally
experimentpath='MANSI_ANALYSIS_20170104\'; 
% experimentpath='MA56-20160211-VitC-NAC-pilot3colorMCF10A_782\';
datadir=([projectpath,experimentpath,'Data2\']);

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

%%% pool data from multiple wells %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rowmat=[1:8];
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
frames=1:120; INC = 1; LOW = 0; EMERG = 2;

cdk2incframesmat=repmat(frames, size(tracesCDK2inc,1), 1);
cdk2incmitosismat=repmat(CDK2incPOM, length(frames), 1)';
cdk2incshiftedframes=cdk2incframesmat-cdk2incmitosismat;
cdk2incshiftedframes=cdk2incshiftedframes./framesperhr ;
tracesCDK2inc_24 = tracesCDK2inc(:,1:120);
punctaCDK2inc_24 = CDK2incPuncta(:,1:120); %JM
INCclass(1:size(punctaCDK2inc_24),1) = INC;

figure, hold on; suptitle('CDK2inc cells')
for i = 1:size(tracesCDK2inc_24,1)
    subplot(2,1,1)
    line(cdk2incshiftedframes(i,:),tracesCDK2inc_24(i,:),'color','b')
    subplot(2,1,2)
    line(cdk2incshiftedframes(i,:),punctaCDK2inc_24(i,:),'color','b')
end
subplot(2,1,1)
xlabel('Time since anaphase (hrs)')
ylabel('CDK2 activity')
subplot(2,1,2)
xlabel('Time since anaphase (hrs)')
ylabel('Number of Puncta')

cdk2lowframesmat=repmat(frames, size(tracesCDK2low,1), 1);
cdk2lowmitosismat=repmat(CDK2lowPOM, length(frames), 1)';
cdk2lowshiftedframes=cdk2lowframesmat-cdk2lowmitosismat;
cdk2lowshiftedframes=cdk2lowshiftedframes./framesperhr ;
tracesCDK2low_24 = tracesCDK2low(:,1:120);
punctaCDK2low_24 = CDK2lowPuncta(:,1:120); %JM
LOWclass(1:size(punctaCDK2low_24),1) = LOW;

figure, hold on; suptitle('CDK2low cells')
for i = 1:size(tracesCDK2low_24,1)
    subplot(2,1,1)
    line(cdk2lowshiftedframes(i,:),tracesCDK2low_24(i,:),'color','r')
    subplot(2,1,2)
    line(cdk2lowshiftedframes(i,:),punctaCDK2low_24(i,:),'color','r')
end
subplot(2,1,1)
xlabel('Time since anaphase (hrs)')
ylabel('CDK2 activity')
subplot(2,1,2)
xlabel('Time since anaphase (hrs)')
ylabel('Number of Puncta')

cdk2emergframesmat=repmat(frames, size(tracesCDK2emerg,1), 1);
cdk2emergmitosismat=repmat(CDK2emergPOM, length(frames), 1)';
cdk2emergshiftedframes=cdk2emergframesmat-cdk2emergmitosismat;
cdk2emergshiftedframes=cdk2emergshiftedframes./framesperhr ;
tracesCDK2emerg_24 = tracesCDK2emerg(:,1:120);
punctaCDK2emerg_24 = CDK2emergPuncta(:,1:120); %JM
EMERGclass(1:size(punctaCDK2emerg_24),1) = EMERG;

figure, hold on; suptitle('CDK2emerge cells')
for i = 1:size(tracesCDK2emerg_24,1)
    subplot(2,1,1)
    line(cdk2emergshiftedframes(i,:),tracesCDK2emerg_24(i,:),'color','g')
    subplot(2,1,2)
    line(cdk2emergshiftedframes(i,:),punctaCDK2emerg_24(i,:),'color','g')
end

subplot(2,1,1)
xlabel('Time since anaphase (hrs)')
ylabel('CDK2 activity')
subplot(2,1,2)
xlabel('Time since anaphase (hrs)')
ylabel('Number of Puncta')

%create a safe place to store all this stuff
timestore       = [cdk2incshiftedframes; cdk2lowshiftedframes; cdk2emergshiftedframes];
CDK2tracestore  = [tracesCDK2inc_24; tracesCDK2low_24; tracesCDK2emerg_24];
POMstore        = [cdk2incmitosismat; cdk2lowmitosismat; cdk2emergmitosismat];
punctastore     = [punctaCDK2inc_24; punctaCDK2low_24; punctaCDK2emerg_24];
classstore      = [INCclass; LOWclass; EMERGclass];

permutation_index = randperm(size(timestore,1),size(timestore,1))';

%here, sh_ will indicate that the "safe" data has been shuffled to prevent
%bias...
sh_POMstore = POMstore(permutation_index,:);
sh_CDK2tracestore = CDK2tracestore(permutation_index,:);
sh_punctastore = punctastore(permutation_index,:);
sh_timestore = timestore(permutation_index,:);
sh_classstore = classstore(permutation_index,:);

%now, sort the puncta in descending order at a given interval
%interval = framesaftermitosis; %right now interval is one value, 20, so that we're looking at a single timepoint...
%interval = 1:20;
num_examine = 0.05*size(sh_POMstore,1);
    
for i = 1:size(sh_POMstore,1)
    pom = find(sh_timestore(i,:) == 0);
    if pom + max(interval) < max(frames)-1
        %rel_puncta(i,1) = sum(sh_punctastore(i,sh_POMstore(i,1) + interval)); %could also just remove the sum if you are only looking at a single timepoint
    rel_puncta(i,1) = mean(sh_punctastore(i,pom + interval)); %could also just remove the sum if you are only looking at a single timepoint
    else
        rel_puncta(i,1) = NaN;
    end
end

[critval,ix] = sort(rel_puncta,'descend'); %left critval in for use as a diagnostic, nothing more. We're interested in ix at this point
valid = find(~isnan(critval));

sorted_puncta = sh_punctastore(valid,:);
sorted_CDK2traces = sh_CDK2tracestore(valid,:);
sorted_time = sh_timestore(valid,:);
sorted_POM = sh_POMstore(valid,:);
sorted_class = sh_classstore(valid,:);

% figure, hold on; suptitle('All cells plotted together') %just to verify that we haven't thrown all our cells out
% for i = 1:size(sorted_CDK2traces,1)
%     if sorted_class(i) == 1
%         color = 'b';
%     elseif sorted_class(i) == 0
%         color = 'r';       
%     elseif sorted_class(i) == 2
%         color = 'g';
%     end
%     line(sorted_time(i,:),sorted_CDK2traces(i,:),'color',color) 
% end
% xlabel('Time since anaphase (hrs)')
% ylabel('CDK2 activity')

% figure, scatter3(valid,sorted_class,critval(valid))
% xlabel('index')
% ylabel('Manual Classificaton')
% zlabel('Number of Puncta')
% 
% figure, dscatter(sorted_class,critval(valid))

%now that we've sorted everything, pull out the top and bottom 5% 
top_puncta = sorted_puncta(1:num_examine,:);
top_CDK2   = sorted_CDK2traces(1:num_examine,:);
top_POM    = sorted_POM(1:num_examine,:);
top_time   = sorted_time(1:num_examine,:);
top_class  = sorted_class(1:num_examine,:);

bottom_puncta = sorted_puncta(size(sorted_POM,1)-num_examine+1:end,:);
bottom_CDK2   = sorted_CDK2traces(size(sorted_POM,1)-num_examine+1:end,:);
bottom_POM    = sorted_POM(size(sorted_POM,1)-num_examine+1:end,:);
bottom_time   = sorted_time(size(sorted_POM,1)-num_examine+1:end,:);
bottom_class  = sorted_class(size(sorted_POM,1)-num_examine+1:end,:);
%Okay, so we have the top and bottom 5% of entries; now we want to plot
%them separately to see if the tops are actually CDK2low and the bottoms
%are actually CDK2inc
%start with the putatively CDK2low, since they'll be the most obvious/easy
%to verify
ctime = [top_time; bottom_time];
c_CDK2 = [top_CDK2; bottom_CDK2];
c_class = [top_class; bottom_class];
c_puncta = [top_puncta; bottom_puncta];

a= figure(); hold on; title(['Colored according to original, manual classification; Interval: ' num2str(interval(1)) ':' num2str(interval(end))])
b=figure(); hold on; title(['High mean puncta cells; Interval: ' num2str(interval(1)) ':' num2str(interval(end))])
c=figure(); hold on; title(['Low mean puncta cells; Interval: ' num2str(interval(1)) ':' num2str(interval(end))])
for i = 1:size(c_class,1)
    if c_class(i) == 1
        color = 'b';
    elseif c_class(i) == 0
        color = 'r';
    else
        color = 'g';
    end
    
    figure(a)
    line(ctime(i,:),c_CDK2(i,:),'color',color);
    
%     if i <= num_examine
%         figure(b)
%         line(ctime(i,:),c_CDK2(i,:),'color','r')
%     else
%         figure(c)
%         line(ctime(i,:),c_CDK2(i,:),'color','b')
%     end
end
for i = 1:size(bottom_class,1)
   figure(b)
   line(top_time(i,:),top_CDK2(i,:),'color','r')
       pom = find(top_time(i,:) == 0);
      %figure(c)
       scatter(top_time(i,pom+interval(1)), top_CDK2(i,pom+interval(1)),'k','filled')
   
   %figure(c)
   %line(bottom_time(i,:),bottom_CDK2(i,:),'color','b')
end


for i = 1:size(sh_POMstore,1)
        pom = find(sh_timestore(i,:) == 0);
    if pom + max(interval) <= max(frames)
        rel_puncta(i,1) = sum(sh_punctastore(i,pom + interval)); %could also just remove the sum if you are only looking at a single timepoint
    %rel_puncta(i,1) = mean(sh_punctastore(i,sh_POMstore(i,1) + interval)); %could also just remove the sum if you are only looking at a single timepoint
    else
        rel_puncta(i,1) = NaN;
    end
end

[critval,ix] = sort(rel_puncta,'descend'); %left critval in for use as a diagnostic, nothing more. We're interested in ix at this point
valid = find(~isnan(critval));

sorted_puncta = sh_punctastore(valid,:);
sorted_CDK2traces = sh_CDK2tracestore(valid,:);
sorted_time = sh_timestore(valid,:);
sorted_POM = sh_POMstore(valid,:);
sorted_class = sh_classstore(valid,:);

% figure, hold on; suptitle('All cells plotted together') %just to verify that we haven't thrown all our cells out
% for i = 1:size(sorted_CDK2traces,1)
%     if sorted_class(i) == 1
%         color = 'b';
%     elseif sorted_class(i) == 0
%         color = 'r';       
%     elseif sorted_class(i) == 2
%         color = 'g';
%     end
%     line(sorted_time(i,:),sorted_CDK2traces(i,:),'color',color) 
% end
% xlabel('Time since anaphase (hrs)')
% ylabel('CDK2 activity')

% figure, scatter3(valid,sorted_class,critval(valid))
% xlabel('index')
% ylabel('Manual Classificaton')
% zlabel('Number of Puncta')
% 
% figure, dscatter(sorted_class,critval(valid))

%now that we've sorted everything, pull out the top and bottom 5% 
top_puncta = sorted_puncta(1:num_examine,:);
top_CDK2   = sorted_CDK2traces(1:num_examine,:);
top_POM    = sorted_POM(1:num_examine,:);
top_time   = sorted_time(1:num_examine,:);
top_class  = sorted_class(1:num_examine,:);

bottom_puncta = sorted_puncta(size(sorted_POM,1)-num_examine+1:end,:);
bottom_CDK2   = sorted_CDK2traces(size(sorted_POM,1)-num_examine+1:end,:);
bottom_POM    = sorted_POM(size(sorted_POM,1)-num_examine+1:end,:);
bottom_time   = sorted_time(size(sorted_POM,1)-num_examine+1:end,:);
bottom_class  = sorted_class(size(sorted_POM,1)-num_examine+1:end,:);
%Okay, so we have the top and bottom 5% of entries; now we want to plot
%them separately to see if the tops are actually CDK2low and the bottoms
%are actually CDK2inc
%start with the putatively CDK2low, since they'll be the most obvious/easy
%to verify
ctime = [top_time; bottom_time];
c_CDK2 = [top_CDK2; bottom_CDK2];
c_class = [top_class; bottom_class];
c_puncta = [top_puncta; bottom_puncta];

a= figure(); hold on; title(['Colored according to original, manual classification; Interval: ' num2str(interval(1)) ':' num2str(interval(end))])
b=figure(); hold on; title(['High sum puncta cells; Interval: ' num2str(interval(1)) ':' num2str(interval(end))])
c=figure(); hold on; title(['Low sum puncta cells; Interval: ' num2str(interval(1)) ':' num2str(interval(end))])
for i = 1:size(c_class,1)
    if c_class(i) == 1
        color = 'b';
    elseif c_class(i) == 0
        color = 'r';
    else
        color = 'g';
    end
    
    figure(a)
    line(ctime(i,:),c_CDK2(i,:),'color',color);
    
%     if i <= num_examine
%         figure(b)
%         line(ctime(i,:),c_CDK2(i,:),'color','r')
%     else
%         figure(c)
%         line(ctime(i,:),c_CDK2(i,:),'color','b')
%     end
end
for i = 1:size(bottom_class,1)
   figure(b)
   line(top_time(i,:),top_CDK2(i,:),'color','r')
   figure(c)
   line(bottom_time(i,:),bottom_CDK2(i,:),'color','b')
end

%% Tomorrow (20170105) try plotting the cells with the same critval together; see if that has any effect...


% if numel(interval) == 1
%     
%     %the following lines extract the indices of top and bottom 5% of the
%     %sorted cells, according to the number of puncta at time interval
% 
%     [critval, ix]      = sort(rel_puncta,'descend');
%          
% elseif numel(interval) > 1
%     for i = 1:size(sh_POMstore,1)
%         sumarray(i,1) = sum(sh_punctastore(i,interval));
%     end
%     [critval,ix] = sort(sumarray,'descend');
%     
% end
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
for pc=1:numtraces
    if ((pc <= length(CDK2incPOM)) && (CDK2incPOM(pc)+ framesaftermitosis < framenum))
        total_traces_after_gating = total_traces_after_gating+1;
        
        if isempty(find(CDK2incPuncta(pc,CDK2incPOM(pc)+6 : CDK2incPOM(pc)+framesaftermitosis)))% +6 stands for 6th frame after POM.
            
            figure(1)
                plotdata=tracesCDK2inc_24(pc,:);
                realframes=find(~isnan(tracesCDK2inc_24(pc,:)));
                plotdata(realframes)=smooth(tracesCDK2inc_24(pc,realframes),smoothingsize);
                tracecolor='b';
                cdk2inc_wo_puncta = cdk2inc_wo_puncta+1;
                traces_wo_puncta = traces_wo_puncta +1;
         else
                figure(2)
                plotdata=tracesCDK2inc_24(pc,:);
                realframes=find(~isnan(tracesCDK2inc_24(pc,:)));
                plotdata(realframes)=smooth(tracesCDK2inc_24(pc,realframes),smoothingsize);
                tracecolor= 'b';
                cdk2inc_with_puncta = cdk2inc_with_puncta + 1;
                traces_with_puncta = traces_with_puncta +1;
        end
        line(cdk2incshiftedframes(pc,:), plotdata, 'color', tracecolor)
    elseif ((pc > length(CDK2incPOM)) && (pc <= length(CDK2incPOM)+ length(CDK2lowPOM)) && (CDK2lowPOM(pc-length(CDK2incPOM))+ framesaftermitosis < framenum))
            ei = pc-length(CDK2incPOM);
            total_traces_after_gating = total_traces_after_gating+1;
        if isempty(find(CDK2lowPuncta(ei,CDK2lowPOM(ei)+6 : CDK2lowPOM(ei)+framesaftermitosis)))% +6 stands for 6th frame after POM.if isempty(find(CDK2lowPuncta(ei,(CDK2lowPOM(ei)+6)):(CDK2lowPOM(ei)+framesaftermitosis)))% +6 stands for 6th frame after POM.
                figure(1)
                plotdata=tracesCDK2low_24(ei,:);
                realframes=find(~isnan(tracesCDK2low_24(ei,:)));
                plotdata(realframes)=smooth(tracesCDK2low_24(ei,realframes),smoothingsize);
                tracecolor='r';
                cdk2low_wo_puncta = cdk2low_wo_puncta+1;
                traces_wo_puncta = traces_wo_puncta +1;
        else
                figure(2)
                plotdata=tracesCDK2low_24(ei,:);
                realframes=find(~isnan(tracesCDK2low_24(ei,:)));
                plotdata(realframes)=smooth(tracesCDK2low_24(ei,realframes),smoothingsize);
                tracecolor= 'r';
                cdk2low_with_puncta = cdk2low_with_puncta + 1;
                traces_with_puncta = traces_with_puncta +1;
        end
        line(cdk2lowshiftedframes(ei,:), plotdata, 'color', tracecolor)
   elseif ((pc > length(CDK2incPOM)+ length(CDK2lowPOM)) && (CDK2emergPOM(pc-length(CDK2incPOM)-length(CDK2lowPOM))+ framesaftermitosis < framenum))
            ai = pc-length(CDK2incPOM)-length(CDK2lowPOM);
            total_traces_after_gating = total_traces_after_gating+1;
        if isempty(find(CDK2emergPuncta(ai,CDK2emergPOM(ai)+6 : CDK2emergPOM(ai)+framesaftermitosis)))% +6 stands for 6th frame after POM.
                figure(1)
                plotdata=tracesCDK2emerg_24(ai,:);
                realframes=find(~isnan(tracesCDK2emerg_24(ai,:)));
                plotdata(realframes)=smooth(tracesCDK2emerg_24(ai,realframes),smoothingsize);
                tracecolor='g';
                cdk2emerg_wo_puncta = cdk2emerg_wo_puncta+1;
                traces_wo_puncta = traces_wo_puncta +1;
         else
                figure(2)
                plotdata=tracesCDK2emerg_24(ai,:);
                realframes=find(~isnan(tracesCDK2emerg_24(ai,:)));
                plotdata(realframes)=smooth(tracesCDK2emerg_24(ai,realframes),smoothingsize);
                tracecolor= 'g';
                cdk2emerg_with_puncta = cdk2emerg_with_puncta + 1;
                traces_with_puncta = traces_with_puncta +1;
         end
       line(cdk2emergshiftedframes(ai,:), plotdata, 'color', tracecolor)
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

for i=1:length(CDK2incPOM)
         
        plotdata=tracesCDK2inc_24(i,:);
        realframes=find(~isnan(tracesCDK2inc_24(i,:)));
        plotdata(realframes)=smooth(tracesCDK2inc_24(i,realframes),smoothingsize);
        figure(3)
        subplot(2,1,1)
        line(cdk2incshiftedframes(i,:),plotdata);
        title('CDK2 inc');
        xlim([-25 25])
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


for i=1:length(CDK2lowPOM)
         
        plotdata=tracesCDK2low_24(i,:);
        realframes=find(~isnan(tracesCDK2low_24(i,:)));
        plotdata(realframes)=smooth(tracesCDK2low_24(i,realframes),smoothingsize);
        figure(4)
        subplot(2,1,1)
        line(cdk2lowshiftedframes(i,:),plotdata);
        xlim([-25 25])
        title('CDK2 low');
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



for i=1:length(CDK2emergPOM)
         
        plotdata=tracesCDK2emerg_24(i,:);
        realframes=find(~isnan(tracesCDK2emerg_24(i,:)));
        plotdata(realframes)=smooth(tracesCDK2emerg_24(i,realframes),smoothingsize);
        figure(5)
        subplot(2,1,1)
        line(cdk2emergshiftedframes(i,:),plotdata);
        title('CDK2 emerg')
        xlim([-25 25])
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

alignedtime1 = alignedtime1./framesperhr;
alignedtime2 = alignedtime2./framesperhr;
alignedtime3 = alignedtime3./framesperhr;

plotdata2_mean = nanmean(incpuncta_datatotal);
plotdata3_mean = nanmean(lowpuncta_datatotal);
plotdata4_mean = nanmean(emergpuncta_datatotal);
plotCdk2emergdata_median = nanmedian(emergdatatotal);
plotCdk2emergdata_mean = nanmean(emergdatatotal);
plotCdk2incdata_median = nanmedian(incdatatotal);
plotCdk2incdata_mean = nanmean(incdatatotal);
plotCdk2lowdata_median = nanmedian(lowdatatotal);

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
% figure(3)
% subplot(2,1,2)
% errorbar(alignedtime1,plotdata2_mean,error1);
% ylim([-0.1 3])
% xlim([-25 25])

figure(3)
subplot(2,1,2)
errorbar(alignedtime1,plotdata2_mean,error1);
ylim([-0.1 3])
xlim([-25 25])

figure(4)
subplot(2,1,2)
errorbar(alignedtime2,plotdata3_mean,error2);
ylim([-0.1 3])
xlim([-25 25])

figure(5)
subplot(2,1,2)
errorbar(alignedtime3,plotdata4_mean,error3);
ylim([-0.1 3])
xlim([-25 25])

figure(6)
errorbar(alignedtime1,plotdata2_mean,error1,'b');hold on;
errorbar(alignedtime2,plotdata3_mean,error2,'r');hold on;
errorbar(alignedtime3,plotdata4_mean,error3,'g');hold off;
%         subplot(2,1,2)

figure(7)
plot(alignedtime1,plotdata2_mean,'b');hold on;
plot(alignedtime2,plotdata3_mean,'r');hold on;
plot(alignedtime3,plotdata4_mean,'g');hold off;
%%% record number of traces %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('%0.0f traces\n',numtraces);

figure(8)
% plot(alignedtime3,plotCdk2emergdata_mean,'k'); hold on;
% yyaxis left
% xlim([-15 15])
% ylim([0 1.6])
plot(alignedtime3,plotdata4_mean,'g');
% yyaxis right
xlim([-15 15])
ylim([0 2.0])

figure(9)
plot(alignedtime1,plotCdk2incdata_mean,'k'); hold on;
plot(alignedtime1,plotdata2_mean,'b');hold off;
xlim([-15 15])
ylim([-0.1 1.8])

figure(10)
plot(alignedtime2,plotCdk2lowdata_median,'k'); hold on;
plot(alignedtime2,plotdata3_mean,'r');hold off;
xlim([-14 14])
ylim([0 2])
