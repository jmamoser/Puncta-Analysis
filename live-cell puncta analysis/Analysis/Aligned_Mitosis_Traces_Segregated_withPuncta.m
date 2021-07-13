close all;
clear all;
%%% filepaths %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
projectpath='\\spencerstorage.int.colorado.edu\LabShare\IXMicroImages-goodNames\Mansi\Microscopy data\Live-cell imaging\'; %path where the processed data will be stored; create a folder called "Data" here
imagepath=  '\\spencerstorage.int.colorado.edu\LabShare\IXMicroImages-goodNames\Mansi\Microscopy data\Live-cell imaging\';  %where the images are stored, generally
% experimentpath='MA55-20160206-mCitp21-p21null-drugs_774\'; %where the raw data for THIS expt is stored; put inames into a new subdirectory here called Raw
experimentpath='MA61-20160419-mChyBP1-drugs_1556\'; 
% experimentpath='MA56-20160211-VitC-NAC-pilot3colorMCF10A_782\';
datadir=([projectpath,experimentpath,'Data\']);

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
rowmat=[8];
colmat=[1];
sitemat=[1];
tracedata=[];
tracestats=[];
motherstats=[];
IFdata=[];
for row=rowmat
    for col=colmat
        for site=sitemat
            %shot=wellnum2str(row,col,site);
            shot=[num2str(row),'_',num2str(col),'_',num2str(site)];
            [tracedatatemp,tracestatstemp,motherstatstemp,IFdatatemp]=gathertracedata_1_puncta_overtime(datadir,shot,motheroption,daughteroption,IFoption);
            tracedata=[tracedata;tracedatatemp];
            tracestats=[tracestats;tracestatstemp];
            motherstats=[motherstats;motherstatstemp];
            IFdata=[IFdata;IFdatatemp];
        end
    end
end
tracedata=tracedata(tracestats(:,1)<=framenum,:,:); %%MM addition
tracestats=tracestats(tracestats(:,1)<=framenum,:); %%MM addition
tracestats(tracestats(:,2)>framenum,2)=framenum; %%MM addition
motherstats=motherstats(tracestats(:,1)<=framenum,:); %%MM addition

numcells = size(tracedata , 1);
 for c=1:numcells
    lengthoftrace(c,1)=find(~isnan(tracedata(c,:,1)),1,'first');
    lengthoftrace(c,2)=find(~isnan(tracedata(c,:,1)),1,'last');
end
lengthoftrace(:,3)=lengthoftrace(:,2)-lengthoftrace(:,1)+1;

%%% gate length %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
minlengthtrace=50;
minlengthmother=30;
badlengths=lengthoftrace(:,3)<minlengthtrace; % | motherstats(:,3)<minlengthmother;




%%% gate CDK2 data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nucchannel=6; cytochannel=8;
nucthreshoption=0;
%0:normalize each trace by percentile of entire trace
%1:normalize by first frame of trace (regardless of whether mitosis or not)
nucthresh=300;    %threshold for nuc intensity according to nucthreshoption
motherthresh=1.2;   %threshold for max DHB ratio of mother. default gating=1.0.  0:ignore
noisethresh=0.2;  %threshold for max positive jump in cyto/nuc ratio
[tracesCDK2,badtracesCDK2]=gate_Cdk2_8_mother(tracedata,nucchannel,cytochannel,tracestats,nucthreshoption,nucthresh,motherthresh,noisethresh,quiescentanalysis);
%%% gate Geminin data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% channelGem=6;
% maxpos=1;      %0:anywhere 1:mothertrace 2:daughterlastframe
% maxthresh=50;  %threshold above which max of each trace must be %50
% minpos=2;      %0:anywhere 1:mothertrace 2:daughtertrace
% minthresh=20; %threshold below which min of each trace must be %50
% [tracesGeminin,badGeminin]=gate_Geminin_9_mother(tracedata,tracestats,motherstats,channelGem,maxthresh,minthresh,maxpos,minpos);

%%% combine gateout %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%badtraces=zeros(size(tracedata,1),1);
badtraces=badlengths | badtracesCDK2;
%%% remove gateout %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tracesCDK2=tracesCDK2(~badtraces,:);
tracedata = tracedata(~badtraces,:,:);

%%% to change all the NaNs in the puncta column to 0
ind = find(isnan(tracedata(:,:,9)));
tracedataB = tracedata(:,:,9);
tracedataB(ind) = 0;

%%% Removing data more than 10 and replacing it with previous puncta value

for i = 1:size(tracedataB,1)
    spikes = find(tracedataB(i,:) > 10);
    if ~isempty(spikes)
      for j = 1:length(spikes)
        tracedataB(i,spikes(j)) = tracedataB(i,spikes(j)-1);
      end
    end
end


tracedata(:,:,9) = tracedataB;

%traces2=traces2(~badtraces,:);
tracestats=tracestats(~badtraces,:);
motherstats=motherstats(~badtraces,:);
%%% normalize traces by max in lineage %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
normoption=0;
%0:normalize each trace by max of its own trace
%1:normalize by median of max of all traces if min>max*0.5
%2:normalize each trace by value at first frame after mitosis
%tracesGeminin=normalizetraces_3(tracesGeminin,tracestats,normoption);

%%% Obtain POI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
POIoption='Mitosis';
switch POIoption
    case 'Mitosis'
        POI=tracestats(:,1);
    case 'R' %requires minlength>=15
        %%% get Points Of Interest %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [POI_Rpoint,badtraces_POICDK2]=GetPointsOfInterest_2('Cdk2_R',tracesCDK2,tracestats,minlengthtrace);
        %[POI_Geminin,badtraces_Geminin]=GetPointsOfInterest_2('Geminin',traces1,tracestats,minlengthtrace);
        %%% gate out badPOI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %badtracesallPOI=zeros(size(POIcalc,1),1); %no bad trace gating
        badtracesallPOI=badtraces_POICDK2;
        POI_Rpoint=POI_Rpoint(~badtracesallPOI);
        tracesCDK2=tracesCDK2(~badtracesallPOI,:);
        tracestats=tracestats(~badtracesallPOI,:);
        motherstats=motherstats(~badtracesallPOI,:);
        %%% define alignment points %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        POI=POI_Rpoint;
end



% %%% for removing traces with mislocalised sensor %%%%%%%
% numtraces=size(tracestats,1);
% temp = 1;
% for cc = 1:numtraces
%    if ~isempty(find(tracesCDK2(cc,(POI(cc)+5)) > 0.9))
%        gated_traces(temp) = cc;
%        temp = temp + 1;
%     end
% end       
% 
% tracesCDK2(gated_traces,:) = [];
% tracedata(gated_traces,:,:) = [];
% tracestats(gated_traces,:) = [];
% motherstats(gated_traces,:) = [];
% POI(gated_traces)= [];

%%% categorize traces %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numtraces=size(tracestats,1);
calldelay=20; %must be >= minlengthtrace
callthresh=0.55;
categorization=ones(numtraces,1)*NaN;
valstore=ones(numtraces,1)*NaN;

for cc=1:numtraces
   calltime=tracestats(cc,1)+calldelay-1; % #frames after mitosis
   if(calltime < size(tracesCDK2,2)-calldelay)
        categorization(cc)=tracesCDK2(cc,calltime)>callthresh;
        valstore(cc)=min(tracesCDK2(cc,tracestats(cc,2)-2:tracestats(cc,2)));
    end
end


Cdk2inc=find(categorization==1);
Cdk2low=find(categorization==0);
%%% align traces %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[alignedtime,aligneddata]=aligntraces_4(tracesCDK2,POI,tracestats,motherstats,daughteroption);
alignedtime=alignedtime/framesperhr;

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

%%% plot traces %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
smoothingsize=5;
ymin=0.2; ymax=2.5;
figure(1), hold on;
figure(2), hold on;
traces_wo_puncta = 0;
traces_with_puncta = 0;
cdk2low_wo_puncta = 0;
cdk2low_with_puncta = 0;
cdk2inc_with_puncta = 0;
cdk2inc_wo_puncta = 0;
total_traces_after_gating = 0;
for pc=1:numtraces
    if (POI(pc)+framesaftermitosis < framenum) 
        total_traces_after_gating = total_traces_after_gating+1;
        if isempty(find(tracedata(pc,POI(pc)+6 :(POI(pc)+framesaftermitosis),9)))% +6 stands for 6th frame after POM.
            figure(1)
            plotdata=aligneddata(pc,:);
            realframes=find(~isnan(aligneddata(pc,:)));
            plotdata(realframes)=smooth(aligneddata(pc,realframes),smoothingsize);
            if ismember(pc,Cdk2inc)
                tracecolor='b';
                cdk2inc_wo_puncta = cdk2inc_wo_puncta+1;
            elseif ismember(pc,Cdk2low)
                tracecolor='r';
                cdk2low_wo_puncta = cdk2low_wo_puncta +1;
            end
            traces_wo_puncta = traces_wo_puncta +1;
        else
            figure(2)
            plotdata=aligneddata(pc,:);
            realframes=find(~isnan(aligneddata(pc,:)));
            plotdata(realframes)=smooth(aligneddata(pc,realframes),smoothingsize);
            if ismember(pc,Cdk2inc)
                tracecolor='b';
                cdk2inc_with_puncta = cdk2inc_with_puncta + 1;
            elseif ismember(pc,Cdk2low)
                tracecolor='r';
                cdk2low_with_puncta = cdk2low_with_puncta +1;
            end
            traces_with_puncta = traces_with_puncta +1;
        end
        line(alignedtime,plotdata,'color',tracecolor);
    end
end

nucarea=tracedata(:, :, 3);  %how to get the nucarea traces that correspond to the selected cdk2 traces?
Cdk2incI=Cdk2inc; %indices of CDK2inc cells
Cdk2lowI=Cdk2low; %indices of CDK2low cells
frameofmitosis=POI;
framessincemitosis=tracestats(:,3);
numberOfPuncta = tracedata(:,:,9);


length(Cdk2inc);
length(Cdk2low);
pct_Cdk2low = ((cdk2low_wo_puncta+cdk2low_with_puncta)/total_traces_after_gating).*100;
pct_cdk2low_wo_puncta = (cdk2low_wo_puncta/traces_wo_puncta)*100;
pct_cdk2low_with_puncta = (cdk2low_with_puncta/traces_with_puncta)*100;
pct_cdk2low_with_puncta_wrt_cdk2low = (cdk2low_with_puncta/(cdk2low_wo_puncta+cdk2low_with_puncta))*100;
pct_cdk2inc_with_puncta_wrt_cdk2inc = (cdk2inc_with_puncta/(cdk2inc_with_puncta+cdk2inc_wo_puncta))*100;
%temp1 = (cdk2low_wo_puncta/length(Cdk2low))*100;
%xlim([-10 15]);
%ylim([ymin ymax]);
xlabel('Time (hrs)');
ylabel('CDK2 Activity');
set(gcf,'color','w','PaperPosition',[0 0 8 6]);
figure(1)
title(['Row: ' num2str(rowmat) '  Col: ' num2str(colmat),'  ,Total Traces: ',num2str(numtraces),'  ,%CDK2low of Total Traces: ' num2str(pct_Cdk2low),' ,Traces w/o puncta: ', num2str(traces_wo_puncta),' ,%CDK2low of Traces w/o puncta: ',num2str(pct_cdk2low_wo_puncta)],'fontsize', 14 );

%title(  ['Row: ' num2str(rowmat) '  Col: ' num2str(colmat) '       Total Traces: ' num2str(numtraces)  '    %CDK2low of Total Traces:  ' num2str(pct_Cdk2low)] , 'fontsize', 16)
saveas(gcf, [datadir, 'Row_' num2str(rowmat) '_Col_' num2str(colmat) '_CDK2alignedToMitosisredblue.fig']);
figure(2)
title(['Row: ' num2str(rowmat) '  Col: ' num2str(colmat),'  ,Total Traces: ',num2str(numtraces),'  ,%CDK2low of Total Traces: ' num2str(pct_Cdk2low),' ,Traces with puncta: ', num2str(traces_with_puncta),' ,%CDK2low of Traces with puncta: ',num2str(pct_cdk2low_with_puncta)],'fontsize', 14 );

%%% save data
save ([datadir, 'Row_' num2str(rowmat) '_Col_' num2str(colmat), '_keyinfo_Cdk2traces'], 'frameofmitosis', 'framessincemitosis', 'tracesCDK2', 'Cdk2incI', 'Cdk2lowI', 'nucarea','numberOfPuncta')

%%% record number of traces %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('%0.0f traces\n',numtraces);