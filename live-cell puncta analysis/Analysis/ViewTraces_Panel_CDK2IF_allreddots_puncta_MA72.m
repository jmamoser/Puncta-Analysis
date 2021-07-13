function ViewTraces_Panel_CDK2IF_allreddots_MA72(row,col,site)
global rawdir maskdir tracedata jitters plotfignum immunoframe nucr channel channelnames IFchannelname

%row=2;col=4;site=1;
row=3;col=3;site=1;
shot=[num2str(row),'_', num2str(col), '_', num2str(site)];
projectpath = '\\spencerstorage.int.colorado.edu\LabShare\IXMicroImages-goodNames\';
imagepath = '\\spencerstorage.int.colorado.edu\LabShare\IXMicroImages-goodNames\';
experimentpath = 'Mansi\MA72-MCF10A-mChyBP1-antioxidants_1861\';
%experimentpath = 'Sara\20150211-Wes10A-H2BTurq-DHBVen-moviepart-take2_197\';
datadir=([projectpath,experimentpath,'PunctaData\']);
rawdir=[imagepath,experimentpath,'Raw\',shot,'_'];
maskdir=[imagepath,experimentpath,'Mask2\',shot,'_'];

immunoframe=121;

%%%%%% Data Settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nucr=12;
framesperhr=5;
drugspike=0/framesperhr;
frames=1:309; % for Mansi's movie
%frames=1:120;
channelnames={'mCherry_' 'YFP_' 'RFP_'};
IFchannelname='Cy5_';
edgemask=1; %0:no masks saved 1: masks saved
%%%%%% User Settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ensemble=0;    %0:one trace per plotn 1:all traces in one plot
selectmode=0;  %View trace images
selectin=[];   %Choose trace IDs (necessary for selectmode)
%selectin = cellid_puncta(:,1);
plotsignal2=0;
plotsignal3=0;
IFoption=0; %0:no IFdata 1:IFdata
tracemax=400; %max number of traces to plot
%%%%%% Axis Settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
y2min=0; y2max=2.2;
y3min=0; y3max=4;
%%% Cell-Cycle Analysis Settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
motheroption=0; %0:no gating 1:mothers 2:no mothers
daughteroption=1; %0:no gating 1:daughters 2:no daughters
quiescentanalysis=0;
%%%%%% Load Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[tracedata,jitters,tracestats,tracestats_og,motherstats,IFdata,IDs]=gathertracedata_4_allreddots_puncta(datadir,shot,motheroption,daughteroption,IFoption);
%[tracedata,jitters,tracestats,tracestats_og,motherstats,IFdata,IDs]=gathertracedata_4_allreddots(datadir,shot,motheroption,daughteroption,IFoption);
numcells=size(tracedata,1);
%%Bug Fix to plot traces that undergo mitosis at the end of the movie
numcells=size(tracedata,1);
for c=1:numcells
    lengthoftrace(c,1)=find(~isnan(tracedata(c,:,1)),1,'first');
    lengthoftrace(c,2)=find(~isnan(tracedata(c,:,1)),1,'last');
end
lengthoftrace(:,3)=lengthoftrace(:,2)-lengthoftrace(:,1)+1;
%%% gate length %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
minlengthtrace=100;
minlengthmother=50;
badlengths=lengthoftrace(:,3)<minlengthtrace; %| motherstats(:,3)<minlengthmother;
%%% set and gate signal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
signalchoice='CDK2';
switch signalchoice
    case 'CDK2'
        y1min=0.3; y1max=2.3; % was 2.8
        nucchannel=6; cytochannel=8; % change this to the tracedata column according to your respective movie
        nucthreshoption=1;
        %0:normalize each trace by percentile of entire trace
        %1:normalize by first frame of trace (regardless of whether mitosis or not)
        nucthresh=300;    %threshold for nuc intensity according to nucthreshoption
        motherthresh=0;   %threshold for max DHB ratio of mother. default gating=1.0.  0:ignore
        noisethresh=0.2;  %threshold for max positive jump in cyto/nuc ratio
        [traces1,badtraces1]=gate_Cdk2_8_mother(tracedata,nucchannel,cytochannel,tracestats,nucthreshoption,nucthresh,motherthresh,noisethresh,quiescentanalysis);
    case 'Geminin'
        y1min=0; y1max=5000;
        channelGem=7;
        maxpos=0;      %0:anywhere 1:firstframe 2:lastframe
        maxthresh=500;  %threshold above which max of each trace must be %50
        minpos=2;      %0:anywhere 1:mothertrace 2:daughtertrace
        minthresh=500; %threshold below which min of each trace must be %50
        [traces1,badtraces1]=gate_Geminin_9_mother(tracedata,tracestats,motherstats,channelGem,maxthresh,minthresh,maxpos,minpos);
    case 'Area'
        y1min=100; y1max=1000;
        traces1=tracedata(:,:,3);
        badtraces1=zeros(numcells,1);
end
%%% Gate IF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
IF=zeros(numcells,1);
IFmin=8; IFmax=14;
if IFoption
    IF=IFdata(:,9); IF(IF<1)=1; IF=log2(IF);
    traces1(:,end)=traces1(:,end-1); %necessary for selectmode
end

%%% combine all gates and remove %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%gating=zeros(length(IDs),1);
gating=badlengths | badtraces1;

traces1=traces1(~gating,:);
IF=IF(~gating,:);
%traces2=traces2(~badtraces,:);
%traces3=traces3(~badtraces,:);
tracestats=tracestats(~gating,:);
pointsOfMitosis = tracestats;
pointsOfMitosis(:,2:4) = [];
%stuff to save:
tracedata=tracedata(~gating,:,:);
gatedTracedata = tracedata; % gated tracedata
tracesCDK2 = traces1; % CDK2 traces
IDs=IDs(~gating,:);
punctaData = tracedata(:,:,9);

%%% Removing puncta data more than 10 and replacing it with previous puncta value

for i = 1:size(punctaData,1)
    spikes = find(punctaData(i,:) > 10);  %%removes puncta values of more than 10
    if ~isempty(spikes)
      for j = 1:length(spikes)
        punctaData(i,spikes(j)) = punctaData(i,spikes(j)-1);
      end
    end
end

% uncomment this if you want to save all the tracedata, which includes centroids



numgated=size(tracestats,1);
if numgated>tracemax
    numgated=tracemax;
end
%%% plotting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set(0,'Units','pixels'); screendims=get(0,'ScreenSize'); screenx=screendims(3); screeny=screendims(4);
%%%%%% Viewer Settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dotsize=8; %default=8 presentation=12
tracewidth=2; %default=2 presentation=4
steps=5;
y1step=round((y1max-y1min)/steps);
y2step=round((y2max-y2min)/steps);
IFstep=round((IFmax-IFmin)/steps);
trace2color=[1 0 0];
if selectmode==0
    drugtime=drugspike;
    drugplot=0;
else
    drugtime=0;
    drugplot=drugspike;
end
if IFoption
    frames=[frames frames(end)+1];
end
xtime=frames/framesperhr;
xtime=xtime-drugtime;
selection=1:numgated;
if ~isempty(selectin)
    selection=find(ismember(IDs,selectin));
end
if ensemble
    figure; hold on;
end

for counter=1:length(selection)
    i=selection(counter);
    if selectmode
        clf;
        set(gcf,'Position',[round(0.6*screenx) round(0.05*screeny) round(0.38*screenx) round(0.38*screeny)]);
        plotfignum=get(gcf,'Number');
    elseif ~ensemble
        figure(ceil(counter/24));
        subaxis(4,6,mod(counter-1,24)+1,'ML',0.02,'MR',0.01,'MT',0.03,'MB',0.03,'SH',0.02); %5x4
    end
    set(gcf,'color','w');
    ysig1=smoothignorenans(traces1(i,:),1);
    if plotsignal2
        ysig2=smoothignorenans(traces2(i,:),3);
        [haxes,hline1,hline2]=plotyy(xtime,ysig1,xtime,ysig2);
        axes(haxes(1));
        axis([xtime(1) xtime(end) y1min y1max]);
        set(hline1,'DisplayName',num2str(i),'color','b','linewidth',tracewidth);
    elseif IFoption
        [haxes,hline1,hline2]=plotyy(xtime,ysig1,xtime(end),IF(i));
        axes(haxes(1));
        axis([xtime(1) xtime(end) y1min y1max]);
        set(hline1,'DisplayName',num2str(i),'color','b','linewidth',tracewidth);
    else
        yyaxis left
        line(xtime,ysig1,'DisplayName',num2str(i),'linewidth',tracewidth);
        axis([xtime(1) xtime(end) y1min y1max]);
        yyaxis right
        line(xtime,punctaData(i,:),'DisplayName',num2str(i),'linewidth',tracewidth);
        axis([xtime(1) xtime(end) y3min y3max]);
    end
    hold on;
    %%%%% mark drug spike and title %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if drugspike
        line([drugplot drugplot],[y1min y1max],'Color','k','linewidth',tracewidth,'linestyle','--');
    end
    %%%%% mark mitosis (only the current cell's birth) %%%%%%%%%%%%%%%%%%%%

    
    cellid_w_mother = i; %  i is the cell index from tracestats
    if ~isnan(pointsOfMitosis(cellid_w_mother,1))
        relmitosis = pointsOfMitosis(cellid_w_mother,1);
        yyaxis left
        plot(xtime(relmitosis),ysig1(relmitosis),'ro','markerfacecolor', 'r','markersize',dotsize);
        hold on;
        c = 1;
        
        while((c <= (size(pointsOfMitosis,2))) && (pointsOfMitosis(cellid_w_mother,c) ~= 0))
            if((c) < size(pointsOfMitosis,2) && (pointsOfMitosis(cellid_w_mother,c+1) ~= 0) )
                relmitosis = pointsOfMitosis(i,c);
                yyaxis left
                plot(xtime(relmitosis),ysig1(relmitosis),'ro','markerfacecolor', 'r','markersize',dotsize);
                hold on;
            end
            c = c+1;
        
        end
    end
%     if ~isnan(tracestats(cellid_w_mother,4))
%         relmitosis = tracestats(cellid_w_mother,1);
%         plot(xtime(relmitosis),ysig1(relmitosis),'ro','markerfacecolor', 'r','markersize',dotsize);
%         hold on;
%         c = 1;
%     
%         while((c <= (size(tracestats,2)-4)) && (tracestats(cellid_w_mother,4+c) ~= 0))
%             if((4+c) < size(tracestats,2) && (tracestats(cellid_w_mother,5+c) ~= 0) )
%                 relmitosis = tracestats(i,4+c);
%                 plot(xtime(relmitosis),ysig1(relmitosis),'ro','markerfacecolor', 'r','markersize',dotsize);
%                 hold on;
%             end
%             c = c+1;
%         
%         end
%     end
%     
    %%% Adjust signal2 settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if plotsignal2
        axes(haxes(2));
        axis([xtime(1) xtime(end) y2min y2max]);
        set(gca,'YAxisLocation','right','YColor',trace2color,'YTick',y2min:ystep2:y2max);
        set(hline2,'DisplayName',num2str(i),'Color',trace2color,'linewidth',tracewidth);
    end
    if IFoption
        axes(haxes(2));
        hold on;
        scatter(xtime(end),IF(i),50,'MarkerFaceColor','w','MarkerEdgeColor',trace2color,'linewidth',tracewidth);
        axis([xtime(1) xtime(end) IFmin IFmax]);
        set(gca,'YAxisLocation','right','YColor',trace2color,'YTick',IFmin:IFstep:IFmax);
    end
    if plotsignal3
        ysig3=smoothignorenans(traces3(i,:),3);
        line(xtime,ysig3,'color','g','DisplayName',num2str(i),'linewidth',tracewidth);
    end
    title(num2str(IDs(i)));
    %%% operate image-viewing mode %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if selectmode
        selectmodedisplay_1(edgemask);
    end
end
fprintf([num2str(numgated),'\n']);
save ([datadir, 'Row_' num2str(row) '_Col_' num2str(col), '_gatedTracedata_redDots'],'gatedTracedata', 'tracesCDK2','IDs','pointsOfMitosis');
%{
set(gcf,'color','w','PaperPosition',[0 0 20 12]);
saveas(gcf,'D:\Downloads\Fig.jpg');
%}