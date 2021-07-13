function Timelapse_puncta_backgroundcorr(row,col,site)
% row=5;col=2;site=1; 

%%% puncta related parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
puncta_sizethreshold=0.2;  %any "puncta" greater than 20% of the nuclear area will be rejected as a puncta
punctasize_lowerlimit=4;  %4 is the best puncta size lowerlimit for 10x 
nuclearradius=12;       %for MCF10A cells at 20X no binning, use nuclearradius=26; at 10X no binning use nuclearradius=13
strelradius=8;       %8 works best.this is the number of pixels around the puncta to take as background, which will get subtracted
getdapimaskradius=3;  %this is a guess at the radius of a typical puncta


%%% set paths %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
projectpath = '\\spencerstorage.int.colorado.edu\LabShare\IXMicroImages-goodNames\';
imagepath = '\\spencerstorage.int.colorado.edu\LabShare\IXMicroImages-goodNames\';
experimentpath = 'Mansi\MA72-MCF10A-mChyBP1-antioxidants_1861\';

biasdir =  'Y:\IXMicroImages-goodNames\Mansi\MA72-MCF10A-mChyBP1-antioxidants_1861\Illumination bias\';
cmosoffsetdir = '\\spencerstorage.int.colorado.edu\LabShare\IXMicroImages-goodNames\Justin\20160229 Illumination Bias Correction\';
load([cmosoffsetdir,'cmosoffset.mat'],'cmosoffset');

nucname = 'CFP-50'; 
sigg1 = 'YFP-50';
sigg2 = 'mCherry-50';
pos = 1;
load([biasdir,nucname,'_',num2str(pos),'.mat']); nucbias=bias; %load illumination bias calc
load([biasdir,sigg1,'_',num2str(pos),'.mat']); YFPbias=bias; %load illumination bias calc
load([biasdir,sigg2,'_',num2str(pos),'.mat']); mCherrybias=bias;

%shot=wellnum2str(row,col,site);
shot=[num2str(row),'_',num2str(col),'_',num2str(site)];
datadir=[projectpath,experimentpath,'PunctaData\'];

rawdir=[imagepath,experimentpath,'Raw\',shot,'_'];
% rawdir=[imagepath,experimentpath,shot,'_'];
maskdir=[imagepath,experimentpath,'Mask2'];

maskwrite=1;
if ~exist(maskdir,'dir') && maskwrite
    mkdir(maskdir);
end
maskdir=[maskdir,'\',shot,'_'];

%%% general settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%
SF=1;EF=309; %1:120
ringcalc=1; %set to zero if ring is unneeded
name1='CFP-50_'; %nuclear channel
name2='YFP-50_'; %this is always the CDK2sensor channel
puncta_channel = 'mCherry-50_';
namenucedge='nucedge_';
%%% segmentation parameters %%%%%%%%%%%%%%%%%%%%%
nucr=12;
debrisarea=100;
boulderarea=1500;
blobthreshold=-0.03;
%%% tracking parameters %%%%%%%%%%%%%%%%%%%%%%%%%
maxjump=nucr*4;
masschangethreshold=0.60;
areachangethreshold=0.60;
daughtervariance=0.10;
trackparams={nucr,maxjump,debrisarea,masschangethreshold,areachangethreshold,daughtervariance};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
frames=SF:EF;
totalframes=numel(frames);
badframes=ones(EF,1)*NaN;
if SF>1
    badframes(1:SF-1)=0;
end
jitters=zeros(EF,2);
blocksize=10000;
maxcellnum=blocksize;
parameternum=9;
tracedata=ones(maxcellnum,EF,parameternum)*NaN;
tracking=ones(maxcellnum,5)*NaN;
timetotal=tic;
[firstgoodindex,blurthreshhigh,blurthreshlow,numthresh,badframes,height,width]=timelapsesetup_4(rawdir,name1,frames,nucr,blobthreshold,debrisarea,badframes,maskwrite);
regheight=1:0.5*height; regwidth=1:0.5*width;
for i=firstgoodindex:totalframes
    f=frames(i); fprintf('frame %0.0f\n',f);
    timeframe=tic;
    %%% read images %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     raw1=double(imread([rawdir,name1,num2str(f),'.tif']));
%     raw2=double(imread([rawdir,name2,num2str(f),'.tif']));
%     cy3=single(imread([rawdir,puncta_channel,num2str(f),'.tif']));
%     punctaimage = cy3;
    
    temp1=double(imread([rawdir,name1,num2str(f),'.tif']));
    lesscmos1 = temp1 - cmosoffset;
    corrimg1 = lesscmos1./nucbias;
    raw1=corrimg1;
    
    temp2=double(imread([rawdir,name2,num2str(f),'.tif']));
    lesscmos2 = temp2 - cmosoffset;
    corrimg2 = lesscmos2./YFPbias;
    raw2=corrimg2;
    
    temp3=double(imread([rawdir,puncta_channel,num2str(f),'.tif']));
    lesscmos3 = temp3 - cmosoffset;
    corrimg3 = lesscmos3./mCherrybias;
    cy3=corrimg3;
    punctaimage = cy3;
    %%% segment nuclei %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if i==firstgoodindex
        firstsegmethod='log';
        switch firstsegmethod
            case 'log'
                nuc_mask=blobdetector_4(log(raw1),nucr,blobthreshold,debrisarea);
            case 'single'
                blurradius=3;
                nuc_mask=threshmask(raw1,blurradius);
                nuc_mask=markershed(nuc_mask,round(nucr*2/3));
            case 'double'
                blurradius=3;
                nuc_mask=threshmask(raw1,blurradius);
                nuc_mask=markershed(nuc_mask,round(nucr*2/3));
                nuc_mask=secondthresh(raw1,blurradius,nuc_mask,boulderarea*2);
        end
        foreground=nuc_mask;
        nuc_mask=bwareaopen(nuc_mask,debrisarea);
        %%% Deflection-Bridging %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        nuc_mask=segmentdeflections_bwboundaries(nuc_mask,nucr,debrisarea);
        eccentricitythresh=0.85;
        nuc_mask=excludelargeandwarped_3(nuc_mask,boulderarea,eccentricitythresh);
    else
        nuc_mask=threshmask(raw1,1);
        %%% calculate jitter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        lastgoodframe=find(badframes==0,1,'last');
        [reljitx,reljity]=registerimages(imfill(extractmask(regheight,regwidth),'holes'),nuc_mask(regheight,regwidth));
        jitters(f,:)=jitters(lastgoodframe,:)+[reljitx,reljity];
        secondsegmethod='apriori';
        switch secondsegmethod
            case 'log'
                nuc_mask=blobdetector_4(log(raw1),nucr,blobthreshold,debrisarea);
            case 'double'
                nuc_mask=threshmask(raw1,blurradius);
                nuc_mask=markershed(nuc_mask,round(nucr*2/3));
                nuc_mask=secondthresh(raw1,blurradius,nuc_mask,boulderarea*2);
            case 'apriori'
               [nuc_mask,marker_mask]=apriori_markermask(nuc_mask,nuc_center,jitters(f,:));
        end
        foreground=nuc_mask;
        nuc_mask=bwareaopen(nuc_mask,debrisarea);
    end
    %%% remove border objects %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nuc_mask=imclearborder(nuc_mask);
    %%% background subtract: masked %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    compression=4;
    nanmask=imdilate(foreground,strel('disk',nucr/2));
    nanmaskcyto=imdilate(foreground,strel('disk',nucr*2));
    blur1=imfilter(raw1,fspecial('disk',3),'symmetric');
    blur2=imfilter(raw2,fspecial('disk',3),'symmetric');
    real1=bgsubmasked_global_2(blur1,nanmask,11,compression,50);
    real2=bgsubmasked_global_2(blur2,nanmaskcyto,3,compression,10);
    %%% extract features %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [nuc_label,numcells]=bwlabel(nuc_mask);
    nuc_info_raw = regionprops(nuc_mask,real1,'Area','Centroid','MeanIntensity','PixelIdxList');
    nuc_info = struct2cell(nuc_info_raw');
    
%     [nuc_label,numcells]=bwlabel(nuc_mask);
%     nuc_info=struct2cell(regionprops(nuc_mask,real1,'Area','Centroid','MeanIntensity')');
     nuc_area=squeeze(cell2mat(nuc_info(1,1,:)));
    %%%%%% detect bad frame %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    mednuc=median(nuc_area);
    if i>firstgoodindex+1 && (numcells<numthresh || mednuc>blurthreshhigh || mednuc<blurthreshlow)   
        fprintf('badframe: frame %0.0f\n',f);
        badframes(f)=1;
        badextractmask=bwmorph(nuc_mask,'remove');
        if maskwrite
%             imwrite(uint16(badextractmask),[maskdir,'\',namenucedge,num2str(f),'.tif']);
%            MA edited the above line to remove backslash after maskdir on
%            1/6/2017
            imwrite(uint16(badextractmask),[maskdir,namenucedge,num2str(f),'.tif']); 
            %imwrite(uint16(real2),[maskdir,'\',name2,num2str(f),'.tif']);
            %imwrite(uint16(real3),[maskdir,'\',name3,num2str(f),'.tif']);
        end
        continue;
    end
    blurthreshhigh=1.1*mednuc;
    blurthreshlow=0.6*mednuc; % was 0.8
    numthresh=0.5*numcells;
    nuc_center=squeeze(cell2mat(nuc_info(2,1,:)))';
    nuc_density=squeeze(cell2mat(nuc_info(4,1,:)));
    %%% calculate masses %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nuc_mass=nuc_density.*nuc_area;
    
    %%% calculate puncta info %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    x = punctainfo_new(nuc_info_raw,uint16(punctaimage), puncta_sizethreshold, punctasize_lowerlimit, strelradius, getdapimaskradius); %x is a struct that is the updated nucxypos with puncta info
    
    %% Additions from IFfunction
    
    puncta1(length(nuc_info_raw),1)= struct();
    
    for j=1:length(nuc_info_raw)
        puncta1(j,1).puncta = 0;
    end

    
    for j=1:length(nuc_info_raw)  %for each cell in the image
       
         puncta1(j,1).puncta = x(j).puncta;   
        
    end
    number_of_puncta = struct2cell(puncta1);
    [nrows, puncta_data] = cellfun(@size, number_of_puncta);
    curdata=[nuc_center(:,1),nuc_center(:,2),nuc_area,nuc_mass,puncta_data'];
    %%% Addition from IF function ends here
    
    
    %curdata=[nuc_center(:,1),nuc_center(:,2),nuc_area,nuc_mass];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if i>firstgoodindex
        nuc_center(:,1)=nuc_center(:,1)+jitters(f,1);
        nuc_center(:,2)=nuc_center(:,2)+jitters(f,2);
        %%% temporarily store values %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        curdata=[nuc_center(:,1),nuc_center(:,2),nuc_area,nuc_mass,puncta_data'];%%% remove puncta_data from here
        debugpackage={extractmask,jitters(lastgoodframe,:),[reljitx,reljity]};
        %%% track & correct merges (update centers, masses and labels) %%%%
        [tracedata,curdata,tracking,nuc_label]=adaptivetrack_9_puncta(f,lastgoodframe,f,tracedata,curdata,tracking,real1,nuc_label,jitters(f,:),trackparams,debugpackage,punctaimage);
        badframes(f)=0;
    end
    %%% visualize extractions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    extractmask=bwmorph(nuc_label,'remove');
    if maskwrite
        imwrite(uint16(extractmask),[maskdir,namenucedge,num2str(f),'.tif']);
        %imwrite(uint16(real2),[maskdir,name2,num2str(f),'.tif']);
        %imwrite(uint16(real3),[maskdir,name3,num2str(f),'.tif']);
    end
    %%% extract features %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    cellid=find(~isnan(curdata(:,1)));
    numlivecells=numel(cellid);
    curdata=curdata(cellid,:);
    nuc_center=curdata(:,[1 2]);
    nuc_area=curdata(:,3);
    nuc_mass=curdata(:,4);
    num_puncta = curdata(:,5);
    nuc_info=regionprops(nuc_label,'PixelIdxList');
    nanvec=ones(numlivecells,1)*NaN; sig1=nanvec; sig2=nanvec;
    for n=1:numlivecells
        cc=cellid(n);
        sig1(n)=mean(real1(nuc_info(cc).PixelIdxList));
        sig2(n)=median(real2(nuc_info(cc).PixelIdxList));
    end
    if ringcalc==1
        ring_label=getcytoring_3(nuc_label,4,real2);
        ring_info=regionprops(ring_label,'PixelIdxList');
        sig2ring_75th=nanvec; sig2ring_fgmedian=nanvec;
        sig2thresh=100;
        for n=1:numlivecells
            cc=cellid(n);
            if cc>numel(ring_info)
                break;
            end
            ring2all=real2(ring_info(cc).PixelIdxList);
            ring2all(ring2all>prctile(ring2all,98))=[];
            sig2ring_75th(n)=prctile(ring2all,75);
            ring2foreground=ring2all(ring2all>sig2thresh);
            if numel(ring2foreground)<100
                 ring2foreground=ring2all;
            end
            if numel(ring2all)>100
                 sig2ring_fgmedian(n)=nanmedian(ring2foreground);
            end
        end
    end
    %%% store data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tracedata(cellid,f,:)=[nuc_center(:,1),nuc_center(:,2),nuc_area,nuc_mass,sig1,sig2,sig2ring_75th,sig2ring_fgmedian,num_puncta];
    if maxcellnum-max(cellid)<blocksize
        tempdata=ones(blocksize,EF,parameternum)*NaN;
        temptrack=ones(blocksize,5)*NaN;
        tracedata=[tracedata;tempdata];
        tracking=[tracking;temptrack];
        maxcellnum=maxcellnum+blocksize;
    end
    clear puncta1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    toc(timeframe);
end
[tracedata,genealogy,jitters]=postprocessing_nolinking(tracedata,cellid,jitters,badframes,tracking,maxcellnum,nucr);
%%% save data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save([datadir,'tracedata_withPuncta_',shot,'.mat'],'tracedata','genealogy','jitters');
toc(timetotal);
clear all; clear mex;

%{
%%% debugging: view images %%%%%%%%%%
extractmask=bwmorph(nuc_mask,'remove');
tempframe=zeros(height,width,3);
tempframe(:,:,1)=imadjust(mat2gray(raw1));
tempframe(:,:,2)=extractmask;;
%tempframe(:,:,3)=marker_mask;
figure,imshow(tempframe);

nuc_info=struct2cell(regionprops(nuc_mask,'Area')');
nuc_area=squeeze(cell2mat(nuc_info(1,1,:)));
hist(nuc_area,100);

anti_mask=bwareaopen(nuc_mask,debrisarea);
temp_mask=nuc_mask-anti_mask;
extractmask=bwmorph(temp_mask,'remove');

anti_mask=bwareaopen(nuc_mask,1000);
extractmask=bwmorph(anti_mask,'remove');
%}