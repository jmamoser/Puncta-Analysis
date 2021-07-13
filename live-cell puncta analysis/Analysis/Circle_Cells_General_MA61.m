%% This function is to track the cell with the given cell id in the last frame
%% Not Accurate for cells in mitosis because 53BP1 is not
%% present in mitosis. The script requires a cellid list to be entered. It
%% then circles the cells in the cellid list and displays their cellid.

%%Parameters to be changed:
% Change the row,col,site,projectpath,imagepath and experimentpath
% variable named lastframe at approx. line #25 needs to be changed to the
% last frame of the movie.
% A cellid_list needs to be provided, which is a list of all the cells you
% want the script to circle.
function Circle_Cells_General(row,col,site,cellid_list)

row=1;col=1;site=1;

%%Setup paths for the image data
projectpath = '\\spencerstorage.int.colorado.edu\LabShare\IXMicroImages-goodNames\Mansi\';
imagepath = '\\spencerstorage.int.colorado.edu\LabShare\IXMicroImages-goodNames\Mansi\';
experimentpath = 'Microscopy data\Live-cell imaging\MA61-20160419-mChyBP1-drugs_1556\';

mkdir([projectpath,experimentpath,'analyzedData_test'])
savingdir=[projectpath,experimentpath,'analyzedData_test\'];
shot=[num2str(row),'_',num2str(col),'_',num2str(site)];
datadir=[projectpath,experimentpath,'Data\'];
%datadir=[projectpath,experimentpath,'Data\'];
rawdir=[imagepath,experimentpath,'Raw\', shot,'_'];
IFoption=0;
lastframe = 100;%Enter the last frame
if(IFoption)
    name1 = 'CFP_stain.tif';
    name2 = 'YFP_stain.tif';
else
    name1 = ['mCherry-50_',num2str(lastframe),'.tif'];
    name2 = ['YFP-50_',num2str(lastframe),'.tif'];
end
figure

I1 = (imread([rawdir,name1])); 
I2 = (imread([rawdir,name2]));
% subplot(1,2,1), imshow(I1,'DisplayRange',[0 1500])%This is the CDK2 image
% subplot(1,2,2), imshow(I2,'DisplayRange',[2000 20000])%This is the 53BP1 image
%load data
load([datadir,'tracedata_withPuncta_LL4_SR8_DM3_',shot,'.mat'],'tracedata','jitters');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%List of the cellids
cellid_list = [872];
%cellid_list = [1378 1557 1847 2035 2162 2295 2322 2474 2506 2835 3136 3291 3343 3558 3580 3583 3603 3604];% All the jagged ones from 1_4_1
%cellid_list = [2402 2648 2938 1303 2250];% All the jagged ones from 2_4_1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:length(cellid_list)
    %prompt1 = 'Please enter the cellid: ';
    %cellid = input(prompt1)
    cellid  = cellid_list(i);
    nuc_centroid = tracedata(cellid,lastframe,[1,2]);
    x1 = nuc_centroid(:,1)-jitters(lastframe,1);
    y1 = nuc_centroid(:,2)-jitters(lastframe,2);
    str = num2str(cellid);
    
    
    pos = [x1+15 y1+15];
   
    value = [cellid];
   
%     figure(1)
    I1 = insertText(I1,pos,value,'textColor','w','BoxColor','black','FontSize',10);
    I1 = insertShape(I1,'circle',[x1 y1 20],'Color','w');
    Gray1 = rgb2gray(I1); % For CDK2 images
     
    I2 = insertText(I2,pos,value,'textColor','w','BoxColor','black','FontSize',10);
    I2 = insertShape(I2,'circle',[x1 y1 30],'Color','w');
    Gray2 = rgb2gray(I2); % For 53BP1 images
    
    subplot(1,2,1), imshow(Gray1)%,'DisplayRange',[0 15000])
    subplot(1,2,2), imshow(Gray2,'DisplayRange',[0 5000])
    hold on;
end
% load([datadir,'tracedata_',shot,'.mat'],'tracedata');
% YFP_trace = tracedata(cellid,:,8)./tracedata(cellid,:,6);
% frames = 1: 120;
% framesperhr = 5;
% x = frames/framesperhr;
% figure
% plot(x,YFP_trace)
% ViewTraces_Panel_CDK2(cellid_list,row,col,site);

