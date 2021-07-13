clear all
close all

projectpath='\\spencerstorage.int.colorado.edu\LabShare\IXMicroImages-goodNames\Mansi\Microscopy data\Live-cell imaging\'; %path where the processed data will be stored; create a folder called "Data" here
experimentpath='MA61-20160419-mChyBP1-drugs_1556\'; 
datadir=([projectpath,experimentpath,'Data\']);
datadir2=([projectpath,experimentpath,'Data2\']);

rowmat=[8];
colmat=[1];
sitemat=[1];
for row=rowmat
    for col=colmat
        for site=sitemat
            %shot=wellnum2str(row,col,site);
            shot=[num2str(row),'_',num2str(col),'_',num2str(site)];
           
        end
    end
end



%load ([datadir, 'Row_' num2str(rowmat) '_Col_' num2str(colmat) 'puncta_keyinfo_v2.mat'])
load ([datadir, 'Row_' num2str(rowmat) '_Col_' num2str(colmat) '_keyinfo_Cdk2traces.mat'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%for CDK2low cells 
j=0;
i=0;
m = 0;
n = 0;
for k=1:length(Cdk2lowI)
    i=Cdk2lowI(k);
    figure(1), hold on
    plot(tracesCDK2(i,:))
    plot(frameofmitosis(i), tracesCDK2(i, frameofmitosis(i)), 'ko', 'markerfacecolor', 'k')
    plot([0 150], [.55 .55], 'r-')
    plot([(frameofmitosis(i)+25) (frameofmitosis(i)+25)], [0.1 1.9], 'k-')
    plot([(frameofmitosis(i)+40) (frameofmitosis(i)+40)], [0.1 1.9], 'k-')
    title(k)
    axis([0 150 .1 1.9])
    
    
    
    reply = input('Classify the data as cdk2 inc,emerging or low: ', 's');  %1: cdk2inc 2:cdk2emerging 3:cdk2low enter:skip
    
    if ~isempty(reply)  %if reply is not empty
        
        if(reply == '1')
            j = j+1;
            cdk2inctraces(j,:)=tracesCDK2(i,:);
            cdk2incnucarea(j,:)=nucarea(i,:);
            cdk2incframeofmitosis(j)=frameofmitosis(i);
            cdk2incframessincemitosis(j)=framessincemitosis(i);
            cdk2incindices(j)=i;
            cdk2incpuncta(j,:)=numberOfPuncta(i,:);
        
        elseif(reply == '2')
            m = m+1;
            cdk2emergtraces(m,:)=tracesCDK2(i,:);
            cdk2emergnucarea(m,:)=nucarea(i,:);
            cdk2emergframeofmitosis(m)=frameofmitosis(i);
            cdk2emergframessincemitosis(m)=framessincemitosis(i);
            cdk2emergindices(m)=i;
            cdk2emergpuncta(m,:)=numberOfPuncta(i,:);
        
        elseif(reply == '3')
            n = n+1;
            cdk2lowtraces(n,:)=tracesCDK2(i,:);
            cdk2lownucarea(n,:)=nucarea(i,:);
            cdk2lowframeofmitosis(n)=frameofmitosis(i);
            cdk2lowframessincemitosis(n)=framessincemitosis(i);
            cdk2lowindices(n)=i;
            cdk2lowpuncta(n,:)=numberOfPuncta(i,:);
        end
    end
    
    close(figure(1))
end
%save ([datadir, 'Row_' num2str(rowmat) '_Col_' num2str(colmat), '_goodcdk2lowcells'],'cdk2lowtraces','cdk2lownucarea','cdk2lowframeofmitosis', 'cdk2lowframessincemitosis', 'cdk2lowindices', 'cdk2lowIF')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%for CDK2inc cells now

%j=0;
i=0;
for k=1:length(Cdk2incI)
    i=Cdk2incI(k);
    figure(1), hold on
    plot(tracesCDK2(i,:))
    plot(frameofmitosis(i), tracesCDK2(i, frameofmitosis(i)), 'ko', 'markerfacecolor', 'k')
    plot([0 150], [.55 .55], 'c-')
    plot([(frameofmitosis(i)+25) (frameofmitosis(i)+25)], [0.1 1.9], 'k-')   %%MA: Plot a line at 5 hrs after frame of anaphase. This is the minimum amount of time the cell needs to be quiescent.
    plot([(frameofmitosis(i)+40) (frameofmitosis(i)+40)], [0.1 1.9], 'k-')  %%MA: Plot a line at 8 hrs after frame of anaphase. This is the max amount of time the cell can be quiescent. Just to get a tighter group of emerging cells for paper figure.
    title(k)
    axis([0 150 .1 1.9])
    
    
    
    reply = input('Classify the data as cdk2 inc,emerging or low: ', 's');  %1: cdk2inc 2:cdk2emerging 3:cdk2low  enter:skip
    
    if ~isempty(reply)  %if reply is not empty
        
        if(reply == '1')
            j=j+1;
            cdk2inctraces(j,:)=tracesCDK2(i,:);
            cdk2incnucarea(j,:)=nucarea(i,:);
            cdk2incframeofmitosis(j)=frameofmitosis(i);
            cdk2incframessincemitosis(j)=framessincemitosis(i);
            cdk2incindices(j)=i;
            cdk2incpuncta(j,:)=numberOfPuncta(i,:);
        elseif(reply == '2')
            m = m+1;
            cdk2emergtraces(m,:)=tracesCDK2(i,:);
            cdk2emergnucarea(m,:)=nucarea(i,:);
            cdk2emergframeofmitosis(m)=frameofmitosis(i);
            cdk2emergframessincemitosis(m)=framessincemitosis(i);
            cdk2emergindices(m)=i;
            cdk2emergpuncta(m,:)=numberOfPuncta(i,:);
        
        elseif(reply == '3')
            n = n+1;
            cdk2lowtraces(n,:)=tracesCDK2(i,:);
            cdk2lownucarea(n,:)=nucarea(i,:);
            cdk2lowframeofmitosis(n)=frameofmitosis(i);
            cdk2lowframessincemitosis(n)=framessincemitosis(i);
            cdk2lowindices(n)=i;
            cdk2lowpuncta(n,:)=numberOfPuncta(i,:);
        end
    end
    
    close(figure(1))
end
save ([datadir2, 'Row_' num2str(rowmat) '_Col_' num2str(colmat), '_allgoodcdk2cells'],'cdk2inctraces', 'cdk2incnucarea','cdk2incframeofmitosis', 'cdk2incframessincemitosis', 'cdk2incindices','cdk2lowtraces','cdk2lownucarea','cdk2lowframeofmitosis', 'cdk2lowframessincemitosis', 'cdk2lowindices','cdk2emergtraces', 'cdk2emergnucarea','cdk2emergframeofmitosis', 'cdk2emergframessincemitosis', 'cdk2emergindices','cdk2incpuncta','cdk2lowpuncta','cdk2emergpuncta')