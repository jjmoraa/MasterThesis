%main v4
%Code by Jose Mora @ UMass
%This one allows reading multiple section files
close all
clear
clc

parallel=0;
if parallel
    if isempty(gcp('nocreate'))
        parpool(5)
    end
end
%% initiation
%JJM: Very important to check these scripts and set your paths correctly!
addNumadPaths
addFrame3DDPath %jjm: This one is self made

%%%%%%%%%%%%%%%%%% Brief Description on input system %%%%%%%%%%%%%%%%%%%%
%JJM:
%The following lines are for reading the inputs for NuMADs blade object.
%The idea here is to prove the input file can come as a different thing
%than an excel, which is sort of NuMADs predilect (?) form. Recently,
%there has been a push to make input files in a common format (wind.io).
%We could definitely use that as input here, required some coding to get
%the information as a blade object
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get the current working directory (current folder)
currentFolder = pwd;
% Get the parent folder (folder above the current folder)
parentFolder = fileparts(currentFolder);

%read inputs
% Open the text file for reading
fid = fopen([parentFolder,'\Inputs\Driver.txt'], 'r');
tsr = fscanf(fid, '%f', 1);fgetl(fid); % Read the tip speed ratio (tsr)
fgetl(fid);
airfoils = textscan(fid, '%f %s', 'HeaderLines', 1);
fclose(fid);
col1=cell2mat(airfoils(:,1));
col2=cellfun(@string, airfoils(:,2), 'UniformOutput', false);
col2 = [col2{:}];
airfoils = table(col1,col2,'VariableNames',{'airfoil no','name'});

% Open the text file for reading
fid = fopen([parentFolder,'\Inputs\BladeSectionsv3.dat'], 'r');%CHECK!
%assign parameters and read comments/discard them
for i=1:4
    fgetl(fid);
end
inputs = textscan(fid, '%f %f %f %f %f %f %f', 'HeaderLines', 2);
fclose(fid); % Close the file
hubRad=4.118878;
wndspeed=10.65;

cases=[10 20 50 100 200];%beam elements can't be smaller than 1m


for count=2
    numel=cases(count);
    referenceCase=struct([]);
    inputs = array2table(cell2mat(inputs),'VariableNames',{'span (r) [m]','not used1','not used2','not used3','twist','chord','airfoil no'});
    
%% reference case run
    [a,ap,amids,apmids,cpvector,ctvector,cmvector,cpmids,ctmids,cmmids,ocp,oct,ocm,dT,dQ,totalpwr,refmoment,blade,D,R,F,L,Ks]=analysis_blocks_v2_00(inputs,hubRad,parentFolder,airfoils,numel);
    
    referenceCase(1).a=a;
    referenceCase(1).ap=ap;
    referenceCase(1).amids=amids;
    referenceCase(1).apmids=apmids;
    referenceCase(1).ocp=ocp;
    referenceCase(1).oct=oct;
    referenceCase(1).ocm=ocm;
    referenceCase(1).cpmids=cpmids;
    referenceCase(1).ctmids=ctmids;
    referenceCase(1).cmmids=cmmids;
    referenceCase(1).dT=dT;
    referenceCase(1).dQ=dQ;
    referenceCase(1).totalpwr=totalpwr;
    referenceCase(1).blade=blade;
    referenceCase(1).D=D;
    referenceCase(1).R=R;
    referenceCase(1).F=F;
    referenceCase(1).L=L;
    referenceCase(1).Ks=Ks;
    %topSecProps(count).secprops=blade.secprops.data;
    %ispans(count).ispan=blade.ispan;
    %Dstory(count).D=D;

    %reference case plots
    figure('Name', 'PreComp Analysis Results');
    data = [blade.secprops.data]; % Concatenate results to form 'data'
    for i = 2:size(data, 2)
        subplot(5, 5, i);
        plot(data(:, 1), data(:, i), 'b-o');
        xlabel(strrep(blade.secprops.labels{1}, '_', '\_'));
        ylabel(strrep(blade.secprops.labels{i}, '_', '\_'));
    end
    
%% jamiesons case

powervector=zeros(50,50,50);
powertable=table;
fprintf('a_in   n_in    p_in\n')
jamiesonsCase=struct([]);

%Future=parallel.FevalFuture;
index=0;
for a_in=33
    for n_in=11 %concave <1 convex
        for p_in=15
            index=index+1;
            fprintf('%4.4f %4.4f %4.4f \n',a_in,n_in,p_in)
            [jinputs]=jamieson(referenceCase,refmoment,wndspeed,a_in/100,n_in/10,p_in/10,hubRad,3,inputs,parentFolder,airfoils);
            
            %topSecProps=struct;
            %ispans=struct;
            %Dstory=struct;
                
                %Future(index)=parfeval(@analysis_blocks,13,jinputs,hubRad,parentFolder,airfoils,numel);
                [a,ap,amids,apmids,cpvector,ctvector,cmvector,cpmids,ctmids,cmmids,ocp,oct,ocm,dT,dQ,totalpwr,moment,blade,D,R,F,L,Ks]=analysis_blocks(jinputs,hubRad,parentFolder,airfoils,numel);
        
                %topSecProps(count).secprops=blade.secprops.data;
                %ispans(count).ispan=blade.ispan;
                %Dstory(count).D=D;
                
                jamiesonsCase(a_in,n_in,p_in).a=a;
                jamiesonsCase(a_in,n_in,p_in).ap=ap;
                jamiesonsCase(a_in,n_in,p_in).amids=amids;
                jamiesonsCase(a_in,n_in,p_in).apmids=apmids;
                jamiesonsCase(a_in,n_in,p_in).cpmids=cpmids;
                jamiesonsCase(a_in,n_in,p_in).ctmids=ctmids;
                jamiesonsCase(a_in,n_in,p_in).cmmids=cmmids;
                jamiesonsCase(a_in,n_in,p_in).ocp=ocp;
                jamiesonsCase(a_in,n_in,p_in).oct=oct;
                jamiesonsCase(a_in,n_in,p_in).ocm=ocm;
                jamiesonsCase(a_in,n_in,p_in).dT=dT;
                jamiesonsCase(a_in,n_in,p_in).dQ=dQ;
                jamiesonsCase(a_in,n_in,p_in).totalpwr=totalpwr;
                jamiesonsCase(a_in,n_in,p_in).blade=blade;
                jamiesonsCase(a_in,n_in,p_in).D=D;
                jamiesonsCase(a_in,n_in,p_in).R=R;
                jamiesonsCase(a_in,n_in,p_in).F=F;
                jamiesonsCase(a_in,n_in,p_in).L=L;
                jamiesonsCase(a_in,n_in,p_in).Ks=Ks;
                newRow={a_in,n_in,p_in,a,ap,ocp,oct,ocm,totalpwr,blade.span(end),D(3,end)};
                %D1=[jamiesonsCase(a_in,n_in,p_in).ocp,jamiesonsCase(a_in,n_in,p_in).oct,jamiesonsCase(a_in,n_in,p_in).ocm,jamiesonsCase(a_in,n_in,p_in).totalpwr,jamiesonsCase(a_in,n_in,p_in).blade.span(end),jamiesonsCase(a_in,n_in,p_in).D(3,end)];
                %D2=[referenceCase.ocp,referenceCase.oct,referenceCase.ocm,referenceCase.totalpwr,referenceCase.blade.span(end),referenceCase.D(3,end)];
                %P=[D1;D2];
                %spider_plot(P,'AxesLabels', {'ocp', 'oct', 'ocm', 'totalpwr', 'R' ,'tipdfl'})
                powertable=[powertable;newRow];
        end
    end
figure()
title('power')
heatmap(powertable,'Var2','Var3','ColorVariable','Var9')
figure()
title('radius')
heatmap(powertable,'Var2','Var3','ColorVariable','Var10')
figure()
title('tip dfl')
heatmap(powertable,'Var2','Var3','ColorVariable','Var11')
end
delete(gcp('nocreate'))
end