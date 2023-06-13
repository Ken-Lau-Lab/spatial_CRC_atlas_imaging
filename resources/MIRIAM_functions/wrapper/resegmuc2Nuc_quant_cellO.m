function resegmuc2Nuc_quant_cell(SlideDir, NucList, MemList, EpiList ,quantify, unsharp)
%Eliot McKinley - 2/23/2017
% No returned output - function saves image of objects and quantification
% of AFremoved markers
%
%SlideDir=directory for slide containing AFRemoved images folder and
%Registered images folder - assumes Round 001 is baseline and all fies are
%.tif - if no folder is supplied- then ui asks for directory
%
%NucList= cell array of nuclear markers to include besides DAPI these must
%match the marker name in AFRemoved images, e.g. {'SOX9' }; if empty cell is
%passed user will be asked to slect markers
%
%MemList= cell array of membrane markers to combine these must match the
%marker name in AFRemoved images,e.g. {'BCAT' 'ECAD' 'NAKATPASE' 'VILLIN'}; if empty cell is
%passed user will be asked to select markers
%
%EpiListList= cell array of epithelila markers to combine these must match the
%marker name in AFRemoved images,e.g. {'BCAT' 'ECAD' 'PCK26' 'MUC2'}; if empty cell is
%passed user will be asked to select markers

%quantify= whether or not to quantify 1= yes, 0=no

%unsharp=option to use unsharp mask for adaptive filters 1= yes 0=no

%% Custom functions needed
%SegDirFormatting-gets files into correct naming format
%main- generates initial nuclear segmentation (mask), supermembrane, and a binary membrane mask file
%tubsegbatch - generates epithelial mask
%Muc2Add- epithelial mask with Muc2 added and cells seperated by membranes
%shortBatch- generates cell mask in cells with high muc2
%ReSegCells- generates resegmented cell mask
%NucCountBatch - resegments cells with multiple nuclei- generating new cellmask and nuclear mask


%% Parse Directory supplied for cell segmentation

%if no directory is supplied ask for user input
if isempty(SlideDir)
    SlideDir=uigetdir;
end

% get formatting for AFRemoved, DAPI, and output directories
[AFRemoved, DAPI, OutDir]=SegDirFormatting(SlideDir);

%get AFremoved images
AFList=dir([AFRemoved '*.tif']);
AFList={AFList.name};
AFList=regexp(AFList, '_AFRemoved_', 'split'); %split characters before and after delimiter
AFList=cat(3,AFList{:}); %reorganize cell
PosList=squeeze(AFList(1,2,:)); %get positions with .tif
AFList=unique(squeeze(AFList(1,1,:))); %get only the Marker names
PosList=unique(regexprep(PosList,'.tif',''));  %get only the positions
OutPos=PosList; %format for output position

%get DAPI images from first round of imaging (assumes Olympus)
DapiList=dir([DAPI '*_dapi.tif']);
DapiList={DapiList.name};
DapiList=strcat(DAPI ,DapiList); %add in path

% Format DAPI images for Cytell based imaging
if isempty(DapiList)
    DapiList=dir([DAPI '*_dapi_*.tif' ]);
    DapiList={DapiList.name};
    DapiList=strcat(DAPI ,DapiList); %add in path
    OutPos=strrep(PosList,'pyr16_spot_',''); %format for output position in cytell images
end

%if no membrane markers are supplied ask for user input from AFRemoved list
if isempty(MemList)
    [a, ~]=listdlg('PromptString', 'Select Membrane Markers', 'SelectionMode', 'multiple', 'ListString', AFList);
    MemList=AFList(a);
end

%if no nuclear markers are supplied ask for user input from AFRemoved list
if isempty(NucList)
    [a, ~]=listdlg('PromptString', 'Select Nuclear Markers', 'SelectionMode', 'multiple', 'ListString', AFList);
    NucList=AFList(a);
end

%if no epithelial marker are supplied ask for user input from AFRemoved list
if isempty(EpiList)
    [a, ~]=listdlg('PromptString', 'Select Markers for Epithelial Mask', 'SelectionMode', 'multiple', 'ListString', AFList);
    NucList=AFList(a);
end

%make sure that the number of DAPI images equals number of poisitons
if length(DapiList) ~= length(PosList)
    error('Dapi Image Mismatch');
end

%print status updates to command line
fprintf(['Segmentation of: ' SlideDir ' ; ' num2str(length(PosList)) ' Positions; Unsharp Mask: ' ])
if unsharp==1
    fprintf('on')
else
    fprintf('off')
end
fprintf('\n')


%% Segmentation and Quantification for each position
for i=1:length(PosList)
    
    fprintf([OutPos{i} ': '])
    
    %% nuclear segmentation and generate SuperMembrane and binary Membrane mask
    
    %check to see if files exist; no: generate files ; yes: read saved images
    if ~exist([OutDir{5} 'NucMask_' OutPos{i} '.png']) || ~exist([OutDir{7} 'SuperMem_' OutPos{i} '.tif']) || ~exist([OutDir{10} 'MemMask_' OutPos{i} '.png'])
        
        fprintf('NucSeg, MemMask, & SuperMem; ')
        
        %format path and filenames for membrane markers
        mem=strcat(AFRemoved, MemList ,['_AFRemoved_' PosList{i} '.tif']);
        
        %format path and filenames for nuclear markers
        nuc=strcat(AFRemoved, NucList ,['_AFRemoved_' PosList{i} '.tif']);
        nuc{end+1}=DapiList{i}; %add DAPI image to marker set
        
        %calculate nuclear mask, supermembrane and membrane mask
        [mask, SuperMem, MemMask]=main(mem, DapiList{i}, nuc, unsharp);
               
        %write images generated from main
        imwrite(uint8(255*(mask>0)), [OutDir{5} 'NucMask_' OutPos{i} '.png'] )
        imwrite(SuperMem, [OutDir{7} 'SuperMem_' OutPos{i} '.tif'] ) %write 16 bit tiff
        imwrite(uint8(255*(MemMask>0)), [OutDir{10} 'MemMask_' OutPos{i} '.png'] ) %write 16 bit tiff
        
    else
        %read files if previously generated
        mask=imread( [OutDir{5} 'NucMask_' OutPos{i} '.png']);
        SuperMem=imread( [OutDir{7} 'SuperMem_' OutPos{i} '.tif']);
        MemMask=im2bw(imread([OutDir{10} 'MemMask_' OutPos{i} '.png']));
    end
    
    %make sure nuclear mask in binary
    mask=logical(mask);
    
    %remove blurred nuclear regions
    mask=mask.*blurimg2_batch(imread(DapiList{i}));
    
    %makes sure data is from Olympus, if not adjust constant for kernal
    %sizes in subsequent steps    
%     s= size(mask);
%     pixadj=1;
%     if s(1)~=2048 || s(2)~=2048 %check if images is from Cytell ~=2048
%         pixadj=3; %adjust for smaller pixel size if Cytell
%     end
    
    
    %% generate epithelial mask
    
    %check to see if files exist; no: generate files ; yes: read saved images
    if ~exist([OutDir{3} 'EpiMask_' OutPos{i} '.png'])
        fprintf('EpiMask; ')
        
        % format path and filenames for epithelial markers
        epi=strcat(AFRemoved, EpiList ,['_AFRemoved_' PosList{i} '.tif']); %get path and filenames for epithelial markers
        
        %generate epithelial mask
        epiMask=tubsegbatch(epi, 10, 0.00001, 0.01);
        
        
        %write image generated from tubesegbatch- forcing into 8 bit binary
        %(0 or 255)
        imwrite(uint8(255*(epiMask>0)), [OutDir{3} 'EpiMask_' OutPos{i} '.png'] ) %write cell borders
    else
        %read file if previously generated
        epiMask=logical(imread([OutDir{3} 'EpiMask_' OutPos{i} '.png']));
    end
    
    %% generate cell (re)segmentation and nuclear segmentation  images
    tic
    %check to see if files exist; no: generate files ; yes: read saved images
    if ~exist([OutDir{12} 'CellSegFinal_' OutPos{i} '.tif']) || ~exist([OutDir{11} 'NucMaskFinal_' OutPos{i} '.png'])
        
        fprintf('CellSeg; ')
        
        %format data for the Muc2 labeling
        if any(strcmp(AFList,'MUC2'))
            mucLoc = strcat(AFRemoved, ['MUC2_AFRemoved_' PosList{i} '.tif']); %get muc2 location
            %generate mask with Muc2 added and cells seperated by membranes
            EpiSep=Muc2Add (mucLoc, epiMask, MemMask, mask);
        elseif  any(strcmp(AFList,'Muc2'))
            mucLoc = strcat(AFRemoved, ['Muc2_AFRemoved_' PosList{i} '.tif']); %get muc2 location
            %generate mask with Muc2 added and cells seperated by membranes
            EpiSep=Muc2Add (mucLoc, epiMask, MemMask, mask);
        else %if Muc2 is not in the data set set MucLoc to empty
            mucLoc={};
        end
        
        %check to see if files exist; no: generate files ; yes: read saved images
        if ~exist([OutDir{1} 'L2_' OutPos{i} '.tif'])
            %if mucLoc is empty, run watershed
            if isempty(mucLoc)
                L2=watershed(imimposemin(double(MemMask),mask));
                imwrite(uint16(L2), [OutDir{1} 'L2_' OutPos{i} '.tif'] ) %write 16 bit tiff
                %if mucLoc exists, run muc2 based reseg
            else
                L2=shortBatchO(EpiSep, MemMask, mucLoc);
                imwrite(uint16(L2), [OutDir{1} 'L2_' OutPos{i} '.tif'] ) %write 16 bit tiff
            end
            
        else
            L2=imread( [OutDir{1} 'L2_' OutPos{i} '.tif']);
        end
        toc
        %check to see if files exist; no: generate files ; yes: read saved images
        if ~exist([OutDir{1} 'CellSeg_' OutPos{i} '.tif'])
            CellSeg= ReSegCells(L2, MemMask);
            imwrite(uint16(CellSeg), [OutDir{1} 'CellSeg_' OutPos{i} '.tif'] )
        else
            CellSeg=imread( [OutDir{1} 'CellSeg_' OutPos{i} '.tif']);
        end
        toc
        
        %check to see if files exist; no: generate files ; yes: read saved images
        if ~exist([OutDir{12} 'CellSegFinal_' OutPos{i} '.tif'])
            %CellSeg=imread( [OutDir{1} 'CellSeg_' OutPos{i} '.tif']);
            SuperMem=imread( [OutDir{7} 'SuperMem_' OutPos{i} '.tif']);
            [WatCellSeg, mask]=NucCountBatch(CellSeg, mask, epiMask, MemMask, mucLoc, imread(DapiList{i}), SuperMem);
            imwrite(uint16(WatCellSeg), [OutDir{12} 'CellSegFinal_' OutPos{i} '.tif'] )
            imwrite(uint8(255*(mask>0)), [OutDir{11} 'NucMaskFinal_' OutPos{i} '.png'] )
        else
            WatCellSeg=imread( [OutDir{12} 'CellSegFinal_' OutPos{i} '.tif']);
            mask=logical(imread([OutDir{11} 'NucMaskFinal_' OutPos{i} '.png']));
        end
    else
        WatCellSeg=imread( [OutDir{12} 'CellSegFinal_' OutPos{i} '.tif']);
        mask=logical(imread([OutDir{11} 'NucMaskFinal_' OutPos{i} '.png']));
     end
    toc
    
    
    
    
    %% Qunatification performed if specified
    if quantify ==1
        
        fprintf('Quant; ')
        
        %check to see if files exist; no: generate files ; yes: read saved files
        if exist([OutDir{13} 'PosStats_' OutPos{i} '.csv']) && exist([OutDir{4} 'Novlp_' OutPos{i} '.png'])
            fprintf('Reading file; ')
            Stats=readtable([OutDir{13} 'PosStats_' OutPos{i} '.csv']);
        else
            % run quantification
            [Stats, NoOvlp]=MxIF_quantify(i, WatCellSeg, AFRemoved, AFList, PosList, mask, MemMask, pixadj, epiMask, OutPos);
            %format data table and write
            names=cell2table(Stats.Properties.VariableNames);
            writetable(names, [OutDir{13} 'PosStats_' OutPos{i} '.csv'], 'WriteVariableNames', false);
            dlmwrite([OutDir{13} 'PosStats_' OutPos{i} '.csv'], table2array(Stats), '-append');
            imwrite(double(NoOvlp), [OutDir{4} 'Novlp_' OutPos{i} '.png'] )
        end
        
        if exist('AllStats','var') %if first position create file
            AllStats=vertcat(AllStats ,Stats);
        else
            AllStats=Stats; %vertically concatenate each positions
        end
        
        fprintf('\n')
        
    end
    
    
end

%after all fields are quantified write to file
if quantify ==1
    
    fprintf('Writing Stats File')
    fprintf('\n')
    
    %format data table and write
    names=cell2table(AllStats.Properties.VariableNames);
    writetable(names, [OutDir{8} 'Stats.csv'], 'WriteVariableNames', false);
    dlmwrite([OutDir{8} 'Stats.csv'], table2array(AllStats), '-append');
    
end

