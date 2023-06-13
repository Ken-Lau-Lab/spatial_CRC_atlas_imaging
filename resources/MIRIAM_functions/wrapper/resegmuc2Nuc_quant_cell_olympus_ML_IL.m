function resegmuc2Nuc_quant_cell_olympus_ML_IL_str(SlideDir ,quantify, start)
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
fprintf(['Segmentation of: ' SlideDir ' ; ' num2str(length(PosList)) ' Positions;\n' ])


%% Segmentation and Quantification for each position
for i=start:length(PosList)
    
    fprintf([OutPos{i} ': '])
    tic
    %make Stacks if they don't exist
    if ~exist([OutDir{15} OutPos{i} '_stack.tif'])
        
        fprintf(['Stack: ' OutPos{i} '\n'])
        %DapiImg=imread(DapiList{i});
        imwrite(imread(DapiList{i}), [OutDir{15} OutPos{i} '_stack.tif'], 'Compression', 'lzw');
        %[s1 s2]=size(DapiImg);
        %stack=zeros(s1, s2, length(AFList)+1);
        %stack(:, :, 1)=DapiImg;
        for j=1:(length(AFList))
            %stack(:,:,j+1)=imread([SlideDir '/AFRemoved/' AFList{j} '_AFRemoved_'  OutPos{i} '.tif']);
            imwrite(imread([SlideDir '/AFRemoved/' AFList{j} '_AFRemoved_'  OutPos{i} '.tif']), [OutDir{15} OutPos{i} '_stack.tif'], 'writemode', 'append', 'Compression', 'lzw');
            
        end
    end
    
    if ~exist([OutDir{14} 'epi_' OutPos{i} '_stack_Probabilities.png'])
        fprintf('No Epithelial Probability File\n')
        continue
    end
    
    if ~exist([OutDir{14} 'mem_' OutPos{i} '_stack_Probabilities.png'])
        fprintf('No Membrane/Nucleus Probability File\n')
        continue
    end
    
    %% nuclear segmentation and generate SuperMembrane and binary Membrane mask
    
    %check to see if files exist; no: generate files ; yes: read saved images
    if ~exist([OutDir{5} 'NucMask_' OutPos{i} '.png']) || ~exist([OutDir{7} 'SuperMem_' OutPos{i} '.tif']) || ~exist([OutDir{10} 'MemMask_' OutPos{i} '.png'])
        
        %fprintf('NucSeg, MemMask, & SuperMem; ')
        
        Probs=imread([OutDir{14} 'mem_' OutPos{i} '_stack_Probabilities.png']);
        
        %mem=Probs(:,:,1);
        %nuc=Probs(:,:,2);
        
        
        
        %imwrite(nuc, '/Users/etmckinley/Dropbox (VUMC)/scan_alpha/TCPS4A/SegQuant/NucMask/mem_009_nuc.png');
        
        %mem=mem>255*.6;
        
        
        
        %format path and filenames for membrane markers
        %mem=strcat(AFRemoved, MemList ,['_AFRemoved_' PosList{i} '.tif']);
        
        %format path and filenames for nuclear markers
        %nuc=strcat(AFRemoved, NucList ,['_AFRemoved_' PosList{i} '.tif']);
        %nuc{end+1}=DapiList{i}; %add DAPI image to marker set
        
        %calculate nuclear mask, supermembrane and membrane mask
        %[mask, SuperMem, MemMask]=main(mem, DapiList{i}, nuc);
        
        %write images generated from main
        mask=uint8(255*(Probs(:,:,2)>255*.6));
        imwrite(mask, [OutDir{5} 'NucMask_' OutPos{i} '.png'] )
        imwrite(Probs(:,:,1), [OutDir{7} 'SuperMem_' OutPos{i} '.tif'] ) %write 16 bit tiff
        imwrite(uint8(255*(Probs(:,:,1)>255*.6)), [OutDir{10} 'MemMask_' OutPos{i} '.png'] ) %write 16 bit tiff
        MemMask=im2bw(imread([OutDir{10} 'MemMask_' OutPos{i} '.png']));
        
    else
        %read files if previously generated
        mask=imread( [OutDir{5} 'NucMask_' OutPos{i} '.png']);
        SuperMem=imread( [OutDir{7} 'SuperMem_' OutPos{i} '.tif']);
        MemMask=im2bw(imread([OutDir{10} 'MemMask_' OutPos{i} '.png']));
    end
    
    %make sure nuclear mask in binary
    mask=logical(mask);
    mask=imfill(mask, 'holes');
    mask=imopen(mask, strel('disk',3));
    %mask=bwulterode(mask);
    
    %remove blurred nuclear regions
    mask=mask.*blurimg2_batch(imread(DapiList{i}));
    
    %makes sure data is from Olympus, if not adjust constant for kernal
    %sizes in subsequent steps
    s= size(mask);
    pixadj=1;
    if s(1)~=2048 || s(2)~=2048 %check if images is from Cytell ~=2048
        pixadj=3; %adjust for smaller pixel size if Cytell
    end
    
    
    %% generate epithelial mask from machine learning
    
    if~ exist([OutDir{3} 'EpiMask_' OutPos{i} '.png'])
        fprintf('EpiMask Processing; ')
        
        epiMask=imread([OutDir{14} 'epi_' OutPos{i} '_stack_Probabilities.png']);
        
        epiMask=ML_probability(epiMask,pixadj*0.01);
        imwrite(uint8(255*(epiMask>0)), [OutDir{3} 'EpiMask_' OutPos{i} '.png'] )
        epiMask=logical(imresize(epiMask, size(mask)));
        
    else
        %read file if previously generated
        epiMask=logical(imread([OutDir{3} 'EpiMask_' OutPos{i} '.png']));
        epiMask=logical(imresize(epiMask, size(mask)));
    end
    
    %check to see if files exist; no: generate files ; yes: read saved images
    %     if ~exist([OutDir{13} 'ML_' OutPos{i} '.jpg'])
    %         fprintf('EpiMask; ')
    %
    %         % format path and filenames for epithelial markers
    %         epi=strcat(AFRemoved, EpiList ,['_AFRemoved_' PosList{i} '.tif']); %get path and filenames for epithelial markers
    %
    %         %generate epithelial mask
    %         R_ML=ML_processing(epi);
    %         G_ML=mat2gray(imadjust(imread(epi{1}),stretchlim(imread(epi{1}),[0 .999])));  %adapthisteq(imread(epi{1}, 'ClipLimit', 0.2));
    %         B_ML=mat2gray(imadjust(imread(DapiList{i}),stretchlim(imread(DapiList{i}),[0 .9999])));
    %         %epiMask=tubsegbatch(epi, pixadj*10, pixadj*0.00001, pixadj*0.01);
    %         ML=(cat(3, R_ML, G_ML, B_ML));
    %         sz=size(ML);
    %         scale=1023/sz(2);
    %         ML=imresize(ML,scale);
    %         imwrite(ML, [OutDir{13} 'ML_' OutPos{i} '.jpg'] )
    %
    %
    %         %write image generated from tubesegbatch- forcing into 8 bit binary
    %         %(0 or 255)
    %         %imwrite(uint8(255*(epiMask>0)), [OutDir{3} 'EpiMask_' OutPos{i} '.png'] ) %write cell borders
    %     else
    %         %read file if previously generated
    %         %epiMask=logical(imread([OutDir{3} 'EpiMask_' OutPos{i} '.png']));
    %     end
    
    %Edges=edge(epiMask,'Sobel');
    %MemMask=imadd(Edges, MemMask);
    % Edges=edge(epiMask,'Sobel');
    %MemMask=bwmorph(MemMask, 'skel',Inf);
    MemMask=bwmorph(MemMask, 'thin',Inf);
    
    
    % MemMask=imadd(Edges, MemMask);
    
    %% generate cell (re)segmentation and nuclear segmentation  images
    %    tic
    %check to see if files exist; no: generate files ; yes: read saved images
    if ~exist([OutDir{12} 'CellSegFinal_' OutPos{i} '.tif']) || ~exist([OutDir{11} 'NucMaskFinal_' OutPos{i} '.png'])
        
        fprintf('CellSeg; ')
        
        mucLoc={};
        %          EpiSep=0;
        %         %format data for the Muc2 labeling
        %         if any(strcmp(AFList,'MUC2'))
        %             mucLoc = strcat(AFRemoved, ['MUC2_AFRemoved_' PosList{i} '.tif']); %get muc2 location
        %             %generate mask with Muc2 added and cells seperated by membranes
        %             EpiSep=Muc2Add (mucLoc, epiMask, MemMask, mask);
        %         elseif  any(strcmp(AFList,'Muc2'))
        %             mucLoc = strcat(AFRemoved, ['Muc2_AFRemoved_' PosList{i} '.tif']); %get muc2 location
        %             %generate mask with Muc2 added and cells seperated by membranes
        %             EpiSep=Muc2Add (mucLoc, epiMask, MemMask, mask);
        %         else %if Muc2 is not in the data set set MucLoc to empty
        %             mucLoc={};
        %         end
        %
        %         if max(EpiSep(:))==0
        %             mucLoc={};
        %         end
        
        
        
        
        %check to see if files exist; no: generate files ; yes: read saved images
        if ~exist([OutDir{1} 'L2_' OutPos{i} '.tif'])
            %if mucLoc is empty, run watershed
            L2=imcomplement (epiMask);
            L2=im2bw(imadd(im2bw(L2),im2bw(MemMask)));
            L2=uint8(L2);
            
            %watseg with nuclei as basins
            L2 = watershed(imimposemin(L2,mask),4);
            
            L2=double(L2);
            %L2=im2bw(L2);
            
            % return only those in tube mask
            L2=L2.*epiMask;
            imwrite(uint16(L2), [OutDir{1} 'L2_' OutPos{i} '.tif'] ) %write 16 bit tiff
            %             if isempty(mucLoc)
            %                 L2=watershed(imimposemin(double(MemMask),mask));
            %                 imwrite(uint16(L2), [OutDir{1} 'L2_' OutPos{i} '.tif'] ) %write 16 bit tiff
            %                 %if mucLoc exists, run muc2 based reseg
            %             else
            %                 L2=shortBatch(EpiSep, MemMask, mucLoc);
            %                 imwrite(uint16(L2), [OutDir{1} 'L2_' OutPos{i} '.tif'] ) %write 16 bit tiff
            %             end
            
        else
            L2=imread( [OutDir{1} 'L2_' OutPos{i} '.tif']);
        end
        %        toc
        %check to see if files exist; no: generate files ; yes: read saved images
        if ~exist([OutDir{1} 'CellSeg_' OutPos{i} '.tif'])
            MemMask=im2bw(imread([OutDir{10} 'MemMask_' OutPos{i} '.png']));
            CellSeg= ReSegCells(L2, MemMask);
            imwrite(uint16(CellSeg), [OutDir{1} 'CellSeg_' OutPos{i} '.tif'] )
        else
            CellSeg=imread( [OutDir{1} 'CellSeg_' OutPos{i} '.tif']);
        end
        %   toc
        
        %check to see if files exist; no: generate files ; yes: read saved images
        if ~exist([OutDir{12} 'CellSegFinal_' OutPos{i} '.tif'])
            CellSeg=imread( [OutDir{1} 'CellSeg_' OutPos{i} '.tif']);
            SuperMem=imread( [OutDir{7} 'SuperMem_' OutPos{i} '.tif']);
            Probs=imread([OutDir{14} 'mem_' OutPos{i} '_stack_Probabilities.png']);
            
            [WatCellSeg, mask]=NucCountBatch(CellSeg, mask, epiMask, MemMask, mucLoc, Probs(:,:,2), SuperMem);
            
            
     
            
           % L=bwlabel(CellSeg,4);
%             L=closesmallholes(CellSeg, 50);
%             L=bwlabel(L,4);
%             
%             %final filering step that can maybe be ignored
%             %mainly look for large objects, and then, first get rif of objects that are very weirdly shapped (circle metric<0.15)
%             %then keep only those that are "long" (eccen) cells.
%             
%             stats = regionprops(L,'Area');
%             
%             % large objects
%             Areas=transpose(cell2mat({stats.Area}));
%             large_ones=find (Areas>2500);
%             large= ismember(L,large_ones).*L;
%             
%             % circle metric
%             bw = bwlabel (large,4);
%             [B] = bwboundaries(bw, 4,'noholes');
%             metric_b=calculate_metric(B,bw);
%             erase=find(metric_b<0.15);
%             remain=find(metric_b>=0.15);
%             erase1=ismember(bw,erase).*bw;
%             remain1=ismember(bw,remain).*bw;
%             
%             % keep long cells
%             L=bwlabel(remain1,4);
%             stats = regionprops(L,'Eccentricity');
%             eccen=transpose(cell2mat({stats.Eccentricity}));
%             erasex=find (eccen<0.5);
%             erase2=ismember(L,erasex).*L;
%             
%             %add the two erases
%             erase1=im2bw(imadd(im2bw(erase1),im2bw(erase2)));
%             
%             %subtract all to erase
%             mem5=im2bw(imsubtract(im2bw(CellSeg),im2bw(erase1)));
%             
%             s=size(mem5);
%             
%             if s(1)~=2048 || s(2)~=2048 %check if images is from Cytell ~=2048
%                 mem5=closesmallholes(mem5, 50); %adjust for smaller pixel size if Cytell
%             else
%                 mem5=closesmallholes(mem5, 20);
%             end
%             
%             %fill in small holes in cells
%             WatCellSeg=bwlabel(mem5,4);
%             
%             
%             
%             
%             
%             % WatCellSeg=CellSeg;
            
            
            WatCellSeg=closesmallholes(WatCellSeg, 10);
            WatCellSeg(epiMask==0)=0;
            WatCellSeg=logical(WatCellSeg);
            WatCellSeg=bwareaopen(WatCellSeg,15);
            WatCellSeg=bwlabel(WatCellSeg,4);
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
    
    
    
    
    %WatCellSeg=watershed(imimposemin(double(MemMask),mask));
    %imwrite(uint16(WatCellSeg), [OutDir{12} 'CellSegFinal_' OutPos{i} '.tif'] )
    % toc
    
    
    
    
    %% Qunatification performed if specified
    if quantify ==1
        
        
        
        %check to see if files exist; no: generate files ; yes: read saved files
        if exist([OutDir{13} 'PosStats_' OutPos{i} '.csv']) && exist([OutDir{4} 'Novlp_' OutPos{i} '.png'])
            %fprintf('Reading file; ')
            %             Stats=readtable([OutDir{13} 'PosStats_' OutPos{i} '.csv']);
        else
            if max(WatCellSeg(:))==0
                fprintf('\n')
                continue
            end
            
            % run quantification
            fprintf('Quant; ')
            [Stats, NoOvlp]=MxIF_quantify(i, WatCellSeg, AFRemoved, AFList, PosList, mask, MemMask, pixadj, epiMask, OutPos);
            %format data table and write
            names=cell2table(Stats.Properties.VariableNames);
            writetable(names, [OutDir{13} 'PosStats_' OutPos{i} '.csv'], 'WriteVariableNames', false);
            dlmwrite([OutDir{13} 'PosStats_' OutPos{i} '.csv'], table2array(Stats), '-append');
            imwrite(double(NoOvlp), [OutDir{4} 'Novlp_' OutPos{i} '.png'] )
        end
        
        %         if exist('AllStats','var') %if first position create file
        %             AllStats=vertcat(AllStats ,Stats);
        %         else
        %             AllStats=Stats; %vertically concatenate each positions
        %         end
        
        
        
    end
    toc
    fprintf('\n')
    
end

%after all fields are quantified write to file
% if quantify ==1
%
%     fprintf('Writing Stats File')
%     fprintf('\n')
%
%     %format data table and write
%     names=cell2table(AllStats.Properties.VariableNames);
%     writetable(names, [OutDir{8} 'Stats.csv'], 'WriteVariableNames', false);
%     dlmwrite([OutDir{8} 'Stats.csv'], table2array(AllStats), '-append');
%
% end

