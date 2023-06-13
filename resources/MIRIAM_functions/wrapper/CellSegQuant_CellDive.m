function CellSegQuant_CellDive(SlideDir, quantify, shape, stroma, tumor, annotate, start, threshold, threshold_file, ilastik_epi_file, ilastik_mem_file, stack_markers)
%Eliot McKinley - 2018-11-26
%eliot.t.mckinley.1@vanderbilt.edu
%
% No returned output - function saves image of objects and quantification
% of AFremoved markers
%
%SlideDir=directory for slide containing AFRemoved images folder and
%Registered images folder - assumes Round 001 is baseline and all fies are
%.tif - if no folder is supplied- then ui asks for directory
%
%quantify= whether or not to quantify 1=yes, 0=no
%
%stroma= whether or not to segment stroma 1= yes, 0-no
%
%start= which image to start processing

%% Custom functions needed
%SegDirFormatting-creates directories to save the output of this function
%blurimg2_batch - finds blurry spots on DAPI images and excludes from analysis
%ML_probability - creates epithelial mask from probability map
%MaskFiltration - filtering for epithelial mask generation
%ReSegCells- generates resegmented cell mask
%NucCountBatch - resegments cells with multiple nuclei- generating new cellmask and nuclear mask
%MxIF_quantify - quantifies marker staining intensity
%stromal_nuclei_segmentation - segments stromal nuclei
%BBoxCalc - get bounding box of cell
%reseg - re-segment cells with internal membranes
%extend - extend line segments on both ends
%intersectfirst - connect line segments
%trimtree - take branches off line segments
%finalconnect_2 - extend from spurs in line
%finalconnect - connect any oustanding endpoint pairs

%% Parse Directory supplied for cell segmentation

pixadj=3;


%if no directory is supplied ask for user input
if isempty(SlideDir)
    SlideDir=uigetdir;
end

% get formatting for AFRemoved, DAPI, and output directories
[AFRemoved, DAPI, OutDir]=SegDirFormatting_COLON_MAP(SlideDir);

%get AFRemoved images
AFList=imageDatastore(AFRemoved);
%AFList=dir([AFRemoved '*.tif']);
AFList={AFList.Files};
AFList=AFList{1,1};
%PosList=AFList;

AF_delimiter='_AFRemoved_'; % the string in between the Marker and position ID
AFList=split(AFList, AF_delimiter);
PosList=AFList(:,2);
AFList=split(AFList(:,1), '/AFRemoved Tiles/');
AFList=AFList(:,2);
AFList=unique(AFList); %get only marker names

PosList=strrep(PosList, ".tif", "");
PosList=unique(PosList);
OutPos=PosList; %format for output position

%get DAPI images from first round of imaging
DapiList=imageDatastore(DAPI);
DapiList={DapiList.Files};
DapiList=DapiList{1,1};

%make sure that the number of DAPI images equals number of poisitons
if length(DapiList) ~= length(PosList)
    error('Dapi Image Mismatch');
end

%get thresholded images
ThreshList=[];
if threshold == 1
    Slide=split(SlideDir, '/scan_alpha/');
    Slide=Slide{2};
    %ThreshList=imageDatastore(join(["/Users/etmckinley/Dropbox (VUMC)/T cell quantification/Images/" Slide "/*threshold.png"], "")).Files;
    ThreshList=imageDatastore(join([threshold_file Slide "/*threshold.png"], "")).Files;
    ThreshList=split(ThreshList, join([Slide "/"], ""));
    ThreshList=ThreshList(:,2);
    ThreshList=split(ThreshList, "_region");
    ThreshList=unique(ThreshList(:,1));
end


%print status updates to command line
fprintf(['Segmentation of: ' SlideDir ' ; ' num2str(length(PosList)) ' Positions;\n' ])



%% Segmentation and Quantification for each position

%Note: you can run this using parfor when Ilastik is not being used to
%generate probability files.
for i=  start: length(PosList)

    
    fprintf([ SlideDir ' ' OutPos{i} ': '])
    tic
    %make Stacks of AFRemoved images and Dapi if they don't exist
    if ~exist([OutDir{15} OutPos{i} '_stack.tif'], 'file')
        %seg_markers=readtable('/Volumes/CellDive/COLON MAP Segmentation/Segmentation markers.csv', 'Delimiter', ',');
        seg_markers = stack_markers;
        %seg_markers=readtable('/Users/etmckinley/Dropbox (VUMC)/Research/Projects/COLON MAP/Cody Slides/Segmentation markers.csv', 'Delimiter', ',');
        %seg_markers=seg_markers.Marker;
        %seg_markers=seg_markers(2:length(seg_markers));
        fprintf(['Stack: ' OutPos{i} '\n'])
        imwrite(imread(DapiList{i}), [OutDir{15} OutPos{i} '_stack.tif'], 'Compression', 'lzw'); %write out DAPI
        for j=1:(length(seg_markers)) %loop through the AFRemoved images and append to tiff stack
            %imwrite(imread([SlideDir '/AFRemoved/' AFList{j} '_AFRemoved_pyr16_spot_'  OutPos{i} '.tif']), [OutDir{15} OutPos{i} '_stack_all.tif'], 'writemode', 'append', 'Compression', 'lzw');
            imwrite(imread([SlideDir '/AFRemoved Tiles/' seg_markers{j} '_AFRemoved_' OutPos{i} '.tif']), [OutDir{15} OutPos{i} '_stack.tif'], 'writemode', 'append', 'Compression', 'lzw');
            
        end
    end
    

        
     % Check for Epithelial Probability file from Ilastik
    if ~exist([OutDir{14} 'epi_' OutPos{i} '_stack_Probabilities.png'], 'file')
        proj_path = ilastik_epi_file;
        if ~exist(proj_path, 'file')
            continue
        end
        file_info= dir([OutDir{15} OutPos{i} '_stack.tif']);
        if file_info.bytes < 2000000 %check to see if stack is too small
            fprintf('Stack too small\n')
            continue
        end
        
        fprintf('Generating Epithelial Probability File in Headless Ilastik\n')
        outPath = [OutDir{14} 'epi_' OutPos{i} '_stack_Probabilities.png'];
        stackPath = [OutDir{15} OutPos{i} '_stack.tif'];
        % note, you need to be pointing to the correct ilastik file below
        % and you can't have any spaces in the paths
        command1 = join(['/Applications/ilastik-1.4.0rc5-OSX.app/Contents/ilastik-release/run_ilastik.sh --headless --project="' proj_path '" --output_format=png --output_filename_format="' outPath '"  "' stackPath '"'],"");
        [ ~, ~ ] =  system(command1);        
        
    end
    
       %Check for Membrane Probability file from Ilastik
    if ~exist([OutDir{14} 'mem_' OutPos{i} '_stack_Probabilities.png'], 'file')
        proj_path = ilastik_mem_file;
        if ~exist(proj_path, 'file')
            continue
        end
        
        
        fprintf('Generating Membrane Probability File in Headless Ilastik\n')
            outPath = [OutDir{14} 'mem_' OutPos{i} '_stack_Probabilities.png'];
            stackPath = [OutDir{15} OutPos{i} '_stack.tif'];
        command1 = join(['/Applications/ilastik-1.4.0rc5-OSX.app/Contents/ilastik-release/run_ilastik.sh --headless --project="' proj_path '" --output_format=png --output_filename_format="' outPath '"  "' stackPath '"'],"");
        [ ~, ~ ] =  system(command1);        
                
    end
    
    
    %% generate epithelial mask from machine learning
    if~ exist([OutDir{3} 'EpiMask_' OutPos{i} '.png'], 'file')
        fprintf('EpiMask Processing; ')
        if ~exist([OutDir{14} 'epi_' OutPos{i} '_stack_Probabilities.png'], 'file')
            continue
        end
        epiMask=imread([OutDir{14} 'epi_' OutPos{i} '_stack_Probabilities.png']);
        if mean(mean(epiMask(:,:,1)))<2 || mean(mean(epiMask(:,:,2)))== 0
            fprintf('Blank Position\n')
            continue
        end
        
        epiMask=ML_probability(epiMask,3*0.01, .45); %create epithelial mask from probability map
        imwrite(uint8(255*(epiMask>0)), [OutDir{3} 'EpiMask_' OutPos{i} '.png'] )
        %epiMask=logical(imresize(epiMask, size(mask)));
        
    else
        %read file if previously generated
        epiMask=logical(imread([OutDir{3} 'EpiMask_' OutPos{i} '.png']));
        if sum(sum(epiMask)) < 15000
            fprintf('Epithelial Area Too Small\n')
            continue
        end
        %epiMask=logical(imresize(epiMask, size(mask)));
    end
    
 
    
    

   
    
    
    
%     if~ exist([OutDir{3} 'EpiMask_' OutPos{i} '.png'], 'file')
%         fprintf('EpiMask Processing; ')
%         
%         epiMask=imread([OutDir{14} 'epi_' OutPos{i} '_stack_Probabilities.png']);
%         if mean(mean(epiMask(:,:,1)))<2 || mean(mean(epiMask(:,:,2)))== 0
%             fprintf('Blank Position\n')
%             continue
%         end
%         
%         epiMask=ML_probability(epiMask,3*0.01, .45); %create epithelial mask from probability map
%         imwrite(uint8(255*(epiMask>0)), [OutDir{3} 'EpiMask_' OutPos{i} '.png'] )
%         %epiMask=logical(imresize(epiMask, size(mask)));
%         
%     else
%         %read file if previously generated
%         epiMask=logical(imread([OutDir{3} 'EpiMask_' OutPos{i} '.png']));
%         if sum(sum(epiMask)) < 15000
%             fprintf('Epithelial Area Too Small\n')
%             continue
%         end
%         %epiMask=logical(imresize(epiMask, size(mask)));
%     end
    
        
        
        %% nuclear segmentation and generate SuperMembrane and binary Membrane mask
        
        if exist([OutDir{14} 'epi_' OutPos{i} '_stack_Probabilities.png'], 'file')
            
            %check to see if files exist; no: generate files ; yes: read saved images
            if ~exist([OutDir{5} 'NucMask_' OutPos{i} '.png'], 'file') || ~exist([OutDir{7} 'SuperMem_' OutPos{i} '.tif'], 'file') || ~exist([OutDir{10} 'MemMask_' OutPos{i} '.png'], 'file')
                %Read in probability image for membrane and nucleus
                if ~exist([OutDir{14} 'mem_' OutPos{i} '_stack_Probabilities.png'], 'file')
                    continue
                end
                Probs=imread([OutDir{14} 'mem_' OutPos{i} '_stack_Probabilities.png']);
                
                
                mask=uint8(255*(Probs(:,:,2)>255*.6)); %set nuclear probability >0.6 as nuclear mask
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
            
            %fill in small holes and smooth
            mask=imfill(mask, 'holes');
            mask=imopen(mask, strel('disk',3));
            
            
            %remove blurred nuclear regions
            mask=mask.*blurimg2_batch(imread(DapiList{i}));
            
            %makes sure data is from Olympus, if not adjust constant for kernal
            %sizes in subsequent steps
            s= size(mask);
            
            
            
            
            %thin the membrane borders prior to initial watershed
            MemMask=bwmorph(MemMask, 'thin',Inf);
            
            %% generate cell (re)segmentation and nuclear segmentation  images
            %check to see if files exist; no: generate files ; yes: read saved images
            if ~exist([OutDir{12} 'CellSegFinal_' OutPos{i} '.tif'], 'file') || ~exist([OutDir{11} 'NucMaskFinal_' OutPos{i} '.png'], 'file')
                
                fprintf('CellSeg; ')
                
                %check to see if files exist; no: generate files ; yes: read saved images
                if ~exist([OutDir{1} 'L2_' OutPos{i} '.tif'], 'file')
                    %if mucLoc is empty, run watershed
                    L2=imcomplement (epiMask);
                    L2=im2bw(imadd((L2),(MemMask)));
                    L2=uint8(L2);
                    
                    %watershed segmentation with nuclei as basins
                    L2 = watershed(imimposemin(L2,mask),4);
                    
                    L2=double(L2);
                    
                    % return only those in the epithelial mask
                    L2=L2.*epiMask;
                    imwrite(uint16(L2), [OutDir{1} 'L2_' OutPos{i} '.tif'] ) %write 16 bit tiff
                    
                    
                else
                    L2=imread( [OutDir{1} 'L2_' OutPos{i} '.tif']);
                end
                
                %check to see if files exist; no: generate files ; yes: read saved images
                if ~exist([OutDir{1} 'CellSeg_' OutPos{i} '.tif'], 'file')
                    MemMask=im2bw(imread([OutDir{10} 'MemMask_' OutPos{i} '.png']));
                    CellSeg= ReSegCells(L2, MemMask);
                    imwrite(uint16(CellSeg), [OutDir{1} 'CellSeg_' OutPos{i} '.tif'] )
                else
                    CellSeg=imread( [OutDir{1} 'CellSeg_' OutPos{i} '.tif']);
                end
                %   toc
                
                %check to see if files exist; no: generate files ; yes: read saved images
                if ~exist([OutDir{12} 'CellSegFinal_' OutPos{i} '.tif'], 'file')
                    CellSeg=imread( [OutDir{1} 'CellSeg_' OutPos{i} '.tif']);
                    SuperMem=imread( [OutDir{7} 'SuperMem_' OutPos{i} '.tif']);
                    Probs=imread([OutDir{14} 'mem_' OutPos{i} '_stack_Probabilities.png']);
                    %check for cells with multiple nuclei and re-segment if they exist
                    [WatCellSeg, mask]=NucCountBatch(CellSeg, mask, epiMask, MemMask, {}, Probs(:,:,2), SuperMem);
                    WatCellSeg=im2bw(WatCellSeg);
                    
                    %fill small crosses
                    %WatCellSeg=closesmallholes(WatCellSeg, 4);
                    filter=zeros(5,5);
                    filter(2,3)=1;
                    filter(3,3)=1;
                    filter(3,2)=1;
                    filter(3,4)=1;
                    filter(4,3)=1;
                    
                    
                    hitmiss=bwhitmiss(WatCellSeg, imcomplement(filter));
                    spots=bwareaopen(hitmiss,2);
                    diff=hitmiss-spots;
                    
                    diff=conv2(diff, filter, 'same');
                    %imshowpair(hitmiss, spots)
                    
                    WatCellSeg=diff+WatCellSeg;
                    
                    
                    
                    %WatCellSeg2= conv2(imcomplement(WatCellSeg),filter, 'same');
                    %set non-epithelial pixels to zero
                    WatCellSeg(epiMask==0)=0;
                    WatCellSeg=logical(WatCellSeg);
                    %opening to clean up segmentation
                    WatCellSeg=bwareaopen(WatCellSeg,15);
                    %label cells
                    
                    WatCellSeg=bwlabel(WatCellSeg,4);
                    %write out data
                    imwrite(uint16(WatCellSeg), [OutDir{12} 'CellSegFinal_' OutPos{i} '.tif'] )
                    imwrite(uint8(255*(mask>0)), [OutDir{11} 'NucMaskFinal_' OutPos{i} '.png'] )
                else
                    %read in data
                    WatCellSeg=imread( [OutDir{12} 'CellSegFinal_' OutPos{i} '.tif']);
                    mask=logical(imread([OutDir{11} 'NucMaskFinal_' OutPos{i} '.png']));
                end
            else
                WatCellSeg=imread( [OutDir{12} 'CellSegFinal_' OutPos{i} '.tif']);
                mask=logical(imread([OutDir{11} 'NucMaskFinal_' OutPos{i} '.png']));
            end
            
            
            
            
            
            
            
            
            
            %% Quantification performed if specified
            
            if quantify ==1   
            
            %check to see if files exist; no: generate files ; yes: read saved files
            if exist([OutDir{13} 'PosStats_' OutPos{i} '.csv'], 'file') && exist([OutDir{4} 'Novlp_' OutPos{i} '.png'], 'file') && annotate==0
                %fprintf('Reading file; ')
                %             Stats=readtable([OutDir{13} 'PosStats_' OutPos{i} '.csv']);
            else
                %skip if the field is empty
                if max(WatCellSeg(:))==0
                    fprintf('\n')
                    continue
                end
                
                % run quantification
                fprintf('Quant; ')
                
                if tumor==0 && annotate ==0
                    
                    if length(unique(WatCellSeg)) <20
                        continue
                    end
                    [Stats, NoOvlp]=MxIF_quantify(i, WatCellSeg, AFRemoved, AFList, PosList, mask, MemMask, [], [], OutPos);
                end
                if tumor==1
                    if length(unique(WatCellSeg)) <20
                        continue
                    end
                    %if ~exist([OutDir{14} 'tum_' OutPos{i} '_stack_Probabilities.png'], 'file')
                    %    fprintf('No Tumor Probability File\n')
                    %    continue
                    %end
                    if  ~exist([OutDir{16} 'TumorMask_' OutPos{i} '.png'], 'file')
                        %[tumorMask,~ , ~]=imread([OutDir{14} 'tum_' OutPos{i} '_stack_Probabilities.png']);
                        tumorMask = false(size(mask,1), size(mask,2)); 
                    else 
                        tumorMask= ( imread([OutDir{16} 'TumorMask_' OutPos{i} '.png']));
                        %tumorMask=uint8(255*(tumorMask(:,:,1)>255*.5)); %set tumor probability >0.5 as tumor mask
                        %tumorMask=ML_probability(tumorMask,pixadj*0.01,.5); %create epithelial mask from probability map
                        %if ~exist([OutDir{16} 'TumorMask_' OutPos{i} '.png'], 'file')
                        %    imwrite(uint8(255*(tumorMask>0)), [OutDir{16} 'TumorMask_' OutPos{i} '.png'] )
                        %end
                     end  
                        [Stats, NoOvlp]=MxIF_quantify(i, WatCellSeg, AFRemoved, AFList, PosList, mask, MemMask, tumorMask, [], OutPos);
                    
                end
                if annotate==1
                 if ~exist([OutDir{18} 'Annotated_' OutPos{i} '.png'], 'file') 
                     if  exist([OutDir{4} 'Novlp_' OutPos{i} '.png'], 'file')
                         continue
                     end
                    [Stats, NoOvlp]=MxIF_quantify(i, WatCellSeg, AFRemoved, AFList, PosList, mask, MemMask, [], [], OutPos);
                    Annotation=table(zeros(height(Stats),1));
                    Annotation.Properties.VariableNames={'Top'}; %rename table variable
                    Stats=[Stats Annotation];
                    Annotation.Properties.VariableNames={'Bottom'}; %rename table variable
                    Stats=[Stats Annotation];
                 end
                 if exist([OutDir{18} 'Annotated_' OutPos{i} '.png'], 'file')
                    annotate_im=imread([OutDir{18} 'Annotated_' OutPos{i} '.png']);
                    [Stats, NoOvlp]=MxIF_quantify(i, WatCellSeg, AFRemoved, AFList, PosList, mask, MemMask, [], annotate_im, OutPos);

                 end
                end
                    
                
                %format data table and write
                %names=cell2table(Stats.Properties.VariableNames);
                writetable(Stats, [OutDir{13} 'PosStats_' OutPos{i} '.csv'], 'WriteVariableNames', true);
                %dlmwrite([OutDir{13} 'PosStats_' OutPos{i} '.csv'], (Stats), '-append');
                imwrite(double(NoOvlp), [OutDir{4} 'Novlp_' OutPos{i} '.png'] )
            end
            
        end
    end
    
    %% Stromal Quantification
    if stroma==1 && quantify ==1
        
        if  ~exist([OutDir{13} 'StrPosStats_' OutPos{i} '.csv'], 'file')
        fprintf('\nStromal Quant; ')
        epiMask=imread([OutDir{14} 'epi_' OutPos{i} '_stack_Probabilities.png']);
        if mean(mean(epiMask(:,:,1)))<2 || mean(mean(epiMask(:,:,2)))== 0
            fprintf('Blank Position\n')
            continue
        end
        
        strMask=ML_str_probability(epiMask,3*0.01, .45); %create epithelial mask from probability map
        
        
        NUCimg=imread(DapiList{i});
        
        %NUCim = imopen(NUCimg,strel('disk',3));
        %NUCopen= imsharpen(NUCimg, 'Radius',5,'Amount',1);
        NUCimg= imguidedfilter(NUCimg, 'NeighborhoodSize', [6 6]);
        %imshow(NUCopen)
        %imshowpair(NUCimg, NUCopen)
        %Fast peak finder
        filt2 = (fspecial('gaussian', 8, 2));
        fastPOS = FastPeakFind(NUCimg, 1, filt2, 2);
        posSPOTS = reshape(fastPOS,2,[]);
        plotSPOTS = zeros(size(NUCimg,1) , size(NUCimg,2));
       
        %Create binary image from the identified peaks
        for Spots=1:length(posSPOTS)
            colSpot = posSPOTS(1,Spots);
            rowSpot = posSPOTS(2,Spots);
            if colSpot > size(NUCimg,2) || rowSpot >  size(NUCimg,1)
                continue
            end
            plotSPOTS(rowSpot,colSpot)= 1; 
        end
           
        NUCspots = imbinarize(plotSPOTS,0.01);
        %NUCspots = NUCspots & tissueMask;
        NUCspots = imdilate(NUCspots,strel('disk',5));
        %NUCblur = imgaussfilt(mat2gray(NUCplots),5);        
        
        %imshowpair(NUCimg, NUCplots)
        NUCspots(strMask==false)=false;
        
        
        water=watershed( imcomplement( NUCspots),4);
        water(strMask == false) =0;
        
        imwrite(uint16(water), [OutDir{12} 'StrCellSegFinal2_' OutPos{i} '.tif'] )
        
        
        if tumor==1
            if length(unique(water)) <20
                continue
            end
            %if ~exist([OutDir{14} 'tum_' OutPos{i} '_stack_Probabilities.png'], 'file')
            %    fprintf('No Tumor Probability File\n')
            %    continue
            %end
            if  ~exist([OutDir{16} 'TumorMask_' OutPos{i} '.png'], 'file')
                %[tumorMask,~ , ~]=imread([OutDir{14} 'tum_' OutPos{i} '_stack_Probabilities.png']);
                tumorMask = false(size(water,1), size(water,2));
            else
                tumorMask= ( imread([OutDir{16} 'TumorMask_' OutPos{i} '.png']));
                %tumorMask=uint8(255*(tumorMask(:,:,1)>255*.5)); %set tumor probability >0.5 as tumor mask
                %tumorMask=ML_probability(tumorMask,pixadj*0.01,.5); %create epithelial mask from probability map
                %if ~exist([OutDir{16} 'TumorMask_' OutPos{i} '.png'], 'file')
                %    imwrite(uint8(255*(tumorMask>0)), [OutDir{16} 'TumorMask_' OutPos{i} '.png'] )
                %end
            end
            %[Stats, NoOvlp]=MxIF_quantify(i, WatCellSeg, AFRemoved, AFList, PosList, mask, MemMask, tumorMask, [], OutPos);
            
        end
        if tumor == 0
            tumorMask = [];
        end
        
         %quantify markers in cells and write out data
            if max(max(water)) >9
                strStats=MxIF_quantify_stroma2(i, water, AFRemoved, AFList, PosList, NUCimg,  tumorMask, OutPos);
                                       
                if ~isempty(strStats)
                
                    writetable(strStats, [OutDir{13} 'StrPosStats_' OutPos{i} '.csv'], 'WriteVariableNames', true);
                    %names=cell2table(strStats.Properties.VariableNames);
                    %writetable(names, [OutDir{13} 'StrPosStats_' OutPos{i} '.csv'], 'WriteVariableNames', false);
                    %dlmwrite([OutDir{13} 'StrPosStats_' OutPos{i} '.csv'], table2array(strStats), '-append');
                    %imwrite(double(strNoOvlp), [OutDir{4} 'StrNovlp_' OutPos{i} '.png'] )
                end
            end
        %end
        
        %imwrite(uint8(255*(stromal_nuclei>0)), [OutDir{11} 'StrNucMaskFinal_' OutPos{i} '.png'] )
        
        %imagesc(water)
        %imwrite( imcomplement(imbinarize(water,0)), 'nuc find test.png')
        
        
        %check to see if stromal nuclei probability map exists
%         if ~exist([OutDir{14} 'str_' OutPos{i} '_stack_Probabilities.png'], 'file')
%             fprintf('\nStromal Quant; ')
%             fprintf('No Stromal Probability File\n')
%             continue
%         end
        %if stats and segmentation already exist, then skip
        %if exist([OutDir{13} 'StrPosStats_' OutPos{i} '.csv'], 'file') && exist([OutDir{4} 'StrNovlp_' OutPos{i} '.png'], 'file')
            
%         else
%             fprintf('\nStromal Quant; ')
%             %run function to segment stromal nuclei
%             stromal_nuclei=stromal_nuclei_segmentation(imread([OutDir{14} 'str_' OutPos{i} '_stack_Probabilities.png']));
%             stromal_nuclei(epiMask==1)=0;
%             stromal_grow=imdilate(stromal_nuclei, strel('disk',3)); %dilate nuclei a bit
%             
%             stromal_label= watershed(imimposemin(uint8(stromal_grow),stromal_nuclei),4); %watershed on dilated cells with nuceli as seed points
%             stromal_label(stromal_grow==0)=0; %set all non-cellular pixels to zero
%             %write out results
%             imwrite(uint16(stromal_label), [OutDir{12} 'StrCellSegFinal_' OutPos{i} '.tif'] )
%             imwrite(uint8(255*(stromal_nuclei>0)), [OutDir{11} 'StrNucMaskFinal_' OutPos{i} '.png'] )
%             
%             %quantify markers in cells and write out data
%             if max(max(stromal_label)) >9
%                 [strStats]=MxIF_quantify_stroma(i, stromal_label, AFRemoved, AFList, PosList, OutPos);
%                 if ~isempty(strStats)
%                 
%                     writetable(strStats, [OutDir{13} 'StrPosStats_' OutPos{i} '.csv'], 'WriteVariableNames', true);
%                     %names=cell2table(strStats.Properties.VariableNames);
%                     %writetable(names, [OutDir{13} 'StrPosStats_' OutPos{i} '.csv'], 'WriteVariableNames', false);
%                     %dlmwrite([OutDir{13} 'StrPosStats_' OutPos{i} '.csv'], table2array(strStats), '-append');
%                     %imwrite(double(strNoOvlp), [OutDir{4} 'StrNovlp_' OutPos{i} '.png'] )
%                 end
%             end
         end
        
    end
    
    
    %% Shape Pre-processing
    if shape==1
        %load the final cell segmentation image, extract the cells, and
        %save to .mat file
        if exist([OutDir{12} 'CellSegFinal_' OutPos{i} '.tif'], 'file') && ~exist([OutDir{17} 'CellShape_' OutPos{i} '.mat' ], 'file')
            fprintf('\nCell Shape Pre-Processing; ')
            CellImages=imread([OutDir{12} 'CellSegFinal_' OutPos{i} '.tif']);
            CellImages=cell_shape_images(CellImages);
            save([OutDir{17} 'CellShape_' OutPos{i} '.mat' ], 'CellImages')
        end
        
        %concatenate all the position cell shapes and
        if i == length(PosList)
            if ~exist([OutDir{17} 'autoencoder.mat'], 'file') && ~exist([OutDir{17} 'encoded_cells.csv'], 'file')
                fprintf('\nTraining Autoencoder; ')
                FileList= dir([OutDir{17} 'CellShape*.mat']);
                
                FileList = fullfile({FileList.folder}.', {FileList.name}.');
                
                %run autoencoder with specified percent of training data
                [autoencoder, trainList]=CellShapeAutoencoder(FileList, 0.2);
                
                save([OutDir{17} 'autoencoder.mat'], 'autoencoder')
                writetable(trainList, [OutDir{17} 'encoded_cells.csv']);
                
            end
            
        end
    end
    fprintf('\n')
end

