function [Stats]=MxIF_quantify_stroma2_thresh(i, WatCellSeg, AFRemoved, AFList, PosList, ThreshList, DAPI, tumorMask, OutPos)
% Funtion to quantify AFRemoved images
%inputs:
%WatCellSeg= cell segmenation
%AFRemoved

%set up cells
        Slide=split(AFRemoved, '/scan_alpha/');
            Slide=Slide{2};
            Slide=split(Slide, '/AFRemoved');
            Slide=Slide{1};
            
        Stats=struct2table(regionprops(WatCellSeg, { 'Centroid' ,'Area'})); %get morphometrics for cell
        C=array2table(Stats.Centroid);
        if size(C,1) <10
            fprintf('Skipping Quantification. No Cell Seg File\n')
            return
        end
        Stats=[C(:,1) C(:,2) Stats(:,1)];
        clear C
        Stats.Properties.VariableNames={'Cell_Centroid_X' 'Cell_Centroid_Y' 'Cell_Area'};
        ID=struct2table(regionprops(WatCellSeg, WatCellSeg, { 'PixelValues'}));       
        ID=table(cellfun(@nanmedian, ID{:,1}));
        ID.Properties.VariableNames={'ID'}; %rename table variable
        Position = repmat({OutPos{i}},height(ID),1); 
        Position=cell2table(Position);
        %Position=array2table(ones(height(ID),1)*str2num(strrep(OutPos{i},'pyr16_spot_','' )));
        Position.Properties.VariableNames={'Pos'}; %rename table variable
        Stats=[ID Position Stats]; %add ID to Stats table
        Stats=sortrows(Stats,1);
        
        AFQuantDAPI= struct2table(regionprops(WatCellSeg, DAPI, { 'PixelValues' }));
        AFQuantDAPI=table(cellfun(@nanmedian, AFQuantDAPI{:,1}));
        AFQuantDAPI.Properties.VariableNames={['Median_Cell_DAPI']};
        AFQuantDAPI=[ID AFQuantDAPI];
        AFQuantDAPI=sortrows(AFQuantDAPI,1);
        AFQuantDAPI=AFQuantDAPI(:,2);
        
        Stats=[Stats AFQuantDAPI];
        %shapes=shape_determination(WatCellSeg);
        
        %s=size(WatCellSeg);
        %NoOvlp=zeros(s(1), s(2), 3);
        
        for j = 1:length(AFList) %quantify each marker
            fprintf([AFList{j} ' '])
            AFim=imread([AFRemoved AFList{j} '_AFRemoved_' PosList{i} '.tif']); %read biomarker image
            %AForig=AFim;
            %Quantify whole cell stats
            %ID=struct2table(regionpropsETM(WatCellSeg, WatCellSeg, { 'MedianIntensity'}));
            AFQuantCell= struct2table(regionprops(WatCellSeg, AFim, { 'PixelValues' }));
            AFQuantCell=table(cellfun(@nanmedian, AFQuantCell{:,1}));
            AFQuantCell.Properties.VariableNames={['Median_Cell_' AFList{j}]};
            AFQuantCell=[ID AFQuantCell];
            AFQuantCell=sortrows(AFQuantCell,1);
            AFQuantCell=AFQuantCell(:,2);
            
            Stats=[Stats AFQuantCell];
            %Stats=[Stats AFQuant(:,2)];
            
            
            
            %quantify nuclear stats
            %AFim=double(AForig); %get only nuclear signal
            %AFim(mask==0)=nan;
            %CellID=WatCellSeg.*uint16(nucmask); %get regions only in nuc
            %imwrite(AFim, [OutDir 'AFNuc_' AFList{j} '_' PosList{i} '.tif'] ) %write 16 bit tiff
            %AFim(AFim==0)=nan;
            %ID= struct2table(regionpropsETM(WatCellSeg, WatCellSeg, { 'MedianIntensity'}));
            %if j==1 %for first marker
                
            %    Area=struct2table(regionprops(WatCellSeg, mask, { 'PixelValues'}));
            %    Area=table(cellfun(@nansum, Area{:,1}));
             %   Area.Properties.VariableNames={'Nuc_Area'};
              %  Area=[ID Area];
             %   Area=sortrows(Area,1);
            %    Stats=[Stats Area(:,2)];
            %end
            
            %AFQuantNuc= struct2table(regionprops(WatCellSeg, AFim, { 'PixelValues' }));
            %AFQuantNuc=table(cellfun(@nanmedian, AFQuantNuc{:,1}));
            %AFQuantNuc.Properties.VariableNames={['Median_Nuc_' AFList{j}]};
            %AFQuantNuc=[ID AFQuantNuc];
            %AFQuantNuc=sortrows(AFQuantNuc,1);
            %AFQuantNuc=AFQuantNuc(:,2);
            %Stats=[Stats AFQuant(:,2)];
            
            
            %quantify cell edge (mem) stats
%             AFim=double(AForig); %get only non-nuclear signal
%             MemMask=WatCellSeg==0;
%             
%             disksize=5;
%             %if pixadj>1
%             %    disksize=ceil(5*pixadj/1.5);
%             %end
%             
%             MemMask=imdilate(MemMask,strel('square',disksize));
%             MemMask(WatCellSeg==0)=0;
%             MemMask(mask==1)=0;
%             AFim(MemMask==0)=nan;
%             %imwrite(AFim, [OutDir 'AFNuc_' AFList{j} '_' PosList{i} '.tif'] ) %write 16 bit tiff
%             %AFim(AFim==0)=nan;
%             %ID= struct2table(regionpropsETM(WatCellSeg, WatCellSeg, { 'MedianIntensity'}));
%             if j==1 %for first marker
%                 
%                 Area=struct2table(regionprops(WatCellSeg, MemMask, { 'PixelValues'}));
%                 Area=table(cellfun(@nansum, Area{:,1}));
%                 Area.Properties.VariableNames={'Mem_Area'};
%                 Area=[ID Area];
%                 Area=sortrows(Area,1);
%                 Stats=[Stats Area(:,2)];
%             end
%             
%             AFQuantMem= struct2table(regionprops(WatCellSeg, AFim, { 'PixelValues' }));
%             AFQuantMem=table(cellfun(@nanmedian, AFQuantMem{:,1}));
%             AFQuantMem.Properties.VariableNames={['Median_Mem_' AFList{j}]};
%             AFQuantMem=[ID AFQuantMem];
%             AFQuantMem=sortrows(AFQuantMem,1);
%             AFQuantMem=AFQuantMem(:,2);
%             %Stats=[Stats AFQuant(:,2)];
            
            
            
            %quantify non nuclear and non mem (cyt) stats
            %AFim=double(AForig); %get only non-nuclear signal
            %CytMask=WatCellSeg>0 & mask==0;
            %AFim(CytMask==0)=nan;
            %imwrite(AFim, [OutDir 'AFNuc_' AFList{j} '_' PosList{i} '.tif'] ) %write 16 bit tiff
            %AFim(AFim==0)=nan;
            %ID= struct2table(regionpropsETM(WatCellSeg, WatCellSeg, { 'MedianIntensity'}));
            %if j==1 %for first marker
                
            %    Area=struct2table(regionprops(WatCellSeg, CytMask, { 'PixelValues'}));
            %    Area=table(cellfun(@nansum, Area{:,1}));
            %    Area.Properties.VariableNames={'Edge_Area'};
            %    Area=[ID Area];
            %    Area=sortrows(Area,1);
            %    Stats=[Stats Area(:,2)];
            %end
            
            %AFQuantCyt= struct2table(regionprops(WatCellSeg, AFim, { 'PixelValues' }));
            %AFQuantCyt=table(cellfun(@nanmedian, AFQuantCyt{:,1}));
            %AFQuantCyt.Properties.VariableNames={['Median_Edge_' AFList{j}]};
            %AFQuantCyt=[ID AFQuantCyt];
            %AFQuantCyt=sortrows(AFQuantCyt,1);
            %AFQuantCyt=AFQuantCyt(:,2);
            
            %Stats=[Stats AFQuantCell AFQuantNuc AFQuantCyt];
            
            
            
            
        end
        %get pixels on edges of watcellseg
        if ~isempty(tumorMask)
                tumQuantCell= struct2table(regionprops(WatCellSeg, tumorMask, { 'PixelValues' }));
                tumQuantCell=table(cellfun(@nanmedian, tumQuantCell{:,1}));
                tumQuantCell.Properties.VariableNames={['Tumor']};
                tumQuantCell=[ID tumQuantCell];
                tumQuantCell=sortrows(tumQuantCell,1);
                tumQuantCell=tumQuantCell(:,2);
                Stats=[Stats tumQuantCell];
        end
        
         if~isempty(ThreshList)
            
            for t=1:length(ThreshList)
            %thresh_image=imread(join(["/Users/etmckinley/Dropbox (VUMC)/T cell quantification/Images/" Slide "/Tiles/" ThreshList{t} "_threshold_" PosList{i} ".png"],""));
            thresh_image=imread(join(["/Users/etmckinley/Dropbox (VUMC)/COLON MAP/Adenomas/THRESHOLDED/" Slide "/Tiles/" ThreshList{t} "_threshold_" PosList{i} ".png"],""));
            ThreshQuant=struct2table(regionprops(WatCellSeg, thresh_image, {'PixelValues'}));
            ThreshQuant=table(cellfun(@nansum, ThreshQuant{:,1}));
            ThreshQuant.Properties.VariableNames={['Area_Thresh_' ThreshList{t}]};
            ThreshQuant=[ID ThreshQuant];
            ThreshQuant=sortrows(ThreshQuant,1);
            ThreshQuant=ThreshQuant(:,2);
            Stats=[Stats ThreshQuant];
            end
        end
        
        
        
    end