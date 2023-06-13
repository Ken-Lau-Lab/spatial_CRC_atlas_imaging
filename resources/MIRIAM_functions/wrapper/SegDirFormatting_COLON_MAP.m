function [AFRemoved, DAPI, OutDir]=SegDirFormatting_TMA(SlideDir)
%SegDirFormatting Create and format data structure of single cell
%segmentation
%Inputs:
%SlideDir= directory for slide containing AFRemoved images folder and
%Registered images folder - assumes Round 001 is baseline and all fies are
%tif format
%
%Outputs
%AFRemoved= string for location of AFRemoved images
%DAPI=string for location of DAPI images
%OutDir= cell array of strings for the output files



AFRemoved=join([SlideDir filesep 'AFRemoved Tiles' filesep],"");

%DAPI=[SlideDir '/RegisteredImages/S002/'];
%DAPI=[SlideDir '/RegisteredImages/S003/'];
DAPI=join([SlideDir filesep 'DAPI Tiles' filesep],"");



% check if various output directories exist and create them if not
if ~exist(join([SlideDir filesep 'SegQuant'],""), 'dir')
    mkdir(join([SlideDir filesep 'SegQuant'],""));
end

if ~exist(join([SlideDir filesep 'SegQuant' filesep 'Stacks'], ""), 'dir')
    mkdir(join([SlideDir filesep 'SegQuant' filesep 'Stacks'], ""));
end


if ~exist(join([SlideDir filesep 'SegQuant' filesep 'CellSeg'], ""), 'dir')
    mkdir(join([SlideDir filesep 'SegQuant' filesep 'CellSeg'], ""));
end

if ~exist(join([SlideDir filesep 'SegQuant' filesep 'CellSegFinal'], ""), 'dir')
    mkdir(join([SlideDir filesep 'SegQuant' filesep 'CellSegFinal'], ""));
end

if ~exist(join([SlideDir filesep 'SegQuant' filesep 'EpiMask'], ""), 'dir')
    mkdir(join([SlideDir filesep 'SegQuant' filesep  'EpiMask'], ""));

end

if ~exist(join([SlideDir filesep 'SegQuant' filesep 'Novlp'], ""), 'dir')
    mkdir(join([SlideDir filesep 'SegQuant' filesep 'Novlp'], ""));

end
if ~exist(join([SlideDir filesep 'SegQuant' filesep 'NucMask'], ""), 'dir')
    mkdir(join([SlideDir filesep 'SegQuant' filesep 'NucMask'], ""));

end

if ~exist(join([SlideDir filesep 'SegQuant' filesep 'SuperMem'], ""), 'dir')
    mkdir(join([SlideDir filesep 'SegQuant' filesep 'SuperMem'], ""));

end

if ~exist(join([SlideDir filesep 'SegQuant' filesep 'MemMask'], ""), 'dir')
    mkdir(join([SlideDir filesep 'SegQuant' filesep 'MemMask'], ""));

end
if ~exist(join([SlideDir filesep 'SegQuant' filesep 'NucMaskFinal'], ""), 'dir')
    mkdir(join([SlideDir filesep 'SegQuant' filesep 'NucMaskFinal'], ""));

end
if ~exist(join([SlideDir filesep 'SegQuant' filesep 'PosStats'], ""), 'dir')
    mkdir(join([SlideDir filesep 'SegQuant' filesep 'PosStats'], ""));

end

if ~exist(join([SlideDir filesep 'SegQuant' filesep 'ML'], ""), 'dir')
    mkdir(join([SlideDir filesep 'SegQuant' filesep 'ML'], ""));

end

if ~exist(join([SlideDir filesep 'SegQuant' filesep 'TumorMask'], ""), 'dir')
    mkdir(join([SlideDir filesep 'SegQuant' filesep 'TumorMask'], ""));
end

if ~exist(join([SlideDir filesep 'SegQuant' filesep 'CellShape'], ""), 'dir')
    mkdir(join([SlideDir filesep 'SegQuant' filesep 'CellShape'], ""));
end


%load the OutDir cell array
OutDir{1}= join([SlideDir filesep 'SegQuant' filesep 'CellSeg' filesep], "");
OutDir{2}= join([SlideDir filesep 'SegQuant' filesep 'CellSegMask' filesep], "");
OutDir{3}= join([SlideDir filesep 'SegQuant' filesep 'EpiMask' filesep], "");
OutDir{4}= join([SlideDir filesep 'SegQuant' filesep 'Novlp' filesep], "");
OutDir{5}= join([SlideDir filesep 'SegQuant' filesep 'NucMask' filesep], "");
OutDir{6}= join([SlideDir filesep 'SegQuant' filesep 'NucSeg' filesep], "");
OutDir{7}= join([SlideDir filesep 'SegQuant' filesep 'SuperMem' filesep], "");
OutDir{8}= join([SlideDir filesep 'SegQuant' filesep], "");
OutDir{9}= join([SlideDir filesep 'SegQuant' filesep 'AFMask' filesep], "");
OutDir{10}= join([SlideDir filesep 'SegQuant' filesep 'MemMask' filesep], "");
OutDir{11}= join([SlideDir filesep 'SegQuant' filesep 'NucMaskFinal' filesep], "");
OutDir{12}= join([SlideDir filesep 'SegQuant' filesep 'CellSegFinal' filesep], "");
OutDir{13}= join([SlideDir filesep 'SegQuant' filesep 'PosStats' filesep], "");
OutDir{14}= join([SlideDir filesep 'SegQuant' filesep 'ML' filesep], "");
OutDir{15}= join([SlideDir filesep 'SegQuant' filesep 'Stacks' filesep], "");
OutDir{16}= join([SlideDir filesep 'SegQuant' filesep 'TumorMask' filesep], "");
OutDir{17}= join([SlideDir filesep 'SegQuant' filesep 'CellShape' filesep], "");
OutDir{18}= join([SlideDir filesep 'SegQuant' filesep 'Annotated' filesep], "");



