% This script initializes and runs the MIRIAM segmentation pipeline
% created by Eliot McKinley
% etmckinley@gmail.com

%% USER INPUTS START %%

%file path to csv of slide names with "slide_id" as the column name
%slide_path = "/Users/etmckinley/Dropbox (VUMC)/Cell Dive Documentation/Example data/slide_examples.csv";
%slide_path = "/Users/etmckinley/Dropbox (VUMC)/COLON MAP/Cody Slides/slide threshold list.csv";
%slide_path = "/Users/etmckinley/Dropbox (VUMC)/COLON MAP/Cody Slides/slide threshold done 2022-12-14.csv";
slide_path = "/Volumes/CellDive/scan_alpha/st_slides_10Jan23.csv";

%file path to location of imaging data (make sure to end with the file
%seperator ("/" of "\")
base = '/Volumes/CellDive/scan_alpha/';
%base ='/Users/etmckinley/Dropbox (VUMC)/Cell Dive Documentation/Example data/';

%% Inputs for Crop function
% 1 or 0 if you want to run the cropping
crop = 0;


%% Inputs for Segmentation Function
% 1 or 0 if you want to run segmentation
segmentation = 1;

% what position you want to start with (typically 1)
start = 1;

% 1 or 0 if you want to quantify your data or just run the segmentation
quantify=1;

% 1 or 0 if you want to segment the stroma
stroma=1;

% 1 or 0 if you have tumor masks in the SegQuant Directory and want to
% use them as a mask
tumor=0;

%location of the full path to the ilastik file from the epithelial/stromal pixel classification
ilastik_epi_file = "../resources/epi.ilp";

%location of the full path to the ilastik file from the membrane/nucleus/cytoplasm pixel classification
ilastik_mem_file ="../resources/mem.ilp";

%location if Ilastik application
ilastik_app_loc= "/Applications/ilastik-1.4.0rc8-OSX.app/Contents/ilastik-release/run_ilastik.sh";

% a string of markers used to create the image stacks, this can be a CSV with DAPI as the first row or
% input directly
seg_markers=readtable('../resources/segmentation_markers.csv', 'Delimiter', ',');
seg_markers=seg_markers.Marker;
seg_markers=seg_markers(2:length(seg_markers));
%seg_markers = {"DAPI", "NAKATPASE"};


% 1 or 0 if you want to use the deep learning shape alorithm (typically
% you do not)
shape=0;


%% USER INPUTS END %%

warning('off','all')


%read in slide name csv
slides=readtable(slide_path, 'Delimiter', ',');
slides=slides.slide_id;
%slides=slides(12:end);


% if you need to filter by a prefix you can use this:
%slides=slides(startsWith(slides, 'HTA'));

%manual input of slides
%slides = {"WD84216"};

%Crop the data
if crop
    cell_dive_crop(base, slides)
end



if segmentation
    for  s= 1:length(slides)
        SlideDir=[base slides{s}];
        CellSegQuant_CellDive_mesmer(SlideDir, quantify, shape, stroma, tumor, start, ilastik_epi_file, ilastik_mem_file, ilastik_app_loc, seg_markers)
    end
end
