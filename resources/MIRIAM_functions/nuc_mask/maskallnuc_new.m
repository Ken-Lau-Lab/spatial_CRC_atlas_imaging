%output: mask = mask using all nuclear markers
%parameters: directory = nuclear dir, lev =50 for nuclear,
%totalmem=totalmembrane mask

function [mask] = maskallnuc_new (directory,lev,totalmem) 
filename=directory;
s=size(filename);
s=s(2);
ps=size(totalmem);
pixadj=1;
if ps(1)~=2048 || ps(2)~=2048 %check if images is from Cytell (>2048
    pixadj=3; %adjust for smaller pixel size if Cytell
end

total=zeros (ps(1),ps(2));
total=im2bw(total);
for i=1:s
    %adaptive thresholding for all files in dir

    
    fox = imread(char(filename(i)));
    fox = uint8( (double(fox) - double(min(fox(:)))) /(double(max(fox(:))) - double(min(fox(:)))) * 255 );
    fox=adapthisteq(fox);
    
    T=adaptthresh(fox,lev);
    I=imbinarize(fox,T);
    I=bwareaopen(I,ceil(50*pixadj));
    
    % do normal segmentation now on individual nuclear markers
    mask = bwlabel(I,4);
    
    %segment this nuclear mask
    output = segprocedure_large(mask,totalmem);
    
    % add to a total mask
    
    %total = comb2nuc(total,output,totalmem);
    total=im2bw(imadd (im2bw(total),im2bw(output)));
    %total=binaryadd (total, output);
end
mask=total;
end
