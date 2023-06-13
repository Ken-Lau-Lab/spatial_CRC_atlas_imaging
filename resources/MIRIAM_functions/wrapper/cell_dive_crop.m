function[] = cell_dive_crop(base, slides)

% This function takes in a csv of slide names and crops AFRemoved and DAPI
% images into approximately 2048 x 2048 pixel tiles from AFRemoved and
% Registered Cell Dive tiffs
% base: base folder with imaging data
% slides: strings of slide names
% AFRemoved: 1 or 0 of whether to process AFRemoved Images or not
% DAPI: 1 or 0 of whether to process DAPI Images or not
% Created by: Eliot McKinley (etmckinley@gmail.com)


%this could be changed to parfor to run quicker
for s=1:length(slides)
    
    AFRemoved =1;
    DAPI = 1;
    folder_dir=[base slides{s} filesep];
    
    image_dir=join([ folder_dir 'AFRemoved'],"");
    save_dir= join([folder_dir 'AFRemoved Tiles' filesep],"");
        
    %change this if you want to use another image round for your DAPI
    dapi_round = "002";

    dapi_dir=join([folder_dir 'RegisteredImages' filesep 'S' dapi_round filesep], '');
    dapi_save_dir=join([ folder_dir 'DAPI Tiles' filesep],"");
    
    if ~exist(save_dir, 'dir')
        mkdir(save_dir);
    end
    
    if ~exist(dapi_save_dir, 'dir')
        mkdir(dapi_save_dir);
    end
       
    if AFRemoved==1
        images=imageDatastore(image_dir);
        images=images.Files;
        images = images(~contains(images, "downsample"));

        for image =1:length(images)
            fprintf(join([images{image} "\n"]))
            info=imfinfo(images{image});
            
            width=info(1).Width;
            height=info(1).Height;
            
            ncol=ceil(width/2048);
            nrow=ceil(height/2048);
            
            y_size=ceil(width/ncol);
            x_size=ceil(height/nrow);
                      
            short=images{image};
            delimiter=join([filesep 'AFRemoved' filesep],"");
            short=regexp(short, delimiter, 'split');
            short=short{2};
            short=strrep(short, '_pyr16', '');
            short=strrep(short, '.tif', '');
            AFList=dir(join([save_dir short '*'],""));
            
            if length(AFList) ~= ncol*nrow
                
                fprintf(join([images{image} ": reading\n"]))
                img=imread(images{image});

                for i=  1:nrow
                    for j=1:ncol
                        filename=strrep(images{image}, '.tif', ['_r_' num2str(i, '%03d') '_c_' num2str(j, '%03d') '.tif' ]);
                        filename=strrep(filename, join([filesep 'AFRemoved' filesep],""), join([filesep 'AFRemoved Tiles' filesep],""));
                        filename=strrep(filename, '_pyr16', '');
                        if ~isfile(filename)
                            
                            if(i<nrow && j<ncol)
                                img_out=img(x_size*(i-1)+1:x_size*i+floor(x_size*.1)+1, y_size*(j-1)+1:y_size*j+floor(y_size*.1)+1,:);
                            end
                            if(i==nrow && j<ncol)
                                img_out=img( height-x_size-x_size*.1:height, y_size*(j-1)+1:y_size*j+y_size*.1+1,:);
                            end
                            if(i<nrow && j==ncol)
                                img_out=img( x_size*(i-1)+1:x_size*i+x_size*.1+1, width-y_size-y_size*.1:width,:);
                            end
                            if(i==nrow && j==ncol)
                                img_out=img( height-x_size-x_size*.1:height, width-y_size-y_size*.1:width,:);
                            end
                            imwrite(img_out, filename, 'Compression', 'lzw');
                        end
                    end
                end
                
            end
        end
    end
    
    
    
    if DAPI==1
        images=imageDatastore(join([dapi_dir  '*_dapi_*'], ""));
        images=images.Files;
        
        for image =1:length(images)
            %fprintf(join([images{image} "\n"]))
            info=imfinfo(images{image});
            
            width=info(1).Width;
            height=info(1).Height;
            
            ncol=ceil(width/2048);
            nrow=ceil(height/2048);
            
            y_size=ceil(width/ncol);
            x_size=ceil(height/nrow);
            
            short=images{image};
            delimiter=join([filesep "S" dapi_round filesep], "");
            short=regexp(short, delimiter, 'split');
            short=short{2};
            short=strrep(short, '_pyr16', '');
            short=strrep(short, '.tif', '');
            AFList=dir(join([dapi_save_dir short '*'],""));
            
            if length(AFList) ~= ncol*nrow
                
                fprintf(join([images{image} ": reading\n"]))
                img=imread(images{image});
                
                for i=1:nrow
                    for j=1:ncol
                        filename=join([dapi_save_dir short '_r_' num2str(i, '%03d') '_c_' num2str(j, '%03d') '.tif' ],"");
                        if ~isfile(filename)
                            
                            if(i<nrow && j<ncol)
                                img_out=img(x_size*(i-1)+1:x_size*i+floor(x_size*.1)+1, y_size*(j-1)+1:y_size*j+floor(y_size*.1)+1,:);
                            end
                            if(i==nrow && j<ncol)
                                img_out=img( height-x_size-x_size*.1:height, y_size*(j-1)+1:y_size*j+y_size*.1+1,:);
                            end
                            if(i<nrow && j==ncol)
                                img_out=img( x_size*(i-1)+1:x_size*i+x_size*.1+1, width-y_size-y_size*.1:width,:);
                            end
                            if(i==nrow && j==ncol)
                                img_out=img( height-x_size-x_size*.1:height, width-y_size-y_size*.1:width,:);
                            end
                            imwrite(img_out, filename, 'Compression', 'lzw');
                        end
                    end
                end
                
            end
        end
    end

    if exist(join([folder_dir "Thresholded"], "" ), 'file')
        thresh_dir=join([folder_dir "Thresholded" filesep], "" );
        save_dir= join([folder_dir "Thresholded Tiles" filesep], "" );
        
        if ~exist(save_dir, 'dir')
            mkdir(save_dir);
        end

        images=imageDatastore(join([thresh_dir  "*threshold*"], ""));
        images=images.Files;
        for image =1:length(images)
            fprintf(join([images{image} "\n"]))
            %info=imfinfo(images{image});
            
            %get size of the original image
            
                      
            short=images{image};
            delimiter=join([filesep 'Thresholded' filesep],"");
            short=regexp(short, delimiter, 'split');
            short=short{2};
            short=strrep(short, '_threshold.png', '');
            short=strrep(short, '_region_', '_THRESHOLDED_region_');
            
            if contains(short, "DAPI")
                dapi_images=imageDatastore(join([dapi_dir  '*_dapi_*'], ""));
                dapi_images=dapi_images.Files;
                info = imfinfo(char(dapi_images(contains(dapi_images, strrep(short, "DAPI_THRESHOLDED_", "")))));
            else
                info = imfinfo(join([image_dir, filesep, strrep(short, "_THRESHOLDED_region_", "_AFRemoved_pyr16_region_"), ".tif"],""));
            end

            


            width=info(1).Width;
            height=info(1).Height;
            
            ncol=ceil(width/2048);
            nrow=ceil(height/2048);
            
            y_size=ceil(width/ncol);
            x_size=ceil(height/nrow);
            
            AFList=dir(join([save_dir short '*'],""));
            
            

            if length(AFList) ~= ncol*nrow
                
                fprintf(join([images{image} ": reading\n"]))
                img=imread(images{image});
                %resize threshold image to the original
                img = imresize(img, [height width]);
                for i=  1:nrow
                    for j=1:ncol
                        
                        filename = join([save_dir short '_r_' num2str(i, '%03d') '_c_' num2str(j, '%03d') '.png'],"");
                        
                        if ~isfile(filename)
                            
                            if(i<nrow && j<ncol)
                                img_out=img(x_size*(i-1)+1:x_size*i+floor(x_size*.1)+1, y_size*(j-1)+1:y_size*j+floor(y_size*.1)+1,:);
                            end
                            if(i==nrow && j<ncol)
                                img_out=img( height-x_size-x_size*.1:height, y_size*(j-1)+1:y_size*j+y_size*.1+1,:);
                            end
                            if(i<nrow && j==ncol)
                                img_out=img( x_size*(i-1)+1:x_size*i+x_size*.1+1, width-y_size-y_size*.1:width,:);
                            end
                            if(i==nrow && j==ncol)
                                img_out=img( height-x_size-x_size*.1:height, width-y_size-y_size*.1:width,:);
                            end
                            imwrite(img_out, filename);
                        end
                    end
                end
                
            end
        end

    end


end
