% Greg Loughnane
% September 23, 2011

%%%% This code works with DREAM.3D version DREAM3D-2011.09.21-Win64 
%%%% (This is written out to 2nd line of output VTK file)

%% Description
% This program down-samples the VTK output file from a microstructure 
% created using DREAM.3D's grain generator. It then writes another VTK file
% back out with the modified grain ID data. The down-sampling origin is
% determined using random translations in the x, y, and z directions.

%% Initialize
clc
clear all
close all
format loose

%% Get input data from Grain Generator .vtk output file

tic % Clock data reading in speed

[filename,pathname] = uigetfile('*.vtk','Choose .vtk File','MultiSelect','on');
file_name = strcat(pathname,filename);
[cube_dimension,xcoord,ycoord,zcoord,grain_ID,xres_old,yres_old,zres_old]=vtkreader(file_name);

toc % Clock data reading in speed

%% Scaling factor for resolution
%vox_input = input('Enter the voxel size scalars to use for downsampling in vector form: ');
%vox_input = [1.335 2.003 4.006 8.011 13.352 40.056]; % Light tail 1/3/5/10/20
%vox_input = [1.502 3.004 6.008 10.014 30.042]; % Heavy tail 1/3/5/10/20/30
vox_input = [10.014 30.042];
depth_input = [1:10];

for vox_scalar = vox_input
    for depth_multiplier = depth_input
%% ***** Define parameter values from input data *****`

%%% Original Structure    
    % Number of points along each dimension
    xpoints = cube_dimension;
    ypoints = cube_dimension;
    zpoints = cube_dimension;
    
    % x,y,z edge lengths (voxels * microns)
    xbox_size = cube_dimension * xres_old;
    ybox_size = cube_dimension * yres_old;
    zbox_size = cube_dimension * zres_old;

%%% Down-Sampled Structure
    % x,y,z resolutions (microns)
    xres_new = xres_old * vox_scalar;
    yres_new = yres_old * vox_scalar;
    zres_new = depth_multiplier * zres_old * vox_scalar;
    
    % Number of points along each dimension
    ds_xpoints = floor(xbox_size / xres_new);
    ds_ypoints = floor(ybox_size / yres_new);
    ds_zpoints = floor(zbox_size / zres_new);


%% ***** Determine new voxel numbers for down-sampled volume *****

    % Preallocate memory with cell arrays
    voxel_new = double(size(ds_xpoints * ds_ypoints * ds_zpoints));

    % Create random numbers for origin translations and output file naming
    % scheme
    rand_no = rand(); % Make this a variable for later on to be written out
    rx = 2*(rand_no - 0.5)*0.499*xres_old;
    ry = 2*(rand() - 0.5) *0.499*yres_old;
    rz = 2*(rand() - 0.5) *0.499*zres_old;

    for i = 1:ds_zpoints
    
        for j = 1:ds_ypoints
        
            for k = 1:ds_xpoints
                
                % Define x,y,z coordinates of the original structure that the
                % down-sampled point is in
                x = (k-1)*xres_new + rx;
                y = (j-1)*yres_new + ry;
                z = (i-1)*zres_new + rz;

                % Compute the column, row and plane of the original structure
                % the down-sampled point is in
                column = floor((x + 0.5 * xres_old) / xres_old);
                row = floor((y + 0.5 * yres_old) / yres_old);
                plane = floor((z + 0.5 * zres_old) /zres_old);

                % Determine down-sampled voxel number
                voxel_new = ((i-1) * ds_xpoints*ds_ypoints) + ((j-1)*ds_xpoints) + k;
                      
                % Determine corresponding original voxel number
                voxel_orig = floor(((plane)*xpoints*ypoints) + ((row)*xpoints) + column);
                        
                % Set grain ID from voxel_orig to the down-sampled points
                ds_grainIDs(voxel_new) = grain_ID(voxel_orig + 1);
         
            end
            
        end
        
    end
    
    % Transpose result
    ds_grainIDs = ds_grainIDs';

%% ***** Define variables for writing back out to VTK *****
   
    tic % Clock file writing out speed
    
    % x,y,z coordinates for down-sampled volume
    xcoord_new = vox_scalar .* xcoord;
    ycoord_new = vox_scalar .* ycoord;
    zcoord_new = vox_scalar .* zcoord;

    % Voxel dimensions of down-sampled cube
    xdim = floor((cube_dimension)/vox_scalar) + 1;
    ydim = floor((cube_dimension)/vox_scalar) + 1;
    zdim = floor((cube_dimension)/vox_scalar) + 1;

    % Number of voxels in down-sampled volume
    cell_data = (xdim - 1)*(ydim - 1)*(zdim - 1);

%% ************ Write data to VTK ************
        
    % Define output file name format
    folder = (exp(1.1)/xres_old)/vox_scalar;
    num = num2str(folder,'%5.0f');
    random = num2str(rand_no*1000,'%5.0f');
    seperator = ('_');
    depth = num2str(depth_multiplier,'%i');
    ds_filename = strcat(num,seperator,depth,seperator,random,seperator,filename);
    
%%%% ********** Write Header **********
    
    fid = fopen(ds_filename,'w');
    fprintf(fid, '# vtk DataFile Version 2.0');
    fprintf(fid, '\n');
    
    fprintf(fid, 'data set from DREAM3D Version');
    fprintf(fid, '\n');
    
    fprintf(fid, 'ASCII');
    fprintf(fid, '\n');
    
    fprintf(fid, 'DATASET RECTILINEAR_GRID');
    fprintf(fid, '\n');
    
    fprintf(fid, 'DIMENSIONS %i %i %i', xdim, ydim, zdim);
    fprintf(fid, '\n');

%%%% ********** Write X-coordinate data **********
    
    fprintf(fid, 'X_COORDINATES %i float', xdim);
    fprintf(fid, '\n');
    xrows = floor((xdim-21)/20);
    
    % # of rows written to VTK cannot be negative
    if xrows < 0
        xrows = 0;
    end
    
    % CASE I: Voxel dimensions do not fill up first written row in VTK file
    if xdim <= 21
        for a = 1:xdim
            fprintf(fid, '%3.6f ', xcoord_new(a));
        end
    end
    
    % CASE II: Voxel dimensions fill up first written row in VTK file
    if xdim > 21
        
        % Write first row
        for a = 1:21
            fprintf(fid, '%3.6f ', xcoord_new(a));
        end
        fprintf(fid, '\n');
        a = a + 1;
        
        % Write rows in between first and last
        while a >= 22
            
            % Case II A: Voxel dimensions do not fill up second written row in
            % VTK file
            if xdim <= 42
                fprintf(fid, '%3.6f ', xcoord_new(a:xdim));
                break
                break
            end
            
            % Case II B: Voxel dimensions do fill up second written row in
            % VTK file
            if xdim >42
                fprintf(fid, '%3.6f ', xcoord_new(a:a + 19));
                fprintf(fid, '\n');
                a = a + 20;
            end
            
            % For Case II B: This requires that if there are not enough
            % elements left to write an entire row of 20 elements, the program
            % moves to the next loop to write the last row
            if a > (xdim - 20)
                break
                break
            end
        end
    end
    
    % The number of complete rows written for x-coordinates cannot be 0
    if xrows ~= 0
        
        % Write last row, which typically has < 20 elements
        for b = (xrows * 20 + 22):xdim
            fprintf(fid, '%3.6f ', xcoord_new(b));
        end
    end
    fprintf(fid, '\n');
    
%%%% ********** Write Y-coordinates **********
    
    fprintf(fid, 'Y_COORDINATES %i float', ydim);
    fprintf(fid, '\n');
    yrows = floor((ydim-21)/20);
    
    % # of rows written to VTK cannot be negative
    if yrows < 0
        yrows = 0;
    end
    
    % CASE I: Voxel dimensions do not fill up first written row in VTK file
    if ydim <= 21
        for a = 1:ydim
            fprintf(fid, '%3.6f ', ycoord_new(a));
        end
    end
    
    % CASE II: Voxel dimensions fill up first written row in VTK file
    if ydim > 21
        
        % Write first row
        for a = 1:21
            fprintf(fid, '%3.6f ', ycoord_new(a));
        end
        fprintf(fid, '\n');
        a = a + 1;
        
        % Write rows in between first and last
        while a >= 22
            
            % Case II A: Voxel dimensions do not fill up second written row in
            % VTK file
            if ydim <= 42
                fprintf(fid, '%3.6f ', ycoord_new(a:ydim));
                break
                break
            end
            
            % Case II B: Voxel dimensions do fill up second written row in
            % VTK file
            if ydim >42
                fprintf(fid, '%3.6f ', ycoord_new(a:a + 19));
                fprintf(fid, '\n');
                a = a + 20;
            end
            
            % For Case II B: This requires that if there are not enough
            % elements left to write an entire row of 20 elements, the program
            % moves to the next loop to write the last row
            if a > (ydim - 20)
                break
                break
            end
        end
    end
    
    % The number of complete rows written for x-coordinates cannot be 0
    if yrows ~= 0
        
        % Write last row, which typically has < 20 elements
        for b = (yrows * 20 + 22):ydim
            fprintf(fid, '%3.6f ', ycoord_new(b));
        end
    end
    fprintf(fid, '\n');
    
%%%% ********** Write Z-coordinates **********
    
    fprintf(fid, 'Z_COORDINATES %i float', zdim);
    fprintf(fid, '\n');
    zrows = floor((zdim-21)/20);
    
    % # of rows written to VTK cannot be negative
    if zrows < 0
        zrows = 0;
    end
    
    % CASE I: Voxel dimensions do not fill up first written row in VTK file
    if zdim <= 21
        for a = 1:zdim
            fprintf(fid, '%3.6f ', zcoord_new(a));
        end
    end
    
    % CASE II: Voxel dimensions fill up first written row in VTK file
    if zdim > 21
        
        % Write first row
        for a = 1:21
            fprintf(fid, '%3.6f ', zcoord_new(a));
        end
        fprintf(fid, '\n');
        a = a + 1;
        
        % Write rows in between first and last
        while a >= 22
            
            % Case II A: Voxel dimensions do not fill up second written row in
            % VTK file
            if zdim <= 42
                fprintf(fid, '%3.6f ', zcoord_new(a:xdim));
                break
                break
            end
            
            % Case II B: Voxel dimensions do fill up second written row in
            % VTK file
            if zdim >42
                fprintf(fid, '%3.6f ', zcoord_new(a:a + 19));
                fprintf(fid, '\n');
                a = a + 20;
            end
            
            % For Case II B: This requires that if there are not enough
            % elements left to write an entire row of 20 elements, the program
            % moves to the next loop to write the last row
            if a > (zdim - 20)
                break
                break
            end
        end
    end
    
    % The number of complete rows written for x-coordinates cannot be 0
    if zrows ~= 0
        
        % Write last row, which typically has < 20 elements
        for b = (zrows * 20 + 22):zdim
            fprintf(fid, '%3.6f ', zcoord_new(b));
        end
    end
    fprintf(fid, '\n');
    
%%%% ********** Write "Cell data" **********
    fprintf(fid, 'CELL_DATA %i', cell_data);
    fprintf(fid, '\n');
    
%%%% ********** Write "Scalars" **********    
    fprintf(fid, 'SCALARS GrainID int 1');
    fprintf(fid, '\n');
    
%%%% ********** Write "Grain ID Lookup Table" **********  
    fprintf(fid, 'LOOKUP_TABLE default');
    fprintf(fid, '\n');
    
%%%% ********** Write Grain IDs **********
    no_grainIDs = length(ds_grainIDs);
    no_rows = floor(no_grainIDs / 20);
    
    a = 1;
    while a <= no_rows * 20
        fprintf(fid, '%i ', ds_grainIDs(a:a + 19));
        fprintf(fid, '\n');
        a = a+20;
    end
    
    for b = no_rows * 20:no_grainIDs
        fprintf(fid, '%i ', ds_grainIDs(b));
    end
    
    fprintf(fid, '\n');
    
%%%% ********** Write "Scalars Surface Voxel **********
    fprintf(fid, 'SCALARS SurfaceVoxel int 1');
    
    toc % Clock file writing out speed
    
    end
end
