% This function is used to get the x, y, and z coordinates and the grain ID
% for each individual voxel from Grain Generator output .vtk files, 
% It can also outputs the cube edge length dimension (voxels).

% The function then transforms the coordinate data that denotes the edges
% of each voxel into nodal coordinate data for ease of downsampling. During
% this step the resolutions defined in DREAM 3D during structure generation
% are also computed and so are used as output for downsampling program.

function [cube_dimension,xcoord,ycoord,zcoord,grain_ID,xres_old,yres_old,zres_old] = vtkdatareader(DREAM3D_FILENAME)

DREAM3D_FILENAME = fopen(DREAM3D_FILENAME,'r'); %open data file

end_of_file = 0; %end of file identifier
ender = 0;
linecount = 0;
cnt = 0;

while end_of_file == 0 && ender == 0
        line1 = fgetl(DREAM3D_FILENAME);
        linecount = linecount+1;

        if linecount == 5
           A = sscanf(line1,'%10c %f %f %f');
           cube_dimension = A(11)-1; %Cube edge dimension
           coord_data_rows = round(cube_dimension/20); %Number of rows coordinate data will be stored in

%     % Define the lines where x-coordinates begin & end
xcoord_low_row = 7;
xcoord_high_row = 7 + coord_data_rows;

           for linecount = xcoord_low_row:xcoord_high_row
               line = fgetl(DREAM3D_FILENAME);
               i  = linecount - (xcoord_low_row - 1);
               X(i) = textscan(line,'%f','MultipleDelimsAsOne',1);
           end

           xcoord = cell2mat(X'); %use transpose of cell arrays

%     % Define the lines where y-coordinates begin & end
ycoord_low_row = xcoord_high_row + 1;
ycoord_high_row = ycoord_low_row + coord_data_rows;

           for linecount = ycoord_low_row:ycoord_high_row
               line = fgetl(DREAM3D_FILENAME);
               i = linecount - (ycoord_low_row - 1);
               Y(i) = textscan(line,'%f','MultipleDelimsAsOne',1);
           end

           ycoord = cell2mat(Y'); %use transpose of cell arrays

%     % Define the lines where z-coordinates begin & end
zcoord_low_row = ycoord_high_row + 1;
zcoord_high_row = zcoord_low_row + coord_data_rows;

           for linecount = zcoord_low_row:zcoord_high_row
               line = fgetl(DREAM3D_FILENAME);
               i = linecount - (zcoord_low_row - 1);
               Z(i) = textscan(line,'%f','MultipleDelimsAsOne',1);
           end

           zcoord = cell2mat(Z'); %use transpose of cell arrays

% Define the lines where grain ID data begins and ends
grain_low_row = zcoord_high_row + 3;
grain_high_row= grain_low_row + (cube_dimension)^3/20;

% Preallocate
grain_ID = zeros(cube_dimension^3,1);

                line = fgetl(DREAM3D_FILENAME);
                line = fgetl(DREAM3D_FILENAME);
                line = fgetl(DREAM3D_FILENAME);

           for linecount = 1:grain_high_row - grain_low_row
               if linecount < grain_high_row - grain_low_row
                    line = fgetl(DREAM3D_FILENAME);
                    grain_ID((linecount*20-19):(linecount*20),1) = (sscanf(line,'%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f '));

               else
                   line = fgetl(DREAM3D_FILENAME);
                   G = textscan(line,'%f');
                   grain_ID_end_row = cell2mat(G'); %use transpose of cell arrays
                   [last_row_elements,~] = size(grain_ID_end_row);
                   grain_ID((linecount*20-19):(linecount*20-19+last_row_elements-1),1) = grain_ID_end_row;
                   ender = 1;
                   
               end
               
           end

        end

        end_of_file = feof(DREAM3D_FILENAME);
end

fclose(DREAM3D_FILENAME);

%% Transform edge coordinates into nodal coordinates

% x,y,z old resolutions
xres_old = abs(xcoord(2) - xcoord(1));
yres_old = abs(ycoord(2) - ycoord(1));
zres_old = abs(zcoord(2) - zcoord(1));

end
