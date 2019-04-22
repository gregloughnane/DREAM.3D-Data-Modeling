function [] = PipelineCreator(PipelineIn,PipelineOut,D3DFileName,CSVFileName,newres,newpnoise,newbnoise)
% 
% PipelineIn is the name of the input .json file, with directory
% PipelineOut is the output pipeline .json file name desired, without directory
% D3DFileName is the output .dream3d file name desired (to be entered in the pipeline file), without directory
% newres is the new value for x, y, z resolution (for isotropic down-sampling)
% newpnoise is the new value for poisson (random) noise
% newbnoise is the new value for boundary noise

fid = fopen(PipelineIn,'r'); %open data file
newfile = fopen(PipelineOut,'w'); %Need to name new file here

end_of_file = 0; %end of file identifier
ender = 0;
linecount = 0;

while end_of_file == 0 && ender == 0
lines = fgetl(fid);
linecount = linecount+1; % Initiate counter
  if    linecount < 453
        ScanPrintLine(lines,newfile);
%% ChangeResolution
    elseif linecount == 453 %%% XRES
        xres = sscanf(lines,'%c');
        nxres = strcat(xres(1,1:17),newres); fprintf(newfile,'%s,\n',nxres);
    elseif linecount == 454 %%% YRES
        yres = sscanf(lines,'%c');
        nyres = strcat(yres(1,1:17),newres); fprintf(newfile,'%s,\n',nyres);
    elseif linecount == 455 %%% ZRES
        zres = sscanf(lines,'%c');
        nzres = strcat(zres(1,1:17),newres); fprintf(newfile,'%s \n',nzres);
    elseif linecount > 455 && linecount < 485
		ScanPrintLine(lines,newfile);        
%% Write DREAM.3D Data File
    elseif linecount == 485
		%%% Output
        d3d = sscanf(lines,'%c')
        nd3d = strcat(d3d(1,1:96),D3DFileName);
        fprintf(newfile,'%s",\n',nd3d);
	elseif linecount > 485 && linecount < 494
		ScanPrintLine(lines,newfile);        
%% Write Feature Data CSV File
    elseif linecount == 494
        %%% Output
        csv = sscanf(lines,'%c');
        ncsv = strcat(csv(1,1:101),CSVFileName);
        fprintf(newfile,'%s",\n',ncsv);
	elseif linecount > 494 && linecount < 533
		ScanPrintLine(lines,newfile);    
%% AddBadData
    elseif linecount == 533 %%% Boundary Noise
        bn = sscanf(lines,'%c'); 
        nbn = strcat(bn(1,1:31),newbnoise); fprintf(newfile,'%s,\n',nbn);
	elseif linecount > 533 && linecount < 543
		ScanPrintLine(lines,newfile);            
	elseif linecount == 543 %%% Poisson Noise
		pn = sscanf(lines,'%c'); 	
		npn = strcat(pn(1,1:30),newpnoise); fprintf(newfile,'%s \n',npn);
    elseif linecount > 543 && linecount < 647	
		ScanPrintLine(lines,newfile);   
		
    end
    
    end_of_file = feof(fid);
    
end
fclose(fid);
fclose(newfile);
end