function[] = ScanPrintLine(lines,newfile)
 
    A = sscanf(lines,'%c'); 
    fprintf(newfile,'%s \n', A);