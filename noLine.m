function nol = noLine(filePath)

fid = fopen(filePath);
fseek(fid, 0, 'eof');
chunksize = ftell(fid);
fseek(fid, 0, 'bof');
ch = fread(fid, chunksize, '*uchar');
nol = sum(ch == sprintf('\n')); % number of lines 
fclose(fid);