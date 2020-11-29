
[File,Path]=uigetfile('*.txt','MultiSelect','off');

 str=sprintf('%s', [Path File]);                            %makes str be the name of file (along with its path)
 file=sprintf('%s',File);                          
 format shortG                                              %change formating so doesn't show 0's for e-11 values. 

 %Need to skip 1st row, b/c it has some other info in it so use import data
 %fnc.
 %data = importdata(str, ' ', 0);              %importdata(FILENAME, DELIM, NHEADERLINES) loads data into a struct from 
                                                %ASCII file FILENAME, reading numeric data starting from line NHEADERLINES+1.
 
 data = importdata(str)
                                                %struct has fields data and textdata
 position = data(:,1);
 fullV(2:73) = data(:,2);
 p_full = data(:,3)/N;  %convert back to scaled form b/c need to calculate nT's!
 n_full = data(:,4)/N;
 
 
pT_full = H_D*(p_full*N/N_HOMO).^(T/T_tD);   %USE .^ for element wise (scalar) power
nT_full = H_A*(n_full*N/N_LUMO).^(T/T_tA);

pT_full = pT_full/N;
nT_full = nT_full/N;


%NOTE: Vbi is a constant for the given system--> depends on energy levels
%of materials
 fullV(1)          = 0.0;
 fullV(num_cell+1) = (Vbi -0.94)/Vt;  %for intial guess: Va = 0, b/c is equil. run, so Vbi-Va --> Vbi
 %0.94, b/c want  solution when 1st itertied to 0.94V


