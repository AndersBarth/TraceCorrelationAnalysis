function save_to_fcsfit(Name1,Name2,Path,FileName,Cor_Times,Cor_Average,Cor_SEM,Cor_Array)
%%% Saves data
Header = ['Correlation file for: ' strrep(fullfile(Path, FileName),'\','\\') ' of Channels ' Name1 ' cross ' Name2];
Counts = [0 0]; % we don't care about the brightness
Valid = true(1,size(Cor_Array,2));
%%% Generates filename
Current_FileName=fullfile(Path,[FileName '_' Name1 '_x_' Name2 '.mcor']);
save(Current_FileName,'Header','Counts','Valid','Cor_Times','Cor_Average','Cor_SEM','Cor_Array');


