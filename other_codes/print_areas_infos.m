function print_areas_infos(Data,iSubject,iSess)


%% make directories for files

display('Creating directories...')
% Data directories %%%%%%%%%%%%%%%%%%%%%%%%%%%
Dir_dataHits = sprintf('Data/Hits/%d_Subject/%d_Sess/',iSubject,iSess)
if ~exist(Dir_dataHits, 'dir')
    mkdir(Dir_dataHits)
end

file_brain_areas_list = sprintf('%s/brain_areas_list.txt',Dir_dataHits);
file_brain_areas_unique_names = sprintf('%s/brain_areas_unique_names.txt',Dir_dataHits);

filePh = fopen(file_brain_areas_list,'w');
fprintf(filePh,'%s\n',Data.RecordPairMRIlabels{:,1});
fclose(filePh);

uniqueBA = unique(Data.RecordPairMRIlabels);
filePh = fopen(file_brain_areas_unique_names,'w');
fprintf(filePh,'%s\n',uniqueBA{:});
fclose(filePh);

% The following is done manually, for each area. It would be nice to find a
% way to go through the structure fields and print the Electrode index in
% each area automatically 

area_ch_indx = sprintf('%s/dPFC_ch_index.txt',Dir_dataHits);
dlmwrite(area_ch_indx,Data.MRIlabels.dPFC.ElecIndx,'delimiter','\n') 



end