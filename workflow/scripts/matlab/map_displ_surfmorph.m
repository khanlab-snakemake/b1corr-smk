%% Segmentation
rois = {'hippocampus'};
n_rois = size(rois,2);

%% Input paths and subjects
base_dir         = '../../output/surfmorph';
labels_dir       = [base_dir,'/labels'];
displacement_dir = [base_dir,'/displacement'];

fid              = fopen([labels_dir,'/participants.tsv']);
subjects         = textscan(fid, '%s', 'HeaderLines', 1);
subjects         = subjects{1,1};
n_subjects       = size(subjects,1);
fclose(fid);

%% Loop over ROIs and subjects to load displacement txt files
for r=1:n_rois
    % Read in surface
    surf     = read_vtk([displacement_dir,'/',subjects{1},'/anat/',subjects{1},'_space-MNI152NLin2009cAsym_label-',rois{r},'_surfmorphvec.vtk']);
    vertices = surf.vertices';
    faces    = surf.faces';
    
    out_dir = [displacement_dir,'/out/',rois{r},'/displacement'];
    if ~exist(out_dir,'dir') mkdir(out_dir); end 

    data = zeros(size(vertices,1),n_subjects/2);

    i=1;
    for s=1:2:n_subjects
        subject = strtok(subjects{s},'_');
        corr  = dlmread([displacement_dir,'/work/',rois{r},'_iter3/',subject,'_corr/',rois{r},'.surf_inout.txt']);
        orig  = dlmread([displacement_dir,'/work/',rois{r},'_iter3/',subject,'_orig/',rois{r},'.surf_inout.txt']);
        data(:,i) = corr-orig;

        % Write single subject VTK file with scalar data
        fname=[out_dir,'/',subject,'_space-MNI152NLin2009cAsym_label-',rois{r},'_surfmorphinout-diff.vtk'];
        write_vtk(fname,'polydata','triangle',vertices(:,1),vertices(:,2),vertices(:,3),faces,'scalars','corr-orig',data(:,i));
        i=i+1;
    end
    
    maastricht_avg = nanmean(data(:,1:32),2);
    new_maastricht_avg = nanmean(data(:,33:34),2);
    london_avg = nanmean(data(:,35:end),2);
    
    fname=[out_dir,'/maastricht-avg_space-MNI152NLin2009cAsym_label-',rois{r},'_surfmorphinout-diff.vtk'];
    write_vtk(fname,'polydata','triangle',vertices(:,1),vertices(:,2),vertices(:,3),faces, ...
        'scalars','maastricht_displacement',maastricht_avg, ...
        'scalars','new_maastricht_displacement',new_maastricht_avg, ...
        'scalars','london_displacement',london_avg);
   
%     site_std = nanstd(data,0,2);        
%     fname=[out_dir,'/',sites{1,s},'-std_space-MNI152NLin2009cAsym_label-',rois{r},'_surfmorphinout-diff.vtk'];
%     write_vtk(fname,'polydata','triangle',vertices(:,1),vertices(:,2),vertices(:,3),faces,'scalars','corr-uncorr',site_std);    
end