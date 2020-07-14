clear variables;
matlabpath(pathdef);

%% Datasets
images  = {'T1-gramag-diff'};
space   = 'MNI152NLin2009cAsym'; %'T1w'

%% ROI
rois   = {'hippocampus'};
n_rois = size(rois,2);

%% Input paths and subjects
base_dir         = '../../output/surfmorph';
labels_dir       = [base_dir,'/labels'];
qmri_dir         = [base_dir,'/qmri'];

fid              = fopen([labels_dir,'/participants.tsv']);
subjects         = textscan(fid, '%s', 'HeaderLines', 1);
subjects         = subjects{1,1};
n_subjects       = size(subjects,1);
fclose(fid);

%% Loop over subjects to load displacement txt files
for r=1:n_rois  
    % Read in surface
    subject = strtok(subjects{1},'_');
    fname   = [qmri_dir,'/',subject,'/anat/',subject,'_corr-orig_space-',space,'_label-',rois{r},'_',images{1},'.vtk']; 
    surf    = read_vtk(fname);
    s       = struct();
    s.tri   = int32(surf.faces');
    s.coord = surf.vertices;
    
    group = [];
    t = []; t_avg = [];
    
    i=1;
    for sub=1:2:n_subjects
        subject  = strtok(subjects{sub},'_');
        site     = strtok(subjects{sub},'-');
        fname    = [qmri_dir,'/',subject,'/anat/',subject,'_corr-orig_space-',space,'_label-',rois{r},'_',images{1},'.vtk'];
        surf     = read_vtk(fname,1);
        t(:,i)   = surf.scalars;
        group{i} = site;
        i = i+1;
    end
    
    maastricht_avg = nanmean(t(:,1:32),2);
    new_maastricht_avg = nanmean(t(:,33:34),2);
    london_avg = nanmean(t(:,35:end),2);    

    %% Stats usting surfstat toolbox
    Group    = term(group);
    Y        = abs(t)';
    contrast = Group.MSTRCHT - Group.SNSX;
    M        = 1 + Group ;
    slm      = SurfStatLinMod( Y, M, s );
    slm      = SurfStatT( slm, contrast );
    [ pval, peak, clus ] = SurfStatP( slm );
    
    fname = [qmri_dir,'/maastricht_gt_london_space-',space,'_label-',rois{r},'_',images{1},'.vtk'];
    write_vtk(fname,'polydata','triangle',surf.vertices(1,:)',surf.vertices(2,:)',surf.vertices(3,:)',surf.faces',...
        'scalars','Maastricht',maastricht_avg, ...
        'scalars','London',london_avg, ... 
        'scalars','Difference',abs(maastricht_avg)-abs(london_avg), ...
        'scalars','T-statistic',slm.t', ...
        'scalars','P-value',pval.P', ...
        'scalars','C-value',pval.C');    
end

    