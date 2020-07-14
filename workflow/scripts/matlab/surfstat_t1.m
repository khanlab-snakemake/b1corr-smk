clear variables;
matlabpath(pathdef);

%% Datasets
sites   = {'MSTRCHT','SNSX'};

%% Input paths and subjects
base       = '../../';
fid        = fopen([base,'output/surfmorph/labels/participants.tsv']);
subjects   = textscan(fid, '%s', 'HeaderLines', 1);
subjects   = subjects{1,1};
n_subjects = size(subjects,1);
fclose(fid);

i=1; group = [];
for sub=1:2:n_subjects
    site     = strtok(subjects{sub},'-');
    if ~(isequal(site, 'NEW'))
        group{i} = site;
        i = i+1;
    end 
end

age=[]; sex=[]; % 0 = M, 1 = F
for s=1:size(sites,2)
    M = readmatrix([base,sites{s},'_subjects.csv']);
    age = [age; M(:,end-1)];
    sex = [sex; M(:,end)];
end

gender=[];
for sub=1:size(sex,1)
    if sex(sub)==0
        gender{sub}='Male';
    else
        gender{sub}='Female';
    end
end

%% Read in surface
fname   = '../../../../templates/gifti/rh.inflated_164k.vtk';
surf    = read_vtk(fname);
s       = struct();
s.tri   = int32(surf.faces');
s.coord = surf.vertices;

%% Read in data
MSTRCHT_orig = readmatrix([base,'output/freesurfer_data/Maastricht/t1_orig_rh.csv']);
SNSX_orig = readmatrix([base,'output/freesurfer_data/London/t1_orig_rh.csv']);
orig = [MSTRCHT_orig SNSX_orig];

MSTRCHT_corr = readmatrix([base,'output/freesurfer_data/Maastricht/t1_corr_rh.csv']);
SNSX_corr = readmatrix([base,'output/freesurfer_data/London/t1_corr_rh.csv']);
corr = [MSTRCHT_corr SNSX_corr];

%% Stats usting surfstat toolbox
Group    = term(group);
Age      = term(age);
Gender   = term(gender);
Y1       = orig';
Y2       = corr';
contrast = Group.MSTRCHT - Group.SNSX;
M        = 1 + Age + Sex + Group;

meant1subj = mean( double( Y2 ), 2 );
SurfStatPlot( Age, meant1subj );

% Do for original data
slm_orig = SurfStatLinMod( Y1, M, s );
slm_orig = SurfStatT( slm_orig, contrast );
[ pval_orig, peak_orig, clus_orig ] = SurfStatP( slm_orig );
tstat_orig = slm_orig.t';

% Do for corrected data
slm_corr = SurfStatLinMod( Y2, M, s );
slm_corr = SurfStatT( slm_corr, contrast );
[ pval_corr, peak_corr, clus_corr ] = SurfStatP( slm_corr );
tstat_corr = slm_corr.t';

resels = SurfStatResels( slm_orig, pval_orig.mask );
treshold = stat_threshold( resels, length(slm_orig.t), 1, slm_orig.df );

% Save to .mat file
save([base,'output/freesurfer_data/stats_output'],'pval_orig','pval_corr','tstat_orig','tstat_corr');

fname = [base,'output/freesurfer_data/stats_t1.vtk'];
write_vtk(fname,'polydata','triangle',surf.vertices(1,:)',surf.vertices(2,:)',surf.vertices(3,:)',surf.faces',...
    'scalars','ORIG-T-statistic',slm_orig.t', ...
    'scalars','ORIG-P-value',pval_orig.P', ...
    'scalars','ORIG-C-value',pval_orig.C', ...
    'scalars','CORR-T-statistic',slm_corr.t', ...
    'scalars','CORR-P-value',pval_corr.P', ...
    'scalars','CORR-C-value',pval_corr.C');