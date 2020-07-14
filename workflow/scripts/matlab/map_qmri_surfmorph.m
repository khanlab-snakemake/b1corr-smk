clear variables;
matlabpath(pathdef);

%% Datasets
timepoints = {'orig','corr'};
images     = {'T1-gramag-diff'};
space      = 'MNI152NLin2009cAsym'; %'T1w'

%% Segmentation
rois   = {'hippocampus'};
n_rois = size(rois,2);

%% Read in subjects
base_dir         = '../../output/surfmorph';
labels_dir       = [base_dir,'/labels'];
displacement_dir = [base_dir,'/displacement'];
qmri_dir         = [base_dir,'/qmri'];

fid              = fopen([labels_dir,'/participants.tsv']);
subjects         = textscan(fid, '%s', 'HeaderLines', 1);
subjects         = subjects{1,1};
n_subjects       = size(subjects,1);
fclose(fid);

%% Sampling distances
distances = linspace(-0.1,0.1,10); %boundary
distances_direction = {'boundary'};

%% Loop over subjects to load displacement txt files
s=1;
for r=1:n_rois
    % Read in template surface
    surf   = read_vtk([displacement_dir,'/',subjects{1},'/anat/',subjects{1},'_space-MNI152NLin2009cAsym_label-',rois{r},'_surfmorphvec.vtk']);
    tmpl_v = surf.vertices';
    tmpl_f = surf.faces';

    data = zeros(size(tmpl_v,1),n_subjects/2);
    for ss=1:2:n_subjects
        subject = strtok(subjects{ss},'_')
        in_qmri_dir = [qmri_dir,'/',subject,'/anat'];           

        % Iterate through different images to map onto surface
        for i=1:size(images,2)
            in_qmri_file = [in_qmri_dir,'/',subject,'_',images{i},'.nii.gz'];
            in_qmri_nii  = load_nifti(in_qmri_file);
            qmri         = in_qmri_nii.vol;

            for j=1:3
                indices{j} = [0:1:size(qmri,j)-1];
            end
            
            [X,Y,Z] = meshgrid(indices{2},indices{1},indices{3});

            % Read subject's surface data
            surf = read_vtk([displacement_dir,'/',subject,'_orig/anat/',subject,'_orig_space-T1w_label-',rois{r},'_surfmorphinout.vtk']);
            v    = surf.vertices';
            f    = surf.faces';
            nbrs = find_neighbors(f,v);
            N    = vertexnormals(f,v);           
            
            % Bring vertex coordinates into voxel space
            vv = [v ones(size(v,1),1)];
            vv = inv(in_qmri_nii.vox2ras)*vv';
            vv = vv';

            % Map data onto vertices
            avg_qmri = interp3(X,Y,Z,qmri,vv(:,2),vv(:,1),vv(:,3));

            % Smooth surface maps using neighbors
            avg_qmri_smooth = vtksmooth(v,nbrs,avg_qmri,1);
            data(:,s) = avg_qmri_smooth;

            % Write single subject VTK file with smoothed data as scalar
            fname=[qmri_dir,'/',subject,'/anat/',subject,'_corr-orig_space-',space,'_label-',rois{r},'_',images{i},'.vtk'];
            vtkwrite(fname,'polydata','triangle',tmpl_v(:,1),tmpl_v(:,2),tmpl_v(:,3),tmpl_f,'scalars',images{i},avg_qmri_smooth);
        end
        s=s+1;
    end
    
%     maastricht_avg     = nanmean(data(:,1:32),2);
%     new_maastricht_avg = nanmean(data(:,33:34),2);
%     london_avg         = nanmean(data(:,35:end),2);
%     
%     fname=[qmri_dir,'/gramag-diff_space-MNI152NLin2009cAsym_label-',rois{r},'_surfmorph.vtk'];
%     write_vtk(fname,'polydata','triangle',tmpl_v(:,1),tmpl_v(:,2),tmpl_v(:,3),tmpl_f, ...
%         'scalars','maastricht',maastricht_avg, ...
%         'scalars','new_maastricht',new_maastricht_avg, ...
%         'scalars','london',london_avg);    
end