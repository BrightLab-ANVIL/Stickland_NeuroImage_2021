
function CVRfiles_LagIndex(fMRI_prefix, input_dir, output_dir)

% Files that should be in the input_dir: index for the max Rsq, CO2-Coef, CO2-Tstat and Mean-Coef across shifts.
% Outputs CO2-Coef, CO2-Tstat and Mean-Coef at the optimum shift for each voxel.

% You need this toolbox in your MATLAB path for this code to work:
% https://www.mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image

% Read inputs (parent_input/subject/task)
index = load_untouch_nii(sprintf('%s/%s_Rsq_MAXid1.nii.gz',input_dir,fMRI_prefix));
CO2coef = load_untouch_nii(sprintf('%s/%s_CO2-Coef.nii.gz',input_dir,fMRI_prefix));
MEANcoef = load_untouch_nii(sprintf('%s/%s_Mean-Coef.nii.gz',input_dir,fMRI_prefix));
CO2tstat = load_untouch_nii(sprintf('%s/%s_CO2-Tstat.nii.gz',input_dir,fMRI_prefix));

index_img = double(index.img);
CO2coef_img = double(CO2coef.img);
MEANcoef_img = double(MEANcoef.img);
CO2tstat_img = double(CO2tstat.img);

index_hdr=index.hdr;
% CO2coef_hdr = CO2coef.hdr;
% MEANcoef_hdr = MEANcoef.hdr;
% CO2tstat_hdr = CO2tstat.hdr;

% Process
x_dim=index_hdr.dime.dim(2);
y_dim=index_hdr.dime.dim(3);
z_dim=index_hdr.dime.dim(4);
CO2coef_out = zeros(x_dim, y_dim, z_dim);
MEANcoef_out = zeros(x_dim, y_dim, z_dim);
CO2tstat_out = zeros(x_dim,y_dim,z_dim); 
    
for x = 1:x_dim
    for y = 1:y_dim
        for z = 1:z_dim
            optlag_index = index_img(x,y,z);
            
            if optlag_index ~= 0 
            CO2coef_Opt = CO2coef_img(x,y,z,optlag_index);
            MEANcoef_Opt = MEANcoef_img(x,y,z,optlag_index);
            CO2tstat_Opt = CO2tstat_img(x,y,z,optlag_index);
            CO2coef_out(x,y,z) = CO2coef_Opt;
            MEANcoef_out(x,y,z) = MEANcoef_Opt;
            CO2tstat_out(x,y,z) = CO2tstat_Opt;
            else
            end
                
        end
    end
end

% Save outputs
nii_A.hdr = index_hdr;
nii_A.img = CO2tstat_out;
save_nii(nii_A, sprintf('%s/%s_CO2-Tstat_OptShift.nii.gz',output_dir,fMRI_prefix));

nii_A.hdr = index_hdr;
nii_A.img = CO2coef_out;
save_nii(nii_A, sprintf('%s/%s_CO2-Coef_OptShift.nii.gz',output_dir,fMRI_prefix));

nii_A.hdr = index_hdr;
nii_A.img = MEANcoef_out;
save_nii(nii_A, sprintf('%s/%s_Mean-Coef_OptShift.nii.gz',output_dir,fMRI_prefix));

exit
