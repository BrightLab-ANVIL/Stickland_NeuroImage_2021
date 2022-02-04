function Shift_hires_CO2regressor(fMRI_file,CO2_file,extra_TRs_half,cut_TRs_start,cut_TRs_end,bulk_shift,save_bulk_hires,fine_shift,save_fine_hires,shift_unit_s,shift_max_s,prefix,output_dir)

% fMRI_file:        Full path to fMRI average time-series (not demeaned).
% CO2_file:         Full path to CO2 regressor, outputted from calc_CO2_regressor (not demeaned).
% extra_TRs_half:   Extra number of TRs you exported within the CO2_file, before/after the scan.
% cut_TRs_start:    Number of TRs to cut from the start of the fMRI input data (this can be 0).
% cut_TRs_end:      Number of TRs to cut from the end of the fMRI input data (this can be 0).
%
% bulk_shift:       Do you want to run bulk-shift? 1 for yes, 0 for no. 
% save_bulk_hires:  1 for yes, 0 for no. By default, the low-res bulk shifted regressor is saved. 
%
% fine_shift:       Do you want to run fine-shift? 1 for yes, 0 for no. 
% save_fine_hires:  1 for yes, 0 for no. By default, the low-res fine shifted regressor is saved. 
% shift_unit_s:     If fine_shift=1. Fine shift unit, in seconds. This should be faster than the new sampling frequency, and shift_unit_s/(1/nf) should be an interger. 
% shift_max_s:      If fine_shift=1. Max shift (in one direction, +ve or -ve). This should be a multiple of shit_unit_s.
%
% prefix:           Output filename for the shifted CO2 regressors. E.g.sub-01_BH_CO2
% output_dir:       Output directory for the shifted CO2 regressors, and plots.
%
% Negative shift means CO2 leads fMRI. Positive shift means CO2 lags fMRI.
% This is inverted later on in the analysis pipeline.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOTE: This code is written based on these inputs: sf of 1000Hz, nf of 40Hz and a TR of 1.2s.
% It should work in principle for other inputs but this has not been tested.
% There are definitely more efficient ways of making theses lagged-regressors :) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sf=1000; % sampling frequency of the CO2_file in Hz
nf=40; %New frequency to down/up sample the CO2/fMRI time-series to, before doing the cross correlation 
%40Hz is pretty arbitrary - just to make the high-res regressor easier to manage because 1000Hz is overkill
TR = 1.2; % Repitition time of your fMRI scan, in seconds.
%% 

%%%%%%%%%%%%%%%%%%%%%%%
% Define/Check Inputs %
%%%%%%%%%%%%%%%%%%%%%%%

%if output_dir does not exist, make it
if exist(sprintf('%s',output_dir),'dir') == 7
else
    mkdir(sprintf('%s',output_dir))
    fprintf(2,'Output Directory does not exist -> Making Output Directory \n \n')
end

extra_TRs_half_samples_sf=extra_TRs_half*TR*sf; %Calculate how many samples correspond to the extra_TRs (one side) at the sf
extra_TRs_half_samples_nf=extra_TRs_half*TR*nf; %Calculate how many samples correspond to the extra_TRs (one side) at the nf
cut_TRs_start_samples_sf=cut_TRs_start*TR*sf; %Calculate how many samples correspond to the TRs to be cut at the sf - start
cut_TRs_start_samples_nf=cut_TRs_start*TR*nf; %Calculate how many samples correspond to the TRs to be cut at the nf - start
cut_TRs_end_samples_sf=cut_TRs_end*TR*sf; %Calculate how many samples correspond to the TRs to be cut at the sf - end
cut_TRs_end_samples_nf=cut_TRs_end*TR*nf; %Calculate how many samples correspond to the TRs to be cut at the nf - end
TR_samples_nf = TR/(1/nf); %Calculate how many samples correspond to the TR at the nf

%Check nf is a multiple of sf
if mod(sf,nf)~=0
    fprintf(2, 'New frequency is not a multiple of sampling frequency. Exiting!\n\n')
    return
else
end

TR_samples_nf=round(TR_samples_nf);
extra_TRs_half_samples_sf=round(extra_TRs_half_samples_sf);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load Data. Cut, demean and resample. %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Load BOLD-fMRI trace, demean, cut and upsample to nf 
fMRI=load(fMRI_file);
n_TRs_total=length(fMRI);
TRstart=cut_TRs_start+1;
TRend=n_TRs_total-cut_TRs_end;
TRstart_string=num2str(TRstart);
TRend_string=num2str(TRend);
fMRI_cut=fMRI(TRstart:TRend); %cut then demean because the first vols may strongly bias the fMRI mean
fMRI_cut_demean=fMRI_cut-mean(fMRI_cut);
n_TRs_process = length(fMRI_cut_demean);
resample_input=round(inv((1/TR)/nf));% E.g. If TR = 1.2 and nf = 40, then (1/1.2)*48=40
fMRI_cut_demean_nf=resample(fMRI_cut_demean,resample_input,1);
fprintf(1,'\n\n The number of TRs in your scan is %d. You are processing %d TRs, from %d to %d. \n\n',n_TRs_total,n_TRs_process,TRstart,TRend);

%Load high-res PETCO2 trace and downsample to nf
CO2=load(CO2_file);
PETCO2=CO2(:,2);
PETCO2_nf = downsample(PETCO2,sf/nf);  %it was PETCO2_nf
 
%Cut PETCO2 to same length as fMRI and downsample to nf
sample_start=extra_TRs_half_samples_sf+cut_TRs_start_samples_sf+1;
sample_end=length(PETCO2)-extra_TRs_half_samples_sf-cut_TRs_end_samples_sf;
sample_start=round(sample_start); 
sample_end=round(sample_end);  
PETCO2_cut=PETCO2(sample_start:sample_end); 
PETCO2_cut_nf = downsample(PETCO2_cut,sf/nf);

%Check that CO2 and fMRI time-series are the same length
if isequal(length(fMRI_cut_demean_nf),length(PETCO2_cut_nf))~=1
    fprintf(1,'\n\nThe fMRI and CO2 inputs for the xcorr are not the same length. Exiting!\n\n');
    return
end

%% 

if bulk_shift == 1
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make the bulk-shifted CO2-regressor %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf(1,'\n** Running cross correlation with the CO2 regressor and fMRI-ts ** \n');

%demean the PETCO2-ts (just for running xcorr - fMRI is already demeaned)
PETCO2_cut_nf_demean=PETCO2_cut_nf - mean(PETCO2_cut_nf); 
%xcorr
[corr, lags] = xcorr(PETCO2_cut_nf_demean,fMRI_cut_demean_nf,round(extra_TRs_half_samples_nf));
[max_corr, max_corr_I] = max(corr); %only positive xcorr
optshift = lags(max_corr_I);
optshift_s = optshift/nf;

%New CO2 with optimized lag (fMRI vols cut at start and end)
PETCO2_nf_shift_cut=PETCO2_nf((extra_TRs_half_samples_nf+cut_TRs_start_samples_nf)+optshift+1:end-((extra_TRs_half_samples_nf+cut_TRs_end_samples_nf)-optshift),1);

if save_bulk_hires==1
    fprintf(1,sprintf('\n** Saving the bulk shifted CO2 regressor at %d Hz ** \n',nf));
    %Output this bulk shifted PETCO2 regressor at high-resolution
    output_bulk_filename=[output_dir '/' prefix '_TR' TRstart_string 'to_TR' TRend_string '_bulkshift_hires.txt'];
    fid=fopen(output_bulk_filename,'wt');
    fprintf(fid,'%f\n',PETCO2_nf_shift_cut');
    fclose(fid);
else
end

% Downsample PETCO2 regressor to TR and then output
TR_samples_nf_index=1;
PETCO2_nf_shift_cut_DS=zeros(n_TRs_process,1);
for x = 1:n_TRs_process
    PETCO2_nf_shift_cut_DS(x,1)=PETCO2_nf_shift_cut(TR_samples_nf_index,1);
    TR_samples_nf_index=TR_samples_nf_index+TR_samples_nf;
end

%save
fprintf(1,'\n** Saving the bulk shifted CO2 regressor at resolution of TR ** \n');
output_bulk_filename=[output_dir '/' prefix '_TR' TRstart_string 'to_TR' TRend_string '_bulkshift.txt'];
fid=fopen(output_bulk_filename,'wt');
fprintf(fid,'%f\n',PETCO2_nf_shift_cut_DS');
fclose(fid);

%Plots

PETCO2_cut_nf_NORM = (PETCO2_cut_nf - min(PETCO2_cut_nf))/(max(PETCO2_cut_nf) - min(PETCO2_cut_nf));
PETCO2_nf_shift_cut_NORM = (PETCO2_nf_shift_cut - min(PETCO2_nf_shift_cut))/(max(PETCO2_nf_shift_cut) - min(PETCO2_nf_shift_cut));
fMRI_cut_demean_nf_NORM = (fMRI_cut_demean_nf - min(fMRI_cut_demean_nf))/(max(fMRI_cut_demean_nf) - min(fMRI_cut_demean_nf));

figure;
subplot(2,1,1)
plot(lags,corr,'k','LineWidth',1)
xlabel('Shift (in samples)')
xlim([min(lags)-5 max(lags)+5])
ylabel('Cross Correlation')
optshift_s_display=num2str(optshift_s);
title(sprintf('Max cross correlation is %0.3g, at a CO2 shift of %s seconds.', max_corr, optshift_s_display))
set(gca,'FontSize',12)
subplot(2,1,2)
plot(PETCO2_cut_nf_NORM,'r','LineWidth',1)
hold on
plot(PETCO2_nf_shift_cut_NORM,'b','LineWidth',1)
hold on
plot(fMRI_cut_demean_nf_NORM,'k','LineWidth',1)
ylabel('Normalized signal intensity')
xlabel('Samples')
xlim([0 length(fMRI_cut_demean_nf_NORM)])
legend('P_E_TCO_2 orig','P_E_TCO_2 shift','BOLD-fMRI')
set(gca,'FontSize',12)

%Save here
output_bulk_filename_plots=[output_dir '/' prefix '_TR' TRstart_string 'to_TR' TRend_string '_bulkshift_PLOTS.png'];
saveas(gcf,output_bulk_filename_plots)

% Make table to save max corr and lag
cell = ({'Prefix' 'MaxCorr' 'OptShift'; prefix, max_corr, optshift_s});
table = cell2table(cell(2,:));
table.Properties.VariableNames = cell(1,:);

% Save table
output_table_filename=[output_dir '/' prefix '_TR' TRstart_string 'to_TR' TRend_string '_bulkshift_TABLE.txt'];
writetable(table,output_table_filename,'Delimiter','\t')

else 
    fprintf(1,'\n** NOT bulk shifting the PETCO2 regressor ** \n');
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make the fine-shifted CO2-regressors %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if fine_shift==1

    if bulk_shift==1
        bulk_optshift=optshift;
        PETCO2_orig=PETCO2_nf_shift_cut;
        output_fine_filename=[output_dir '/' prefix '_TR' TRstart_string 'to_TR' TRend_string '_fineshift.txt'];
        output_fine_filename_hr=[output_dir '/' prefix '_TR' TRstart_string 'to_TR' TRend_string '_fineshift_hires.txt'];
    elseif bulk_shift==0
        bulk_optshift=0;
        PETCO2_orig=PETCO2_cut_nf;
        output_fine_filename=[output_dir '/' prefix '_TR' TRstart_string 'to_TR' TRend_string '_fineshift_NoBulk.txt'];
        output_fine_filename_hr=[output_dir '/' prefix '_TR' TRstart_string 'to_TR' TRend_string '_fineshift_NoBulk_hires.txt'];
    else
    end

    %Check there is enough data to shift
    if shift_max_s > extra_TRs_half*TR
       fprintf(2, 'Not enough extra data exported for a shift max of %d. Exiting!\n\n',shift_max_s)
       return
    else
    end

    %Check the shift unit is larger than the sample frequency
    if shift_unit_s < 1/nf
       fprintf(2, 'Cannot shift by %d with the new frequency of %dHz. Exiting!\n\n',shift_unit_s,nf)
       return
    else
    end

    %Check that the shift unit is a multiple of shift max
    n_shifts = shift_max_s/shift_unit_s;
    %Check if this is a whole number
     if rem(n_shifts,1) ~= 0
        fprintf(2, 'Change shift_unit_s or shift_max_s so that shift_max_s/shift_unit_s is an interger. \n \n')
        return
     else
     end

    %Calculate how many samples corresponds to the shift_unit (at the nf)
     shift_samples = shift_unit_s/(1/nf);
     shift_samples = round(shift_samples); 
    %Check if this is a whole number
     if rem(shift_samples,1) ~= 0
        fprintf(2, 'Change shift_unit so it is a multiple of new frequency (1/%d=interger) or change new frequency. \n \n',nf)
        return
     else
     end

    %%MAKE THE FINE-SHIFTED PETCO2 REGRESSORS%%

    fprintf(1,'\n** Fine shifting the PETCO2 regressor ** \n');

    %Time stamp of the start and end samples of the PETCO2 regressor, based on the bulk shift results
    start_bulkshift=(extra_TRs_half_samples_nf+cut_TRs_start_samples_nf)+bulk_optshift+1;
    end_bulkshift=(extra_TRs_half_samples_nf+cut_TRs_end_samples_nf)-bulk_optshift;

    %Make empty matrix to output shifted CO2 regressors
    PETCO2_shifted = zeros(length(PETCO2_orig),(n_shifts*2)+1);
    %Put the bulk shifted data in the centre
    PETCO2_shifted(:,n_shifts+1) = PETCO2_orig(:,1);

    %Shift the CO2 regressor (negative)
    for x = 1:n_shifts
        n_shift_samples = shift_samples*x;
        SHIFTED = PETCO2_nf(start_bulkshift-n_shift_samples:end-(end_bulkshift+n_shift_samples),1);
        PETCO2_shifted(:,n_shifts-(x-1)) = SHIFTED(:,:);
    end

    %Shift the CO2 regressor (positive)
    for x = 1:n_shifts
        n_shift_samples = shift_samples*x;
        SHIFTED = PETCO2_nf(start_bulkshift+n_shift_samples:end-(end_bulkshift-n_shift_samples),1);
        PETCO2_shifted(:,n_shifts+(x+1)) = SHIFTED(:,:);
    end


    if save_fine_hires==1
    fprintf(1,sprintf('\n** Saving the fine shifted CO2 regressor at %d Hz ** \n',nf));
    %Output this fine shifted PETCO2 regressor at high-resolution
    save(output_fine_filename_hr, 'PETCO2_shifted', '-ascii', '-double', '-tabs');
    else
    end

    %Downsample PETCO2 regressor to TR and then output
    TR_samples_nf_index=1;
    PETCO2_shifted_DS=zeros(n_TRs_process,(n_shifts*2)+1);
    for x = 1:n_TRs_process
        PETCO2_shifted_DS(x,:)=PETCO2_shifted(TR_samples_nf_index,:);
        TR_samples_nf_index=TR_samples_nf_index+TR_samples_nf;
    end

    figure;
    for x = 1:((n_shifts*2)+1)
        plot(PETCO2_shifted_DS(:,x))
        hold on
    end
    title('Fine-shifted CO2 regressors (*Plot not saved*)')
    
    %Save
    fprintf(1,'\n** Saving the fine shifted CO2 regressor at resolution of TR ** \n');
    save(output_fine_filename, 'PETCO2_shifted_DS', '-ascii', '-double', '-tabs');

else
     fprintf(1,'\n** NOT fine shifting the PETCO2 regressor ** \n');
end

end
