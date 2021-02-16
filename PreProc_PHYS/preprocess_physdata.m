function preprocess_physdata(inputfile_txt, time_col,vol_col,sf,TR, n_TRs, prefix, output_dir,extra_TRs,inputfile_json)
%% Usage:preprocess_physdata(inputfile_txt, time_col,vol_col,sf,TR, n_TRs, prefix, output_dir,extra_TRs, inputfile_json)
%
% inputfile_txt -----> full path to physiological text file
% time_col ----------> column number in text file that indicates time
% vol_col -----------> column number in text file that indicates volume triggers
% sf ----------------> sampling frequency of input text file in Hz
% TR ----------------> TR of the fMRI aquistion in seconds.
% n_TRs -------------> Total number of TRs (volumes) in the fMRI aquisition.
% prefix ------------> prefix for output file e.g. 'subjectID-scanID'
% output_dir --------> full path to output directory
%
% extra_TRs ---------> OPTIONAL. Extra data to export before & after scan, in units of TR.
%                      Write '0' for no extra data.
% inputfile_json ----> OPTIONAL. Full path to json file (slice timing info) for creating slice triggers.
%                      Write '0' if you do not want slice triggers created.

%
% This function will read in a text file exported from LabChart, and
% save out a text file corresponding to the specific scan to be analyzed.
% For the inputfile_txt, it expects one column to be time and another
% column to be triggers sent from the MRI scanner. The other columns can vary.
% The code binarizes the scanner trigger values, and only codes the onset of
% the trigger in the outputted text file.
%
% THERE IS THE OPTION TO MAKE NEW COLUMNS:
%
% 1. Volume triggers within the extra data exported before & after the scan
% 2. Slice triggers based on the original volume triggers
% 3. Slice triggers based on the original & new volume triggers
%
% If extra_data > 0 and slice triggers are requested:
% - New vol col including the extra data is in [last col - 2]
% - Slice triggers based on original vol col is in [last col - 1]
% - Slice triggers based on new vol col is in [last col]
%
% If extra_data > 0 and slice triggers are not requested:
% - New vol col including the extra data is in [last col]

% If extra_data = 0 and slice triggers are requested:
% - Slice triggers based on original vol col is in [last col]

%% DO CHECKS/SAVES

if nargin ~= 10
    help preprocess_physdata
    return
end

%if output_dir does not end with a '/' add it
if strcmp(output_dir(end),'/') == 0
   output_dir = strcat(output_dir,'/');
else
end

%if output_dir does not exist, make it
if exist(sprintf('%s',output_dir),'dir') == 7
else
    mkdir(sprintf('%s',output_dir))
    fprintf(2,'Output Directory does not exist -> Making Output Directory \n \n')
end

%print inputs to command window
fprintf(1,'\n \n');
fprintf(1,'The following inputs were given to this function:\n \n');
fprintf(1,'Path and name of physiological text file -------------> %s \n\n',inputfile_txt);
fprintf(1,'Time column in physiological text file ---------------> %d \n\n',time_col);
fprintf(1,'Volume column in physiological text file -------------> %d \n\n',vol_col);
fprintf(1,'Sampling frequency of physiological text file (Hz) ---> %d \n\n',sf);
fprintf(1,'Prefix used for all file outputs ---------------------> %s \n\n',prefix);
fprintf(1,'Location of output directory -------------------------> %s \n\n',output_dir);
fprintf(1,'Extra data to export at start & end (number of TRs) --> %d \n\n',extra_TRs);
fprintf(1,'TR of fMRI scan (secs) -------------------------------> %d \n\n',TR);
fprintf(1,'Path and name of json file ---------------------------> %s \n\n',inputfile_json);

%% START

%define number of subplots
if isequal(extra_TRs,0) == 0 && isequal(inputfile_json,0) == 0
    n_subplots = 4;
elseif isequal(extra_TRs,0) == 0 && isequal(inputfile_json,0) == 1
    n_subplots = 2;
elseif isequal(extra_TRs,0) == 1 && isequal(inputfile_json,0) == 0
    n_subplots = 2;
elseif isequal(extra_TRs,0) == 1 && isequal(inputfile_json,0) == 1
    n_subplots = 1;
else
end

%express extra data in seconds
extra_data_secs=extra_TRs*TR;

%load the data file and define some variables
all_data = load(inputfile_txt); %if this doesn't work try'importdata'
vol = all_data(:,vol_col);
time = all_data(:,time_col);

%binarise trigger values and only code onset of trigger
vol_data = all_data(:,vol_col);
vol_data(vol_data <=2) = 0;
vol_data(vol_data >2) = 1;

vol_data_onset = vol_data;

for x = 2:size(vol_data,1)
    if vol_data(x-1,1)~= 0
    vol_data_onset(x,1) = 0;
    else
    vol_data_onset(x,1) = vol_data(x,1);
    end
end

all_data(:,vol_col) = vol_data_onset;

%Find start and end times for the scan
figure('units','normalized','outerposition',[0 0 1 1]); %make the figure the same size as my screen
plot(time,vol);
title('First click: just before the scan.        Second click: just after the scan.','FontSize',30)
[x,~]=ginput(2);
close %close the gui
pause(0.1);

start_time_GUI = x(1);
finish_time_GUI = x(2);

if start_time_GUI > finish_time_GUI
   fprintf(2,'The start time of your scan appears to be after your end time - you probably clicked incorrectly! \n')
   fprintf(2,'Exiting! \n \n')
   return
else
end

all_data_vols = all_data(all_data(:,2)>0);
[~,start_time] = min(abs(all_data_vols(:,1) - start_time_GUI)); %closest time stamp to where the user clicked, only in rows that have a volume trigger
[~,finish_time] = min(abs(all_data_vols(:,1) - finish_time_GUI)); %closet time stamp to where the user clicked, only in rows that have a volume trigger
start_vol_index = find(all_data(:,1) == all_data_vols(start_time)); %index of this time stamp in all_data
finish_vol_index = find(all_data(:,1) == all_data_vols(finish_time));%index of this time stamp in all_data

%extract the specific data from original file
scan_data = all_data(start_vol_index:start_vol_index+(n_TRs*TR*sf)-1,:);

%checks
n_TRs_volcol = sum(scan_data(:,vol_col));
if isequal(n_TRs_volcol,n_TRs) ~=1
fprintf(1, '** Number of volume triggers (TRs) counted = %d. This does not match what you inputed (%d) ** \n \n',n_TRs_volcol, n_TRs)
return
end

scan_data_samples=length(scan_data);
n_TRs_samples = round((scan_data_samples/sf)/TR); %%%%%%%%%%%%%%%%%%%%%%% WORK OUT A BETTER WAY TO SOLVE THIS
n_TRs_inputted=n_TRs;
if isequal(n_TRs_samples,n_TRs_inputted) ~=1
fprintf(1, '** Number of samples extracted in units of TR = (%d/%d)/%d. This does not match what you inputed (%d) ** \n \n',scan_data_samples,sf,TR, n_TRs)
return
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If requested, extract extra phys data before & after the scan. This creates extra volume triggers in a new column %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isequal(extra_TRs,0) == 0

    old_start=start_vol_index;
    old_finish=start_vol_index+(n_TRs*TR*sf)-1;
    scan_data = all_data(old_start - (extra_data_secs * sf):old_finish + (extra_data_secs * sf),:);

    %checks
    scan_data_samples=length(scan_data);
    n_TRs_samples = round((scan_data_samples/sf)/TR); %%%%%%%%%%%%%%%%%%%%%%% WORK OUT A BETTER WAY TO SOLVE THIS
    n_TR_samples_requested=n_TRs+extra_TRs+extra_TRs;
    if isequal(n_TRs_samples,n_TR_samples_requested) ~=1
        fprintf(1, '** With extra_TRs: Number of samples extracted in units of TR = (%d/%d)/%d. This does not match what you requested (%d+%d+%d) ** \n \n',scan_data_samples,sf,TR, n_TRs,extra_TRs,extra_TRs)
        return
    end

    vol_extra = scan_data(:,vol_col);
    scan_data_1 = [scan_data vol_extra];
    [~,colnum] = size(scan_data);
    vol_extra_col = colnum + 1;

    [index,~] = find(scan_data_1(:,vol_extra_col) == 1);
    index_start = index(1);
    index_end = index(end);

    %Find time spacing between volume triggers (should be the TR)
    spacing1 = abs(index(1) - index(2));
    spacing2 = abs(index(2) - index(3));
    spacing3 = abs(index(end-1) - index(end));
    spacing4 = abs(index(end-2) - index(end-1));

    if isequal(spacing1,spacing2,spacing3,spacing4) == 0
       fprintf(2,'The spacing between volume triggers is not equal. \n \n')
       fprintf(2,'Exiting! \n \n')
       return
    else
    end

    spacing_ms = spacing1;
    spacing_seconds = spacing1/1000;

    if isequal(spacing_seconds,TR) == 0
       fprintf(2,'The TR you inputted does not match what is in your data. \n \n')
       fprint(2, 'The extra volume triggers added may be wrong. \n \n')
       fprintf(2,'Exiting! \n \n')
    end

    %append vol triggers to start of data
    for x = 1:extra_TRs
        scan_data_1(index_start-spacing_ms,vol_extra_col) = 1;
        [index,~] = find(scan_data_1(:,vol_extra_col) == 1);
        index_start = index(1);
    end


    %append vol triggers to end of data
    for x = 1:extra_TRs
        scan_data_1(index_end+spacing_ms,vol_extra_col) = 1;
        [index,~] = find(scan_data_1(:,vol_extra_col) == 1);
        index_end = index(end);
    end

    %checks
    n_TRs_volcol=length(find(scan_data_1(:,vol_extra_col) == 1));
    n_TR_samples_requested=n_TRs+extra_TRs+extra_TRs;
    if isequal(n_TRs_volcol,n_TR_samples_requested) ~=1
        fprintf(1, '** With extra_TRs: Number of volume triggers in the new vol_col (%d) does not match what you requested (%d) ** \n \n',n_TRs_volcol,n_TR_samples_requested)
        return
    end

    %Plot data extracted and extra triggers
    figure;
    subplot(n_subplots,1,1)
    plot(time,vol, 'Color', [65/255 105/255 225/255])
    rectangle('Position',[time(start_vol_index - (extra_data_secs*sf)) min(vol)-0.2 length(time(start_vol_index - (extra_data_secs*sf)):time(finish_vol_index + ((extra_data_secs+TR)*sf))) max(vol)+0.2],'EdgeColor','k','LineWidth',2)
    title('Data extracted. If this overlaps with a different scan try changing the extra-TRs input','FontSize',15)

    vols = scan_data(:,vol_col);
    vols(vols(:,1) == 1) = 1.5; %make vol col higher than slice col to stand out

    subplot(n_subplots,1,2)
    plot(scan_data_1(:,time_col),scan_data_1(:,vol_extra_col),'Color', [65/255 105/255 225/255])
    ylim([0 1.55])
    hold on
    plot(scan_data(:,time_col),vols, 'Color', [65/255 105/255 225/255])
    ylim([0 1.55])
    title(sprintf('Volume Triggers - Should be %d extra at the start & end (plotted shorter)',extra_TRs),'FontSize',15)

else
    scan_data_1 = scan_data;

    %Plot data extracted
    figure;
    subplot(n_subplots,1,1)
    plot(time,vol, 'Color', [65/255 105/255 225/255])
    rectangle('Position',[time(start_vol_index) min(vol)-0.2 length(time(start_vol_index):time(start_vol_index+(n_TRs*TR*sf)-1)) max(vol)+0.2],'EdgeColor','k','LineWidth',2)
    title('Data extracted. If this overlaps with a different scan try changing the extra-TRs input','FontSize',15)

end

%% If requested, create pseudo slice triggers and save as new column(s)

if isequal(inputfile_json,0) == 0

    %load the json file to get slice timings
    json_matlab = jsondecode(fileread(inputfile_json));
    Slice_Timing = json_matlab.SliceTiming;
    n_slices = size(Slice_Timing,1);

    json_tr = json_matlab.RepetitionTime;
    if isequal(json_tr,TR) == 0
       fprintf(2,'The TR you inputted does not match what is in the json file. \n \n')
       fprint(2, 'The extra volume triggers added may be wrong. \n \n')
       fprintf(2,'Exiting! \n \n')
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Create slice triggers based on the original vol col %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    vol_col_no = vol_col;

    %load the phsy file and make an extra column
    make_slice_trig = zeros(size(scan_data_1,1),1);
    scan_data_2 = [scan_data_1 make_slice_trig];
    [~,colnum] = size(scan_data_1);
    slice_trig_col = colnum + 1;

    %add in slice triggers - based on the volume triggers recorded in the data and slice timing in json file
   [vol_rows,~] = find(scan_data_1(:,vol_col_no)==1);
    n_vols = size(vol_rows,1);

    biggest_diff = 0;
    for volume = 1:n_vols

    vol_index = vol_rows(volume,1);
    vol_start = scan_data_2(vol_index,time_col);%time stamp of vol start

    %add the rest of the slice triggers
    for slice = 1:n_slices

        slice_trig_time_N = vol_start + Slice_Timing(slice);

        [diff,slice_trigN_time_row] = min(abs(scan_data_2(:,time_col) - slice_trig_time_N)); %find the time point in the data that is closest to the slice timing

        if diff > biggest_diff
            biggest_diff = diff;
        else
        end

        scan_data_2(slice_trigN_time_row,slice_trig_col) = 1; %max(data_new(:,vol_col));
    end

    end

    expected_diff_R = round((1/sf)/2,4);
    biggest_diff_R = round(biggest_diff,4);

    if biggest_diff_R <= expected_diff_R % biggest error in slice timing is within expected range, for this sampling frequency
    else
    fprintf(2,'** Original vol col: Biggest error in slice timings is %d - This is NOT expected with a sampling frequency of %d Hz ** \n \n',biggest_diff,sf)
    fprintf(2,'Exiting! \n \n')
    return
    end

    %Plot result
    vols = scan_data_2(:,vol_col_no);
    vols(vols(:,1) == 1) = 1.5; %make vol col higher than slice col to stand out

    if isequal(extra_TRs,0) == 0
        subplot_index = 3;
    else
        subplot_index = 2;
    end

    subplot(n_subplots,1,subplot_index)
    plot(scan_data_2(:,time_col),vols, 'Color', [65/255 105/255 225/255])
    ylim([0 1.55])
    hold on
    plot(scan_data_2(:,time_col),scan_data_2(:,slice_trig_col), 'Color', [255/255 140/255 0/255])
    ylim([0 1.55])
    title(sprintf('Slice triggers (orange) for each original vol trigger. There should be %d/MBfactor per vol trigger',n_slices),'FontSize',15)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Create slice triggers based on the new vol col - if created %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if isequal(extra_TRs,0) == 0

        vol_col_no = vol_extra_col;

        %load the phsy file and make an extra column
        make_slice_trig = zeros(size(scan_data_2,1),1);
        scan_data_3 = [scan_data_2 make_slice_trig];
        [~,colnum] = size(scan_data_2);
        slice_trig_col = colnum + 1;

        %add in slice triggers - based on the volume triggers recorded in the data and slice timing in json file
        [vol_rows,~] = find(scan_data_2(:,vol_col_no)==1);
        n_vols = size(vol_rows,1);

        biggest_diff = 0;
        for volume = 1:n_vols

        vol_index = vol_rows(volume,1);
        vol_start = scan_data_3(vol_index,time_col);%time stamp of vol start

        %add the rest of the slice triggers
        for slice = 1:n_slices

            slice_trig_time_N = vol_start + Slice_Timing(slice);

            [diff,slice_trigN_time_row] = min(abs(scan_data_3(:,time_col) - slice_trig_time_N)); %find the time point in the data that is closest to the slice timing

            if diff > biggest_diff
                biggest_diff = diff;
            else
            end

            scan_data_3(slice_trigN_time_row,slice_trig_col) = 1; %max(data_new(:,vol_col));
        end

        end

        expected_diff_R = round((1/sf)/2,4);
        biggest_diff_R = round(biggest_diff,4);

        if biggest_diff_R <= expected_diff_R % biggest error in slice timing is within expected range, for this sampling frequency
        else
        fprintf(2,'** New vol col: Biggest error in slice timings is %d - This is NOT expected with a sampling frequency of %d Hz ** \n \n',biggest_diff,sf)
        fprintf(2,'Exiting! \n \n')
        return
        end

        %Plot result
        vols = scan_data_3(:,vol_col_no);
        vols(vols(:,1) == 1) = 1.5; %make vol col higher than slice col to stand out

        subplot(n_subplots,1,subplot_index+1)
        plot(scan_data_3(:,time_col),vols, 'Color', [65/255 105/255 225/255])
        ylim([0 1.55])
        hold on
        plot(scan_data_3(:,time_col),scan_data_3(:,slice_trig_col), 'Color', [255/255 140/255 0/255])
        ylim([0 1.55])
        title(sprintf('Slice triggers (orange) for each ORIG+NEW vol trigger. There should be %d/MBfactor per vol trigger',n_slices),'FontSize',15)

    else
        scan_data_3 = scan_data_2;
    end

else
    scan_data_3 = scan_data_1;
end


%% Save the file in .txt format

%define naming of output file
if isequal(extra_TRs,0) == 1 && isequal(inputfile_json,0) == 1
    out_filename_save = strcat(output_dir, prefix);
elseif isequal(extra_TRs,0) == 0 && isequal(inputfile_json,0) == 1
    out_filename_save = strcat(output_dir, prefix, sprintf('_%dTRPadding',extra_TRs));
elseif isequal(extra_TRs,0) == 0 && isequal(inputfile_json,0) == 0
    out_filename_save = strcat(output_dir, prefix, sprintf('_%dTRPadding_SliceTrig',extra_TRs));
elseif isequal(extra_TRs,0) == 1 && isequal(inputfile_json,0) == 0
    out_filename_save = strcat(output_dir, prefix, '_SliceTrig');
else
end

T = table(scan_data_3); % put the data in table form (easier to then export)
writetable(T, out_filename_save,'Delimiter','\t','WriteVariableNames',0) % write data to text file using tab delimiters and no headers

end
