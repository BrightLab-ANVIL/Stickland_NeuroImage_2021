function calc_CO2_regressor(filename,fs, vol_col, CO2_col, output_dir, prefix,HRFconv,output_hires,td,demean,peak_check)

% Usage:
%
%    calc_CO2_regressor(filename,fs, vol_col, CO2_col, output_dir, prefix,HRFconv,output_hires,td,demean,peak_check)
%
%    filename ----------> full path to physiological text file
%    fs ----------------> sampling frequency in Hz
%    vol_col -----------> column number in text file that indicates volume triggers
%    CO2_col -----------> column number in text file that indicates CO2 trace 
%    output_dir --------> full path the output directory
%    prefix ------------> prefix for all the output files created
%
%    Put '1' for yes, '0' for no:
%
%    HRFconv -----------> convolve the CO2 regressor with the HRF and output to separate file
%    output_hires ------> output high-resolution regressors (at fs)
%    td ----------------> output a temporal derivative of CO2 regressor
%    demean ------------> demean the regressors, and their temporal derivaties, before writing out
%    peak_check --------> turn on plots to check peak detection accuracy (RECOMMENDED).
%                         PLEASE NOTE: if there is a '_peakcheck' file in the output directory, this
%                         will be re-loaded automatically. If this is not desired, delete this old file.
%
% Data in the input file should be stored in columns as text. The program expects a column of volume triggers to be included also (0 = no trigger, 1 = trigger).
% This script reads in the physiology recording and outputs an end-tidal CO2 regressor.

% Last Updated by Rachael Stickland Jan 2021
% Based on code from Kevin Murphy (CUBRIC, Cardiff University, 2013/2014)
% This is a simplified subset of a larger script, that can also create ETO2, heart-rate, respiration, and RETROICOR regressors. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin ~= 11
   help calc_CO2_regressor
   return
end

fprintf(1,'\n\nRunning the calc_CO2_regressor script with ');

fprintf(1,'\tInput physiology file: %s\n\tOutput prefix: %s\n\tOutput dir: %s\n',filename,prefix,output_dir);
fprintf(1,'The following control constants are being used:\n ');
fprintf(1,'\tfs = %d\t\t-\tSampling rate of input physiology file in Hz\n',fs);
fprintf(1,'\tvol_col = %d\t\t-\tColumn in physiology file containing volume triggers\n',vol_col);
fprintf(1,'\tCO2_col = %d\t\t-\tColumn in physiology file containing end-tidal CO2 trace\n',CO2_col);

if (HRFconv)
    fprintf(1,'CO2 regressor will be convolved with a HRF. \n');
end

if (td)
    fprintf(1,'Temporal derivative of the CO2 regressor will be saved. \n');
end

if (demean)
	fprintf(1,'CO2 regressor (and temporal derivatives, if requested) will be demeaned. \n');
end

fprintf(1,'\nRemember to check these constants carefully!!!\n\n');

CO2_td=td;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Load file and extract relevant columns %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (exist(filename,'file')==0)
	fprintf(1,'\nError! Input physiology file %s does not exist. Exiting\n\n',filename);
	return
end


fprintf(1,'Loading physiology file\n');
A = load(filename);
vol = A(:,vol_col);
ivols = find(vol==1);

CO2 = A(:,CO2_col);

clear A;

num_vols = length(ivols);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% CO2 Peak Detection Algorithm %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf(1,'Detecting end-tidal peaks\n');

frame = round(fs/5); if (round(frame./2) == (frame./2)); frame = frame + 1; end
CO2filt = sgolayfilt(CO2,1,frame);

% Calculate where we are going to assume the subject starts breathing out
x = ((diff(CO2filt)-max(diff(CO2filt))/4)>0).*(diff(CO2filt)-max(diff(CO2filt))/4);
indzero = find(x==0);
indrise = find(x(indzero(1:end-1)+1)>0);
imaxstart = indzero(indrise)+1;

% Calc position of max between the start of one positive lobe and the next
imax = [];
[m ind] = max(CO2filt(1:imaxstart(2)));
imax(1) = ind;
for j=2:length(imaxstart)-1;
    [m ind] = max(CO2filt(imaxstart(j):imaxstart(j+1)));
    imax(j) = imaxstart(j) + ind;
end
[m ind] = max(CO2filt(imaxstart(end):length(CO2filt)));
imax(length(imaxstart)) = imaxstart(end) + ind;
imax = imax'-1;


% Remove spurious values where the max is the same as the value at imaxstart
iremove = find(imaxstart(2:end)==imax(1:end-1));
imax(iremove) = [];

% Remove short breaths where end-tidal CO2 is not reached with a sliding window technique
iremove = [];
if (imax(1) == 0); imax(1) = 1; end
for i = 6:length(imax)-5;
    m = mean([CO2filt(imax(i-5:i-1));CO2filt(imax(i+1:i+5))]);
    l = min(CO2filt(imax(i-5):imax(i+5)));
    if (CO2filt(imax(i)) < (0.5 * (m + l))); iremove = [iremove;i]; end
end
imax(iremove) = [];

% Check for peaks that are within 2s of each other and choose only the highest peak
iremove = [];
i = 1;
while (i < length(imax));
    count = 0;
    j = 1;
    check_inds = [];
    too_long = 0;
    while ((too_long == 0) && (imax(i+j) - imax(i+j-1)) < 2 * fs)
        check_inds = [check_inds;i+j-1];
        j = j+1;
        count = 1;
        if ((i+j) > length(imax)); too_long = 1; end
    end

    if (count == 1)
        check_inds = [check_inds;i+j-1];
        [m ind] = max(CO2filt(imax(check_inds)));
        check_inds(ind) = [];
        iremove = check_inds;
        imax(iremove) = [];
        imaxstart(iremove) = [];
    else
        i = i+1;
    end
end

if (peak_check)

    %%% Read in an existing peakcheck file if it exists in the output directory
    if exist(strcat(output_dir, prefix,'_CO2_peakcheck.txt'),'file') == 2
        peakcheck_file=load(strcat(output_dir, prefix,'_CO2_peakcheck.txt'));
        volpeaks_file=peakcheck_file(:,2);
        imax=find(volpeaks_file==1);
        fprintf(1,'\n ---- CO2 peakcheck file was found in the output directory. Loading in these peaks! ----\n\n');
        pause(2)
    else
    end

    %%% Check peak detection results
    x = zeros(length(CO2filt),1);
    x(imax) = max(CO2filt);
    fig = figure;
set(fig,'units','normalized','outerposition',[0 0 1 1]);
    set(fig,'Toolbar','figure');
    plot(CO2filt)
    title('CHANGE CO2 PEAKS','FontSize',16)
    hold on
    red = plot(x,'r');

    temp_maxenv = [];
    temp_maxenv(1:imax(1)) = CO2filt(imax(1));
    for j = 1:length(imax)-1;
        temp_maxenv(imax(j):imax(j+1)) = ...
            CO2filt(imax(j)) + (CO2filt(imax(j+1))-CO2filt(imax(j))) ...
                    * (((imax(j):imax(j+1))-imax(j)) / (imax(j+1) - imax(j)));
    end
    temp_maxenv(imax(end):length(CO2filt))=CO2filt(imax(end));
    black = plot(temp_maxenv,'k');

    but1 = uicontrol('Style','togglebutton','String','Remove Point','units','normalized','Position',[0.01 0.65 0.1 0.05]);
    but2 = uicontrol('Style','togglebutton','String','Insert Point','units','normalized','Position',[0.01 0.6 0.1 0.05]);
    but3 = uicontrol('Style','togglebutton','String','Insert Coord','units','normalized','Position',[0.01 0.55 0.1 0.05]);
    but4 = uicontrol('Style','togglebutton','String','DONE!','units','normalized','Position',[0.01 0.45 0.1 0.05]);

    loop = 1;
    while (loop == 1);
        while ((get(but1,'Value') == 0) && (get(but2,'Value') == 0) && (get(but3,'Value') == 0) && (get(but4,'Value') == 0))
            pause(0.01)
        end
        if (get(but1,'Value') == 1)
            h = impoint;
            pos = getPosition(h);
            delete(h);
            remove_ind = knnsearch(imax,pos(1));
            imax(remove_ind) = [];
            x = zeros(length(CO2filt),1);
            x(imax) = max(CO2filt);
            delete(red);
            red = plot(x,'r');
            delete(black);
            temp_maxenv = [];
            temp_maxenv(1:imax(1)) = CO2filt(imax(1));
            for j = 1:length(imax)-1;
                temp_maxenv(imax(j):imax(j+1)) = ...
                    CO2filt(imax(j)) + (CO2filt(imax(j+1))-CO2filt(imax(j))) ...
                            * (((imax(j):imax(j+1))-imax(j)) / (imax(j+1) - imax(j)));
            end
            temp_maxenv(imax(end):length(CO2filt))=CO2filt(imax(end));
            black = plot(temp_maxenv,'k');
            delete(but1);
            but1 = uicontrol('Style','togglebutton','String','Remove Point','units','normalized','Position',[0.01 0.65 0.1 0.05]);
        elseif (get(but2,'Value') == 1)
            h = impoint;
            pos = getPosition(h);
            delete(h);
            insert_ind = round(pos(1));
            [temp_max temp_ind] = max(CO2filt(insert_ind-100:insert_ind+100));
            insert_ind = insert_ind-100 + temp_ind;
            imax = sort([imax; insert_ind]);
            x = zeros(length(CO2filt),1);
            x(imax) = max(CO2filt);
            delete(red);
            red = plot(x,'r');
            delete(black);
            temp_maxenv = [];
            temp_maxenv(1:imax(1)) = CO2filt(imax(1));
            for j = 1:length(imax)-1;
                temp_maxenv(imax(j):imax(j+1)) = ...
                    CO2filt(imax(j)) + (CO2filt(imax(j+1))-CO2filt(imax(j))) ...
                            * (((imax(j):imax(j+1))-imax(j)) / (imax(j+1) - imax(j)));
            end
            temp_maxenv(imax(end):length(CO2filt))=CO2filt(imax(end));
            black = plot(temp_maxenv,'k');
            delete(but2);
            but2 = uicontrol('Style','togglebutton','String','Insert Point','units','normalized','Position',[0.01 0.6 0.1 0.05]);
        elseif (get(but3,'Value') == 1)
            h = impoint;
            pos = getPosition(h);
            delete(h);
            insert_ind = round(pos(1));
            imax = sort([imax; insert_ind]);
            CO2filt(insert_ind) = pos(2);
            x = zeros(length(CO2filt),1);
            x(imax) = max(CO2filt);
            delete(red);
            red = plot(x,'r');
            delete(black);
            temp_maxenv = [];
            temp_maxenv(1:imax(1)) = CO2filt(imax(1));
            for j = 1:length(imax)-1;
                temp_maxenv(imax(j):imax(j+1)) = ...
                    CO2filt(imax(j)) + (CO2filt(imax(j+1))-CO2filt(imax(j))) ...
                            * (((imax(j):imax(j+1))-imax(j)) / (imax(j+1) - imax(j)));
            end
            temp_maxenv(imax(end):length(CO2filt))=CO2filt(imax(end));
            black = plot(temp_maxenv,'k');
            delete(but3);
            but3 = uicontrol('Style','togglebutton','String','Insert Coord','units','normalized','Position',[0.01 0.55 0.1 0.05]);
        elseif (get(but4,'Value') == 1)
            loop = 0;
            delete(but1);
            delete(but2);
            delete(but3);
            delete(but4);
            title('End-tidal CO2 peaks','FontSize',16)
        end
    end

    pause(0.1)

    output_filename = strcat(output_dir, prefix,'_CO2_peakcheck.txt');
    fprintf(1,'\tWriting CO2 peaks to output file: %s\n',output_filename);
    fid=fopen(output_filename,'wt');
    volpeaks = zeros(length(CO2filt),1);
    volpeaks(imax) = 1;
        fprintf(fid,'%f\t%d\n',[CO2filt volpeaks]');
end

iCO2peaks = imax;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Calculate and write out CO2 regressor %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(1,'Calculating CO2 regressor\n');

if (exist('iCO2peaks') == 0)
    fprintf(1,'\n\nOutput of CO2 peak detection algorithm does not exist\nExiting!\n\n');
    return
end

maxenv = [];
maxenv(1:iCO2peaks(1)) = CO2filt(iCO2peaks(1));
for j = 1:length(iCO2peaks)-1;
    maxenv(iCO2peaks(j):iCO2peaks(j+1)) = ...
                CO2filt(iCO2peaks(j)) + (CO2filt(iCO2peaks(j+1))-CO2filt(iCO2peaks(j))) ...
                            * (((iCO2peaks(j):iCO2peaks(j+1))-iCO2peaks(j)) / (iCO2peaks(j+1) - iCO2peaks(j)));
end
maxenv(iCO2peaks(end):length(CO2filt))=CO2filt(iCO2peaks(end));
maxenv = maxenv';

CO2regressor = maxenv(ivols);

output_filename = strcat(output_dir, prefix,'_CO2.txt');
fprintf(1,'\tWriting CO2 regressor to output file: %s\n',output_filename);
fid=fopen(output_filename,'wt');
if (CO2_td)
    CO2regressor_td = diff(CO2); CO2regressor_td = CO2regressor_td(ivols);
    if (demean); CO2regressor = CO2regressor - mean(CO2regressor); CO2regressor_td = CO2regressor_td - mean(CO2regressor_td); end
        fprintf(fid,'%f\t%f\n',[CO2regressor CO2regressor_td]');
else
    if (demean); CO2regressor = CO2regressor - mean(CO2regressor); end
        fprintf(fid,'%f\n',CO2regressor');
end
fclose(fid);

if (HRFconv == 1)
    fprintf(1,'\tConvolving CO2 regressor with HRF\n');
    t = 0:1/fs:25;
    HRF = exp(-t) .* ((0.00833333 .* t .^ 5) - (1.27e-13 .* t .^ 15)) ;
    CO2conv = conv([maxenv(1)*ones(60*fs,1);maxenv],HRF)/sum(HRF);
    CO2conv = CO2conv(60*fs+1:60*fs+length(maxenv));
    CO2conv=rescale(CO2conv, min(maxenv),max(maxenv)); %scale to the same units as unconvolved regressor
    CO2convregressor = CO2conv(ivols);
    if (CO2_td); CO2convregressor_td = diff(CO2conv); CO2convregressor_td = CO2convregressor_td(ivols); end

    output_filename = strcat(output_dir, prefix,'_CO2_HRFconv.txt');
    fprintf(1,'\tWriting CO2 HRFconv regressor to output file: %s\n',output_filename);
    fid=fopen(output_filename,'wt');
    if (CO2_td)
        if (demean); CO2convregressor = CO2convregressor - mean(CO2convregressor); CO2convregressor_td = CO2convregressor_td - mean(CO2convregressor_td); end
            fprintf(fid,'%f\t%f\n',[CO2convregressor CO2convregressor_td]');
    else
        if (demean); CO2convregressor = CO2convregressor - mean(CO2convregressor); end
            fprintf(fid,'%f\n',CO2convregressor');
    end
    fclose(fid);
end

if (output_hires)
    output_filename = strcat(output_dir, prefix,'_CO2_hires.txt');
    fprintf(1,'\tWriting high resolution CO2 regressor to output file: %s\n',output_filename);
    fid=fopen(output_filename,'wt');
    if (CO2_td)
        if (HRFconv)
            CO2td = diff(CO2conv); CO2td(end+1) = CO2td(end);
            if (demean); maxenv = maxenv - mean(maxenv); CO2conv = CO2conv - mean(CO2conv(1:length(maxenv))); CO2td = CO2td - mean(CO2td(1:length(maxenv))); end
                fprintf(fid,'%f\t%f\t%f\t%d\n',[maxenv CO2conv CO2td vol]');
        else
            CO2td = diff(maxenv); CO2td(end+1) = CO2td(end);
            if (demean); maxenv = maxenv - mean(maxenv); CO2td = CO2td - mean(CO2td); end
                fprintf(fid,'%f\t%f\t%d\n',[maxenv CO2td vol]');
        end
    else
        if (HRFconv)
            if (demean); maxenv = maxenv - mean(maxenv); CO2conv = CO2conv - mean(CO2conv(1:length(maxenv)));end
                    fprintf(fid,'%f\t%f\t%d\n',[maxenv CO2conv vol]');
        else
            if (demean); maxenv = maxenv - mean(maxenv); end
                fprintf(fid,'%f\t%d\n',[maxenv vol]');
        end
    end
    fclose(fid);
end


fprintf(1,'\nFinished!\n\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% End of calc_CO2_regressor script %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% The following subfunctions are included so that this file will  %%%%%%%%
%%%%%%%% work as a standalone matlab script. This functions have been    %%%%%%%%
%%%%%%%% borrowed from Jimmy Shen (jimmy@rotman-baycrest.on.ca). See     %%%%%%%%
%%%%%%%% http://www.rotman-baycrest.on.ca/~jimmy/NIFTI/ for more details %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (c) 2009, Jimmy Shen                                                %
% All rights reserved.                                                          %
%                                                                               %
% Redistribution and use in source and binary forms, with or without            %
% modification, are permitted provided that the following conditions are        %
% met:                                                                          %
%                                                                               %
%     * Redistributions of source code must retain the above copyright          %
%       notice, this list of conditions and the following disclaimer.           %
%     * Redistributions in binary form must reproduce the above copyright       %
%       notice, this list of conditions and the following disclaimer in         %
%       the documentation and/or other materials provided with the distribution %
%                                                                               %
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"   %
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE     %                                                                               %
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE    %
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE      %
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR           %
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF          %
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS      %
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN       %                                                                              %
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)       %
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE    %
% POSSIBILITY OF SUCH DAMAGE.                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function nii = load_untouch_nii(filename, img_idx, dim5_idx, dim6_idx, dim7_idx, ...
			old_RGB, slice_idx)

   if ~exist('filename','var')
      error('Usage: nii = load_untouch_nii(filename, [img_idx], [dim5_idx], [dim6_idx], [dim7_idx], [old_RGB], [slice_idx])');
   end

   if ~exist('img_idx','var') | isempty(img_idx)
      img_idx = [];
   end

   if ~exist('dim5_idx','var') | isempty(dim5_idx)
      dim5_idx = [];
   end

   if ~exist('dim6_idx','var') | isempty(dim6_idx)
      dim6_idx = [];
   end

   if ~exist('dim7_idx','var') | isempty(dim7_idx)
      dim7_idx = [];
   end

   if ~exist('old_RGB','var') | isempty(old_RGB)
      old_RGB = 0;
   end

   if ~exist('slice_idx','var') | isempty(slice_idx)
      slice_idx = [];
   end

   %  Read the dataset header
   %
   [nii.hdr,nii.filetype,nii.fileprefix,nii.machine] = load_nii_hdr(filename);

   if nii.filetype == 0
      nii.hdr = load_untouch0_nii_hdr(nii.fileprefix,nii.machine);
      nii.ext = [];
   else
      nii.hdr = load_untouch_nii_hdr(nii.fileprefix,nii.machine,nii.filetype);

      %  Read the header extension
      %
      nii.ext = load_nii_ext(filename);
   end

   %  Read the dataset body
   %
   [nii.img,nii.hdr] = load_untouch_nii_img(nii.hdr,nii.filetype,nii.fileprefix, ...
		nii.machine,img_idx,dim5_idx,dim6_idx,dim7_idx,old_RGB,slice_idx);

   %  Perform some of sform/qform transform
   %
%   nii = xform_nii(nii, tolerance, preferredForm);

   nii.untouch = 1;

   return					% load_untouch_nii

%---------------------------------------------------------------------

function [hdr, filetype, fileprefix, machine] = load_nii_hdr(fileprefix)

   if ~exist('fileprefix','var'),
      error('Usage: [hdr, filetype, fileprefix, machine] = load_nii_hdr(filename)');
   end

   machine = 'ieee-le';
   new_ext = 0;

   if findstr('.nii',fileprefix)
      new_ext = 1;
      fileprefix = strrep(fileprefix,'.nii','');
   end

   if findstr('.hdr',fileprefix)
      fileprefix = strrep(fileprefix,'.hdr','');
   end

   if findstr('.img',fileprefix)
      fileprefix = strrep(fileprefix,'.img','');
   end

   if new_ext
      fn = sprintf('%s.nii',fileprefix);

      if ~exist(fn)
         msg = sprintf('Cannot find file "%s.nii".', fileprefix);
         error(msg);
      end
   else
      fn = sprintf('%s.hdr',fileprefix);

      if ~exist(fn)
         msg = sprintf('Cannot find file "%s.hdr".', fileprefix);
         error(msg);
      end
   end

   fid = fopen(fn,'r',machine);

   if fid < 0,
      msg = sprintf('Cannot open file %s.',fn);
      error(msg);
   else
      fseek(fid,0,'bof');

      if fread(fid,1,'int32') == 348
         hdr = read_header(fid);
         fclose(fid);
      else
         fclose(fid);

         %  first try reading the opposite endian to 'machine'
         %
         switch machine,
         case 'ieee-le', machine = 'ieee-be';
         case 'ieee-be', machine = 'ieee-le';
         end

         fid = fopen(fn,'r',machine);

         if fid < 0,
            msg = sprintf('Cannot open file %s.',fn);
            error(msg);
         else
            fseek(fid,0,'bof');

            if fread(fid,1,'int32') ~= 348

               %  Now throw an error
               %
               msg = sprintf('File "%s" is corrupted.',fn);
               error(msg);
            end

            hdr = read_header(fid);
            fclose(fid);
         end
      end
   end

   if strcmp(hdr.hist.magic, 'n+1')
      filetype = 2;
   elseif strcmp(hdr.hist.magic, 'ni1')
      filetype = 1;
   else
      filetype = 0;
   end

   return					% load_nii_hdr


%---------------------------------------------------------------------
function [ dsr ] = read_header(fid)

        %  Original header structures
	%  struct dsr
	%       {
	%       struct header_key hk;            /*   0 +  40       */
	%       struct image_dimension dime;     /*  40 + 108       */
	%       struct data_history hist;        /* 148 + 200       */
	%       };                               /* total= 348 bytes*/

    dsr.hk   = header_key(fid);
    dsr.dime = image_dimension(fid);
    dsr.hist = data_history(fid);

    %  For Analyze data format
    %
    if ~strcmp(dsr.hist.magic, 'n+1') & ~strcmp(dsr.hist.magic, 'ni1')
        dsr.hist.qform_code = 0;
        dsr.hist.sform_code = 0;
    end

    return					% read_header


%---------------------------------------------------------------------
function [ hk ] = header_key(fid)

    fseek(fid,0,'bof');

	%  Original header structures
	%  struct header_key                     /* header key      */
	%       {                                /* off + size      */
	%       int sizeof_hdr                   /*  0 +  4         */
	%       char data_type[10];              /*  4 + 10         */
	%       char db_name[18];                /* 14 + 18         */
	%       int extents;                     /* 32 +  4         */
	%       short int session_error;         /* 36 +  2         */
	%       char regular;                    /* 38 +  1         */
	%       char dim_info;   % char hkey_un0;        /* 39 +  1 */
	%       };                               /* total=40 bytes  */
	%
	% int sizeof_header   Should be 348.
	% char regular        Must be 'r' to indicate that all images and
	%                     volumes are the same size.

    v6 = version;
    if str2num(v6(1))<6
       directchar = '*char';
    else
       directchar = 'uchar=>char';
    end

    hk.sizeof_hdr    = fread(fid, 1,'int32')';	% should be 348!
    hk.data_type     = deblank(fread(fid,10,directchar)');
    hk.db_name       = deblank(fread(fid,18,directchar)');
    hk.extents       = fread(fid, 1,'int32')';
    hk.session_error = fread(fid, 1,'int16')';
    hk.regular       = fread(fid, 1,directchar)';
    hk.dim_info      = fread(fid, 1,'uchar')';

    return					% header_key


%---------------------------------------------------------------------
function [ dime ] = image_dimension(fid)

	%  Original header structures
	%  struct image_dimension
	%       {                                /* off + size      */
	%       short int dim[8];                /* 0 + 16          */
        %       /*
        %           dim[0]      Number of dimensions in database; usually 4.
        %           dim[1]      Image X dimension;  number of *pixels* in an image row.
        %           dim[2]      Image Y dimension;  number of *pixel rows* in slice.
        %           dim[3]      Volume Z dimension; number of *slices* in a volume.
        %           dim[4]      Time points; number of volumes in database
        %       */
	%       float intent_p1;   % char vox_units[4];   /* 16 + 4       */
	%       float intent_p2;   % char cal_units[8];   /* 20 + 4       */
	%       float intent_p3;   % char cal_units[8];   /* 24 + 4       */
	%       short int intent_code;   % short int unused1;   /* 28 + 2 */
	%       short int datatype;              /* 30 + 2          */
	%       short int bitpix;                /* 32 + 2          */
	%       short int slice_start;   % short int dim_un0;   /* 34 + 2 */
	%       float pixdim[8];                 /* 36 + 32         */
	%	/*
	%		pixdim[] specifies the voxel dimensions:
	%		pixdim[1] - voxel width, mm
	%		pixdim[2] - voxel height, mm
	%		pixdim[3] - slice thickness, mm
	%		pixdim[4] - volume timing, in msec
	%					..etc
	%	*/
	%       float vox_offset;                /* 68 + 4          */
	%       float scl_slope;   % float roi_scale;     /* 72 + 4 */
	%       float scl_inter;   % float funused1;      /* 76 + 4 */
	%       short slice_end;   % float funused2;      /* 80 + 2 */
	%       char slice_code;   % float funused2;      /* 82 + 1 */
	%       char xyzt_units;   % float funused2;      /* 83 + 1 */
	%       float cal_max;                   /* 84 + 4          */
	%       float cal_min;                   /* 88 + 4          */
	%       float slice_duration;   % int compressed; /* 92 + 4 */
	%       float toffset;   % int verified;          /* 96 + 4 */
	%       int glmax;                       /* 100 + 4         */
	%       int glmin;                       /* 104 + 4         */
	%       };                               /* total=108 bytes */

    dime.dim        = fread(fid,8,'int16')';
    dime.intent_p1  = fread(fid,1,'float32')';
    dime.intent_p2  = fread(fid,1,'float32')';
    dime.intent_p3  = fread(fid,1,'float32')';
    dime.intent_code = fread(fid,1,'int16')';
    dime.datatype   = fread(fid,1,'int16')';
    dime.bitpix     = fread(fid,1,'int16')';
    dime.slice_start = fread(fid,1,'int16')';
    dime.pixdim     = fread(fid,8,'float32')';
    dime.vox_offset = fread(fid,1,'float32')';
    dime.scl_slope  = fread(fid,1,'float32')';
    dime.scl_inter  = fread(fid,1,'float32')';
    dime.slice_end  = fread(fid,1,'int16')';
    dime.slice_code = fread(fid,1,'uchar')';
    dime.xyzt_units = fread(fid,1,'uchar')';
    dime.cal_max    = fread(fid,1,'float32')';
    dime.cal_min    = fread(fid,1,'float32')';
    dime.slice_duration = fread(fid,1,'float32')';
    dime.toffset    = fread(fid,1,'float32')';
    dime.glmax      = fread(fid,1,'int32')';
    dime.glmin      = fread(fid,1,'int32')';

    return					% image_dimension


%---------------------------------------------------------------------
function [ hist ] = data_history(fid)

	%  Original header structures
	%  struct data_history
	%       {                                /* off + size      */
	%       char descrip[80];                /* 0 + 80          */
	%       char aux_file[24];               /* 80 + 24         */
	%       short int qform_code;            /* 104 + 2         */
	%       short int sform_code;            /* 106 + 2         */
	%       float quatern_b;                 /* 108 + 4         */
	%       float quatern_c;                 /* 112 + 4         */
	%       float quatern_d;                 /* 116 + 4         */
	%       float qoffset_x;                 /* 120 + 4         */
	%       float qoffset_y;                 /* 124 + 4         */
	%       float qoffset_z;                 /* 128 + 4         */
	%       float srow_x[4];                 /* 132 + 16        */
	%       float srow_y[4];                 /* 148 + 16        */
	%       float srow_z[4];                 /* 164 + 16        */
	%       char intent_name[16];            /* 180 + 16        */
	%       char magic[4];   % int smin;     /* 196 + 4         */
	%       };                               /* total=200 bytes */

    v6 = version;
    if str2num(v6(1))<6
       directchar = '*char';
    else
       directchar = 'uchar=>char';
    end

    hist.descrip     = deblank(fread(fid,80,directchar)');
    hist.aux_file    = deblank(fread(fid,24,directchar)');
    hist.qform_code  = fread(fid,1,'int16')';
    hist.sform_code  = fread(fid,1,'int16')';
    hist.quatern_b   = fread(fid,1,'float32')';
    hist.quatern_c   = fread(fid,1,'float32')';
    hist.quatern_d   = fread(fid,1,'float32')';
    hist.qoffset_x   = fread(fid,1,'float32')';
    hist.qoffset_y   = fread(fid,1,'float32')';
    hist.qoffset_z   = fread(fid,1,'float32')';
    hist.srow_x      = fread(fid,4,'float32')';
    hist.srow_y      = fread(fid,4,'float32')';
    hist.srow_z      = fread(fid,4,'float32')';
    hist.intent_name = deblank(fread(fid,16,directchar)');
    hist.magic       = deblank(fread(fid,4,directchar)');

    fseek(fid,253,'bof');
    hist.originator  = fread(fid, 5,'int16')';

    return					% data_history


%---------------------------------------------------------------------

function hdr = load_untouch_nii_hdr(fileprefix, machine, filetype)

   if filetype == 2
      fn = sprintf('%s.nii',fileprefix);

      if ~exist(fn)
         msg = sprintf('Cannot find file "%s.nii".', fileprefix);
         error(msg);
      end
   else
      fn = sprintf('%s.hdr',fileprefix);

      if ~exist(fn)
         msg = sprintf('Cannot find file "%s.hdr".', fileprefix);
         error(msg);
      end
   end

   fid = fopen(fn,'r',machine);

   if fid < 0,
      msg = sprintf('Cannot open file %s.',fn);
      error(msg);
   else
      fseek(fid,0,'bof');
      hdr = read_header1(fid);
      fclose(fid);
   end

   return					% load_nii_hdr


%---------------------------------------------------------------------
function [ dsr ] = read_header1(fid)

        %  Original header structures
	%  struct dsr
	%       {
	%       struct header_key hk;            /*   0 +  40       */
	%       struct image_dimension dime;     /*  40 + 108       */
	%       struct data_history hist;        /* 148 + 200       */
	%       };                               /* total= 348 bytes*/

    dsr.hk   = header_key1(fid);
    dsr.dime = image_dimension1(fid);
    dsr.hist = data_history1(fid);

    %  For Analyze data format
    %
    if ~strcmp(dsr.hist.magic, 'n+1') & ~strcmp(dsr.hist.magic, 'ni1')
        dsr.hist.qform_code = 0;
        dsr.hist.sform_code = 0;
    end

    return					% read_header


%---------------------------------------------------------------------
function [ hk ] = header_key1(fid)

    fseek(fid,0,'bof');

	%  Original header structures
	%  struct header_key                     /* header key      */
	%       {                                /* off + size      */
	%       int sizeof_hdr                   /*  0 +  4         */
	%       char data_type[10];              /*  4 + 10         */
	%       char db_name[18];                /* 14 + 18         */
	%       int extents;                     /* 32 +  4         */
	%       short int session_error;         /* 36 +  2         */
	%       char regular;                    /* 38 +  1         */
	%       char dim_info;   % char hkey_un0;        /* 39 +  1 */
	%       };                               /* total=40 bytes  */
	%
	% int sizeof_header   Should be 348.
	% char regular        Must be 'r' to indicate that all images and
	%                     volumes are the same size.

    v6 = version;
    if str2num(v6(1))<6
       directchar = '*char';
    else
       directchar = 'uchar=>char';
    end

    hk.sizeof_hdr    = fread(fid, 1,'int32')';	% should be 348!
    hk.data_type     = deblank(fread(fid,10,directchar)');
    hk.db_name       = deblank(fread(fid,18,directchar)');
    hk.extents       = fread(fid, 1,'int32')';
    hk.session_error = fread(fid, 1,'int16')';
    hk.regular       = fread(fid, 1,directchar)';
    hk.dim_info      = fread(fid, 1,'uchar')';

    return					% header_key


%---------------------------------------------------------------------
function [ dime ] = image_dimension1(fid)

	%  Original header structures
	%  struct image_dimension
	%       {                                /* off + size      */
	%       short int dim[8];                /* 0 + 16          */
        %       /*
        %           dim[0]      Number of dimensions in database; usually 4.
        %           dim[1]      Image X dimension;  number of *pixels* in an image row.
        %           dim[2]      Image Y dimension;  number of *pixel rows* in slice.
        %           dim[3]      Volume Z dimension; number of *slices* in a volume.
        %           dim[4]      Time points; number of volumes in database
        %       */
	%       float intent_p1;   % char vox_units[4];   /* 16 + 4       */
	%       float intent_p2;   % char cal_units[8];   /* 20 + 4       */
	%       float intent_p3;   % char cal_units[8];   /* 24 + 4       */
	%       short int intent_code;   % short int unused1;   /* 28 + 2 */
	%       short int datatype;              /* 30 + 2          */
	%       short int bitpix;                /* 32 + 2          */
	%       short int slice_start;   % short int dim_un0;   /* 34 + 2 */
	%       float pixdim[8];                 /* 36 + 32         */
	%	/*
	%		pixdim[] specifies the voxel dimensions:
	%		pixdim[1] - voxel width, mm
	%		pixdim[2] - voxel height, mm
	%		pixdim[3] - slice thickness, mm
	%		pixdim[4] - volume timing, in msec
	%					..etc
	%	*/
	%       float vox_offset;                /* 68 + 4          */
	%       float scl_slope;   % float roi_scale;     /* 72 + 4 */
	%       float scl_inter;   % float funused1;      /* 76 + 4 */
	%       short slice_end;   % float funused2;      /* 80 + 2 */
	%       char slice_code;   % float funused2;      /* 82 + 1 */
	%       char xyzt_units;   % float funused2;      /* 83 + 1 */
	%       float cal_max;                   /* 84 + 4          */
	%       float cal_min;                   /* 88 + 4          */
	%       float slice_duration;   % int compressed; /* 92 + 4 */
	%       float toffset;   % int verified;          /* 96 + 4 */
	%       int glmax;                       /* 100 + 4         */
	%       int glmin;                       /* 104 + 4         */
	%       };                               /* total=108 bytes */

    dime.dim        = fread(fid,8,'int16')';
    dime.intent_p1  = fread(fid,1,'float32')';
    dime.intent_p2  = fread(fid,1,'float32')';
    dime.intent_p3  = fread(fid,1,'float32')';
    dime.intent_code = fread(fid,1,'int16')';
    dime.datatype   = fread(fid,1,'int16')';
    dime.bitpix     = fread(fid,1,'int16')';
    dime.slice_start = fread(fid,1,'int16')';
    dime.pixdim     = fread(fid,8,'float32')';
    dime.vox_offset = fread(fid,1,'float32')';
    dime.scl_slope  = fread(fid,1,'float32')';
    dime.scl_inter  = fread(fid,1,'float32')';
    dime.slice_end  = fread(fid,1,'int16')';
    dime.slice_code = fread(fid,1,'uchar')';
    dime.xyzt_units = fread(fid,1,'uchar')';
    dime.cal_max    = fread(fid,1,'float32')';
    dime.cal_min    = fread(fid,1,'float32')';
    dime.slice_duration = fread(fid,1,'float32')';
    dime.toffset    = fread(fid,1,'float32')';
    dime.glmax      = fread(fid,1,'int32')';
    dime.glmin      = fread(fid,1,'int32')';

    return					% image_dimension


%---------------------------------------------------------------------
function [ hist ] = data_history1(fid)

	%  Original header structures
	%  struct data_history
	%       {                                /* off + size      */
	%       char descrip[80];                /* 0 + 80          */
	%       char aux_file[24];               /* 80 + 24         */
	%       short int qform_code;            /* 104 + 2         */
	%       short int sform_code;            /* 106 + 2         */
	%       float quatern_b;                 /* 108 + 4         */
	%       float quatern_c;                 /* 112 + 4         */
	%       float quatern_d;                 /* 116 + 4         */
	%       float qoffset_x;                 /* 120 + 4         */
	%       float qoffset_y;                 /* 124 + 4         */
	%       float qoffset_z;                 /* 128 + 4         */
	%       float srow_x[4];                 /* 132 + 16        */
	%       float srow_y[4];                 /* 148 + 16        */
	%       float srow_z[4];                 /* 164 + 16        */
	%       char intent_name[16];            /* 180 + 16        */
	%       char magic[4];   % int smin;     /* 196 + 4         */
	%       };                               /* total=200 bytes */

    v6 = version;
    if str2num(v6(1))<6
       directchar = '*char';
    else
       directchar = 'uchar=>char';
    end

    hist.descrip     = deblank(fread(fid,80,directchar)');
    hist.aux_file    = deblank(fread(fid,24,directchar)');
    hist.qform_code  = fread(fid,1,'int16')';
    hist.sform_code  = fread(fid,1,'int16')';
    hist.quatern_b   = fread(fid,1,'float32')';
    hist.quatern_c   = fread(fid,1,'float32')';
    hist.quatern_d   = fread(fid,1,'float32')';
    hist.qoffset_x   = fread(fid,1,'float32')';
    hist.qoffset_y   = fread(fid,1,'float32')';
    hist.qoffset_z   = fread(fid,1,'float32')';
    hist.srow_x      = fread(fid,4,'float32')';
    hist.srow_y      = fread(fid,4,'float32')';
    hist.srow_z      = fread(fid,4,'float32')';
    hist.intent_name = deblank(fread(fid,16,directchar)');
    hist.magic       = deblank(fread(fid,4,directchar)');

    return					% data_history

%---------------------------------------------------------------------

function ext = load_nii_ext(fileprefix)

   if ~exist('fileprefix','var'),
      error('Usage: ext = load_nii_ext(filename)');
   end

   machine = 'ieee-le';
   new_ext = 0;

   if findstr('.nii',fileprefix)
      new_ext = 1;
      fileprefix = strrep(fileprefix,'.nii','');
   end

   if findstr('.hdr',fileprefix)
      fileprefix = strrep(fileprefix,'.hdr','');
   end

   if findstr('.img',fileprefix)
      fileprefix = strrep(fileprefix,'.img','');
   end

   if new_ext
      fn = sprintf('%s.nii',fileprefix);

      if ~exist(fn)
         msg = sprintf('Cannot find file "%s.nii".', fileprefix);
         error(msg);
      end
   else
      fn = sprintf('%s.hdr',fileprefix);

      if ~exist(fn)
         msg = sprintf('Cannot find file "%s.hdr".', fileprefix);
         error(msg);
      end
   end

   fid = fopen(fn,'r',machine);
   vox_offset = 0;

   if fid < 0,
      msg = sprintf('Cannot open file %s.',fn);
      error(msg);
   else
      fseek(fid,0,'bof');

      if fread(fid,1,'int32') == 348
         if new_ext
            fseek(fid,108,'bof');
            vox_offset = fread(fid,1,'float32');
         end

         ext = read_extension(fid, vox_offset);
         fclose(fid);
      else
         fclose(fid);

         %  first try reading the opposite endian to 'machine'
         %
         switch machine,
         case 'ieee-le', machine = 'ieee-be';
         case 'ieee-be', machine = 'ieee-le';
         end

         fid = fopen(fn,'r',machine);

         if fid < 0,
            msg = sprintf('Cannot open file %s.',fn);
            error(msg);
         else
            fseek(fid,0,'bof');

            if fread(fid,1,'int32') ~= 348

               %  Now throw an error
               %
               msg = sprintf('File "%s" is corrupted.',fn);
               error(msg);
            end

            if new_ext
               fseek(fid,108,'bof');
               vox_offset = fread(fid,1,'float32');
            end

            ext = read_extension(fid, vox_offset);
            fclose(fid);
         end
      end
   end

   return                                       % load_nii_ext


%---------------------------------------------------------------------
function ext = read_extension(fid, vox_offset)

   ext = [];

   if vox_offset
      end_of_ext = vox_offset;
   else
      fseek(fid, 0, 'eof');
      end_of_ext = ftell(fid);
   end

   if end_of_ext > 352
      fseek(fid, 348, 'bof');
      ext.extension = fread(fid,4)';
   end

   if isempty(ext) | ext.extension(1) == 0
      ext = [];
      return;
   end

   i = 1;

   while(ftell(fid) < end_of_ext)
      ext.section(i).esize = fread(fid,1,'int32');
      ext.section(i).ecode = fread(fid,1,'int32');
      ext.section(i).edata = char(fread(fid,ext.section(i).esize-8)');
      i = i + 1;
   end

   ext.num_ext = length(ext.section);

   return                                               % read_extension
%---------------------------------------------------------------------

function [img,hdr] = load_untouch_nii_img(hdr,filetype,fileprefix,machine,img_idx,dim5_idx,dim6_idx,dim7_idx,old_RGB,slice_idx)

   if ~exist('hdr','var') | ~exist('filetype','var') | ~exist('fileprefix','var') | ~exist('machine','var')
      error('Usage: [img,hdr] = load_nii_img(hdr,filetype,fileprefix,machine,[img_idx],[dim5_idx],[dim6_idx],[dim7_idx],[old_RGB],[slice_idx]);');
   end

   if ~exist('img_idx','var') | isempty(img_idx) | hdr.dime.dim(5)<1
      img_idx = [];
   end

   if ~exist('dim5_idx','var') | isempty(dim5_idx) | hdr.dime.dim(6)<1
      dim5_idx = [];
   end

   if ~exist('dim6_idx','var') | isempty(dim6_idx) | hdr.dime.dim(7)<1
      dim6_idx = [];
   end

   if ~exist('dim7_idx','var') | isempty(dim7_idx) | hdr.dime.dim(8)<1
      dim7_idx = [];
   end

   if ~exist('old_RGB','var') | isempty(old_RGB)
      old_RGB = 0;
   end

   if ~exist('slice_idx','var') | isempty(slice_idx) | hdr.dime.dim(4)<1
      slice_idx = [];
   end

   %  check img_idx
   %
   if ~isempty(img_idx) & ~isnumeric(img_idx)
      error('"img_idx" should be a numerical array.');
   end

   if length(unique(img_idx)) ~= length(img_idx)
      error('Duplicate image index in "img_idx"');
   end

   if ~isempty(img_idx) & (min(img_idx) < 1 | max(img_idx) > hdr.dime.dim(5))
      max_range = hdr.dime.dim(5);

      if max_range == 1
         error(['"img_idx" should be 1.']);
      else
         range = ['1 ' num2str(max_range)];
         error(['"img_idx" should be an integer within the range of [' range '].']);
      end
   end

   %  check dim5_idx
   %
   if ~isempty(dim5_idx) & ~isnumeric(dim5_idx)
      error('"dim5_idx" should be a numerical array.');
   end

   if length(unique(dim5_idx)) ~= length(dim5_idx)
      error('Duplicate index in "dim5_idx"');
   end

   if ~isempty(dim5_idx) & (min(dim5_idx) < 1 | max(dim5_idx) > hdr.dime.dim(6))
      max_range = hdr.dime.dim(6);

      if max_range == 1
         error(['"dim5_idx" should be 1.']);
      else
         range = ['1 ' num2str(max_range)];
         error(['"dim5_idx" should be an integer within the range of [' range '].']);
      end
   end

   %  check dim6_idx
   %
   if ~isempty(dim6_idx) & ~isnumeric(dim6_idx)
      error('"dim6_idx" should be a numerical array.');
   end

   if length(unique(dim6_idx)) ~= length(dim6_idx)
      error('Duplicate index in "dim6_idx"');
   end

   if ~isempty(dim6_idx) & (min(dim6_idx) < 1 | max(dim6_idx) > hdr.dime.dim(7))
      max_range = hdr.dime.dim(7);

      if max_range == 1
         error(['"dim6_idx" should be 1.']);
      else
         range = ['1 ' num2str(max_range)];
         error(['"dim6_idx" should be an integer within the range of [' range '].']);
      end
   end

   %  check dim7_idx
   %
   if ~isempty(dim7_idx) & ~isnumeric(dim7_idx)
      error('"dim7_idx" should be a numerical array.');
   end

   if length(unique(dim7_idx)) ~= length(dim7_idx)
      error('Duplicate index in "dim7_idx"');
   end

   if ~isempty(dim7_idx) & (min(dim7_idx) < 1 | max(dim7_idx) > hdr.dime.dim(8))
      max_range = hdr.dime.dim(8);

      if max_range == 1
         error(['"dim7_idx" should be 1.']);
      else
         range = ['1 ' num2str(max_range)];
         error(['"dim7_idx" should be an integer within the range of [' range '].']);
      end
   end

   %  check slice_idx
   %
   if ~isempty(slice_idx) & ~isnumeric(slice_idx)
      error('"slice_idx" should be a numerical array.');
   end

   if length(unique(slice_idx)) ~= length(slice_idx)
      error('Duplicate index in "slice_idx"');
   end

   if ~isempty(slice_idx) & (min(slice_idx) < 1 | max(slice_idx) > hdr.dime.dim(4))
      max_range = hdr.dime.dim(4);

      if max_range == 1
         error(['"slice_idx" should be 1.']);
      else
         range = ['1 ' num2str(max_range)];
         error(['"slice_idx" should be an integer within the range of [' range '].']);
      end
   end

   [img,hdr] = read_image(hdr,filetype,fileprefix,machine,img_idx,dim5_idx,dim6_idx,dim7_idx,old_RGB,slice_idx);

   return					% load_nii_img


%---------------------------------------------------------------------
function [img,hdr] = read_image(hdr,filetype,fileprefix,machine,img_idx,dim5_idx,dim6_idx,dim7_idx,old_RGB,slice_idx)

   switch filetype
   case {0, 1}
      fn = [fileprefix '.img'];
   case 2
      fn = [fileprefix '.nii'];
   end

   fid = fopen(fn,'r',machine);

   if fid < 0,
      msg = sprintf('Cannot open file %s.',fn);
      error(msg);
   end

   %  Set bitpix according to datatype
   %
   %  /*Acceptable values for datatype are*/
   %
   %     0 None                     (Unknown bit per voxel) % DT_NONE, DT_UNKNOWN
   %     1 Binary                         (ubit1, bitpix=1) % DT_BINARY
   %     2 Unsigned char         (uchar or uint8, bitpix=8) % DT_UINT8, NIFTI_TYPE_UINT8
   %     4 Signed short                  (int16, bitpix=16) % DT_INT16, NIFTI_TYPE_INT16
   %     8 Signed integer                (int32, bitpix=32) % DT_INT32, NIFTI_TYPE_INT32
   %    16 Floating point    (single or float32, bitpix=32) % DT_FLOAT32, NIFTI_TYPE_FLOAT32
   %    32 Complex, 2 float32      (Use float32, bitpix=64) % DT_COMPLEX64, NIFTI_TYPE_COMPLEX64
   %    64 Double precision  (double or float64, bitpix=64) % DT_FLOAT64, NIFTI_TYPE_FLOAT64
   %   128 uint8 RGB                 (Use uint8, bitpix=24) % DT_RGB24, NIFTI_TYPE_RGB24
   %   256 Signed char            (schar or int8, bitpix=8) % DT_INT8, NIFTI_TYPE_INT8
   %   511 Single RGB              (Use float32, bitpix=96) % DT_RGB96, NIFTI_TYPE_RGB96
   %   512 Unsigned short               (uint16, bitpix=16) % DT_UNINT16, NIFTI_TYPE_UNINT16
   %   768 Unsigned integer             (uint32, bitpix=32) % DT_UNINT32, NIFTI_TYPE_UNINT32
   %  1024 Signed long long              (int64, bitpix=64) % DT_INT64, NIFTI_TYPE_INT64
   %  1280 Unsigned long long           (uint64, bitpix=64) % DT_UINT64, NIFTI_TYPE_UINT64
   %  1536 Long double, float128  (Unsupported, bitpix=128) % DT_FLOAT128, NIFTI_TYPE_FLOAT128
   %  1792 Complex128, 2 float64  (Use float64, bitpix=128) % DT_COMPLEX128, NIFTI_TYPE_COMPLEX128
   %  2048 Complex256, 2 float128 (Unsupported, bitpix=256) % DT_COMPLEX128, NIFTI_TYPE_COMPLEX128
   %
   switch hdr.dime.datatype
   case   1,
      hdr.dime.bitpix = 1;  precision = 'ubit1';
   case   2,
      hdr.dime.bitpix = 8;  precision = 'uint8';
   case   4,
      hdr.dime.bitpix = 16; precision = 'int16';
   case   8,
      hdr.dime.bitpix = 32; precision = 'int32';
   case  16,
      hdr.dime.bitpix = 32; precision = 'float32';
   case  32,
      hdr.dime.bitpix = 64; precision = 'float32';
   case  64,
      hdr.dime.bitpix = 64; precision = 'float64';
   case 128,
      hdr.dime.bitpix = 24; precision = 'uint8';
   case 256
      hdr.dime.bitpix = 8;  precision = 'int8';
   case 511
      hdr.dime.bitpix = 96; precision = 'float32';
   case 512
      hdr.dime.bitpix = 16; precision = 'uint16';
   case 768
      hdr.dime.bitpix = 32; precision = 'uint32';
   case 1024
      hdr.dime.bitpix = 64; precision = 'int64';
   case 1280
      hdr.dime.bitpix = 64; precision = 'uint64';
   case 1792,
      hdr.dime.bitpix = 128; precision = 'float64';
   otherwise
      error('This datatype is not supported');
   end

   tmp = hdr.dime.dim(2:end);
   tmp(find(tmp < 1)) = 1;
   hdr.dime.dim(2:end) = tmp;

   %  move pointer to the start of image block
   %
   switch filetype
   case {0, 1}
      fseek(fid, 0, 'bof');
   case 2
      fseek(fid, hdr.dime.vox_offset, 'bof');
   end

   %  Load whole image block for old Analyze format or binary image;
   %  otherwise, load images that are specified in img_idx, dim5_idx,
   %  dim6_idx, and dim7_idx
   %
   %  For binary image, we have to read all because pos can not be
   %  seeked in bit and can not be calculated the way below.
   %
   if hdr.dime.datatype == 1 | isequal(hdr.dime.dim(4:8),ones(1,5)) | ...
	(isempty(img_idx) & isempty(dim5_idx) & isempty(dim6_idx) & isempty(dim7_idx) & isempty(slice_idx))

      %  For each frame, precision of value will be read
      %  in img_siz times, where img_siz is only the
      %  dimension size of an image, not the byte storage
      %  size of an image.
      %
      img_siz = prod(hdr.dime.dim(2:8));

      %  For complex float32 or complex float64, voxel values
      %  include [real, imag]
      %
      if hdr.dime.datatype == 32 | hdr.dime.datatype == 1792
         img_siz = img_siz * 2;
      end

      %MPH: For RGB24, voxel values include 3 separate color planes
      %
      if hdr.dime.datatype == 128 | hdr.dime.datatype == 511
	 img_siz = img_siz * 3;
      end

      img = fread(fid, img_siz, sprintf('*%s',precision));

      d1 = hdr.dime.dim(2);
      d2 = hdr.dime.dim(3);
      d3 = hdr.dime.dim(4);
      d4 = hdr.dime.dim(5);
      d5 = hdr.dime.dim(6);
      d6 = hdr.dime.dim(7);
      d7 = hdr.dime.dim(8);

      if isempty(slice_idx)
         slice_idx = 1:d3;
      end

      if isempty(img_idx)
         img_idx = 1:d4;
      end

      if isempty(dim5_idx)
         dim5_idx = 1:d5;
      end

      if isempty(dim6_idx)
         dim6_idx = 1:d6;
      end

      if isempty(dim7_idx)
         dim7_idx = 1:d7;
      end
   else

      d1 = hdr.dime.dim(2);
      d2 = hdr.dime.dim(3);
      d3 = hdr.dime.dim(4);
      d4 = hdr.dime.dim(5);
      d5 = hdr.dime.dim(6);
      d6 = hdr.dime.dim(7);
      d7 = hdr.dime.dim(8);

      if isempty(slice_idx)
         slice_idx = 1:d3;
      end

      if isempty(img_idx)
         img_idx = 1:d4;
      end

      if isempty(dim5_idx)
         dim5_idx = 1:d5;
      end

      if isempty(dim6_idx)
         dim6_idx = 1:d6;
      end

      if isempty(dim7_idx)
         dim7_idx = 1:d7;
      end

      %ROMAN: begin
      roman = 1;
      if(roman)

         %  compute size of one slice
         %
         img_siz = prod(hdr.dime.dim(2:3));

         %  For complex float32 or complex float64, voxel values
         %  include [real, imag]
         %
         if hdr.dime.datatype == 32 | hdr.dime.datatype == 1792
            img_siz = img_siz * 2;
         end

         %MPH: For RGB24, voxel values include 3 separate color planes
         %
         if hdr.dime.datatype == 128 | hdr.dime.datatype == 511
            img_siz = img_siz * 3;
         end

         % preallocate img
         img = zeros(img_siz, length(slice_idx)*length(img_idx)*length(dim5_idx)*length(dim6_idx)*length(dim7_idx) );
         currentIndex = 1;
      else
        img = [];
      end; %if(roman)
      % ROMAN: end

      for i7=1:length(dim7_idx)
         for i6=1:length(dim6_idx)
            for i5=1:length(dim5_idx)
               for t=1:length(img_idx)
               for s=1:length(slice_idx)

                  %  Position is seeked in bytes. To convert dimension size
                  %  to byte storage size, hdr.dime.bitpix/8 will be
                  %  applied.
                  %
                  pos = sub2ind([d1 d2 d3 d4 d5 d6 d7], 1, 1, slice_idx(s), ...
			                    img_idx(t), dim5_idx(i5),dim6_idx(i6),dim7_idx(i7)) -1;
                  pos = pos * hdr.dime.bitpix/8;

                  % ROMAN: begin
                  if(roman)
                      % do nothing
                  else
                     img_siz = prod(hdr.dime.dim(2:3));

                     %  For complex float32 or complex float64, voxel values
                     %  include [real, imag]
                     %
                     if hdr.dime.datatype == 32 | hdr.dime.datatype == 1792
                        img_siz = img_siz * 2;
                     end

                     %MPH: For RGB24, voxel values include 3 separate color planes
                     %
                     if hdr.dime.datatype == 128 | hdr.dime.datatype == 511
                        img_siz = img_siz * 3;
                     end
                  end; % if (roman)
                  % ROMAN: end

                  if filetype == 2
                     fseek(fid, pos + hdr.dime.vox_offset, 'bof');
                  else
                     fseek(fid, pos, 'bof');
                  end

                  %  For each frame, fread will read precision of value
                  %  in img_siz times
                  %
                  % ROMAN: begin
                  if(roman)
                     img(:,currentIndex) = fread(fid, img_siz, sprintf('*%s',precision));
                     currentIndex = currentIndex +1;
                  else
                     img = [img fread(fid, img_siz, sprintf('*%s',precision))];
                  end; %if(roman)
                  % ROMAN: end

               end
               end
            end
         end
      end
   end

   %  For complex float32 or complex float64, voxel values
   %  include [real, imag]
   %
   if hdr.dime.datatype == 32 | hdr.dime.datatype == 1792
      img = reshape(img, [2, length(img)/2]);
      img = complex(img(1,:)', img(2,:)');
   end

   fclose(fid);

   %  Update the global min and max values
   %
   hdr.dime.glmax = double(max(img(:)));
   hdr.dime.glmin = double(min(img(:)));

   %  old_RGB treat RGB slice by slice, now it is treated voxel by voxel
   %
   if old_RGB & hdr.dime.datatype == 128 & hdr.dime.bitpix == 24
      % remove squeeze
      img = (reshape(img, [hdr.dime.dim(2:3) 3 length(slice_idx) length(img_idx) length(dim5_idx) length(dim6_idx) length(dim7_idx)]));
      img = permute(img, [1 2 4 3 5 6 7 8]);
   elseif hdr.dime.datatype == 128 & hdr.dime.bitpix == 24
      % remove squeeze
      img = (reshape(img, [3 hdr.dime.dim(2:3) length(slice_idx) length(img_idx) length(dim5_idx) length(dim6_idx) length(dim7_idx)]));
      img = permute(img, [2 3 4 1 5 6 7 8]);
   elseif hdr.dime.datatype == 511 & hdr.dime.bitpix == 96
      img = double(img(:));
      img = single((img - min(img))/(max(img) - min(img)));
      % remove squeeze
      img = (reshape(img, [3 hdr.dime.dim(2:3) length(slice_idx) length(img_idx) length(dim5_idx) length(dim6_idx) length(dim7_idx)]));
      img = permute(img, [2 3 4 1 5 6 7 8]);
   else
      % remove squeeze
      img = (reshape(img, [hdr.dime.dim(2:3) length(slice_idx) length(img_idx) length(dim5_idx) length(dim6_idx) length(dim7_idx)]));
   end

   if ~isempty(slice_idx)
      hdr.dime.dim(4) = length(slice_idx);
   end

   if ~isempty(img_idx)
      hdr.dime.dim(5) = length(img_idx);
   end

   if ~isempty(dim5_idx)
      hdr.dime.dim(6) = length(dim5_idx);
   end

   if ~isempty(dim6_idx)
      hdr.dime.dim(7) = length(dim6_idx);
   end

   if ~isempty(dim7_idx)
      hdr.dime.dim(8) = length(dim7_idx);
   end

   return						% read_image

%------------------------------------------------
function save_untouch_nii(nii, filename)

   if ~exist('nii','var') | isempty(nii) | ~isfield(nii,'hdr') | ...
	~isfield(nii,'img') | ~exist('filename','var') | isempty(filename)

      error('Usage: save_untouch_nii(nii, filename)');
   end

   if ~isfield(nii,'untouch') | nii.untouch == 0
      error('Usage: please use ''save_nii.m'' for the modified structure.');
   end

   if isfield(nii.hdr.hist,'magic') & strcmp(nii.hdr.hist.magic(1:3),'ni1')
      filetype = 1;
   elseif isfield(nii.hdr.hist,'magic') & strcmp(nii.hdr.hist.magic(1:3),'n+1')
      filetype = 2;
   else
      filetype = 0;
   end

   [p,f] = fileparts(filename);
   fileprefix = fullfile(p, f);

   write_nii(nii, filetype, fileprefix);

%   %  So earlier versions of SPM can also open it with correct originator
 %  %
  % if filetype == 0
   %   M=[[diag(nii.hdr.dime.pixdim(2:4)) -[nii.hdr.hist.originator(1:3).*nii.hdr.dime.pixdim(2:4)]'];[0 0 0 1]];
    %  save(fileprefix, 'M');
%   elseif filetype == 1
 %     M=[];
  %    save(fileprefix, 'M');
   %end

   return					% save_untouch_nii


%-----------------------------------------------------------------------------------
function write_nii(nii, filetype, fileprefix)

   hdr = nii.hdr;

   if isfield(nii,'ext') & ~isempty(nii.ext)
      ext = nii.ext;
      [ext, esize_total] = verify_nii_ext(ext);
   else
      ext = [];
   end

   switch double(hdr.dime.datatype),
   case   1,
      hdr.dime.bitpix = int16(1 ); precision = 'ubit1';
   case   2,
      hdr.dime.bitpix = int16(8 ); precision = 'uint8';
   case   4,
      hdr.dime.bitpix = int16(16); precision = 'int16';
   case   8,
      hdr.dime.bitpix = int16(32); precision = 'int32';
   case  16,
      hdr.dime.bitpix = int16(32); precision = 'float32';
   case  32,
      hdr.dime.bitpix = int16(64); precision = 'float32';
   case  64,
      hdr.dime.bitpix = int16(64); precision = 'float64';
   case 128,
      hdr.dime.bitpix = int16(24); precision = 'uint8';
   case 256
      hdr.dime.bitpix = int16(8 ); precision = 'int8';
   case 512
      hdr.dime.bitpix = int16(16); precision = 'uint16';
   case 768
      hdr.dime.bitpix = int16(32); precision = 'uint32';
   case 1024
      hdr.dime.bitpix = int16(64); precision = 'int64';
   case 1280
      hdr.dime.bitpix = int16(64); precision = 'uint64';
   case 1792,
      hdr.dime.bitpix = int16(128); precision = 'float64';
   otherwise
      error('This datatype is not supported');
   end

%   hdr.dime.glmax = round(double(max(nii.img(:))));
 %  hdr.dime.glmin = round(double(min(nii.img(:))));

   if filetype == 2
      fid = fopen(sprintf('%s.nii',fileprefix),'w');

      if fid < 0,
         msg = sprintf('Cannot open file %s.nii.',fileprefix);
         error(msg);
      end

      hdr.dime.vox_offset = 352;

      if ~isempty(ext)
         hdr.dime.vox_offset = hdr.dime.vox_offset + esize_total;
      end

      hdr.hist.magic = 'n+1';
      save_untouch_nii_hdr(hdr, fid);

      if ~isempty(ext)
         save_nii_ext(ext, fid);
      end
   elseif filetype == 1
      fid = fopen(sprintf('%s.hdr',fileprefix),'w');

      if fid < 0,
         msg = sprintf('Cannot open file %s.hdr.',fileprefix);
         error(msg);
      end

      hdr.dime.vox_offset = 0;
      hdr.hist.magic = 'ni1';
      save_untouch_nii_hdr(hdr, fid);

      if ~isempty(ext)
         save_nii_ext(ext, fid);
      end

      fclose(fid);
      fid = fopen(sprintf('%s.img',fileprefix),'w');
   else
      fid = fopen(sprintf('%s.hdr',fileprefix),'w');

      if fid < 0,
         msg = sprintf('Cannot open file %s.hdr.',fileprefix);
         error(msg);
      end

      save_untouch0_nii_hdr(hdr, fid);

      fclose(fid);
      fid = fopen(sprintf('%s.img',fileprefix),'w');
   end

   ScanDim = double(hdr.dime.dim(5));		% t
   SliceDim = double(hdr.dime.dim(4));		% z
   RowDim   = double(hdr.dime.dim(3));		% y
   PixelDim = double(hdr.dime.dim(2));		% x
   SliceSz  = double(hdr.dime.pixdim(4));
   RowSz    = double(hdr.dime.pixdim(3));
   PixelSz  = double(hdr.dime.pixdim(2));

   x = 1:PixelDim;

   if filetype == 2 & isempty(ext)
      skip_bytes = double(hdr.dime.vox_offset) - 348;
   else
      skip_bytes = 0;
   end

   if double(hdr.dime.datatype) == 128

      %  RGB planes are expected to be in the 4th dimension of nii.img
      %
      if(size(nii.img,4)~=3)
         error(['The NII structure does not appear to have 3 RGB color planes in the 4th dimension']);
      end

      nii.img = permute(nii.img, [4 1 2 3 5 6 7 8]);
   end

   %  For complex float32 or complex float64, voxel values
   %  include [real, imag]
   %
   if hdr.dime.datatype == 32 | hdr.dime.datatype == 1792
      real_img = real(nii.img(:))';
      nii.img = imag(nii.img(:))';
      nii.img = [real_img; nii.img];
   end

   if skip_bytes
      fwrite(fid, zeros(1,skip_bytes), 'uint8');
   end

   fwrite(fid, nii.img, precision);
%   fwrite(fid, nii.img, precision, skip_bytes);        % error using skip
   fclose(fid);

   return;					% write_nii




%------------------------------------------------
function save_untouch_nii_hdr(hdr, fid)

   if ~isequal(hdr.hk.sizeof_hdr,348),
      error('hdr.hk.sizeof_hdr must be 348.');
   end

   write_header(hdr, fid);

   return;					% save_nii_hdr


%---------------------------------------------------------------------
function write_header(hdr, fid)

        %  Original header structures
	%  struct dsr				/* dsr = hdr */
	%       {
	%       struct header_key hk;            /*   0 +  40       */
	%       struct image_dimension dime;     /*  40 + 108       */
	%       struct data_history hist;        /* 148 + 200       */
	%       };                               /* total= 348 bytes*/

   header_key2(fid, hdr.hk);
   image_dimension2(fid, hdr.dime);
   data_history2(fid, hdr.hist);

   %  check the file size is 348 bytes
   %
   fbytes = ftell(fid);

   if ~isequal(fbytes,348),
      msg = sprintf('Header size is not 348 bytes.');
      warning(msg);
   end

   return;					% write_header


%---------------------------------------------------------------------
function header_key2(fid, hk)

   fseek(fid,0,'bof');

	%  Original header structures
	%  struct header_key                      /* header key      */
	%       {                                /* off + size      */
	%       int sizeof_hdr                   /*  0 +  4         */
	%       char data_type[10];              /*  4 + 10         */
	%       char db_name[18];                /* 14 + 18         */
	%       int extents;                     /* 32 +  4         */
	%       short int session_error;         /* 36 +  2         */
	%       char regular;                    /* 38 +  1         */
	%       char dim_info;   % char hkey_un0;        /* 39 +  1 */
	%       };                               /* total=40 bytes  */

   fwrite(fid, hk.sizeof_hdr(1),    'int32');	% must be 348.

   % data_type = sprintf('%-10s',hk.data_type);	% ensure it is 10 chars from left
   % fwrite(fid, data_type(1:10), 'uchar');
   pad = zeros(1, 10-length(hk.data_type));
   hk.data_type = [hk.data_type  char(pad)];
   fwrite(fid, hk.data_type(1:10), 'uchar');

   % db_name   = sprintf('%-18s', hk.db_name);	% ensure it is 18 chars from left
   % fwrite(fid, db_name(1:18), 'uchar');
   pad = zeros(1, 18-length(hk.db_name));
   hk.db_name = [hk.db_name  char(pad)];
   fwrite(fid, hk.db_name(1:18), 'uchar');

   fwrite(fid, hk.extents(1),       'int32');
   fwrite(fid, hk.session_error(1), 'int16');
   fwrite(fid, hk.regular(1),       'uchar');	% might be uint8

   % fwrite(fid, hk.hkey_un0(1),    'uchar');
   % fwrite(fid, hk.hkey_un0(1),    'uint8');
   fwrite(fid, hk.dim_info(1),      'uchar');

   return;					% header_key


%---------------------------------------------------------------------
function image_dimension2(fid, dime)

	%  Original header structures
	%  struct image_dimension
	%       {                                /* off + size      */
	%       short int dim[8];                /* 0 + 16          */
	%       float intent_p1;   % char vox_units[4];   /* 16 + 4       */
	%       float intent_p2;   % char cal_units[8];   /* 20 + 4       */
	%       float intent_p3;   % char cal_units[8];   /* 24 + 4       */
	%       short int intent_code;   % short int unused1;   /* 28 + 2 */
	%       short int datatype;              /* 30 + 2          */
	%       short int bitpix;                /* 32 + 2          */
	%       short int slice_start;   % short int dim_un0;   /* 34 + 2 */
	%       float pixdim[8];                 /* 36 + 32         */
	%			/*
	%				pixdim[] specifies the voxel dimensions:
	%				pixdim[1] - voxel width
	%				pixdim[2] - voxel height
	%				pixdim[3] - interslice distance
	%				pixdim[4] - volume timing, in msec
	%					..etc
	%			*/
	%       float vox_offset;                /* 68 + 4          */
	%       float scl_slope;   % float roi_scale;     /* 72 + 4 */
	%       float scl_inter;   % float funused1;      /* 76 + 4 */
	%       short slice_end;   % float funused2;      /* 80 + 2 */
	%       char slice_code;   % float funused2;      /* 82 + 1 */
	%       char xyzt_units;   % float funused2;      /* 83 + 1 */
	%       float cal_max;                   /* 84 + 4          */
	%       float cal_min;                   /* 88 + 4          */
	%       float slice_duration;   % int compressed; /* 92 + 4 */
	%       float toffset;   % int verified;          /* 96 + 4 */
	%       int glmax;                       /* 100 + 4         */
	%       int glmin;                       /* 104 + 4         */
	%       };                               /* total=108 bytes */

   fwrite(fid, dime.dim(1:8),        'int16');
   fwrite(fid, dime.intent_p1(1),  'float32');
   fwrite(fid, dime.intent_p2(1),  'float32');
   fwrite(fid, dime.intent_p3(1),  'float32');
   fwrite(fid, dime.intent_code(1),  'int16');
   fwrite(fid, dime.datatype(1),     'int16');
   fwrite(fid, dime.bitpix(1),       'int16');
   fwrite(fid, dime.slice_start(1),  'int16');
   fwrite(fid, dime.pixdim(1:8),   'float32');
   fwrite(fid, dime.vox_offset(1), 'float32');
   fwrite(fid, dime.scl_slope(1),  'float32');
   fwrite(fid, dime.scl_inter(1),  'float32');
   fwrite(fid, dime.slice_end(1),    'int16');
   fwrite(fid, dime.slice_code(1),   'uchar');
   fwrite(fid, dime.xyzt_units(1),   'uchar');
   fwrite(fid, dime.cal_max(1),    'float32');
   fwrite(fid, dime.cal_min(1),    'float32');
   fwrite(fid, dime.slice_duration(1), 'float32');
   fwrite(fid, dime.toffset(1),    'float32');
   fwrite(fid, dime.glmax(1),        'int32');
   fwrite(fid, dime.glmin(1),        'int32');

   return;					% image_dimension


%---------------------------------------------------------------------
function data_history2(fid, hist)

	% Original header structures
	%struct data_history
	%       {                                /* off + size      */
	%       char descrip[80];                /* 0 + 80          */
	%       char aux_file[24];               /* 80 + 24         */
	%       short int qform_code;            /* 104 + 2         */
	%       short int sform_code;            /* 106 + 2         */
	%       float quatern_b;                 /* 108 + 4         */
	%       float quatern_c;                 /* 112 + 4         */
	%       float quatern_d;                 /* 116 + 4         */
	%       float qoffset_x;                 /* 120 + 4         */
	%       float qoffset_y;                 /* 124 + 4         */
	%       float qoffset_z;                 /* 128 + 4         */
	%       float srow_x[4];                 /* 132 + 16        */
	%       float srow_y[4];                 /* 148 + 16        */
	%       float srow_z[4];                 /* 164 + 16        */
	%       char intent_name[16];            /* 180 + 16        */
	%       char magic[4];   % int smin;     /* 196 + 4         */
	%       };                               /* total=200 bytes */

   % descrip     = sprintf('%-80s', hist.descrip);     % 80 chars from left
   % fwrite(fid, descrip(1:80),    'uchar');
   pad = zeros(1, 80-length(hist.descrip));
   hist.descrip = [hist.descrip  char(pad)];
   fwrite(fid, hist.descrip(1:80), 'uchar');

   % aux_file    = sprintf('%-24s', hist.aux_file);    % 24 chars from left
   % fwrite(fid, aux_file(1:24),   'uchar');
   pad = zeros(1, 24-length(hist.aux_file));
   hist.aux_file = [hist.aux_file  char(pad)];
   fwrite(fid, hist.aux_file(1:24), 'uchar');

   fwrite(fid, hist.qform_code,    'int16');
   fwrite(fid, hist.sform_code,    'int16');
   fwrite(fid, hist.quatern_b,   'float32');
   fwrite(fid, hist.quatern_c,   'float32');
   fwrite(fid, hist.quatern_d,   'float32');
   fwrite(fid, hist.qoffset_x,   'float32');
   fwrite(fid, hist.qoffset_y,   'float32');
   fwrite(fid, hist.qoffset_z,   'float32');
   fwrite(fid, hist.srow_x(1:4), 'float32');
   fwrite(fid, hist.srow_y(1:4), 'float32');
   fwrite(fid, hist.srow_z(1:4), 'float32');

   % intent_name = sprintf('%-16s', hist.intent_name);	% 16 chars from left
   % fwrite(fid, intent_name(1:16),    'uchar');
   pad = zeros(1, 16-length(hist.intent_name));
   hist.intent_name = [hist.intent_name  char(pad)];
   fwrite(fid, hist.intent_name(1:16), 'uchar');

   % magic	= sprintf('%-4s', hist.magic);		% 4 chars from left
   % fwrite(fid, magic(1:4),           'uchar');
   pad = zeros(1, 4-length(hist.magic));
   hist.magic = [hist.magic  char(pad)];
   fwrite(fid, hist.magic(1:4),        'uchar');

   return;					% data_history

%------------------------------------------------
function [ext, esize_total] = verify_nii_ext(ext)

   if ~isfield(ext, 'section')
      error('Incorrect NIFTI header extension structure.');
   elseif ~isfield(ext, 'num_ext')
      ext.num_ext = length(ext.section);
   elseif ~isfield(ext, 'extension')
      ext.extension = [1 0 0 0];
   end

   esize_total = 0;

   for i=1:ext.num_ext
      if ~isfield(ext.section(i), 'ecode') | ~isfield(ext.section(i), 'edata')
         error('Incorrect NIFTI header extension structure.');
      end

      ext.section(i).esize = ceil((length(ext.section(i).edata)+8)/16)*16;
      ext.section(i).edata = ...
	[ext.section(i).edata ...
	 zeros(1,ext.section(i).esize-length(ext.section(i).edata)-8)];
      esize_total = esize_total + ext.section(i).esize;
   end

   return                                       % verify_nii_ext

%------------------------------------------------
function save_nii_ext(ext, fid)

   if ~exist('ext','var') | ~exist('fid','var')
      error('Usage: save_nii_ext(ext, fid)');
   end

   if ~isfield(ext,'extension') | ~isfield(ext,'section') | ~isfield(ext,'num_ext')
      error('Wrong header extension');
   end

   write_ext(ext, fid);

   return;                                      % save_nii_ext


%---------------------------------------------------------------------
function write_ext(ext, fid)

   fwrite(fid, ext.extension, 'uchar');

   for i=1:ext.num_ext
      fwrite(fid, ext.section(i).esize, 'int32');
      fwrite(fid, ext.section(i).ecode, 'int32');
      fwrite(fid, ext.section(i).edata, 'uchar');
   end

   return;                                      % write_ext

%------------------------------------------------
