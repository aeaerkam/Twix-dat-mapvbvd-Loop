%% mapVBVD - a short Tutorial
% Philipp Ehses <philipp.ehses@tuebingen.mpg.de>
% Created: September 2016

% loop

%% Parse twix file
% The basic idea of mapVBVD is to divide reading the data into two parts:
% First, parsing the twix file and saving all relevant information from the
% 'measurement data headers' (mdh) for each acquired line of k-space 
% together with their position in the file (byte position).
% Second, read the data into memory after choosing which data should be 
% read (which slice/s, which repetition/s, which line/s, ...). This is
% similar to memory-mapping, hence the name.
%
% Parsing of the data can be achieved in several ways:

% 1) open dialog to choose file from
% twix = mapVBVD();

% 2) use full filename (including path if necessary)
% twix = mapVBVD('meas_MID00501_FID24867_t2_tse_tra_192_4mm_2ave.dat');

% 3) use measurement id (you need to 'cd' to the containing folder)
% close all;
% clear all;

 % modded 2023 
addpath ..
names=cell(1,10e3);
imgs = cell(size(names)); 
twixes=cell(size(names));
addpath(genpath('..'))


 plotFig=0;  % Flag to plot all intermediate figures
 


if ~exist("firstRun")
firstRun = 1;
end

if firstRun
 MeasurementList = {};
 firstRun = 0;
end

% cd ("dat23May\")   

for MID = [721]
% for MID=1457:1462 
% for MID=1610:1612  % Faulty distorted supply modulation block, pulseq TR=100
% for MID=1383:1394

%  for MID=[1383:1394]
%  for MID=[2145:2146]  %5V 10V localizer New Coil 
 %  for MID=[1385 1392 1394]
% for MID=[7:10]
%  for MID=[2068:2125]
%  for MID=[976:1024]  March 9  

% for MID = [721]
%  for MID = [1024]

% for MID=1386 % 15V 4ms
%     for MID=1391 % 15V 2ms
try
  twix = mapVBVD(MID);  % try catch in case of any errors
catch ME
%    if (strcmp(ME.identifier,'MATLAB:catenate:dimensionMismatch'))
      msg = ['Error at MID: ' num2str(MID)];
        causeException = MException('MATLAB:myCode:dimensions',msg);
        ME = addCause(ME,causeException);
		disp(msg);
%    end
%    rethrow(ME)
end 


% filename=MID;
% 
% if ~exist('filename','var') || isempty(filename)
%     info = 'Please select binary file to read';
%     [fname,pathname]=uigetfile('*.dat',info);
%     if isempty(pathname)
%         return
%     end
%     filename=[pathname fname];
% else
%     if ischar(filename)
%         % assume that complete path is given
%         if  ~strcmpi(filename(end-3:end),'.dat');
%             filename=[filename '.dat'];   %% adds filetype ending to file
%         end
%     else
%         % filename not a string, so assume that it is the MeasID
%         measID   = filename;
%         filelist = dir('*.dat');
%         filesfound = 0;
%         for k=1:numel(filelist)
%             if regexp(filelist(k).name,['^meas_MID0*' num2str(measID) '_.*\.dat'])==1
%                 if filesfound == 0
%                     filename = filelist(k).name;
%                 end
%                 filesfound = filesfound+1;
%             end
%         end
%         if filesfound == 0
%             error(['File with meas. id ' num2str(measID) ' not found.']);
%         elseif filesfound > 1
%             disp(['Multiple files with meas. id ' num2str(measID) ' found. Choosing first occurence.']);
%         end
%     end
% end

% % add absolute path, when no path is given
% [pathstr, name, ext] = fileparts(filename);
% 
% if isempty(pathstr)
%     pathstr  = pwd;
%     filename = fullfile(pathstr, [name ext]);
% end



%% First look at 'twix' 
% There are two cells in twix (this is a multi-raid VD/VE file).

 twix;
  

%%
% The first 'measurement' includes the noise data and its header, the 
% second holds the image data (as well as other data such as phasecor,
% refscan,...).
% Note that twix data from VA and VB systems will never show multiple
% measurements and 'twix' will simply be a struct with the contents of
% twix{2}.
% 
% twix{:}

%%
% In this tutorial we do not care for noise data, so let's get rid of the
% first measurement:
% 
% twix = twix{2};


%% A closer look
% Now let's look at the object that holds the image information
twix.image

%%
% As you can see there is a lot of information, it seems that we acquired
% three slices, two averages and 11 segments. Here's the matrix size that
% would be created if we would call twix.image():
twix.image.dataSize(1:11)

%% Making use of flags
% Let's get rid of the average dimension by telling the twix object to
% automatically average them during the read procedure. In addition, we
% also want to get rid of the oversampling that is almost always performed
% in read direction (2 x oversampling = 2 x larger FoV in read):
twix.image.flagDoAverage = false;
twix.image.flagRemoveOS  = false;

%%
% Btw, flagDoAverage is a shortcut to flagAverageDim, and in fact we can
% average every of the 16 dimensions by setting the corresponding entry in
% flagAverageDim to true;
%
% Now let's see how the two flags changed the 'size' of our data:

twix.image.dataSize(1:11);

%%
% The average dimension is gone and the 'Col'(read)-dimension is only half
% its previous size. This will greatly reduce our memory footprint...
%
% There are a bunch of our flags - most of them are just shortcuts to
% flagAverageDim. 'flagRampSampRegrid' may be of interest for people
% working with EPI data as it tells the twix-object to automatically regrid
% ramp-sampled data during the read procedure.

%% First look at our data
% Now, let's look at different segments of our data. What are segments?
% Segments indicate how the acquisition of the data was 'segmented', e.g.
% in an EPI or TSE acquisition. This data is from a TSE and the number of
% segments is equal to the turbo factor. So each of this segments
% corresponds to a different echo in our TSE echo train.

%%
% Let's read in the second slice of our data. There are ways to shorten the
% following expression, but this is probably the safest way:

%                 Col Cha Lin Par Sli Ave Phs Eco Rep Set Seg 
% data = twix.image(  :,  :,  :,  1,  2,  1,  1,  1,  1,  1, :);
% 
% %%
% % Now we look at four segments for one coil:
% subplot(2,2,1);
% imagesc(abs(squeeze(data(:,1,:, 1))).^0.2);
% subplot(2,2,2);
% imagesc(abs(squeeze(data(:,1,:, 4))).^0.2);
% subplot(2,2,3);
% imagesc(abs(squeeze(data(:,1,:, 8))).^0.2);
% subplot(2,2,4);
% imagesc(abs(squeeze(data(:,1,:,11))).^0.2);

%% 
% As you can see, each segment holds a continuous 'segment' of k-space.
% There are lots of zeros in the data, so our data matrix is much larger
% than it needs to be. The segment information can be important for phase
% correction but let's omit it in order to reduce the size of our data:

twix.image.flagIgnoreSeg = true;

% Let's look again at the data size to see that the segments dim. is gone:
twix.image.dataSize(1:11);


%%
% By the way, I got quite a few emails from people that had problems with
% their data not fitting into their computer's memory. This flag as well as
% the other averaging flag can help a lot!
%
% Let's read in the data again. This time we can read in all slices at
% once, our matrix is much smaller now:
data = twix.image();
size(data)

%% 
% And here's the k-space for the first coil and the first slice:
if plotFig==1
figure(MID)
ax(1)=subplot(1,2,1); 
% imagesc(abs(squeeze(data(:,1,:,1))).^0.2); colorbar;
imagesc(abs(squeeze(data(:,1,:,1)))); colorbar; title(['Linear scale k-space Data of MID: ' num2str(MID)])

figure(MID+1000)
axx(1)=subplot(1,2,1); 
% imagesc(abs(squeeze(data(:,1,:,1))).^0.2); colorbar;
imagesc(angle(squeeze(data(:,1,:,1)))); colorbar; title(['Phase of k-space Data of MID: ' num2str(MID)])
else
end
%%
%  data = squeeze(twix.image(:,:,:,1,sli));
% test=(ifftshift(ifft(fftshift(data,f),[],f),f));
% figure
%  imagesc(abs(test))
%% Let's reconstruct the data
% Now, we're ready to reconstruct the data. For large datasets it can
% be very important to never read in all the data into memory but only read
% in as much data as is required for each reconstruction step. In this 
% example, we can reconstruct each slice separately in order to reduce our
% memory footprint.

%
% The final matrix size of our image (single precision is enough):
os  = 1; % oversampling factor
img = zeros([twix.image.NCol/os, twix.image.NLin, twix.image.NSli], 'single');
imgph = zeros([twix.image.NCol/os, twix.image.NLin, twix.image.NSli], 'single');
for sli = 1:twix.image.NSli
    % read in the data for slice 'sli'
    data = twix.image(:,:,:,1,sli);
    
    % fft in col and lin dimension:
    fft_dims = [1 3];
    for f = fft_dims
        data = ifftshift(ifft(fftshift(data,f),[],f),f);
    end
    
    % sum-of-square coil combination:
    img(:,:,sli) = squeeze(sqrt(sum(abs(data).^2,2)));
    imgph(:,:,sli) = squeeze(angle(sum((data),2)));
end
imgs{MID}=img;
twixes{MID}=twix;

%%
% And here's the result for the three slices:

scrsz = get(groot,'ScreenSize');
if plotFig==1
figure(MID)
ax(2)=subplot(1,2,2);
% figure('Position',[1 scrsz(4)/4 scrsz(3)/2 scrsz(4)/4])
% imagesc(flipud(img(:,:)), [0, 0.7*max(img(:))]), axis image off; colorbar;
imagesc(flipud(img(:,:))), axis image off; colorbar;
colormap gray; %caxis([0 3e-7])
title(twix.hdr.Dicom.tProtocolName,'Interpreter', 'none')
else
end
%% !!!!!!  COMMENT THIS IF YOU WANT QUICK PLOTTING WITHOUT INTERRUPTIONS
% uicontrol('Style', 'pushbutton', 'Callback', 'uiresume(gcbf)','String','Continue');
% title('Move and scale the ellipse. Press Continue when you are ready...');
% h = imellipse(gca,[20 20 10 10]);
% uiwait(gcf);
% setResizable(h,0);
% BW = createMask(h);
% disp(['THE AVERAGE VALUE on your ELLIPSE IS::: '  num2str(sum(sum(img(BW)))/length(BW))])
%%
%phase
if plotFig==1
figure(MID+1000)
axx(2)=subplot(1,2,2); 
imagesc(flipud(imgph(:,:))), axis image off; colorbar;
title(twix.hdr.Dicom.tProtocolName,'Interpreter', 'none')
else
end
 names{MID}=twix.hdr.Dicom.tProtocolName;
% colormap(ax(2),gray)
%%
%
%

%% Appendix: Header information
% In case you are looking for header information, the twix-structure also
% contains a 'hdr' entry. It contains all the header information that was
% recognized by read_twix_hdr.m. It's a very basic script, so unfortunately
% some header information is lost.

twix.hdr

%%
% There are different types of headers with a lot of redundant data. I
% usually use the 'MeasYaps' structure which includes most of the stuff
% that you can access in IDEA via 'MrProt' (also, same structure):
%
% A few examples:

TR = twix.hdr.MeasYaps.alTR{1} ; %us
TE = twix.hdr.MeasYaps.alTE{1}  ;%us

% sli_thickness = twix.hdr.MeasYaps.sSliceArray.asSlice{1}.dThickness
% sli_position  = twix.hdr.MeasYaps.sSliceArray.asSlice{1}.sPosition
% 
% turbo_factor = twix.hdr.MeasYaps.sFastImaging.lTurboFactor
% 

figure(128)
title('middle cutlines')
if MID==721
    cut = 64
else
    cut = 128
end

aa=flipud(img(:,:));
MidLine=aa(cut,:);

meanbackground = mean(MidLine(1:43)) %  
peakval = max(MidLine)
SNRloose = peakval/meanbackground

o.MID = MID;
o.name = twix.hdr.Dicom.tProtocolName;
o.meanbackground = meanbackground;
o.SNRloose = SNRloose;
o.MidLine = MidLine;
o.twix = twix;
o.img = img;
o.imgph = imgph;
o.sli = sli;
o.aa = aa;
o.data = reshape(data,[size(data,1) size(data,3)]);
plot(MidLine, 'LineWidth',1.5, 'DisplayName',twix.hdr.Dicom.tProtocolName)
hold on; legend(gca,'show','Interpreter', 'none')
grid on; grid minor;

MeasurementList{end+1} =  o


end




 


% arrangefigures()





