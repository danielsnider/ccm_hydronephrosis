set(0,'DefaultFigureWindowStyle','docked');
addpath(genpath('functions'));

base_path = 'C:\Users\danie\Documents\CCM\LaurenErdman\Left-Right-Sagittal-Transverse-Only';
folders = dir(base_path);
ResultTable = table();


plugin_name = 'Hydronepherosis';
thresh_param = Inf;
thresh_param = 60;
thresh_param = 40;
debug_level = 'All';
debug_level = 'Result Only';
threshold_smooth_param = 7;
watershed_smooth_param = 15;
seeds = false;
boarder_clear = false;
boarder_clear = 0;
intensity_filter_param = 19;
% intensity_filter_param = 30;
max_area = Inf;
min_area = 100;

save_figs = true;
save_mag = '-native';

date_str = datestr(now,'yyyymmddTHHMMSS');
% save_path_prefix = sprintf('plots/%s_%s/',date_str, strrep(folder_name,'\','__'));
save_path_prefix = sprintf('%s/',date_str);
mkdir(['plots/' save_path_prefix]);


for folder = folders'
  if ~folder.isdir | strcmp(folder.name,'.') | strcmp(folder.name,'..')
    continue
  end

  % folder_name = 'D2219006_3\1.2.840.113663.1500.1.193471476.2.1.20061211.141402.453'
  % full_path = fullfile(base_path, folder_name);
  folder_name = folder.name;
  full_path = fullfile(folder.folder, folder.name);
  files = dir(sprintf('%s/*.dcm', full_path));
  if isempty(files)
    % instead of finding images, we found a folder in a folder, now lets look inside that folder
    files = dir(full_path);
    for file = files'  % should be '.', '..', and a directory containing images
      if file.isdir
        folder = file;
      end
    end
  end
  full_path = fullfile(folder.folder, folder.name);
  files = dir(sprintf('%s/*.dcm', full_path));

  for file=files'
    known_types = {'rt-sag', 'rt-trans', 'lt-sag', 'lt-trans'};
    contains_known_type = [contains(file.name,'rt-sag'), contains(file.name,'rt-trans'), contains(file.name,'lt-sag'), contains(file.name,'lt-trans')];
    if sum(contains_known_type)==0
      continue  % only process known image types
    end
    kidney_type = known_types{find(contains_known_type)};

    close all;
    file_path = fullfile(file.folder, file.name)
    orig_img = dicomread(file_path);
    img = orig_img;
    img_size_y = size(img,1);
    img_size_x = size(img,2);

    % Skip if three dimensions because that means a weird overlay is happening
    if ndims(orig_img)==3
      %warning(sprintf('Skipping image because it has a color overlay: %s', file.name));
      %continue;
      img=rgb2gray(img);
    end

    figure;imshow(img);

    % Reduce effect of crosshairs (reducing to more middle intensity)
    img(img==226) = 128;

    % Crop
%     img = img(128:128+580,173:173+709);
%     if ismember(debug_level,{'All'})
%       f = figure(886); clf; set(f,'name','crop','NumberTitle', 'off');
%       imshow(img,[]);
%     end

    % Smooth
    img_smooth = imgaussfilt(img,threshold_smooth_param);
    if ismember(debug_level,{'All'})
      f = figure(886); clf; set(f,'name','smooth for threshold','NumberTitle', 'off');
      imshow(img_smooth,[]);
    end


    % threshold
    img_thresh = img_smooth < thresh_param;
    if ismember(debug_level,{'All'})
      f = figure(885); clf; set(f,'name','threshold','NumberTitle', 'off');
      imshow(img_thresh,[]);
    end

    % remove seeds outside of our img mask
    if ~isequal(seeds,false)
      seeds(img_thresh==0)=0;

      % Debug with plot
      if ismember(debug_level,{'All'})
        [X Y] = find(seeds);
        f = figure(826); clf; set(f,'name','input seeds','NumberTitle', 'off')
        imshow(img,[]);
        hold on;
        plot(Y,X,'or','markersize',2,'markerfacecolor','r')
      end
    end

    if isequal(watershed_smooth_param,false)
      img_smooth2 = img;
    else
      img_smooth2 = imgaussfilt(img,watershed_smooth_param);
    end
    if ismember(debug_level,{'All'})
      f = figure(889); clf; set(f,'name','smooth for watershed','NumberTitle', 'off');
      imshow(img_smooth2,[]);
    end


    %% Watershed
    if ~isequal(seeds,false)
      img_min = imimposemin(max(img_smooth2(:))-img_smooth2,seeds); % set locations of seeds to be -Inf as per  matlab's watershed
    else
      % img_min = max(img_smooth2(:))-img_smooth2;
      img_min = img_smooth2;
    end

    if ismember(debug_level,{'All'})
      f = figure(564); clf; set(f,'name','imimposemin','NumberTitle', 'off')
      imshow(img_min,[]);
    end

    img_ws = watershed(img_min);
    if ismember(debug_level,{'All'})
      f = figure(562); clf; set(f,'name','watershed','NumberTitle', 'off')
      imshow(img_ws,[]);
    end

    img_ws(img_thresh==0)=0; % remove areas that aren't in our img mask
    if ismember(debug_level,{'All'})
      f = figure(561); clf; set(f,'name','watershed & threshold','NumberTitle', 'off')
      imshow(img_ws,[]);
    end

    % Clear cells touching the boarder
    if isnumeric(boarder_clear)
      if isequal(boarder_clear,0)
        bordercleared_img = imclearborder(img_ws);
      else
        % Clear objects touching boarder too much (ex. too much could be 1/4 of perimeter)
        bordercleared_img = img_ws;
        for idx=1:max(img_ws(:))
          single_object = img_ws == idx;
          [x y] = find(single_object);
          count_edge_touches = ismember([x; y], [1 size(single_object,1), size(single_object,2)]);
          count_perim = bwperim(single_object);
          % If the object touches the edge for more than 1/5 the length of the perimeter, delete it
          if sum(count_edge_touches) > sum(count_perim(:)) / (100 / boarder_clear)
            bordercleared_img(bordercleared_img==idx)=0; % delete this object
          end
        end
      end
      if ismember(debug_level,{'All'})
        f = figure(511); clf; set(f,'name','imclearborder','NumberTitle', 'off')
        imshow(bordercleared_img,[]);
      end
    else
      bordercleared_img = img_ws;
    end

    % Fill holes
    filled_img = imfill(bordercleared_img,'holes');
    if ismember(debug_level,{'All'})
      f = figure(512); clf; set(f,'name','imfill','NumberTitle', 'off')
      imshow(filled_img,[]);
    end

    % Remove objects that are too small or too large
    filled_img = bwlabel(filled_img); stats=regionprops('table', filled_img, img, 'MeanIntensity','Centroid','Area');
    size_filtered_im = filled_img;
    reject_idx = find(stats.Area > max_area | stats.Area < min_area);
    size_filtered_im(ismember(size_filtered_im,reject_idx))=0;
    if ismember(debug_level,{'All'})
      f = figure(5100); clf; set(f,'name','size_filt','NumberTitle', 'off')
      imshow(size_filtered_im,[]);
    end

    % Keep only central centroids
    centroid_filtered_im = size_filtered_im;
    if ~isempty(stats)
      reject_idx = find(stats.Centroid(:,2) < img_size_y * 0.25 | stats.Centroid(:,2) > img_size_y * 0.65 | stats.Centroid(:,1) < img_size_x * 0.29 | stats.Centroid(:,1) > img_size_x * 0.65 ); % 29 becouse there is info on the right, 25 because there is a lot of detail at the top. Any larger numbers and too many objects would be filtered
      centroid_filtered_im(ismember(centroid_filtered_im,reject_idx))=0;
      if ismember(debug_level,{'All'})
        f = figure(5012); clf; set(f,'name','centroid_filt','NumberTitle', 'off')
        imshow(centroid_filtered_im,[]);
      end
    end

    % Keep only dark blobs
    intensity_filtered_im = centroid_filtered_im;
    int_reject_idx = find(stats.MeanIntensity>intensity_filter_param);
    intensity_filtered_im(ismember(intensity_filtered_im,int_reject_idx))=0;
    if ismember(debug_level,{'All'})
      f = figure(5192); clf; set(f,'name','int_filt','NumberTitle', 'off')
      imshow(intensity_filtered_im,[]);
    end

    % Label result
    labelled_img_all = bwlabel(centroid_filtered_im);
    labelled_img_disease = bwlabel(intensity_filtered_im);

    % Visualization
    if ismember(debug_level,{'All','Result Only','Result With Seeds'})
      f = figure(743); clf; set(f,'name',[plugin_name ' Result'],'NumberTitle', 'off')
      % Display original image
      % Cast img as double, had issues with 32bit
      img8 =  orig_img; %im2uint8(double(img));
      if min(img8(:)) < prctile(img8(:),99.5)
          min_max = [min(img8(:)) prctile(img8(:),99.5)];
      else
          min_max = [];
      end
      imshow(img8,[min_max]);
      hold on
      % Display color overlay (grey all objects)
      labelled_perim = imdilate(bwperim(labelled_img_all),strel('disk',1));
      labelled_rgb = label2rgb(uint32(labelled_perim), [.7 .7 .7], [1 1 1]);
      himage = imshow(im2uint8(labelled_rgb),[min_max]);
      himage.AlphaData = labelled_perim*1;
      % Display color overlay (colored possible hydronephrosis)
      labelled_perim = imdilate(bwlabel(bwperim(labelled_img_disease)),strel('disk',1));
      labelled_rgb = label2rgb(uint32(labelled_perim), 'jet', [1 1 1], 'shuffle');
      himage = imshow(im2uint8(labelled_rgb),[min_max]);
      himage.AlphaData = labelled_perim*1;
      if ismember(debug_level,{'All','Result With Seeds'})
        if ~isequal(seeds,false)
          seeds(aligned_labelled_img<1)=0;
          % Display red dots for seeds
          [xm,ym]=find(seeds);
          hold on
          plot(ym,xm,'or','markersize',2,'markerfacecolor','r','markeredgecolor','r')
        end
      end
      hold off
      if save_figs
        pause(0.1)
        fig_name = sprintf('plots/%s/%s_%s.png',save_path_prefix,folder_name,file.name);
        export_fig(fig_name,save_mag);
      end
    end

    % Save Stats (Metadata)
    ThisResult = table();
    ThisResult.PatientMDR = {folder_name};
    ThisResult.KidneyType = {kidney_type};
    % Save Stats (All)
    ThisResult.NumObjectsAll = max(max(labelled_img_all));
    ThisResult.TotalAreaAll = sum(sum(labelled_img_all>0));
    ThisResult.MeanIntensityAll = mean(img(labelled_img_all>0));
    % Save Stats (Desease)
    ThisResult.NumObjectsSuspect = max(max(labelled_img_disease));
    ThisResult.TotalAreaSuspect = sum(sum(labelled_img_disease>0));
    ThisResult.MeanIntensitySuspect = mean(img(labelled_img_disease>0));
    % Save Stats (Per Object)
    if ThisResult.NumObjectsAll > 0
      stats = regionprops('table',labelled_img_all,img,'Area','Centroid','Eccentricity','MajorAxisLength','MinorAxisLength','Orientation','Solidity','MeanIntensity');
      ThisResult.Areas = {regexprep(num2str(stats.Area'),'\s+',';')};
      ThisResult.MeanIntensities = {regexprep(num2str(stats.MeanIntensity'),'\s+',';')};
      ThisResult.CentroidX = {regexprep(num2str(stats.Centroid(:,1)'),'\s+',';')};
      ThisResult.CentroidY = {regexprep(num2str(stats.Centroid(:,2)'),'\s+',';')};
      ThisResult.MajorAxisLength = {regexprep(num2str(stats.MajorAxisLength'),'\s+',';')};
      ThisResult.MinorAxisLength = {regexprep(num2str(stats.MinorAxisLength'),'\s+',';')};
      ThisResult.Orientation = {regexprep(num2str(stats.Orientation'),'\s+',';')};
      ThisResult.Solidity = {regexprep(num2str(stats.Solidity'),'\s+',';')};
      ThisResult.Eccentricity = {regexprep(num2str(stats.Eccentricity'),'\s+',';')};
    else
      ThisResult.Areas = {NaN};
      ThisResult.MeanIntensities = {NaN};
      ThisResult.CentroidX = {NaN};
      ThisResult.CentroidY = {NaN};
      ThisResult.MajorAxisLength = {NaN};
      ThisResult.MinorAxisLength = {NaN};
      ThisResult.Orientation = {NaN};
      ThisResult.Solidity = {NaN};
      ThisResult.Eccentricity = {NaN};
    end

    % Store
    ResultTable(height(ResultTable)+1,:) = ThisResult;
       
    break
  end
end
% Save to disk

save_name = sprintf('plots/%s/results.csv',save_path_prefix);
writetable(ResultTable,save_name);

msgbox('Done!')
