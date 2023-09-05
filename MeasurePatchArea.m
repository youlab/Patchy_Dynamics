% Author: Emrah Simsek
% This code was used for patch area measurements.


% This code contains substantial adaptations 
% from previously unpublished code used for "Simsek et al. ISME (2022)",
% which was used to measure the bacterial colony sizes and fluorescence signals in the corresponding paper. 

% This code contains substantial adaptations 
% from "Karig et al. PNAS (2018), Turing_Spot_Analysis2.m",
% which was used to analyze the statistics of the spots from Turing
% pattern images in the corresponding paper. 

% Input: images from experiments
% Output: Statistics of the spots, including the max, min, mean, and the total
% number of the spots 
% Filter is applied to discard tiny spots below a threshold 

clear; 
close all;


filename = 'fileNAMEgoesHERE';
mkdir('Results')

%read the image
    FigX2 = imread(filename,'tif');
    
%image crop
%     I = imcrop(FigX2);
    I = FigX2;
    FigX2 = rgb2gray(I); 
    
%show the image and convert it to black-and-white with right threshold
    %figure(1),
    %subplot(1,2,1),
    figure(1),
    imshow(FigX2);  
    
    imwrite(FigX2,[filename '_greytest.tif'])
    
    moddivide= 30; bgval= mode(round(FigX2/moddivide), 'all')*moddivide; 
    
    Fig_bs = FigX2 - bgval;
    
    figure(2)
    imshow(Fig_bs);  

%     imwrite(Fig_bs, [filename '_bstest.tif'])
    
    % thresholding
    Th = 20; % threshold for imaging conversion (im2bw) to black-and-white. The threshold above is very critical.
    BW = Fig_bs > Th; 
    
    % Define a structuring element
    sep = 0; SE = strel('disk', sep, 4);

    % use MATLAB' s built-in imclose function for filling in the gaps as desired
    % using the structuting element created above
    BW = imclose(BW, SE);
    BW = imfill(BW,'holes');
    BW = imclearborder(BW);


    figure(3)
    imshow(BW)

%     imwrite(BW,[filename '_bwtest.tif'])

    % find the edges of the binary image
    BW_edge= edge(BW);

    % find the perimeter of the edges of the binary image
    EM_p = bwperim(BW_edge); 
    
    % for display purposes, overlay the phase contrast image with the mask
    % perimeter
    C = imfuse(FigX2, EM_p);

    % display the new binary mask created above corrected by removing small 
    % objects for the current inverted phase contrast image
    fig3= figure(4); 
    imshowpair(I, C, 'montage')
    title([sprintf(filename) '--Original cropped image (left), Masked image (right)'], 'Interpreter', 'None');
    set(gca,'FontSize', 18)
       
    Fig_bw = BW;

    %show the domains using pseudocolor image
    PseudoFig = label2rgb(Fig_bw, @spring, 'c', 'shuffle'); 
    figure(5)
    imshow(PseudoFig,'InitialMagnification', 'fit');  

    %Analyze domain characters
    %Count and label the domains
    [labeledImage, numObjects] = bwlabel(Fig_bw,8);
    %numObjects; % Count all distinct objects in the image.
    %disp(['original numObjects= ' num2str(numObjects)]);

    domaindata = regionprops(labeledImage,'basic'); %which returns: area, centroid, and 'bounding box'

    %filter small domains
    para_mini_ditected_area = 12; %minimal spot sizes to be detected (areas smaller than that will be filtered out)
    thres_area = para_mini_ditected_area;
%     idx = find([domaindata.Area] > thres_area );
%     Judge_domain = ismember(labeledImage, idx);

    % filter coalesced patches to be analyzed manually later
    para_max_ditected_area = 2000; %minimal spot sizes to be detected (areas smaller than that will be filtered out)
    up_thres_area = para_max_ditected_area;
    idx = find([domaindata.Area] > thres_area & [domaindata.Area] < up_thres_area);
    Judge_domain = ismember(labeledImage, idx);

    %show the filtered image
    Re_labeledImg = labeledImage .* Judge_domain ;
    
    fig4 = figure(6);
    imshow(Re_labeledImg,'InitialMagnification', 'fit');
%     saveas(fig4, 'Filtered_Binary.jpg');

    %Analyze the filtered image
    domaindataF = regionprops(Re_labeledImg,'basic'); %which returns: area, centroid, and 'bounding box'

    alldomainsF = [domaindataF.Area];
       
    %Collect all of the actual (non-zero)areas:
    alldomainsFR=alldomainsF(alldomainsF~=0); %get all the non zero values

    % convert from pixels to standard units
    %Nomalize the image according to the scale bar
    mm_per_pix = 0.085784314; %1.8872; 
    alldomainsFRS = alldomainsFR * mm_per_pix * mm_per_pix;

    
    
    num_growing_areaFRS = length(alldomainsFRS);  % Find area number for the filtered image.
    mean_growing_areaFRS = mean(alldomainsFRS);  % Find the mean area size for the filtered image.
    std_growing_areaFRS = std(alldomainsFRS);  % Find the std of area size for the filtered image.
    max_growing_areaFRS = max(alldomainsFRS);  % Find the max area size for the filtered image.
    min_growing_areaFRS = min(alldomainsFRS);  % Find the min area size for the filtered image.

%     % Here manually add the number of empty patches
%     num_vertices = 593;
%     
%     num_empty = 593-286; %num_vertices - num_growing_areaFRS;
%     
%     alldomainsFRS_corrected = [alldomainsFRS zeros(1, num_empty)];
%     
%     disp(['Filtered, Corrected, Real Scaled']);
%     num_areaFRS=length(alldomainsFRS_corrected);  % Find area number for the filtered image.
%     disp(['num_area(filRS)= ' num2str(num_areaFRS)]);
%     mean_areaFRS=mean(alldomainsFRS_corrected);  % Find the mean area size for the filtered image.
%     disp(['mean_area(filRS)= ' num2str(mean_areaFRS)]);
%     max_areaFRS = max(alldomainsFRS_corrected);  % Find the maximal area in the filtered image.
%     disp(['max_area(filRS)= ' num2str(max_areaFRS)]);
%     min_areaFRS = min(alldomainsFRS_corrected);  % Find the minimal area in the filtered image.
%     disp(['min_area(filRS)= ' num2str(min_areaFRS)]);
%     med_areaFRS = median(alldomainsFRS_corrected);  % Find the median area in the filtered image.
%     disp(['med_area(filRS)= ' num2str(med_areaFRS)]);
%     std_areaFRS = std(alldomainsFRS_corrected);  % Find the standard deviation area in the filtered image.
%     disp(['std_area(filRS)= ' num2str(std_areaFRS)]);
%     TotalArea = size(FigX2, 1) * size(FigX2, 2) * mm_per_pix*mm_per_pix; % Find the total area analyzed in square milimeters
%     disp(['TotalArea(mm^2)= ' num2str(TotalArea)]);
%     SpotDensity = num_areaFRS/TotalArea; % Find the spot density (i.e., number of spots per unit area analyzed in 1/square milimeters
%     disp(['SpotDenst(1/mm^2) ' num2str(SpotDensity)]);
%     CV_areaFRS = std_areaFRS/mean_areaFRS;  % Find the coefficient of variation of patch area in the filtered image.
%     disp(['CV_area(filRS)= ' num2str(CV_areaFRS)]);
%     
    bw2 = bwperim(Re_labeledImg);
    
%     bw3 = imdilate(bw2, strel('disk', 1));
    
    M = imoverlay(FigX2, bw2, 'red');  
    
    fig5 = figure(7);
    imshow(M);
    
%     imwrite(fig5, [filename '_MaskMontage.jpg'])
    saveas(fig5, [filename '_MaskMontage.jpg']);
    
%     output = [num_vertices; num_growing_areaFRS; 100*(num_vertices-num_empty)/num_vertices; mean_growing_areaFRS; std_growing_areaFRS; std_growing_areaFRS/mean_growing_areaFRS; max_growing_areaFRS; min_growing_areaFRS; mean_areaFRS; std_areaFRS; std_areaFRS/mean_areaFRS; sum(alldomainsFRS)];

    [ycdf,xcdf] = cdfcalc(alldomainsFRS);
    xccdf = xcdf;
    yccdf = 1-ycdf(1:end-1);
    
%     fig6 = figure(8);
%     subplot(2,2,1)
%         % imshow(Re_labeledImg,'InitialMagnification', 'fit');
%         montage({FigX2, Re_labeledImg}); 
%         title('Left: Original image, Right: Processed binary image')
%         set(gca, 'FontSize', 24)
%     subplot(2,2,2)
%         histogram(alldomainsFRS)
%         xlabel('Patch area, a (mm^2)')
%         ylabel('Count')
%         title('Histogram')
%         set(gca, 'FontSize', 24)
%     subplot(2,2,3)
%         semilogy(xcdf, yccdf, 'bo', 'LineWidth', 3);
%         xlabel('Patch area, a (mm^2)')
%         ylabel('P(A>a)')
%         title('Semilog scale')
%         set(gca, 'FontSize', 24)
%     subplot(2,2,4)
%         loglog(xcdf, yccdf, 'ro', 'LineWidth', 3);
%         xlabel('Patch area, a (mm^2)')
%         ylabel('P(A>a)')
%         title('log-log scale')
%         set(gca, 'FontSize', 24)
% %     saveas(fig6, 'PatchSizeDistribution.png')

XY =[domaindataF.Centroid];

Xpos = XY(1:2:end)';
Ypos = XY(2:2:end)';

A = (mm_per_pix^2)*[domaindataF.Area]'; % patch areas is mm^2

results = [Xpos, Ypos, [domaindataF.Area]', A];
results(any(isnan(results),2),:) = [];

stat = [length(results(:,4)), mean(results(:,4)), max(results(:,4)), min(results(:,4)), std(results(:,4)), std(results(:,4))/mean(results(:,4))];


colNames1 = {'X', 'Y', 'Area_pixels', 'Area_mm2'};
colNames2 = {'NumberOfPatches', 'MeanPatchArea', 'MaxPatchArea', 'MinPatchArea', 'SDPatchArea', 'CVPatchArea'};

sTable1 = array2table(results, 'VariableNames', colNames1);
sTable2 = array2table(stat, 'VariableNames', colNames2);

writetable(sTable1, [pwd '/Results/' filename '_ResultsToBeManuallyFinalized.csv'])
writetable(sTable2, [pwd '/Results/' filename '_ConvertedStatsToBeManuallyFinalized.csv'])

% writetable(sTable1, [pwd '/Results/' filename '_ResultsFinal.csv'])
% writetable(sTable2, [pwd '/Results/' filename '_ConvertedStatsFinal.csv'])