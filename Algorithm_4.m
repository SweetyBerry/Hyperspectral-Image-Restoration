clear all; close all; clc
%=================================================================================
% !!!!!!!!!!!!!!!    Stupid Algorithm Warning    !!!!!!!!!!!!!
%       This algorithm may take you about 7 to 10 minutes.
% You may try grab some food or check your instagram while running this program.
%=================================================================================
%% parameter setting
N=3;

%% Load Data
load Nevada.mat

%% Generate True A & S (using Nevada.mat)
[A_tru, S_tru, time_tru] = HyperCSI(reshape(X,150*150,183)',N);
X_reconstructed_true = A_tru*S_tru;
S_tru_3d = reshape(S_tru',[150,150,N]);

%% Create mask
corrupted_band_size = [2,1,1,1,2,1,1,1,1,2]*2; %stripes size
corrupted_band = [10,24,40,56,73,80,88,125,132,140];  %location where the stripes start
corrupted_band_end = corrupted_band + corrupted_band_size - 1 ; %location where the stripes end


%% Add mask && Cropping
X_corrupted = X;
X_croped = X;

% Graph with black stripes
for i=1:size(corrupted_band,2)
    X_corrupted(:,corrupted_band(i):corrupted_band_end(i),:) = 0;
end

% Graph without black stripes but shrinked width
X_croped(:,[corrupted_band(1):corrupted_band_end(1),corrupted_band(2):corrupted_band_end(2), ...
    corrupted_band(3):corrupted_band_end(3),corrupted_band(4):corrupted_band_end(4), ...
    corrupted_band(5):corrupted_band_end(5),corrupted_band(6):corrupted_band_end(6), ...
    corrupted_band(7):corrupted_band_end(7),corrupted_band(8):corrupted_band_end(8) ...
    corrupted_band(9):corrupted_band_end(9),corrupted_band(10):corrupted_band_end(10)],:) = [];

width = 150-sum(corrupted_band_size);
X_reshaped = reshape(X_croped,[width*150,183])';

%% Run HyperCSI
tic
[A_est, S_est, time] = HyperCSI(X_reshaped,N);
X_reconstructed_damaged = A_est*S_est;

%% Reconstruct(extend) S
S_3d = reshape(S_est',[150,width,N]);
S_3d_extend = zeros(150,150,N);

%---------------------------------Marking damaged section---------------------------------
%Caculate sections size
section_size = zeros(1,size(corrupted_band_size,2)+1);
for i=1:size(corrupted_band_size,2)+1
   if i==1
       section_size(1) = corrupted_band(1) - 1;
   elseif i==size(corrupted_band_size,2)+1
       section_size(i) = 150-corrupted_band_end(i-1);
   else
       section_size(i) = corrupted_band(i) - corrupted_band_end(i-1) - 1;
   end
end

%---------------------------------Adding Zeros(Stripes) back to S---------------------------------
pointer = 1;
for i=1:size(corrupted_band,2)+1
   if i==1
       S_3d_extend(:,1:(corrupted_band(i)-1),:) =  S_3d(:,1:section_size(i),:);
   elseif i==size(corrupted_band,2)+1
       S_3d_extend(:,(corrupted_band_end(i-1)+1):150,:) =  S_3d(:,pointer:pointer+section_size(i)-1,:);
   else
       S_3d_extend(:,(corrupted_band_end(i-1)+1):corrupted_band(i)-1,:) =  S_3d(:,pointer:pointer+section_size(i)-1,:);
   end
   pointer = pointer + section_size(i);
end

%---------------------------------Marking location where Algorithem should restore---------------------------------
even_or_odd = mod( corrupted_band_size , 2 );
times = floor(corrupted_band_size/2);
start = zeros(2,2,2);

start_l_1 = (corrupted_band(1)-1:corrupted_band(1)+corrupted_band_size(1)/2-2);
start_r_1 = (corrupted_band(1)+corrupted_band_size(1)/2:corrupted_band_end(1));
start_l_2 = (corrupted_band(2)-1:corrupted_band(2)+corrupted_band_size(2)/2-2);
start_r_2 = (corrupted_band(2)+corrupted_band_size(2)/2:corrupted_band_end(2));
start_l_3 = (corrupted_band(3)-1:corrupted_band(3)+corrupted_band_size(3)/2-2);
start_r_3 = (corrupted_band(3)+corrupted_band_size(3)/2:corrupted_band_end(3));
start_l_4 = (corrupted_band(4)-1:corrupted_band(4)+corrupted_band_size(4)/2-2);
start_r_4 = (corrupted_band(4)+corrupted_band_size(4)/2:corrupted_band_end(4));
start_l_5 = (corrupted_band(5)-1:corrupted_band(5)+corrupted_band_size(5)/2-2);
start_r_5 = (corrupted_band(5)+corrupted_band_size(5)/2:corrupted_band_end(5));
start_l_6 = (corrupted_band(6)-1:corrupted_band(6)+corrupted_band_size(6)/2-2);
start_r_6 = (corrupted_band(6)+corrupted_band_size(6)/2:corrupted_band_end(6));
start_l_7 = (corrupted_band(7)-1:corrupted_band(7)+corrupted_band_size(7)/2-2);
start_r_7 = (corrupted_band(7)+corrupted_band_size(7)/2:corrupted_band_end(7));
start_l_8 = (corrupted_band(8)-1:corrupted_band(8)+corrupted_band_size(8)/2-2);
start_r_8 = (corrupted_band(8)+corrupted_band_size(8)/2:corrupted_band_end(8));
start_l_9 = (corrupted_band(9)-1:corrupted_band(9)+corrupted_band_size(9)/2-2);
start_r_9 = (corrupted_band(9)+corrupted_band_size(9)/2:corrupted_band_end(9));
start_l_10 = (corrupted_band(10)-1:corrupted_band(10)+corrupted_band_size(10)/2-2);
start_r_10 = (corrupted_band(10)+corrupted_band_size(10)/2:corrupted_band_end(10));

start_l = [start_l_1,start_l_2,start_l_3,start_l_4,start_l_5,start_l_6,start_l_7,start_l_8,start_l_9,start_l_10] -1;
start_r = [start_r_1,start_r_2,start_r_3,start_r_4,start_r_5,start_r_6,start_r_7,start_r_8,start_r_9,start_r_10];
start_r = flip(start_r);

%---------------------------------Marking location where kernel should run---------------------------------
ok_section = [1:148];
for i=1:size(corrupted_band,2)
    ignore_section = corrupted_band(i)-3:corrupted_band_end(i)+1;
    ok_section = setdiff(ok_section,ignore_section);
end

%---------------------------------Restoring the left parts of srtipes---------------------------------
%left part
kernal = zeros(2,3);
corrupted_window = zeros(2,3);
error_map_window = ones(149,148);
for i=1:N   %observed dimention
    for j=start_l %cottupted window's column
        for k=1:150-1   %cottupted window's row
            corrupted_window = S_3d_extend(k:k+1,j:j+2,i);
            for m = 1:149  %kernal's row
                for n = ok_section %kernal's column
                    kernal = S_3d_extend(m:m+1,n:n+2,i);
                    error_map_window(m,n) = compute_error(kernal,corrupted_window,1);
                end
            end
            minimum = min(min(error_map_window));
            [r_index,c_index ]= find(error_map_window==minimum);
            kernal = S_3d_extend(r_index:r_index+1,c_index:c_index+2,i);
            element1 = kernal(1,3);
            element2 = kernal(2,3);
            S_3d_extend(k,j+2,i) = element1;
            S_3d_extend(k+1,j+2,i) = element2;
            fprintf('(Left side) N:%d / Load:%d & %d / into:(%d,%d) (%d,%d) \n', i, element1, element2 ,k,j+2 ,k+1,j+2);
        end
    end
end

%---------------------------------Restoring the right parts of srtipes---------------------------------
%right part
kernal = zeros(2,3);
corrupted_window = zeros(2,3);
error_map_window = ones(149,148);
for i=1:N   %observed dimention
    for j=start_r %cottupted window's column
        for k=1:150-1   %cottupted window's row
            corrupted_window = S_3d_extend(k:k+1,j:j+2,i);
            for m = 1:149  %kernal's row
                for n = ok_section %kernal's column
                    kernal = S_3d_extend(m:m+1,n:n+2,i);
                    error_map_window(m,n) = compute_error(kernal,corrupted_window,2);
                end
            end
            minimum = min(min(error_map_window));
            [r_index,c_index ]= find(error_map_window==minimum);
            kernal = S_3d_extend(r_index:r_index+1,c_index:c_index+2,i);
            element1 = kernal(1,1);
            element2 = kernal(2,1);
            S_3d_extend(k,j,i) = element1;
            S_3d_extend(k+1,j,i) = element2;
            fprintf('(Right side) N:%d / Load:%d & %d / into:(%d,%d) (%d,%d)\n', i,element1, element2 ,k,j,k+1,j);
        end
    end
end
toc

%% Caculate Performance
S_reconsruct_flatten = reshape(S_3d_extend,150*150,N)';
X_reconstructed = A_est*S_reconsruct_flatten;
error_map = X_reconstructed - reshape(X,150*150,183)';
error_map = error_map.^2;
X_gd_reshape = reshape(X,150*150,183)';
error = sqrt( sum( sum(error_map,2)./ sum(X_gd_reshape,2) )/183)*100;
fprintf('Error:%d (percent)\n', error);
missing_rate = sum(corrupted_band_size,'all')/150*100;
fprintf('Missing rate:%d (percent)\n', missing_rate);

%% Plot
f = figure;  
f.Position = [110 110 1700 800];

% original graph (undamaged)
subplot(3,4,1);
original_rgb = zeros(150,150,3);
original_rgb(:,:,1) = X(:,:,1) / mean(X(:,:,1),'all') * 255;
original_rgb(:,:,2) = X(:,:,2) / mean(X(:,:,2),'all') * 140;
original_rgb(:,:,3) = X(:,:,3) / mean(X(:,:,3),'all') * 10;
imshow(uint8(original_rgb));title('original rgb');

% corrupted graph (with stripes of 0s)
subplot(3,4,5);
corrupted_rgb = zeros(150,150,3);
corrupted_rgb(:,:,1) = X_corrupted(:,:,1) / mean(X_corrupted(:,:,1),'all') * 255;
corrupted_rgb(:,:,2) = X_corrupted(:,:,2) / mean(X_corrupted(:,:,2),'all') * 140;
corrupted_rgb(:,:,3) = X_corrupted(:,:,3) / mean(X_corrupted(:,:,3),'all') * 10;
imshow(uint8(corrupted_rgb));title('corruputed rgb');

% cropped graph (stripes removed)
subplot(3,4,9);
cropped_rgb(:,:,1:3) = X_croped(:,:,1:3);
imshow(uint8(cropped_rgb*255));title('cropped rgb');

% S1 (quantity of 1st matter)
subplot(3,4,2);
imshow(S_3d_extend(:,:,1));title('S1 reconstuct');

% S2 (quantity of 2nd matter)
subplot(3,4,6);
imshow(S_3d_extend(:,:,2));title('S2 reconstuct');

% S3 (quantity of 3rd matter)
subplot(3,4,10);
S_3 = S_3d_extend(:,:,3);
imshow(S_3d_extend(:,:,3));title('S3 reconstuct');

% A1 (identity of 1st matter)
subplot(3,4,3);
plot(A_est(:,1)); title('A1 est');
axis([1 183 0 1.5]);

% A2 (identity of 2nd matter)
subplot(3,4,7);
plot(A_est(:,2)); title('A2 est');
axis([1 183 0 1.5]);

% A3 (identity of 3rd matter)
subplot(3,4,11);
plot(A_est(:,3)); title('A3 est');
axis([1 183 0 1.5]);

% X_reconstructed
subplot(3,4,8);
X_reconsructed_3d = reshape(X_reconstructed',[150,150,183]);
reconstruct_rgb = zeros(150,150,3);
reconstruct_rgb(:,:,1) = X_reconsructed_3d(:,:,1) / mean(X_reconsructed_3d(:,:,1),'all') * 255;
reconstruct_rgb(:,:,2) = X_reconsructed_3d(:,:,2) / mean(X_reconsructed_3d(:,:,2),'all') * 140;
reconstruct_rgb(:,:,3) = X_reconsructed_3d(:,:,3) / mean(X_reconsructed_3d(:,:,3),'all') * 10;
imshow(uint8(reconstruct_rgb));title('X Reconstructed');

% Reconstructed graph (A*S (cropped))
subplot(3,4,4);
X_reconstructed_damaged = reshape(X_reconstructed_damaged',[150,width,183]);
imshow(uint8(X_reconstructed_damaged(:,:,1:3)*255));title('Reconstructed damaged rgb');

%% Save mat
Y = X_reconsructed_3d;
save('Data_4.mat','Y');

%% function
function [error] = compute_error(K,C,direction)
if direction==1
    E = K(1:2,1:2) - C(1:2,1:2);
    E = E.^2;
else
    E = K(1:2,2:3) - C(1:2,2:3);
    E = E.^2;
end
error = sum(E, 'all');
end