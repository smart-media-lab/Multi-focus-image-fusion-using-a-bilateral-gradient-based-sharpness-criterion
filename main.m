function fusion_OC_2011_main.m
%
% This is a demo program of the paper J. Tian, L. Chen, L. Ma and W. Yu,
% "Multi-focus image fusion using a bilateral gradient-based sharpness criterion," 
% Optics Communications, Vol. 284, Jan. 2011, pp. 80-87.
%


close all; clear all; clc;

% Load two images with different focus levels
MA = double(imread(['clock_A.bmp']));
MB = double(imread(['clock_B.bmp']));

% Parameter setting
param.win = 5;  %equation (16)
param.alpha = 1; %equation (16)
param.beta = 0.5; %equation (16)

% Proposed sharpness measure see equation (16)
sharp_strength_A = func_sharpness_measuring_strength(MA);
sharp_strength_B = func_sharpness_measuring_strength(MB);
sharp_phase_A = func_sharpness_measuring_phase(MA, param);
sharp_phase_B = func_sharpness_measuring_phase(MB, param);

% Perform fusion
weight_A = sharp_strength_A.^param.alpha.*sharp_phase_A.^param.beta;
weight_B = sharp_strength_B.^param.alpha.*sharp_phase_B.^param.beta;

mask1 = (weight_A>=weight_B);
mask2 = (weight_A<weight_B);
MF = MA.*mask1 + MB.*mask2;

% Write the output image
imwrite(uint8(MF), ['clock_fused.bmp'], 'bmp');

% Calculate the performance using two criterions
fprintf('Mutural Information is %.2f\n', func_evaluate_mutural_information(MA, MB, MF, 256));
fprintf('Spatial Frequency is %.2f\n', func_evaluate_spatial_frequency(MA, MB, MF));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Inner Function %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function result = func_sharpness_measuring_strength(img)
% Sharpness measure using gradient strength (see equation (14))

dx = [-1 0 1; -2 0 2; -1 0 1];
dy = dx';
Ix = conv2(img, dx, 'same');
Iy = conv2(img, dy, 'same');
Ix2 = Ix.^2;
Iy2 = Iy.^2;
Ixy = Ix.*Iy;
aa = Ix2;
bb = Ixy;
cc = Ixy;
dd = Iy2;
Eig11 = (aa+dd)./2 + sqrt(((aa+dd).*(aa+dd))./4 + bb.*cc-aa.*dd);
Eig22 = (aa+dd)./2 - sqrt(((aa+dd).*(aa+dd))./4 + bb.*cc-aa.*dd);
Eig11(isnan(Eig11))=0;
Eig22(isnan(Eig22))=0;
result = Eig11;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Inner Function %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sharpness measure using phase coherence (see equation (15))

function result = func_sharpness_measuring_phase(img, param)
win_size = param.win;
img = padarray(img,[(win_size-1)/2 (win_size-1)/2],'replicate','both');
dx = [-1 0 1; -2 0 2; -1 0 1];
dy = dx';
Ix = conv2(img, dx, 'same');
Iy = conv2(img, dy, 'same');
theta = atan2(Iy,Ix);
for i=(win_size-1)/2+1:size(img,1)-(win_size-1)/2
    for j=(win_size-1)/2+1:size(img,2)-(win_size-1)/2
        temp = theta(i-(win_size-1)/2:i+(win_size-1)/2,j-(win_size-1)/2:j+(win_size-1)/2);        
        temp = temp-mean(mean(temp));
        result(i-(win_size-1)/2,j-(win_size-1)/2) = sum(sum(cos(temp)));
    end
end
result = -result;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Inner Function %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Evaluatoin of image fusion algorithm using mutural information
function mutural_informationR = func_evaluate_mutural_information(grey_matrixA,grey_matrixB,grey_matrixF,grey_level)

HA=entropy_fusion(grey_matrixA,grey_level);
HB=entropy_fusion(grey_matrixB,grey_level);
HF=entropy_fusion(grey_matrixF,grey_level);
HFA=Hab(grey_matrixF,grey_matrixA,grey_level);
HFB=Hab(grey_matrixF,grey_matrixB,grey_level);
MIFA=HA+HF-HFA;
MIFB=HB+HF-HFB;
mutural_informationR=MIFA+MIFB;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Inner Function %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Evaluatoin of image fusion algorithm using spatial frequency
function entropyR=entropy_fusion(grey_matrix,grey_level)
[row,column]=size(grey_matrix);
total=row*column;

counter=zeros(1,grey_level);
grey_matrix=grey_matrix+1;
for i=1:row
    for j=1:column
        indexx= grey_matrix(i,j);
        counter(indexx)=counter(indexx)+1;
    end
end
total= sum(counter(:));
index = find(counter~=0);
 p = counter/total;
entropyR= sum(sum(-p(index).*log2(p(index))));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Inner Function %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function HabR=Hab(grey_matrixA,grey_matrixB,grey_level)
% HabR=Hab(grey_matrixA,grey_matrixB,grey_level)
% compute mutural information of the image
% grey_matrixA , grey_matrixB,grey_matrixF are grey values of imageA,imageB and fusion image
% grey_level is the grayscale degree of image
% --------- 
[row,column]=size(grey_matrixA);
counter = zeros(256,256);
grey_matrixA=grey_matrixA+1;
grey_matrixB=grey_matrixB+1;
for i=1:row
    for j=1:column
        indexx = grey_matrixA(i,j);
        indexy = grey_matrixB(i,j);
        counter(indexx,indexy) = counter(indexx,indexy)+1;
    end
end
total= sum(counter(:));
index = find(counter~=0);
p = counter/total;
HabR = sum(sum(-p(index).*log2(p(index))));        

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Inner Function %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function SF = func_evaluate_spatial_frequency(MA, MB, MF)
RF = diff(MF,1,1);
RF = sqrt(mean(mean(RF.^2)));
CF = diff(MF,1,2);
CF = sqrt(mean(mean(CF.^2)));
SF = sqrt(RF^2+CF^2);

