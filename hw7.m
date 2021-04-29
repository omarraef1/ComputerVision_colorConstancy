function hw7()
% task: infer or correct chromaticity
% provide triples: origina, corrected, canonical

close all;
format compact;

%canonicalImg(351, 122, 1) = 250;
%canonicalImg(351, 122, 2) = 250;
%canonicalImg(351, 122, 3) = 250;

canonicalImg = imread('macbeth_syl-50MR16Q.tif');
figure(1),
imshow(canonicalImg);
datacursormode on
%whos('canonicalImg')

% part 1
% average and scale illuminant
% 110, 340 to 135, 365

illSet = canonicalImg(340:365,110:135,:);
whos('illSet')
rSet = illSet(:,:,1);
gSet = illSet(:,:,2);
bSet = illSet(:,:,3);

meanR = mean(mean(rSet));
meanG = mean(mean(gSet));
meanB = mean(mean(bSet));
meanR
meanG
meanB
% scale factor = 1.87816
sf = 1.87816;
nuR = meanR*sf;
nuG = meanG*sf;
nuB = meanB*sf;
nuR
nuG
nuB
% max number (B value) is now 250

%part 2: same as part 1

img2 = imread('macbeth_solux-4100.tif');
figure(2),
imshow(img2);
datacursormode on

% average and scale illuminant
% 110, 340 to 135, 365

illSet2 = img2(340:365,110:135,:);
whos('illSet')
rSet2 = illSet2(:,:,1);
gSet2 = illSet2(:,:,2);
bSet2 = illSet2(:,:,3);

meanR2 = mean(mean(rSet2));
meanG2 = mean(mean(gSet2));
meanB2 = mean(mean(bSet2));
meanR2
meanG2
meanB2
% scale factor = 1.7304
sf2 = 1.7304;
nuR2 = meanR2*sf2;
nuG2 = meanG2*sf2;
nuB2 = meanB2*sf2;
nuR2
nuG2
nuB2

% max number (B value) is now 250

% part 3 angular error between 1 and 2
firstRGB = [nuR nuG nuB];
secondRGB = [nuR2 nuG2 nuB2];
dotP = dot(firstRGB, secondRGB);
m1 = sqrt(nuR^2 + nuG^2 + nuB^2);
m2 = sqrt(nuR2^2 + nuG2^2 + nuB2^2);
prod = dotP/(m1*m2);
angle = acosd(prod);
angle


%THIS COMMENT STUB STATES THAT 
%THIS CODE IS THE PROPERTY OF OMAR R.G. (UofA Student)

% part 4 lec 15 diagonal model

% try with zeros then ones then other equation
img3 = imread('macbeth_solux-4100.tif');
ratioMat = [nuR/nuR2 nuG/nuG2 nuB/nuB2];
for i = 1:637
    for j = 1:468
        img3(j,i,1) = img2(j,i,1)*ratioMat(1);
        img3(j,i,2) = img2(j,i,2)*ratioMat(2);
        img3(j,i,3) = img2(j,i,3)*ratioMat(3);
    end
end
figure(3),
imshow(img3);
datacursormode on

mxrcanon = max(max(canonicalImg(:,:,1)));
mxgcanon = max(max(canonicalImg(:,:,2)));
mxbcanon = max(max(canonicalImg(:,:,3)));
mxrcanon
mxgcanon
mxbcanon


mxrblu = max(max(img2(:,:,1)));
mxgblu = max(max(img2(:,:,2)));
mxbblu = max(max(img2(:,:,3)));
mxrblu
mxgblu
mxbblu

mxrcorr = max(max(img3(:,:,1)));
mxgcorr = max(max(img3(:,:,2)));
mxbcorr = max(max(img3(:,:,3)));
mxrcorr
mxgcorr
mxbcorr

% after scaling
scaledCan = imread('macbeth_syl-50MR16Q.tif');
scaledBlu = imread('macbeth_solux-4100.tif');
scaledCorr = img3;
for i = 1:637
    for j = 1:468
        scaledCan(j,i,1) = canonicalImg(j,i,1) * 1.4368;
        scaledCan(j,i,2) = canonicalImg(j,i,2) * 1.4368;
        scaledCan(j,i,3) = canonicalImg(j,i,3) * 1.4368;
        
        scaledBlu(j,i,1) = img2(j,i,1) * 1.4205;
        scaledBlu(j,i,2) = img2(j,i,2) * 1.4205;
        scaledBlu(j,i,3) = img2(j,i,3) * 1.4205;
        
        scaledCorr(j,i,1) = img3(j,i,1) * 1.2755;
        scaledCorr(j,i,2) = img3(j,i,2) * 1.2755;
        scaledCorr(j,i,3) = img3(j,i,3) * 1.2755;
        
    end
end

scaledImgs = [scaledBlu scaledCorr scaledCan];
figure(4),
imshow(scaledImgs);
datacursormode on

% part 5 ignore dark pixels with r+g+b < 10
%A) rms between original (blu) and canonical
%B) rms between correction and canonical

%A)original (blu) and canonical
bluR = scaledBlu(:,:,1)./((scaledBlu(:,:,1)+scaledBlu(:,:,2)+scaledBlu(:,:,3)));
bluG = scaledBlu(:,:,2)./((scaledBlu(:,:,1)+scaledBlu(:,:,2)+scaledBlu(:,:,3)));

canR = scaledCan(:,:,1)./((scaledCan(:,:,1)+scaledCan(:,:,2)+scaledCan(:,:,3)));
canG = scaledCan(:,:,2)./((scaledCan(:,:,1)+scaledCan(:,:,2)+scaledCan(:,:,3)));

newbluR = [];
newbluG = [];

newcanR = [];
newcanG = [];

for i = 1:637
    for j = 1:468
        if((scaledBlu(j,i,1)+scaledBlu(j,i,2)+scaledBlu(j,i,3)) >= 10)
            newbluR = [newbluR; scaledBlu(j,i,1)./((scaledBlu(j,i,1)+scaledBlu(j,i,2)+scaledBlu(j,i,3)))];
            newbluG = [newbluG; scaledBlu(j,i,2)./((scaledBlu(j,i,1)+scaledBlu(j,i,2)+scaledBlu(j,i,3)))];
            newcanR = [newcanR; scaledCan(j,i,1)./((scaledCan(j,i,1)+scaledCan(j,i,2)+scaledCan(j,i,3)))];
            newcanG = [newcanG; scaledCan(j,i,2)./((scaledCan(j,i,1)+scaledCan(j,i,2)+scaledCan(j,i,3)))];
        
        end
    end
end
FinalbluR = newbluR;
FinalbluG = newbluG;
FinalcanR = newcanR;
FinalcanG = newcanG;

diffR = double(FinalcanR-FinalbluR);
diffRS = diffR.^2;
sumR = sum(diffRS);
rmsr = sqrt(sum(sumR/(468*637)));
rmsr

diffg = double(FinalcanG-FinalbluG);
diffgS = diffg.^2;
sumg = sum(diffgS);
rmsg = sqrt(sum(sumg/(468*637)));
rmsg

%B)correction and canonical
corrR = scaledCorr(:,:,1)./((scaledCorr(:,:,1)+scaledCorr(:,:,2)+scaledCorr(:,:,3)));
corrG = scaledCorr(:,:,2)./((scaledCorr(:,:,1)+scaledCorr(:,:,2)+scaledCorr(:,:,3)));

newcorrR = [];
newcorrG = [];

newcanR = [];
newcanG = [];

for i = 1:637
    for j = 1:468
        if((scaledCorr(j,i,1)+scaledCorr(j,i,2)+scaledCorr(j,i,3)) >= 10)
            newcorrR = [newcorrR; scaledCorr(j,i,1)./((scaledCorr(j,i,1)+scaledCorr(j,i,2)+scaledCorr(j,i,3)))];
            newcorrG = [newcorrG; scaledCorr(j,i,2)./((scaledCorr(j,i,1)+scaledCorr(j,i,2)+scaledCorr(j,i,3)))];
            newcanR = [newcanR; scaledCan(j,i,1)./((scaledCan(j,i,1)+scaledCan(j,i,2)+scaledCan(j,i,3)))];
            newcanG = [newcanG; scaledCan(j,i,2)./((scaledCan(j,i,1)+scaledCan(j,i,2)+scaledCan(j,i,3)))];
        
        end
    end
end
FinalcorrR = newcorrR;
FinalcorrG = newcorrG;
Final2canR = newcanR;
Final2canG = newcanG;

diffR2 = double(Final2canR-FinalcorrR);
diffRS2 = diffR2.^2;
sumR2 = sum(diffRS2);
rmsr2 = sqrt(sum(sumR2/(468*637)));
rmsr2

diffg2 = double(Final2canG-FinalcorrG);
diffgS2 = diffg2.^2;
sumg2 = sum(diffgS2);
rmsg2 = sqrt(sum(sumg2/(468*637)));
rmsg2

%part 6: 
% A) MaxRGB estimate light color
% B) angular errors

% A)
apples = imread('apples2_solux-4100.tif');
ball = imread('ball_solux-4100.tif');
blocks = imread('blocks1_solux-4100.tif');

mxrapple = double(max(max(apples(:,:,1))));
mxgapple = double(max(max(apples(:,:,2))));
mxbapple = double(max(max(apples(:,:,3))));
mxrapple
mxgapple
mxbapple

mxrball = double(max(max(ball(:,:,1))));
mxgball = double(max(max(ball(:,:,2))));
mxbball = double(max(max(ball(:,:,3))));
mxrball
mxgball
mxbball

mxrblocks = double(max(max(blocks(:,:,1))));
mxgblocks = double(max(max(blocks(:,:,2))));
mxbblocks = double(max(max(blocks(:,:,3))));
mxrblocks
mxgblocks
mxbblocks

% B)
applesRGB = [mxrapple mxgapple mxbapple];
ballRGB = [mxrball mxgball mxbball];
blocksRGB = [mxrblocks mxgblocks mxbblocks];

%apple angular error
dotPap = dot(applesRGB, secondRGB);
m1ap = sqrt(mxrapple^2 + mxgapple^2 + mxbapple^2);
m2ap = sqrt(nuR2^2 + nuG2^2 + nuB2^2);
prodap = dotPap/(m1ap*m2ap);
angleApple = acosd(prodap);
angleApple
%ball angular error
dotPball = dot(ballRGB, secondRGB);
m1ball = sqrt(mxrball^2 + mxgball^2 + mxbball^2);
m2ball = sqrt(nuR2^2 + nuG2^2 + nuB2^2);
prodball = dotPball/(m1ball*m2ball);
angleBall = acosd(prodball);
angleBall
%blocks angular error
dotPblocks = dot(blocksRGB, secondRGB);
m1blocks = sqrt(mxrblocks^2 + mxgblocks^2 + mxbblocks^2);
m2blocks = sqrt(nuR2^2 + nuG2^2 + nuB2^2);
prodblocks = dotPblocks/(m1blocks*m2blocks);
angleBlocks = acosd(prodblocks);
angleBlocks

% part 7:
%canonicals
applesCan = imread('apples2_syl-50MR16Q.tif');
ballCan = imread('ball_syl-50MR16Q.tif');
blocksCan = imread('blocks1_syl-50MR16Q.tif');
whos('applesCan')
whos('ballCan')
whos('blocksCan')

mxrapple
mxgapple
mxbapple

mxrball
mxgball
mxbball

mxrblocks
mxgblocks
mxbblocks

ratioMatApple = [nuR/mxrapple nuG/mxgapple nuB/mxbapple];
ratioMatBall = [nuR/mxrball nuG/mxgball nuB/mxbball];
ratioMatBlock = [nuR/mxrblocks nuG/mxgblocks nuB/mxbblocks];

applesCorr = imread('apples2_solux-4100.tif');
ballCorr = imread('ball_solux-4100.tif');
blocksCorr = imread('blocks1_solux-4100.tif');

for i = 1:637
    for j = 1:468
        applesCorr(j,i,1) = apples(j,i,1)*ratioMatApple(1);
        applesCorr(j,i,2) = apples(j,i,2)*ratioMatApple(2);
        applesCorr(j,i,3) = apples(j,i,3)*ratioMatApple(3);
        
        ballCorr(j,i,1) = ball(j,i,1)*ratioMatBall(1);
        ballCorr(j,i,2) = ball(j,i,2)*ratioMatBall(2);
        ballCorr(j,i,3) = ball(j,i,3)*ratioMatBall(3);
        
        blocksCorr(j,i,1) = blocks(j,i,1)*ratioMatBlock(1);
        blocksCorr(j,i,2) = blocks(j,i,2)*ratioMatBlock(2);
        blocksCorr(j,i,3) = blocks(j,i,3)*ratioMatBlock(3);
    end
end

scaledNine = [apples applesCorr applesCan; ball ballCorr ballCan; blocks blocksCorr blocksCan];
figure(5),
imshow(scaledNine);

%rms errors of images (like part 5)
% Ignore dark pixels with r + g + b < 10

%______APPLES____%

%A)original (apples) and canonical
orApR = apples(:,:,1)./((apples(:,:,1)+apples(:,:,2)+apples(:,:,3)));
orApG = apples(:,:,2)./((apples(:,:,1)+apples(:,:,2)+apples(:,:,3)));

canApR = applesCan(:,:,1)./((applesCan(:,:,1)+applesCan(:,:,2)+applesCan(:,:,3)));
canApG = applesCan(:,:,2)./((applesCan(:,:,1)+applesCan(:,:,2)+applesCan(:,:,3)));

neworApR = [];
neworApG = [];

newcanApR = [];
newcanApG = [];

for i = 1:637
    for j = 1:468
        if((apples(j,i,1)+apples(j,i,2)+apples(j,i,3)) >= 10)
            neworApR = [neworApR; apples(j,i,1)./((apples(j,i,1)+apples(j,i,2)+apples(j,i,3)))];
            neworApG = [neworApG; apples(j,i,2)./((apples(j,i,1)+apples(j,i,2)+apples(j,i,3)))];
            newcanApR = [newcanApR; applesCan(j,i,1)./((applesCan(j,i,1)+applesCan(j,i,2)+applesCan(j,i,3)))];
            newcanApG = [newcanApG; applesCan(j,i,2)./((applesCan(j,i,1)+applesCan(j,i,2)+applesCan(j,i,3)))];
        
        end
    end
end
FinalorApR = neworApR;
FinalorApG = neworApG;
FinalcanApR = newcanApR;
FinalcanApG = newcanApG;

diffRApp1 = double(FinalcanApR-FinalorApR);
diffRSApp1 = diffRApp1.^2;
sumRApp1 = sum(diffRSApp1);
rmsrApp1 = sqrt(sum(sumRApp1/(468*637)));
rmsrApp1

diffgApp1 = double(FinalcanApG-FinalorApG);
diffgSApp1 = diffgApp1.^2;
sumgApp1 = sum(diffgSApp1);
rmsgApp1 = sqrt(sum(sumgApp1/(468*637)));
rmsgApp1

%B)correction apples and canonical
corrApR = applesCorr(:,:,1)./((applesCorr(:,:,1)+applesCorr(:,:,2)+applesCorr(:,:,3)));
corrApG = applesCorr(:,:,2)./((applesCorr(:,:,1)+applesCorr(:,:,2)+applesCorr(:,:,3)));

newcorrApR = [];
newcorrApG = [];

newcanApR = [];
newcanApG = [];

for i = 1:637
    for j = 1:468
        if((applesCorr(j,i,1)+applesCorr(j,i,2)+applesCorr(j,i,3)) >= 10)
            newcorrApR = [newcorrApR; applesCorr(j,i,1)./((applesCorr(j,i,1)+applesCorr(j,i,2)+applesCorr(j,i,3)))];
            newcorrApG = [newcorrApG; applesCorr(j,i,2)./((applesCorr(j,i,1)+applesCorr(j,i,2)+applesCorr(j,i,3)))];
            newcanApR = [newcanApR; applesCan(j,i,1)./((applesCan(j,i,1)+applesCan(j,i,2)+applesCan(j,i,3)))];
            newcanApG = [newcanApG; applesCan(j,i,2)./((applesCan(j,i,1)+applesCan(j,i,2)+applesCan(j,i,3)))];
        
        end
    end
end
FinalcorrApR = newcorrApR;
FinalcorrApG = newcorrApG;
Final2canApR = newcanApR;
Final2canApG = newcanApG;

diffR2App1 = double(Final2canApR-FinalcorrApR);
diffRS2App1 = diffR2App1.^2;
sumR2App1 = sum(diffRS2App1);
rmsr2App1 = sqrt(sum(sumR2App1/(468*637)));
rmsr2App1

diffg2App1 = double(Final2canApG-FinalcorrApG);
diffgS2App1 = diffg2App1.^2;
sumg2App1 = sum(diffgS2App1);
rmsg2App1 = sqrt(sum(sumg2App1/(468*637)));
rmsg2App1

%______APPLES_END____%

%______BALL____%

%A)original (ball) and canonical
orballR = ball(:,:,1)./((ball(:,:,1)+ball(:,:,2)+ball(:,:,3)));
orballG = ball(:,:,2)./((ball(:,:,1)+ball(:,:,2)+ball(:,:,3)));

canballR = ballCan(:,:,1)./((ballCan(:,:,1)+ballCan(:,:,2)+ballCan(:,:,3)));
canballG = ballCan(:,:,2)./((ballCan(:,:,1)+ballCan(:,:,2)+ballCan(:,:,3)));

neworballR = [];
neworballG = [];

newcanballR = [];
newcanballG = [];

for i = 1:637
    for j = 1:468
        if((ball(j,i,1)+ball(j,i,2)+ball(j,i,3)) >= 10)
            neworballR = [neworballR; ball(j,i,1)./((ball(j,i,1)+ball(j,i,2)+ball(j,i,3)))];
            neworballG = [neworballG; ball(j,i,2)./((ball(j,i,1)+ball(j,i,2)+ball(j,i,3)))];
            newcanballR = [newcanballR; ballCan(j,i,1)./((ballCan(j,i,1)+ballCan(j,i,2)+ballCan(j,i,3)))];
            newcanballG = [newcanballG; ballCan(j,i,2)./((ballCan(j,i,1)+ballCan(j,i,2)+ballCan(j,i,3)))];
        
        end
    end
end
FinalorballR = neworballR;
FinalorballG = neworballG;
FinalcanballR = newcanballR;
FinalcanballG = newcanballG;

diffRBall1 = double(FinalcanballR-FinalorballR);
diffRSBall1 = diffRBall1.^2;
sumRBall1 = sum(diffRSBall1);
rmsrBall1 = sqrt(sum(sumRBall1/(468*637)));
rmsrBall1

diffgBall1 = double(FinalcanballG-FinalorballG);
diffgSBall1 = diffgBall1.^2;
sumgBall1 = sum(diffgSBall1);
rmsgBall1 = sqrt(sum(sumgBall1/(468*637)));
rmsgBall1

%B)correction ball and canonical
corrballR = ballCorr(:,:,1)./((ballCorr(:,:,1)+ballCorr(:,:,2)+ballCorr(:,:,3)));
corrballG = ballCorr(:,:,2)./((ballCorr(:,:,1)+ballCorr(:,:,2)+ballCorr(:,:,3)));

newcorrballR = [];
newcorrballG = [];

newcanballR = [];
newcanballG = [];

for i = 1:637
    for j = 1:468
        if((ballCorr(j,i,1)+ballCorr(j,i,2)+ballCorr(j,i,3)) >= 10)
            newcorrballR = [newcorrballR; ballCorr(j,i,1)./((ballCorr(j,i,1)+ballCorr(j,i,2)+ballCorr(j,i,3)))];
            newcorrballG = [newcorrballG; ballCorr(j,i,2)./((ballCorr(j,i,1)+ballCorr(j,i,2)+ballCorr(j,i,3)))];
            newcanballR = [newcanballR; ballCan(j,i,1)./((ballCan(j,i,1)+ballCan(j,i,2)+ballCan(j,i,3)))];
            newcanballG = [newcanballG; ballCan(j,i,2)./((ballCan(j,i,1)+ballCan(j,i,2)+ballCan(j,i,3)))];
        
        end
    end
end
FinalcorrballR = newcorrballR;
FinalcorrballG = newcorrballG;
Final2canballR = newcanballR;
Final2canballG = newcanballG;

diffR2ball1 = double(Final2canballR-FinalcorrballR);
diffRS2ball1 = diffR2ball1.^2;
sumR2ball1 = sum(diffRS2ball1);
rmsr2Ball1 = sqrt(sum(sumR2ball1/(468*637)));
rmsr2Ball1

diffg2ball1 = double(Final2canballG-FinalcorrballG);
diffgS2ball1 = diffg2ball1.^2;
sumg2ball1 = sum(diffgS2ball1);
rmsg2Ball1 = sqrt(sum(sumg2ball1/(468*637)));
rmsg2Ball1

%______BALL_END____%

%______Blocks____%

%A)original (blocks) and canonical
orblocksR = blocks(:,:,1)./((blocks(:,:,1)+blocks(:,:,2)+blocks(:,:,3)));
orblocksG = blocks(:,:,2)./((blocks(:,:,1)+blocks(:,:,2)+blocks(:,:,3)));

canblocksR = blocksCan(:,:,1)./((blocksCan(:,:,1)+blocksCan(:,:,2)+blocksCan(:,:,3)));
canblocksG = blocksCan(:,:,2)./((blocksCan(:,:,1)+blocksCan(:,:,2)+blocksCan(:,:,3)));

neworblocksR = [];
neworblocksG = [];

newcanblocksR = [];
newcanblocksG = [];

for i = 1:637
    for j = 1:468
        if((blocks(j,i,1)+blocks(j,i,2)+blocks(j,i,3)) >= 10)
            neworblocksR = [neworblocksR; blocks(j,i,1)./((blocks(j,i,1)+blocks(j,i,2)+blocks(j,i,3)))];
            neworblocksG = [neworblocksG; blocks(j,i,2)./((blocks(j,i,1)+blocks(j,i,2)+blocks(j,i,3)))];
            newcanblocksR = [newcanblocksR; blocksCan(j,i,1)./((blocksCan(j,i,1)+blocksCan(j,i,2)+blocksCan(j,i,3)))];
            newcanblocksG = [newcanblocksG; blocksCan(j,i,2)./((blocksCan(j,i,1)+blocksCan(j,i,2)+blocksCan(j,i,3)))];
        
        end
    end
end
FinalorblocksR = neworblocksR;
FinalorblocksG = neworblocksG;
FinalcanblocksR = newcanblocksR;
FinalcanblocksG = newcanblocksG;

diffRblocks1 = double(FinalcanblocksR-FinalorblocksR);
diffRSblocks1 = diffRblocks1.^2;
sumRblocks1 = sum(diffRSblocks1);
rmsrBlocks1 = sqrt(sum(sumRblocks1/(468*637)));
rmsrBlocks1

diffgblocks1 = double(FinalcanblocksG-FinalorblocksG);
diffgSblocks1 = diffgblocks1.^2;
sumgblocks1 = sum(diffgSblocks1);
rmsgBlocks1 = sqrt(sum(sumgblocks1/(468*637)));
rmsgBlocks1

%B)correction blocks and canonical
corrblocksR = blocksCorr(:,:,1)./((blocksCorr(:,:,1)+blocksCorr(:,:,2)+blocksCorr(:,:,3)));
corrblocksG = blocksCorr(:,:,2)./((blocksCorr(:,:,1)+blocksCorr(:,:,2)+blocksCorr(:,:,3)));

newcorrblocksR = [];
newcorrblocksG = [];

newcanblocksR = [];
newcanblocksG = [];

for i = 1:637
    for j = 1:468
        if((blocksCorr(j,i,1)+blocksCorr(j,i,2)+blocksCorr(j,i,3)) >= 10)
            newcorrblocksR = [newcorrblocksR; blocksCorr(j,i,1)./((blocksCorr(j,i,1)+blocksCorr(j,i,2)+blocksCorr(j,i,3)))];
            newcorrblocksG = [newcorrblocksG; blocksCorr(j,i,2)./((blocksCorr(j,i,1)+blocksCorr(j,i,2)+blocksCorr(j,i,3)))];
            newcanblocksR = [newcanblocksR; blocksCan(j,i,1)./((blocksCan(j,i,1)+blocksCan(j,i,2)+blocksCan(j,i,3)))];
            newcanblocksG = [newcanblocksG; blocksCan(j,i,2)./((blocksCan(j,i,1)+blocksCan(j,i,2)+blocksCan(j,i,3)))];
        
        end
    end
end
FinalcorrblocksR = newcorrblocksR;
FinalcorrblocksG = newcorrblocksG;
Final2canblocksR = newcanblocksR;
Final2canblocksG = newcanblocksG;

diffR2blocks1 = double(Final2canblocksR-FinalcorrblocksR);
diffRS2blocks1 = diffR2blocks1.^2;
sumR2blocks1 = sum(diffRS2blocks1);
rmsr2Blocks1 = sqrt(sum(sumR2blocks1/(468*637)));
rmsr2Blocks1

diffg2blocks1 = double(Final2canblocksG-FinalcorrblocksG);
diffgS2blocks1 = diffg2blocks1.^2;
sumg2blocks1 = sum(diffgS2blocks1);
rmsg2Blocks1 = sqrt(sum(sumg2blocks1/(468*637)));
rmsg2Blocks1

%______Blocks_END____%

%Is there good agreement between________________________________________________________________________%%
%the general trend of the (r,g) RMS error 
%and the angular error computed in problem 6)


%part 8: gray world on apples ball and blocks
rAppSet = apples(:,:,1);
gAppSet = apples(:,:,2);
bAppSet = apples(:,:,3);
meanRapp = 2*double(mean(mean(rAppSet)));
meanGapp = 2*double(mean(mean(gAppSet)));
meanBapp = 2*double(mean(mean(bAppSet)));
meanRapp
meanGapp
meanBapp

rballSet = ball(:,:,1);
gballSet = ball(:,:,2);
bballSet = ball(:,:,3);
meanRball = 2*double(mean(mean(rballSet)));
meanGball = 2*double(mean(mean(gballSet)));
meanBball = 2*double(mean(mean(bballSet)));
meanRball
meanGball
meanBball

rblocksSet = blocks(:,:,1);
gblocksSet = blocks(:,:,2);
bblocksSet = blocks(:,:,3);
meanRblocks = 2*double(mean(mean(rblocksSet)));
meanGblocks = 2*double(mean(mean(gblocksSet)));
meanBblocks = 2*double(mean(mean(bblocksSet)));
meanRblocks
meanGblocks
meanBblocks


% Angular Errors (like part 6)
applesRGB2 = [meanRapp meanGapp meanBapp];
ballRGB2 = [meanRball meanGball meanBball];
blocksRGB2 = [meanRblocks meanGblocks meanBblocks];

%apple angular error
dotPap2 = dot(applesRGB2, secondRGB);
m1ap2 = sqrt(meanRapp^2 + meanGapp^2 + meanBapp^2);
m2ap2 = sqrt(nuR2^2 + nuG2^2 + nuB2^2);
prodap2 = dotPap2/(m1ap2*m2ap2);
angleApple2 = acosd(prodap2);
angleApple2
%ball angular error
dotPball2 = dot(ballRGB2, secondRGB);
m1ball2 = sqrt(meanRball^2 + meanGball^2 + meanBball^2);
m2ball2 = sqrt(nuR2^2 + nuG2^2 + nuB2^2);
prodball2 = dotPball2/(m1ball2*m2ball2);
angleBall2 = acosd(prodball2);
angleBall2
%blocks angular error
dotPblocks2 = dot(blocksRGB2, secondRGB);
m1blocks2 = sqrt(meanRblocks^2 + meanGblocks^2 + meanBblocks^2);
m2blocks2 = sqrt(nuR2^2 + nuG2^2 + nuB2^2);
prodblocks2 = dotPblocks2/(m1blocks2*m2blocks2);
angleBlocks2 = acosd(prodblocks2);
angleBlocks2

% RMS Errors (r,g) (from part 7)


ratioMatApple2 = [nuR/meanRapp nuG/meanGapp nuB/meanBapp];
ratioMatBall2 = [nuR/meanRball nuG/meanGball nuB/meanBball];
ratioMatBlock2 = [nuR/meanRblocks nuG/meanGblocks nuB/meanBblocks];

applesCorr2 = imread('apples2_solux-4100.tif');
ballCorr2 = imread('ball_solux-4100.tif');
blocksCorr2 = imread('blocks1_solux-4100.tif');

for i = 1:637
    for j = 1:468
        applesCorr2(j,i,1) = apples(j,i,1)*ratioMatApple2(1);
        applesCorr2(j,i,2) = apples(j,i,2)*ratioMatApple2(2);
        applesCorr2(j,i,3) = apples(j,i,3)*ratioMatApple2(3);
        
        ballCorr2(j,i,1) = ball(j,i,1)*ratioMatBall2(1);
        ballCorr2(j,i,2) = ball(j,i,2)*ratioMatBall2(2);
        ballCorr2(j,i,3) = ball(j,i,3)*ratioMatBall2(3);
        
        blocksCorr2(j,i,1) = blocks(j,i,1)*ratioMatBlock2(1);
        blocksCorr2(j,i,2) = blocks(j,i,2)*ratioMatBlock2(2);
        blocksCorr2(j,i,3) = blocks(j,i,3)*ratioMatBlock2(3);
    end
end

scaledNine2 = [apples applesCorr2 applesCan; ball ballCorr2 ballCan; blocks blocksCorr2 blocksCan];
figure(6),
imshow(scaledNine2);


%Apples 
%B)correction apples and canonical
corrApRgray = applesCorr2(:,:,1)./((applesCorr2(:,:,1)+applesCorr2(:,:,2)+applesCorr2(:,:,3)));
corrApGgray = applesCorr2(:,:,2)./((applesCorr2(:,:,1)+applesCorr2(:,:,2)+applesCorr2(:,:,3)));

newcorrApRgray = [];
newcorrApGgray = [];

newcanApRgray = [];
newcanApGgray = [];

for i = 1:637
    for j = 1:468
        if((applesCorr2(j,i,1)+applesCorr2(j,i,2)+applesCorr2(j,i,3)) >= 10)
            newcorrApRgray = [newcorrApRgray; applesCorr2(j,i,1)./((applesCorr2(j,i,1)+applesCorr2(j,i,2)+applesCorr2(j,i,3)))];
            newcorrApGgray = [newcorrApGgray; applesCorr2(j,i,2)./((applesCorr2(j,i,1)+applesCorr2(j,i,2)+applesCorr2(j,i,3)))];
            newcanApRgray = [newcanApRgray; applesCan(j,i,1)./((applesCan(j,i,1)+applesCan(j,i,2)+applesCan(j,i,3)))];
            newcanApGgray = [newcanApGgray; applesCan(j,i,2)./((applesCan(j,i,1)+applesCan(j,i,2)+applesCan(j,i,3)))];
        
        end
    end
end
FinalcorrApRgray = newcorrApRgray;
FinalcorrApGgray = newcorrApGgray;
Final2canApRgray = newcanApRgray;
Final2canApGgray = newcanApGgray;

diffR2App1gray = double(Final2canApRgray-FinalcorrApRgray);
diffRS2App1gray = diffR2App1gray.^2;
sumR2App1gray = sum(diffRS2App1gray);
rmsr2App1gray = sqrt(sum(sumR2App1gray/(468*637)));
rmsr2App1gray

diffg2App1gray = double(Final2canApGgray-FinalcorrApGgray);
diffgS2App1gray = diffg2App1gray.^2;
sumg2App1gray = sum(diffgS2App1gray);
rmsg2App1gray = sqrt(sum(sumg2App1gray/(468*637)));
rmsg2App1gray
%______Apples end____%

%______BALL____%
%B)correction ball and canonical
corrballRgray = ballCorr2(:,:,1)./((ballCorr2(:,:,1)+ballCorr2(:,:,2)+ballCorr2(:,:,3)));
corrballGgray = ballCorr2(:,:,2)./((ballCorr2(:,:,1)+ballCorr2(:,:,2)+ballCorr2(:,:,3)));

newcorrballRgray = [];
newcorrballGgray = [];

newcanballRgray = [];
newcanballGgray = [];

for i = 1:637
    for j = 1:468
        if((ballCorr2(j,i,1)+ballCorr2(j,i,2)+ballCorr2(j,i,3)) >= 10)
            newcorrballRgray = [newcorrballRgray; ballCorr2(j,i,1)./((ballCorr2(j,i,1)+ballCorr2(j,i,2)+ballCorr2(j,i,3)))];
            newcorrballGgray = [newcorrballGgray; ballCorr2(j,i,2)./((ballCorr2(j,i,1)+ballCorr2(j,i,2)+ballCorr2(j,i,3)))];
            newcanballRgray = [newcanballRgray; ballCan(j,i,1)./((ballCan(j,i,1)+ballCan(j,i,2)+ballCan(j,i,3)))];
            newcanballGgray = [newcanballGgray; ballCan(j,i,2)./((ballCan(j,i,1)+ballCan(j,i,2)+ballCan(j,i,3)))];
        
        end
    end
end
FinalcorrballRgray = newcorrballRgray;
FinalcorrballGgray = newcorrballGgray;
Final2canballRgray = newcanballRgray;
Final2canballGgray = newcanballGgray;

diffR2ball1gray = double(Final2canballRgray-FinalcorrballRgray);
diffRS2ball1gray = diffR2ball1gray.^2;
sumR2ball1gray = sum(diffRS2ball1gray);
rmsr2Ball1gray = sqrt(sum(sumR2ball1gray/(468*637)));
rmsr2Ball1gray

diffg2ball1gray = double(Final2canballGgray-FinalcorrballGgray);
diffgS2ball1gray = diffg2ball1gray.^2;
sumg2ball1gray = sum(diffgS2ball1gray);
rmsg2Ball1gray = sqrt(sum(sumg2ball1gray/(468*637)));
rmsg2Ball1gray

%______BALL_END____%

%______Blocks____%
%B)correction blocks and canonical
corrblocksRgray = blocksCorr2(:,:,1)./((blocksCorr2(:,:,1)+blocksCorr2(:,:,2)+blocksCorr2(:,:,3)));
corrblocksGgray = blocksCorr2(:,:,2)./((blocksCorr2(:,:,1)+blocksCorr2(:,:,2)+blocksCorr2(:,:,3)));

newcorrblocksRgray = [];
newcorrblocksGgray = [];

newcanblocksRgray = [];
newcanblocksGgray = [];

for i = 1:637
    for j = 1:468
        if((blocksCorr2(j,i,1)+blocksCorr2(j,i,2)+blocksCorr2(j,i,3)) >= 10)
            newcorrblocksRgray = [newcorrblocksRgray; blocksCorr2(j,i,1)./((blocksCorr2(j,i,1)+blocksCorr2(j,i,2)+blocksCorr2(j,i,3)))];
            newcorrblocksGgray = [newcorrblocksGgray; blocksCorr2(j,i,2)./((blocksCorr2(j,i,1)+blocksCorr2(j,i,2)+blocksCorr2(j,i,3)))];
            newcanblocksRgray = [newcanblocksRgray; blocksCan(j,i,1)./((blocksCan(j,i,1)+blocksCan(j,i,2)+blocksCan(j,i,3)))];
            newcanblocksGgray = [newcanblocksGgray; blocksCan(j,i,2)./((blocksCan(j,i,1)+blocksCan(j,i,2)+blocksCan(j,i,3)))];
        
        end
    end
end
FinalcorrblocksRgray = newcorrblocksRgray;
FinalcorrblocksGgray = newcorrblocksGgray;
Final2canblocksRgray = newcanblocksRgray;
Final2canblocksGgray = newcanblocksGgray;

diffR2blocks1gray = double(Final2canblocksRgray-FinalcorrblocksRgray);
diffRS2blocks1gray = diffR2blocks1gray.^2;
sumR2blocks1gray = sum(diffRS2blocks1gray);
rmsr2Blocks1gray = sqrt(sum(sumR2blocks1gray/(468*637)));
rmsr2Blocks1gray

diffg2blocks1gray = double(Final2canblocksGgray-FinalcorrblocksGgray);
diffgS2blocks1gray = diffg2blocks1gray.^2;
sumg2blocks1gray = sum(diffgS2blocks1gray);
rmsg2Blocks1gray = sqrt(sum(sumg2blocks1gray/(468*637)));
rmsg2Blocks1gray

%______Blocks_END____%


%which method is working better________MAXRGB BECAUSE GRAY WORLD LOOKS LIKE OVEREXPOSED FILM (LIGHT)_________________________________%%



end