clear ;close all ;clc;
warning off
imagePath='.\data\911\1.jpg';
[numOfGrains,numOfBud] = calcGrainNum(imagePath);
germinationRate=numOfBud/numOfGrains;