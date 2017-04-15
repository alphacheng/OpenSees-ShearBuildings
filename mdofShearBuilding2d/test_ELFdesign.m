% test_ELFdesign.m
% Units: kslug, kip, ft, sec

clear; close all; clc;
load('ground_motions.mat');

%% Define Building
nStories = 3;
bldg = mdofShearBuilding2d(nStories);
bldg.echoOpenSeesOutput = false;
bldg.deleteFilesAfterAnalysis = false;
bldg.g = 32.2;

bldg.seismicDesignCategory = 'Dmax';
bldg.responseModificationCoefficient = 8;
bldg.deflectionAmplificationFactor = 5.5;
bldg.overstrengthFactor = 3;
bldg.importanceFactor = 1;

bldg.storyHeight = [20 15 15];
storyDL = [0.080 0.080 0.030];
storyArea = 90*90*ones(1,nStories);

bldg.storyMass = (storyDL .* storyArea)/bldg.g;

results = bldg.ELFdesign();
