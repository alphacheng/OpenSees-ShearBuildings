%############################# test_Strongback.m ##############################%
%                                                                              %
% Script for testing the functioning of the Strongback model.                  %
%                                                                              %
%                                                                              %
%                                                                              %
%                                                                              %
%                                                                              %
%##############################################################################%

clear all; close all; clc; %#ok<CLALL>

%################################# Definition #################################%
nStories = 3;
bldg = Strongback(nStories);
bldg.seismicDesignCategory = 'Dmax';

%----------------------------- Units and constants ----------------------------%
bldg.units.force = 'kip';
bldg.units.mass  = 'kslug';
bldg.units.length= 'ft';
bldg.units.time  = 'sec';

bldg.g = 32.2;       % Acceleration due to gravity

%------------------------------ Story definition ------------------------------%
storyHeight = 15;
firstHeight = 20;
storyDL     = 0.080;
roofDL      = 0.030;
storyArea   = 90*90;

bldg.storyHeight    = ones(1,nStories)*storyHeight;
bldg.storyHeight(1) = firstHeight;

storyDL        = ones(1,nStories)*storyDL;
storyDL(end)   = roofDL;
bldg.storyMass = (storyDL*storyArea)/bldg.g;

bldg.storySpringDefinition = {
    'uniaxialMaterial Bilin 1 3847.50 0.03 0.03 461.700 -461.700 10 10 10 10 1 1 1 1 1.00 1.00 1.500 1.500 0.3 0.3 52.4 52.4 1 1 0'
    'uniaxialMaterial Bilin 2 3738.81 0.03 0.03 336.493 -336.493 10 10 10 10 1 1 1 1 0.75 0.75 1.125 1.125 0.3 0.3 39.3 39.3 1 1 0'
    'uniaxialMaterial Bilin 3 1304.24 0.03 0.03 117.381 -117.381 10 10 10 10 1 1 1 1 0.75 0.75 1.125 1.125 0.3 0.3 39.3 39.3 1 1 0'
};
bldg.storyTrussDefinition = {
    'uniaxialMaterial Elastic 11 769500'
    'uniaxialMaterial Elastic 12 445500'
    'uniaxialMaterial Elastic 13 121500'
};
bldg.strongbackDefinition.Area    = 1;
bldg.strongbackDefinition.Modulus = 1e6;
bldg.strongbackDefinition.Inertia = 1e6;

%----------------------------------- Options ----------------------------------%
bldg.echoOpenSeesOutput = false;
bldg.deleteFilesAfterAnalysis = false;

bldg.optionsPushover.maxDrift = sum(bldg.storyHeight);
bldg.optionsPushover.test.print = 0;

%################################# Run tests! #################################%

F = bldg.pushoverForceDistribution();
pushover = bldg.pushover(F,'TargetPostPeakRatio',0.75);
