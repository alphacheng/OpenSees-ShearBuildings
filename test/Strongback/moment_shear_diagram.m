%########################### moment_shear_diagram.m ###########################%
%                                                                              %
% Script for generating shear and moment diagrams.                             %
%                                                                              %
%                                                                              %
%                                                                              %
%                                                                              %
%                                                                              %
%##############################################################################%

% clear all; close all; clc; %#ok<CLALL>

function moment_shear_diagram(nStories, EI)

%################################# Definition #################################%
% nStories = 6;
bldg = Strongback(nStories);
bldg.seismicDesignCategory = 'Dmax';

%----------------------------- Units and constants ----------------------------%
bldg.units.force = 'kip';
bldg.units.mass  = 'kslug';
bldg.units.length= 'ft';
bldg.units.time  = 'sec';

bldg.g = 32.2;       % Acceleration due to gravity

bldg.seismicDesignCategory = 'Dmax';
bldg.respModCoeff = 8;
bldg.deflAmplFact = 8;
bldg.overstrengthFactor = 3;
bldg.impFactor = 1;

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

springGivens.as       =  0.03;  % strain hardening ratio
springGivens.Lambda_S = 10.00;  % Cyclic deterioration parameter - strength
springGivens.Lambda_K = 10.00;  % Cyclic deterioration parameter - stiffness
springGivens.c_S      =  1.00;  % rate of deterioration - strength
springGivens.c_K      =  1.00;  % rate of deterioration - stiffness
springGivens.Res      =  0.30;  % residual strength ratio (relative to yield)
springGivens.D        =  1.00;  % rate of cyclic deterioration
springGivens.nFactor  =  0.00;  % elastic stiffness amplification factor
springGivens.C_yc     =  0.80;  % ratio of yield strength to capping strength
springGivens.C_upc    = 20.00;  % ratio of ultimate deflection to u_y + u_p + u_pc
springGivens.ad       =  0.10;  % deterioration stiffness ratio -- higher values mean faster deterioration
springGivens.includePDelta = false;

springGivens.stiffnessSafety = 1.0;
springGivens.strengthSafety  = 1.0;

springGivens.enforceMinimumStiffness = false;
springGivens.enforceMinimumStrength = false;
springGivens.minimumRatio = 0.7;

ELF = equivalentLateralForceAnalysis(bldg);
spring = bldg.springDesign(ELF, springGivens);
bldg.storySpringDefinition = {spring.definition}';

targetTrussDeformation = 0.002;  % Ratio of story height
bldg.storyTrussDefinition = cell(nStories,1);
trussModulus = (cumsum(bldg.storyMass,'reverse')*bldg.g)/targetTrussDeformation;
for i = 1:nStories
    bldg.storyTrussDefinition{i} = sprintf('uniaxialMaterial Elastic %i %g',i+10,trussModulus(i));
end

% EI = 1e6;
bldg.strongbackDefinition.Area    = 1;
bldg.strongbackDefinition.Modulus = sqrt(EI);
bldg.strongbackDefinition.Inertia = sqrt(EI);

%----------------------------------- Options ----------------------------------%
bldg.echoOpenSeesOutput = false;
bldg.deleteFilesAfterAnalysis = false;
% analysisType = 'responseHistory';

paths = fieldnames(pathOf);
for i = 1:length(paths)
    bldg.pathOf.(paths{i}) = pathOf.(paths{i});
end
bldg.pathOf.tclfunctions = '/home/petertalley/Github/OpenSees-ShearBuildings/lib';

bldg.optionsPushover.maxDrift = sum(bldg.storyHeight);
bldg.optionsPushover.test.print = 0;

bldg.optionsIDA.dt_min = 0;
bldg.optionsIDA.nMotions = 14;

gm_mat = '../ground_motions.mat';
gmIndex_start = 1;
gmIndex_end = min(bldg.optionsIDA.nMotions + gmIndex_start - 1, 44);

%################################ Run analysis ################################%

load(gm_mat)

for analysisType = {'responseHistory', 'pushover'}
    momentShearEnvelopes(bldg, analysisType{:}, ground_motions, [gmIndex_start, gmIndex_end]);
end

end
