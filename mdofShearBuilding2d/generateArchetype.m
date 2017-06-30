% generateArchetype.m

function bldg = generateArchetype(nStories,options)

bldg = mdofShearBuilding2d(nStories);

bldg.verbose = options.verbose;
bldg.echoOpenSeesOutput = options.echoOpenSeesOutput;
bldg.deleteFilesAfterAnalysis = options.deleteFilesAfterAnalysis;

bldg.units       = options.units;
bldg.g           = options.g;

bldg.storyHeight = ones(1,nStories)*options.storyHeight;
bldg.storyHeight(1) = options.firstHeight;

storyDL   = ones(1,nStories)*options.storyDL;
storyDL(end) = options.roofDL;
storyArea = ones(1,nStories)*options.storyArea;
bldg.storyMass   = (storyDL .* storyArea)/options.g;

bldg.pushover_stepSize = options.pushover_stepSize;
bldg.pushover_maxDrift = options.pushover_maxDrift;

bldg.damping_ModeA  = options.damping_ModeA;
bldg.damping_ModeB  = options.damping_ModeB;
bldg.damping_RatioA = options.damping_RatioA;
bldg.damping_RatioB = options.damping_RatioB;

bldg.optionsPushover        = options.pushover;
bldg.optionsResponseHistory = options.responseHistory;
bldg.optionsIDA             = options.IDA;

bldg.seismicDesignCategory  = options.seismicDesignCategory;
bldg.respModCoeff           = options.respModCoeff;
bldg.deflAmplFact           = options.deflAmplFact;
bldg.overstrengthFactor     = options.overstrengthFactor;
bldg.impFactor              = options.impFactor;

end
