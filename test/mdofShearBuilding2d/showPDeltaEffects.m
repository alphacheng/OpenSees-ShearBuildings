% showPDeltaEffects.m

test_options
nStories = 3;

springGivens.includePDelta = true;
bldg_with_pdelta = generateArchetype(nStories,options);
ELF = bldg_with_pdelta.ELFanalysis();
spring_with_pdelta = bldg_with_pdelta.springDesign(ELF,springGivens);
storySpringDefinition = cell(nStories,1);
for i = 1:nStories
    storySpringDefinition{i} = spring_with_pdelta(i).definition;
end
bldg_with_pdelta.storySpringDefinition = storySpringDefinition;

fig = bldg_with_pdelta.plotBackboneCurves(spring_with_pdelta);
grarray = fig.Children;
ax = grarray(2);
lines = ax.Children;
data  = cell(nStories,2);
color = cell(nStories,1);
for i = 1:nStories
    data{i,1} = lines(i).XData;
    data{i,2} = lines(i).YData;
    color{i}  = lines(i).Color;
end
close(fig)

springGivens.includePDelta = false;
bldg_without_pdelta = generateArchetype(nStories,options);
ELF = bldg_without_pdelta.ELFanalysis();
spring_without_pdelta = bldg_without_pdelta.springDesign(ELF,springGivens);
storySpringDefinition = cell(nStories,1);
for i = 1:nStories
    storySpringDefinition{i} = spring_without_pdelta(i).definition;
end
bldg_without_pdelta.storySpringDefinition = storySpringDefinition;

fig2 = bldg_without_pdelta.plotBackboneCurves(spring_without_pdelta);
grarray = fig2.Children;
ax = grarray(2);
for i = nStories:-1:1
    plot(ax,data{i,1},data{i,2},'--','Color',color{i})
end
