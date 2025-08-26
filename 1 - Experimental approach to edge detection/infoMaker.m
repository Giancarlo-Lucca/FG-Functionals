
clear Experiment;

externalHD = 0;

if (strcmp(OS_MODE,'linux'))
    ROOT_FOLDER='/home/user/';
    if (~exist(ROOT_FOLDER,'dir'))
    	ROOT_FOLDER='/home/user/Research/';
    end
elseif (strcmp(OS_MODE,'mac'))
    externalHD = 0;
    ROOT_FOLDER='/Users/user/Research/';
elseif(strcmp(OS_MODE,'win') || strcmp(OS_MODE,'windows'))
    ROOT_FOLDER='D:/Research/';
else
    error('Error at infoMaker.m> Unknown OS_MODE [%s]\n',OS_MODE);
end


%
% Paths
%

Experiment.set='originalRGB';
Experiment.experimentName = 'FG-func'

Experiment.databaseName = 'BSDS5full';
Experiment.originalFolder = 'originalRGB';%originalGS, originalRGB


% Experiment.orgDir=[ROOT_FOLDER 'Resources/Images/BSDS5full/originalGS/' Experiment.set '/'];
% Experiment.gtBnDir=[ROOT_FOLDER 'Resources/Images/BSDS5full/gtruthBN/'];

if (strcmp(Experiment.databaseName,'BSDS5full'))
    Experiment.orgDir=[ROOT_FOLDER 'Resources/Images/BSDS5full/' Experiment.originalFolder '/' Experiment.set '/'];
    Experiment.gtBnDir=[ROOT_FOLDER 'Resources/Images/BSDS5full/gtruthBN/' Experiment.set '/'];
    Experiment.researchDataDir=[ROOT_FOLDER 'ResearchData/' Experiment.experimentName '/'];
elseif (strcmp(Experiment.databaseName,'BSDS5test'))
    Experiment.orgDir=[ROOT_FOLDER 'Resources/Images/BSDS5test/test/'];
    Experiment.gtBnDir=[ROOT_FOLDER 'Resources/Images/BSDS5test/testGT/'];
    Experiment.researchDataDir=[ROOT_FOLDER 'ResearchData/' Experiment.experimentName '/'];
else
    error('Wrong Experiment.databaseName at infoMaker.m');
end

if externalHD
    Experiment.resultBaseDir=['/Volumes/BACKUP/' 'ResearchData/' Experiment.experimentName '/' Experiment.set];
else
    Experiment.resultBaseDir=[ROOT_FOLDER 'ResearchData/' Experiment.experimentName '/' Experiment.set];
end

Experiment.smDir = [Experiment.resultBaseDir '/smooIm/'];
Experiment.featImDir = [Experiment.resultBaseDir '/featIm/'];
Experiment.bdryDir = [Experiment.resultBaseDir '/bn/'];
Experiment.cpDir = [Experiment.resultBaseDir '/cp/'];
Experiment.texDir = [Experiment.resultBaseDir '/plots/'];


if (~exist(Experiment.smDir,'dir'))
    mkdir(Experiment.smDir);
end
if (~exist(Experiment.featImDir,'dir'))
    mkdir(Experiment.featImDir);
end
if (~exist(Experiment.bdryDir,'dir'))
    mkdir(Experiment.bdryDir);
end
if (~exist(Experiment.cpDir,'dir'))
    mkdir(Experiment.cpDir);
end
if (~exist(Experiment.texDir,'dir'))
    mkdir(Experiment.texDir);
end


%
% Parameters
%

%
% Files
%

% Experiment.RFmodelDir=strcat(ROOT_FOLDER,'research/extern/structured-edges/models');

Experiment.imagesFrom=1;
Experiment.imagesTo=200;

Experiment.smPrefix='sm-';%Smoothing
Experiment.ftPrefix='ft-';%Feature image
Experiment.bdryPrefix='bdry-';%Boundary image
Experiment.cpPrefix='cp-';%Comparison

Experiment.imageExt='png';
if strcmp(Experiment.originalFolder,'originalRGB')
    Experiment.orgImagesExt='jpg';
else
    Experiment.orgImagesExt='png';%jpg
end
Experiment.bnPrefix='bn';
Experiment.cpPrefix='cp';
Experiment.dataExt='mat';

Experiment.forceProcessMaker=0;
Experiment.forceSmMaker=0;
Experiment.forceFtMaker=0;
Experiment.forceBnMaker=0;
Experiment.forceCpMaker=0;
Experiment.saveBnImages=0;

Experiment.forceFtNorm = 0;

Experiment.writeImages=0;
Experiment.writeImagesDir=0;
if Experiment.writeImages
    Experiment.writeImagesDir=0;
end

%
% Smoothing
%
1;%
Experiment.smoothingMethod = {{'gauss', 1},...
                              {'gauss', 2},...
                              };
Experiment.gauss.sigma=[2];%2
Experiment.mshift.spatialSupport=5;
Experiment.mshift.tonalSupport=25;
Experiment.mshift.stopCondition=0.05;

Experiment.grav.iterations = 50;%10, 30, 50
Experiment.grav.minDistInfFactor = 0.02;%0.05, 0.02
Experiment.grav.gConst = 0.05;%0.05 0.5
Experiment.grav.colorFactor = 70;%20, 70
Experiment.grav.colorMetric = 'euc';
Experiment.grav.posMetric = 'euc';

%
% Choquet
%


Experiment.config.featureMethod = {'canny','FGsugeno'}; % 'dXCF', 'dCF', 'dXCF1F2', 'sugeno', 'FGsugeno', 'canny','fuzzyM','ged','SF'

if strcmp(Experiment.config.featureMethod{1},'dCF')
    Experiment.config.feat.dCF.Dis = {'d0'}; % 'd0','d1','d2','d3','d4','d5'
    Experiment.config.feat.dCF.Gen = {'choquet'};%{'choquet','C_F','O_B','F_BPC','hamacher'};%{'T_M','choquet','T_L','hamacher','O_B','O_mM','O_alpha','O_Div','GM','HM','Sin','O_RS','C_F','C_L','F_GL','F_BPC','F_BD1','FNA','FNA2'};% 'hamacher','FNA'
elseif strcmp(Experiment.config.featureMethod{1},'dXCF')
    Experiment.config.feat.dCF.Dis = {'d0','d1','d2','d3','d4','d5'}; % 'd0','d1','d2','d3','d4','d5'
    Experiment.config.feat.dCF.Gen = {'T_M','prod','lukasiewicz','hamacker','OB','OmM','ODiv','CF','CL'};%{'T_M','prod','lukasiewicz','hamacker','OB','OmM','ODiv','CF','CL'}%{'choquet','C_F','O_B','F_BPC','hamacher'};%{'T_M','choquet','T_L','hamacher','O_B','O_mM','O_alpha','O_Div','GM','HM','Sin','O_RS','C_F','C_L','F_GL','F_BPC','F_BD1','FNA','FNA2'};% 'hamacher','FNA'
elseif strcmp(Experiment.config.featureMethod{1},'dXCF1F2')
    Experiment.config.feat.dCF.Dis = {'d0'}; % 'd0','d1','d2','d3','d4','d5'
    Experiment.config.feat.dCF.Gen = {{'prod', 'prod'},{'T_M', 'T_M '},{'hamacker', 'hamacker'},{'GL', 'T_M'},{'GM', 'FBPC'}};%{'T_M','prod','T_L','hamacher','O_B','O_mM','O_alpha','O_Div','GM','HM','Sin','O_RS','C_F','C_L','F_GL','F_BPC','F_BD1','FNA','FNA2'};% 'hamacher','FNA'
elseif strcmp(Experiment.config.featureMethod{1},'sugeno')
    Experiment.config.feat.dCF.Dis = {''};
    Experiment.config.feat.dCF.Gen = {''};
elseif strcmp(Experiment.config.featureMethod{1},'FGsugeno')
    Experiment.config.feat.dCF.Dis = {''};
    Experiment.config.feat.dCF.Gen = {{'sum', 'hamacker'},...
                                      {'sum', 'T_M'},...
                                      {'sum', 'prod'},...
                                      {'max', 'hamacker'},...
                                      {'max', 'T_M'},...
                                      {'max', 'prod'},...
                                      {'Fi', 'hamacker'},...
                                      {'Fi', 'T_M'},...
                                      {'Fi', 'prod'},...
                                      {'DivFi', 'hamacker'},...
                                      {'DivFi', 'T_M'},...
                                      {'DivFi', 'prod'}};
end
Experiment.config.feat.dCF.F1 = {'prod'};
Experiment.config.feat.dCF.F2 = {'prod'};
Experiment.config.feat.canny.sigma =  2.25;
Experiment.config.feat.fuzzyM = {'T_nM','I_KD','SS'};
Experiment.config.feat.ged.F = {'S_M','S_P'};
Experiment.config.feat.ged.wSize = 3;
Experiment.config.feat.ged.kmode = 'global';

Experiment.config.NumClassicMethods = 5  ;

Experiment.tamW=[3];% 3, 5, 7
Experiment.colorAgg = 'max';


Experiment.config.measure = {'power'};
Experiment.params.measure.power.useSelected_q = 0;
Experiment.params.measure.power.q = 1;%0.1:0.1:1;
Experiment.params.measure.power.selected_q = [1];
if Experiment.params.measure.power.useSelected_q
    Experiment.params.measure.power.q = Experiment.params.measure.power.selected_q;
end
Experiment.params.measure.power.type = {''};%'adapVec-T_Max'};%'adaptative-max', 'adaptative-min', 'adaptative-mean', ...'adapVec-hamacker'};Experiment.params.measure.tam = Experiment.tamW.^2 - 1;

if (ismember('power',Experiment.config.measure))

    for idxQ = 1:length(Experiment.params.measure.power.q)
        mName{1,idxQ} = ['power-' sigma2name(Experiment.params.measure.power.q(idxQ))];
    end
    
    idxPower = find(not(cellfun('isempty', strfind(Experiment.config.measure, 'power'))));

    Experiment.config.measure = [Experiment.config.measure(idxPower) Experiment.config.measure(1:end ~= idxPower)];
    Experiment.config.measureComplete = [mName Experiment.config.measure(1:end ~= idxPower)];
else
    Experiment.config.measureComplete = Experiment.config.measure;
end

if any(ismember(Experiment.config.featureMethod,{'dCF'}))
    Experiment.numRes = length(Experiment.smoothingMethod)*length(Experiment.params.measure.power.q)*length(Experiment.config.feat.dCF.Dis)*length(Experiment.config.feat.dCF.Gen);
end
Experiment.p=inf;%0.5; % 0.3

%
% ColorMap for feateure output
%

Experiment.dtDiffColorMap=createColorMap([0.9, 0.9, 0.9],...
                                         [0.14, 0.12, 0.1],...
                                         [0.92, 0.37, 0],...
                                         256);

Experiment.colorMap = makeColorMap([0; 70;  100],...
                                   ['#010101'; '#085a3f';  '#ede6df'],...
                                   256, 1);

%
% Matching parameters
%


Experiment.matching = 'EJMbCM-F'; % {'DistbCM-F','DilbCM-F','EJMbCM-F','csaBCM-F'}
Experiment.matchingTolerance = 0.025;%25





