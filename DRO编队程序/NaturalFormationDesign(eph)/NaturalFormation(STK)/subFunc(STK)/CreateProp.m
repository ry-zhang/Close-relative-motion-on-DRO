function CreateProp(scenario,auxSFF)
% 创建一个新的Propagator，用户可自定义
astroComp = scenario.ComponentDirectory.GetComponents('eComponentAstrogator');
propComp = astroComp.GetFolder('Propagators');
% 复制HPOP,重命名
propComp.DuplicateComponent('Earth HPOP Default v10', auxSFF.PropagatorName);
hpopComp = propComp.Item(auxSFF.PropagatorName);
propFunComp = hpopComp.PropagatorFunctions;
try
    propFunComp.Remove('Jacchia-Roberts');
catch
end
propMoon = propFunComp.Item('Moon');
propMoon.SetModeType('eVAThirdBodyModeGravityField');
propMoon.Mode.Degree = 8;
propMoon.Mode.Order = 8;
propFunComp.Add('Mercury');
propFunComp.Add('Venus');
propFunComp.Add('Mars');
propFunComp.Add('Jupiter System');
propFunComp.Add('Saturn System');

% %% 设置地球引力
% propCentralBody = propFunComp.Item('Gravitational Force');
% propCentralBody.GravityFilename = 'STKData\CentralBodies\Earth\GGM01C.grv';
% propCentralBody.Degree = 70;
% propCentralBody.Order = 70;
% propCentralBody.UseSecularVariations = 0;
% % Ocean Tides
% propCentralBody.UseOceanTides = 1;
% propCentralBody.OceanTideMaxDegree = 30;
% propCentralBody.OceanTideMaxOrder = 30;
% propCentralBody.SolidTideType = 'eSolidTidePermanent'; % 'eSolidTideFull'
% % propCentralBody.TruncateSolidTides = 1; % only valid for eSolidTideFull type
% % propCentralBody.IncludeTimeDependentSolidTides = 0; % only valid for eSolidTideFull type
% 
% %% 设置月球引力
% propMoon = propFunComp.Item('Moon');
% propMoon.SetModeType('eVAThirdBodyModeGravityField');
% propMoon.EphemSource = 'eVAEphemSourceDEFile'; % 'eVAEphemSourceCbFile'
% propMoon.Mode.GravityFilename = 'STKData\CentralBodies\Moon\GL0660B.grv';
% propMoon.Mode.Degree = 20;
% propMoon.Mode.Order = 20;
% 
% %% 设置太阳引力
% propSun = propFunComp.Item('Sun');
% propSun.SetModeType('eVAThirdBodyModePointMass');
% propSun.EphemSource = 'eVAEphemSourceDEFile'; % 'eVAEphemSourceCbFile'
% propSun.EphemSource = 'eVAEphemSourceDEFile'; % 'eVAEphemSourceCbFile'
% propSun.Mode.GravSource = 'eVAGravParamSourceDEFile'; % 'eVAGravParamSourceCbFile'
% % propSun.Mode.Mu = 132712440041.939;
% 
% %% 设置太阳光压
% propSRP = propFunComp.Item('Spherical SRP');
% propSRP.EclipsingBodies.RemoveAll;
% % propSRP.EclipsingBodies.Add('Moon');
% % propSRP.UseSunCbFileValues = 0;
% % propSRP.Luminosity = 3.839e+26;
% % propSRP.MeanFlux = 1365.07787607266;
% 
% %% 设置大气密度
% try
%     propFunComp.Remove('Jacchia-Roberts');
% catch
% end
% propFunComp.Add('NRLMSISE 2000');
% propDensity = propFunComp.Item('NRLMSISE 2000');
% propDensity.AtmosDataSource = 'eVAAtmosDataSourceFile'; %'eVAAtmosDataSourceConstant'
% propDensity.AtmosDataFilename = 'C:\ProgramData\AGI\STK 11 (x64)\DynamicEarthData\SpaceWeather-v1.2.txt'; % 文件路径
% propDensity.AtmosDataGeoMagneticFluxSource = 'eVAGeoMagneticFluxSourceAp'; %'eVAGeoMagneticFluxSourceKp'

