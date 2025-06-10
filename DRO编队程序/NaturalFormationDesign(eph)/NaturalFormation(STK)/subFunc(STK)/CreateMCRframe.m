function CreateMCRframe(scenario,root)
% 建立月心地月旋转系

%% 建立Moon
planetMoon = scenario.Children.New('ePlanet', 'Moon');
planetMoon.CommonTasks.SetPositionSourceCentralBody(...
    'Moon', 'eEphemJPLDE');

%% 建立MCEMR坐标系
MoonVgt = root.CentralBodies.Moon.Vgt;
axeEMR = MoonVgt.Axes.Factory.Create(...
    'EMR', 'Earth-Moon Rotation', 'eCrdnAxesTypeAlignedAndConstrained');
axeEMR.AlignmentDirection.AssignXYZ(-1, 0, 0);
axeEMR.AlignmentReferenceVector.SetVector(MoonVgt.Vectors.Item('Earth'));
axeEMR.ConstraintDirection.AssignXYZ(0, 0, 1);
axeEMR.ConstraintReferenceVector.SetVector(MoonVgt.Vectors.Item('Orbit_Normal'));
MCEMR = MoonVgt.Systems.Factory.Create('MCEMR','Moon-centered EMR','eCrdnSystemTypeAssembled');
MCEMR.OriginPoint.SetPoint(MoonVgt.Points.Item('Center'));
MCEMR.ReferenceAxes.SetAxes(axeEMR);