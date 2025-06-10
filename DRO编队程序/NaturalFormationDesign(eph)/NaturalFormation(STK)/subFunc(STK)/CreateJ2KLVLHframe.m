function CreateJ2KLVLHframe(auxSAT)
% 建立地心J2000中的LVLH坐标系

for ii_index = 1:length(auxSAT)
    SatVgt = auxSAT(ii_index).satA.sat.Vgt;
    axeJ2KLVLH = SatVgt.Axes.Factory.Create(...
        'J2KLVLH', 'J2000 LVLH', 'eCrdnAxesTypeTrajectory');
    axeJ2KLVLH.ReferenceSystem.SetPath('CentralBody/Earth J2000');
    axeJ2KLVLH.TrajectoryAxesType = 'eCrdnTrajectoryAxesLVLH';
    
    J2KLVLH_SC = SatVgt.Systems.Factory.Create('J2KLVLH_SC','Earth-centered J2000 LVLH','eCrdnSystemTypeAssembled');
    J2KLVLH_SC.OriginPoint.SetPoint(SatVgt.Points.Item('Center'));
    J2KLVLH_SC.ReferenceAxes.SetAxes(axeJ2KLVLH);
end