function CreateMCRLVLHframe(auxSAT)
% 建立月心地月旋转系中的LVLH坐标系

for ii_index = 1:length(auxSAT)
    SatVgt = auxSAT(ii_index).satA.sat.Vgt;
    axeEMRLVLH = SatVgt.Axes.Factory.Create(...
        'MCEMRLVLH', 'Earth-Moon Rotation LVLH', 'eCrdnAxesTypeTrajectory');
    axeEMRLVLH.ReferenceSystem.SetPath('CentralBody/Moon MCEMR');
    axeEMRLVLH.TrajectoryAxesType = 'eCrdnTrajectoryAxesLVLH';
    
    MCEMRLVLH_SC = SatVgt.Systems.Factory.Create('MCEMRLVLH_SC','Earth-Moon Rotation LVLH','eCrdnSystemTypeAssembled');
    MCEMRLVLH_SC.OriginPoint.SetPoint(SatVgt.Points.Item('Center'));
    MCEMRLVLH_SC.ReferenceAxes.SetAxes(axeEMRLVLH);
end