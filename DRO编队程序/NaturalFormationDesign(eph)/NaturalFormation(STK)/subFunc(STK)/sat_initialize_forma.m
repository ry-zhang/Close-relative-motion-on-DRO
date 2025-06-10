function sat_Interf = sat_initialize_forma(scenario,rv0,t0_sec,dt1_sec,dt2_sec,satName,PropagatorName)

sat_Interf.sat = scenario.Children.New('eSatellite',satName);
sat_Interf.sat.SetPropagatorType('ePropagatorAstrogator');
satMSC = sat_Interf.sat.Propagator.MainSequence;
satMSC.RemoveAll;

% insert segments
sat_Interf.InitialState = satMSC.Insert('eVASegmentTypeInitialstate','InitialState','-');
sat_Interf.Propagate1 = satMSC.Insert('eVASegmentTypePropagate','Propagate1','-');
sat_Interf.Propagate2 = satMSC.Insert('eVASegmentTypePropagate','Propagate2','-');

% set the initial state
sat_Interf.InitialState.OrbitEpoch = t0_sec;
sat_Interf.InitialState.CoordSystemName = 'CentralBody/Earth J2000';
sat_Interf.InitialState.SetElementType('eVAElementTypeCartesian');
set(sat_Interf.InitialState.Element, 'X', rv0(1),'Y',rv0(2),...
    'Z',rv0(3),'Vx',rv0(4),'Vy',rv0(5),'Vz',rv0(6));

% set the propergate1 segements
sat_Interf.Propagate1.PropagatorName = PropagatorName;
sat_Interf.Propagate1.StoppingConditions.Item('Duration').Properties.Trip = dt1_sec;
% set the propergate2 segements,to achieve the theta_tar
sat_Interf.Propagate2.PropagatorName = PropagatorName;
sat_Interf.Propagate2.StoppingConditions.Item('Duration').Properties.Trip = dt2_sec;

% run the astrogator
sat_Interf.sat.Propagator.RunMCS;