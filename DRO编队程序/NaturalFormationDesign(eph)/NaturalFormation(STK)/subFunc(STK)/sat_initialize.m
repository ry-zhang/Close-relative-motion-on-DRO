function sat_Interf = sat_initialize(scenario,rv0,CoordSystem,t0_sec,dt_sec,satName,PropagatorName)

sat_Interf.sat = scenario.Children.New('eSatellite',satName);
sat_Interf.sat.SetPropagatorType('ePropagatorAstrogator');
satMSC = sat_Interf.sat.Propagator.MainSequence;
satMSC.RemoveAll;

% insert segments
sat_Interf.InitialState = satMSC.Insert('eVASegmentTypeInitialstate','InitialState','-');
sat_Interf.Propagate = satMSC.Insert('eVASegmentTypePropagate','Propagate','-');

% set the initial state
sat_Interf.InitialState.OrbitEpoch = t0_sec;
sat_Interf.InitialState.CoordSystemName = CoordSystem;
sat_Interf.InitialState.SetElementType('eVAElementTypeCartesian');
set(sat_Interf.InitialState.Element, 'X', rv0(1),'Y',rv0(2),...
    'Z',rv0(3),'Vx',rv0(4),'Vy',rv0(5),'Vz',rv0(6));

% set the propergate segements
sat_Interf.Propagate.PropagatorName = PropagatorName;
sat_Interf.Propagate.StoppingConditions.Item('Duration').Properties.Trip = dt_sec;

% run the astrogator
sat_Interf.sat.Propagator.RunMCS;