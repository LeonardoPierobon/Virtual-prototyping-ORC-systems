within ;
package Design "Package for the component design"
  import Modelica.Math.*;
  import Modelica.SIunits.*;
  import Modelica.Constants.*;
  package Components "Library with the design of the components "
    package HEX "Heat exchangers"

      model shell_and_tube
        "Shell and tube heat exchanger where the hot fluid flows on the shell and enters from the top. The cold fluid enters at the bottom."
        //THE WORKING FLUIDS
        replaceable package Medium_s = Media.OneRandomOrganicFluid constrainedby
          Modelica.Media.Interfaces.PartialMedium "Medium model shell" annotation(choicesAllMatching = true);
        replaceable package Medium_t = Media.OneRandomOrganicFluid constrainedby
          Modelica.Media.Interfaces.PartialMedium "Medium model tubes" annotation(choicesAllMatching = true);
        replaceable Materials.void_material Material_t "Material model shell"  annotation(choicesAllMatching = true, Dialog(tab="Geometry"));
        replaceable Materials.void_material Material_s "Material model shell"  annotation(choicesAllMatching = true, Dialog(tab="Geometry"));
        parameter Boolean use_dp = false
          "Use pressure drop for outlet pressure if true";
        parameter Boolean offdesign = false "Off design mode on if true" annotation (Dialog(tab="Off-design"));
        parameter Modelica.SIunits.Length N_tubes_od = 1
          "Number of tubes fixed for off-design"
                                               annotation (Dialog(tab="Off-design"));

        //GEOMETRY OF THE HEAT EXCHANGER
        parameter Modelica.SIunits.Length Dhyd "Hydraulic diameter" annotation (Dialog(tab="Geometry"));
        parameter Modelica.SIunits.Length thick_t "Tube thickness" annotation (Dialog(tab="Geometry"));
        parameter Modelica.SIunits.Length thick_s "Shell thickness" annotation (Dialog(tab="Geometry"));
        parameter Modelica.SIunits.Length l "Tube lenght" annotation (Dialog(tab="Geometry"));
        parameter Real pitch_f "Tube pitch" annotation (Dialog(tab="Geometry"));
        parameter Integer layout "Tube layout, 1 = triangular, 2 = squared" annotation (Dialog(tab="Geometry"));
        parameter Integer N_passes = 2 "Number of tube passes" annotation (Dialog(tab="Geometry"));
        parameter Integer N_baffles = 4 "Number of baffles" annotation (Dialog(tab="Geometry"));
        parameter Real pin_s = 1
          "Pin for the heat flow in the shell -1 cold fluid else hot";
        parameter Real pin_t = -1
          "Pin for the heat flow in the tubes -1 cold fluid else hot";
        parameter Integer N_baffles_d = N_baffles + 1
          "The number of discretized baffles";
        final parameter Integer Ncell = N_baffles_d*N_passes
          "Number of cell elements";
        parameter Modelica.SIunits.SpecificEnthalpy h_s_in_start
          "Inlet specific enthalpy start value hot side" annotation (Dialog(tab="Start"));
        parameter Modelica.SIunits.SpecificEnthalpy h_s_out_start
          "Outlet specific enthalpy start value hot side" annotation (Dialog(tab="Start"));
        parameter Modelica.SIunits.SpecificEnthalpy h_t_in_start
          "Inlet specific enthalpy start value cold side" annotation (Dialog(tab="Start"));
        parameter Modelica.SIunits.SpecificEnthalpy h_t_out_start
          "Outlet specific enthalpy start value cold side" annotation (Dialog(tab="Start"));
        parameter Real N_tubes_start
          "Number of tubes in the bundle start value"                            annotation (Dialog(tab="Start"));
        parameter Modelica.SIunits.CoefficientOfHeatTransfer ht_t_f1 = 3e3
          "Tube fouling heat transfer coefficient" annotation (Dialog(tab="Geometry"));
        parameter Modelica.SIunits.CoefficientOfHeatTransfer ht_s_f1 = 5e3
          "Shell fouling heat transfer coefficient" annotation (Dialog(tab="Geometry"));
        Modelica.SIunits.Length d_s "Shell diameter";
        Modelica.SIunits.Length l_b "Baffle lenght";
        Real N_tubes(start = N_tubes_start, fixed = true)
          "Number of tubes in the bundle";
        Real N_t_p_p "Number of tubes per pass";
        Real bs_f "Baffle spacing in percent";
        Modelica.SIunits.CoefficientOfHeatTransfer U
          "Global heat transfer coefficient";
        Modelica.SIunits.Area A "Heat transfer area";
        Modelica.SIunits.Mass W_dry "Dry weight of the heat exchanger";
        Modelica.SIunits.Mass W_fluids "Weight of the fluids";
        Modelica.SIunits.Mass W_wet "Wet weight of the heat exchanger";
        Real PEC "Purchases equipment cost";
        Modelica.SIunits.HeatFlowRate Q "Heat rate";

        //Heat transfer and pressure drop correlations
        replaceable Heat_transfer.Tubes.Sieder_Tate hT_tube(Dhyd=Dhyd,
        eta_wall=bundle.state[1].eta*ones(Ncell))
          constrainedby Heat_transfer.Tubes.Base_classes.base_ht(redeclare
            package Medium =
                     Medium_t,
          Ncell = Ncell, state = bundle.state, mdot = tube_in.m_flow/N_t_p_p,
          Aflow = bundle.Aflow) annotation(choicesAllMatching = true);

        replaceable Pressure_drops.Tubes.Frank dp_tube(Dhyd = Dhyd, l = l/N_baffles_d,
        heads = 2.5/Ncell)
          constrainedby Pressure_drops.Tubes.Base_classes.base_dp(redeclare
            package Medium =
                     Medium_t,
          Ncell = Ncell, state = bundle.state, mdot = tube_in.m_flow/N_t_p_p,
          Aflow = bundle.Aflow) annotation(choicesAllMatching = true);

        replaceable Heat_transfer.Shell.single_phase_Kern hT_shell(
        Dhyd_o = Dhyd + 2*thick_t, layout = layout, pitch_f = pitch_f)
          constrainedby Heat_transfer.Shell.Base_classes.base_ht(redeclare
            package Medium =
                     Medium_s,
          Ncell = Ncell, state = shell.state, mdot = shell_in.m_flow, Aflow = shell.Aflow) annotation(choicesAllMatching = true);

        replaceable Pressure_drops.Shell.single_phase_Johnston dp_shell(
          l = l/N_baffles_d/N_passes,
          d_s = d_s,
          Dhyd = hT_shell.d_s_eq,
          l_b = l_b,
          eta_wall = shell.state[1].eta*ones(Ncell)) constrainedby
          Pressure_drops.Shell.Base_classes.base_dp(
          redeclare package Medium = Medium_s,
          Ncell = Ncell,
          state = shell.state,
          mdot = shell_in.m_flow,
          Aflow = shell.Aflow) annotation (choicesAllMatching=true);

        //Defining the model for the bundle clearance
        replaceable function bundle_clearance =
        Miscellanea.Shell_clearance.base_clearance annotation(choicesAllMatching = true);

        //Defining the model for the cost
        replaceable function cost = Miscellanea.Cost.base_cost  annotation(choicesAllMatching = true);

        //Definiing the tubes and the shell
        Objects.tube_bundle
        bundle(redeclare package Medium = Medium_t,
        h_in = tube_in.h,
        h_out = tube_out.h,
        h_out_start=h_t_out_start,
        h_in_start=h_t_in_start,
        Ncell = Ncell, Aflow = 0.25*pi*Dhyd^2, mdot = tube_in.m_flow,
        Dhyd = Dhyd, thick = thick_t, lambda = Material_t.lambda,
        rho = Material_t.rho, N_tubes = N_tubes,N_passes = N_passes, layout = layout,
        pin = pin_t, pitch_f = pitch_f);

        Objects.tube shell(
          redeclare package Medium = Medium_s,
          h_in=shell_in.h,
          h_out_start=h_s_out_start,
          h_in_start=h_s_in_start,
          h_out=shell_out.h,
          Ncell=Ncell,
          Aflow=(1 - 1/pitch_f)*d_s*l_b,
          mdot=shell_in.m_flow,
          Dhyd=d_s,
          thick=thick_s,
          lambda=Material_s.lambda,
          rho=Material_s.rho,
          pin=pin_s);

        Miscellanea.check_velocity check_shell(redeclare package Medium = Medium_s,
                    T = shell.state[1].T, umin = min(hT_shell.u),
                    umax = max(hT_shell.u), geometry = "shell", op_p = shell_in.p);

        Miscellanea.check_velocity check_tube(redeclare package Medium = Medium_t,
                    T = bundle.state[1].T, umin = min(hT_tube.u),
                    umax = max(hT_tube.u), geometry = "tube", op_p = tube_in.p);

      protected
         parameter Modelica.SIunits.CoefficientOfHeatTransfer ht_t_f[Ncell]=
          ones(Ncell)*ht_t_f1 "Tube fouling heat transfer coefficient (array)";
         parameter Modelica.SIunits.CoefficientOfHeatTransfer ht_s_f[Ncell]=
         ones(Ncell)*ht_s_f1 "Shell fouling heat transfer coefficient (array)";
         Integer k[Ncell] "Index for the shell to tube topology";
         Modelica.SIunits.ThermalConductance kA_tot[Ncell]
          "Overall heat transfer coefficient (array)";
         Modelica.SIunits.ThermalConductance htA_shell[Ncell]
          "Shell heat transfer coefficient (array)";
         Modelica.SIunits.ThermalConductance htA_tubes[Ncell]
          "Tube heat transfer coefficient (array)";
         Modelica.SIunits.ThermalConductance G_wall[Ncell]
          "Wall heat transfer coefficient (array)";

      public
        Nodes.Node_out shell_out(redeclare package Medium = Medium_s)
          "Outlet node shell side"
                                  annotation (Placement(
              transformation(extent={{-108,-73},{-88,-53}}), iconTransformation(
                extent={{-56,-106},{-44,-94}})));
        Nodes.Node_in tube_in(redeclare package Medium = Medium_t)
          "Inlet node tube side"  annotation (Placement(
              transformation(extent={{-104,-103},{-84,-83}}), iconTransformation(
                extent={{-106,-36},{-94,-24}})));
        Nodes.Node_out tube_out(redeclare package Medium = Medium_t)
          "Outlet node tube side" annotation (Placement(
              transformation(extent={{-108,-74},{-88,-54}}), iconTransformation(
                extent={{-106,24},{-94,37}})));
        Nodes.Node_in shell_in(redeclare package Medium = Medium_s)
          "Inlet node shell side"  annotation (Placement(
              transformation(extent={{-104,-103},{-84,-83}}), iconTransformation(
                extent={{64,94},{76,106}})));
      equation

          //Mass balance
          shell_in.m_flow = shell_out.m_flow;
          tube_in.m_flow  = tube_out.m_flow;

          //Pressure balance
          if use_dp then
            shell_out.p   = shell_in.p - dp_shell.dp_tot;
            tube_out.p    = tube_in.p - dp_tube.dp_tot;
          else
            shell_in.p    = shell_out.p;
            tube_in.p     = tube_out.p;
          end if;

          //Energy balance
          Q                = pin_s*shell_in.m_flow*(shell_in.h - shell_out.h);
          Q                = pin_t*tube_in.m_flow*(tube_in.h - tube_out.h);

          //Fixed length if off-design mode is active
          if offdesign then
            N_tubes        = N_tubes_od;
          end if;

          //Set boundary conditions at the inlet and outelt
          bundle.p_in        = tube_in.p;
          shell.p_in         = shell_in.p;
          bundle.h[1]        = tube_in.h;
          shell.h[1]         = shell_in.h;
          shell.h[Ncell + 1] = shell_out.h;

          //Things get tricky here. We need to get the index of the shell cell seen by the tube cell
          for j in 1:N_passes loop
               for i in 1+(j-1)*N_baffles_d:2:j*N_baffles_d loop
                  if mod(j,2) <> 0 then
                    if mod(N_baffles_d,2) <> 0 then
                      k[i]   = Ncell - N_passes*(i-(1+(j-1)*N_baffles_d)) - (j-1);
                    else
                      k[i]   = Ncell - N_passes*(i-(1+(j-1)*N_baffles_d)) - (N_passes
                      - j);
                    end if;
                  else
                      k[i]   = (N_passes + 1 - j) + N_passes*(i-(1+(j-1)*N_baffles_d));
                  end if;
               end for;
           end for;

           for j in 1:N_passes loop
               for i in 2+(j-1)*N_baffles_d:2:j*N_baffles_d loop
                  if mod(j,2) <> 0 then
                    if mod(N_baffles_d,2) <> 0 then
                      k[i]   = Ncell - (2*N_passes - 1) + (j-1) - N_passes*(i-(2+(j-1)
                      *N_baffles_d));
                    else
                      k[i]   = Ncell - N_passes - (j-1) - N_passes*(i-(2+(j-1)
                      *N_baffles_d));
                    end if;
                  else
                      k[i]   = N_passes + j + N_passes*(i-(2+(j-1)*N_baffles_d));
                  end if;
               end for;
           end for;

          //A sperate loop for the heat transfer and energy balances
          for i in 1:Ncell loop
            G_wall[i]          = bundle.At*bundle.G_wall/bundle.Dhyd;
            htA_tubes[i]       = bundle.At/(1/hT_tube.ht[i] + 1/ht_t_f[i]);
            htA_shell[k[i]]    = shell.At/(1/hT_shell.ht[k[i]] + 1/ht_s_f[k[i]]);
            kA_tot[i]          = 1/(1/htA_tubes[i] + 1/G_wall[i] + 1/htA_shell[k[i]]);
            bundle.qdot[i]     = shell.qdot[k[i]];
            bundle.qdot[i]     = kA_tot[i]*(pin_s*shell.state[k[i]].T + pin_t*
            bundle.state[i].T);

            //Fluid weight calculation
            bundle.W_fluids[i] = 0.25*pi*Dhyd^2*l/N_baffles_d/N_passes*N_tubes*
            bundle.state[i].d;
            shell.W_fluids[i]  = 0.25*pi*(shell.Dhyd^2 - N_tubes*bundle.Dhyd_o^2)*l/
            N_baffles_d/N_passes*shell.state[i].d;

            //Check second principle of Thermodynamics
            assert(pin_s*shell.state[k[i]].T > pin_t*bundle.state[i].T, "Second principle of Thermodynamics not respected", AssertionLevel.warning);
          end for;

          //Geometric calculations
          d_s                     = bundle.d_b + bundle_clearance(bundle.d_b);
          l_b                     = l/(N_baffles + 1);
          bs_f                    = 1e2*l_b/d_s
          "Baffle spacing as a percentage og the shell diameter 20 - 100 per cent";

          bundle.At     = N_t_p_p*pi*l/N_baffles_d*bundle.Dhyd
          "I see this heat transfer area if I am inside the tubes";
          N_t_p_p       = N_tubes/N_passes "Number of tubes per pass";
          shell.At      = N_t_p_p*pi*l/N_baffles_d*bundle.Dhyd_o
          "I see this heat transfer area if I am outside the tubes";

          //Area, global heat transfer coefficient and corrected DMTL
          A              = N_tubes*pi*l*bundle.Dhyd_o;
          U              = sum(kA_tot)/A;

          //Weight calculation
          bundle.W_dry   = 0.25*pi*(bundle.Dhyd_o^2 - bundle.Dhyd^2)*l*N_tubes
                           *bundle.rho;
          shell.W_dry    = 0.25*pi*(shell.Dhyd_o^2 - shell.Dhyd^2)*l*shell.rho;
          W_dry          = 1.25*bundle.W_dry + 1.3*shell.W_dry;
          W_fluids       = sum(bundle.W_fluids) + sum(shell.W_fluids);
          W_wet          = W_dry + W_fluids;

          //Cost calculation
          PEC            = cost(A);

        annotation (experiment(
            __Dymola_NumberOfIntervals=1,
            Tolerance=1e-006,
            __Dymola_Algorithm="Dassl"),
            __Dymola_experimentSetupOutput,
          Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,100}}),
              graphics={
              Rectangle(
                extent={{-80,50},{90,-50}},
                lineColor={0,0,0},
                fillColor={255,255,255},
                fillPattern=FillPattern.Solid),
              Rectangle(
                extent={{-60,-50},{-40,-70}},
                lineColor={0,0,0},
                fillColor={255,255,255},
                fillPattern=FillPattern.Solid),
              Rectangle(
                extent={{60,70},{80,50}},
                lineColor={0,0,0},
                fillColor={255,255,255},
                fillPattern=FillPattern.Solid),
              Line(
                points={{70,100},{70,56}},
                color={255,0,0},
                smooth=Smooth.None),
              Line(
                points={{2,-2},{-2,2}},
                color={255,0,0},
                smooth=Smooth.None,
                origin={72,58},
                rotation=90),
              Line(
                points={{2,2},{-2,-2}},
                color={255,0,0},
                smooth=Smooth.None,
                origin={68,58},
                rotation=90),
              Line(
                points={{-50,-56},{-50,-100}},
                color={255,0,0},
                smooth=Smooth.None),
              Line(
                points={{2,-2},{-2,2}},
                color={255,0,0},
                smooth=Smooth.None,
                origin={-48,-64},
                rotation=90),
              Line(
                points={{2,2},{-2,-2}},
                color={255,0,0},
                smooth=Smooth.None,
                origin={-52,-64},
                rotation=90),
              Rectangle(extent={{-100,-30},{60,-30}},lineColor={85,170,255}),
              Rectangle(extent={{-100,30},{60,30}},lineColor={85,170,255}),
              Line(
                points={{60,30},{60,-30}},
                color={85,170,255},
                smooth=Smooth.None),
              Line(
                points={{-20,-30},{-24,-26}},
                color={85,170,255},
                smooth=Smooth.None),
              Line(
                points={{-20,-30},{-24,-34}},
                color={85,170,255},
                smooth=Smooth.None),
              Line(
                points={{-20,30},{-16,34}},
                color={85,170,255},
                smooth=Smooth.None),
              Line(
                points={{-20,30},{-16,26}},
                color={85,170,255},
                smooth=Smooth.None)}));
      end shell_and_tube;

      model Flat_plate "Flat plate heat exchanger"
        replaceable package Medium_hot = Media.OneRandomOrganicFluid constrainedby
          Modelica.Media.Interfaces.PartialMedium "Medium model hot cells" annotation(choicesAllMatching = true);
        replaceable package Medium_cold = Media.OneRandomOrganicFluid constrainedby
          Modelica.Media.Interfaces.PartialMedium "Medium model cold cells" annotation(choicesAllMatching = true);
        replaceable Materials.void_material material
          "Material model for the plate"
          annotation(choicesAllMatching = true, Dialog(tab="Geometry"));
        final parameter Integer N_ch = tpg_hot.N_passes*N_ch_p
          "Total number of channels" annotation (Dialog(tab="Geometry"));
        parameter Integer N_ch_p = 3 "Number of channels per pass" annotation (Dialog(tab="Geometry"));
        parameter Integer N_plates = 2*N_ch + 1 "Number of plates" annotation (Dialog(tab="Geometry"));
        parameter Integer N_cell_pc = 3 "Number of cells per channel";
        parameter Boolean use_dp = false
          "Use pressure drop for outlet pressure if true";
        parameter Boolean offdesign = false "Off design mode on if true" annotation (Dialog(tab="Off-design"));
        parameter Modelica.SIunits.Length thick "Plate thickness" annotation (Dialog(tab="Geometry"));
        parameter Modelica.SIunits.Length b "Flow thickness" annotation (Dialog(tab="Geometry"));
        parameter Modelica.SIunits.Length w "Width" annotation (Dialog(tab="Geometry"));
        parameter Real X "Wave number" annotation (Dialog(tab="Geometry"));
        parameter Modelica.SIunits.Angle beta "Plate inclination angle" annotation (Dialog(tab="Geometry"));
        parameter Modelica.SIunits.Length d_pt "Port diameter" annotation (Dialog(tab="Geometry"));
        parameter Modelica.SIunits.CoefficientOfHeatTransfer ht_hot_f1
          "Fouling heat hot side"  annotation(Dialog(tab="Geometry"));
        parameter Modelica.SIunits.CoefficientOfHeatTransfer ht_cold_f1
          "Fouling heat cold side"  annotation(Dialog(tab="Geometry"));
        parameter Modelica.SIunits.Length l_start = 1 "Length" annotation (Dialog(tab="Start"));
        parameter Modelica.SIunits.Length l_od = 1
          "Length fixed for off-design"                                          annotation (Dialog(tab="Off-design"));
        parameter Modelica.SIunits.SpecificEnthalpy h_hot_in_start
          "Inlet specific enthalpy start value hot side" annotation (Dialog(tab="Start"));
        parameter Modelica.SIunits.SpecificEnthalpy h_hot_out_start
          "Outlet specific enthalpy start value hot side" annotation (Dialog(tab="Start"));
        parameter Modelica.SIunits.SpecificEnthalpy h_cold_in_start
          "Inlet specific enthalpy start value cold side" annotation (Dialog(tab="Start"));
        parameter Modelica.SIunits.SpecificEnthalpy h_cold_out_start
          "Outlet specific enthalpy start value cold side" annotation (Dialog(tab="Start"));
        Modelica.SIunits.Length l(start = l_start, fixed = true) "Length";
        Modelica.SIunits.CoefficientOfHeatTransfer U
          "Global heat transfer coefficient";
        Modelica.SIunits.Area A "Heat transfer area";
        Modelica.SIunits.Area At(start = l_start*w/N_cell_pc)
          "Area of one cell";
        Modelica.SIunits.Temp_C  DTML_tilde
          "Logarithmic mean temperature difference corrected";
        Modelica.SIunits.Mass W_dry "Dry weight of the heat exchanger";
        Modelica.SIunits.Mass W_fluids "Weight of the fluids";
        Modelica.SIunits.Mass W_wet "Wet weight of the heat exchanger";
        Real PEC "Purchased equipment cost";
        Modelica.SIunits.HeatFlowRate Q "Heat rate";

        //Plate model
        Objects.plate plate(redeclare package Medium_hot = Medium_hot,
        redeclare package Medium_cold = Medium_cold,
        h_hot_in = node_h_in.h, h_hot_out = node_h_out.h, h_cold_in = node_c_in.h,
        h_cold_out = node_c_out.h, b = b, N_cell_pc = N_cell_pc, X = X,
        N_ch = N_ch, mdot_hot = tpg_hot.mdot_p, pin_hot = tpg_hot.pin,
        mdot_cold = tpg_cold.mdot_p, p_hot_in = node_h_in.p, p_cold_in = node_c_in.p,
        pin_cold = tpg_cold.pin, h_hot_start = tpg_hot.h_start,
        h_cold_start = tpg_cold.h_start);

        //Topologies
        replaceable Miscellanea.topology_PHE.parallel tpg_hot
        constrainedby Miscellanea.topology_PHE.Base_classes.base_top(
        redeclare package Medium = Medium_hot,
          N_cell_pc = N_cell_pc,
          N_ch = N_ch,
          mdot = node_h_in.m_flow,
          h = plate.h_hot,
          h_in = node_h_in.h,
          h_out = node_h_out.h,
          h_in_start = h_hot_in_start,
          h_out_start = h_hot_out_start,
          p_in = node_h_in.p,
          boundary=false,
          p_out = node_h_out.p) annotation(choicesAllMatching = true);

        //Heat transfer model for the hot side
        replaceable Heat_transfer.Plates.single_phase_Martin ht_hot(beta=beta)
        constrainedby Heat_transfer.Plates.Base_classes.base_ht(
        redeclare package Medium = Medium_hot,
          N_cell_pc=N_cell_pc,
          N_ch=N_ch,
          mdot=tpg_hot.mdot_p,
          state=plate.state_hot,
          Aflow=w*b,
          Dhyd=plate.Dhyd) annotation (choicesAllMatching=true);

        //Pressure drop model for the hot side
        replaceable Pressure_drops.Plates.single_phase_Martin dp_hot(beta=beta)
          constrainedby Pressure_drops.Plates.Base_classes.base_dp(
          redeclare package Medium = Medium_hot,
          N_cell_pc=N_cell_pc,
          N_ch=N_ch,
          N_ch_p=N_ch_p,
          N_passes=tpg_hot.N_passes,
          state_p=tpg_hot.state_p,
          mdot=tpg_hot.mdot_p,
          state=plate.state_hot,
          Aflow=w*b,
          Dhyd=plate.Dhyd,
          l=l/N_cell_pc,
          A_pt=0.25*Modelica.Constants.pi*d_pt^2) annotation (choicesAllMatching=true);

        replaceable Miscellanea.topology_PHE.parallel tpg_cold
        constrainedby Miscellanea.topology_PHE.Base_classes.base_top(
        redeclare package Medium = Medium_cold,
          N_cell_pc = N_cell_pc,
          N_ch = N_ch,
          mdot = node_c_in.m_flow,
          h = plate.h_cold,
          h_in = node_c_in.h,
          h_out = node_c_out.h,
          p_in = node_c_in.p,
          h_in_start = h_cold_in_start,
          h_out_start = h_cold_out_start,
          boundary=true,
          p_out = node_c_out.p) annotation(choicesAllMatching = true);

        //Heat transfer model for the cold side
        replaceable Heat_transfer.Plates.single_phase_Martin ht_cold(beta=beta)
        constrainedby Heat_transfer.Plates.Base_classes.base_ht(
        redeclare package Medium = Medium_cold,
          N_cell_pc=N_cell_pc,
          N_ch=N_ch,
          mdot=tpg_cold.mdot_p,
          state=plate.state_cold,
          Aflow=w*b,
          Dhyd=plate.Dhyd) annotation (choicesAllMatching=true);

        //Pressure drop model for the cold side
        replaceable Pressure_drops.Plates.single_phase_Martin dp_cold(beta=beta)
          constrainedby Pressure_drops.Plates.Base_classes.base_dp(
          redeclare package Medium = Medium_cold,
          N_cell_pc=N_cell_pc,
          N_ch=N_ch,
          N_ch_p=N_ch_p,
          N_passes=tpg_cold.N_passes,
          state_p=tpg_cold.state_p,
          mdot=tpg_cold.mdot_p,
          state=plate.state_cold,
          Aflow=w*b,
          Dhyd=plate.Dhyd,
          l=l/N_cell_pc,
          A_pt=0.25*Modelica.Constants.pi*d_pt^2) annotation (choicesAllMatching=true);

        //Defining the model for the cost
        replaceable function cost = Miscellanea.Cost.base_cost annotation(choicesAllMatching = true);

        replaceable Miscellanea.check_velocity check_hot(redeclare package
            Medium =
        Medium_hot, T = plate.state_hot[1, 1].T, umin = min(ht_hot.u), umax = max(ht_hot.u),
        geometry = "tube", op_p = node_h_in.p);

        replaceable Miscellanea.check_velocity check_cold(redeclare package
            Medium =
        Medium_cold, T = plate.state_cold[1, 1].T, umin = min(ht_cold.u), umax = max(ht_cold.u),
        geometry = "tube", op_p = node_c_in.p);

      protected
        parameter Modelica.SIunits.CoefficientOfHeatTransfer ht_hot_f[N_ch, N_cell_pc]=
          fill(ht_hot_f1, N_ch, N_cell_pc)
          "Fouling heat transfer coefficient (array) hot side";
        parameter Modelica.SIunits.CoefficientOfHeatTransfer ht_cold_f[N_ch, N_cell_pc]=
          fill(ht_cold_f1, N_ch, N_cell_pc)
          "Fouling heat transfer coefficient (array) cold side";
         Integer k1[N_ch] "Index for the flat plate topology columns";
         Integer k2[N_cell_pc] "Index for the flat plate topology rows";
         Modelica.SIunits.CoefficientOfHeatTransfer h_wall
          "Thermal conductance of each wall cell";
         Modelica.SIunits.CoefficientOfHeatTransfer G_tot[2, N_ch, N_cell_pc]
          "Global thermal conductance of each cell";
      public
        Nodes.Node_out node_h_out(redeclare package Medium = Medium_hot)
          "Outlet node hot side"
                                annotation (Placement(
              transformation(extent={{-108,-73},{-88,-53}}), iconTransformation(
                extent={{-106,-66},{-94,-54}})));
        Nodes.Node_in node_c_in(redeclare package Medium = Medium_cold)
          "Inlet node cold side"  annotation (Placement(
              transformation(extent={{-104,-103},{-84,-83}}), iconTransformation(
                extent={{-106,-97},{-94,-85}})));
        Nodes.Node_out node_c_out(redeclare package Medium = Medium_cold)
          "Outlet node cold side" annotation (Placement(
              transformation(extent={{-108,-74},{-88,-54}}), iconTransformation(
                extent={{94,43},{106,56}})));
        Nodes.Node_in node_h_in(redeclare package Medium = Medium_hot)
          "Inlet node hot side"  annotation (Placement(
              transformation(extent={{-104,-103},{-84,-83}}), iconTransformation(
                extent={{94,74},{106,86}})));

      equation
            //Mass balance
            node_h_in.m_flow = node_h_out.m_flow;
            node_c_in.m_flow = node_c_out.m_flow;

            //Pressure balance
            if use_dp then
              node_h_out.p      = node_h_in.p - dp_hot.dp_tot;
              node_c_out.p      = node_c_in.p - dp_cold.dp_tot;
            else
              node_h_in.p      = node_h_out.p;
              node_c_in.p      = node_c_out.p;
            end if;

            //Energy balance
            Q                = node_h_in.m_flow*(node_h_in.h - node_h_out.h);
            Q                = node_c_in.m_flow*(node_c_out.h - node_c_in.h);

            //Fixed length if off-design mode is active
            if offdesign then
              l        = l_od;
            end if;

            //Component equations
            for j in 1:N_ch loop
              k1[j] = N_ch + 1 - j;
              for i in 1:N_cell_pc loop

                //Heat transfer coefficient
                G_tot[1, j, i]           = 1/(1/ht_hot.ht[j, i] + 1/ht_hot_f[j, i] +
                1/ht_cold.ht[k1[j], k2[i]] + 1/ht_cold_f[k1[j], k2[i]] + 1/h_wall);

                //The heat rate of the metal wall follows the hot side convention
                plate.qdot_wall[1, j, i] = At*G_tot[1, j, i]*(plate.state_hot[j, i].T
                - plate.state_cold[k1[j], k2[i]].T);
                plate.qdot_hot[j, i]     = plate.qdot_wall[1, j, i] +
                plate.qdot_wall[2, j, i];

                //Weight of the fluids
                plate.W_fluids[j, i]     = At*b*(plate.state_hot[j, i].d +
                plate.state_cold[j, i].d);

                //First row hot side and last row cold side have only one heat flux
                if (j == 1) then
                  k2[i]                         = N_cell_pc + 1 - i;
                  plate.qdot_wall[1, j, i]      = plate.qdot_cold[k1[j], k2[i]];
                  plate.qdot_wall[2, j, i]      = At*G_tot[2, j, i]* (
                  plate.state_hot[j, i].T - plate.state_cold[k1[j] - 1, k2[i]].T);
                  G_tot[2, j, i]                = 1/(1/ht_hot.ht[j, i] +
                  1/ht_hot_f[j, i] + 1/ht_cold.ht[k1[j] - 1, k2[i]] +
                  1/ht_cold_f[k1[j] - 1, k2[i]] + 1/h_wall);
                  //Check second principle of Thermodynamics
                  assert(plate.state_hot[j, i].T > plate.state_cold[k1[j] - 1,
                  k2[i]].T, "Second principle of Thermodynamics not respected",
                  AssertionLevel.warning);
                elseif (j == N_ch) then
                  G_tot[2, j, i]                = 0;
                  plate.qdot_wall[2, j, i]      = 0;
                  plate.qdot_cold[k1[j], k2[i]] = plate.qdot_wall[1, j, i] +
                  plate.qdot_wall[2, j - 1, i];
                else
                  plate.qdot_wall[2, j, i]      =  At*G_tot[2, j, i]*(
                  plate.state_hot[j, i].T - plate.state_cold[k1[j] - 1, k2[i]].T);
                  plate.qdot_cold[k1[j], k2[i]] = plate.qdot_wall[1, j, i] +
                  plate.qdot_wall[2, j - 1, i];
                  G_tot[2, j, i]                = 1/(1/ht_hot.ht[j, i] +
                  1/ht_hot_f[j, i] + 1/ht_cold.ht[k1[j], k2[i]] +
                  1/ht_cold_f[k1[j], k2[i]] + 1/h_wall);
                  //Check second principle of Thermodynamics
                  assert(plate.state_hot[j, i].T > plate.state_cold[k1[j] - 1,
                  k2[i]].T, "Second principle of Thermodynamics not respected",
                  AssertionLevel.warning);
                end if;

                //Check second principle of Thermodynamics
                 assert(plate.state_hot[j, i].T > plate.state_cold[k1[j], k2[i]].T,
                 "Second principle of Thermodynamics not respected",
                 AssertionLevel.warning);
              end for;
            end for;

            //Area and thermal conductance of the metal wall
            At             = l*w/N_cell_pc;
            h_wall         = material.lambda/thick;

            //Weight calculation
            plate.W_dry    = l*w*thick*material.rho;
            W_dry          = N_plates*plate.W_dry;
            W_fluids       = sum(plate.W_fluids);
            W_wet          = W_dry + W_fluids;

            //Area, global heat transfer coefficient and corrected DMTL
            A              = 2*l*w*N_ch;
            U              = At*sum(G_tot)/A;
            Q              = DTML_tilde*U*A;

            //Cost calculation
            PEC            = cost(A);

        annotation (experiment(Tolerance=1e-006), __Dymola_experimentSetupOutput,
          Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},
                  {100,100}}),
              graphics={
              Rectangle(
                extent={{-80,93},{80,-99}},
                lineColor={255,255,255},
                fillColor={255,255,255},
                fillPattern=FillPattern.Solid),
              Line(
                points={{-60,81},{100,81}},
                color={255,0,0},
                smooth=Smooth.None),
              Line(
                points={{-60,81},{-60,-59}},
                color={255,0,0},
                smooth=Smooth.None),
              Line(
                points={{20,81},{20,-59}},
                color={255,0,0},
                smooth=Smooth.None),
              Line(
                points={{-100,-59},{20,-59}},
                color={255,0,0},
                smooth=Smooth.None),
              Line(
                points={{-20,50},{100,50}},
                color={0,128,255},
                smooth=Smooth.None),
              Line(
                points={{-100,-90},{60,-90}},
                color={0,128,255},
                smooth=Smooth.None),
              Line(
                points={{-66,7},{-60,1},{-54,7}},
                color={255,0,0},
                smooth=Smooth.None),
              Line(
                points={{14,7},{20,1},{26,7}},
                color={255,0,0},
                smooth=Smooth.None),
              Line(
                points={{-20,50},{-20,-90}},
                color={0,128,255},
                smooth=Smooth.None),
              Line(
                points={{6,3},{0,-3},{-6,3}},
                color={0,128,255},
                smooth=Smooth.None,
                origin={-20,-2},
                rotation=180),
              Line(
                points={{60,50},{60,-90}},
                color={0,128,255},
                smooth=Smooth.None),
              Line(
                points={{6,3},{0,-3},{-6,3}},
                color={0,128,255},
                smooth=Smooth.None,
                origin={60,-2},
                rotation=180),
              Rectangle(
                extent={{-42,47},{-38,-53}},
                lineColor={0,0,0},
                fillColor={255,255,255},
                fillPattern=FillPattern.VerticalCylinder),
              Rectangle(
                extent={{-2,47},{2,-53}},
                lineColor={0,0,0},
                fillColor={255,255,255},
                fillPattern=FillPattern.VerticalCylinder),
              Rectangle(
                extent={{38,47},{42,-53}},
                lineColor={0,0,0},
                fillColor={255,255,255},
                fillPattern=FillPattern.VerticalCylinder)}));
      end Flat_plate;
      annotation (Icon(graphics={
            Ellipse(extent={{-80,80},{80,-80}}, lineColor={0,0,0},
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid),
            Line(
              points={{-40,40},{40,-40}},
              color={0,0,0},
              smooth=Smooth.None),
            Line(
              points={{-80,0},{-60,0},{-40,40}},
              color={0,0,0},
              smooth=Smooth.None),
            Line(
              points={{80,0},{60,0},{40,-40}},
              color={0,0,0},
              smooth=Smooth.None)}));
    end HEX;

    package Nodes "Node package"

      connector hp_mdot "A node with a thermodynamic state"
        replaceable package Medium = Test.Media.OneRandomOrganicFluid constrainedby
          Modelica.Media.Interfaces.PartialMedium "Medium model hot cells" annotation(choicesAllMatching = true);
        Medium.MassFlowRate m_flow;
        Modelica.SIunits.SpecificEnthalpy h "Specific enthalpy";
        Modelica.SIunits.AbsolutePressure p "Pressure";
      end hp_mdot;

      connector Node_in "Fluid connector with filled icon"
        extends hp_mdot;
        annotation(Diagram(coordinateSystem(preserveAspectRatio = true, extent = {{-100, -100}, {100, 100}}), graphics={  Ellipse(extent = {{-100, 100}, {100, -100}}, lineColor=
                    {0,127,255},                                                                                                    fillColor=
                    {255,255,255},
                  fillPattern=FillPattern.Solid),                                                                                                    Ellipse(extent = {{-100, 100}, {100, -100}}, lineColor=
                    {0,0,0}),                                                                                                    Text(extent={{
                    -100,194},{100,100}},                                                                                                    textString = "%name", lineColor=
                    {0,0,0})}),                                                                                                    Icon(coordinateSystem(preserveAspectRatio = true, extent = {{-100, -100}, {100, 100}}), graphics={  Ellipse(extent = {{-100, 100}, {100, -100}},
                pattern=LinePattern.None,
                lineColor={0,0,0}),                                                                                                    Ellipse(extent = {{-100, 100}, {100, -100}}, lineColor=
                    {0,0,0},
                fillColor={255,255,255},
                fillPattern=FillPattern.Solid)}),                                                                                                    Documentation(info = "<html>Modelica.Media.Examples.Tests.Components.FluidPort_a
                                                                                                </html>"));
      end Node_in;

      connector Node_out "Fluid connector with filled icon"
        extends hp_mdot;
        annotation(Diagram(coordinateSystem(preserveAspectRatio = true, extent = {{-100, -100}, {100, 100}}), graphics={  Ellipse(extent = {{-100, 100}, {100, -100}}, lineColor=
                    {0,127,255},                                                                                                    fillColor=
                    {0,0,0},
                  fillPattern=FillPattern.Solid),                                                                                                    Ellipse(extent = {{-100, 100}, {100, -100}}, lineColor=
                    {0,0,0}),                                                                                                    Text(extent={{
                    -94,192},{106,98}},                                                                                                    textString = "%name", lineColor=
                    {0,0,0})}),                                                                                                    Icon(coordinateSystem(preserveAspectRatio = true, extent = {{-100, -100}, {100, 100}}), graphics={  Ellipse(extent = {{-100, 100}, {100, -100}},
                pattern=LinePattern.None,
                lineColor={0,0,0}),                                                                                                    Ellipse(extent = {{-100, 100}, {100, -100}}, lineColor=
                    {0,0,0},
                fillColor={0,0,0},
                fillPattern=FillPattern.Solid)}),                                                                                                    Documentation(info = "<html>Modelica.Media.Examples.Tests.Components.FluidPort_a
                                                                                                </html>"));
      end Node_out;

      annotation (Icon(graphics={
            Ellipse(lineColor={255,0,0},
              pattern=LinePattern.None,
              extent={{-100,-100},{100,100}},
              fillColor={255,0,0},
              fillPattern=FillPattern.Sphere)}));
    end Nodes;

    package Flags "Components to impose or initialize variables"

      model START
        "Flag to initialize a condition related to a process variable"
        replaceable package Medium = Media.OneRandomOrganicFluid constrainedby
          Modelica.Media.Interfaces.PartialMedium "Medium model hot cells" annotation(choicesAllMatching = true);
        parameter Modelica.SIunits.AbsolutePressure p = 1e5 "Pressure" annotation (Dialog(tab="Start"));
        parameter Modelica.SIunits.SpecificEnthalpy h = 1e5 "Enthalpy" annotation (Dialog(tab="Start"));
        parameter Modelica.SIunits.MassFlowRate m_flow = 1 "Mass flow" annotation (Dialog(tab="Start"));
        CycleTempo.Components.Nodes.Node_in node(
          redeclare package Medium = Medium,
          p(start=p),
          h(start=h),
          m_flow(start=m_flow)) annotation (Placement(transformation(extent={{90,-10},
                  {110,10}}), iconTransformation(extent={{90,-10},{110,10}})));

        annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
                  -100},{100,100}}), graphics), Icon(coordinateSystem(
                preserveAspectRatio=false, extent={{-100,-100},{100,100}}), graphics={
              Rectangle(
              extent={{-40,-50},{60,50}},
              lineColor={0,0,0},
              fillColor={0,255,0},
              fillPattern=FillPattern.Solid), Line(
              points={{101,0},{60,0}},
              color={0,0,0},
              smooth=Smooth.None),
              Text(
                extent={{0,20},{22,-20}},
                lineColor={0,0,0},
                textString="S")}));
      end START;

      model ADDCO "Flag to impose and visualize the process variables"
        replaceable package Medium = Media.OneRandomOrganicFluid constrainedby
          Modelica.Media.Interfaces.PartialMedium "Medium model hot cells" annotation(choicesAllMatching = true);
        parameter Boolean visible = true "if true display mode on";
        parameter Modelica.SIunits.Temperature T = 298.15 "Temperature" annotation (Dialog(tab="Addco"));
        parameter Modelica.SIunits.AbsolutePressure p = 1e5 "Pressure" annotation (Dialog(tab="Addco"));
        parameter Modelica.SIunits.SpecificEnthalpy h = 1e5 "Enthalpy" annotation (Dialog(tab="Addco"));
        parameter Modelica.SIunits.MassFlowRate m_flow = 1 "Mass flow" annotation (Dialog(tab="Addco"));
        parameter Boolean use_T = false "True if temperature is given"  annotation (Dialog(tab="Addco"));
        parameter Boolean use_p = false "True if pressure is given"  annotation (Dialog(tab="Addco"));
        parameter Boolean use_h = false "True if enthalpy is given"  annotation (Dialog(tab="Addco"));
        parameter Boolean use_m_flow = false "True if mass flow is given"  annotation (Dialog(tab="Addco"));
        Medium.ThermodynamicState state "Thermodynamic state";
        CycleTempo.Components.Nodes.Node_in node(redeclare package Medium = Medium)
          annotation (Placement(transformation(extent={{91,-10},{111,10}}),
              iconTransformation(extent={{91,-10},{111,10}})));
      equation
        //Component equations
        state = Medium.setState_ph(node.p, node.h);

        //Boundary equations
        if use_T then
          node.h =  Medium.specificEnthalpy(Medium.setState_pT(node.p, T));
        end if;
        if use_p then
          p      = node.p;
        end if;
        if use_h then
          h      = node.h;
        end if;
        if use_m_flow then
          m_flow = node.m_flow;
        end if;
        annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
                  -100},{100,100}}), graphics), Icon(coordinateSystem(
                preserveAspectRatio=false, extent={{-100,-100},{100,100}}), graphics={
              Rectangle(extent={{-280,122},{-40,-122}}, lineColor={0,0,0},
                fillColor={255,255,255},
                fillPattern=FillPattern.Solid),
              Text(visible=visible,
              extent={{-269,-150},{-179,-30}},
              textString=DynamicSelect("0.0", String(Modelica.SIunits.Conversions.to_bar(node.p), significantDigits=2, format=  ".1f"))),
              Text(visible=visible,
              extent={{-268,27},{-178,147}},
              textString=DynamicSelect("0.0", String(node.m_flow, significantDigits=2, format=  ".2f"))),
              Text(visible=visible,
              extent={{-268,-68},{-179,62}},
              textString=DynamicSelect("0.0", String(Modelica.SIunits.Conversions.to_degC(state.T), significantDigits=2, format=  ".1f"))),
              Text(visible=visible,
              extent={{-140,-68},{-51,62}},
              textString=DynamicSelect("°C", String("°C"))),
              Text(visible=visible,
              extent={{-140,-149},{-51,-19}},
              textString=DynamicSelect("bar", String("bar"))),
              Text(visible=visible,
              extent={{-140,16},{-51,146}},
              textString=DynamicSelect("kg/s", String("kg/s"))),
              Rectangle(
              extent={{-40,-50},{60,50}},
              lineColor={0,0,0},
              fillColor={0,128,255},
              fillPattern=FillPattern.Solid), Line(
              points={{101,0},{60,0}},
              color={0,0,0},
              smooth=Smooth.None),
              Text(
                extent={{0,20},{22,-20}},
                lineColor={0,0,0},
                textString="A"),                                                                  Text(lineColor=
                    {0,0,0},                                                                                                    extent={{
                    -112,27.5},{112,-27.5}},                                                                                                textString = "%name",
                origin={-162,150.5},
                rotation=0)}));
      end ADDCO;

    end Flags;
  end Components;

  package Materials "Package containing the properties of different materials"
    extends Modelica.Icons.MaterialPropertiesPackage;
    class S_AISI_1010 "S_AISI_1010"
        extends Design.Materials.void_material(
        materialName="S_AISI_1010",
        lambda=65.2,
        rho=7850);
    end S_AISI_1010;

    class SS_AISI_410 "SS_AISI_410"
        extends Design.Materials.void_material(
        materialName="SS_AISI_410",
        lambda=24.9,
        rho=7740);
    end SS_AISI_410;

    class SS_AISI_304 "SS_AISI_304"
        extends Design.Materials.void_material(
        materialName="SS_AISI_304",
        lambda=18.3866,
        rho=8030);
    end SS_AISI_304;

    class S_AISI_1040 "S_AISI_1040"
        extends Design.Materials.void_material(
        materialName="S_AISI_1040",
        lambda=50.7,
        rho=7845);
    end S_AISI_1040;

    class void_material "Void material"
        extends Modelica.Icons.MaterialProperty;
        parameter String materialName = "Random";
        parameter Modelica.SIunits.ThermalConductivity lambda = 50;
        parameter Modelica.SIunits.Density rho = 8000;
    end void_material;

    class T_R50250 "Titanium alloy R50250"
        extends Design.Materials.void_material(
        materialName="T_R50250",
        lambda=20.8,
        rho=4510);
    end T_R50250;
  end Materials;

  package Miscellanea "Package containing useful functions"
    function log_mean_delta_T "logarithmic mean temperature difference"
      input Real T1;
      input Real T2;
      input Real t1;
      input Real t2;
      output Real LMTD;

    algorithm
      LMTD :=((T1 - t2) - (T2 - t1))/log((T1 - t2)/(T2 - t1));
    end log_mean_delta_T;

    function bundle_diameter "Function to calculate the bundle diameter"

      parameter Integer K_ind[8] = {1, 2, 0, 3, 0, 4, 0, 5};
      parameter Real K[2,5] = [0.319, 0.249, 0.175, 0.0743, 0.0365;
                          0.215, 0.156, 0.158, 0.0402, 0.0331];
      parameter Real n[2,5] = [2.142, 2.207, 2.285, 2.4990, 2.6750;
                          2.207, 2.291, 2.263, 2.6170, 2.6430];
      input Real N_tubes;
      input Integer N_passes;
      input Real d_eq;
      input Integer layout;
      output Real d_b;

    algorithm
      d_b :=d_eq*(N_tubes/K[layout, K_ind[N_passes]])^(1/n[layout, K_ind[N_passes]]);
    end bundle_diameter;

    class check_velocity
      "This class is used to check the velocities inside heat exchangers"
      replaceable package Medium = Design.Media.OneRandomOrganicFluid
        "Medium model";
      input Modelica.SIunits.Velocity umax "Maximum velocity";
      input Modelica.SIunits.Velocity umin "Minimum velocity";
      parameter String geometry = "tube"
        "Type of geometry. E.g. shell, tubes, etc...";
      input Modelica.SIunits.AbsolutePressure op_p
        "Operating pressure. If the fluid is a vapour then the limit depends on the pressure";
      input Modelica.SIunits.Temperature T_sat = Medium.saturationTemperature(op_p)
        "Saturation temperature";
      input Modelica.SIunits.Temperature T "Fluid temperature";

    equation
      if (T < T_sat) then
        if (geometry == "tube") then
          assert(umax < 4*2, "Velocity twice higher than 4 m/s", AssertionLevel.warning);
          assert(umin > 1*0.5,"Velocity twice lower than 0.8 m/s", AssertionLevel.warning);
        elseif (geometry == "shell") then
          assert(umax < 1.2*2, "Velocity twice higher than 1.2 m/s", AssertionLevel.warning);
          assert(umin > 0.3*0.5,"Velocity twice lower than 0.3 m/s", AssertionLevel.warning);
        end if;
      elseif (T > T_sat) then
        if (op_p > 1e6) then
          assert(umax < 70*2, "Velocity twice higher than 70 m/s", AssertionLevel.warning);
          assert(umin > 50*0.5,"Velocity twice lower than 50 m/s", AssertionLevel.warning);
        elseif (op_p < 1e4) then
          assert(umax < 10*2, "Velocity twice higher than 10 m/s", AssertionLevel.warning);
          assert(umin > 5*0.5,"Velocity twice lower than 5 m/s", AssertionLevel.warning);
        else
          assert(umax < 30*2, "Velocity twice higher than 30 m/s", AssertionLevel.warning);
          assert(umin > 10*0.5,"Velocity twice lower than 10 m/s", AssertionLevel.warning);
        end if;
      end if;
    end check_velocity;

    package Cost "Functions to evaluate the cost of heat exchangers "
      function CS_Hall "Carbon steel heat exchanger Hall et al."
         extends base_cost(a = {30800, 750, 0.81});
      end CS_Hall;

      function base_cost
        "Base class for the cost calculation of heat exchangers"
        input Modelica.SIunits.Area A "heat transfer area";
        output Real PEC "Purchased equipment cost";
      protected
        parameter Real a[3] = {30800, 1644, 0.81};
      algorithm
          PEC := a[1] + a[2]*A^a[3];
      end base_cost;

      function SS_Hall "Stainless steel heat exchanger Hall et al."
         extends base_cost(a = {30800, 1644, 0.81});
      end SS_Hall;

      function FP_Rafferty
        "Stainless steel flat plate heat exchanger Rafferty et al."
         extends base_cost(a = {1378.1781, 10.8962, 1});
      end FP_Rafferty;
      annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
                -100},{100,100}}), graphics), Icon(graphics={
                                 Text(
              extent={{-127,-13},{128,-208}},
              lineColor={0,0,0},
              textString="$
")}));
    end Cost;

    package Shell_clearance
      "Functions to evaluate the shell clearance based on the diameter of the tube bundle "
      function Fixed_Utube "Fixed_U_tube shape clearance shell-bundle diameter"
         extends base_clearance(a = {10, 10, 0.2});
      end Fixed_Utube;

      function OPH "Outside packed head clearance shell-bundle diameter"
         extends base_clearance(a = {38, 0, 0.2});
      end OPH;

      function SRFH "Split ring floating head clearance shell-bundle diameter"
         extends base_clearance(a = {50, 28, 0.2});
      end SRFH;

      function PTFH
        "Pull through floating head clearance shell-bundle diameter"
         extends base_clearance(a = {88, 11, 0.2});
      // algorithm
      //   clearance :=1e-3*(88 + 11*(d_b - 0.2));
      end PTFH;

      function base_clearance "Base class for clearance shell-bundle diameter"
        input Modelica.SIunits.Length d_b "bundle diameter";
        output Modelica.SIunits.Length clearance
          "Shell to bundle diameter clearance";
      protected
        parameter Real a[3] = {10, 0, 0};
      algorithm
          clearance :=1e-3*(a[1] + a[2]*(d_b - a[3]));
      end base_clearance;
      annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{
                -100,-100},{100,100}}), graphics={
            Line(
              points={{-40,0},{40,0}},
              color={0,0,0},
              smooth=Smooth.None),
            Line(
              points={{-40,-20},{-40,20}},
              color={0,0,0},
              smooth=Smooth.None),
            Line(
              points={{40,-20},{40,20}},
              color={0,0,0},
              smooth=Smooth.None)}));
    end Shell_clearance;

    package numbers "Package containing adimensional numbers"
      function Prandtl "Prandtl number"
        input Real cp;
        input Real eta;
        input Real lambda;
        output Real Pr;

      algorithm
        Pr := cp*eta/lambda;
      end Prandtl;

      function Nusselt "Nusselt number"
        input Real ht;
        input Real d_eq;
        input Real lambda;
        output Real Nu;

      algorithm
        Nu := ht*d_eq/lambda;
      end Nusselt;

      function Reynolds "Reynolds number"
        input Real u;
        input Real d;
        input Real mu;
        input Real d_eq;
        output Real Re;

      algorithm
        Re := u*d*d_eq/mu;
      end Reynolds;
    end numbers;

    package topology_PHE "Topologies for the flat plate heat exchanger"
      class series "Series"
        extends Base_classes.base_top;
        final parameter Integer N_passes = 1 "Number of passes";
        input Modelica.SIunits.MassFlowRate mdot_p = mdot
          "Mass flow rate plate channel";
         parameter Modelica.SIunits.SpecificEnthalpy h_start[N_ch, N_cell_pc + 1]=
         fill(linspace(h_in_start, h_out_start, N_cell_pc + 1), N_ch)
          "Start value for the hot stream enthalpy matrix";
        Real pin[N_ch, N_cell_pc] "Pin for the heat flow";
        Medium.ThermodynamicState state_p[1]
          "Thermodynamic states at the inlet port";

      equation
        for j in 1:N_ch loop
          for i in 1:N_cell_pc loop
            if (mod(j,2) <> 0) then
              pin[j, i]  = N_ch;
            else
              pin[j, i]  = -N_ch;
            end if;
            //Boundary conditions
            if (i == 1) then
              if (j == 1) then
                h_in                 = h[1, 1];
              else
                if (mod(j,2) <> 0) then
                  h[j, 1]              = h[j - 1, 1];
                else
                  h[j, N_cell_pc + 1]  = h[j - 1, N_cell_pc + 1];
                end if;
              end if;
            end if;
          end for;
        end for;

        if boundary then
          h_out      = h[N_ch, N_cell_pc + 1];
        end if;

        state_p[1] = Medium.setState_ph(p_in, h_in);

      end series;

      class parallel "Parallel"
        extends Base_classes.base_top;
        final parameter Integer N_passes = 1 "Number of passes";
        input Modelica.SIunits.MassFlowRate mdot_p = mdot/N_ch
          "Mass flow rate plate channel";
         parameter Modelica.SIunits.SpecificEnthalpy h_start[N_ch, N_cell_pc + 1]=
         fill(linspace(h_in_start, h_out_start, N_cell_pc + 1), N_ch)
          "Start value for the hot stream enthalpy matrix";
        Real pin[N_ch, N_cell_pc] "Pin for the heat flow";
        Medium.ThermodynamicState state_p[1]
          "Thermodynamic states at the inlet port";
      equation
        for j in 1:N_ch loop
          for i in 1:N_cell_pc loop
            pin[j, i]  = 1;
            //Boundary conditions
            if (i == 1) then
              h[j, i]       = h_in;
            end if;
          end for;
        end for;

        if boundary then
          h_out      = sum(h[:, N_cell_pc + 1])/N_ch;
        end if;

        state_p[1] = Medium.setState_ph(p_in, h_in);

      end parallel;

      class two_pass_two_pass "Two pass two pass configuration"
        extends Base_classes.base_top;
        input Modelica.SIunits.MassFlowRate mdot_p = mdot/N_ch_p
          "Mass flow rate plate channel";
        final parameter Integer N_passes = 2 "Number of passes";
        parameter Integer N_ch_p( start = 3) "Number of channels per pass";
        parameter Integer stype "Stream type";
        parameter Modelica.SIunits.SpecificEnthalpy
        h_start[N_ch, N_cell_pc + 1] = (if stype == 1 then
        [fill(linspace(h_out_start, 0.5*(h_out_start + h_in_start), N_cell_pc + 1), N_ch_p);
        fill(linspace(h_in_start, 0.5*(h_out_start + h_in_start), N_cell_pc + 1), N_ch_p)] else
        [fill(linspace(0.5*(h_out_start + h_in_start), h_out_start, N_cell_pc + 1), N_ch_p);
        fill(linspace(0.5*(h_out_start + h_in_start), h_in_start, N_cell_pc + 1), N_ch_p)])
          "Start value for the hot stream enthalpy matrix";
        Real pin[N_ch, N_cell_pc] "Pin for the heat flow";
        Medium.ThermodynamicState state_p[2]
          "Thermodynamic states at the ports";
      equation
        if (stype == 1) then
          for j in 1:N_ch loop
            for i in 1:N_cell_pc loop
              if (j <= N_ch_p) then
                pin[j, i]             = -1;
                if (i == N_cell_pc) then
                  h[j, i + 1] = sum(h[N_ch_p + 1:N_ch, N_cell_pc + 1])/N_ch_p;
                end if;
              else
                pin[j, i]             = 1;
                if (i == 1) then
                  h[j, i]             = h_in;
                end if;
              end if;
            end for;
          end for;
          if boundary then
            h_out = sum(h[1:N_ch_p, 1])/N_ch_p;
          end if;
          state_p[1] = Medium.setState_ph(p_in, sum(h[N_ch_p + 1:N_ch, N_cell_pc + 1]
          /N_ch_p));
          state_p[2] = Medium.setState_ph(p_in, h_in);
        else
          for j in 1:N_ch loop
            for i in 1:N_cell_pc loop
              if (j <= N_ch_p) then
                pin[j, i]             = 1;
                if (i == 1) then
                  h[j, i] = sum(h[N_ch_p + 1:N_ch, 1])/N_ch_p;
                end if;
              else
                pin[j, i]             = -1;
                if (i == N_cell_pc) then
                  h[j, i + 1]         = h_in;
                end if;
              end if;
            end for;
          end for;
          if boundary then
            h_out = sum(h[1:N_ch_p, N_cell_pc + 1])/N_ch_p;
          end if;
          state_p[1] = Medium.setState_ph(p_in, sum(h[N_ch_p + 1:N_ch, 1])/N_ch_p);
          state_p[2] = Medium.setState_ph(p_in, h_in);
        end if;

      end two_pass_two_pass;

      package Base_classes "Base classes for the topology of the PHE"
        class base_top "Base class topology PHE"
          replaceable package Medium = Media.OneRandomOrganicFluid
            "Medium model";
          parameter Integer N_ch(start = 6) "Total number of channels";
          parameter Integer N_cell_pc = 3 "Number of cells per channel";
          parameter Boolean boundary "Boundary condition for the enthalpy";
          parameter Modelica.SIunits.SpecificEnthalpy h_in_start
            "Inlet specific enthalpy start value";
          parameter Modelica.SIunits.SpecificEnthalpy h_out_start
            "Outlet specific enthalpy start value";
          input Modelica.SIunits.SpecificEnthalpy h_in
            "Inlet specific enthalpy";
          input Modelica.SIunits.SpecificEnthalpy h_out
            "Outlet specific enthalpy";
          input Modelica.SIunits.AbsolutePressure p_in "Inlet pressure";
          input Modelica.SIunits.AbsolutePressure p_out "Outlet pressure";
          input Modelica.SIunits.MassFlowRate mdot "Mass flow rate";
          input Modelica.SIunits.SpecificEnthalpy h[N_ch, N_cell_pc + 1]
            "Stream enthalpy matrix";
        end base_top;
      end Base_classes;
    end topology_PHE;
  annotation (Icon(graphics={
          Rectangle(
            lineColor={200,200,200},
            fillColor={248,248,248},
            fillPattern=FillPattern.HorizontalCylinder,
            extent={{-100,-100},{100,100}},
            radius=25.0),
          Rectangle(
            lineColor={128,128,128},
            fillPattern=FillPattern.None,
            extent={{-100,-100},{100,100}},
            radius=25.0),
          Ellipse(
            origin={10,10},
            fillColor={76,76,76},
            pattern=LinePattern.None,
            fillPattern=FillPattern.Solid,
            extent={{-80.0,-80.0},{-20.0,-20.0}}),
          Ellipse(
            origin={10,10},
            pattern=LinePattern.None,
            fillPattern=FillPattern.Solid,
            extent={{0.0,-80.0},{60.0,-20.0}}),
          Ellipse(
            origin={10,10},
            fillColor={128,128,128},
            pattern=LinePattern.None,
            fillPattern=FillPattern.Solid,
            extent={{0.0,0.0},{60.0,60.0}}),
          Ellipse(
            origin={10,10},
            lineColor={128,128,128},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid,
            extent={{-80.0,0.0},{-20.0,60.0}})}));
  end Miscellanea;

  package Objects "Package containing all the objects of the VIP"

    class tube "I am a tube and I contain all the my relevant informations"
        replaceable package Medium = Media.OneRandomOrganicFluid constrainedby
        Modelica.Media.Interfaces.PartialMedium "Medium model" annotation(choicesAllMatching = true);
        parameter Integer Ncell(start=3) "Number of cell elements";
        parameter Modelica.SIunits.Length thick "Thickness";
        parameter Modelica.SIunits.ThermalConductivity lambda
        "Thermal conductivity of the wall";
        parameter Modelica.SIunits.Density rho "Density of the wall";
        parameter Real pin "This pin identifies if the fluid is hot or cold";
        parameter Modelica.SIunits.SpecificEnthalpy h_in_start
        "Inlet specific enthalpy start value hot side";
        parameter Modelica.SIunits.SpecificEnthalpy h_out_start
        "Outlet specific enthalpy start value hot side";
        input Modelica.SIunits.SpecificEnthalpy h_in "Inlet specific enthalpy";
        input Modelica.SIunits.SpecificEnthalpy h_out
        "Outlet specific enthalpy";
        input Modelica.SIunits.MassFlowRate mdot "Mass flow rate";
        input Modelica.SIunits.Length Dhyd "Hydraulic diameter";
        input Modelica.SIunits.Area Aflow "Flow area";
        Modelica.SIunits.Area At "Heat transfer area";
        Modelica.SIunits.Length Dhyd_o "Outer hydraulic diameter";
        Modelica.SIunits.ThermalConductivity G_wall
        "Thermal conductivity of the wall";
        Modelica.SIunits.AbsolutePressure p_in "Inlet pressure";
        Medium.ThermodynamicState state[Ncell]
        "Thermodynamic states of the cells";
        Modelica.SIunits.SpecificEnthalpy h[Ncell + 1](
         start=linspace(h_in_start, h_out_start, Ncell + 1))
        "Specific enthalpies";
        Modelica.SIunits.HeatFlowRate qdot[Ncell] "Heat rate";
        Modelica.SIunits.Mass W_dry "Dry weight of one element";
        Modelica.SIunits.Mass W_fluids[Ncell] "Fluids weight of one element";

    equation
        //Trivial calculations
        Dhyd_o = Dhyd + 2*thick;
        G_wall = 2*lambda/log(Dhyd_o/Dhyd);

        for i in 1:Ncell loop
          qdot[i]           = pin*mdot*(h[i] - h[i + 1]);
          state[i]          = Medium.setState_ph(p_in, 0.5*(h[i] + h[i + 1]));
        end for;
      annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},
                {100,100}}), graphics={
            Rectangle(extent={{-80,-40},{60,-40}}, lineColor={85,170,255}),
            Rectangle(extent={{-80,20},{60,20}}, lineColor={85,170,255}),
            Line(
              points={{60,20},{60,-40}},
              color={85,170,255},
              smooth=Smooth.None),
            Line(
              points={{-20,-40},{-24,-36}},
              color={85,170,255},
              smooth=Smooth.None),
            Line(
              points={{-20,-40},{-24,-44}},
              color={85,170,255},
              smooth=Smooth.None),
            Line(
              points={{-20,20},{-16,24}},
              color={85,170,255},
              smooth=Smooth.None),
            Line(
              points={{-20,20},{-16,16}},
              color={85,170,255},
              smooth=Smooth.None)}));
    end tube;

    class tube_bundle
      "I am a tube bundle and I contain all the my relevant information"
        extends Objects.tube;
        parameter Integer N_passes "Number of tube passes";
        parameter Integer layout "Tube layout, 1 = triangular, 2 = squared";
        parameter Real pitch_f "Tube pitch";
        input Real N_tubes "Number of tubes in the bundle";
        Modelica.SIunits.Length d_b "Bundle diameter";
    equation
        d_b  =Miscellanea.bundle_diameter(
            N_tubes,
            N_passes,
            Dhyd_o,
            layout) "calculate the tube bundle";
    end tube_bundle;

    class plate "I am a plate and I contain all the my relevant informations"
        replaceable package Medium_hot = Media.OneRandomOrganicFluid constrainedby
        Modelica.Media.Interfaces.PartialMedium "Medium model hot cells";
        replaceable package Medium_cold = Media.OneRandomOrganicFluid constrainedby
        Modelica.Media.Interfaces.PartialMedium "Medium model cold cells";
        parameter Integer N_ch(start = 6) "Total number of channels";
        parameter Integer N_cell_pc(start = 3) "Number of cells per channels";
        parameter Modelica.SIunits.Length b "Plate thickness";
        parameter Real X "Wavenumber";
        input Modelica.SIunits.SpecificEnthalpy h_hot_in
        "Inlet specific enthalpy hot side";
        input Modelica.SIunits.SpecificEnthalpy h_hot_out
        "Outlet specific enthalpy hot side";
        input Modelica.SIunits.SpecificEnthalpy h_cold_in
        "Inlet specific enthalpy cold side";
        input Modelica.SIunits.SpecificEnthalpy h_cold_out
        "Outlet specific enthalpy cold side";
        final parameter Modelica.SIunits.Length Dhyd = 2*b/phi
        "Hydraulic diameter";
        input Modelica.SIunits.AbsolutePressure p_hot_in
        "Inlet pressure hot cells";
        input Modelica.SIunits.AbsolutePressure p_cold_in
        "Inlet pressure cold cells";
        parameter Modelica.SIunits.SpecificEnthalpy h_hot_start[N_ch, N_cell_pc + 1]
        "Hot stream temperature matrix";
        parameter Modelica.SIunits.SpecificEnthalpy h_cold_start[N_ch, N_cell_pc + 1]
        "Cold stream temperature matrix";
        input Modelica.SIunits.MassFlowRate mdot_hot "Mass flow rate hot cells";
        input Modelica.SIunits.MassFlowRate mdot_cold
        "Mass flow rate cold cells";
        input Real pin_hot[N_ch, N_cell_pc] "Pin hot cells";
        input Real pin_cold[N_ch, N_cell_pc] "Pin cold cells";
        Medium_hot.ThermodynamicState state_hot[N_ch, N_cell_pc]
        "Thermodynamic states hot cells";
        Medium_cold.ThermodynamicState state_cold[N_ch, N_cell_pc]
        "Thermodynamic states cold cells";
        Modelica.SIunits.HeatFlowRate qdot_hot[N_ch, N_cell_pc]
        "Heat rate hot cells";
        Modelica.SIunits.HeatFlowRate qdot_cold[N_ch, N_cell_pc]
        "Heat rate cold cells";
        Modelica.SIunits.HeatFlowRate qdot_wall[2, N_ch, N_cell_pc]
        "Heat rate metal wall";
        Modelica.SIunits.SpecificEnthalpy h_hot[N_ch, N_cell_pc + 1](
         start = h_hot_start) "Hot stream enthalpy matrix";
        Modelica.SIunits.SpecificEnthalpy h_cold[N_ch, N_cell_pc + 1](
         start = h_cold_start) "Cold stream enthalpy matrix";
        Modelica.SIunits.Mass W_dry "Dry weight of the plate";
        Modelica.SIunits.Mass W_fluids[N_ch, N_cell_pc] "Fluid weights";
    protected
        final parameter Real phi = (1 + sqrt(1 + X^2) + 4*sqrt(1 + 0.5*X^2))/6
        "Corrugation pitch";

    equation
          for j in 1:N_ch loop
            for i in 1:N_cell_pc loop
              qdot_hot[j, i]     = pin_hot[j, i]*mdot_hot*(h_hot[j, i] -
              h_hot[j, i + 1]);
              qdot_cold[j, i]    = pin_cold[j, i]*mdot_cold*(h_cold[j, i + 1] -
              h_cold[j, i]);
              state_hot[j, i]    = Medium_hot.setState_ph(p_hot_in,
                0.5*(h_hot[j, i] + h_hot[j, i + 1]));
              state_cold[j, i]   = Medium_cold.setState_ph(p_cold_in,
                0.5*(h_cold[j, i] + h_cold[j, i + 1]));
            end for;
        end for;
      annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},
                {100,100}}), graphics={Polygon(
              points={{-28,-80},{32,-60},{32,80},{-28,60},{-28,-80}},
              lineColor={0,0,0},
              smooth=Smooth.None,
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid)}));
    end plate;
  end Objects;

  package Pressure_drops "A package containing pressure drops correlations"
    package Tubes "heat transfer correlations in tubes"
      class Frank "Frank correlation for tubes"
          extends Pressure_drops.Tubes.Base_classes.base_dp;
          output Modelica.SIunits.AbsolutePressure dp[Ncell]
          "Pressure drops cells";
          output Modelica.SIunits.AbsolutePressure dp_tot
          "Pressure drops tubes";
          input Modelica.SIunits.Length l "Lenght";
          input Modelica.SIunits.Length Dhyd "Hydraulic diameter";
          input Modelica.SIunits.AbsolutePressure heads
          "Number of velocity heads";
          Modelica.SIunits.ReynoldsNumber Re[Ncell] "Reynolds number";
          Modelica.SIunits.Velocity u[Ncell](start=ones(Ncell)) "Velocity";
          Real csi[Ncell] "Friction factor";

      equation
        for i in 1:Ncell loop
          u[i]     = mdot/state[i].d/Aflow;
          Re[i]    = Miscellanea.numbers.Reynolds(
                  u[i],
                  state[i].d,
                  state[i].eta,
                  Dhyd);
          if (Re[i] < 8e2) then
            csi[i] = 10^0.94244* Re[i]^(-1.03935);
          else
            csi[i] = 10^(-1.4249925) * Re[i]^(-0.225699);
          end if;
          dp[i]             =  0.5*(8*csi[i]*l/Dhyd + heads)*state[i].d*u[i]^2
            "Local pressure drops";
        end for;

        dp_tot = sum(dp) "Total pressure drops";

      end Frank;

      package Base_classes "Base classes for the pressure drops"
        class base_dp "Basic pressure drop correlation"
            replaceable package Medium = Media.OneRandomOrganicFluid constrainedby
            Modelica.Media.Interfaces.PartialMedium "Medium model" annotation(choicesAllMatching = true);
            parameter Integer Ncell(start=3) "Number of cell elements";
            input Modelica.SIunits.Area Aflow "Cross-sectional area";
            input Modelica.SIunits.MassFlowRate mdot "Mass flow rate";
            input Medium.ThermodynamicState state[Ncell];
        end base_dp;
      end Base_classes;
      annotation (Icon(graphics={
            Rectangle(extent={{-80,-40},{60,-40}}, lineColor={85,170,255}),
            Rectangle(extent={{-80,20},{60,20}}, lineColor={85,170,255}),
            Line(
              points={{60,20},{60,-40}},
              color={85,170,255},
              smooth=Smooth.None),
            Line(
              points={{-20,-40},{-24,-36}},
              color={85,170,255},
              smooth=Smooth.None),
            Line(
              points={{-20,-40},{-24,-44}},
              color={85,170,255},
              smooth=Smooth.None),
            Line(
              points={{-20,20},{-16,24}},
              color={85,170,255},
              smooth=Smooth.None),
            Line(
              points={{-20,20},{-16,16}},
              color={85,170,255},
              smooth=Smooth.None)}));
    end Tubes;

    package Shell "heat transfer correlations in shells"

      class single_phase_Kern
        "Kern correlation for shell single phase to single phase"
          extends Pressure_drops.Shell.Base_classes.base_Kern;
          Modelica.SIunits.AbsolutePressure dp[Ncell] "Pressure drops cells";
          Modelica.SIunits.AbsolutePressure dp_tot "Pressure drops tubes";
          input Modelica.SIunits.DynamicViscosity eta_wall[Ncell]
          "Viscosity of the fluid at the wall";
          Modelica.SIunits.ReynoldsNumber Re[Ncell](start=10e6*ones(Ncell))
          "Reynolds number";
          Modelica.SIunits.Velocity u[Ncell](start=ones(Ncell)) "Velocity";
          Real csi[Ncell] "Friction factor";

      equation
        for i in 1:Ncell loop
          u[i]     = mdot/state[i].d/Aflow;
          Re[i]    = Miscellanea.numbers.Reynolds(
                  u[i],
                  state[i].d,
                  state[i].eta,
                  Dhyd);
          if (Re[i] < 3e2) then
            csi[i] = 18.4810 * Re[i]^(-0.9344);
          else
            csi[i] = 0.2581 * Re[i]^(-0.1759);
          end if;
          dp[i]    =  4*csi[i]*d_s/Dhyd*l/l_b*state[i].d*u[i]^2;
        end for;

        dp_tot = sum(dp);

      end single_phase_Kern;

      class condensation_Kern "Kern correlation for shell condensation"
          extends Pressure_drops.Shell.Base_classes.base_Kern;
          parameter Real X "Factor for pressure drop in condensation";
          input Modelica.SIunits.DynamicViscosity eta_wall[Ncell]
          "exponent for the viscosity correction";
          input Modelica.SIunits.AbsolutePressure p_in "Inlet pressure";
          Modelica.SIunits.ReynoldsNumber Re[Ncell](start=10e6*ones(Ncell))
          "Reynolds number";
          Modelica.SIunits.Velocity u[Ncell](start=ones(Ncell)) "Velocity";
          Real csi[Ncell] "Friction factor";
          Modelica.SIunits.AbsolutePressure dp[Ncell] "Pressure drops cells";
          Modelica.SIunits.AbsolutePressure dp_tot "Pressure drops tubes";
          Medium.SaturationProperties sat "Saturation properties";
          Medium.ThermodynamicState state_l
          "Thermodynamic state in saturated liquid";
          Medium.ThermodynamicState state_v
          "Thermodynamic state in saturated vapor";

      equation
          sat      = Medium.setSat_p(p_in);
          state_l  = Medium.setState_ph(p_in, sat.hl);
          state_v  = Medium.setState_ph(p_in, sat.hv);
        for i in 1:Ncell loop
          u[i]     = mdot/state[i].d/Aflow "tube velocity";
          Re[i]    = Miscellanea.numbers.Reynolds(
                  u[i],
                  state[i].d,
                  state[i].eta,
                  Dhyd);
          if (Re[i] < 3e2) then
            csi[i] = 18.4810 * Re[i]^(-0.9344);
          else
            csi[i] = 0.2581 * Re[i]^(-0.1759);
          end if;
          if (state[i].h < state_l.h or state[i].h > state_v.h) then
            dp[i]  =  4*csi[i]*d_s/Dhyd*l/l_b*state[i].d*u[i]^2;
          else
            dp[i]  =  4*X*csi[i]*d_s/Dhyd*l/l_b*state[i].d*u[i]^2;
          end if;
        end for;

        dp_tot = sum(dp);

      end condensation_Kern;

      class single_phase_Johnston "Wills and Johnston correlation for shell"
          extends Pressure_drops.Shell.Base_classes.base_Kern;
          output Modelica.SIunits.AbsolutePressure dp[Ncell]
          "Pressure drops tubes";
          output Modelica.SIunits.AbsolutePressure dp_tot
          "Pressure drops tubes";
          input Modelica.SIunits.DynamicViscosity eta_wall[Ncell]
          "exponent for the viscosity correction";
          Modelica.SIunits.ReynoldsNumber Re[Ncell](start=10e6*ones(Ncell))
          "Reynolds number tubes (average)";
          Modelica.SIunits.Velocity u[Ncell](start=ones(Ncell))
          "Velocity (average)";
          Real csi[Ncell] "Friction factor";

      equation
        for i in 1:Ncell loop
          u[i]   = mdot/state[i].d/Aflow;
          Re[i]  = Miscellanea.numbers.Reynolds(
                  u[i],
                  state[i].d,
                  state[i].eta,
                  Dhyd);
          assert(Re[i] > 4e2, "Reynolds number is lower than 1e2 to use Wills_Johnston", AssertionLevel.warning);
          csi[i] = 1.7789*Re[i]^(-0.195868);
          dp[i]  = 0.5*csi[i]*d_s/Dhyd*l/l_b*(eta_wall[i]/state[i].eta)^0.14*state[i].d*u[i]^2;
        end for;

        dp_tot = sum(dp);

      end single_phase_Johnston;

      class single_phase_Bell_Delaware
        "Bell Delaware correlation for shell side for single phase to single phase"
          extends Pressure_drops.Shell.Base_classes.base_Bell_Delaware;
          input Modelica.SIunits.DynamicViscosity eta_wall[Ncell]
          "exponent for the viscosity correction";
          Modelica.SIunits.AbsolutePressure dp_c[Ncell]
          "Pressure drops cross flow zone";
          Modelica.SIunits.AbsolutePressure dp_w[Ncell]
          "Pressure drops window zone";
          Modelica.SIunits.AbsolutePressure dp_e[Ncell]
          "Pressure drops end wall zone";
          Modelica.SIunits.AbsolutePressure dp_w_tot
          "Pressure drops window zone total";
          Modelica.SIunits.AbsolutePressure dp_c_tot
          "Pressure drops cross flow zone total";
          Modelica.SIunits.AbsolutePressure dp_e_tot
          "Pressure drops end wall zone total";
          Modelica.SIunits.AbsolutePressure dp_tot "Pressure drops total";
          Modelica.SIunits.ReynoldsNumber Re[Ncell]
          "Reynolds number shell (average)";
          Modelica.SIunits.Velocity u[Ncell](start=ones(Ncell))
          "Velocities inside the shell";
          Modelica.SIunits.Velocity u_w[Ncell](start=ones(Ncell))
          "Velocities inside the shell";
          Real csi[Ncell] "Friction factor";

      protected
           Real b1[Ncell];
           Real b2[Ncell];
           Real b[Ncell];
      equation

        for i in 1:Ncell loop

          u[i]    = mdot/state[i].d/S_m;
          u_w[i]  = mdot/state[i].d/S_w;
          Re[i]   = Miscellanea.numbers.Reynolds(
                  u[i],
                  state[i].d,
                  state[i].eta,
                  Dhyd_o);

          if (Re[i] <= 1e1) then
            if (layout == 1) then
              b1[i]  = 48;
              b2[i]  = -1;
            else
              b1[i]  = 35;
              b2[i]  = -1;
            end if;
          elseif (Re[i] > 1e1 and Re[i] <= 1e2) then
            if (layout == 1) then
              b1[i]  = 45.1;
              b2[i]  = -0.973;
            else
              b1[i]  = 32.1;
              b2[i]  = -0.963;
            end if;
          elseif (Re[i] > 1e2 and Re[i] <= 1e3) then
            if (layout == 1) then
              b1[i]  = 4.57;
              b2[i]  = -0.476;
            else
              b1[i]  = 6.09;
              b2[i]  = -0.602;
            end if;
          elseif (Re[i] > 1e3 and Re[i] < 1e4) then
            if (layout == 1) then
              b1[i]  = 0.486;
              b2[i]  = -0.152;
            else
              b1[i]  = 0.0815;
              b2[i]  = 0.022;
            end if;
          else
            if (layout == 1) then
              b1[i]  = 0.372;
              b2[i]  = -0.123;
            else
              b1[i]  = 0.391;
              b2[i]  = -0.148;
            end if;
          end if;

          if (layout == 1) then
            b[i]   = 7/(1 + 0.14*Re[i]^0.5);
          else
            b[i]   = 6.3/(1 + 0.14*Re[i]^0.378);
          end if;

          csi[i] = b1[i]*(1.33/pitch_f)^b[i]*Re[i]^b2[i];

          if (i > Ncell - N_passes) then
            dp_w[i]  = 0;
          else
            dp_w[i]  = (2 + 0.6*N_tcw)*mdot^2/(2*S_w*S_m*state[i].d)*(N_baffles
            /(N_passes*N_baffles_d));
          end if;

          if (i <= N_passes or i > Ncell - N_passes) then
            dp_c[i] = 0;
            dp_e[i] = 2*(N_tcc + N_tcw)*csi[i]*(state[i].d*u[i]^2/2)*(N_baffles
            /(N_passes*N_baffles_d));
          else
            dp_e[i] = 0;
            dp_c[i] = ff*2*N_tcc*csi[i]*(state[i].d*u[i]^2/2)*(N_baffles
            /(N_passes*N_baffles_d));
          end if;

        end for;

        dp_w_tot = sum(dp_w)*R_l;
        dp_c_tot = sum(dp_c)*R_b*R_l;
        dp_e_tot = sum(dp_e)*R_b;

        dp_tot = dp_w_tot + dp_c_tot + dp_e_tot;

      end single_phase_Bell_Delaware;

      class condensation_Bell_Delaware
        "Bell Delaware correlation for shell side for condensation"
          extends Pressure_drops.Shell.Base_classes.base_Bell_Delaware;
          input Modelica.SIunits.AbsolutePressure p_in "Inlet pressure";
          Modelica.SIunits.AbsolutePressure dp_c[Ncell]
          "Pressure drops cross flow zone";
          Modelica.SIunits.AbsolutePressure dp_w[Ncell]
          "Pressure drops window zone";
          Modelica.SIunits.AbsolutePressure dp_e[Ncell]
          "Pressure drops end wall zone";
          Modelica.SIunits.AbsolutePressure dp_w_tot
          "Pressure drops window zone total";
          Modelica.SIunits.AbsolutePressure dp_c_tot
          "Pressure drops cross flow zone total";
          Modelica.SIunits.AbsolutePressure dp_e_tot
          "Pressure drops end wall zone total";
          Modelica.SIunits.AbsolutePressure dp_tot "Pressure drops total";
          Modelica.SIunits.ReynoldsNumber Re[Ncell] "Reynolds number";
          Modelica.SIunits.ReynoldsNumber Re_l(start=1e5)
          "Reynolds number if all liquid phase";
          Modelica.SIunits.Velocity u[Ncell](start=ones(Ncell))
          "Velocities inside the shell";
          Modelica.SIunits.Velocity u_w[Ncell](start=ones(Ncell))
          "Velocities inside the shell";
          Modelica.SIunits.Velocity u_l(start=1) "Velocity if all liquid phase";
          Real csi[Ncell] "Friction factor";
          Real csi_l "Friction factor if all liquid phase";
          Real xq[Ncell] "Vapor quality";
          Real phi_lo[Ncell] "Two phase multiplier cross flow";
          Real phi_lo_w[Ncell]
          "Two phase multiplier baffle window, horizontal cut";
          Real Y "Parameter cross flow";
          Real Y_w "Parameter baffle window, horizontal cut";
          Medium.SaturationProperties sat "Saturation properties";
          Medium.ThermodynamicState state_l
          "Thermodynamic state in saturated liquid";
          Medium.ThermodynamicState state_v
          "Thermodynamic state in saturated vapor";
      protected
           Real b1[Ncell];
           Real b2[Ncell];
           Real b[Ncell];
           Real b1_l;
           Real b2_l;
           Real b_l;
      equation
        sat       = Medium.setSat_p(p_in);
        state_l   = Medium.setState_ph(p_in, sat.hl);
        state_v   = Medium.setState_ph(p_in, sat.hv);
        u_l       = mdot/state_l.d/S_m;
        Y         = sqrt(state_l.d/state_v.d)*(state_v.eta/state_l.eta)^0.23;
        Y_w       = sqrt(state_l.d/state_v.d);
        Re_l      = Miscellanea.numbers.Reynolds(u_l, state_l.d, state_l.eta, Dhyd_o);

          if (Re_l <= 1e1) then
            if (layout == 1) then
              b1_l  = 48;
              b2_l  = -1;
            else
              b1_l  = 35;
              b2_l  = -1;
            end if;
          elseif (Re_l > 1e1 and Re_l <= 1e2) then
            if (layout == 1) then
              b1_l  = 45;
              b2_l  = -0.973;
            else
              b1_l  = 32.1;
              b2_l  = -0.963;
            end if;
          elseif (Re_l > 1e2 and Re_l <= 1e3) then
            if (layout == 1) then
              b1_l  = 4.57;
              b2_l  = -0.476;
            else
              b1_l  = 6.09;
              b2_l  = -0.602;
            end if;
          elseif (Re_l > 1e3 and Re_l < 1e4) then
            if (layout == 1) then
              b1_l  = 0.486;
              b2_l  = -0.152;
            else
              b1_l  = 0.0815;
              b2_l  = 0.022;
            end if;
          else
            if (layout == 1) then
              b1_l  = 0.372;
              b2_l  = -0.123;
            else
              b1_l  = 0.391;
              b2_l  = -0.148;
            end if;
          end if;

          if (layout == 1) then
            b_l   = 7/(1 + 0.14*Re_l^0.5);
          else
            b_l   = 6.3/(1 + 0.14*Re_l^0.378);
          end if;

          csi_l = b1_l*(1.33/pitch_f)^b_l*Re_l^b2_l;

        for i in 1:Ncell loop

          u[i]    = mdot/state[i].d/S_m;
          u_w[i]  = mdot/state[i].d/S_w;
          Re[i]   = Miscellanea.numbers.Reynolds(
                  u[i],
                  state[i].d,
                  state[i].eta,
                  Dhyd_o);

          if (Re[i] <= 1e1) then
            if (layout == 1) then
              b1[i]  = 48;
              b2[i]  = -1;
            else
              b1[i]  = 35;
              b2[i]  = -1;
            end if;
          elseif (Re[i] > 1e1 and Re[i] <= 1e2) then
            if (layout == 1) then
              b1[i]  = 45;
              b2[i]  = -0.973;
            else
              b1[i]  = 32.1;
              b2[i]  = -0.963;
            end if;
          elseif (Re[i] > 1e2 and Re[i] <= 1e3) then
            if (layout == 1) then
              b1[i]  = 4.57;
              b2[i]  = -0.476;
            else
              b1[i]  = 6.09;
              b2[i]  = -0.602;
            end if;
          elseif (Re[i] > 1e3 and Re[i] < 1e4) then
            if (layout == 1) then
              b1[i]  = 0.486;
              b2[i]  = -0.152;
            else
              b1[i]  = 0.0815;
              b2[i]  = 0.022;
            end if;
          else
            if (layout == 1) then
              b1[i]  = 0.372;
              b2[i]  = -0.123;
            else
              b1[i]  = 0.391;
              b2[i]  = -0.148;
            end if;
          end if;

          if (layout == 1) then
            b[i]   = 7/(1 + 0.14*Re[i]^0.5);
          else
            b[i]   = 6.3/(1 + 0.14*Re[i]^0.378);
          end if;

          csi[i] = b1[i]*(1.33/pitch_f)^b[i]*Re[i]^b2[i];

          if (state[i].h < state_l.h or state[i].h > state_v.h) then

            if (i > Ncell - N_passes) then
              dp_w[i]    = 0;
            else
              dp_w[i]    = (2 + 0.6*N_tcw)*mdot^2/(2*S_w*S_m*state[i].d)*(N_baffles/
              (N_passes*N_baffles_d));
            end if;

            if (i <= N_passes or i > Ncell - N_passes) then
              dp_c[i]   = 0;
              dp_e[i]   = 2*(N_tcc + N_tcw)*csi[i]*(state[i].d*u[i]^2/2)*(N_baffles/
              (N_passes*N_baffles_d));
            else
              dp_e[i]   = 0;
              dp_c[i]   = ff*2*N_tcc*csi[i]*(state[i].d*u[i]^2/2)*(N_baffles/
              (N_passes*N_baffles_d));
            end if;

            xq[i]       = 1;
            phi_lo[i]   = 1;
            phi_lo_w[i] = 1;

          else

            xq[i]       = (state[i].h - state_l.h)/(state_v.h - state_l.h);
            phi_lo[i]   = 1 + (Y^2 - 1)*(0.75*(xq[i]*(1 - xq[i]))^0.77 + xq[i]^1.54);
            phi_lo_w[i] = 1 + (Y_w^2 - 1)*(2*xq[i]/(Y_w + 1));

            if (i > Ncell - N_passes) then
              dp_w[i]   = 0;
            else
              dp_w[i]   = phi_lo_w[i]*(2 + 0.6*N_tcw)*mdot^2/(2*S_w*S_m*state_l.d)*
              (N_baffles/(N_passes*N_baffles_d));
            end if;

            if (i <= N_passes or i > Ncell - N_passes) then
              dp_c[i]   = 0;
              dp_e[i]   = phi_lo[i]*(N_tcc + N_tcw)*csi_l*(state_l.d*u_l^2)*
              (N_baffles/(N_passes*N_baffles_d));
            else
              dp_e[i]   = 0;
              dp_c[i]   = phi_lo[i]*ff*N_tcc*csi_l*(state_l.d*u_l^2)*(N_baffles/
              (N_passes*N_baffles_d));
            end if;
          end if;
        end for;

        dp_w_tot = sum(dp_w)*R_l;
        dp_c_tot = sum(dp_c)*R_b*R_l;
        dp_e_tot = sum(dp_e)*R_b;

        dp_tot = dp_w_tot + dp_c_tot + dp_e_tot;

      end condensation_Bell_Delaware;

      package Base_classes "Base classes for the pressure drops"
        class base_dp "Basic pressure drop correlation"
            replaceable package Medium = Media.OneRandomOrganicFluid constrainedby
            Modelica.Media.Interfaces.PartialMedium "Medium model" annotation(choicesAllMatching = true);
            parameter Integer Ncell(start=3) "Number of cell elements";
            input Modelica.SIunits.Area Aflow "Cross-sectional area";
            input Modelica.SIunits.MassFlowRate mdot "Mass flow rate";
            input Medium.ThermodynamicState state[Ncell];
        end base_dp;

        class base_Bell_Delaware
          "Bell Delaware correlation for shell side base class"
            extends Pressure_drops.Shell.Base_classes.base_dp;
            input Modelica.SIunits.Length Dhyd "Hydraulic diameter";
            input Modelica.SIunits.Length Dhyd_o "Outer hydraulic diameter";
            parameter Real b_cut "Baffle cut";
            parameter Real ttb = 0.8e-3 "tube to baffle clearance";
            parameter Real bts = 4e-3 "tube to baffle clearance";
            parameter Real pitch_f "Tube pitch";
            parameter Real N_ss
            "The number of sealing strips (pairs) in one baffle spacing";
            parameter Integer layout "Tube layout 1 = triangular, 2 = squared";
            parameter Integer N_passes "Number of tube passes";
            parameter Real ff "Fouling factor (=1 if clean)";
            parameter Integer N_baffles "Number of baffles";
            parameter Integer N_baffles_d
            "The number of discretization volumes is always N_baffles + 1";
            input Modelica.SIunits.Length d_s "Shell diameter";
            input Modelica.SIunits.Length d_b "Bundle diameter";
            input Real N_tubes "Number of tubes in the bundle";
            input Modelica.SIunits.Length l_b "Baffle lenght";
            Modelica.SIunits.Length tpv "vertical tube pitch";
            Real N_tcc
            "Number of constrictions crossed=number of tube rows between the baffle tips";
            Real N_tcw
            "The effective number of tube rows crossed in the baffle window";
            Real F_sbp
            "The ratio of the bypass area, Sb, to the overall crossflow area";
            Real F_w "Fraction of number of tubes in the baffle window";
            Real R_l
            "Correction factors for baffle leakage effects for pressure drops";
            Real R_b
            "Correction factors for bundle bypass effects for pressure drops";
            Modelica.SIunits.Area S_w "Window area";
            Modelica.SIunits.Area S_tb "The tube-to-baffle leakage area";
            Modelica.SIunits.Area S_sb "The shell-to-baffle leakage area";
            Modelica.SIunits.Area S_m "Reference normal area";
            Modelica.SIunits.Area S_wg "The gross window flow area";
            Modelica.SIunits.Area S_wt
            "The segmental baffle window area occupied by the tubes";
            Modelica.SIunits.Area S_b
            "The bypass area between the shell and the tube bundle within one baffle";
            Real r_s "Correlational parameter";
            Real r_lm "Correlational parameter";
            Modelica.SIunits.Angle teta_ds "The centriangle of baffle cut";
            Modelica.SIunits.Angle teta_ctl
            "The upper centriangle of baffle cut";

        equation
          //Tube row correction factor F_n
          if (layout == 1) then
            tpv = 0.866*pitch_f*Dhyd_o;
          else
            tpv = pitch_f*Dhyd_o;
          end if;

          teta_ctl = 2*acos(d_s/(d_b - Dhyd_o)*(1 - 2*b_cut));
          teta_ds  = 2*acos(1 - 2*b_cut);
          N_tcc    = d_s*(1 - 2*b_cut)/tpv;
          S_m      = l_b*(d_s - d_b + (d_b - Dhyd_o)*(1 - 1/pitch_f));
          S_wg     = 0.25*pi*d_s^2*(teta_ds/(2*pi) - sin(teta_ds)/(2*pi));
          F_w      = teta_ctl/(2*pi) - sin(teta_ctl)/(2*pi);
          S_wt     = 0.25*N_tubes*F_w*pi*Dhyd_o^2;

          S_tb     = 0.25*pi*((ttb + Dhyd_o)^2 - Dhyd_o^2)*N_tubes*(1 - F_w);
          S_sb     = 0.5*bts*d_s*(1 - teta_ds/(2*pi));
          r_s      = S_sb/(S_tb + S_sb);
          r_lm     = (S_tb + S_sb)/S_m;
          R_l      = exp(-1.33*(1 + r_s)*r_lm^(-0.15*(1 + r_s) + 0.8));

          S_b      = l_b*(d_s - d_b);
          F_sbp    = S_b/S_m;
          R_b      = exp(-3.7*F_sbp*(1 - (2/N_ss)^(1/3)));

          S_w      = S_wg - S_wt;
          N_tcw    = 0.8*(d_s*b_cut - 0.5*(d_s - (d_b - Dhyd_o)))/tpv;

        end base_Bell_Delaware;

        class base_Kern "Kern correlation for shell base class"
            extends Pressure_drops.Shell.Base_classes.base_dp;
            parameter Modelica.SIunits.Length l "Lenght";
            input Modelica.SIunits.Length Dhyd "Hydraulic diameter";
            input Modelica.SIunits.Length d_s "Shell diameter";
            input Modelica.SIunits.Length l_b "Baffle length";
        end base_Kern;
      end Base_classes;
      annotation (Icon(graphics={
            Rectangle(
              extent={{-80,40},{80,-40}},
              lineColor={0,0,0}),
            Rectangle(
              extent={{-70,-40},{-50,-60}},
              lineColor={0,0,0}),
            Rectangle(
              extent={{50,60},{70,40}},
              lineColor={0,0,0}),
            Line(
              points={{60,70},{60,46}},
              color={255,0,0},
              smooth=Smooth.None),
            Line(
              points={{2,-2},{-2,2}},
              color={255,0,0},
              smooth=Smooth.None,
              origin={62,48},
              rotation=90),
            Line(
              points={{2,2},{-2,-2}},
              color={255,0,0},
              smooth=Smooth.None,
              origin={58,48},
              rotation=90),
            Line(
              points={{-60,-42},{-60,-66}},
              color={255,0,0},
              smooth=Smooth.None),
            Line(
              points={{2,-2},{-2,2}},
              color={255,0,0},
              smooth=Smooth.None,
              origin={-58,-64},
              rotation=90),
            Line(
              points={{2,2},{-2,-2}},
              color={255,0,0},
              smooth=Smooth.None,
              origin={-62,-64},
              rotation=90)}));
    end Shell;

    package Plates
      class single_phase_Martin
        "Martin pressure drop correlation for single phase"
          extends Base_classes.base_dp;
          parameter Modelica.SIunits.Angle beta "Plate inclination angle";
          Modelica.SIunits.ReynoldsNumber Re[N_ch, N_cell_pc] "Reynolds number";
          Modelica.SIunits.Velocity u[N_ch, N_cell_pc] "Velocity";
          Modelica.SIunits.Velocity u_pt[N_passes] "Port velocity";
          Modelica.SIunits.AbsolutePressure dp[N_ch, N_cell_pc]
          "Pressure drops cells";
          Modelica.SIunits.AbsolutePressure dp_pt[N_passes]
          "Pressure drops ports";
          Modelica.SIunits.AbsolutePressure dp_plates[N_ch]
          "Pressure drops plates";
          Modelica.SIunits.AbsolutePressure dp_tot "Total pressure drops";

      protected
          Real csi[N_ch, N_cell_pc] "Friction factor";
          Real csi0[N_ch, N_cell_pc] "Friction factor 0";
          Real csi1[N_ch, N_cell_pc] "Friction factor 1";
      equation
          for j in 1:N_ch loop
            dp_plates[j]   = sum(dp[j,:]);
            for i in 1:N_cell_pc loop
              u[j, i]      = mdot/state[j, i].d/Aflow;
              Re[j, i]     =Design.Miscellanea.numbers.Reynolds(
                    u[j, i],
                    state[j, i].d,
                    state[j, i].eta,
                    Dhyd);
              if (Re[j, i] < 2e3) then
                csi0[j, i] = 16/Re[j, i];
                csi1[j, i] = 149.25/Re[j, i] + 0.9625;
              else
                csi0[j, i] = 1/(1.56*log(Re[j, i]) - 3)^2;
                csi1[j, i] = 9.75/Re[j, i]^0.289;
              end if;
              csi[j, i]    = 1/((1 - cos(beta))/sqrt(3.8*csi1[j, i]) + cos(beta)/
              sqrt(0.045*tan(beta) + 0.09*sin(beta) + csi0[j, i]/cos(beta)))^2;
              dp[j, i]     = 2*csi[j, i]*(l/Dhyd)*state[j, i].d*u[j, i]^2;
            end for;
          end for;

          for j in 1:N_passes loop
            u_pt[j]     = N_ch_p*mdot/(state_p[j].d*A_pt);
            dp_pt[j]    = 0.5*1.3*state_p[j].d*u_pt[j]^2;
          end for;

          dp_tot        = sum(dp_plates)/N_ch_p + sum(dp_pt);

      end single_phase_Martin;

      class evaporation_Martin
        "Lockhart-Martinelli pressure drop correlation for evaporation heat transfer. Martin for single phase."
          extends Base_classes.base_dp;
          parameter Modelica.SIunits.Angle beta "Plate inclination angle";
          input Modelica.SIunits.AbsolutePressure p_in "Inlet pressure";
          Modelica.SIunits.ReynoldsNumber Re[N_ch, N_cell_pc] "Reynolds number";
          Modelica.SIunits.ReynoldsNumber Re_l[N_ch, N_cell_pc]
          "Reynolds number liquid phase";
          Modelica.SIunits.ReynoldsNumber Re_v[N_ch, N_cell_pc]
          "Reynolds number vapor phase";
          Modelica.SIunits.Velocity u[N_ch, N_cell_pc] "Velocity";
          Modelica.SIunits.Velocity u_pt[N_passes] "Port velocity";
          Real xq[N_ch, N_cell_pc] "Vapor quality";
          Modelica.SIunits.AbsolutePressure dp[N_ch, N_cell_pc]
          "Pressure drops cells";
          Modelica.SIunits.AbsolutePressure dp_pt[N_passes]
          "Pressure drops ports";
          Modelica.SIunits.AbsolutePressure dp_plates[N_ch]
          "Pressure drops plates";
          Modelica.SIunits.AbsolutePressure dp_tot "Total pressure drops";
          Medium.SaturationProperties sat "Saturation properties";
          Medium.ThermodynamicState state_l
          "Thermodynamic state in saturated liquid";
          Medium.ThermodynamicState state_v
          "Thermodynamic state in saturated vapor";
      protected
          Modelica.SIunits.AbsolutePressure dp_l[N_ch, N_cell_pc]
          "Pressure drops liquid phase";
          Modelica.SIunits.AbsolutePressure dp_v[N_ch, N_cell_pc]
          "Pressure drops vapor phase";
          Real X_tt[N_ch, N_cell_pc] "Martinelli parameter";
          Real phi_l[N_ch, N_cell_pc] "Multiplier liquid phase";
          Real phi_v[N_ch, N_cell_pc] "Multiplier vapor phase";
          Modelica.SIunits.Velocity u_l[N_ch, N_cell_pc]
          "Velocity liquid phase";
          Modelica.SIunits.Velocity u_v[N_ch, N_cell_pc] "Velocity vapor phase";
          Real csi[N_ch, N_cell_pc] "Friction factor";
          Real csi_l[N_ch, N_cell_pc] "Friction factor liquid phase";
          Real csi_v[N_ch, N_cell_pc] "Friction factor vapor phase";
          Real csi0[N_ch, N_cell_pc] "Friction factor 0";
          Real csi1[N_ch, N_cell_pc] "Friction factor 1";
          Real C[N_ch, N_cell_pc] "Constant for the two phase multiplier";
      equation

          sat     = Medium.setSat_p(p_in);
          state_l = Medium.setState_ph(p_in, sat.hl);
          state_v = Medium.setState_ph(p_in, sat.hv);

          for j in 1:N_ch loop
            dp_plates[j]   = sum(dp[j,:]);
            for i in 1:N_cell_pc loop
              u[j, i]      = mdot/state[j, i].d/Aflow;
              Re[j, i]     =Design.Miscellanea.numbers.Reynolds(
                    u[j, i],
                    state[j, i].d,
                    state[j, i].eta,
                    Dhyd);
              if (Re[j, i] < 2e3) then
                csi0[j, i] = 16/Re[j, i];
                csi1[j, i] = 149.25/Re[j, i] + 0.9625;
              else
                csi0[j, i] = 1/(1.56*log(Re[j, i]) - 3)^2;
                csi1[j, i] = 9.75/Re[j, i]^0.289;
              end if;
              csi[j, i]    = 1/((1 - cos(beta))/sqrt(3.8*csi1[j, i]) + cos(beta)/
              sqrt(0.045*tan(beta) + 0.09*sin(beta) + csi0[j, i]/cos(beta)))^2;
              if (state[j, i].h < state_l.h or state[j, i].h > state_v.h) then
                xq[j, i]     = 0;
                u_l[j, i]    = 0;
                u_v[j, i]    = 0;
                Re_l[j, i]   = 0;
                Re_v[j, i]   = 0;
                X_tt[j, i]   = 0;
                C[j, i]      = 0;
                phi_l[j, i]  = 0;
                phi_v[j, i]  = 0;
                dp_l[j, i]   = 0;
                dp_v[j, i]   = 0;
                csi_l[j, i]  = 0;
                csi_v[j, i]  = 0;
                dp[j, i]     = 2*csi[j, i]*(l/Dhyd)*state[j, i].d*u[j, i]^2;
              else
                xq[j, i]     = (state[j, i].h - state_l.h)/(state_v.h - state_l.h);
                u_l[j, i]    = mdot*(1 - xq[j, i])/state_l.d/Aflow;
                u_v[j, i]    = mdot*xq[j, i]/state_v.d/Aflow;
                Re_l[j, i]   =Design.Miscellanea.numbers.Reynolds(
                      u_l[j, i],
                      state_l.d,
                      state_l.eta,
                      Dhyd);
                Re_v[j, i]   =Design.Miscellanea.numbers.Reynolds(
                      u_v[j, i],
                      state_v.d,
                      state_v.eta,
                      Dhyd);
                if (Re_l[j, i] <= 1e3) then
                  csi_l[j, i] = 16/Re_l[j, i];
                elseif (Re_l[j, i] > 1e3 and Re_l[j, i] <= 2e3) then
                  csi_l[j, i] = (1e-3*Re_l[j, i] - 1)*(16/Re_l[j, i]) + (2 - 1e-3*
                  Re_l[j, i])*4.6e-2/Re_l[j, i]^0.2;
                else
                  csi_l[j, i] = 4.6e-2/Re_l[j, i]^0.2;
                end if;
                if (Re_v[j, i] <= 1e3) then
                  csi_v[j, i] = 16/Re_v[j, i];
                elseif (Re_v[j, i] > 1e3 and Re_v[j, i] <= 2e3) then
                  csi_v[j, i] = (1e-3*Re_v[j, i] - 1)*(16/Re_v[j, i]) + (2 - 1e-3*
                  Re_v[j, i])*4.6e-2/Re_v[j, i]^0.2;
                else
                  csi_v[j, i] = 4.6e-2/Re_v[j, i]^0.2;
                end if;
                dp_l[j, i]   = 2*csi_l[j, i]*(l/Dhyd)*state_l.d*u_l[j, i]^2;
                dp_v[j, i]   = 2*csi_v[j, i]*(l/Dhyd)*state_v.d*u_v[j, i]^2;
                X_tt[j, i]   = sqrt(dp_l[j, i]/dp_v[j, i]);
                if (Re_l[j, i] > 1.5e3 and Re_v[j, i] > 1.5e3) then
                  C[j, i]    = 20;
                elseif (Re_l[j, i] <= 1.5e3 and Re_v[j, i] > 1.5e3) then
                  C[j, i]    = 12;
                elseif (Re_l[j, i] > 1.5e3 and Re_v[j, i] <= 1.5e3) then
                  C[j, i]    = 10;
                else
                  C[j, i]    = 5;
                end if;
                phi_l[j, i]  = 1 + C[j, i]/X_tt[j, i] + 1/X_tt[j, i]^2;
                phi_v[j, i]  = 1 + C[j, i]*X_tt[j, i] + X_tt[j, i]^2;
                if (dp_l[j, i]*phi_l[j, i] > dp_v[j, i]*phi_v[j, i]) then
                  dp[j, i] = dp_v[j, i]*phi_v[j, i];
                else
                  dp[j, i] = dp_l[j, i]*phi_l[j, i];
                end if;
              end if;
            end for;
          end for;

          for j in 1:N_passes loop
            u_pt[j]     = N_ch_p*mdot/(state_p[j].d*A_pt);
            dp_pt[j]    = 0.5*1.3*state_p[j].d*u_pt[j]^2;
          end for;

          dp_tot        = sum(dp_plates)/N_ch_p + sum(dp_pt);

      end evaporation_Martin;

      class single_phase_Coulson
        "Coulson pressure drop correlation for single phase"
          extends Base_classes.base_dp;
          Modelica.SIunits.ReynoldsNumber Re[N_ch, N_cell_pc] "Reynolds number";
          Modelica.SIunits.Velocity u[N_ch, N_cell_pc] "Velocity";
          Modelica.SIunits.Velocity u_pt[N_passes] "Port velocity";
          Modelica.SIunits.AbsolutePressure dp[N_ch, N_cell_pc]
          "Pressure drops cells";
          Modelica.SIunits.AbsolutePressure dp_pt[N_passes]
          "Pressure drops ports";
          Modelica.SIunits.AbsolutePressure dp_plates[N_ch]
          "Pressure drops plates";
          Modelica.SIunits.AbsolutePressure dp_tot "Total pressure drops";
      protected
          Real csi[N_ch, N_cell_pc] "Friction factor";
      equation

          for j in 1:N_ch loop
            dp_plates[j]   = sum(dp[j,:]);
            for i in 1:N_cell_pc loop
              u[j, i]   = mdot/state[j, i].d/Aflow;
              Re[j, i]  =Design.Miscellanea.numbers.Reynolds(
                    u[j, i],
                    state[j, i].d,
                    state[j, i].eta,
                    Dhyd);
              csi[j, i] = 0.6/Re[j, i]^0.3;
              dp[j, i]  = 4*csi[j, i]*(l/Dhyd)*state[j, i].d*u[j, i]^2;
            end for;
          end for;

          for j in 1:N_passes loop
            u_pt[j]     = N_ch_p*mdot/(state_p[j].d*A_pt);
            dp_pt[j]    = 0.5*1.3*state_p[j].d*u_pt[j]^2;
          end for;

          dp_tot        = sum(dp_plates)/N_ch_p + sum(dp_pt);

      end single_phase_Coulson;

      package Base_classes "Base classes for the heat transfer"
        class base_dp "Basic pressure drop correlation"
            replaceable package Medium = Media.OneRandomOrganicFluid constrainedby
            Modelica.Media.Interfaces.PartialMedium "Medium model" annotation(choicesAllMatching = true);
            parameter Integer N_cell_pc( start = 3)
            "Number of cells per channel";
            parameter Integer N_ch(start = 3) "Total number of channels";
            parameter Modelica.SIunits.Area A_pt "Port area";
            parameter Modelica.SIunits.Length Dhyd "Hydraulic diameter";
            parameter Integer N_ch_p(start = 3) "Number of channels per pass";
            parameter Integer N_passes(start = 1) "Number of passes";
            input Modelica.SIunits.MassFlowRate mdot "Mass flow rate";
            input Medium.ThermodynamicState state_p[N_passes] "Index matrix";
            input Modelica.SIunits.Area Aflow "Cross-sectional area";
            input Modelica.SIunits.Length l "Length";
            input Medium.ThermodynamicState state[N_ch, N_cell_pc]
            "Thermodynamic states";
        end base_dp;
      end Base_classes;
      annotation (Icon(graphics={      Polygon(
              points={{-32,-80},{28,-60},{28,80},{-32,60},{-32,-80}},
              lineColor={0,0,0},
              smooth=Smooth.None,
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid)}));
    end Plates;
  annotation (Icon(graphics={
        Rectangle(
          extent={{70,88},{90,-88}},
          lineColor={95,95,95},
          fillColor={215,215,215},
          fillPattern=FillPattern.Solid),
        Polygon(
          points={{-20,-2},{-20,20},{20,0},{-20,-20},{-20,-2}},
          lineColor={0,128,255},
          smooth=Smooth.None,
          fillColor={0,128,255},
          fillPattern=FillPattern.Solid,
          origin={33,52},
          rotation=90),
        Rectangle(
          extent={{-50,6},{50,-6}},
          lineColor={0,128,255},
          fillColor={0,128,255},
          fillPattern=FillPattern.Solid,
          origin={33,-18},
          rotation=90),
        Polygon(
          points={{-20,-2},{-20,20},{20,0},{-20,-20},{-20,-2}},
          lineColor={0,128,255},
          smooth=Smooth.None,
          fillColor={0,128,255},
          fillPattern=FillPattern.Solid,
          origin={-41,52},
          rotation=90),
        Rectangle(
          extent={{-50,6},{50,-6}},
          lineColor={0,128,255},
          fillColor={0,128,255},
          fillPattern=FillPattern.Solid,
          origin={-41,-18},
          rotation=90)}));
  end Pressure_drops;

  package Heat_transfer "A package containing heat transfer correlations"
    package Tubes "heat transfer correlations in tubes"
      class Dittus_Boelter "Dittus Boelter correlation for tubes"
          extends Base_classes.base_ht;
          Modelica.SIunits.NusseltNumber Nu[Ncell] "Nusselt number";
          Modelica.SIunits.CoefficientOfHeatTransfer ht[Ncell]
          "Heat transfer coefficient";
          parameter Modelica.SIunits.Length Dhyd "Hydraulic diameter";
          parameter Real alfa "exponent for the Prandtl number";
          Modelica.SIunits.ReynoldsNumber Re[Ncell] "Reynolds number";
          Modelica.SIunits.PrandtlNumber Pr[Ncell] "Prandtl number";
          Modelica.SIunits.Velocity u[Ncell] "Velocity";
      equation
            for i in 1:Ncell loop
              u[i]   = mdot/state[i].d/Aflow;
              Re[i]  =Miscellanea.numbers.Reynolds(
                  u[i],
                  state[i].d,
                  state[i].eta,
                  Dhyd);
              Pr[i]  =Miscellanea.numbers.Prandtl(
                  state[i].cp,
                  state[i].eta,
                  state[i].lambda);
              assert(Re[i] > 1e4, "Reynolds number is lower than 1e4 to use Dittus and Boelter",
                AssertionLevel.warning);
              Nu[i]  =  2.3e-2*Re[i]^0.8*Pr[i]^alfa;
              ht[i]  =Miscellanea.numbers.Nusselt(
                  Nu[i],
                  state[i].lambda,
                  Dhyd);
            end for;

      end Dittus_Boelter;

      class Sieder_Tate "Sieder Tate correlation for tubes"
          extends Heat_transfer.Tubes.Base_classes.base_ht;
          Modelica.SIunits.NusseltNumber Nu[Ncell] "Nusselt number";
          Modelica.SIunits.CoefficientOfHeatTransfer ht[Ncell]
          "Heat transfer coefficient";
          parameter Modelica.SIunits.Length Dhyd "Hydraulic diameter";
          input Modelica.SIunits.DynamicViscosity eta_wall[Ncell]
          "exponent for the viscosity correction";
          Modelica.SIunits.ReynoldsNumber Re[Ncell] "Reynolds number";
          Modelica.SIunits.PrandtlNumber Pr[Ncell] "Prandtl number";
          Modelica.SIunits.Velocity u[Ncell] "Velocity";
      equation
            for i in 1:Ncell loop
              u[i]  = mdot/state[i].d/Aflow;
              Re[i]  =Miscellanea.numbers.Reynolds(
                  u[i],
                  state[i].d,
                  state[i].eta,
                  Dhyd);
              Pr[i]  =Miscellanea.numbers.Prandtl(
                  state[i].cp,
                  state[i].eta,
                  state[i].lambda);
              assert(Re[i] > 1e4,"Reynolds number is lower than 1e4 to use Sieder and Tate",
                AssertionLevel.warning);
              assert(Pr[i] > 0.6, "Prandtl number is lower than 0.6 to be use Sieder and Tate",
                AssertionLevel.warning);
              Nu[i] = 2.3e-2*Re[i]^0.8*Pr[i]^(1.0/3)*(eta_wall[i]
              /state[i].eta)^0.14;
              ht[i]  =Miscellanea.numbers.Nusselt(
                  Nu[i],
                  state[i].lambda,
                  Dhyd);
            end for;

      end Sieder_Tate;

      class Gnielinski "Gnielinski correlation for tubes"
          extends Heat_transfer.Tubes.Base_classes.base_ht;
          parameter Modelica.SIunits.Length l "Lenght";
          parameter Modelica.SIunits.Length Dhyd "Hydraulic diameter";
          Modelica.SIunits.ReynoldsNumber Re[Ncell] "Reynolds number";
          Modelica.SIunits.PrandtlNumber Pr[Ncell] "Prandtl number";
          Modelica.SIunits.Velocity u[Ncell] "Velocity";
          Modelica.SIunits.NusseltNumber Nu[Ncell] "Nusselt number";
          Modelica.SIunits.CoefficientOfHeatTransfer ht[Ncell]
          "Heat transfer coefficient";
      protected
          Real csi[Ncell] "Friction factor";
      equation
            for i in 1:Ncell loop
              u[i]   = mdot/state[i].d/Aflow;
              Re[i]  =Miscellanea.numbers.Reynolds(
                  u[i],
                  state[i].d,
                  state[i].eta,
                  Dhyd);
              Pr[i]  =Miscellanea.numbers.Prandtl(
                  state[i].cp,
                  state[i].eta,
                  state[i].lambda);
              csi[i] = 1/(0.78*log(Re[i]) - 1.5)^2;
              Nu[i]  = ((csi[i]/8)*Re[i]*Pr[i])/(1 + 12.7*sqrt(csi[i]/8)*(Pr[i]^
              (2.0/3) - 1))*(1 +  (Dhyd/l)^(2.0/3));
              ht[i]  =Miscellanea.numbers.Nusselt(
                  Nu[i],
                  state[i].lambda,
                  Dhyd);
            end for;
      end Gnielinski;

      class EagleFerguson
        "Eagle-Ferguson heat transfer correlation (only for liquid water)"
          extends Heat_transfer.Tubes.Base_classes.base_ht;
          Modelica.SIunits.NusseltNumber Nu[Ncell] "Nusselt number tubes ";
          parameter Modelica.SIunits.Length Dhyd "Hydraulic diameter";
          Modelica.SIunits.Velocity u[Ncell] "Velocity";
          Modelica.SIunits.CoefficientOfHeatTransfer ht[Ncell]
          "Heat transfer coefficient tubes";

      equation
            for i in 1:Ncell loop
              u[i]  = mdot/state[i].d/Aflow;
              Nu[i] =  4.2e3*Dhyd/state[i].lambda*(1.35 + 2e-2*(state[i].T - 273.15))
              *u[i]^0.8/(1e3*Dhyd)^0.2;
              ht[i]  =Miscellanea.numbers.Nusselt(
                  Nu[i],
                  state[i].lambda,
                  Dhyd);
            end for;
      end EagleFerguson;

      package Base_classes "Base classes for the heat transfer"
        class base_ht "Basic heat transfer correlation"
            replaceable package Medium = Media.OneRandomOrganicFluid constrainedby
            Modelica.Media.Interfaces.PartialMedium "Medium model" annotation(choicesAllMatching = true);
            parameter Integer Ncell(start=3) "Number of cell elements";
            input Modelica.SIunits.Area Aflow
            "Cross-sectional area (single tube)";
            input Modelica.SIunits.MassFlowRate mdot "Mass flow rate";
            input Medium.ThermodynamicState state[Ncell] "Thermodynamic state";
        end base_ht;
      end Base_classes;
      annotation (Icon(graphics={
            Rectangle(extent={{-80,-40},{60,-40}}, lineColor={85,170,255}),
            Rectangle(extent={{-80,20},{60,20}}, lineColor={85,170,255}),
            Line(
              points={{60,20},{60,-40}},
              color={85,170,255},
              smooth=Smooth.None),
            Line(
              points={{-20,-40},{-24,-36}},
              color={85,170,255},
              smooth=Smooth.None),
            Line(
              points={{-20,-40},{-24,-44}},
              color={85,170,255},
              smooth=Smooth.None),
            Line(
              points={{-20,20},{-16,24}},
              color={85,170,255},
              smooth=Smooth.None),
            Line(
              points={{-20,20},{-16,16}},
              color={85,170,255},
              smooth=Smooth.None)}));
    end Tubes;

    package Shell "heat transfer correlations in shells"

      class single_phase_Kern
        "Kern method for shell side only single phase to single phase"
          extends Heat_transfer.Shell.Base_classes.base_Kern;
          Modelica.SIunits.NusseltNumber Nu[Ncell] "Nusselt number";
          Modelica.SIunits.CoefficientOfHeatTransfer ht[Ncell]
          "Heat transfer coefficient";
          Modelica.SIunits.ReynoldsNumber Re[Ncell] "Reynolds number";
          Modelica.SIunits.PrandtlNumber Pr[Ncell] "Prandtl number ";
          Modelica.SIunits.Length d_s_eq "Equivalent shell diameter";
          Modelica.SIunits.Velocity u[Ncell] "Velocity";

      equation
          for i in 1:Ncell loop
            u[i]  = mdot/state[i].d/Aflow;
            Re[i] =Miscellanea.numbers.Reynolds(
                  u[i],
                  state[i].d,
                  state[i].eta,
                  d_s_eq);
            Pr[i] =Miscellanea.numbers.Prandtl(
                  state[i].cp,
                  state[i].eta,
                  state[i].lambda);
            Nu[i] = 0.6246*Re[i]^(-0.4989)*Re[i]*Pr[i]^(1/3);
            ht[i] =Miscellanea.numbers.Nusselt(
                  Nu[i],
                  state[i].lambda,
                  d_s_eq);
          end for;
      end single_phase_Kern;

      class condensation_Kern
        "Kern correlation for shell side for condensation outside tube bundles"
          extends Heat_transfer.Shell.Base_classes.base_Kern;
          parameter Boolean shear_vapor = true
          "= true, if shear vapor effect is considered";
          input Modelica.SIunits.Length l "Lenght";
          input Real N_tubes "Number of tubes in the bundle";
          input Modelica.SIunits.AbsolutePressure p_in "Inlet pressure";
          Modelica.SIunits.CoefficientOfHeatTransfer ht[Ncell]
          "Heat transfer coefficient";
          Modelica.SIunits.CoefficientOfHeatTransfer ht_l
          "Heat transfer coefficient if all liquid phase";
          Modelica.SIunits.NusseltNumber ht_g[Ncell]
          "Nusselt number gravity-controlled";
          Modelica.SIunits.NusseltNumber ht_sh[Ncell]
          "Nusselt number shear-controlled condensation";
          Modelica.SIunits.Length d_s_eq "Equivalent shell diameter";
          Modelica.SIunits.ReynoldsNumber Re[Ncell] "Reynolds number";
          Modelica.SIunits.PrandtlNumber Pr[Ncell] "Prandtl number";
          Modelica.SIunits.ReynoldsNumber Re_l
          "Reynolds number if all liquid phase";
          Modelica.SIunits.PrandtlNumber Pr_l
          "Prandtl number if all liquid phase";
          Modelica.SIunits.Velocity u[Ncell] "Velocity";
          Modelica.SIunits.Velocity u_l "Velocity if all liquid phase";
          Real X_tt[Ncell] "Martinelli parameter";
          Real xq[Ncell] "Vapor quality";
          Medium.SaturationProperties sat "Saturation properties";
          Medium.ThermodynamicState state_l
          "Thermodynamic state in saturated liquid";
          Medium.ThermodynamicState state_v
          "Thermodynamic state in saturated vapor";

      equation
          sat     = Medium.setSat_p(p_in);
          state_l = Medium.setState_ph(p_in, sat.hl);
          state_v = Medium.setState_ph(p_in, sat.hv);
          u_l     = mdot/state_l.d/Aflow;
          Re_l    = Miscellanea.numbers.Reynolds(u_l, state_l.d, state_l.eta, d_s_eq);
          Pr_l    = Miscellanea.numbers.Prandtl(state_l.cp, state_l.eta,
          state_l.lambda);
          ht_l  = 0.6246*Re_l^(-0.4989)*Re_l*Pr_l^(1/3)*state_l.lambda/d_s_eq;

          for i in 1:Ncell loop

            u[i]  = mdot/state[i].d/Aflow "tube velocity";
            Re[i] =Miscellanea.numbers.Reynolds(
                  u[i],
                  state[i].d,
                  state[i].eta,
                  d_s_eq);
            Pr[i] =Miscellanea.numbers.Prandtl(
                  state[i].cp,
                  state[i].eta,
                  state[i].lambda);

            if (state[i].h < state_l.h or state[i].h > state_v.h) then

              ht_g[i]  = 0.6246*Re[i]^(-0.4989)*Re[i]*Pr[i]^(1/3)*state[i].lambda/d_s_eq;
              xq[i]    = 1;
              X_tt[i]  = 0;
              ht_sh[i] = ht_g[i];

            else

              xq[i]    = (state[i].h - state_l.h)/(state_v.h - state_l.h);
              X_tt[i]  = (1/xq[i] - 1)^0.9*(state_v.d/state_l.d)^0.5*(state_l.eta
              /state_v.eta)^0.1;
              ht_sh[i] = 1.26*(1/X_tt[i])^0.78*ht_l;
              ht_g[i]  = 0.95*state_l.lambda*(g_n*l*N_tubes^(2/3)*state_l.d*(state_l.d
              - state_v.d)/(state_l.eta*mdot))^(1/3);

            end if;

            if shear_vapor then
              ht[i]    = sqrt(ht_g[i]^2 + ht_sh[i]^2);
            else
              ht[i]    = ht_g[i];
            end if;

          end for;

      end condensation_Kern;

      class single_phase_Bell_Delaware
        "Bell method for shell side only single phase to single phase"
          extends Heat_transfer.Shell.Base_classes.base_Bell_Delaware;
          Modelica.SIunits.CoefficientOfHeatTransfer ht[Ncell]
          "Heat transfer coefficient";
          Modelica.SIunits.ReynoldsNumber Re[Ncell] "Reynolds number";
          Modelica.SIunits.PrandtlNumber Pr[Ncell] "Prandtl number";
          Modelica.SIunits.Velocity u[Ncell] "Velocity";
      protected
          Real csi[Ncell] "Friction factor";
          Real a1[Ncell];
          Real a2[Ncell];
          Real a[Ncell];
      equation

        for i in 1:Ncell loop

          u[i]  = mdot/state[i].d/S_m;
          Re[i] =Miscellanea.numbers.Reynolds(
                  u[i],
                  state[i].d,
                  state[i].eta,
                  Dhyd_o);
          Pr[i] =Miscellanea.numbers.Prandtl(
                  state[i].cp,
                  state[i].eta,
                  state[i].lambda);

          if (Re[i] <= 1e1) then
            if (layout == 1) then
              a1[i]  = 1.4;
              a2[i]  = -0.667;
            else
              a1[i]  = 0.97;
              a2[i]  = -0.667;
            end if;
          elseif (Re[i] > 1e1 and Re[i] <= 1e2) then
            if (layout == 1) then
              a1[i]  = 1.36;
              a2[i]  = -0.657;
            else
              a1[i]  = 0.9;
              a2[i]  = -0.631;
            end if;
          elseif (Re[i] > 1e2 and Re[i] <= 1e3) then
            if (layout == 1) then
              a1[i]  = 0.593;
              a2[i]  = -0.477;
            else
              a1[i]  = 0.408;
              a2[i]  = -0.460;
            end if;
          elseif (Re[i] > 1e3 and Re[i] < 1e4) then
            if (layout == 1) then
              a1[i]  = 0.321;
              a2[i]  = -0.388;
            else
              a1[i]  = 0.107;
              a2[i]  = -0.266;
            end if;
          else
            if (layout == 1) then
              a1[i]  = 0.321;
              a2[i]  = -0.388;
            else
              a1[i]  = 0.370;
              a2[i]  = -0.395;
            end if;
          end if;

          if (layout == 1) then
            a[i]   = 1.450/(1 + 0.14*Re[i]^0.519);
          else
            a[i]   = 1.187/(1 + 0.14*Re[i]^0.370);
          end if;

          csi[i] = a1[i]*(1.33/pitch_f)^a[i]*Re[i]^a2[i];
          ht[i]  = J_l*J_b*J_c*csi[i]*state[i].cp*mdot/(S_m*Pr[i]^(2/3));

        end for;

      end single_phase_Bell_Delaware;

      class condensation_Bell_Delaware
        "Bell method for shell side for condensation outside tube bundles"
          extends Heat_transfer.Shell.Base_classes.base_Bell_Delaware;
          parameter Boolean shear_vapor = true
          "= true, if shear vapor effect is considered";
          input Modelica.SIunits.AbsolutePressure p_in "Inlet pressure";
          input Modelica.SIunits.Length l "Lenght";
          input Real N_tubes "Number of tubes in the bundle";
          Modelica.SIunits.CoefficientOfHeatTransfer ht[Ncell]
          "Heat transfer coefficient";
          Modelica.SIunits.CoefficientOfHeatTransfer ht_l
          "Heat transfer coefficient if all liquid phase";
          Modelica.SIunits.NusseltNumber ht_g[Ncell]
          "Nusselt number gravity-controlled";
          Modelica.SIunits.NusseltNumber ht_sh[Ncell](start = 1e4*ones(Ncell),fixed = true)
          "Nusselt number shear-controlled condensation";
          Modelica.SIunits.ReynoldsNumber Re[Ncell] "Reynolds number";
          Modelica.SIunits.PrandtlNumber Pr[Ncell] "Prandtl number";
          Modelica.SIunits.ReynoldsNumber Re_l
          "Reynolds number if all liquid phase";
          Modelica.SIunits.PrandtlNumber Pr_l
          "Prandtl number if all liquid phase";
          Modelica.SIunits.Velocity u[Ncell] "Velocity";
          Modelica.SIunits.Velocity u_l "Velocity if all liquid phase";
          Real X_tt[Ncell] "Martinelli parameter";
          Real xq[Ncell] "Vapor quality";
          Medium.SaturationProperties sat "Saturation properties";
          Medium.ThermodynamicState state_l
          "Thermodynamic state in saturated liquid";
          Medium.ThermodynamicState state_v
          "Thermodynamic state in saturated vapor";

      protected
          Real csi[Ncell] "Friction factor";
          Real a1[Ncell];
          Real a2[Ncell];
          Real a[Ncell];
          Real csi_l "Friction factor if all liquid phase";
          Real a1_l;
          Real a2_l;
          Real a_l;
      equation
        sat     = Medium.setSat_p(p_in);
        state_l = Medium.setState_ph(p_in, sat.hl);
        state_v = Medium.setState_ph(p_in, sat.hv);
        u_l     = mdot/state_l.d/S_m;
        Re_l    = Miscellanea.numbers.Reynolds(u_l, state_l.d, state_l.eta, Dhyd_o);
        Pr_l    = Miscellanea.numbers.Prandtl(state_l.cp, state_l.eta,
        state_l.lambda);
        if (Re_l <= 1e1) then
          if (layout == 1) then
            a1_l  = 1.4;
            a2_l  = -0.667;
          else
            a1_l  = 0.97;
            a2_l  = -0.667;
          end if;
        elseif (Re_l > 1e1 and Re_l <= 1e2) then
          if (layout == 1) then
            a1_l  = 1.36;
            a2_l  = -0.657;
          else
            a1_l  = 0.9;
            a2_l  = -0.631;
          end if;
        elseif (Re_l > 1e2 and Re_l <= 1e3) then
          if (layout == 1) then
            a1_l  = 0.593;
            a2_l  = -0.477;
          else
            a1_l  = 0.408;
            a2_l  = -0.460;
          end if;
        elseif (Re_l > 1e3 and Re_l < 1e4) then
          if (layout == 1) then
            a1_l  = 0.321;
            a2_l  = -0.388;
          else
            a1_l  = 0.107;
            a2_l  = -0.266;
          end if;
        else
          if (layout == 1) then
            a1_l  = 0.321;
            a2_l  = -0.388;
          else
            a1_l  = 0.370;
            a2_l  = -0.395;
          end if;
        end if;

        if (layout == 1) then
          a_l   = 1.450/(1 + 0.14*Re_l^0.519);
        else
          a_l   = 1.187/(1 + 0.14*Re_l^0.370);
        end if;

        csi_l = a1_l*(1.33/pitch_f)^a_l*Re_l^a2_l;
        ht_l  = J_l*J_b*J_c*csi_l*state_l.cp*mdot/(S_m*Pr_l^(2/3));

        for i in 1:Ncell loop

          u[i]  = mdot/state[i].d/S_m;
          Re[i] =Miscellanea.numbers.Reynolds(
                  u[i],
                  state[i].d,
                  state[i].eta,
                  Dhyd_o);
          Pr[i] =Miscellanea.numbers.Prandtl(
                  state[i].cp,
                  state[i].eta,
                  state[i].lambda);

          if (Re[i] <= 1e1) then
            if (layout == 1) then
              a1[i]  = 1.4;
              a2[i]  = -0.667;
            else
              a1[i]  = 0.97;
              a2[i]  = -0.667;
            end if;
          elseif (Re[i] > 1e1 and Re[i] <= 1e2) then
            if (layout == 1) then
              a1[i]  = 1.36;
              a2[i]  = -0.657;
            else
              a1[i]  = 0.9;
              a2[i]  = -0.631;
            end if;
          elseif (Re[i] > 1e2 and Re[i] <= 1e3) then
            if (layout == 1) then
              a1[i]  = 0.593;
              a2[i]  = -0.477;
            else
              a1[i]  = 0.408;
              a2[i]  = -0.460;
            end if;
          elseif (Re[i] > 1e3 and Re[i] < 1e4) then
            if (layout == 1) then
              a1[i]  = 0.321;
              a2[i]  = -0.388;
            else
              a1[i]  = 0.107;
              a2[i]  = -0.266;
            end if;
          else
            if (layout == 1) then
              a1[i]  = 0.321;
              a2[i]  = -0.388;
            else
              a1[i]  = 0.370;
              a2[i]  = -0.395;
            end if;
          end if;

          if (layout == 1) then
            a[i]   = 1.450/(1 + 0.14*Re[i]^0.519);
          else
            a[i]   = 1.187/(1 + 0.14*Re[i]^0.370);
          end if;

          csi[i] = a1[i]*(1.33/pitch_f)^a[i]*Re[i]^a2[i];

          if (state[i].h < state_l.h or state[i].h > state_v.h) then

            ht_g[i]  = J_l*J_b*J_c*csi[i]*state[i].cp*mdot/(S_m*Pr[i]^(2/3));
            xq[i]    = 1;
            X_tt[i]  = 0;
            ht_sh[i] = ht_g[i];

          else

            xq[i]    = (state[i].h - state_l.h)/(state_v.h - state_l.h);
            X_tt[i]  = (1/xq[i] - 1)^0.9*(state_v.d/state_l.d)^0.5*(state_l.eta
            /state_v.eta)^0.1;
            ht_sh[i] = 1.26*(1/X_tt[i])^0.78*ht_l;
            ht_g[i]  = 0.95*state_l.lambda*(g_n*l*N_tubes^(2/3)*state_l.d*(state_l.d
              - state_v.d)/(state_l.eta*mdot))^(1/3);

          end if;

          if shear_vapor then
            ht[i]    = sqrt(ht_g[i]^2 + ht_sh[i]^2);
          else
            ht[i]    = ht_g[i];
          end if;

        end for;

      end condensation_Bell_Delaware;

      package Base_classes "Base classes for the heat transfer"
        class base_Kern "Kern correlation for shell side base class"
            extends Heat_transfer.Shell.Base_classes.base_ht;
            input Modelica.SIunits.Length Dhyd_o "Outer hydraulic diameter";
            parameter Real pitch_f "Tube pitch";
            parameter Integer layout = 1
            "Tube layout 1 = triangular, 2 = squared";
            Modelica.SIunits.Length d_s_eq "Equivalent shell diameter";

        equation
            if layout == 1 then
              d_s_eq  = 1.1*Dhyd_o*(pitch_f^2 - 0.917);
            elseif layout == 2 then
              d_s_eq  = 1.27*Dhyd_o*(pitch_f^2 - 0.785);
            else
              d_s_eq  = 1;
            end if;
        end base_Kern;

        class base_Bell_Delaware "Bell Delaware correlation for shell side"
            extends Heat_transfer.Shell.Base_classes.base_ht;
            input Modelica.SIunits.Length Dhyd "Hydraulic diameter";
            input Modelica.SIunits.Length Dhyd_o "Outer hydraulic diameter";
            parameter Real b_cut "Baffle cut";
            parameter Real ttb = 8e-4 "Tube to baffle clearance";
            parameter Real bts = 4.8e-3 "Baffle to shell clearance";
            parameter Real pitch_f "Tube pitch";
            parameter Real N_ss
            "The number of sealing strips (pairs) in one baffle spacing";
            parameter Integer layout "Tube layout, 1 = triangular, 2 = squared";
            parameter Integer N_passes "Number of tube passes";
            parameter Real ff = 1 "Fouling factor (=1 if clean)";
            input Modelica.SIunits.Length d_s "Shell diameter";
            input Modelica.SIunits.Length d_b "Bundle diameter";
            input Real N_tubes "Number of tubes in the bundle";
            input Modelica.SIunits.Length l_b "Baffle length";
            Modelica.SIunits.Length tpv "Vertical tube pitch";
            Real N_tcc "Number of tube rows between the baffle tips";
            Real N_tcw
            "The effective number of tube rows crossed in the baffle window";
            Real F_sbp
            "The ratio of the bypass area to the overall crossflow area";
            Real F_w "Fraction of number of tubes in the baffle window";
            Real F_c "Fraction of number of tubes in pure crossflow";
            Real J_c "Segmental baffle window correction factor";
            Real J_l
            "Correction factors for baffle leakage effects for heat transfer";
            Real J_b
            "Correction factors for bundle bypass effects for heat transfer";
            Modelica.SIunits.Area S_w "Window area";
            Modelica.SIunits.Area S_tb "The tube-to-baffle leakage area";
            Modelica.SIunits.Area S_sb "The shell-to-baffle leakage area";
            Modelica.SIunits.Area S_m "Reference normal area";
            Modelica.SIunits.Area S_wg "The gross window flow area";
            Modelica.SIunits.Area S_wt
            "The segmental baffle window area occupied by the tubes";
            Modelica.SIunits.Area S_b
            "The bypass area between the shell and the tube bundle within one baffle";
            Real r_s "Correlational parameter";
            Real r_lm "Correlational parameter";
            Modelica.SIunits.Angle teta_ds "The centriangle of baffle cut";
            Modelica.SIunits.Angle teta_ctl
            "The upper centriangle of baffle cut";

        equation
          if (layout == 1) then
            tpv = 0.866*pitch_f*Dhyd_o;
          else
            tpv = pitch_f*Dhyd_o;
          end if;

          teta_ctl = 2*acos(d_s/(d_b - Dhyd_o)*(1 - 2*b_cut));
          teta_ds  = 2*acos(1 - 2*b_cut);
          N_tcc    = d_s*(1 - 2*b_cut)/tpv;
          S_m      = l_b*(d_s - d_b + (d_b - Dhyd_o)*(1 - 1/pitch_f));
          S_wg     = 0.25*pi*d_s^2*(teta_ds/(2*pi) - sin(teta_ds)/(2*pi));
          F_w      = teta_ctl/(2*pi) - sin(teta_ctl)/(2*pi);
          S_wt     = 0.25*N_tubes*F_w*pi*Dhyd_o^2;

          S_tb     = 0.25*pi*((ttb + Dhyd_o)^2 - Dhyd_o^2)*N_tubes*(1 - F_w);
          S_sb     = 0.5*bts*d_s*(1 - teta_ds/(2*pi));
          r_s      = S_sb/(S_tb + S_sb);
          r_lm     = (S_tb + S_sb)/S_m;
          J_l      = 0.44*(1 - r_s) + (1 - 0.44*(1 - r_s))*exp(-2.2*r_lm);

          S_b      = l_b*(d_s - d_b);
          F_sbp    = S_b/S_m;
          J_b      = exp(-1.35*F_sbp*(1 - (2/N_ss)^(1/3)));

          F_c      = 1 - 2*F_w;
          J_c      = F_c + 0.54*(1 - F_c)^0.345;
          S_w      = S_wg - S_wt;
          N_tcw    = 0.8*(d_s*b_cut - 0.5*(d_s - (d_b - Dhyd_o)))/tpv;

        end base_Bell_Delaware;

        class base_ht "Basic heat transfer correlation"
            replaceable package Medium = Media.OneRandomOrganicFluid constrainedby
            Modelica.Media.Interfaces.PartialMedium "Medium model" annotation(choicesAllMatching = true);
            parameter Integer Ncell(start=3) "Number of cell elements";
            input Modelica.SIunits.Area Aflow
            "Cross-sectional area (single tube)";
            input Modelica.SIunits.MassFlowRate mdot "Mass flow rate";
            input Medium.ThermodynamicState state[Ncell];
        end base_ht;
      end Base_classes;
      annotation (Icon(graphics={
            Rectangle(
              extent={{-80,40},{80,-40}},
              lineColor={0,0,0}),
            Rectangle(
              extent={{-70,-40},{-50,-60}},
              lineColor={0,0,0}),
            Rectangle(
              extent={{50,60},{70,40}},
              lineColor={0,0,0}),
            Line(
              points={{60,70},{60,46}},
              color={255,0,0},
              smooth=Smooth.None),
            Line(
              points={{2,-2},{-2,2}},
              color={255,0,0},
              smooth=Smooth.None,
              origin={62,48},
              rotation=90),
            Line(
              points={{2,2},{-2,-2}},
              color={255,0,0},
              smooth=Smooth.None,
              origin={58,48},
              rotation=90),
            Line(
              points={{-60,-42},{-60,-66}},
              color={255,0,0},
              smooth=Smooth.None),
            Line(
              points={{2,-2},{-2,2}},
              color={255,0,0},
              smooth=Smooth.None,
              origin={-58,-64},
              rotation=90),
            Line(
              points={{2,2},{-2,-2}},
              color={255,0,0},
              smooth=Smooth.None,
              origin={-62,-64},
              rotation=90)}));
    end Shell;

    package Plates
      class single_phase_Martin
        "Martin heat transfer correlation for single phase heat transfer"
          extends Base_classes.base_ht;
          parameter Modelica.SIunits.Angle beta "Plate inclination angle";
          Modelica.SIunits.ReynoldsNumber Re[N_ch, N_cell_pc] "Reynolds number";
          Modelica.SIunits.PrandtlNumber Pr[N_ch, N_cell_pc] "Prandtl number";
          Modelica.SIunits.Velocity u[N_ch, N_cell_pc] "Velocity";
          Modelica.SIunits.NusseltNumber Nu[N_ch, N_cell_pc] "Nusselt number";
          Modelica.SIunits.CoefficientOfHeatTransfer ht[N_ch, N_cell_pc]
          "Heat transfer coefficient";
      protected
          Real csi[N_ch, N_cell_pc] "Friction factor";
          Real csi0[N_ch, N_cell_pc] "Friction factor 0";
          Real csi1[N_ch, N_cell_pc] "Friction factor 1";
      equation
          for j in 1:N_ch loop
            for i in 1:N_cell_pc loop
              u[j, i]   = mdot/state[j, i].d/Aflow;
              Re[j, i]  =Miscellanea.numbers.Reynolds(
                    u[j, i],
                    state[j, i].d,
                    state[j, i].eta,
                    Dhyd);
              Pr[j, i]  =Miscellanea.numbers.Prandtl(
                    state[j, i].cp,
                    state[j, i].eta,
                    state[j, i].lambda);
              if (Re[j, i] < 2e3) then
                csi0[j, i] = 16/Re[j, i];
                csi1[j, i] = 149.25/Re[j, i] + 0.9625;
              else
                csi0[j, i] = 1/(1.56*log(Re[j, i]) - 3)^2;
                csi1[j, i] = 9.75/Re[j, i]^0.289;
              end if;
              csi[j, i]       = 1/((1 - cos(beta))/sqrt(3.8*csi1[j, i]) + cos(beta)/
              sqrt(0.045*tan(beta) + 0.09*sin(beta) + csi0[j, i]/cos(beta)))^2;
              Nu[j, i]  = 0.205*Pr[j, i]^(1/3)*(csi[j, i]*Re[j, i]^2*sin(2*beta))
              ^0.374;
              ht[j, i]  =Miscellanea.numbers.Nusselt(
                    Nu[j, i],
                    state[j, i].lambda,
                    Dhyd);
            end for;
          end for;
      end single_phase_Martin;

      class evaporation_Coulson
        "Cooper heat transfer correlation for evaporation heat transfer. Coulson for single phase."
          extends Base_classes.base_ht;
          parameter Medium.ThermodynamicState state_l
          "Thermodynamic state in saturated liquid";
          parameter Medium.ThermodynamicState state_v
          "Thermodynamic state in saturated vapor";
          parameter Modelica.SIunits.Length R_p = 1e-6
          "Relative roughness of the surface";
          input Modelica.SIunits.Area At "Width";
          parameter Modelica.SIunits.MolarMass M = Medium.fluidConstants[1].molarMass
          "Molar mass of the fluid";
          parameter Real pc = Medium.fluidConstants[1].criticalPressure
          "Critical pressure";
          parameter Modelica.SIunits.HeatFlux qdot_tilde_start
          "Heat flux start";
          input Modelica.SIunits.AbsolutePressure p_in "Inlet pressure";
          input Modelica.SIunits.Length l "Length";
          input Modelica.SIunits.HeatFlowRate qdot[N_ch, N_cell_pc] "Heat rate";
          Modelica.SIunits.ReynoldsNumber Re[N_ch, N_cell_pc] "Reynolds number";
          Modelica.SIunits.PrandtlNumber Pr[N_ch, N_cell_pc] "Prandtl number";
          Modelica.SIunits.Velocity u[N_ch, N_cell_pc] "Velocity";
          Modelica.SIunits.CoefficientOfHeatTransfer ht[N_ch, N_cell_pc]
          "Heat transfer coefficient";
          Modelica.SIunits.HeatFlux qdot_tilde[N_ch, N_cell_pc](
           start = fill(qdot_tilde_start, N_ch, N_cell_pc)) "Heat flux";

      equation
          for j in 1:N_ch loop
            for i in 1:N_cell_pc loop
              qdot_tilde[j, i] = qdot[j, i]/(2*At);
              u[j, i]          = mdot/state[j, i].d/Aflow;
              Re[j, i]         =Miscellanea.numbers.Reynolds(
                    u[j, i],
                    state[j, i].d,
                    state[j, i].eta,
                    Dhyd);
              Pr[j, i]         =Miscellanea.numbers.Prandtl(
                    state[j, i].cp,
                    state[j, i].eta,
                    state[j, i].lambda);
              if (state[j, i].h < state_l.h or state[j, i].h > state_v.h) then
                ht[j, i]       = 0.26*Pr[j, i]^0.4*Re[j, i]^0.65*state[j, i].lambda/Dhyd;
              else
                ht[j, i]      = 55*(p_in/pc)^(0.12 - 0.2*log10(R_p))/(-log10(
                p_in/pc))^0.55*abs(qdot_tilde[j, i])^0.67/sqrt(1e3*M);
              end if;
            end for;
          end for;
      end evaporation_Coulson;

      class evaporation_Martin
        "Cooper heat transfer correlation for evaporation heat transfer. Martin for single phase."
          extends Base_classes.base_ht;
          parameter Modelica.SIunits.Angle beta "Plate inclination angle";
          parameter Modelica.SIunits.Length R_p = 1e-6
          "Relative surface roughness";
          input Modelica.SIunits.Area At "Plate cell area";
          parameter Modelica.SIunits.MolarMass M = Medium.fluidConstants[1].molarMass
          "Molar mass of the fluid";
          parameter Real pc = Medium.fluidConstants[1].criticalPressure
          "Critical pressure";
          parameter Modelica.SIunits.HeatFlux qdot_tilde_start
          "Heat flux start";
          input Modelica.SIunits.AbsolutePressure p_in "Inlet pressure";
          input Modelica.SIunits.HeatFlowRate qdot[N_ch, N_cell_pc] "Heat rate";
          Modelica.SIunits.ReynoldsNumber Re[N_ch, N_cell_pc] "Reynolds number";
          Modelica.SIunits.PrandtlNumber Pr[N_ch, N_cell_pc] "Prandtl number";
          Modelica.SIunits.Velocity u[N_ch, N_cell_pc] "Velocity";
          Modelica.SIunits.CoefficientOfHeatTransfer ht[N_ch, N_cell_pc]
          "Heat transfer coefficient";
          Modelica.SIunits.HeatFlux qdot_tilde[N_ch, N_cell_pc](
            start = fill(qdot_tilde_start, N_ch, N_cell_pc)) "Heat flux";
          Medium.SaturationProperties sat "Saturation properties";
          Medium.ThermodynamicState state_l
          "Thermodynamic state in saturated liquid";
          Medium.ThermodynamicState state_v
          "Thermodynamic state in saturated vapor";
      protected
           Real csi[N_ch, N_cell_pc] "Friction factor";
           Real csi0[N_ch, N_cell_pc] "Friction factor 0";
           Real csi1[N_ch, N_cell_pc] "Friction factor 1";
      equation
          sat     = Medium.setSat_p(p_in);
          state_l = Medium.setState_ph(p_in, sat.hl);
          state_v = Medium.setState_ph(p_in, sat.hv);

          for j in 1:N_ch loop
            for i in 1:N_cell_pc loop
              qdot_tilde[j, i] = qdot[j, i]/(2*At);
              u[j, i]          = mdot/state[j, i].d/Aflow;
              Re[j, i]         = Miscellanea.numbers.Reynolds(
                    u[j, i],
                    state[j, i].d,
                    state[j, i].eta,
                    Dhyd);
              Pr[j, i]         = Miscellanea.numbers.Prandtl(
                    state[j, i].cp,
                    state[j, i].eta,
                    state[j, i].lambda);
              if (Re[j, i] < 2e3) then
                csi0[j, i] = 16/Re[j, i];
                csi1[j, i] = 149.25/Re[j, i] + 0.9625;
              else
                csi0[j, i] = 1/(1.56*log(Re[j, i]) - 3)^2;
                csi1[j, i] = 9.75/Re[j, i]^0.289;
              end if;
              csi[j, i]       = 1/((1 - cos(beta))/sqrt(3.8*csi1[j, i]) + cos(beta)/
              sqrt(0.045*tan(beta) + 0.09*sin(beta) + csi0[j, i]/cos(beta)))^2;
              if (state[j, i].h < state_l.h or state[j, i].h > state_v.h) then
                ht[j, i]      = 0.205*Pr[j, i]^(1/3)*(csi[j, i]*Re[j, i]^2*sin(2*beta))
                ^0.374*state[j, i].lambda/Dhyd;
              else
                 ht[j, i]      = 55*(p_in/pc)^(0.12 - 0.2*log10(R_p))/(-log10(
                 p_in/pc))^0.55*abs(qdot_tilde[j, i])^0.67/sqrt(1e3*M);
              end if;
            end for;
          end for;
      end evaporation_Martin;

      class single_phase_Coulson
        "Coulson heat transfer correlation for single phase heat transfer"
          extends Base_classes.base_ht;
          Modelica.SIunits.ReynoldsNumber Re[N_ch, N_cell_pc] "Reynolds number";
          Modelica.SIunits.PrandtlNumber Pr[N_ch, N_cell_pc] "Prandtl number";
          Modelica.SIunits.Velocity u[N_ch, N_cell_pc] "Velocity";
          Modelica.SIunits.NusseltNumber Nu[N_ch, N_cell_pc] "Nusselt number";
          Modelica.SIunits.CoefficientOfHeatTransfer ht[N_ch, N_cell_pc]
          "Heat transfer coefficient";

      equation
          for j in 1:N_ch loop
            for i in 1:N_cell_pc loop
              u[j, i]   = mdot/state[j, i].d/Aflow;
              Re[j, i]  =Miscellanea.numbers.Reynolds(
                    u[j, i],
                    state[j, i].d,
                    state[j, i].eta,
                    Dhyd);
              Pr[j, i]  =Miscellanea.numbers.Prandtl(
                    state[j, i].cp,
                    state[j, i].eta,
                    state[j, i].lambda);
              Nu[j, i]  = 0.26*Pr[j, i]^0.4*Re[j, i]^0.65;
              ht[j, i]  =Miscellanea.numbers.Nusselt(
                    Nu[j, i],
                    state[j, i].lambda,
                    Dhyd);
            end for;
          end for;
      end single_phase_Coulson;

      package Base_classes "Base classes for the heat transfer"
        class base_ht "Basic heat transfer correlation"
            replaceable package Medium = Media.OneRandomOrganicFluid constrainedby
            Modelica.Media.Interfaces.PartialMedium "Medium model" annotation(choicesAllMatching = true);
            parameter Integer N_ch(start = 3) "Total number of channels";
            parameter Integer N_cell_pc( start = 3)
            "Number of cells per channel";
            parameter Modelica.SIunits.Length Dhyd "Hydraulic diameter";
            input Modelica.SIunits.Area Aflow "Cross-sectional area";
            input Modelica.SIunits.MassFlowRate mdot "Mass flow rate";
            input Medium.ThermodynamicState state[N_ch, N_cell_pc]
            "Thermodynamic states";
        end base_ht;
      end Base_classes;
      annotation (Icon(graphics={      Polygon(
              points={{-30,-80},{30,-60},{30,80},{-30,60},{-30,-80}},
              lineColor={0,0,0},
              smooth=Smooth.None,
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid)}));
    end Plates;
  annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,
            -100},{100,100}}), graphics={
        Rectangle(
          extent={{60,88},{80,-88}},
          lineColor={95,95,95},
          fillColor={215,215,215},
          fillPattern=FillPattern.Solid),
        Polygon(
          points={{6,58},{6,80},{46,60},{6,40},{6,58}},
          lineColor={255,0,0},
          smooth=Smooth.None,
          fillColor={255,0,0},
          fillPattern=FillPattern.Solid),
        Rectangle(
          extent={{-94,66},{6,54}},
          lineColor={255,0,0},
          fillColor={255,0,0},
          fillPattern=FillPattern.Solid),
        Polygon(
          points={{6,-42},{6,-20},{46,-40},{6,-60},{6,-42}},
          lineColor={255,0,0},
          smooth=Smooth.None,
          fillColor={255,0,0},
          fillPattern=FillPattern.Solid),
        Rectangle(
          extent={{-94,-34},{6,-46}},
          lineColor={255,0,0},
          fillColor={255,0,0},
          fillPattern=FillPattern.Solid)}));
  end Heat_transfer;

  package Media "Fluid library"
    package OneRandomOrganicFluid
      "Change the name to something more appropriate"
        extends ExternalMedia.Media.FluidPropMedium(
        mediumName = "Name of the fluid for documentation purposes",
        libraryName = "FluidProp.RefProp",
        substanceNames = {"cyclopentane"},
        ThermoStates = Modelica.Media.Interfaces.PartialMedium.Choices.IndependentVariables.ph);
    end OneRandomOrganicFluid;

  annotation (Icon(graphics={
          Rectangle(
            lineColor={200,200,200},
            fillColor={248,248,248},
            fillPattern=FillPattern.HorizontalCylinder,
            extent={{-100,-100},{100,100}},
            radius=25.0),
          Rectangle(
            lineColor={128,128,128},
            fillPattern=FillPattern.None,
            extent={{-100,-100},{100,100}},
            radius=25.0),
          Line(
            points={{-76,-80},{-62,-30},{-32,40},{4,66},{48,66},{73,45},{62,-8},
              {48,-50},{38,-80}},
            color={64,64,64},
            smooth=Smooth.Bezier),
          Line(
            points={{-40,20},{68,20}},
            color={175,175,175},
            smooth=Smooth.None),
          Line(
            points={{-40,20},{-44,88},{-44,88}},
            color={175,175,175},
            smooth=Smooth.None),
          Line(
            points={{68,20},{86,-58}},
            color={175,175,175},
            smooth=Smooth.None),
          Line(
            points={{-60,-28},{56,-28}},
            color={175,175,175},
            smooth=Smooth.None),
          Line(
            points={{-60,-28},{-74,84},{-74,84}},
            color={175,175,175},
            smooth=Smooth.None),
          Line(
            points={{56,-28},{70,-80}},
            color={175,175,175},
            smooth=Smooth.None),
          Line(
            points={{-76,-80},{38,-80}},
            color={175,175,175},
            smooth=Smooth.None),
          Line(
            points={{-76,-80},{-94,-16},{-94,-16}},
            color={175,175,175},
            smooth=Smooth.None)}));
  end Media;
  annotation (uses(Modelica(version="3.2.1")));
end Design;
