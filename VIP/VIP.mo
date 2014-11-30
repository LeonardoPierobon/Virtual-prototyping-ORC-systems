within ;
package VIP "I am a package for the Virtual Prototyping Environment"
  package Heat_transfer "A package containing heat transfer correlations"
    package Tubes "heat transfer correlations in tubes"
      class Dittus_Boelter "Dittus Boelter correlation for tubes"
          replaceable package Medium = VIP.Media.OneRandomOrganicFluid
          "Medium model";
          parameter Integer Ncell(start=3) "Number of cell elements";
          parameter Modelica.SIunits.Length l "Lenght (single tube)";
          output Modelica.SIunits.NusseltNumber Nu[Ncell]
          "Nusselt number tubes ";
          output Modelica.SIunits.CoefficientOfHeatTransfer ht[Ncell]
          "Heat transfer coefficient tubes";
          input Modelica.SIunits.Length Dhyd "Hydraulic Diameter (single tube)";
          input Medium.ThermodynamicState state[Ncell];
          input Modelica.SIunits.ReynoldsNumber Re[Ncell]
          "Reynolds number tubes (average)";
          input Real alfa "exponent for the Prandtl number";
          input Modelica.SIunits.PrandtlNumber Pr[Ncell] "Prandtl number tubes";
      equation
            for i in 1:Ncell loop
              assert(Re[i] > 1e4, "Reynolds number is lower than 1e4 to use Dittus and Boelter", AssertionLevel.warning);
              Nu[i]             =  2.3e-2*Re[i]^0.8*Pr[i]^alfa
            "Nusselt number tubes";
              ht[i]             = Useful_functions.Pure_numbers.Nusselt(
                                                       Nu[i], state[i].lambda, Dhyd)
            "Heat transfer coefficient tube side";
            end for;

      end Dittus_Boelter;

      class Sieder_Tate "Sieder Tate correlation for tubes"
          replaceable package Medium = VIP.Media.OneRandomOrganicFluid
          "Medium model";
          parameter Integer Ncell(start=3) "Number of cell elements";
          parameter Modelica.SIunits.Length l "Lenght (single tube)";
          output Modelica.SIunits.NusseltNumber Nu[Ncell]
          "Nusselt number tubes ";
          output Modelica.SIunits.CoefficientOfHeatTransfer ht[Ncell]
          "Heat transfer coefficient tubes";
          input Modelica.SIunits.Length Dhyd "Hydraulic Diameter (single tube)";
          input Medium.ThermodynamicState state[Ncell];
          input Modelica.SIunits.ReynoldsNumber Re[Ncell]
          "Reynolds number tubes (average)";
          input Modelica.SIunits.DynamicViscosity eta_wall[Ncell]
          "exponent for the Prandtl number";
          input Modelica.SIunits.PrandtlNumber Pr[Ncell] "Prandtl number tubes";
      equation
            for i in 1:Ncell loop
              assert(Re[i] > 1e4, "Reynolds number is lower than 1e4 to use Sieder and Tate", AssertionLevel.warning);
              assert(Pr[i] > 0.6,  "Prandtl number is lower than 0.6 to be used in Sieder and Tate", AssertionLevel.warning);
              Nu[i]             = 2.3e-2*Re[i]^0.8*Pr[i]^(1.0/3)*(eta_wall[i]/state[i].eta)^0.14
            "Nusselt number tubes";
              ht[i]             = Useful_functions.Pure_numbers.Nusselt(
                                                       Nu[i], state[i].lambda, Dhyd)
            "Heat transfer coefficient tube side";
            end for;

      end Sieder_Tate;

      class Gnielinski "Gnielinski correlation for tubes"
          replaceable package Medium = VIP.Media.OneRandomOrganicFluid
          "Medium model";
          parameter Integer Ncell(start=3) "Number of cell elements";
          parameter Modelica.SIunits.Length l "Lenght (single tube)";
          output Modelica.SIunits.NusseltNumber Nu[Ncell]
          "Nusselt number tubes ";
          output Modelica.SIunits.CoefficientOfHeatTransfer ht[Ncell]
          "Heat transfer coefficient tubes";
          input Modelica.SIunits.Length Dhyd "Hydraulic Diameter (single tube)";
          input Medium.ThermodynamicState state[Ncell];
          input Modelica.SIunits.ReynoldsNumber Re[Ncell]
          "Reynolds number tubes (average)";
          input Modelica.SIunits.PrandtlNumber Pr[Ncell] "Prandtl number tubes";
          output Real csi[Ncell] "Friction factor";
      equation
            for i in 1:Ncell loop
              csi[i]            = 1/(0.78*log(Re[i]) - 1.5)^2 "Friction factor";
              Nu[i]             =  ((csi[i]/8)*Re[i]*Pr[i])/(1 + 12.7*sqrt(csi[i]/8)*(Pr[i]^(2.0/3) - 1))*(1 +  (Dhyd/l)^(2.0/3))
            "Nusselt number tubes";
              ht[i]             = Useful_functions.Pure_numbers.Nusselt(
                                                       Nu[i], state[i].lambda, Dhyd)
            "Heat transfer coefficient tube side";
            end for;

      end Gnielinski;

      class EagleFerguson
        "Eagle-Ferguson heat transfer correlation (only for liquid water)"
          replaceable package Medium = VIP.Media.OneRandomOrganicFluid
          "Medium model";
          parameter Integer Ncell(start=3) "Number of cell elements";
          parameter Modelica.SIunits.Length l "Lenght (single tube)";
          output Modelica.SIunits.NusseltNumber Nu[Ncell]
          "Nusselt number tubes ";
          input Modelica.SIunits.Length Dhyd "Hydraulic Diameter (single tube)";
          input Modelica.SIunits.Velocity u[Ncell](start=ones(Ncell))
          "Velocity (average)";
          Modelica.SIunits.CoefficientOfHeatTransfer ht[Ncell]
          "Heat transfer coefficient tubes";
          Medium.ThermodynamicState state[Ncell];

      equation
            for i in 1:Ncell loop
              Nu[i]             =  4.2e3*Dhyd/state[i].lambda*(1.35 + 2e-2*(state[i].T - 273.15))*u[i]^0.8/(1e3*Dhyd)^0.2
            "Nusselt number tubes";
              ht[i]             = Useful_functions.Pure_numbers.Nusselt(
                                                       Nu[i], state[i].lambda, Dhyd)
            "Heat transfer coefficient tube side";
            end for;
      end EagleFerguson;

      class basic_ht "Basic heat transfer correlation"
          replaceable package Medium = VIP.Media.OneRandomOrganicFluid
          "Medium model";
          parameter Integer Ncell(start=3) "Number of cell elements";
          parameter Modelica.SIunits.Length l "Lenght (single tube)";
          input Modelica.SIunits.Length Dhyd "Hydraulic Diameter (single tube)";
          input Medium.ThermodynamicState state[Ncell];
      end basic_ht;
    end Tubes;

    package Shell "heat transfer correlations in shells"

      class Kern "Kern correlation for shell side"
          replaceable package Medium = VIP.Media.OneRandomOrganicFluid
          "Medium model";
          parameter Integer Ncell(start=3) "Number of cell elements";
          parameter Modelica.SIunits.Length l "Lenght (single tube)";
          output Modelica.SIunits.NusseltNumber Nu[Ncell]
          "Nusselt number tubes ";
          output Modelica.SIunits.CoefficientOfHeatTransfer ht[Ncell]
          "Heat transfer coefficient tubes";
          input Modelica.SIunits.Length Dhyd "Hydraulic Diameter (single tube)";
          input Medium.ThermodynamicState state[Ncell];
          input Modelica.SIunits.ReynoldsNumber Re[Ncell](start=10e6*ones(Ncell))
          "Reynolds number tubes (average)";
          input Modelica.SIunits.PrandtlNumber Pr[Ncell] "Prandtl number tubes";
      protected
          parameter Real x[2] = {1e5, 1e1};
          parameter Real y[2] = {2e-3, 1.98e-1};
      equation
            for i in 1:Ncell loop
              Nu[i]             = 10^log10(2e-3/1e5^(log10(y[2]/y[1])/log10(x[2]/x[1])))*Re[i]^(log10(y[2]/y[1])/log10(x[2]/x[1]))*Re[i]*Pr[i]^(1/3)
            "Nusselt number tubes";
              ht[i]             = Useful_functions.Pure_numbers.Nusselt(
                                                       Nu[i], state[i].lambda, Dhyd)
            "Heat transfer coefficient tube side";
            end for;
      end Kern;
    end Shell;
  end Heat_transfer;

  package Pressure_drops "A package containing pressure drops correlations"
    package Tubes "heat transfer correlations in tubes"
      class Frank "Frank correlation for tubes"
          replaceable package Medium = VIP.Media.OneRandomOrganicFluid
          "Medium model";
          parameter Integer Ncell(start=3) "Number of cell elements";
          parameter Modelica.SIunits.Length l "Lenght (single tube)";
          output Modelica.SIunits.AbsolutePressure dp[Ncell]
          "Pressure drops tubes";
          output Modelica.SIunits.AbsolutePressure dp_tot
          "Pressure drops tubes";
          input Modelica.SIunits.Length Dhyd "Hydraulic Diameter (single tube)";
          input Medium.ThermodynamicState state[Ncell];
          input Modelica.SIunits.ReynoldsNumber Re[Ncell]
          "Reynolds number tubes (average)";
          input Modelica.SIunits.Velocity u[Ncell](start=ones(Ncell))
          "Velocity (average)";
          input Modelica.SIunits.AbsolutePressure heads(start=2.5)
          "Number of velocity heads";
          output Real csi[Ncell] "Friction factor";

      equation
        for i in 1:Ncell loop
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

      class basic_dp "Basic pressure drop correlation"
          replaceable package Medium = VIP.Media.OneRandomOrganicFluid
          "Medium model";
          parameter Integer Ncell(start=3) "Number of cell elements";
          parameter Modelica.SIunits.Length l "Lenght (single tube)";
          input Modelica.SIunits.Length Dhyd "Hydraulic Diameter (single tube)";
          input Medium.ThermodynamicState state[Ncell];
      end basic_dp;
    end Tubes;

    package Shell "heat transfer correlations in shells"

      class Kern "Kern correlation for shell"
          replaceable package Medium = VIP.Media.OneRandomOrganicFluid
          "Medium model";
          parameter Integer Ncell(start=3) "Number of cell elements";
          parameter Modelica.SIunits.Length l "Lenght (single tube)";
          output Modelica.SIunits.AbsolutePressure dp[Ncell]
          "Pressure drops tubes";
          output Modelica.SIunits.AbsolutePressure dp_tot
          "Pressure drops tubes";
          input Modelica.SIunits.Length Dhyd "Hydraulic Diameter (single tube)";
          input Medium.ThermodynamicState state[Ncell];
          input Modelica.SIunits.ReynoldsNumber Re[Ncell](start=10e6*ones(Ncell))
          "Reynolds number tubes (average)";
          input Modelica.SIunits.Velocity u[Ncell](start=ones(Ncell))
          "Velocity (average)";
          input Modelica.SIunits.Length d_s "Shell diameter";
          input Modelica.SIunits.Length l_b "Baffle length";
          output Real csi[Ncell] "Friction factor";
      protected
          parameter Real x[2] = {3e2, 10};
          parameter Real y[2] = {1e-1, 2.4};
      equation

        for i in 1:Ncell loop

          if (Re[i] < 3e2) then
            csi[i] = 10^log10(2.5e-1/1e2^(log10(2.4/1e-1)/log10(10/3e2))) * Re[i]^(log10(2.4/1e-1)/log10(10/3e2));
          else
            csi[i] = 10^log10(4e-2/4e4^(log10(1e-1/2.4e-2)/log10(3e2/1e6))) * Re[i]^(log10(1e-1/2.4e-2)/log10(3e2/1e6));
          end if;
          dp[i]             =  4*csi[i]*d_s/Dhyd*l/l_b*state[i].d*u[i]^2
            "Local pressure drops";
        end for;

        dp_tot = sum(dp) "Total pressure drops";

      end Kern;
    end Shell;
  end Pressure_drops;

  package Objects "Package containing all the objects of the VIP"

    class tube_bundle
      "I am a tube bundle and I contain all the my relevant informations"
        extends VIP.Objects.tubes;
        replaceable function bundle_clearance =
          VIP.Useful_functions.Shell.Fixed_Utube;
        parameter Integer N_passes "Number of tube passes";
        parameter Integer layout "Tube layout 1 = triangular, 2 = squared";
        parameter Real pitch_f = 1.25
        "Tube pitch as a fraction of the outer tube diameter";
        input Real N_tubes "Number of tubes in the bundle";
        Modelica.SIunits.Length d_b "Bundle Diameter";
        Modelica.SIunits.Length clearance(start=68e-3)
        "Bundle - shell clearance";
    equation
        d_b       = Useful_functions.Shell.bundle_diameter(N_tubes, N_passes, Dhyd_o, layout)
        "calculate the tube bundle";
        clearance = bundle_clearance(d_b) "calculate the bundle clearance";

    end tube_bundle;

    class tubes "I am a tube and I contain all the my relevant informations"
        import constants = Modelica.Constants "to import pi-greco";
        replaceable package Medium = VIP.Media.OneRandomOrganicFluid
        "Medium model";
        replaceable VIP.Heat_transfer.Tubes.basic_ht heatTransfer constrainedby
        VIP.Heat_transfer.Tubes.basic_ht(Medium = Medium, Ncell = Ncell, l=l, state=state, Dhyd = Dhyd);
        replaceable VIP.Pressure_drops.Tubes.basic_dp pressureDrop constrainedby
        VIP.Pressure_drops.Tubes.basic_dp(Medium = Medium, Ncell = Ncell, l=l, state=state, Dhyd = Dhyd);
        parameter Integer Ncell(start=3) "Number of cell elements";
        parameter Modelica.SIunits.Length thick "Thickness of the tube";
        parameter Modelica.SIunits.Length l "Lenght (single tube)";
        parameter Modelica.SIunits.ThermalConductivity lambda_wall = 50
        "Thermal conductivity of the tube wall";
        input Modelica.SIunits.Area Aflow "Cross-sectional area (single tube)";
        input Modelica.SIunits.Length Dhyd "Hydraulic Diameter (single tube)";
        Modelica.SIunits.Area At "Heat transfer area associated to the tubes";
        Modelica.SIunits.Length Dhyd_o "Outer Diameter (single tube)";
        Modelica.SIunits.ThermalConductivity G_wall;
        Modelica.SIunits.AbsolutePressure p_in "Inlet pressure tube side";
        Modelica.SIunits.Velocity u[Ncell](start=ones(Ncell))
        "Velocity (average)";
        Modelica.SIunits.MassFlowRate m(start=1) "Mass flow shell side";
        Medium.ThermodynamicState state[Ncell];
        Modelica.SIunits.ReynoldsNumber Re[Ncell]
        "Reynolds number tubes (average)";
        Modelica.SIunits.PrandtlNumber Pr[Ncell] "Prandtl number tubes";
        Modelica.SIunits.SpecificEnthalpy h[2*Ncell]
        "Tube stream 2*Ncell enthalpies";
    equation
        Dhyd_o = Dhyd + 2*thick;
        G_wall = 2*lambda_wall/log(Dhyd_o/Dhyd);

        for i in 1:Ncell loop
          u[i]              = m/state[i].d/Aflow "tube velocity";
          state[i]          = Medium.setState_ph(p_in, (h[2*i - 1] + h[2*i])/2)
          "thermodynamic states cold cells";
          Re[i]             = Useful_functions.Pure_numbers.Reynolds(
                                                    u[i], state[i].d,  state[i].eta, Dhyd)
          "Reynolds number tubes";
          Pr[i]             = Useful_functions.Pure_numbers.Prandtl(
                                                   state[i].cp, state[i].eta, state[i].lambda)
          "Prandtl number tubes";
        end for;
      annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},
                {100,100}}), graphics={
            Line(
              points={{-60,20},{60,20}},
              color={0,0,0},
              smooth=Smooth.None),
            Line(
              points={{-60,-40},{60,-40}},
              color={0,0,0},
              smooth=Smooth.None),
            Ellipse(extent={{54,20},{64,-40}}, lineColor={0,0,0}),
            Ellipse(extent={{-66,20},{-56,-40}}, lineColor={0,0,0})}));
    end tubes;

    class cell "I am a cell element for the discretization of heat exchangers "
       parameter Integer Ncell "Number of cell elements";
       Modelica.SIunits.MassFlowRate mdot_hot "Hot mass flows";
       Modelica.SIunits.MassFlowRate mdot_cold "Cold mass flows";
       Modelica.SIunits.SpecificEnthalpy h_hot[2*Ncell]
        "Hot stream 2*Ncell enthalpies";
       Modelica.SIunits.SpecificEnthalpy h_cold[2*Ncell]
        "Cold stream 2*Ncell enthalpies";
       Integer pin_hot[Ncell]
        "Pins to tell the direction of heat of the hot stream";
       Integer pin_cold[Ncell]
        "Pins to tell the direction of heat of the cold stream";
       Modelica.SIunits.Temperature  T_hot[Ncell]
        "Temperatures of the shell cells";
       Modelica.SIunits.Temperature  T_cold[Ncell]
        "Temperatures of the tube cells";
       Modelica.SIunits.HeatFlowRate qdot[Ncell] "Heat rates";
       Modelica.SIunits.ThermalConductance htA_1[Ncell]
        "Thermal conductance shell side";
       Modelica.SIunits.ThermalConductance htA_2[Ncell]
        "Thermal conductance tube side";
       Modelica.SIunits.ThermalConductance G_wall[Ncell]
        "Thermal conductance of the tube wall";

    equation
      for j in 1:Ncell loop
           qdot[j]   = pin_hot[j]*mdot_hot*(h_hot[2*j - 1] - h_hot[2*j])
          "heat balance shell side";
           qdot[j]   = pin_cold[j]*mdot_cold*(h_cold[2*j] - h_cold[2*j - 1])
          "heat balance tube side";
           qdot[j]   = 1/(1/htA_1[j] + 1/G_wall[j] + 1/htA_2[j])*(T_hot[j] - T_cold[j])
          "heat transfer equation";
           assert(T_hot[j] > T_cold[j], "Check II principle of Thermodynamics", AssertionLevel.warning);
           //Here I get the sign of the heat across the tube cells
      end for;

      annotation (experiment(__Dymola_NumberOfIntervals=1, Tolerance=0.001),
          __Dymola_experimentSetupOutput,
        Icon(graphics={
            Rectangle(extent={{-100,100},{100,-100}}, lineColor={0,0,0}),
            Ellipse(
              extent={{-94,-6},{-106,6}},
              lineColor={0,0,0},
              fillPattern=FillPattern.Solid),
            Ellipse(
              extent={{106,-6},{94,6}},
              lineColor={0,0,0},
              fillPattern=FillPattern.Solid),
            Ellipse(
              extent={{6,94},{-6,106}},
              lineColor={0,0,0},
              fillPattern=FillPattern.Solid),
            Ellipse(
              extent={{6,-106},{-6,-94}},
              lineColor={0,0,0},
              fillPattern=FillPattern.Solid),
            Line(
              points={{-40,0},{40,0},{40,0}},
              color={0,0,0},
              smooth=Smooth.None),
            Line(
              points={{32,8},{40,0},{32,-8}},
              color={0,0,0},
              smooth=Smooth.None),
            Line(
              points={{0,40},{0,-40},{0,-40},{0,-40}},
              color={0,0,0},
              smooth=Smooth.None),
            Line(
              points={{0,-40},{8,-32}},
              color={0,0,0},
              smooth=Smooth.None),
            Line(
              points={{0,-40},{-8,-32}},
              color={0,0,0},
              smooth=Smooth.None)}));
    end cell;
  end Objects;

  package Useful_functions "Package containing useful functions"
    function log_mean_delta_T "logarithmic mean temperature difference"
      input Real T1;
      input Real T2;
      input Real t1;
      input Real t2;
      output Real LMTD;

    algorithm
      LMTD :=((T1 - t2) - (T2 - t1))/log((T1 - t2)/(T2 - t1));
    end log_mean_delta_T;

    package Shell
      "Functions to evaluate the shell clearance based on the diameter of the tube bundle "
      function Fixed_Utube "Fixed_U_tube shape clearance shell-bundle diameter"

        input Real d_b "bundle diameter";
        output Real clearance "bundle clearance";

      algorithm
        clearance :=1e-3*(10 + 10*(d_b - 0.2));
      end Fixed_Utube;

      function OPH "Outside packed head clearance shell-bundle diameter"

        input Real d_b "bundle diameter";
        output Real clearance "bundle clearance";

      algorithm
        clearance :=38e-3;
      end OPH;

      function SRFH "Split ring floating head clearance shell-bundle diameter"

        input Real d_b "bundle diameter";
        output Real clearance "bundle clearance";

      algorithm
        clearance :=1e-3*(50 + 28*(d_b - 0.2));
      end SRFH;

      function PTFH
        "Pull through floating head clearance shell-bundle diameter"

        input Real d_b "bundle diameter";
        output Real clearance "bundle clearance";

      algorithm
        clearance :=1e-3*(88 + 11*(d_b - 0.2));
      end PTFH;

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
    end Shell;

    package Pure_numbers "Package containing adimensional numbers"
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
    end Pure_numbers;
  end Useful_functions;

  package Media "I am the packag containing the definition of the fluids"

    package OneRandomOrganicFluid
      "Change the name to something more appropriate"
        extends ExternalMedia.Media.FluidPropMedium(
        mediumName = "Name of the fluid for documentation purposes",
        libraryName = "FluidProp.RefProp",
        substanceNames = {"cyclopentane"},
        ThermoStates = Modelica.Media.Interfaces.PartialMedium.Choices.IndependentVariables.ph);
    end OneRandomOrganicFluid;

    package Methanol_CoolProp "CoolProp model of Methanol"
      extends ExternalMedia.Media.CoolPropMedium(
        mediumName = "methanol",
        substanceNames = {"methanol|enable_TTSE=1|enable_EXTTP=1"},
        ThermoStates = Modelica.Media.Interfaces.PartialMedium.Choices.IndependentVariables.ph);
    end Methanol_CoolProp;

    package Water_Refprop "Water"
        extends ExternalMedia.Media.FluidPropMedium(
        mediumName = "water",
        libraryName = "FluidProp.Refprop",
        substanceNames = {"water"},
        ThermoStates = Modelica.Media.Interfaces.PartialMedium.Choices.IndependentVariables.ph);
    end Water_Refprop;

    package Water_CoolProp "CoolProp model of Water"
      extends ExternalMedia.Media.CoolPropMedium(
        mediumName = "Water",
        substanceNames = {"Water|enable_TTSE=1|enable_EXTTP=1"},
        ThermoStates = Modelica.Media.Interfaces.PartialMedium.Choices.IndependentVariables.ph);
    end Water_CoolProp;
  end Media;

  package Topologies
    "This package contains different topologies for heat exchangers"

    class shell_tube
      "Shell and tube topology. The shell fluid enters at the top and the tube enters at the bottom"
      parameter Integer N_baffles_d(start=3);
      parameter Integer Ncell(start=6);
      parameter Integer N_passes(start=2) "Number of tube passes";
      parameter Modelica.SIunits.SpecificEnthalpy h_s_in;
      parameter Modelica.SIunits.SpecificEnthalpy h_t_in;
      parameter Modelica.SIunits.SpecificEnthalpy h_s_out;
      output Modelica.SIunits.SpecificEnthalpy h_tube[2*Ncell]
        "Tube enthalpies";
      output Modelica.SIunits.SpecificEnthalpy h_shell[2*Ncell]
        "Shell enthalpies";
      output Integer pin_tube[Ncell]
        "Pins to tell the direction of heat of the tubes";
      output Integer pin_shell[Ncell]
        "Pins to tell the direction of heat of the shell";

    equation
          h_shell[2*N_baffles_d - 1]            = h_s_in;
          h_shell[2*Ncell - 2*N_baffles_d + 2]  = h_s_out;

          if (N_passes <> 1) then
            h_tube[2*Ncell - 2*N_baffles_d + 1] = h_t_in;
          else
            h_tube[2*N_baffles_d]               = h_t_in;
          end if;

          //Connecting the cells on the tube side (ok confirmed for all config.)
           for j in 2*N_baffles_d:4*N_baffles_d:2*Ncell - 2*N_baffles_d loop
               h_tube[j] = h_tube[j + 2*N_baffles_d];
           end for;

           for j in 2*N_baffles_d + 1:4*N_baffles_d:2*Ncell - 2*N_baffles_d loop
               h_tube[j] = h_tube[j + 2*N_baffles_d];
           end for;

          //Connecting the cells on the shell side (ok confirmed for all config.)
          for i in 1:N_passes loop
            for j in 2 + 2*N_baffles_d*(i-1):2:2*N_baffles_d*i-1 loop
              h_tube[j] = h_tube[j + 1];
              if (i <> N_passes) then
                h_shell[j] = h_shell[j + 2*N_baffles_d - 1];
              end if;
            end for;
            if (i <> N_passes) then
              h_shell[2*N_baffles_d*i] = h_shell[2*N_baffles_d*i + 2*N_baffles_d - 1];
            end if;
          end for;

           for j in 2*Ncell:-4:2*(Ncell - N_baffles_d) + 3 loop
             h_shell[j] = h_shell[j - 2];
           end for;

          for j in 2*N_baffles_d-3:-4:2 loop
            h_shell[j] = h_shell[j - 2];
          end for;

          for i in 1:N_passes loop
            for j in 1 + N_baffles_d*(i-1):N_baffles_d*i loop
              if mod(i,2) == 0 then
                pin_tube[j] = 1;
              else
                pin_tube[j] = -1;
              end if;
            end for;
          end for;
    //       for j in 1:2*N_baffles_d:Ncell loop
    //         for i in 1:N_baffles_d loop
    //             pin_tube[j + i - 1]               = -1;
    //             pin_tube[j + i - 1 + N_baffles_d] = 1;
    //         end for;
    //       end for;

          for i in N_baffles_d:-2:1 loop
            for j in N_passes:-1:1 loop
                pin_shell[i+N_baffles_d*(j-1)]   = 1;
            end for;
          end for;

          for i in N_baffles_d-1:-2:1 loop
            for j in N_passes:-1:1 loop
                pin_shell[i+N_baffles_d*(j-1)]   = -1;
            end for;
          end for;

    end shell_tube;
  end Topologies;

  package Components "Library with the design of the components "
    package Shell_and_tubes "Shell and tube heat exchangers"

      model shell_and_tube
        "U_type heat exchanger where the hot fluid flows on the shell and enters from the top. The cold fluid enters at the bottom."
        import constants = Modelica.Constants;
        //THE WORKING FLUIDS
        package Medium_s = VIP.Media.Methanol_CoolProp "Medium model";
        package Medium_t = VIP.Media.Water_CoolProp "Medium model";

        //GEOMETRY OF THE HEAT EXCHANGER
        parameter Modelica.SIunits.Length Dhyd = 16e-3
          "Hydraulic Diameter (single tube)";
        parameter Modelica.SIunits.Length thick = 2e-3 "Thickness of the tube";
        parameter Modelica.SIunits.Length l = 4.83 "Lenght (single tube)";
        parameter Real pitch_f = 1.25
          "Tube pitch as a fraction of the outer tube diameter";
        parameter Modelica.SIunits.ThermalConductivity lambda_wall = 50
          "Thermal conductivity of the tube wall";
        parameter Integer layout = 1 "Tube layout 1 = triangular, 2 = squared";
        parameter Integer N_passes = 2 "Number of tube passes";
        parameter Integer N_baffles = 27 "Number of baffles";
        parameter Integer N_baffles_d = 9
          "Number of baffles which are actually discretized (it should be N_baffle/2 or /3";
        parameter Integer Ncell = N_baffles_d*N_passes
          "Number of cell elements";
        parameter Integer allocation = 1
          "allocation = 1 hot fluid on the shell side";
        Modelica.SIunits.HeatFlowRate qtot "Heat flow rate";
        Modelica.SIunits.Temperature t_s_out;
        Modelica.SIunits.Temperature t_t_out;
        Modelica.SIunits.Length d_s "Shell diameter";
        Modelica.SIunits.Length d_s_eq "Equivalent shell diameter";
        Modelica.SIunits.Length l_b "Baffle lenght";
        Real N_tubes(start = 600) "Number of tubes in the bundle";
        Real N_t_p_p "Number of tubes per pass";
        Real bs_f "Baffle spacing as a fraction of the shell diameter";

        VIP.Objects.cell     mycells(Ncell=Ncell);

        VIP.Topologies.shell_tube topology(Ncell=Ncell,N_baffles_d=N_baffles_d,N_passes=N_passes,h_s_in=h_s_in, h_s_out=h_s_out, h_t_in=h_t_in,
                                             h_tube(start=linspace(h_t_in, h_t_out, 2*Ncell)), h_shell(start=linspace(h_s_in, h_s_out, 2*Ncell)));
        VIP.Objects.tube_bundle
                    mytubes(redeclare package Medium = Medium_t,
                    h(start=linspace(h_t_in, h_t_out, 2*Ncell)),
                    redeclare VIP.Heat_transfer.Tubes.Sieder_Tate heatTransfer(Re=mytubes.Re, Pr=mytubes.Pr, eta_wall=mytubes.state[1].eta*ones(Ncell)),
                    redeclare VIP.Pressure_drops.Tubes.Frank pressureDrop(l=l/N_baffles_d, Re=mytubes.Re, u=mytubes.u, heads = 2.5/Ncell),
                    redeclare function bundle_clearance =
                    VIP.Useful_functions.Shell.SRFH,
                    Ncell=N_baffles_d*N_passes, Aflow = 0.25*constants.pi*Dhyd^2,
                    Dhyd=Dhyd, l=l, thick=thick, N_tubes=N_tubes, N_passes=N_passes, layout=layout);

        VIP.Objects.tubes
                    myshell(redeclare package Medium = Medium_s,
                    h(start=linspace(h_s_in, h_s_out, 2*Ncell)),
                    redeclare VIP.Heat_transfer.Shell.Kern heatTransfer(Re=myshell.Re, Pr=myshell.Pr),
                    redeclare VIP.Pressure_drops.Shell.Kern pressureDrop(l=l/N_baffles_d/N_passes, d_s=d_s, Dhyd=d_s_eq, Re=mytubes.Re, u=myshell.u, l_b = l_b),
                    Ncell=N_baffles_d*N_passes, Aflow = (1 - 1/pitch_f)*d_s*l_b,
                    Dhyd=d_s_eq, l=l, thick=thick);

      protected
        parameter Modelica.SIunits.CoefficientOfHeatTransfer ht_t_f[Ncell] = ones(Ncell)*3e3
          "fouling heat transfer coefficient (tube side)";
        parameter Modelica.SIunits.CoefficientOfHeatTransfer ht_s_f[Ncell] = ones(Ncell)*5e3
          "fouling heat transfer coefficient (tube side)";
         parameter Modelica.SIunits.MassFlowRate m_s = 1e5/3.6e3
          "Mass flow shell side";
         parameter Medium_s.SaturationProperties state_sat_s = Medium_s.setSat_T(368.15)
          "Thermodynamic state at the inlet of the hot side";
         parameter Modelica.SIunits.AbsolutePressure p_s_in = state_sat_s.psat
          "Inlet pressure shell side";
         parameter Modelica.SIunits.SpecificEnthalpy h_s_in= Medium_s.bubbleEnthalpy(state_sat_s)
          "Inlet specific enthalpy shell side";
         parameter Modelica.SIunits.AbsolutePressure p_s_out = p_s_in
          "Outlet pressure shell side";
         parameter Modelica.SIunits.SpecificEnthalpy h_s_out = Medium_s.specificEnthalpy_pT(p_s_out, 40 + 273.15)
          "Outlet specific enthalpy shell side";

         parameter Modelica.SIunits.MassFlowRate m_t = 68.9
          "Mass flow tube side";
         parameter Modelica.SIunits.AbsolutePressure p_t_in = 1.0e5
          "Inlet pressure tube side";
         parameter Modelica.SIunits.SpecificEnthalpy h_t_in = Medium_t.specificEnthalpy_pT(p_t_in, 25 + 273.15)
          "Inlet specific enthalpy tube side";
         parameter Modelica.SIunits.SpecificEnthalpy h_t_out = Medium_t.specificEnthalpy_pT(p_t_out, 40 + 273.15)
          "Outlet specific enthalpy tube side";
         parameter Modelica.SIunits.AbsolutePressure p_t_out = 1.0e5
          "Outlet pressure tube side";

      equation
          mytubes.p_in = p_t_in;
          mytubes.m    = m_t/N_t_p_p;
          myshell.p_in = p_s_in;
          myshell.m    = m_s;

          qtot       = m_s*(h_s_in - h_s_out);

          if (allocation == 1) then
            mycells.mdot_cold    = m_t;
            mycells.mdot_hot     = m_s;
          else
             mycells.mdot_cold    = m_s;
             mycells.mdot_hot     = m_t;
          end if;

          for j in 1:2*Ncell loop
            myshell.h[j] = topology.h_shell[j];
            mytubes.h[j] = topology.h_tube[j];
             if (allocation == 1) then
              mycells.h_cold[j]    = topology.h_tube[j];
              mycells.h_hot[j]     = topology.h_shell[j];
             else
               mycells.h_cold[j]    = topology.h_shell[j];
               mycells.h_hot[j]     = topology.h_tube[j];
             end if;
          end for;

          for j in 1:Ncell loop
            mycells.T_hot[j]  = myshell.state[j].T;
            mycells.T_cold[j] = mytubes.state[j].T;
            mycells.htA_1[j]    = myshell.At/(1/myshell.heatTransfer.ht[j] + 1/ht_s_f[j]);
            mycells.htA_2[j]    = mytubes.At/(1/mytubes.heatTransfer.ht[j] + 1/ht_t_f[j]);
            mycells.G_wall[j]   = mytubes.At*mytubes.G_wall/mytubes.Dhyd;
             if (allocation == 1) then
              mycells.pin_cold[j]    = topology.pin_tube[j];
              mycells.pin_hot[j]     = topology.pin_shell[j];
             else
               mycells.pin_cold[j]    = topology.pin_shell[j];
               mycells.pin_hot[j]     = topology.pin_tube[j];
             end if;
          end for;

          d_s                     = mytubes.d_b + mytubes.clearance
          "Shell diameter";
          l_b                     = l/N_baffles "Baffle spacing";
          bs_f                    = 1e2*l_b/d_s
          "Baffle spacing as a percentage og the shell diameter 20 - 100 per cent";

          //Equivalent shell diameter
          if layout == 1 then
            d_s_eq  = 1.1*mytubes.Dhyd_o*(pitch_f^2 - 0.917);
          elseif layout == 2 then
            d_s_eq  = 1.27*mytubes.Dhyd_o*(pitch_f^2 - 0.785);
          else
            d_s_eq  = 1;
          end if;

          mytubes.At     = N_t_p_p*constants.pi*l/N_baffles_d*mytubes.Dhyd
          "I see this heat transfer area if I am inside the tubes";
          N_t_p_p        = floor(N_tubes/N_passes) "Number of tubes per pass";
          myshell.At     = N_t_p_p*constants.pi*l/N_baffles_d*mytubes.Dhyd_o
          "I see this heat transfer area if I am outside the tubes";

           //The outlet temperatures just to check the outlet conditions
           t_s_out         = Medium_s.temperature_ph(p_s_out, myshell.h[Ncell + 2]);
           t_t_out         = Medium_t.temperature_ph(p_t_out, mytubes.h[1]);

        annotation (experiment(
            __Dymola_NumberOfIntervals=1,
            Tolerance=0.001,
            __Dymola_Algorithm="Dassl"),
            __Dymola_experimentSetupOutput,
          Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},
                  {100,100}}),
              graphics={
              Rectangle(extent={{-60,-20},{66,-20}}, lineColor={85,170,255}),
              Rectangle(extent={{-60,20},{66,20}}, lineColor={85,170,255}),
              Rectangle(
                extent={{-40,30},{80,-30}},
                lineColor={0,0,0},
                lineThickness=0.5),
              Rectangle(
                extent={{-36,-30},{-24,-42}},
                lineColor={0,0,0},
                lineThickness=0.5),
              Rectangle(
                extent={{64,42},{76,30}},
                lineColor={0,0,0},
                lineThickness=0.5),
              Line(
                points={{66,20},{66,-20}},
                color={85,170,255},
                smooth=Smooth.None),
              Line(
                points={{0,-20},{-4,-16}},
                color={85,170,255},
                smooth=Smooth.None),
              Line(
                points={{0,-20},{-4,-24}},
                color={85,170,255},
                smooth=Smooth.None),
              Line(
                points={{0,20},{4,24}},
                color={85,170,255},
                smooth=Smooth.None),
              Line(
                points={{0,20},{4,16}},
                color={85,170,255},
                smooth=Smooth.None),
              Line(
                points={{70,60},{70,36}},
                color={255,0,0},
                smooth=Smooth.None),
              Line(
                points={{2,-2},{-2,2}},
                color={255,0,0},
                smooth=Smooth.None,
                origin={72,38},
                rotation=90),
              Line(
                points={{2,2},{-2,-2}},
                color={255,0,0},
                smooth=Smooth.None,
                origin={68,38},
                rotation=90),
              Line(
                points={{-30,-32},{-30,-56}},
                color={255,0,0},
                smooth=Smooth.None),
              Line(
                points={{2,-2},{-2,2}},
                color={255,0,0},
                smooth=Smooth.None,
                origin={-28,-54},
                rotation=90),
              Line(
                points={{2,2},{-2,-2}},
                color={255,0,0},
                smooth=Smooth.None,
                origin={-32,-54},
                rotation=90)}));
      end shell_and_tube;
    end Shell_and_tubes;
  end Components;

  package Tests "I test the components here"
    model validation_Aspen "Verification with the results given by Aspen"

      VIP.Components.Shell_and_tubes.shell_and_tube U_type__hot_up_cold_down(
        N_baffles=24,
        l=4.32,
        N_tubes(fixed=true, start=800),
        N_baffles_d=48)
        annotation (Placement(transformation(extent={{-108,-74},{88,66}})));
    end validation_Aspen;

    model validation_Coulson "Verification with the results given by Coulson"

      VIP.Components.Shell_and_tubes.shell_and_tube U_type__hot_up_cold_down(N_tubes(
            fixed=true, start=800), N_baffles_d=9)
        annotation (Placement(transformation(extent={{-108,-74},{88,66}})));
    end validation_Coulson;

    model cell_method_VDI_example
      "Cell method verified with an example from VDI C1 pag. 48"
      import constants = Modelica.Constants;
      parameter Integer Ncell=4 "Number of cell elements";
      parameter Modelica.SIunits.Temp_C t_s_in = 100
        "Inlet temperature of the hot source";
      parameter Modelica.SIunits.Temp_C t_t_in = 20
        "Inlet temperature of the cold source";
      parameter Modelica.SIunits.CoefficientOfHeatTransfer kA_tot = 4749/Ncell;
      parameter Modelica.SIunits.HeatCapacity C_hot = 3.5e3
        "Capacity of the hot stream";
      parameter Modelica.SIunits.HeatCapacity C_cold = 3.5e3
        "Capacity of the cold stream";
      Modelica.SIunits.HeatFlowRate qdot[Ncell] "Heat rate";
      Modelica.SIunits.Temperature t_hot[2*Ncell](start=linspace(t_s_in,60,2*Ncell))
        "Hot stream linear distribution with a guess";
      Modelica.SIunits.Temperature t_cold[2*Ncell](start=linspace(t_t_in,60,2*Ncell))
        "Cold stream linear distribution with a guess";
      Integer pin_hot[Ncell]
        "A pin to tell the direction of heat of the hot stream";
      Integer pin_cold[Ncell]
        "A pin to tell the direction of heat of the cold stream";
      Real adim_hot[2*Ncell] "Adimensional temperatures of the hot stream";
      Real adim_cold[2*Ncell] "Adimensional temperatures of the cold stream";
      Integer index[2*Ncell] "Index for plotting purposes";

    equation
        //Energy balances and heat transfer across the wall for each cell
        for j in 1:Ncell loop
           if (pin_hot[j] == 1) then
             qdot[j]   = C_hot*(t_hot[2*j - 1] - t_hot[2*j]);
           else
             qdot[j]   = C_hot*(t_hot[2*j] - t_hot[2*j - 1]);
           end if;
           if (pin_cold[j] == 1) then
             qdot[j]   = C_cold*(t_cold[2*j] - t_cold[2*j - 1]);
           else
             qdot[j]   = C_cold*(t_cold[2*j -1] - t_cold[2*j]);
           end if;
           qdot[j]   = kA_tot*((t_hot[2*j - 1] + t_hot[2*j])/2 - (t_cold[2*j - 1] + t_cold[2*j])/2);
         end for;

         //Boundary conditions
         t_hot[3]  = t_s_in;
         t_cold[5] = t_t_in;

         //The flow path in the tubes
         t_cold[2]  = t_cold[3];
         t_cold[4]  = t_cold[8];
         t_cold[7] = t_cold[6];
         pin_cold[1] = 0;
         pin_cold[2] = 0;
         pin_cold[3] = 1;
         pin_cold[4] = 1;

         //The flow path in the shell
         t_hot[4] = t_hot[7];
         t_hot[6] = t_hot[8];
         t_hot[5] = t_hot[2];
         pin_hot[1] = 0;
         pin_hot[3] = 0;
         pin_hot[2] = 1;
         pin_hot[4] = 1;

         //Calculate adimensional temperatures
         for i in 1:2*Ncell loop
           index[i] = i;
           adim_hot[i] = (t_hot[i] - t_t_in)/(t_s_in - t_t_in);
           adim_cold[i] = (t_cold[i] - t_t_in)/(t_s_in - t_t_in);
         end for;

      annotation (experiment(__Dymola_NumberOfIntervals=1, Tolerance=1e-009),
          __Dymola_experimentSetupOutput,
        Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,100}}),
            graphics));
    end cell_method_VDI_example;
  end Tests;
  annotation (uses(Modelica(version="3.2")), Icon(coordinateSystem(
          preserveAspectRatio=false, extent={{-100,-100},{100,100}}), graphics));
end VIP;
