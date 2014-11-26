within ;
package VIP "I am a package for the Virtual Prototyping Environment"
  package Heat_transfer "A package containing heat transfer correlations"
    package Tubes "heat transfer correlations in tubes"
      class Dittus_Boelter "Dittus Boelter correlation for tubes"
          replaceable package Medium = DynamicORCDesign.Media.Methanol_CoolProp
          "Medium model";
          parameter Integer Ncell(start=3) "Number of cell elements";
          parameter Modelica.SIunits.Length l "Lenght (single tube)";
          input Modelica.SIunits.Length Dhyd "Hydraulic Diameter (single tube)";
          input Modelica.SIunits.Area Aflow
          "Cross-sectional area (single tube)";
          Modelica.SIunits.CoefficientOfHeatTransfer ht[Ncell]
          "Heat transfer coefficient tubes";
          Medium.ThermodynamicState state[Ncell];
          Modelica.SIunits.ReynoldsNumber Re[Ncell]
          "Reynolds number tubes (average)";
          input Real alfa "exponent for the Prandtl number";
          Modelica.SIunits.PrandtlNumber Pr[Ncell] "Prandtl number tubes";
          Modelica.SIunits.NusseltNumber Nu[Ncell] "Nusselt number tubes ";
          Modelica.SIunits.SpecificEnthalpy h[2*Ncell]
          "Tube stream 2*Ncell enthalpies";
          Modelica.SIunits.AbsolutePressure p_in "Inlet pressure tube side";
          Modelica.SIunits.Velocity u[Ncell](start=ones(Ncell))
          "Velocity (average)";
          Modelica.SIunits.MassFlowRate m(start=1) "Mass flow shell side";

      equation
            for i in 1:Ncell loop
              u[i]              = m/state[i].d/Aflow "tube velocity";
              state[i]          = Medium.setState_ph(p_in, (h[2*i - 1] + h[2*i])/2)
            "thermodynamic states cold cells";
              Re[i]             = Pure_numbers.Reynolds(
                                                 u[i], state[i].d,  state[i].eta, Dhyd)
            "Reynolds number tubes";
              assert(Re[i] > 1e4, "Reynolds number is lower than 1e4 to use Dittus and Boelter", AssertionLevel.warning);
              Pr[i]             = Pure_numbers.Prandtl(
                                                state[i].cp, state[i].eta, state[i].lambda)
            "Prandtl number tubes";
              Nu[i]             =  2.3e-2*Re[i]^0.8*Pr[i]^alfa
            "Nusselt number tubes";
              ht[i]             = Pure_numbers.Nusselt(
                                                Nu[i], state[i].lambda, Dhyd)
            "Heat transfer coefficient tube side";
            end for;

      end Dittus_Boelter;

      class Sieder_Tate "Sieder Tate correlation for tubes"
          replaceable package Medium = DynamicORCDesign.Media.Methanol_CoolProp
          "Medium model";
          parameter Integer Ncell(start=3) "Number of cell elements";
          parameter Modelica.SIunits.Length l "Lenght (single tube)";
          input Modelica.SIunits.Length Dhyd "Hydraulic Diameter (single tube)";
          input Modelica.SIunits.Area Aflow
          "Cross-sectional area (single tube)";
          Modelica.SIunits.CoefficientOfHeatTransfer ht[Ncell]
          "Heat transfer coefficient tubes";
          Medium.ThermodynamicState state[Ncell];
          Modelica.SIunits.ReynoldsNumber Re[Ncell]
          "Reynolds number tubes (average)";
          input Modelica.SIunits.DynamicViscosity eta_wall[Ncell]
          "exponent for the Prandtl number";
          Modelica.SIunits.PrandtlNumber Pr[Ncell] "Prandtl number tubes";
          Modelica.SIunits.NusseltNumber Nu[Ncell] "Nusselt number tubes ";
          Modelica.SIunits.SpecificEnthalpy h[2*Ncell]
          "Tube stream 2*Ncell enthalpies";
          Modelica.SIunits.AbsolutePressure p_in "Inlet pressure tube side";
          Modelica.SIunits.Velocity u[Ncell](start=ones(Ncell))
          "Velocity (average)";
          Modelica.SIunits.MassFlowRate m(start=1) "Mass flow shell side";

      equation
            for i in 1:Ncell loop
              u[i]              = m/state[i].d/Aflow "tube velocity";
              state[i]          = Medium.setState_ph(p_in, (h[2*i - 1] + h[2*i])/2)
            "thermodynamic states cold cells";
              Re[i]             = Pure_numbers.Reynolds(
                                                 u[i], state[i].d,  state[i].eta, Dhyd)
            "Reynolds number tubes";
              assert(Re[i] > 1e4, "Reynolds number is lower than 1e4 to use Sieder and Tate", AssertionLevel.warning);
              assert(Pr[i] > 0.6,  "Prandtl number is lower than 0.6 to be used in Sieder and Tate", AssertionLevel.warning);
              Pr[i]             = Pure_numbers.Prandtl(
                                                state[i].cp, state[i].eta, state[i].lambda)
            "Prandtl number tubes";
              Nu[i]             = 2.3e-2*Re[i]^0.8*Pr[i]^(1.0/3)*(eta_wall[i]/state[i].eta)^0.14
            "Nusselt number tubes";
              ht[i]             = Pure_numbers.Nusselt(
                                                Nu[i], state[i].lambda, Dhyd)
            "Heat transfer coefficient tube side";
            end for;

      end Sieder_Tate;

      class Gnielinski "Gnielinski correlation for tubes"
          replaceable package Medium = DynamicORCDesign.Media.Methanol_CoolProp
          "Medium model";
          parameter Integer Ncell(start=3) "Number of cell elements";
          parameter Modelica.SIunits.Length l "Lenght (single tube)";
          input Modelica.SIunits.Length Dhyd "Hydraulic Diameter (single tube)";
          input Modelica.SIunits.Area Aflow
          "Cross-sectional area (single tube)";
          Modelica.SIunits.CoefficientOfHeatTransfer ht[Ncell]
          "Heat transfer coefficient tubes";
          Medium.ThermodynamicState state[Ncell];
          Modelica.SIunits.ReynoldsNumber Re[Ncell]
          "Reynolds number tubes (average)";
          Real csi[Ncell] "Friction factor";
          Modelica.SIunits.PrandtlNumber Pr[Ncell] "Prandtl number tubes";
          Modelica.SIunits.NusseltNumber Nu[Ncell] "Nusselt number tubes ";
          Modelica.SIunits.SpecificEnthalpy h[2*Ncell]
          "Tube stream 2*Ncell enthalpies";
          Modelica.SIunits.AbsolutePressure p_in "Inlet pressure tube side";
          Modelica.SIunits.Velocity u[Ncell](start=ones(Ncell))
          "Velocity (average)";
          Modelica.SIunits.MassFlowRate m(start=1) "Mass flow shell side";

      equation
            for i in 1:Ncell loop
              u[i]              = m/state[i].d/Aflow "tube velocity";
              state[i]          = Medium.setState_ph(p_in, (h[2*i - 1] + h[2*i])/2)
            "thermodynamic states cold cells";
              Re[i]             = Pure_numbers.Reynolds(
                                                 u[i], state[i].d,  state[i].eta, Dhyd)
            "Reynolds number tubes";
              assert(Re[i] > 1e4 or Re[i] < 1e6,  "Reynolds number is lower than 1e4 or greater than 1e6 to be used in Gnielinski", AssertionLevel.warning);
              assert(Pr[i] > 0.7 or Pr[i] < 1e3,  "Prandtl number is lower than 0.7 or greater than 1e3 to be used in Gnielinski", AssertionLevel.warning);

              Pr[i]             = Pure_numbers.Prandtl(
                                                state[i].cp, state[i].eta, state[i].lambda)
            "Prandtl number tubes";
              csi[i]            = 1/(0.78*log(Re[i]) - 1.5)^2 "Friction factor";
              Nu[i]             =  ((csi[i]/8)*Re[i]*Pr[i])/(1 + 12.7*sqrt(csi[i]/8)*(Pr[i]^(2.0/3) - 1))*(1 +  (Dhyd/l)^(2.0/3))
            "Nusselt number tubes";
              ht[i]             = Pure_numbers.Nusselt(
                                                Nu[i], state[i].lambda, Dhyd)
            "Heat transfer coefficient tube side";
            end for;

      end Gnielinski;

      class EagleFerguson
        "Eagle-Ferguson heat transfer correlation (only for liquid water)"
          replaceable package Medium = DynamicORCDesign.Media.Methanol_CoolProp
          "Medium model";
          parameter Integer Ncell(start=3) "Number of cell elements";
          parameter Modelica.SIunits.Length l "Lenght (single tube)";
          input Modelica.SIunits.Length Dhyd "Hydraulic Diameter (single tube)";
          input Modelica.SIunits.Area Aflow
          "Cross-sectional area (single tube)";
          Modelica.SIunits.CoefficientOfHeatTransfer ht[Ncell]
          "Heat transfer coefficient tubes";
          Medium.ThermodynamicState state[Ncell];
          Modelica.SIunits.NusseltNumber Nu[Ncell] "Nusselt number tubes ";
          Modelica.SIunits.SpecificEnthalpy h[2*Ncell]
          "Tube stream 2*Ncell enthalpies";
          Modelica.SIunits.AbsolutePressure p_in "Inlet pressure tube side";
          Modelica.SIunits.Velocity u[Ncell](start=ones(Ncell))
          "Velocity (average)";
          Modelica.SIunits.MassFlowRate m(start=1) "Mass flow shell side";

      equation
            for i in 1:Ncell loop
              u[i]              = m/state[i].d/Aflow "tube velocity";
              state[i]          = Medium.setState_ph(p_in, (h[2*i - 1] + h[2*i])/2)
            "thermodynamic states cold cells";
              Nu[i]             =  4.2e3*Dhyd/state[i].lambda*(1.35 + 2e-2*(state[i].T - 273.15))*u[i]^0.8/(1e3*Dhyd)^0.2
            "Nusselt number tubes";
              ht[i]             = Pure_numbers.Nusselt(
                                                Nu[i], state[i].lambda, Dhyd)
            "Heat transfer coefficient tube side";
            end for;
      end EagleFerguson;
    end Tubes;

    package Shell "heat transfer correlations in shells"

      class Kern_shell "Kern correlation for shell side"
          replaceable package Medium = DynamicORCDesign.Media.Methanol_CoolProp
          "Medium model";
          parameter Integer Ncell(start=3) "Number of cell elements";
          parameter Modelica.SIunits.Length l "Lenght (single tube)";
          parameter Real x[2] = {1e5, 1e1};
          parameter Real y[2] = {2e-3, 1.98e-1};
          input Modelica.SIunits.Length Dhyd "Hydraulic Diameter (single tube)";
          input Modelica.SIunits.Area Aflow
          "Cross-sectional area (single tube)";
          Modelica.SIunits.CoefficientOfHeatTransfer ht[Ncell]
          "Heat transfer coefficient tubes";
          Medium.ThermodynamicState state[Ncell];
          Modelica.SIunits.ReynoldsNumber Re[Ncell]
          "Reynolds number tubes (average)";
          Modelica.SIunits.PrandtlNumber Pr[Ncell] "Prandtl number tubes";
          Modelica.SIunits.NusseltNumber Nu[Ncell] "Nusselt number tubes ";
          Modelica.SIunits.SpecificEnthalpy h[2*Ncell]
          "Tube stream 2*Ncell enthalpies";
          Modelica.SIunits.AbsolutePressure p_in "Inlet pressure tube side";
          Modelica.SIunits.Velocity u[Ncell](start=ones(Ncell))
          "Velocity (average)";
          Modelica.SIunits.MassFlowRate m(start=1) "Mass flow shell side";

      equation
            for i in 1:Ncell loop
              u[i]              = m/state[i].d/Aflow "tube velocity";
              state[i]          = Medium.setState_ph(p_in, (h[2*i - 1] + h[2*i])/2)
            "thermodynamic states cold cells";
              Re[i]             = Pure_numbers.Reynolds(
                                                 u[i], state[i].d,  state[i].eta, Dhyd)
            "Reynolds number tubes";
              assert(Re[i] > 1e4, "Reynolds number is lower than 1e4 to use Dittus and Boelter", AssertionLevel.warning);
              Pr[i]             = Pure_numbers.Prandtl(
                                                state[i].cp, state[i].eta, state[i].lambda)
            "Prandtl number tubes";
              Nu[i]             = 10^log10(2e-3/1e5^(log10(y[2]/y[1])/log10(x[2]/x[1])))*Re[i]^(log10(y[2]/y[1])/log10(x[2]/x[1]))*Re[i]*Pr[i]^(1/3)
            "Nusselt number tubes";
              ht[i]             = Pure_numbers.Nusselt(
                                                Nu[i], state[i].lambda, Dhyd)
            "Heat transfer coefficient tube side";
            end for;
      end Kern_shell;
    end Shell;
  end Heat_transfer;

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

  package Objects "Package containing all the objects of the VIP"
    class cell "I am a cell element for the discretization of heat exchangers "
       parameter Integer Ncell(start=3) "Number of cell elements";
       Modelica.SIunits.HeatFlowRate qdot[Ncell] "Heat rates";
       Modelica.SIunits.Temperature T_hot[Ncell](start=ones(Ncell)*400);
       Modelica.SIunits.Temperature T_cold[Ncell](start=ones(Ncell)*300);
       Modelica.SIunits.SpecificEnthalpy h_hot[2*Ncell]
        "Hot stream 2*Ncell enthalpies";
       Modelica.SIunits.SpecificEnthalpy h_cold[2*Ncell]
        "Cold stream 2*Ncell enthalpies";
       Modelica.SIunits.AbsolutePressure p_hot
        "Hot pressure for now no pressure drops";
       Modelica.SIunits.AbsolutePressure p_cold
        "Cold pressure for now no pressure drops";
       Modelica.SIunits.MassFlowRate mdot_hot "Hot mass flows";
       Modelica.SIunits.MassFlowRate mdot_cold "Cold mass flows";
       Integer pin_hot[Ncell]
        "Pins to tell the direction of heat of the hot stream";
       Integer pin_cold[Ncell]
        "Pins to tell the direction of heat of the cold stream";
       Modelica.SIunits.ThermalConductance htA_hot[Ncell]
        "Thermal conductance hot side";
       Modelica.SIunits.ThermalConductance htA_cold[Ncell]
        "Thermal conductance cold side";
       Modelica.SIunits.ThermalConductance G_wall[Ncell](start=ones(Ncell)*1e9)
        "Thermal conductance of the tube wall";

    equation
      for j in 1:Ncell loop
           qdot[j]   = pin_hot[j]*mdot_hot*(h_hot[2*j - 1] - h_hot[2*j])
          "heat balance hot side";
           qdot[j]   = pin_cold[j]*mdot_cold*(h_cold[2*j] - h_cold[2*j - 1])
          "heat balance cold side";
           qdot[j]   = 1/(1/htA_hot[j] + 1/G_wall[j] + 1/htA_cold[j])*(T_hot[j] - T_cold[j])
          "heat transfer equation";
           assert(T_hot[j] > T_cold[j], "Check II principle of Thermodynamics", AssertionLevel.warning);
           //Here I get the sign of the heat across the cold cells
      end for;

      annotation (experiment(__Dymola_NumberOfIntervals=1, Tolerance=0.001),
          __Dymola_experimentSetupOutput);
    end cell;

    class tube_bundle
      "I am a tube bundle and I contain all the my relevant informations"
        extends VIP.Objects.tubes;
        parameter Integer N_passes "Number of tube passes";
        parameter Integer layout "Tube layout 1 = triangular, 2 = squared";
        parameter Modelica.SIunits.Length clearance = 68e-3
        "Bundle - shell clearance";
        parameter Real pitch_f = 1.25
        "Tube pitch as a fraction of the outer tube diameter";
        input Real N_tubes "Number of tubes in the bundle";
        Modelica.SIunits.Length d_b "Bundle Diameter";
    equation
        d_b      = Useful_functions.bundle_diameter(
                                         N_tubes, N_passes, Dhyd_o, layout);
    end tube_bundle;

    class tubes "I am a tube and I contain all the my relevant informations"
        import constants = Modelica.Constants "to import pi-greco";
        //replaceable HeatTransfer.IdealHeatTransfer heatTransfer
        replaceable VIP.Heat_transfer.Tubes.EagleFerguson
                                        heatTransfer(Ncell = Ncell, Dhyd = Dhyd, l=l, Aflow=Aflow);
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
    equation
        Dhyd_o = Dhyd + 2*thick;
        G_wall = 2*lambda_wall/log(Dhyd_o/Dhyd);
    end tubes;
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
  end Useful_functions;

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

  model shell_and_tube_1D_Coulson
    "Cell method to design a shell and tube heat exchanger (design)"
    import constants = Modelica.Constants;
    //THE WORKING FLUIDS
    package Medium_s = DynamicORCDesign.Media.Methanol_CoolProp "Medium model";
    package Medium_t = DynamicORCDesign.Media.Water_CoolProp "Medium model";

    //GEOMETRY OF THE HEAT EXCHANGER
    parameter Modelica.SIunits.Length Dhyd = 16e-3
      "Hydraulic Diameter (single tube)";
    parameter Modelica.SIunits.Length thick = 2e-3 "Thickness of the tube";
    parameter Modelica.SIunits.Length l = 4.83 "Lenght (single tube)";
    parameter Modelica.SIunits.Length clearance = 68e-3
      "Bundle - shell clearance";
    parameter Real pitch_f = 1.25
      "Tube pitch as a fraction of the outer tube diameter";
    parameter Modelica.SIunits.ThermalConductivity lambda_wall = 50
      "Thermal conductivity of the tube wall";
    parameter Integer layout = 1 "Tube layout 1 = triangular, 2 = squared";
    parameter Integer N_passes = 2 "Number of tube passes";
    parameter Integer N_baffles = 27 "Number of baffles";
    parameter Integer N_baffles_d = 9
      "Number of baffles which are actually discretized (it should be N_baffle/2 or /3";
    parameter Integer Ncell = N_baffles_d*N_passes "Number of cell elements";

    Modelica.SIunits.HeatFlowRate qtot "Heat flow rate";
    Modelica.SIunits.Temperature t_s_out;
    Modelica.SIunits.Temperature t_t_out;
    Modelica.SIunits.Length d_s "Shell diameter";
    Modelica.SIunits.Length d_s_eq "Equivalent shell diameter";
    Modelica.SIunits.Length l_b "Baffle lenght";
    Real N_tubes(start = 900) "Number of tubes in the bundle";
    Real N_t_p_p(start = N_tubes/N_passes) "Number of tubes per pass";
    Real bs_f(start = 20) "Baffle spacing as a fraction of the shell diameter";

    VIP.Objects.cell
               mycells(Ncell=Ncell);
    VIP.Objects.tube_bundle
                      mytubes(redeclare VIP.Heat_transfer.Tubes.Dittus_Boelter
                                                             heatTransfer(redeclare
          package Medium = Medium_t,h(start=linspace(h_t_in, h_t_out, 2*Ncell)),alfa = 0.4),
                         Ncell=N_baffles_d*N_passes, Aflow = 0.25*constants.pi*Dhyd^2,
                         Dhyd=Dhyd, l=l, thick=thick, N_tubes=N_tubes, N_passes=N_passes, layout=layout);

    VIP.Objects.tubes
                myshell(redeclare VIP.Heat_transfer.Shell.Kern_shell
                                                   heatTransfer(redeclare
          package Medium =
                   Medium_s, h(start=linspace(h_s_in, h_s_out, 2*Ncell))),
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

    parameter Modelica.SIunits.MassFlowRate m_t = 68.9 "Mass flow tube side";
    parameter Modelica.SIunits.AbsolutePressure p_t_in = 1.0e5
      "Inlet pressure tube side";
    parameter Modelica.SIunits.SpecificEnthalpy h_t_in = Medium_t.specificEnthalpy_pT(p_t_in, 25 + 273.15)
      "Inlet specific enthalpy tube side";
    parameter Modelica.SIunits.SpecificEnthalpy h_t_out = Medium_t.specificEnthalpy_pT(p_t_out, 40 + 273.15)
      "Outlet specific enthalpy tube side";
    parameter Modelica.SIunits.AbsolutePressure p_t_out = 1.0e5
      "Outlet pressure tube side";

  equation
      //Boundary conditions when projecting the heat exchanger
      mycells.h_hot[Ncell - 1]  = h_s_in;
      mycells.h_cold[Ncell + 1] = h_t_in;
      mycells.h_hot[Ncell + 2]  = h_s_out;
      mycells.mdot_cold         = m_t;
      mycells.p_cold            = p_t_in;
      mycells.mdot_hot          = m_s;
      mycells.p_hot             = p_s_in;
      qtot                      = m_s*(h_s_in - h_s_out);

      mytubes.heatTransfer.p_in = p_t_in;
      mytubes.heatTransfer.m = m_t/N_t_p_p;

      myshell.heatTransfer.p_in = p_s_in;
      myshell.heatTransfer.m = m_s;

      for j in 1:2*Ncell loop
        myshell.heatTransfer.h[j] = mycells.h_hot[j];
        mytubes.heatTransfer.h[j] = mycells.h_cold[j];
      end for;

      for j in N_baffles_d:-2:1 loop
          mycells.pin_hot[j] = 1;
          mycells.pin_hot[j + N_baffles_d] = 1;
      end for;
      for j in N_baffles_d-1:-2:1 loop
          mycells.pin_hot[j] = -1;
          mycells.pin_hot[j + N_baffles_d] = -1;
      end for;

      //Connecting the cells on the cold side
      for j in 2:2:2*Ncell-1 loop
        if (j == Ncell) then
          mycells.h_cold[Ncell] = mycells.h_cold[2*Ncell];
        else
          mycells.h_cold[j] = mycells.h_cold[j + 1];
        end if;
      end for;

      //Connecting the cells on the hot side
      for j in 2:2:Ncell loop
        mycells.h_hot[j] = mycells.h_hot[j + Ncell - 1];
      end for;

      for j in 2*Ncell:-4:Ncell+4 loop
        mycells.h_hot[j] = mycells.h_hot[j - 2];
      end for;

      for j in Ncell - 3:-4:2 loop
        mycells.h_hot[j] = mycells.h_hot[j - 2];
      end for;

      for j in 1:Ncell loop
        mycells.T_hot[j]   = myshell.heatTransfer.state[j].T;
        mycells.T_cold[j]  = mytubes.heatTransfer.state[j].T;
        if (j <= N_baffles_d) then
           mycells.pin_cold[j] = -1;
         else
           mycells.pin_cold[j] = 1;
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

      for i in 1:Ncell loop
        //All the heat transfer conductances here
        mycells.htA_hot[i]    = myshell.At/(1/myshell.heatTransfer.ht[i] + 1/ht_s_f[i]);
        mycells.htA_cold[i]   = mytubes.At/(1/mytubes.heatTransfer.ht[i] + 1/ht_t_f[i]);
        mycells.G_wall[i]     = mytubes.At*mytubes.G_wall/mytubes.Dhyd;
      end for;

        //The outlet temperatures just to check the outlet conditions
        t_s_out         = Medium_s.temperature_ph(p_s_out, mycells.h_hot[Ncell + 2]);
        t_t_out         = Medium_t.temperature_ph(p_t_out, mycells.h_cold[1]);

    annotation (experiment(
        __Dymola_NumberOfIntervals=1,
        Tolerance=0.001,
        __Dymola_Algorithm="Dassl"),
        __Dymola_experimentSetupOutput,
      Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,100}}),
          graphics={Rectangle(extent={{-32,40},{44,-44}}, lineColor={0,0,255})}));
  end shell_and_tube_1D_Coulson;
  annotation (uses(Modelica(version="3.2")));
end VIP;
