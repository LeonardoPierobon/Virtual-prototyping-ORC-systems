within ;
package VIP "I am a package for the Virtual Prototyping Environment"
  import Modelica.Math.*;
  import Modelica.SIunits.*;
  import Modelica.Constants.*;

  package Icons "Icons for the virtual prototpying environment"
    partial model shell_tube "Shell and tube icon"
      extends VIP.Icons.tube;
      extends VIP.Icons.shell;
      annotation (Icon(graphics={
            Rectangle(extent={{-80,-10},{46,-10}}, lineColor={85,170,255},
              lineThickness=0.5),
            Rectangle(extent={{-80,30},{46,30}}, lineColor={85,170,255},
              lineThickness=0.5),
            Rectangle(
              extent={{-60,40},{60,-20}},
              lineColor={0,0,0},
              lineThickness=1),
            Rectangle(
              extent={{-56,-20},{-44,-32}},
              lineColor={0,0,0},
              lineThickness=1),
            Rectangle(
              extent={{44,52},{56,40}},
              lineColor={0,0,0},
              lineThickness=1),
            Line(
              points={{46,30},{46,-10}},
              color={85,170,255},
              smooth=Smooth.None,
              thickness=0.5),
            Line(
              points={{-20,-10},{-24,-6}},
              color={85,170,255},
              smooth=Smooth.None,
              thickness=0.5),
            Line(
              points={{-20,-10},{-24,-14}},
              color={85,170,255},
              smooth=Smooth.None,
              thickness=0.5),
            Line(
              points={{-20,30},{-16,34}},
              color={85,170,255},
              smooth=Smooth.None,
              thickness=0.5),
            Line(
              points={{-20,30},{-16,26}},
              color={85,170,255},
              smooth=Smooth.None,
              thickness=0.5),
            Line(
              points={{50,70},{50,46}},
              color={255,0,0},
              smooth=Smooth.None,
              thickness=0.5),
            Line(
              points={{2,-2},{-2,2}},
              color={255,0,0},
              smooth=Smooth.None,
              origin={52,48},
              rotation=90,
              thickness=0.5),
            Line(
              points={{2,2},{-2,-2}},
              color={255,0,0},
              smooth=Smooth.None,
              origin={48,48},
              rotation=90,
              thickness=0.5),
            Line(
              points={{-50,-22},{-50,-46}},
              color={255,0,0},
              smooth=Smooth.None,
              thickness=0.5),
            Line(
              points={{2,-2},{-2,2}},
              color={255,0,0},
              smooth=Smooth.None,
              origin={-48,-44},
              rotation=90,
              thickness=0.5),
            Line(
              points={{2,2},{-2,-2}},
              color={255,0,0},
              smooth=Smooth.None,
              origin={-52,-44},
              rotation=90,
              thickness=0.5)}));
    end shell_tube;

    partial model tube "tube icon"

      annotation (Icon(graphics={
            Rectangle(extent={{-80,-10},{46,-10}}, lineColor={85,170,255},
              lineThickness=0.5),
            Rectangle(extent={{-80,30},{46,30}}, lineColor={85,170,255},
              lineThickness=0.5),
            Line(
              points={{46,30},{46,-10}},
              color={85,170,255},
              smooth=Smooth.None,
              thickness=0.5),
            Line(
              points={{-20,-10},{-24,-6}},
              color={85,170,255},
              smooth=Smooth.None,
              thickness=0.5),
            Line(
              points={{-20,-10},{-24,-14}},
              color={85,170,255},
              smooth=Smooth.None,
              thickness=0.5),
            Line(
              points={{-20,30},{-16,34}},
              color={85,170,255},
              smooth=Smooth.None,
              thickness=0.5),
            Line(
              points={{-20,30},{-16,26}},
              color={85,170,255},
              smooth=Smooth.None,
              thickness=0.5)}));
    end tube;

    partial model shell "Shell icon"

      annotation (Icon(graphics={
            Rectangle(
              extent={{-60,40},{60,-20}},
              lineColor={0,0,0},
              lineThickness=1),
            Rectangle(
              extent={{-56,-20},{-44,-32}},
              lineColor={0,0,0},
              lineThickness=1),
            Rectangle(
              extent={{44,52},{56,40}},
              lineColor={0,0,0},
              lineThickness=1),
            Line(
              points={{50,70},{50,46}},
              color={255,0,0},
              smooth=Smooth.None,
              thickness=0.5),
            Line(
              points={{2,-2},{-2,2}},
              color={255,0,0},
              smooth=Smooth.None,
              origin={52,48},
              rotation=90,
              thickness=0.5),
            Line(
              points={{2,2},{-2,-2}},
              color={255,0,0},
              smooth=Smooth.None,
              origin={48,48},
              rotation=90,
              thickness=0.5),
            Line(
              points={{-50,-22},{-50,-46}},
              color={255,0,0},
              smooth=Smooth.None,
              thickness=0.5),
            Line(
              points={{2,-2},{-2,2}},
              color={255,0,0},
              smooth=Smooth.None,
              origin={-48,-44},
              rotation=90,
              thickness=0.5),
            Line(
              points={{2,2},{-2,-2}},
              color={255,0,0},
              smooth=Smooth.None,
              origin={-52,-44},
              rotation=90,
              thickness=0.5)}));
    end shell;
  end Icons;

  package Heat_transfer "A package containing heat transfer correlations"
    package Tubes "heat transfer correlations in tubes"
      class Dittus_Boelter "Dittus Boelter correlation for tubes"
          extends VIP.Heat_transfer.Tubes.basic_ht;
          output Modelica.SIunits.NusseltNumber Nu[Ncell]
          "Nusselt number tubes ";
          output Modelica.SIunits.CoefficientOfHeatTransfer ht[Ncell]
          "Heat transfer coefficient tubes";
          input Modelica.SIunits.Length Dhyd "Hydraulic Diameter (single tube)";
          input Real alfa "exponent for the Prandtl number";
          Modelica.SIunits.ReynoldsNumber Re[Ncell]
          "Reynolds number tubes (average)";
          Modelica.SIunits.PrandtlNumber Pr[Ncell] "Prandtl number tubes";
          Modelica.SIunits.Velocity u[Ncell](start=ones(Ncell))
          "Velocities inside the tubes";
      equation
            for i in 1:Ncell loop
              u[i]              = mdot/state[i].d/Aflow "tube velocity";
              Re[i]             = Useful_functions.Pure_numbers.Reynolds(u[i], state[i].d,  state[i].eta, Dhyd)
            "Reynolds number tubes";
              Pr[i]             = Useful_functions.Pure_numbers.Prandtl(state[i].cp, state[i].eta, state[i].lambda)
            "Prandtl number tubes";
              assert(Re[i] > 1e4, "Reynolds number is lower than 1e4 to use Dittus and Boelter", AssertionLevel.warning);
              Nu[i]             =  2.3e-2*Re[i]^0.8*Pr[i]^alfa
            "Nusselt number tubes";
              ht[i]             = Useful_functions.Pure_numbers.Nusselt(
                                                       Nu[i], state[i].lambda, Dhyd)
            "Heat transfer coefficient tube side";
            end for;

      end Dittus_Boelter;

      class Sieder_Tate "Sieder Tate correlation for tubes"
          extends VIP.Heat_transfer.Tubes.basic_ht;
          output Modelica.SIunits.NusseltNumber Nu[Ncell]
          "Nusselt number tubes ";
          output Modelica.SIunits.CoefficientOfHeatTransfer ht[Ncell]
          "Heat transfer coefficient tubes";
          input Modelica.SIunits.Length Dhyd "Hydraulic Diameter (single tube)";
          input Modelica.SIunits.DynamicViscosity eta_wall[Ncell]
          "exponent for the viscosity correction";
          Modelica.SIunits.ReynoldsNumber Re[Ncell]
          "Reynolds number tubes (average)";
          Modelica.SIunits.PrandtlNumber Pr[Ncell] "Prandtl number tubes";
          Modelica.SIunits.Velocity u[Ncell](start=ones(Ncell))
          "Velocities inside the tubes";
      equation
            for i in 1:Ncell loop
              u[i]              = mdot/state[i].d/Aflow "tube velocity";
              Re[i]             = Useful_functions.Pure_numbers.Reynolds(u[i], state[i].d,  state[i].eta, Dhyd)
            "Reynolds number tubes";
              Pr[i]             = Useful_functions.Pure_numbers.Prandtl(state[i].cp, state[i].eta, state[i].lambda)
            "Prandtl number tubes";
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
          extends VIP.Heat_transfer.Tubes.basic_ht;
          parameter Modelica.SIunits.Length l "Lenght (single tube)";
          output Modelica.SIunits.NusseltNumber Nu[Ncell]
          "Nusselt number tubes ";
          output Modelica.SIunits.CoefficientOfHeatTransfer ht[Ncell]
          "Heat transfer coefficient tubes";
          input Modelica.SIunits.Length Dhyd "Hydraulic Diameter (single tube)";
          Modelica.SIunits.ReynoldsNumber Re[Ncell]
          "Reynolds number tubes (average)";
          Modelica.SIunits.PrandtlNumber Pr[Ncell] "Prandtl number tubes";
          Modelica.SIunits.Velocity u[Ncell](start=ones(Ncell))
          "Velocities inside the tubes";
          Real csi[Ncell] "Friction factor";
      equation
            for i in 1:Ncell loop
              u[i]              = mdot/state[i].d/Aflow "tube velocity";
              Re[i]             = Useful_functions.Pure_numbers.Reynolds(u[i], state[i].d,  state[i].eta, Dhyd)
            "Reynolds number tubes";
              Pr[i]             = Useful_functions.Pure_numbers.Prandtl(state[i].cp, state[i].eta, state[i].lambda)
            "Prandtl number tubes";
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
          extends VIP.Heat_transfer.Tubes.basic_ht;
          output Modelica.SIunits.NusseltNumber Nu[Ncell]
          "Nusselt number tubes ";
          input Modelica.SIunits.Length Dhyd "Hydraulic Diameter (single tube)";
          Modelica.SIunits.Velocity u[Ncell](start=ones(Ncell))
          "Velocity (average)";
          Modelica.SIunits.CoefficientOfHeatTransfer ht[Ncell]
          "Heat transfer coefficient tubes";

      equation
            for i in 1:Ncell loop
              u[i]              = mdot/state[i].d/Aflow "tube velocity";
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
          input Medium.ThermodynamicState state[Ncell];
          input Modelica.SIunits.Area Aflow
          "Cross-sectional area (single tube)";
          input Modelica.SIunits.MassFlowRate mdot "Mass flow rate";
      end basic_ht;
    end Tubes;

    package Shell "heat transfer correlations in shells"

      class Kern "Kern correlation for shell side"
          extends VIP.Heat_transfer.Tubes.basic_ht;
          output Modelica.SIunits.NusseltNumber Nu[Ncell]
          "Nusselt number tubes ";
          output Modelica.SIunits.CoefficientOfHeatTransfer ht[Ncell]
          "Heat transfer coefficient tubes";
          parameter Modelica.SIunits.Length Dhyd_o
          "Outer Diameter (single tube)";
          parameter Real pitch_f
          "Tube pitch as a fraction of the outer tube diameter";
          parameter Integer layout = 1
          "Tube layout 1 = triangular, 2 = squared";
          Modelica.SIunits.ReynoldsNumber Re[Ncell](start=10e6*ones(Ncell))
          "Reynolds number tubes (average)";
          Modelica.SIunits.Length d_s_eq "Equivalent shell diameter";
          Modelica.SIunits.PrandtlNumber Pr[Ncell] "Prandtl number tubes";
          Modelica.SIunits.Velocity u[Ncell](start=ones(Ncell))
          "Velocities inside the tubes";

      equation
          //Equivalent shell diameter
          if layout == 1 then
            d_s_eq  = 1.1*Dhyd_o*(pitch_f^2 - 0.917);
          elseif layout == 2 then
            d_s_eq  = 1.27*Dhyd_o*(pitch_f^2 - 0.785);
          else
            d_s_eq  = 1;
          end if;

          for i in 1:Ncell loop
            u[i]              = mdot/state[i].d/Aflow "tube velocity";
            Re[i]             = Useful_functions.Pure_numbers.Reynolds(u[i], state[i].d,  state[i].eta, d_s_eq)
            "Reynolds number tubes";
            Pr[i]             = Useful_functions.Pure_numbers.Prandtl(state[i].cp, state[i].eta, state[i].lambda)
            "Prandtl number tubes";
            Nu[i]             = 10^log10(2e-3/1e5^(-0.25*log10(1.98/2e-2)))*Re[i]^(-0.25*log10(1.98/2e-2))*Re[i]*Pr[i]^(1/3)
            "Nusselt number tubes";
            ht[i]             = Useful_functions.Pure_numbers.Nusselt(Nu[i], state[i].lambda, d_s_eq)
            "Heat transfer coefficient tube side";
            end for;
      end Kern;


      class Bell_Delaware "Bell Delaware correlation for shell side"
          extends VIP.Heat_transfer.Tubes.basic_ht;
          output Modelica.SIunits.CoefficientOfHeatTransfer ht[Ncell]
          "Heat transfer coefficient tubes";
          input Modelica.SIunits.Length Dhyd "Hydraulic Diameter (single tube)";
          input Modelica.SIunits.Length Dhyd_o "Outer Diameter (single tube)";
          parameter Real b_cut "Baffle cut";
          final parameter Real ttb = 8e-4
          "tube to baffle clearance (from Standards)";
          final parameter Real bts = 4.8e-3
          "baffle to shell clearance (from Standards)";
          parameter Real pitch_f
          "Tube pitch as a fraction of the outer tube diameter";
          parameter Real N_ss
          "The number of sealing strips (pairs) in one baffle spacing";
          parameter Integer layout "Tube layout 1 = triangular, 2 = squared";
          parameter Integer N_passes "Number of tube passes";
          parameter Real ff "Fouling factor (=1 if clean)";
          input Modelica.SIunits.Length d_s "shell diameter";
          input Modelica.SIunits.Length d_b "bundle diameter";
          input Real N_tubes "Number of tubes in the bundle";
          input Modelica.SIunits.Length l_b "Baffle lenght";
          input Modelica.SIunits.DynamicViscosity eta_wall[Ncell]
          "exponent for the viscosity correction";
          Modelica.SIunits.ReynoldsNumber Re[Ncell]
          "Reynolds number shell (average)";
          Modelica.SIunits.PrandtlNumber Pr[Ncell] "Prandtl number tubes";
          Modelica.SIunits.Velocity u[Ncell](start=ones(Ncell))
          "Velocities inside the shell";
          Real csi[Ncell] "Friction factor";
          Modelica.SIunits.Length tpv "vertical tube pitch";
          Real N_tcc
          "Number of constrictions crossed=number of tube rows between the baffle tips";
          Real N_tcw
          "The effective number of tube rows crossed in the baffle window";
          Real F_sbp
          "The ratio of the bypass area, Sb, to the overall crossflow area";
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
          Modelica.SIunits.Angle teta_ctl "The upper centriangle of baffle cut";

      protected
           Real a1[Ncell];
           Real a2[Ncell];
           Real a[Ncell];
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
        J_l      = 0.44*(1 - r_s) + (1 - 0.44*(1 - r_s))*exp(-2.2*r_lm);

        S_b      = l_b*(d_s - d_b);
        F_sbp    = S_b/S_m;
        J_b      = exp(-1.35*F_sbp*(1 - (2/N_ss)^(1/3)));

        F_c      = 1 - 2*F_w;
        J_c      = F_c + 0.54*(1 - F_c)^0.345;
        S_w      = S_wg - S_wt;
        N_tcw    = 0.8*(d_s*b_cut - 0.5*(d_s - (d_b - Dhyd_o)))/tpv;

        for i in 1:Ncell loop

          u[i]    = mdot/state[i].d/S_m;
          Re[i]   = Useful_functions.Pure_numbers.Reynolds(u[i], state[i].d,  state[i].eta, Dhyd_o);
          Pr[i]   = Useful_functions.Pure_numbers.Prandtl(state[i].cp, state[i].eta, state[i].lambda);

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

          csi[i] = a1[i]*(1.33*pitch_f)^a[i]*Re[i]^a2[i];
          ht[i]  = J_l*J_b*J_c*csi[i]*state[i].cp*mdot/(S_m*Pr[i]^(2/3));
        end for;

      end Bell_Delaware;
    end Shell;
  end Heat_transfer;

  package Pressure_drops "A package containing pressure drops correlations"
    package Tubes "heat transfer correlations in tubes"
      class Frank "Frank correlation for tubes"
          extends VIP.Heat_transfer.Tubes.basic_ht;
          output Modelica.SIunits.AbsolutePressure dp[Ncell]
          "Pressure drops cells";
          output Modelica.SIunits.AbsolutePressure dp_tot
          "Pressure drops tubes";
          input Modelica.SIunits.Length l "Lenght (single tube)";
          input Modelica.SIunits.Length Dhyd "Hydraulic Diameter (single tube)";
          input Modelica.SIunits.AbsolutePressure heads(start=2.5)
          "Number of velocity heads";
          Modelica.SIunits.ReynoldsNumber Re[Ncell]
          "Reynolds number tubes (average)";
          Modelica.SIunits.Velocity u[Ncell](start=ones(Ncell))
          "Velocity (average)";
          Real csi[Ncell] "Friction factor";

      equation
        for i in 1:Ncell loop
          u[i]     = mdot/state[i].d/Aflow "tube velocity";
          Re[i]    = Useful_functions.Pure_numbers.Reynolds(u[i], state[i].d,  state[i].eta, Dhyd)
            "Reynolds number tubes";
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
          input Medium.ThermodynamicState state[Ncell];
      end basic_dp;
    end Tubes;

    package Shell "heat transfer correlations in shells"

      class Kern "Kern correlation for shell"
          extends VIP.Heat_transfer.Tubes.basic_ht;
          Modelica.SIunits.AbsolutePressure dp[Ncell] "Pressure drops cells";
          Modelica.SIunits.AbsolutePressure dp_tot "Pressure drops tubes";
          input Modelica.SIunits.Length l "Lenght (single tube)";
          input Modelica.SIunits.Length Dhyd "Hydraulic Diameter (single tube)";
          input Modelica.SIunits.Length d_s "Shell diameter";
          input Modelica.SIunits.Length l_b "Baffle length";
          input Modelica.SIunits.DynamicViscosity eta_wall[Ncell]
          "exponent for the viscosity correction";
          Modelica.SIunits.ReynoldsNumber Re[Ncell](start=10e6*ones(Ncell))
          "Reynolds number tubes (average)";
          Modelica.SIunits.Velocity u[Ncell](start=ones(Ncell))
          "Velocity (average)";
          Real csi[Ncell] "Friction factor";

      equation
        for i in 1:Ncell loop
          u[i]     = mdot/state[i].d/Aflow "tube velocity";
          Re[i]    = Useful_functions.Pure_numbers.Reynolds(u[i], state[i].d,  state[i].eta, Dhyd)
            "Reynolds number tubes";
          if (Re[i] < 3e2) then
            csi[i] = 10^log10(2.5e-1/1e2^(log10(2.4/1e-1)/log10(10/3e2))) * Re[i]^(log10(2.4/1e-1)/log10(10/3e2));
          else
            csi[i] = 10^log10(4e-2/4e4^(log10(1e-1/2.4e-2)/log10(3e2/1e6))) * Re[i]^(log10(1e-1/2.4e-2)/log10(3e2/1e6));
          end if;
          dp[i]             =  4*csi[i]*d_s/Dhyd*l/l_b*(eta_wall[i]/state[i].eta)^0.14*state[i].d*u[i]^2
            "Local pressure drops";
        end for;

        dp_tot = sum(dp) "Total pressure drops";

      end Kern;

      class Wills_Johnston "Wills and Johnston correlation for shell"
          extends VIP.Heat_transfer.Tubes.basic_ht;
          output Modelica.SIunits.AbsolutePressure dp[Ncell]
          "Pressure drops tubes";
          output Modelica.SIunits.AbsolutePressure dp_tot
          "Pressure drops tubes";
          input Modelica.SIunits.Length l "Lenght (single tube)";
          input Modelica.SIunits.Length Dhyd "Hydraulic Diameter (single tube)";
          input Modelica.SIunits.Length d_s "Shell diameter";
          input Modelica.SIunits.Length l_b "Baffle length";
          input Modelica.SIunits.DynamicViscosity eta_wall[Ncell]
          "exponent for the viscosity correction";
          Modelica.SIunits.ReynoldsNumber Re[Ncell](start=10e6*ones(Ncell))
          "Reynolds number tubes (average)";
          Modelica.SIunits.Velocity u[Ncell](start=ones(Ncell))
          "Velocity (average)";
          Real csi[Ncell] "Friction factor";

      equation
        for i in 1:Ncell loop
          u[i]     = mdot/state[i].d/Aflow "tube velocity";
          Re[i]    = Useful_functions.Pure_numbers.Reynolds(u[i], state[i].d,  state[i].eta, Dhyd)
            "Reynolds number tubes";
          assert(Re[i] > 4e2, "Reynolds number is lower than 1e2 to use Wills_Johnston", AssertionLevel.warning);
          csi[i] = 1.7789*Re[i]^(-0.195868);
          dp[i]  =  0.5*csi[i]*d_s/Dhyd*l/l_b*(eta_wall[i]/state[i].eta)^0.14*state[i].d*u[i]^2
            "Local pressure drops";
        end for;

        dp_tot = sum(dp) "Total pressure drops";

      end Wills_Johnston;


      class Bell_Delaware "Bell Delaware correlation for shell side"
          extends VIP.Heat_transfer.Tubes.basic_ht;
          input Modelica.SIunits.Length Dhyd "Hydraulic Diameter (single tube)";
          input Modelica.SIunits.Length Dhyd_o "Outer Diameter (single tube)";
          parameter Real b_cut "Baffle cut";
          final parameter Real ttb = 8e-4
          "tube to baffle clearance (from Standards)";
          final parameter Real bts = 4.8e-3
          "baffle to shell clearance (from Standards)";
          parameter Real pitch_f
          "Tube pitch as a fraction of the outer tube diameter";
          parameter Real N_ss
          "The number of sealing strips (pairs) in one baffle spacing";
          parameter Integer layout "Tube layout 1 = triangular, 2 = squared";
          parameter Integer N_passes "Number of tube passes";
          parameter Real ff "Fouling factor (=1 if clean)";
          input Modelica.SIunits.Length d_s "shell diameter";
          input Modelica.SIunits.Length d_b "bundle diameter";
          input Real N_tubes "Number of tubes in the bundle";
          input Modelica.SIunits.Length l_b "Baffle lenght";
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
          Modelica.SIunits.ReynoldsNumber Re[Ncell](start=10e6*ones(Ncell))
          "Reynolds number shell (average)";
          Modelica.SIunits.Velocity u[Ncell](start=ones(Ncell))
          "Velocities inside the shell";
          Real csi[Ncell] "Friction factor";
          Modelica.SIunits.Length tpv "vertical tube pitch";
          Real N_tcc
          "Number of constrictions crossed=number of tube rows between the baffle tips";
          Real N_tcw
          "The effective number of tube rows crossed in the baffle window";
          Real F_sbp
          "The ratio of the bypass area, Sb, to the overall crossflow area";
          Real F_w "Fraction of number of tubes in the baffle window";
          Real R_l
          "Correction factors for baffle leakage effects for heat transfer";
          Real R_b
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
          Modelica.SIunits.Angle teta_ctl "The upper centriangle of baffle cut";

      protected
           Real b1[Ncell];
           Real b2[Ncell];
           Real b[Ncell];
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
        R_l      = exp(-1.33*(1 + r_s))*r_lm^(-0.15*(1 + r_s) + 0.8);

        S_b      = l_b*(d_s - d_b);
        F_sbp    = S_b/S_m;
        R_b      = exp(-4.5*F_sbp*(1 - (2/N_ss)^(1/3)));

        S_w      = S_wg - S_wt;
        N_tcw    = 0.8*(d_s*b_cut - 0.5*(d_s - (d_b - Dhyd_o)))/tpv;

        for i in 1:Ncell loop

          u[i]    = mdot/state[i].d/S_m;
          Re[i]   = Useful_functions.Pure_numbers.Reynolds(u[i], state[i].d,  state[i].eta, Dhyd_o);

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

          csi[i] = b1[i]*(1.33*pitch_f)^b[i]*Re[i]^b2[i];

          if (i > Ncell - N_passes) then
            dp_w[i]  = 0;
          else
            dp_w[i]  = (2 + 0.6*N_tcw)*mdot^2/(2*S_w*S_m*state[i].d)/N_passes;
          end if;

          if (i <= N_passes or i > Ncell - N_passes) then
            dp_c[i] = 0;
            dp_e[i] = 2*(N_tcc + N_tcw)*csi[i]*(state[i].d*u[i]^2/2)/N_passes;
          else
            dp_e[i] = 0;
            dp_c[i] = ff*2*N_tcc*csi[i]*(state[i].d*u[i]^2/2)/N_passes;
          end if;

        end for;

        dp_w_tot = sum(dp_w)*R_l;
        dp_c_tot = sum(dp_c)*R_b*R_l;
        dp_e_tot = sum(dp_e)*R_b;

        dp_tot = dp_w_tot + dp_c_tot + dp_e_tot "Total pressure drops";

      end Bell_Delaware;
    end Shell;
  end Pressure_drops;

  package Objects "Package containing all the objects of the VIP"

    class tubes "I am a tube and I contain all the my relevant informations"
        extends VIP.Icons.tube;
        replaceable package Medium = VIP.Media.OneRandomOrganicFluid
        "Medium model";
        parameter Integer Ncell(start=3) "Number of cell elements";
        parameter Modelica.SIunits.Length thick "Thickness of the tube";
        parameter Modelica.SIunits.ThermalConductivity lambda_wall = 50
        "Thermal conductivity of the tube wall";
        parameter Modelica.SIunits.SpecificEnthalpy h_in
        "Inlet specific enthalpy shell side";
        parameter Modelica.SIunits.SpecificEnthalpy h_out
        "Outlet specific enthalpy shell side";
        parameter Modelica.SIunits.MassFlowRate mdot "Mass flow rate total";
        parameter Real pin "This pin identifies if the fluid is hot or cold";
        input Modelica.SIunits.Area Aflow "Cross-sectional area (single tube)";
        input Modelica.SIunits.Length Dhyd "Hydraulic Diameter (single tube)";
        input Modelica.SIunits.MassFlowRate mdot_pt
        "Mass flow rate per tube!!!!";
        Modelica.SIunits.Area At "Heat transfer area associated to the tubes";
        Modelica.SIunits.Length Dhyd_o "Outer Diameter (single tube)";
        Modelica.SIunits.ThermalConductivity G_wall;
        Modelica.SIunits.AbsolutePressure p_in "Inlet pressure tube side";
        Medium.ThermodynamicState state[Ncell]
        "Thermodynamic states of the cells";
        Modelica.SIunits.SpecificEnthalpy h[Ncell + 1](start=linspace(h_in, h_out, Ncell + 1))
        "Tube stream Ncell+1 enthalpies";
        Modelica.SIunits.HeatFlowRate qdot[Ncell] "Heat rate of each cell";

    equation
        //Trivial calculations
        Dhyd_o = Dhyd + 2*thick;
        G_wall = 2*lambda_wall/log(Dhyd_o/Dhyd);

        for i in 1:Ncell loop
          qdot[i]           = pin*mdot*(h[i] - h[i + 1]);
          state[i]          = Medium.setState_ph(p_in, (h[i] + h[i + 1])/2);
        end for;
      annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},
                {100,100}}), graphics));
    end tubes;

    class tube_bundle
      "I am a tube bundle and I contain all the my relevant information"
        extends VIP.Objects.tubes;
        parameter Integer N_passes "Number of tube passes";
        parameter Integer layout "Tube layout 1 = triangular, 2 = squared";
        parameter Real pitch_f
        "Tube pitch as a fraction of the outer tube diameter";
        input Real N_tubes "Number of tubes in the bundle";
        input Modelica.SIunits.Length clearance "Bundle - shell clearance";
        Modelica.SIunits.Length d_b "Bundle Diameter";
    equation

        d_b       = Useful_functions.Shell.bundle_diameter(N_tubes, N_passes, Dhyd_o, layout)
        "calculate the tube bundle";

    end tube_bundle;
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
        mediumName = "REFPROP-methanol",
        substanceNames = {"REFPROP-methanol"},
        ThermoStates = Modelica.Media.Interfaces.PartialMedium.Choices.IndependentVariables.ph);
    end Methanol_CoolProp;

    package Water_RefProp "Water"
        extends ExternalMedia.Media.FluidPropMedium(
        mediumName = "water",
        libraryName = "FluidProp.Refprop",
        substanceNames = {"water"},
        ThermoStates = Modelica.Media.Interfaces.PartialMedium.Choices.IndependentVariables.ph);
    end Water_RefProp;

    package Water_CoolProp "CoolProp model of Water"
      extends ExternalMedia.Media.CoolPropMedium(
        mediumName = "Water",
        substanceNames = {"Water|enable_TTSE=1|enable_EXTTP=1"},
        ThermoStates = Modelica.Media.Interfaces.PartialMedium.Choices.IndependentVariables.ph);
    end Water_CoolProp;

    package Air_CoolProp "CoolProp model of Air"
      extends ExternalMedia.Media.CoolPropMedium(
        mediumName = "Water",
        substanceNames = {"Water|enable_TTSE=1|enable_EXTTP=1"},
        ThermoStates = Modelica.Media.Interfaces.PartialMedium.Choices.IndependentVariables.ph);
    end Air_CoolProp;

    package Methanol_RefProp "methanol"
        extends ExternalMedia.Media.FluidPropMedium(
        mediumName = "methanol",
        libraryName = "FluidProp.Refprop",
        substanceNames = {"methanol"},
        ThermoStates = Modelica.Media.Interfaces.PartialMedium.Choices.IndependentVariables.ph);
    end Methanol_RefProp;
  end Media;

  package Components "Library with the design of the components "
    package HEX "Heat exchangers"

      model shell_and_tube "Shell and tube heat exchanger where the hot fluid flows on the shell and enters from the top. The cold fluid enters at the bottom.
  It can model 1 tube pass or 2*N tube passes."
        extends VIP.Icons.shell_tube;

        //THE WORKING FLUIDS
        replaceable package Medium_s = VIP.Media.Methanol_CoolProp constrainedby
          Modelica.Media.Interfaces.PartialMedium "Medium model" annotation(choicesAllMatching = true);
        replaceable package Medium_t = VIP.Media.Water_CoolProp               constrainedby
          Modelica.Media.Interfaces.PartialMedium "Medium model" annotation(choicesAllMatching = true);

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
        parameter Integer N_baffles = 26 "Number of baffles";
        parameter Modelica.SIunits.CoefficientOfHeatTransfer U_guess = 600
          "Guess value for the global heat transfer coefficient";
        final parameter Integer N_baffles_d = N_baffles + 1
          "The number of discretization volumes is always N_baffles + 1";
        final parameter Integer Ncell = N_baffles_d*N_passes
          "Number of cell elements";
        final parameter Modelica.SIunits.HeatFlowRate qdot = m_s*abs(h_s_in - h_s_out)
          "Heat flow rate";
        final parameter Modelica.SIunits.Temp_C  DTML = VIP.Useful_functions.log_mean_delta_T(t_s_in, t_s_out, t_t_in, t_t_out)
          "Logarithmic mean temperature difference";
        Modelica.SIunits.Length d_s "Shell diameter";
        Modelica.SIunits.Length l_b "Baffle lenght";
        Real N_tubes( start = qdot/(DTML*U_guess)/(pi*l*(Dhyd + 2*thick)))
          "Number of tubes in the bundle";
        Real N_t_p_p "Number of tubes per pass";
        Real bs_f "Baffle spacing as a fraction of the shell diameter";

        //Heat transfer and pressure drop correlations
        replaceable VIP.Heat_transfer.Tubes.Sieder_Tate hT_tube(Dhyd=Dhyd, eta_wall=mytubes.state[1].eta*ones(Ncell))
           constrainedby VIP.Heat_transfer.Tubes.basic_ht(Medium = Medium_t, Ncell = Ncell, state=mytubes.state, mdot=m_t/N_t_p_p, Aflow=mytubes.Aflow) annotation(choicesAllMatching = true);
        replaceable VIP.Pressure_drops.Tubes.Frank dp_tube(Dhyd=Dhyd,l=l/N_baffles_d, heads = 2.5/Ncell)
          constrainedby VIP.Heat_transfer.Tubes.basic_ht(Medium = Medium_t, Ncell = Ncell, state=mytubes.state, mdot=m_t/N_t_p_p, Aflow=mytubes.Aflow) annotation(choicesAllMatching = true);

        replaceable VIP.Heat_transfer.Shell.Kern hT_shell(Dhyd_o=Dhyd+2*thick, layout=layout, pitch_f = pitch_f)
          constrainedby VIP.Heat_transfer.Tubes.basic_ht(Medium = Medium_s, Ncell = Ncell, state=myshell.state, mdot=m_s, Aflow=myshell.Aflow) annotation(choicesAllMatching = true);
        replaceable VIP.Pressure_drops.Shell.Wills_Johnston dp_shell(l=l/N_baffles_d/N_passes, d_s=d_s, Dhyd=hT_shell.d_s_eq, l_b = l_b,
                      eta_wall=myshell.state[1].eta*ones(Ncell))
         constrainedby VIP.Heat_transfer.Tubes.basic_ht(Medium = Medium_s, Ncell = Ncell, state=myshell.state, mdot=m_s, Aflow=myshell.Aflow) annotation(choicesAllMatching = true);

        //Defining the model for the bundle clearance
        replaceable function bundle_clearance =
            VIP.Useful_functions.Shell.Fixed_Utube;

        //Definiing the tubes and the shell
        VIP.Objects.tube_bundle
                    mytubes(redeclare package Medium = Medium_t,
                    h_in = h_t_in,
                    h_out = h_t_out,
                    clearance = bundle_clearance(mytubes.d_b),
                    Ncell=N_baffles_d*N_passes, Aflow = 0.25*pi*Dhyd^2, mdot = m_t,
                    Dhyd=Dhyd, thick = thick, lambda_wall = lambda_wall, N_tubes=N_tubes, mdot_pt = m_t/N_t_p_p,
                    N_passes=N_passes, layout=layout, pin=-1, pitch_f=pitch_f);

        VIP.Objects.tubes
                    myshell(redeclare package Medium = Medium_s,
                    h_in = h_s_in,
                    h_out = h_s_out,
                    Ncell=N_baffles_d*N_passes, Aflow = (1 - 1/pitch_f)*d_s*l_b, mdot = m_s, mdot_pt = m_s,
                    Dhyd=d_s, thick = thick, lambda_wall = lambda_wall, pin=1);

         //Boundary conditions at inlet and outlet
         parameter Modelica.SIunits.MassFlowRate m_s = 27.78
          "Mass flow shell side";
         parameter Modelica.SIunits.SpecificEnthalpy h_s_in = 4.25006e+08
          "Inlet specific enthalpy shell side";
         parameter Modelica.SIunits.SpecificEnthalpy h_s_out= 4.24849e+08
          "Outlet specific enthalpy shell side";
         parameter Modelica.SIunits.AbsolutePressure p_s_in = 301435
          "Inlet pressure shell side";
         parameter Modelica.SIunits.AbsolutePressure p_s_out = 301435
          "Outlet pressure shell side";
         parameter Modelica.SIunits.MassFlowRate m_t = 68.9
          "Mass flow tube side";
         parameter Modelica.SIunits.SpecificEnthalpy h_t_in = Medium_t.specificEnthalpy_pT(1e5,298.15)
          "Inlet specific enthalpy tube side";
         parameter Modelica.SIunits.SpecificEnthalpy h_t_out = Medium_t.specificEnthalpy_pT(1e5,313.15)
          "Outlet specific enthalpy tube side";
         parameter Modelica.SIunits.AbsolutePressure p_t_in = 1e5
          "Inlet pressure tube side";
         parameter Modelica.SIunits.AbsolutePressure p_t_out = 1e5
          "Outlet pressure tube side";
         final parameter Modelica.SIunits.Temperature  t_s_in = Medium_s.temperature_ph(p_s_in, h_s_in)
          "Inlet temperature shell side";
         final parameter Modelica.SIunits.Temperature  t_s_out = Medium_s.temperature_ph(p_s_out, h_s_out)
          "Outlet temperature shell side";
         final parameter Modelica.SIunits.Temperature  t_t_in = Medium_t.temperature_ph(p_t_in, h_t_in)
          "Inlet temperature tube side";
         final parameter Modelica.SIunits.Temperature  t_t_out = Medium_t.temperature_ph(p_t_out, h_t_out)
          "Outlet temperature tube side";
         parameter Modelica.SIunits.CoefficientOfHeatTransfer ht_t_f[Ncell] = ones(Ncell)*3e3
          "fouling heat transfer coefficient (tube side)";
         parameter Modelica.SIunits.CoefficientOfHeatTransfer ht_s_f[Ncell] = ones(Ncell)*5e3
          "fouling heat transfer coefficient (tube side)";

      protected
         Integer k[Ncell]
          "This index points at the array elements of the shell from the perspective of the tubes";
         Modelica.SIunits.ThermalConductance kA_tot[Ncell]
          "Overall thermal conductance";
         Modelica.SIunits.ThermalConductance htA_shell[Ncell]
          "Thermal conductance shell side";
         Modelica.SIunits.ThermalConductance htA_tubes[Ncell]
          "Thermal conductance tube side";
         Modelica.SIunits.ThermalConductance G_wall[Ncell]
          "Thermal conductance of the tube wall";

      equation
          //Set boundary conditions at the inlet and outelt
          mytubes.p_in = p_t_in;
          myshell.p_in = p_s_in;
          mytubes.h[1]         = h_t_in;
          myshell.h[1]         = h_s_in;
          myshell.h[Ncell + 1] = h_s_out;

          //Things get tricky here. We need to get the index of the shell cell seen by the tube cell
          for j in 1:N_passes loop
               for i in 1+(j-1)*N_baffles_d:2:j*N_baffles_d loop
                  if mod(j,2) <> 0 then
                    if mod(N_baffles_d,2) <> 0 then
                      k[i]            = Ncell - N_passes*(i-(1+(j-1)*N_baffles_d)) - (j-1);
                    else
                      k[i]            = Ncell - N_passes*(i-(1+(j-1)*N_baffles_d)) - (N_passes - j);
                    end if;
                  else
                      k[i]            = (N_passes + 1 - j) + N_passes*(i-(1+(j-1)*N_baffles_d));
                  end if;
               end for;
           end for;

           for j in 1:N_passes loop
               for i in 2+(j-1)*N_baffles_d:2:j*N_baffles_d loop
                  if mod(j,2) <> 0 then
                    if mod(N_baffles_d,2) <> 0 then
                      k[i]            = Ncell - (2*N_passes - 1) + (j-1) - N_passes*(i-(2+(j-1)*N_baffles_d));
                    else
                      k[i]            = Ncell - N_passes + (j-1) - N_passes*(i-(2+(j-1)*N_baffles_d));
                    end if;
                  else
                      k[i]            = N_passes + j + N_passes*(i-(2+(j-1)*N_baffles_d));
                  end if;
               end for;
           end for;

          //A sperate loop for the heat transfer and energy balances
          for i in 1:Ncell loop
            G_wall[i]       = mytubes.At*mytubes.G_wall/mytubes.Dhyd;
            htA_tubes[i]    = mytubes.At/(1/hT_tube.ht[i] + 1/ht_t_f[i]);
            htA_shell[k[i]] = myshell.At/(1/hT_shell.ht[k[i]] + 1/ht_s_f[k[i]]);
            kA_tot[i]       = 1/(1/htA_tubes[i] + 1/G_wall[i] + 1/htA_shell[k[i]]);
            mytubes.qdot[i] = myshell.qdot[k[i]];
            mytubes.qdot[i] = kA_tot[i]*(myshell.state[k[i]].T - mytubes.state[i].T);
          end for;

          //Geometric calculations
          d_s                     = mytubes.d_b + mytubes.clearance;
          l_b                     = l/(N_baffles+1);
          bs_f                    = 1e2*l_b/d_s
          "Baffle spacing as a percentage og the shell diameter 20 - 100 per cent";

          mytubes.At     = N_t_p_p*pi*l/N_baffles_d*mytubes.Dhyd
          "I see this heat transfer area if I am inside the tubes";
          N_t_p_p        = N_tubes/N_passes "Number of tubes per pass";
          myshell.At     = N_t_p_p*pi*l/N_baffles_d*mytubes.Dhyd_o
          "I see this heat transfer area if I am outside the tubes";

        annotation (experiment(
            __Dymola_NumberOfIntervals=1,
            Tolerance=1e-006,
            __Dymola_Algorithm="Dassl"),
            __Dymola_experimentSetupOutput,
          Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},
                  {100,100}}),
              graphics));
      end shell_and_tube;
    end HEX;
  end Components;

  package Tests "I test the components here"

    model Coulson_Kern
      "Verification with the results given by Coulson et al. using the Kern method"

      replaceable package Medium_s = VIP.Media.Methanol_CoolProp "Medium model";
      replaceable package Medium_t = VIP.Media.Water_CoolProp "Medium model";
      parameter Modelica.SIunits.SpecificEnthalpy h_s_in= 4.25006e+08
        "Inlet specific enthalpy shell side";
      parameter Modelica.SIunits.SpecificEnthalpy h_s_out= 4.24849e+08
        "Outlet specific enthalpy shell side";
      parameter Modelica.SIunits.SpecificEnthalpy h_t_in = Medium_t.specificEnthalpy_pT(1e5,298.15)
        "Inlet specific enthalpy tube side";
      parameter Modelica.SIunits.SpecificEnthalpy h_t_out = Medium_t.specificEnthalpy_pT(1e5,313.15)
        "Outlet specific enthalpy tube side";

      Components.HEX.shell_and_tube             shell_and_tube(
        Dhyd=16e-3,
        thick=2e-3,
        l=4.83,
        pitch_f=1.25,
        lambda_wall=50,
        layout=1,
        N_passes=2,
        N_baffles=27,
        h_s_in=h_s_in,
        h_s_out=h_s_out,
        h_t_in=h_t_in,
        h_t_out=h_t_out,
        m_s=1e2/3.6,
        m_t=68.9,
        p_s_in=301435,
        p_s_out=301435,
        p_t_in=100000,
        p_t_out=100000)
        annotation (Placement(transformation(extent={{-80,-68},{68,58}})));
     annotation (Placement(transformation(extent={{-108,-74},{88,66}})),
        experiment,
        __Dymola_experimentSetupOutput);
    end Coulson_Kern;

    model Coulson_Bell_Delaware
      "Verification with the results given by Coulson et al. using the Bell Delaware method."

      replaceable package Medium_s = VIP.Media.Methanol_CoolProp "Medium model";
      replaceable package Medium_t = VIP.Media.Water_CoolProp "Medium model";
      parameter Modelica.SIunits.SpecificEnthalpy h_s_in= Medium_s.specificEnthalpy_pT(3.2e5,368.15)
        "Inlet specific enthalpy shell side";
      parameter Modelica.SIunits.SpecificEnthalpy h_s_out= Medium_s.specificEnthalpy_pT(3.2e5,313.15)
        "Outlet specific enthalpy shell side";
      parameter Modelica.SIunits.SpecificEnthalpy h_t_in = Medium_t.specificEnthalpy_pT(1e5,298.15)
        "Inlet specific enthalpy tube side";
      parameter Modelica.SIunits.SpecificEnthalpy h_t_out = Medium_t.specificEnthalpy_pT(1e5,313.15)
        "Outlet specific enthalpy tube side";

      Components.HEX.shell_and_tube  shell_and_tube(
        Dhyd=16e-3,
        thick=2e-3,
        l=4.83,
        pitch_f=1.25,
        lambda_wall=50,
        layout=1,
        N_passes=2,
        N_baffles=12,
        h_s_in=h_s_in,
        h_s_out=h_s_out,
        h_t_in=h_t_in,
        h_t_out=h_t_out,
        m_s=1e2/3.6,
        m_t=68.9,
        redeclare VIP.Heat_transfer.Shell.Bell_Delaware   hT_shell(pitch_f = shell_and_tube.pitch_f, Dhyd_o = shell_and_tube.mytubes.Dhyd_o,
        N_tubes = shell_and_tube.N_tubes, d_s = shell_and_tube.d_s, d_b = shell_and_tube.mytubes.d_b, Dhyd = shell_and_tube.d_s,
        layout = shell_and_tube.layout, l_b = shell_and_tube.l_b, N_ss = 5,
        eta_wall=shell_and_tube.myshell.state[1].eta*ones(shell_and_tube.Ncell),
        N_passes = shell_and_tube.N_passes, ff =  1.0, b_cut = 0.25),
        redeclare VIP.Pressure_drops.Shell.Bell_Delaware   dp_shell(pitch_f = shell_and_tube.pitch_f, Dhyd_o = shell_and_tube.mytubes.Dhyd_o,
        N_tubes = shell_and_tube.N_tubes, d_s = shell_and_tube.d_s, d_b = shell_and_tube.mytubes.d_b, Dhyd = shell_and_tube.d_s,
        layout = shell_and_tube.layout, l_b = shell_and_tube.l_b, N_ss = 5,
        eta_wall=shell_and_tube.myshell.state[1].eta*ones(shell_and_tube.Ncell),
        N_passes = shell_and_tube.N_passes, ff =  1.0, b_cut = 0.25),
        redeclare function bundle_clearance =
          VIP.Useful_functions.Shell.SRFH,
        p_s_in=301435,
        p_s_out=301435,
        p_t_in=100000,
        p_t_out=100000)
        annotation (Placement(transformation(extent={{-80,-68},{68,58}})));

     annotation (Placement(transformation(extent={{-108,-74},{88,66}})));
    end Coulson_Bell_Delaware;


    model Aspen_Bell_Delaware
      "Verification with the results given by Aspen using the Bell Delaware method."

      replaceable package Medium_s = VIP.Media.Methanol_CoolProp "Medium model";
      replaceable package Medium_t = VIP.Media.Water_CoolProp "Medium model";
      parameter Modelica.SIunits.SpecificEnthalpy h_s_in= Medium_s.specificEnthalpy_pT(3.2e5,368.15)
        "Inlet specific enthalpy shell side";
      parameter Modelica.SIunits.SpecificEnthalpy h_s_out= Medium_s.specificEnthalpy_pT(3.2e5,313.15)
        "Outlet specific enthalpy shell side";
      parameter Modelica.SIunits.SpecificEnthalpy h_t_in = Medium_t.specificEnthalpy_pT(1e5,298.15)
        "Inlet specific enthalpy tube side";
      parameter Modelica.SIunits.SpecificEnthalpy h_t_out = Medium_t.specificEnthalpy_pT(1e5,313.15)
        "Outlet specific enthalpy tube side";

       Components.HEX.shell_and_tube  shell_and_tube(
         redeclare package Medium_s = Medium_s,
         redeclare package Medium_t = Medium_t,
         Dhyd=16e-3,
         thick=2e-3,
         l=4.65,
         pitch_f=1.25,
         lambda_wall=50,
         layout=1,
         N_passes=2,
         N_baffles=24,
         h_s_in=h_s_in,
         h_s_out=h_s_out,
         h_t_in=h_t_in,
         h_t_out=h_t_out,
         m_s=1e2/3.6,
         m_t=68.9,
         redeclare VIP.Heat_transfer.Shell.Bell_Delaware hT_shell(pitch_f = shell_and_tube.pitch_f, Dhyd_o = shell_and_tube.mytubes.Dhyd_o,
         N_tubes = shell_and_tube.N_tubes, d_s = shell_and_tube.d_s, d_b = shell_and_tube.mytubes.d_b, Dhyd = shell_and_tube.d_s,
         layout = shell_and_tube.layout, l_b = shell_and_tube.l_b, N_ss = 5,
         eta_wall=shell_and_tube.myshell.state[1].eta*ones(shell_and_tube.Ncell),
         N_passes = shell_and_tube.N_passes, ff =  1.0, b_cut = 0.1588),
         redeclare VIP.Pressure_drops.Shell.Bell_Delaware dp_shell(pitch_f = shell_and_tube.pitch_f, Dhyd_o = shell_and_tube.mytubes.Dhyd_o,
         N_tubes = shell_and_tube.N_tubes, d_s = shell_and_tube.d_s, d_b = shell_and_tube.mytubes.d_b, Dhyd = shell_and_tube.d_s,
         layout = shell_and_tube.layout, l_b = shell_and_tube.l_b, N_ss = 5,
         eta_wall=shell_and_tube.myshell.state[1].eta*ones(shell_and_tube.Ncell),
         N_passes = shell_and_tube.N_passes, ff =  1.0, b_cut = 0.1588),
         p_s_in=301435,
         p_s_out=301435,
         p_t_in=100000,
         p_t_out=100000)
        annotation (Placement(transformation(extent={{-80,-68},{68,58}})));

     annotation (Placement(transformation(extent={{-108,-74},{88,66}})));
    end Aspen_Bell_Delaware;


    model cell_method_VDI_example
      "Cell method verified with an example from VDI C1 pag. 48"
      import constants = Modelica.Constants;
      parameter Integer N_passes=2 "Number of cell passes";
      parameter Integer N_baffles=12 "Number of cell baffles";
      parameter Integer Ncell=N_passes*N_baffles "Number of cell elements";
      parameter Modelica.SIunits.Temp_C t_s_in = 100
        "Inlet temperature of the hot source";
      parameter Modelica.SIunits.Temp_C t_t_in = 20
        "Inlet temperature of the cold source";
      parameter Modelica.SIunits.CoefficientOfHeatTransfer kA_tot = 4749/Ncell;
      parameter Modelica.SIunits.HeatCapacity C_hot = 3.5e3
        "Capacity of the hot stream";
      parameter Modelica.SIunits.HeatCapacity C_cold = 3.5e3
        "Capacity of the cold stream";
    //  Modelica.SIunits.HeatFlowRate qdot[Ncell-1] "Heat rate";
      Modelica.SIunits.HeatFlowRate qdot_hot[Ncell] "Heat rate hot side";
      Modelica.SIunits.HeatFlowRate qdot_cold[Ncell] "Heat rate cold side";
      Modelica.SIunits.Temp_C t_hot[Ncell+1](start=linspace(t_s_in,40,Ncell+1))
        "Hot stream linear distribution with a guess";
      Modelica.SIunits.Temp_C t_cold[Ncell+1](start=linspace(t_t_in,40,Ncell+1))
        "Cold stream linear distribution with a guess";
      Modelica.SIunits.Temp_C t_av_hot[Ncell](start=linspace(t_s_in,40,Ncell))
        "Hot stream linear distribution with a guess";
      Modelica.SIunits.Temp_C t_av_cold[Ncell](start=linspace(t_t_in,40,Ncell))
        "Cold stream linear distribution with a guess";
       Real adim_hot[Ncell] "Adimensional temperatures of the hot stream";
       Real adim_cold[Ncell] "Adimensional temperatures of the cold stream";

    equation
        for j in 1:Ncell loop
          qdot_hot[j]  = C_hot*(t_hot[j] - t_hot[j+1]);
          t_av_hot[j]  = 0.5*(t_hot[j+1] + t_hot[j]);
          qdot_cold[j] = C_cold*(t_cold[j+1] - t_cold[j]);
          t_av_cold[j] = 0.5*(t_cold[j+1] + t_cold[j]);
          adim_hot[j]  = (t_hot[j] - t_t_in)/(t_s_in - t_t_in);
          adim_cold[j] = (t_cold[j] - t_t_in)/(t_s_in - t_t_in);
        end for;

         for j in 1:N_passes loop
             for i in 1+(j-1)*N_baffles:2:j*N_baffles loop
                if mod(j,2) <> 0 then
                  if mod(N_baffles,2) <> 0 then
                    qdot_cold[i] = qdot_hot[Ncell - N_passes*(i-(1+(j-1)*N_baffles)) - (j-1)];
                    qdot_cold[i] = kA_tot*(t_av_hot[Ncell - N_passes*(i-(1+(j-1)*N_baffles)) - (j-1)] - t_av_cold[i]);
                  else
                    qdot_cold[i] = qdot_hot[Ncell - N_passes*(i-(1+(j-1)*N_baffles)) - (N_passes - j)];
                    qdot_cold[i] = kA_tot*(t_av_hot[Ncell - N_passes*(i-(1+(j-1)*N_baffles)) - (N_passes - j)] - t_av_cold[i]);
                  end if;
                else
                    qdot_cold[i] = qdot_hot[(N_passes + 1 - j) + N_passes*(i-(1+(j-1)*N_baffles))];
                    qdot_cold[i] = kA_tot*(t_av_hot[(N_passes + 1 - j) + N_passes*(i-(1+(j-1)*N_baffles))] - t_av_cold[i]);
                end if;
             end for;
         end for;

         for j in 1:N_passes loop
             for i in 2+(j-1)*N_baffles:2:j*N_baffles loop
                if mod(j,2) <> 0 then
                  if mod(N_baffles,2) <> 0 then
                    qdot_cold[i] = qdot_hot[Ncell - (2*N_passes - 1) + (j-1) - N_passes*(i-(2+(j-1)*N_baffles))];
                    qdot_cold[i] = kA_tot*(t_av_hot[Ncell - (2*N_passes - 1) + (j-1) - N_passes*(i-(2+(j-1)*N_baffles))] - t_av_cold[i]);
                  else
                    qdot_cold[i] = qdot_hot[Ncell - N_passes + (j-1) - N_passes*(i-(2+(j-1)*N_baffles))];
                    qdot_cold[i] = kA_tot*(t_av_hot[Ncell - N_passes + (j-1) - N_passes*(i-(2+(j-1)*N_baffles))] - t_av_cold[i]);
                  end if;
                else
                    qdot_cold[i] = qdot_hot[N_passes + j + N_passes*(i-(2+(j-1)*N_baffles))];
                    qdot_cold[i] = kA_tot*(t_av_hot[N_passes + j + N_passes*(i-(2+(j-1)*N_baffles))] - t_av_cold[i]);
                end if;
             end for;
         end for;

         //Boundary conditions
         t_hot[1]  = t_s_in;
         t_cold[1] = t_t_in;

      annotation (experiment(__Dymola_NumberOfIntervals=1, Tolerance=1e-009),
          __Dymola_experimentSetupOutput,
        Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,100}}),
            graphics));
    end cell_method_VDI_example;

    model Grossman_Bell_Delaware
      "Verification with the results given by Grossman using the Bell Delaware method."

      replaceable package Medium_s = VIP.Media.Methanol_CoolProp "Medium model";
      replaceable package Medium_t = VIP.Media.Water_CoolProp "Medium model";
      parameter Modelica.SIunits.SpecificEnthalpy h_s_in= Medium_s.specificEnthalpy_pT(3.2e5,368.15)
        "Inlet specific enthalpy shell side";
      parameter Modelica.SIunits.SpecificEnthalpy h_s_out= Medium_s.specificEnthalpy_pT(3.2e5,313.15)
        "Outlet specific enthalpy shell side";
      parameter Modelica.SIunits.SpecificEnthalpy h_t_in = Medium_t.specificEnthalpy_pT(1e5,298.15)
        "Inlet specific enthalpy tube side";
      parameter Modelica.SIunits.SpecificEnthalpy h_t_out = Medium_t.specificEnthalpy_pT(1e5,313.15)
        "Outlet specific enthalpy tube side";

       Components.HEX.shell_and_tube  shell_and_tube(
         redeclare package Medium_s = Medium_s,
         redeclare package Medium_t = Medium_t,
         Dhyd=12.6e-3,
         thick=3.3e-3,
         l=4.88,
         pitch_f=1.25,
         lambda_wall=50,
         layout=2,
         N_passes=2,
         N_baffles=8,
         h_s_in=h_s_in,
         h_s_out=h_s_out,
         h_t_in=h_t_in,
         h_t_out=h_t_out,
         m_s=1e2/3.6,
         m_t=68.9,
         redeclare VIP.Heat_transfer.Shell.Bell_Delaware   hT_shell(pitch_f = shell_and_tube.pitch_f, Dhyd_o = shell_and_tube.mytubes.Dhyd_o,
         N_tubes = shell_and_tube.N_tubes, d_s = shell_and_tube.d_s, d_b = shell_and_tube.mytubes.d_b, Dhyd = shell_and_tube.d_s,
         layout = shell_and_tube.layout, l_b = shell_and_tube.l_b, N_ss = 6,
         eta_wall=shell_and_tube.myshell.state[1].eta*ones(shell_and_tube.Ncell),
         N_passes = shell_and_tube.N_passes, ff =  1.0, b_cut = 0.25),
         redeclare VIP.Pressure_drops.Shell.Bell_Delaware   dp_shell(pitch_f = shell_and_tube.pitch_f, Dhyd_o = shell_and_tube.mytubes.Dhyd_o,
         N_tubes = shell_and_tube.N_tubes, d_s = shell_and_tube.d_s, d_b = shell_and_tube.mytubes.d_b, Dhyd = shell_and_tube.d_s,
         layout = shell_and_tube.layout, l_b = shell_and_tube.l_b, N_ss = 6,
         eta_wall=shell_and_tube.myshell.state[1].eta*ones(shell_and_tube.Ncell),
         N_passes = shell_and_tube.N_passes, ff =  1.0, b_cut = 0.25),
         p_s_in=301435,
         p_s_out=301435,
         p_t_in=100000,
         p_t_out=100000)
        annotation (Placement(transformation(extent={{-80,-68},{68,58}})));

     annotation (Placement(transformation(extent={{-108,-74},{88,66}})));
    end Grossman_Bell_Delaware;
  end Tests;

  annotation (uses(Modelica(version="3.2")), Icon(coordinateSystem(
          preserveAspectRatio=false, extent={{-100,-100},{100,100}}), graphics));
end VIP;
