within ;
package VIP "I am a package for the Virtual Prototyping Environment"
  import Modelica.Math.*;
  import Modelica.SIunits.*;
  import Modelica.Constants.*;

  package Icons "Icons for the virtual prototpying environment"
    partial model shell_tube "Shell and tube icon"

      annotation (Icon(graphics={
            Rectangle(
              extent={{-80,60},{90,-40}},
              lineColor={0,0,0}),
            Rectangle(
              extent={{-60,-40},{-40,-60}},
              lineColor={0,0,0}),
            Rectangle(
              extent={{60,80},{80,60}},
              lineColor={0,0,0}),
            Line(
              points={{70,90},{70,66}},
              color={255,0,0},
              smooth=Smooth.None),
            Line(
              points={{2,-2},{-2,2}},
              color={255,0,0},
              smooth=Smooth.None,
              origin={72,68},
              rotation=90),
            Line(
              points={{2,2},{-2,-2}},
              color={255,0,0},
              smooth=Smooth.None,
              origin={68,68},
              rotation=90),
            Line(
              points={{-50,-42},{-50,-66}},
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
            Rectangle(extent={{-90,-20},{60,-20}}, lineColor={85,170,255}),
            Rectangle(extent={{-90,40},{60,40}}, lineColor={85,170,255}),
            Line(
              points={{60,40},{60,-20}},
              color={85,170,255},
              smooth=Smooth.None),
            Line(
              points={{-20,-20},{-24,-16}},
              color={85,170,255},
              smooth=Smooth.None),
            Line(
              points={{-20,-20},{-24,-24}},
              color={85,170,255},
              smooth=Smooth.None),
            Line(
              points={{-20,40},{-16,44}},
              color={85,170,255},
              smooth=Smooth.None),
            Line(
              points={{-20,40},{-16,36}},
              color={85,170,255},
              smooth=Smooth.None)}));
    end shell_tube;

    partial package tube "tube icon"

      annotation (Icon(graphics={
            Rectangle(extent={{-80,-20},{60,-20}}, lineColor={85,170,255}),
            Rectangle(extent={{-80,40},{60,40}}, lineColor={85,170,255}),
            Line(
              points={{60,40},{60,-20}},
              color={85,170,255},
              smooth=Smooth.None),
            Line(
              points={{-20,-20},{-24,-16}},
              color={85,170,255},
              smooth=Smooth.None),
            Line(
              points={{-20,-20},{-24,-24}},
              color={85,170,255},
              smooth=Smooth.None),
            Line(
              points={{-20,40},{-16,44}},
              color={85,170,255},
              smooth=Smooth.None),
            Line(
              points={{-20,40},{-16,36}},
              color={85,170,255},
              smooth=Smooth.None)}));
    end tube;

    partial package shell "Shell icon"

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
    end shell;

    partial package clearance "clearance icon"

      annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},
                {100,100}}),
                       graphics={Ellipse(extent={{-20,20},{20,-20}}, lineColor={0,0,
                  0}), Ellipse(extent={{60,-60},{-60,60}}, lineColor={0,0,0})}));
    end clearance;

    partial package cost "cost icon"

      annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},
                {100,100}}),
                       graphics={Text(
              extent={{-126,-24},{129,-219}},
              lineColor={0,0,0},
              textString="$
")}));
    end cost;

    partial package HEX "heat exchanger icon"

      annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,
                -100},{100,100}}),
                       graphics={
            Ellipse(extent={{-80,80},{80,-80}}, lineColor={0,0,0}),
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
  end Icons;

  package Heat_transfer "A package containing heat transfer correlations"
    package Tubes "heat transfer correlations in tubes"
      extends VIP.Icons.tube;
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
              Re[i]             = Miscellanea.Pure_numbers.Reynolds(     u[i], state[i].d,  state[i].eta, Dhyd)
            "Reynolds number tubes";
              Pr[i]             = Miscellanea.Pure_numbers.Prandtl(     state[i].cp, state[i].eta, state[i].lambda)
            "Prandtl number tubes";
              assert(Re[i] > 1e4, "Reynolds number is lower than 1e4 to use Dittus and Boelter", AssertionLevel.warning);
              Nu[i]             =  2.3e-2*Re[i]^0.8*Pr[i]^alfa
            "Nusselt number tubes";
              ht[i]             = Miscellanea.Pure_numbers.Nusselt(
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
              Re[i]             = Miscellanea.Pure_numbers.Reynolds(     u[i], state[i].d,  state[i].eta, Dhyd)
            "Reynolds number tubes";
              Pr[i]             = Miscellanea.Pure_numbers.Prandtl(     state[i].cp, state[i].eta, state[i].lambda)
            "Prandtl number tubes";
              assert(Re[i] > 1e4, "Reynolds number is lower than 1e4 to use Sieder and Tate", AssertionLevel.warning);
              assert(Pr[i] > 0.6,  "Prandtl number is lower than 0.6 to be used in Sieder and Tate", AssertionLevel.warning);
              Nu[i]             = 2.3e-2*Re[i]^0.8*Pr[i]^(1.0/3)*(eta_wall[i]/state[i].eta)^0.14
            "Nusselt number tubes";
              ht[i]             = Miscellanea.Pure_numbers.Nusselt(
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
              Re[i]             = Miscellanea.Pure_numbers.Reynolds(     u[i], state[i].d,  state[i].eta, Dhyd)
            "Reynolds number tubes";
              Pr[i]             = Miscellanea.Pure_numbers.Prandtl(     state[i].cp, state[i].eta, state[i].lambda)
            "Prandtl number tubes";
              csi[i]            = 1/(0.78*log(Re[i]) - 1.5)^2 "Friction factor";
              Nu[i]             =  ((csi[i]/8)*Re[i]*Pr[i])/(1 + 12.7*sqrt(csi[i]/8)*(Pr[i]^(2.0/3) - 1))*(1 +  (Dhyd/l)^(2.0/3))
            "Nusselt number tubes";
              ht[i]             = Miscellanea.Pure_numbers.Nusselt(
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
              ht[i]             = Miscellanea.Pure_numbers.Nusselt(
                                                       Nu[i], state[i].lambda, Dhyd)
            "Heat transfer coefficient tube side";
            end for;
      end EagleFerguson;

      class basic_ht "Basic heat transfer correlation"
          replaceable package Medium = VIP.Media.OneRandomOrganicFluid
          "Medium model";
          parameter Integer Ncell(start=3) "Number of cell elements";
          input Modelica.SIunits.Area Aflow
          "Cross-sectional area (single tube)";
          input Modelica.SIunits.MassFlowRate mdot "Mass flow rate";
          input Medium.ThermodynamicState state[Ncell];
      end basic_ht;
    end Tubes;

    package Shell "heat transfer correlations in shells"
      extends VIP.Icons.shell;

      class Kern "Kern correlation for shell side"
          extends VIP.Heat_transfer.Shell.basic_ht;
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
            Re[i]             = Miscellanea.Pure_numbers.Reynolds(     u[i], state[i].d,  state[i].eta, d_s_eq)
            "Reynolds number tubes";
            Pr[i]             = Miscellanea.Pure_numbers.Prandtl(     state[i].cp, state[i].eta, state[i].lambda)
            "Prandtl number tubes";
            Nu[i]             = 10^log10(2e-3/1e5^(-0.25*log10(1.98/2e-2)))*Re[i]^(-0.25*log10(1.98/2e-2))*Re[i]*Pr[i]^(1/3)
            "Nusselt number tubes";
            ht[i]             = Miscellanea.Pure_numbers.Nusselt(     Nu[i], state[i].lambda, d_s_eq)
            "Heat transfer coefficient tube side";
            end for;
      end Kern;

      class Bell_Delaware "Bell Delaware correlation for shell side"
          extends VIP.Heat_transfer.Shell.basic_ht;
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
          Re[i]   = Miscellanea.Pure_numbers.Reynolds(     u[i], state[i].d,  state[i].eta, Dhyd_o);
          Pr[i]   = Miscellanea.Pure_numbers.Prandtl(     state[i].cp, state[i].eta, state[i].lambda);

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

      class basic_ht "Basic heat transfer correlation"
          replaceable package Medium = VIP.Media.OneRandomOrganicFluid
          "Medium model";
          parameter Integer Ncell(start=3) "Number of cell elements";
          input Modelica.SIunits.Area Aflow
          "Cross-sectional area (single tube)";
          input Modelica.SIunits.MassFlowRate mdot "Mass flow rate";
          input Medium.ThermodynamicState state[Ncell];
      end basic_ht;
    end Shell;
  end Heat_transfer;

  package Pressure_drops "A package containing pressure drops correlations"
    package Tubes "heat transfer correlations in tubes"
        extends VIP.Icons.tube;
      class Frank "Frank correlation for tubes"
          extends VIP.Pressure_drops.Tubes.basic_dp;
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
          Re[i]    = Miscellanea.Pure_numbers.Reynolds(     u[i], state[i].d,  state[i].eta, Dhyd)
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
          input Modelica.SIunits.Area Aflow
          "Cross-sectional area (single tube)";
          input Modelica.SIunits.MassFlowRate mdot "Mass flow rate";
          input Medium.ThermodynamicState state[Ncell];
      end basic_dp;
    end Tubes;

    package Shell "heat transfer correlations in shells"
      extends VIP.Icons.shell;
      class Kern "Kern correlation for shell"
          extends VIP.Pressure_drops.Shell.basic_dp;
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
          Re[i]    = Miscellanea.Pure_numbers.Reynolds(     u[i], state[i].d,  state[i].eta, Dhyd)
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
          extends VIP.Pressure_drops.Shell.basic_dp;
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
          Re[i]    = Miscellanea.Pure_numbers.Reynolds(     u[i], state[i].d,  state[i].eta, Dhyd)
            "Reynolds number tubes";
          assert(Re[i] > 4e2, "Reynolds number is lower than 1e2 to use Wills_Johnston", AssertionLevel.warning);
          csi[i] = 1.7789*Re[i]^(-0.195868);
          dp[i]  =  0.5*csi[i]*d_s/Dhyd*l/l_b*(eta_wall[i]/state[i].eta)^0.14*state[i].d*u[i]^2
            "Local pressure drops";
        end for;

        dp_tot = sum(dp) "Total pressure drops";

      end Wills_Johnston;

      class Bell_Delaware "Bell Delaware correlation for shell side"
          extends VIP.Pressure_drops.Shell.basic_dp;
          input Modelica.SIunits.Length Dhyd "Hydraulic Diameter (single tube)";
          input Modelica.SIunits.Length Dhyd_o "Outer Diameter (single tube)";
          parameter Real b_cut "Baffle cut";
          final parameter Real ttb = 0.8e-3
          "tube to baffle clearance (from Standards)";
          final parameter Real bts = 4e-3
          "tube to baffle clearance (from Standards)";
          parameter Real pitch_f
          "Tube pitch as a fraction of the outer tube diameter";
          parameter Real N_ss
          "The number of sealing strips (pairs) in one baffle spacing";
          parameter Integer layout "Tube layout 1 = triangular, 2 = squared";
          parameter Integer N_passes "Number of tube passes";
          parameter Real ff "Fouling factor (=1 if clean)";
          parameter Integer N_baffles "Number of baffles";
          parameter Integer N_baffles_d
          "The number of discretization volumes is always N_baffles + 1";
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
          Modelica.SIunits.ReynoldsNumber Re[Ncell]
          "Reynolds number shell (average)";
          Modelica.SIunits.Velocity u[Ncell](start=ones(Ncell))
          "Velocities inside the shell";
          Modelica.SIunits.Velocity u_w[Ncell](start=ones(Ncell))
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
        R_l      = exp(-1.33*(1 + r_s)*r_lm^(-0.15*(1 + r_s) + 0.8));

        S_b      = l_b*(d_s - d_b);
        F_sbp    = S_b/S_m;
        R_b      = exp(-3.7*F_sbp*(1 - (2/N_ss)^(1/3)));

        S_w      = S_wg - S_wt;
        N_tcw    = 0.8*(d_s*b_cut - 0.5*(d_s - (d_b - Dhyd_o)))/tpv;

        for i in 1:Ncell loop

          u[i]    = mdot/state[i].d/S_m;
          u_w[i]  = mdot/state[i].d/S_w;
          Re[i]   = Miscellanea.Pure_numbers.Reynolds(     u[i], state[i].d,  state[i].eta, Dhyd_o);

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
            dp_w[i]  = (2 + 0.6*N_tcw)*mdot^2/(2*S_w*S_m*state[i].d)*(N_baffles/(N_passes*N_baffles_d));
          end if;

          if (i <= N_passes or i > Ncell - N_passes) then
            dp_c[i] = 0;
            dp_e[i] = 2*(N_tcc + N_tcw)*csi[i]*(state[i].d*u[i]^2/2)*(N_baffles/(N_passes*N_baffles_d));
          else
            dp_e[i] = 0;
            dp_c[i] = ff*2*N_tcc*csi[i]*(state[i].d*u[i]^2/2)*(N_baffles/(N_passes*N_baffles_d));
          end if;

        end for;

        dp_w_tot = sum(dp_w)*R_l;
        dp_c_tot = sum(dp_c)*R_b*R_l;
        dp_e_tot = sum(dp_e)*R_b;

        dp_tot = dp_w_tot + dp_c_tot + dp_e_tot "Total pressure drops";

      end Bell_Delaware;

      class basic_dp "Basic pressure drop correlation"
          replaceable package Medium = VIP.Media.OneRandomOrganicFluid
          "Medium model";
          parameter Integer Ncell(start=3) "Number of cell elements";
          input Modelica.SIunits.Area Aflow
          "Cross-sectional area (single tube)";
          input Modelica.SIunits.MassFlowRate mdot "Mass flow rate";
          input Medium.ThermodynamicState state[Ncell];
      end basic_dp;
    end Shell;
  end Pressure_drops;

  package Objects "Package containing all the objects of the VIP"

    class tubes "I am a tube and I contain all the my relevant informations"
        replaceable package Medium = VIP.Media.OneRandomOrganicFluid
        "Medium model";
        parameter Integer Ncell(start=3) "Number of cell elements";
        parameter Modelica.SIunits.Length thick "Thickness of the tube";
        parameter Modelica.SIunits.ThermalConductivity lambda
        "Thermal conductivity of the tube wall";
        parameter Modelica.SIunits.Density rho "Density of the tube wall";
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
        Modelica.SIunits.Mass W_dry "Dry weight of one tube";
        Modelica.SIunits.Mass W_wet[Ncell] "Dry weight of one tube";

    equation
        //Trivial calculations
        Dhyd_o = Dhyd + 2*thick;
        G_wall = 2*lambda/log(Dhyd_o/Dhyd);

        for i in 1:Ncell loop
          qdot[i]           = pin*mdot*(h[i] - h[i + 1]);
          state[i]          = Medium.setState_ph(p_in, (h[i] + h[i + 1])/2);
        end for;
      annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},
                {100,100}}), graphics={
            Rectangle(extent={{-70,0},{56,0}},     lineColor={85,170,255}),
            Rectangle(extent={{-70,40},{56,40}}, lineColor={85,170,255}),
            Line(
              points={{56,40},{56,0}},
              color={85,170,255},
              smooth=Smooth.None),
            Line(
              points={{-10,0},{-14,4}},
              color={85,170,255},
              smooth=Smooth.None),
            Line(
              points={{-10,0},{-14,-4}},
              color={85,170,255},
              smooth=Smooth.None),
            Line(
              points={{-10,40},{-6,44}},
              color={85,170,255},
              smooth=Smooth.None),
            Line(
              points={{-10,40},{-6,36}},
              color={85,170,255},
              smooth=Smooth.None)}));
    end tubes;

    class tube_bundle
      "I am a tube bundle and I contain all the my relevant information"
        extends VIP.Objects.tubes;
        parameter Integer N_passes "Number of tube passes";
        parameter Integer layout "Tube layout 1 = triangular, 2 = squared";
        parameter Real pitch_f
        "Tube pitch as a fraction of the outer tube diameter";
        input Real N_tubes "Number of tubes in the bundle";
    //     input Modelica.SIunits.Length clearance "Bundle - shell clearance";
        Modelica.SIunits.Length d_b "Bundle Diameter";
    equation

        d_b       = Miscellanea.bundle_diameter(     N_tubes, N_passes, Dhyd_o, layout)
        "calculate the tube bundle";

    end tube_bundle;
  end Objects;

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

    package Shell_clearance
      "Functions to evaluate the shell clearance based on the diameter of the tube bundle "
      extends VIP.Icons.clearance;
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
    end Shell_clearance;

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
      replaceable package Medium = VIP.Media.OneRandomOrganicFluid
        "Medium model";
      input Modelica.SIunits.Velocity umax "Maximum velocity";
      input Modelica.SIunits.Velocity umin "Minimum velocity";
      parameter String geometry = "tube"
        "Type of geometry. E.g. shell, tubes, etc...";
      parameter Modelica.SIunits.AbsolutePressure op_p
        "Operating pressure. If the fluid is a vapour then the limit depends on the pressure";
      final parameter Modelica.SIunits.Temperature T_sat = Medium.saturationTemperature(op_p)
        "Saturation temperature";
      input Modelica.SIunits.Temperature T "Fluid temperature";

    equation
      if (T < T_sat) then
        if (geometry == "tube") then
          assert(umax < 4*2, "Velocity in the tubes twice higher than 4 m/s", AssertionLevel.warning);
          assert(umin > 1*0.5,"Velocity in the tubes twice lower than 0.8 m/s", AssertionLevel.warning);
        elseif (geometry == "shell") then
          assert(umax < 1.2*2, "Velocity in the shell twice higher than 1.2 m/s", AssertionLevel.warning);
          assert(umin > 0.3*0.5,"Velocity in the shell twice lower than 0.3 m/s", AssertionLevel.warning);
        end if;
      elseif (T > T_sat) then
        if (op_p > 1e6) then
          assert(umax < 70*2, "Velocity in the tubes twice higher than 70 m/s", AssertionLevel.warning);
          assert(umin > 50*0.5,"Velocity in the tubes twice lower than 50 m/s", AssertionLevel.warning);
        elseif (op_p < 1e4) then
          assert(umax < 10*2, "Velocity in the shell twice higher than 10 m/s", AssertionLevel.warning);
          assert(umin > 5*0.5,"Velocity in the shell twice lower than 5 m/s", AssertionLevel.warning);
        else
          assert(umax < 30*2, "Velocity in the shell twice higher than 30 m/s", AssertionLevel.warning);
          assert(umin > 10*0.5,"Velocity in the shell twice lower than 10 m/s", AssertionLevel.warning);
        end if;
      end if;
    end check_velocity;


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

    package Cost "Functions to evaluate the cost of heat exchangers "
      extends VIP.Icons.cost;
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
      annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
                -100},{100,100}}), graphics), Icon(graphics));
    end Cost;
  end Miscellanea;

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

    package Cyclopentane_CoolProp "CoolProp model of Cyclopentane"
      extends ExternalMedia.Media.CoolPropMedium(
        mediumName = "Cyclopentane",
        substanceNames = {"Cyclopentane|enable_TTSE=1|enable_EXTTP=1"},
        ThermoStates = Modelica.Media.Interfaces.PartialMedium.Choices.IndependentVariables.ph);
    end Cyclopentane_CoolProp;
  end Media;

  package Materials "Package containing the properties of different materials"
    extends Modelica.Icons.MaterialPropertiesPackage;
    class S_AISI_1010 "S_AISI_1010"
        extends VIP.Materials.void_material( materialName = "S_AISI_1010",
        lambda = 65.2,
        rho = 7850);
    end S_AISI_1010;

    class SS_AISI_410 "SS_AISI_410"
        extends VIP.Materials.void_material(
                                       materialName = "SS_AISI_410",
        lambda = 24.9, rho = 7740);
    end SS_AISI_410;

    class S_AISI_1040 "S_AISI_1040"
        extends VIP.Materials.void_material(materialName = "S_AISI_1040",
        lambda = 50.7,
        rho = 7845);
    end S_AISI_1040;

    class void_material "Void material"
        extends Modelica.Icons.MaterialProperty;
        parameter String materialName = "Random";
        parameter Modelica.SIunits.ThermalConductivity lambda = 50;
        parameter Modelica.SIunits.Density rho = 8000;
    end void_material;
  end Materials;

  package Components "Library with the design of the components "
    package HEX "Heat exchangers"
      extends VIP.Icons.HEX;
      model shell_and_tube
        "Shell and tube heat exchanger where the hot fluid flows on the shell and enters from the top. The cold fluid enters at the bottom."
        extends VIP.Icons.shell_tube;

        //THE WORKING FLUIDS
        replaceable package Medium_s = VIP.Media.OneRandomOrganicFluid constrainedby
          Modelica.Media.Interfaces.PartialMedium "Medium model shell" annotation(choicesAllMatching = true);
        replaceable package Medium_t = VIP.Media.OneRandomOrganicFluid constrainedby
          Modelica.Media.Interfaces.PartialMedium "Medium model tubes" annotation(choicesAllMatching = true);
        replaceable VIP.Materials.void_material
                                           Material_t "Material model shell"        annotation(choicesAllMatching = true);
        replaceable VIP.Materials.void_material
                                           Material_s "Material model shell"         annotation(choicesAllMatching = true);

        //GEOMETRY OF THE HEAT EXCHANGER
         parameter Modelica.SIunits.Length Dhyd
          "Hydraulic Diameter (single tube)";
         parameter Modelica.SIunits.Length thick_t "Thickness of the tube";
         parameter Modelica.SIunits.Length thick_s "Thickness of the tube";
         parameter Modelica.SIunits.Length l "Lenght (single tube)";
         parameter Real pitch_f
          "Tube pitch as a fraction of the outer tube diameter";
         parameter Integer layout "Tube layout 1 = triangular, 2 = squared";
         parameter Integer N_passes = 2 "Number of tube passes";
         parameter Integer N_baffles = 4 "Number of baffles";
         parameter Integer N_heads
          "Number of heads (maximum = 2 front and rear";
         parameter Real pin_s
          "Pin for the heat flow in the shell -1 cold fluid else hot";
         parameter Real pin_t
          "Pin for the heat flow in the tubes -1 cold fluid else hot";
         parameter Modelica.SIunits.CoefficientOfHeatTransfer U_guess = 600
          "Guess value for the global heat transfer coefficient";
         parameter Integer N_baffles_d = N_baffles + 1
          "The number of discretization volumes is always N_baffles + 1";
         final parameter Integer Ncell = N_baffles_d*N_passes
          "Number of cell elements";
         final parameter Modelica.SIunits.HeatFlowRate qdot = m_s*abs(h_s_in - h_s_out)
          "Heat flow rate";
         parameter Modelica.SIunits.Temp_C  DTML = VIP.Miscellanea.log_mean_delta_T(t_s_in, t_s_out, t_t_in, t_t_out)
          "Logarithmic mean temperature difference";
         //Boundary conditions at inlet and outlet
         parameter Modelica.SIunits.MassFlowRate m_s "Mass flow shell side";
         parameter Modelica.SIunits.SpecificEnthalpy h_s_in
          "Inlet specific enthalpy shell side";
         parameter Modelica.SIunits.SpecificEnthalpy h_s_out
          "Outlet specific enthalpy shell side";
         parameter Modelica.SIunits.AbsolutePressure p_s_in
          "Inlet pressure shell side";
         parameter Modelica.SIunits.AbsolutePressure p_s_out
          "Outlet pressure shell side";
         parameter Modelica.SIunits.MassFlowRate m_t "Mass flow tube side";
         parameter Modelica.SIunits.SpecificEnthalpy h_t_in
          "Inlet specific enthalpy tube side";
         parameter Modelica.SIunits.SpecificEnthalpy h_t_out
          "Outlet specific enthalpy tube side";
         parameter Modelica.SIunits.AbsolutePressure p_t_in
          "Inlet pressure tube side";
         parameter Modelica.SIunits.AbsolutePressure p_t_out
          "Outlet pressure tube side";
         final parameter Modelica.SIunits.Temperature  t_s_in = Medium_s.temperature_ph(p_s_in, h_s_in)
          "Inlet temperature shell side";
         final parameter Modelica.SIunits.Temperature  t_s_out = Medium_s.temperature_ph(p_s_out, h_s_out)
          "Outlet temperature shell side";
         final parameter Modelica.SIunits.Temperature  t_t_in = Medium_t.temperature_ph(p_t_in, h_t_in)
          "Inlet temperature tube side";
         final parameter Modelica.SIunits.Temperature  t_t_out = Medium_t.temperature_ph(p_t_out, h_t_out)
          "Outlet temperature tube side";
         parameter Modelica.SIunits.CoefficientOfHeatTransfer ht_t_f1 = 3e3
          "fouling heat transfer coefficient (tube side)";
         parameter Modelica.SIunits.CoefficientOfHeatTransfer ht_s_f1 = 5e3
          "fouling heat transfer coefficient (shell side)";
         Modelica.SIunits.Length d_s "Shell diameter";
         Modelica.SIunits.Length l_b "Baffle lenght";
         Real N_tubes( start = qdot/(DTML*U_guess)/(pi*l*(Dhyd + 2*thick_t)))
          "Number of tubes in the bundle";
         Real N_t_p_p "Number of tubes per pass";
         Real bs_f "Baffle spacing as a fraction of the shell diameter";
         Modelica.SIunits.CoefficientOfHeatTransfer U
          "Global heat transfer coefficient";
         Modelica.SIunits.Area A "Heat transfer area";
         Modelica.SIunits.Temp_C  DTML_tilde
          "Logarithmic mean temperature difference corrected";
         Modelica.SIunits.Mass W_dry "Dry weight of the heat exchanger";
         Modelica.SIunits.Mass W_fluids "Weight of the fluids";
         Modelica.SIunits.Mass W_wet "Wet weight of the heat exchanger";
         Real PEC "Purchases equipment cost";

        //Heat transfer and pressure drop correlations
        replaceable VIP.Heat_transfer.Tubes.Sieder_Tate hT_tube(Dhyd=Dhyd, eta_wall=bundle.state[1].eta*ones(Ncell))
           constrainedby VIP.Heat_transfer.Tubes.basic_ht(Medium = Medium_t, Ncell = Ncell, state=bundle.state, mdot=m_t/N_t_p_p, Aflow=bundle.Aflow) annotation(choicesAllMatching = true);
        replaceable VIP.Pressure_drops.Tubes.Frank dp_tube(Dhyd=Dhyd,l=l/N_baffles_d, heads = 2.5/Ncell)
          constrainedby VIP.Pressure_drops.Tubes.basic_dp(Medium = Medium_t, Ncell = Ncell, state=bundle.state, mdot=m_t/N_t_p_p, Aflow=bundle.Aflow) annotation(choicesAllMatching = true);

        replaceable VIP.Heat_transfer.Shell.Kern hT_shell(Dhyd_o=Dhyd+2*thick_t, layout=layout, pitch_f = pitch_f)
          constrainedby VIP.Heat_transfer.Shell.basic_ht(Medium = Medium_s, Ncell = Ncell, state=shell.state, mdot=m_s, Aflow=shell.Aflow) annotation(choicesAllMatching = true);
        replaceable VIP.Pressure_drops.Shell.Wills_Johnston dp_shell(l=l/N_baffles_d/N_passes, d_s=d_s, Dhyd=hT_shell.d_s_eq, l_b = l_b,
                      eta_wall=shell.state[1].eta*ones(Ncell))
         constrainedby VIP.Pressure_drops.Shell.basic_dp(Medium = Medium_s, Ncell = Ncell, state=shell.state, mdot=m_s, Aflow=shell.Aflow) annotation(choicesAllMatching = true);

        //Defining the model for the bundle clearance
        replaceable function bundle_clearance =
            VIP.Miscellanea.Shell_clearance.base_clearance
                                              annotation(choicesAllMatching = true);

        //Defining the model for the cost
        replaceable function cost =
            Miscellanea.Cost.base_cost        annotation(choicesAllMatching = true);

        //Definiing the tubes and the shell
        VIP.Objects.tube_bundle
                    bundle(redeclare package Medium = Medium_t,
                    h_in = h_t_in,
                    h_out = h_t_out,
                    Ncell = Ncell, Aflow = 0.25*pi*Dhyd^2, mdot = m_t,
                    Dhyd = Dhyd, thick = thick_t, lambda = Material_t.lambda,
                    rho = Material_t.rho, N_tubes = N_tubes,
                    mdot_pt = m_t/N_t_p_p,N_passes = N_passes, layout = layout,
                    pin = pin_t, pitch_f = pitch_f);

        VIP.Objects.tubes
                    shell(redeclare package Medium = Medium_s,
                    h_in = h_s_in,
                    h_out = h_s_out,
                    Ncell = Ncell, Aflow = (1 - 1/pitch_f)*d_s*l_b, mdot = m_s,
                    mdot_pt = m_s, Dhyd = d_s, thick = thick_s,
                    lambda = Material_s.lambda, rho = Material_s.rho, pin = pin_s);

        VIP.Miscellanea.check_velocity check_shell(redeclare package Medium = Medium_s,
                    T = shell.state[1].T, umin = min(hT_shell.u),
                    umax = max(hT_shell.u), geometry = "shell", op_p = p_s_in);

        VIP.Miscellanea.check_velocity check_tube(redeclare package Medium = Medium_t,
                    T = bundle.state[1].T, umin = min(hT_tube.u),
                    umax = max(hT_tube.u), geometry = "tube", op_p = p_t_in);

      protected
         parameter Modelica.SIunits.CoefficientOfHeatTransfer ht_t_f[Ncell] = ones(Ncell)*ht_t_f1
          "fouling heat transfer coefficient (tube side)";
         parameter Modelica.SIunits.CoefficientOfHeatTransfer ht_s_f[Ncell] = ones(Ncell)*ht_s_f1
          "fouling heat transfer coefficient (shell side)";
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
          bundle.p_in         = p_t_in;
          shell.p_in         = p_s_in;
          bundle.h[1]         = h_t_in;
          shell.h[1]         = h_s_in;
          shell.h[Ncell + 1] = h_s_out;

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
            G_wall[i]       = bundle.At*bundle.G_wall/bundle.Dhyd;
            htA_tubes[i]    = bundle.At/(1/hT_tube.ht[i] + 1/ht_t_f[i]);
            htA_shell[k[i]] = shell.At/(1/hT_shell.ht[k[i]] + 1/ht_s_f[k[i]]);
            kA_tot[i]       = 1/(1/htA_tubes[i] + 1/G_wall[i] + 1/htA_shell[k[i]]);
            bundle.qdot[i]  = shell.qdot[k[i]];
            bundle.qdot[i]  = kA_tot[i]*(pin_s*shell.state[k[i]].T + pin_t*bundle.state[i].T);

            //Fluid weight calculation
            bundle.W_wet[i] = 0.25*pi*Dhyd^2*l/N_baffles_d/N_passes*N_tubes*bundle.state[i].d;
            shell.W_wet[i]  = 0.25*pi*(shell.Dhyd^2 - N_tubes*bundle.Dhyd_o^2)*l/N_baffles_d/N_passes*shell.state[i].d;

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
          N_t_p_p        = N_tubes/N_passes "Number of tubes per pass";
          shell.At     = N_t_p_p*pi*l/N_baffles_d*bundle.Dhyd_o
          "I see this heat transfer area if I am outside the tubes";

          //Area, global heat transfer coefficient and corrected DMTL
          A              = N_tubes*pi*l*bundle.Dhyd_o;
          U              = sum(kA_tot)/A;
          DTML_tilde     = qdot/(U*A);

          //Weight calculation
          bundle.W_dry = 0.25*pi*(bundle.Dhyd_o^2 - bundle.Dhyd^2)*l*N_tubes
                           *bundle.rho;
          shell.W_dry  = 0.25*pi*(shell.Dhyd_o^2 - shell.Dhyd^2)*l*shell.rho;
          W_dry          = 1.25*bundle.W_dry + 1.3*N_heads*shell.W_dry;
          W_fluids       = sum(bundle.W_wet) + sum(shell.W_wet);
          W_wet          = W_dry + W_fluids;

          //Cost calculation
          PEC            = cost(A);

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
  extends Modelica.Icons.ExamplesPackage;
    model Coulson_Kern
      "Verification with the results given by Coulson et al. using the Kern method"

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

      Components.HEX.shell_and_tube             shell_and_tube(
        Dhyd=16e-3,
        l=4.83,
        pitch_f=1.25,
        layout=1,
        N_passes=2,
        N_baffles=27,
        h_s_in=h_s_in,
        h_s_out=h_s_out,
        h_t_in=h_t_in,
        h_t_out=h_t_out,
        m_s=1e2/3.6,
        m_t=68.9,
        redeclare package Medium_s = VIP.Media.Methanol_CoolProp,
        redeclare package Medium_t = VIP.Media.Water_CoolProp,
        redeclare VIP.Materials.S_AISI_1040 Material_t,
        redeclare VIP.Materials.S_AISI_1040 Material_s,
        redeclare function bundle_clearance =
            VIP.Miscellanea.Shell_clearance.SRFH,
        pin_s=1,
        pin_t=-1,
        thick_t=2e-3,
        thick_s=10e-3,
        N_heads=1,
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
      "Verification with the results given by Coulson using the Bell Delaware method."

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

       Components.HEX.shell_and_tube  st(
         redeclare package Medium_s = Medium_s,
         redeclare package Medium_t = Medium_t,
         Dhyd =        16e-3,
         l =           4.83,
         pitch_f =     1.25,
         layout =      1,
         N_passes =    2,
         N_baffles =   12,
         N_baffles_d = 25,
         pin_t =       -1,
         pin_s =       1,
         h_s_in =      h_s_in,
         h_s_out =     h_s_out,
         h_t_in =      h_t_in,
         h_t_out =     h_t_out,
         m_s =         1e2/3.6,
         m_t =         68.9,
         redeclare VIP.Heat_transfer.Shell.Bell_Delaware hT_shell(pitch_f = st.pitch_f, Dhyd_o = st.bundle.Dhyd_o,
         N_tubes = st.N_tubes, d_s = st.d_s, d_b = st.bundle.d_b, Dhyd = st.d_s,
         layout = st.layout, l_b = st.l_b, N_ss = 5,
         eta_wall=st.shell.state[1].eta*ones(st.Ncell),
         N_passes = st.N_passes, ff =  1.0, b_cut = 0.25),
         redeclare VIP.Pressure_drops.Shell.Bell_Delaware dp_shell(pitch_f = st.pitch_f, Dhyd_o = st.bundle.Dhyd_o,
         N_tubes = st.N_tubes, d_s = st.d_s, d_b = st.bundle.d_b, Dhyd = st.d_s,
         layout = st.layout, l_b = st.l_b, N_ss = 5,
         eta_wall=st.shell.state[1].eta*ones(st.Ncell),
         N_passes = st.N_passes, N_baffles_d = st.N_baffles_d,
         N_baffles = st.N_baffles, ff =  1.0, b_cut = 0.25),
        redeclare VIP.Materials.S_AISI_1040 Material_t,
        redeclare VIP.Materials.S_AISI_1040 Material_s,
        redeclare function bundle_clearance =
            VIP.Miscellanea.Shell_clearance.SRFH,
        thick_t=2e-3,
        thick_s=10e-3,
        N_heads=1,
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

       Components.HEX.shell_and_tube  st(
         redeclare package Medium_s = Medium_s,
         redeclare package Medium_t = Medium_t,
         Dhyd =        16e-3,
         thick_t =     2e-3,
         thick_s =     10e-3,
         l =           4.65,
         pitch_f =     1.25,
         layout =      1,
         N_passes =    2,
         N_baffles =   24,
         N_baffles_d = 25,
         pin_t =       -1,
         pin_s =       1,
         h_s_in =      h_s_in,
         h_s_out =     h_s_out,
         h_t_in =      h_t_in,
         h_t_out =     h_t_out,
         m_s =         1e2/3.6,
         m_t =         68.9,
         redeclare VIP.Heat_transfer.Shell.Bell_Delaware hT_shell(pitch_f = st.pitch_f, Dhyd_o = st.bundle.Dhyd_o,
         N_tubes = st.N_tubes, d_s = st.d_s, d_b = st.bundle.d_b, Dhyd = st.d_s,
         layout = st.layout, l_b = st.l_b, N_ss = 5,
         eta_wall=st.shell.state[1].eta*ones(st.Ncell),
         N_passes = st.N_passes, ff =  1.0, b_cut = 0.1588),
         redeclare VIP.Pressure_drops.Shell.Bell_Delaware dp_shell(pitch_f = st.pitch_f, Dhyd_o = st.bundle.Dhyd_o,
         N_tubes = st.N_tubes, d_s = st.d_s, d_b = st.bundle.d_b, Dhyd = st.d_s,
         layout = st.layout, l_b = st.l_b, N_ss = 5,
         eta_wall=st.shell.state[1].eta*ones(st.Ncell),
         N_passes = st.N_passes, N_baffles_d = st.N_baffles_d,
         N_baffles = st.N_baffles, ff =  1.0, b_cut = 0.1588),
        redeclare VIP.Materials.S_AISI_1040 Material_t,
        redeclare VIP.Materials.S_AISI_1040 Material_s,
        redeclare function bundle_clearance =
            VIP.Miscellanea.Shell_clearance.Fixed_Utube,
        N_heads=1,
        p_s_in=301435,
        p_s_out=301435,
        p_t_in=100000,
        p_t_out=100000,
        redeclare function cost = VIP.Miscellanea.Cost.CS_Hall)
        annotation (Placement(transformation(extent={{-80,-68},{68,58}})));

     annotation (Placement(transformation(extent={{-108,-74},{88,66}})), Icon(
            coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{
                100,100}}), graphics));
    end Aspen_Bell_Delaware;

    model Aspen_Bell_Delaware_upsidedown
      "Verification with the results given by Aspen using the Bell Delaware method. We switch the fluid allocation here."

      replaceable package Medium_t = VIP.Media.Methanol_CoolProp "Medium model";
      replaceable package Medium_s = VIP.Media.Water_CoolProp "Medium model";
      parameter Modelica.SIunits.SpecificEnthalpy h_t_in= Medium_t.specificEnthalpy_pT(3.2e5,368.15)
        "Inlet specific enthalpy shell side";
      parameter Modelica.SIunits.SpecificEnthalpy h_t_out= Medium_t.specificEnthalpy_pT(3.2e5,313.15)
        "Outlet specific enthalpy shell side";
      parameter Modelica.SIunits.SpecificEnthalpy h_s_in = Medium_s.specificEnthalpy_pT(1e5,298.15)
        "Inlet specific enthalpy tube side";
      parameter Modelica.SIunits.SpecificEnthalpy h_s_out = Medium_s.specificEnthalpy_pT(1e5,313.15)
        "Outlet specific enthalpy tube side";

       Components.HEX.shell_and_tube  st(
         redeclare package Medium_s = Medium_s,
         redeclare package Medium_t = Medium_t,
         pin_t = 1,
         pin_s = -1,
         Dhyd = 16e-3,
         l = 4.65,
         pitch_f=1.25,
         layout = 1,
         N_passes = 4,
         N_baffles = 24,
         h_s_in = h_s_in,
         h_s_out = h_s_out,
         h_t_in = h_t_in,
         h_t_out = h_t_out,
         m_t = 1e2/3.6,
         m_s = 68.9,
         DTML = 50,
         redeclare VIP.Heat_transfer.Shell.Bell_Delaware hT_shell(pitch_f = st.pitch_f, Dhyd_o = st.bundle.Dhyd_o,
         N_tubes = st.N_tubes, d_s = st.d_s, d_b = st.bundle.d_b, Dhyd = st.d_s,
         layout = st.layout, l_b = st.l_b, N_ss = 5,
         eta_wall=st.shell.state[1].eta*ones(st.Ncell),
         N_passes = st.N_passes, ff =  1.0, b_cut = 0.1588),
         redeclare VIP.Pressure_drops.Shell.Bell_Delaware dp_shell(pitch_f = st.pitch_f, Dhyd_o = st.bundle.Dhyd_o,
         N_tubes = st.N_tubes, d_s = st.d_s, d_b = st.bundle.d_b, Dhyd = st.d_s,
         layout = st.layout, l_b = st.l_b, N_ss = 5,
         eta_wall=st.shell.state[1].eta*ones(st.Ncell),
         N_passes = st.N_passes, N_baffles_d = st.N_baffles_d,
         N_baffles = st.N_baffles, ff =  1.0, b_cut = 0.1588),
        redeclare VIP.Materials.S_AISI_1040 Material_t,
        redeclare VIP.Materials.S_AISI_1040 Material_s,
        redeclare function bundle_clearance =
            VIP.Miscellanea.Shell_clearance.Fixed_Utube,
        thick_t=2e-3,
        thick_s=10e-3,
        N_heads=1,
        p_s_in=100000,
        p_s_out=100000,
        p_t_in=301435,
        p_t_out=301435)
        annotation (Placement(transformation(extent={{-80,-68},{68,58}})));

     annotation (Placement(transformation(extent={{-108,-74},{88,66}})));
    end Aspen_Bell_Delaware_upsidedown;

    model DYNDES_Bell_Delaware
      "Verification with the results given by the DYNDES simulation tool using the Bell Delaware method."

      replaceable package Medium_s = VIP.Media.Cyclopentane_CoolProp
        "Medium model";
      replaceable package Medium_t = VIP.Media.Cyclopentane_CoolProp
        "Medium model";
      parameter Modelica.SIunits.SpecificEnthalpy h_s_in= Medium_s.specificEnthalpy_pT(1.1e5,403.15)
        "Inlet specific enthalpy shell side";
      parameter Modelica.SIunits.SpecificEnthalpy h_s_out= Medium_s.specificEnthalpy_pT(1.1e5,357.15)
        "Outlet specific enthalpy shell side";
      parameter Modelica.SIunits.SpecificEnthalpy h_t_in = Medium_t.specificEnthalpy_pT(3e6,323.15)
        "Inlet specific enthalpy tube side";
      parameter Modelica.SIunits.SpecificEnthalpy h_t_out = Medium_t.specificEnthalpy_pT(3e6,323.15)
        "Outlet specific enthalpy tube side";

      Components.HEX.shell_and_tube  st(
        redeclare package Medium_s = Medium_s,
        redeclare package Medium_t = Medium_t,
        Dhyd=16e-3,
        l=7.32,
        pitch_f=1.25,
        layout=1,
        N_passes=2,
        N_baffles=3,
        N_baffles_d=4,
        pin_t =       -1,
        pin_s =       1,
        h_s_in=h_s_in,
        h_s_out=h_s_out,
        h_t_in=h_t_in,
        h_t_out=h_t_out,
        m_s = 40,
        m_t = 40,
        redeclare VIP.Heat_transfer.Shell.Bell_Delaware   hT_shell(pitch_f = st.pitch_f, Dhyd_o = st.bundle.Dhyd_o,
        N_tubes = st.N_tubes, d_s = st.d_s, d_b = st.bundle.d_b, Dhyd = st.d_s,
        layout = st.layout, l_b = st.l_b, N_ss = 5,
        eta_wall=st.shell.state[1].eta*ones(st.Ncell),
        N_passes = st.N_passes, ff =  1.0, b_cut = 0.25),
        redeclare VIP.Pressure_drops.Shell.Bell_Delaware   dp_shell(pitch_f = st.pitch_f, Dhyd_o = st.bundle.Dhyd_o,
        N_tubes = st.N_tubes, d_s = st.d_s, d_b = st.bundle.d_b, Dhyd = st.d_s,
        layout = st.layout, l_b = st.l_b, N_ss = 5,
        eta_wall=st.shell.state[1].eta*ones(st.Ncell),
        N_baffles = st.N_baffles, N_baffles_d = st.N_baffles_d,
        N_passes = st.N_passes, ff =  1.0, b_cut = 0.25),
        redeclare function bundle_clearance =
          VIP.Miscellanea.Shell_clearance.SRFH,
        redeclare VIP.Materials.S_AISI_1040 Material_t,
        redeclare VIP.Materials.S_AISI_1040 Material_s,
        thick_t=2e-3,
        thick_s=10e-3,
        N_heads=1,
        p_s_in=110000,
        p_s_out=110000,
        p_t_in=3000000,
        p_t_out=3000000)
        annotation (Placement(transformation(extent={{-80,-68},{68,58}})));

     annotation (Placement(transformation(extent={{-108,-74},{88,66}})));
    end DYNDES_Bell_Delaware;

    model Optimization "Shell and tube model for external optmization"

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

      Components.HEX.shell_and_tube             shell_and_tube(
        Dhyd=16e-3,
        l=4.83,
        pitch_f=1.25,
        layout=1,
        N_passes=2,
        N_baffles=27,
        h_s_in=h_s_in,
        h_s_out=h_s_out,
        h_t_in=h_t_in,
        h_t_out=h_t_out,
        m_s=1e2/3.6,
        m_t=68.9,
        redeclare package Medium_s = VIP.Media.Methanol_CoolProp,
        redeclare package Medium_t = VIP.Media.Water_CoolProp,
        redeclare VIP.Materials.S_AISI_1040 Material_t,
        redeclare VIP.Materials.S_AISI_1040 Material_s,
        redeclare function bundle_clearance =
            VIP.Miscellanea.Shell_clearance.SRFH,
        pin_s=1,
        pin_t=-1,
        thick_t=2e-3,
        thick_s=10e-3,
        N_heads=1,
        p_s_in=301435,
        p_s_out=301435,
        p_t_in=100000,
        p_t_out=100000)
        annotation (Placement(transformation(extent={{-80,-68},{68,58}})));

     annotation (Placement(transformation(extent={{-108,-74},{88,66}})),
        experiment,
        __Dymola_experimentSetupOutput);
    end Optimization;

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
  end Tests;

  annotation (uses(Modelica(version="3.2.1")),
                                             Icon(coordinateSystem(
          preserveAspectRatio=false, extent={{-100,-100},{100,100}}), graphics));
end VIP;
