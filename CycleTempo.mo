within ;
package CycleTempo "Cycle Tempo 2_0"
  import Modelica.Constants.*;
  package Components "Component library"
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
        annotation(Diagram(coordinateSystem(preserveAspectRatio = true, extent = {{-100, -100}, {100, 100}}), graphics={  Ellipse(extent=  {{-100, 100}, {100, -100}}, lineColor=
                    {0,127,255},                                                                                                    fillColor=
                    {0,0,0},
                  fillPattern=FillPattern.Solid),                                                                                                    Ellipse(extent=  {{-100, 100}, {100, -100}}, lineColor=
                    {0,0,0}),                                                                                                    Text(extent={{
                    -94,192},{106,98}},                                                                                                    textString=  "%name", lineColor=
                    {0,0,0})}),                                                                                                    Icon(coordinateSystem(preserveAspectRatio = true, extent = {{-100, -100}, {100, 100}}), graphics={  Ellipse(extent=  {{-100, 100}, {100, -100}},
                pattern=LinePattern.None,
                lineColor={0,0,0}),                                                                                                    Ellipse(extent=  {{-100, 100}, {100, -100}}, lineColor=
                    {0,0,0},
                fillColor={0,0,0},
                fillPattern=FillPattern.Solid)}),                                                                                                    Documentation(info = "<html>Modelica.Media.Examples.Tests.Components.FluidPort_a
                                                                                                </html>"));
      end Node_out;

      connector power "A node with a thermodynamic state"
        Modelica.SIunits.Power W "Power";
      end power;

      connector terminal "Electric terminal"
        extends power;
        annotation(Diagram(coordinateSystem(preserveAspectRatio = true, extent = {{-100, -100}, {100, 100}}), graphics={         Text(extent={{
                    -98,202},{102,108}},                                                                                                    textString=  "%name", lineColor=
                    {255,128,0}), Rectangle(
                extent={{-100,100},{100,-100}},
                lineColor={0,0,0},
                fillColor={255,255,0},
                fillPattern=FillPattern.Solid)}),                                                                                                    Icon(coordinateSystem(preserveAspectRatio=false,  extent={{-100,
                  -100},{100,100}}),                                                                                                    graphics={
                Rectangle(
                extent={{-100,100},{100,-100}},
                lineColor={0,0,0},
                fillColor={255,255,0},
                fillPattern=FillPattern.Solid)}),                                                                                                    Documentation(info = "<html>Modelica.Media.Examples.Tests.Components.FluidPort_a
                                                                                                </html>"));
      end terminal;
      annotation (Icon(graphics={
            Ellipse(lineColor={255,0,0},
              pattern=LinePattern.None,
              extent={{-100,-100},{100,100}},
              fillColor={255,0,0},
              fillPattern=FillPattern.Sphere)}));
    end Nodes;

    package Turbomachinery "Pump and turbine models"
      model Pump
        replaceable package Medium = Test.Media.OneRandomOrganicFluid constrainedby
          Modelica.Media.Interfaces.PartialMedium "Medium model hot cells" annotation(choicesAllMatching = true);
        parameter Real eta_is "Isentropic efficiency";
        parameter Real eta_m "Mechanical efficiency";
        parameter Boolean activated = false "Off design mode on if true" annotation (Dialog(tab="Off-design"));
        parameter Boolean constant_n = false
          "Constant rotational speed if true"                                    annotation (Dialog(tab="Off-design"));
        parameter Modelica.SIunits.Height head_d = 1 "Head at design" annotation (Dialog(tab="Off-design"));
        parameter Modelica.SIunits.Conversions.NonSIunits.AngularVelocity_rpm n_d = 3e3
          "rotational speed at design" annotation (Dialog(tab="Off-design"));
        parameter Modelica.SIunits.VolumeFlowRate V_flow_d = 1
          "Volume flow rate at design" annotation (Dialog(tab="Off-design"));
        Modelica.SIunits.Power W "Power consumption";
        Modelica.SIunits.SpecificEnthalpy h_is
          "Outlet isentropic specific enthalpy";
        Modelica.SIunits.SpecificEntropy s "Inlet specific entropy";

        replaceable
          CycleTempo.Components.Turbomachinery.Part_load.Pump.Exponential
         off_design constrainedby
          CycleTempo.Components.Turbomachinery.Part_load.Pump.Base_classes.base(
          redeclare package Medium = Medium,
          final offdesign=activated,
          final constant_n=constant_n,
          final node_in=node_in,
          final node_out=node_out,
          final head_d=head_d,
          final V_flow_d=V_flow_d,
          final n_d=n_d) annotation (Dialog(tab="Off-design"));

        Nodes.Node_in node_in(redeclare package Medium = Medium) "Inlet node"
          annotation (Placement(transformation(extent={{-108,-10},{-88,10}}),
              iconTransformation(extent={{-108,-10},{-88,10}})));
        Nodes.Node_out node_out(redeclare package Medium = Medium)
          "Outlet node"
          annotation (Placement(transformation(extent={{90,-10},{110,10}}),
              iconTransformation(extent={{90,-10},{110,10}})));
        Nodes.terminal terminal annotation (Placement(transformation(extent={{-9,90},{
                  11,110}}),  iconTransformation(extent={{-11,-110},{10,-90}})));
      equation

        //Mass balance
        node_in.m_flow = node_out.m_flow;

        //Energy balance
        W              = node_in.m_flow*(node_out.h - node_in.h)/eta_m;
        W              = terminal.W;

        //Component equations
        s              = Medium.specificEntropy(Medium.setState_phX(node_in.p, node_in.h));
        h_is           = Medium.specificEnthalpy(Medium.setState_psX(node_out.p, s));
        node_out.h     = (h_is - node_in.h)/eta_is + node_in.h;

        annotation(Documentation(info = "<HTML>
                                                                                              <p>This is a cycle model for pumps.
                                                                                              </html>"),  Icon(coordinateSystem(extent={{-100,
                  -100},{100,100}},                                                                                                    preserveAspectRatio=false,  initialScale = 0.1, grid = {2, 2}), graphics={                       Ellipse(
                                                                                                  extent={{
                    -60,60},{60,-60}},                                                                                                    endAngle=  360,
                lineColor={0,0,0},
                fillColor={255,255,255},
                fillPattern=FillPattern.Solid),                                                                                                    Polygon(origin={
                    -20.716,0.7156},                                                                                                    rotation=  180, fillColor=
                    {255,255,255},
                  fillPattern=FillPattern.HorizontalCylinder,                                                                                                    points={{
                    19.284,42.7156},{19.284,-41.2844},{-80.7156,0.715564},{
                    19.284,42.7156}},
                lineColor={0,0,0}),                                                                                               Text(lineColor=
                    {0,0,0},                                                                                                    extent={{
                    -100,15},{100,-15}},
                origin={0,114},
                rotation=0,
                textString="%name"),
              Line(
                points={{60,0},{100,0}},
                color={0,0,0},
                smooth=Smooth.None),
              Line(
                points={{-100,0},{-60,0}},
                color={0,0,0},
                smooth=Smooth.None),
              Line(
                points={{0,-100},{0,-80},{0,-60}},
                color={0,0,0},
                smooth=Smooth.None)}),
          Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,
                  100}}), graphics));
      end Pump;

      model Turbine "Turbine model"

        replaceable package Medium = Test.Media.OneRandomOrganicFluid constrainedby
          Modelica.Media.Interfaces.PartialMedium "Medium model hot cells" annotation(choicesAllMatching = true);
        parameter Real eta_is "Isentropic efficiency";
        parameter Real eta_m "Mechanical efficiency";
        parameter Boolean activated = false "Off design mode on if true" annotation (Dialog(tab="Off-design"));
        Modelica.SIunits.Power W "Power production";
        Modelica.SIunits.SpecificEnthalpy h_is
          "Outlet isentropic specific enthalpy";
        Modelica.SIunits.SpecificEntropy s "Inlet specific entropy";
        replaceable
          CycleTempo.Components.Turbomachinery.Part_load.Turbine.Stodola
         off_design constrainedby
          CycleTempo.Components.Turbomachinery.Part_load.Turbine.Base_classes.base(
          redeclare package Medium = Medium,
          offdesign=activated,
          node_in=node_in,
          node_out=node_out) annotation (choicesAllMatching = true, Dialog(tab="Off-design"));

        Nodes.Node_in node_in(redeclare package Medium = Medium) "Inlet node"
          annotation (Placement(transformation(extent={{-50,90},{-30,110}}),
              iconTransformation(extent={{-50,90},{-30,110}})));
        Nodes.Node_out node_out(redeclare package Medium = Medium)
          "Outlet node"
          annotation (Placement(transformation(extent={{4,-130},{24,-110}}),
              iconTransformation(extent={{4,-130},{24,-110}})));
        Nodes.terminal terminal annotation (Placement(transformation(extent={{91,-10},
                  {111,10}}), iconTransformation(extent={{89,-10},{110,10}})));
      equation

        //Mass balance
        node_in.m_flow = node_out.m_flow;

        //Energy balance
        W              = node_in.m_flow*(node_in.h - node_out.h)*eta_m;
        W              = terminal.W;

        //Component equations
        s              = Medium.specificEntropy(Medium.setState_phX(node_in.p, node_in.h));
        h_is           = Medium.specificEnthalpy(Medium.setState_psX(node_out.p, s));
        node_out.h     = node_in.h - (node_in.h - h_is)*eta_is;
        connect(node_out, node_out) annotation (Line(
            points={{14,-120},{42,-120},{42,-140},{44,-140},{44,-120},{14,
                -120}},
            color={0,127,255},
            smooth=Smooth.None));
        annotation(Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,
                  -100},{100,100}}),
                        graphics={                                                                                                    Rectangle(extent={{
                    14,10},{100,-10}},                                                                                                  lineColor=
                    {0,0,0},
                fillColor={255,255,255},
                fillPattern=FillPattern.Solid),                                                                                                    Polygon(points={{
                    -40,28},{-40,-26},{14,-80},{14,80},{-40,28}},                                                                                                    lineColor=
                    {0,0,0},
                fillColor={255,255,255},
                fillPattern=FillPattern.Solid),                                                   Text(lineColor=
                    {0,0,0},                                                                                                    extent={{
                    -100,15},{100,-15}},                                                                                                    textString=  "%name",
                origin={48,114},
                rotation=0),
              Line(
                points={{14,-80},{14,-120}},
                color={0,0,0},
                smooth=Smooth.None),
              Line(
                points={{-40,28},{-40,100}},
                color={0,0,0},
                smooth=Smooth.None)}),                                                                                            Diagram(
              coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},
                  {100,100}}),                                                                                                    graphics), Documentation(info = "<html>
                                                                                                   <p>This is a cycle model for turbine
                                                                                                   </html>"));
      end Turbine;

      model Compressor "Compressor model"
        replaceable package Medium = Test.Media.OneRandomOrganicFluid constrainedby
          Modelica.Media.Interfaces.PartialMedium "Medium model hot cells" annotation(choicesAllMatching = true);
        parameter Real eta_is "Isentropic efficiency";
        parameter Real eta_m "Mechanical efficiency";
        Modelica.SIunits.Power W "Power consumption";
        Modelica.SIunits.SpecificEnthalpy h_is
          "Outlet isentropic specific enthalpy";
        Modelica.SIunits.SpecificEntropy s "Inlet specific entropy";
        Nodes.Node_in node_in(redeclare package Medium = Medium) "Inlet node"
          annotation (Placement(transformation(extent={{-50,90},{-30,110}}),
              iconTransformation(extent={{14,-120},{34,-100}})));
        Nodes.Node_out node_out(redeclare package Medium = Medium)
          "Outlet node"
          annotation (Placement(transformation(extent={{4,-130},{24,-110}}),
              iconTransformation(extent={{-40,100},{-20,120}})));
        Nodes.terminal terminal annotation (Placement(transformation(extent={{91,-10},
                  {111,10}}), iconTransformation(extent={{99,0},{120,20}})));
      equation

        //Mass balance
        node_in.m_flow = node_out.m_flow;

        //Energy balance
        W              = node_in.m_flow*(node_out.h - node_in.h)/eta_m;
        W              = terminal.W;

        //Component equations
        s              = Medium.specificEntropy(Medium.setState_phX(node_in.p, node_in.h));
        h_is           = Medium.specificEnthalpy(Medium.setState_psX(node_out.p, s));
        node_out.h     = (h_is - node_in.h)/eta_is + node_in.h;

        annotation(Documentation(info = "<HTML>
                                                                                              <p>This is a cycle model for pumps.
                                                                                              </html>"),  Icon(coordinateSystem(extent={{-100,
                  -100},{100,100}},                                                                                                    preserveAspectRatio=false,  initialScale = 0.1, grid = {2, 2}), graphics={
                                                                                                    Rectangle(extent={{
                    24,20},{110,0}},                                                                                                    lineColor=
                    {0,0,0},
                fillColor={255,255,255},
                fillPattern=FillPattern.Solid),                                                                                                    Polygon(points={{
                    -30,38},{-30,-16},{24,-70},{24,90},{-30,38}},                                                                                                    lineColor=
                    {0,0,0},
                fillColor={255,255,255},
                fillPattern=FillPattern.Solid),                                                   Text(lineColor=
                    {0,0,0},                                                                                                    extent={{
                    -100,15},{100,-15}},                                                                                                    textString=  "%name",
                origin={58,124},
                rotation=0),
              Line(
                points={{24,-70},{24,-110}},
                color={0,0,0},
                smooth=Smooth.None),
              Line(
                points={{-30,38},{-30,110}},
                color={0,0,0},
                smooth=Smooth.None)}),
          Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,
                  100}}), graphics));
      end Compressor;

      package Part_load "Curves to model the components at off-design"
        package Turbine "Part load models for turbines"
          class de_Laval "de Laval nozzle (chocked)"
            extends Base_classes.base;
            parameter Modelica.SIunits.Area A_d = 1
              "Total nozzle area at design";
            parameter Modelica.SIunits.AbsolutePressure p_is_t_start = 1e5
              "Total nozzle area at design";
            Modelica.SIunits.Velocity c "Speed of sound at the throat";
            Modelica.SIunits.Area A "Total nozzle area";
            Medium.ThermodynamicState state_in "Inlet thermodynamic state";
            Medium.ThermodynamicState state_t "Throat thermodynamic state";
            Modelica.SIunits.SpecificEnthalpy h_is_t
              "Isentropic specific enthalpy at the throat";
            Modelica.SIunits.AbsolutePressure p_is_t(start=p_is_t_start)
              "Isentropic specific enthalpy at the throat";
          equation
            state_in = Medium.setState_ph(node_in.p, node_in.h);
            h_is_t   = Medium.specificEnthalpy_ps(p_is_t, state_in.s);
            state_t  = Medium.setState_ps(p_is_t, state_in.s);
            c        = Medium.velocityOfSound(state_t);
            h_is_t   = state_in.h - 0.5*c^2;
            node_in.m_flow = A*state_t.d*c;

            if offdesign then
              A = A_d;
            end if;

          end de_Laval;

          class Stodola "Stodola part load model"
            extends Base_classes.base;
            parameter Real CT_d = 1 "Stodola constant user defined";
            Real CT "Stodola constant calculated";
            Medium.ThermodynamicState state_in "Inlet thermodynamic state";
          equation
            state_in = Medium.setState_phX(node_in.p, node_in.h);
            CT       = node_in.m_flow*Design.Miscellanea.sqrtReg(state_in.T)/
            Design.Miscellanea.sqrtReg(node_in.p^2 - node_out.p^2);

            if offdesign then
              CT = CT_d;
            end if;

          end Stodola;

          package Base_classes "Base classes for the turbine"
            class base "Basic part load model"
              replaceable package Medium = Test.Media.OneRandomOrganicFluid constrainedby
                Modelica.Media.Interfaces.PartialMedium "Medium model" annotation(choicesAllMatching = true);
              parameter Boolean offdesign = false "Off design mode on if true";
              input Nodes.Node_in node_in(redeclare package Medium = Medium);
              input Nodes.Node_out node_out(redeclare package Medium = Medium);
            end base;
          end Base_classes;
        end Turbine;

        package Pump "Part load models for pumps"
          class Exponential
            "Head versus volume flow rate part load model. The curve is exponential."
            extends Base_classes.base;
            parameter Real c[2] = {2.462215552,  -0.53791904}
              "Default coefficients for head vs. volumetric flow curve";
            Modelica.SIunits.Conversions.NonSIunits.AngularVelocity_rpm n(start=n_d)
              "rotational speed at off-design";
          equation
            if offdesign then
              head = head_d*(c[1] + c[2]*exp(V_flow/V_flow_d))*(n/n_d)^2;
              if constant_n then
                n  = n_d;
              end if;
            else
              n    = n_d;
            end if;

          end Exponential;

          package Base_classes "Base classes for the pump"
            class base "Basic part load model"
              replaceable package Medium = Test.Media.OneRandomOrganicFluid constrainedby
                Modelica.Media.Interfaces.PartialMedium "Medium model" annotation(choicesAllMatching = true);
              parameter Boolean offdesign = false "Off design mode on if true";
              parameter Boolean constant_n = false
                "Constant rotational speed if true";
              parameter Modelica.SIunits.Height head_d "Head at design";
              parameter
                Modelica.SIunits.Conversions.NonSIunits.AngularVelocity_rpm         n_d
                "rotational speed at design";
              parameter Modelica.SIunits.VolumeFlowRate V_flow_d
                "Volume flow rate at design";
              input Nodes.Node_in node_in(redeclare package Medium = Medium);
              input Nodes.Node_out node_out(redeclare package Medium = Medium);
              Modelica.SIunits.Height head(start=head_d) "Head at off-design";
              Modelica.SIunits.VolumeFlowRate V_flow(start=V_flow_d)
                "Volume flow rate at off-design";
              Medium.ThermodynamicState state_in
                "Thermodynamic state at the in";
            equation
              state_in = Medium.setState_phX(node_in.p, node_in.h);
              V_flow   = node_in.m_flow/state_in.d;
              head     = (node_out.h - node_in.h)/g_n;
            end base;
          end Base_classes;
        end Pump;
      end Part_load;
    end Turbomachinery;

    package HEX "Heat exchangers"
      model Evaporator
        "Evaporator with three zones preheater, evaporator, superheater"
        extends CycleTempo.Components.HEX.Heat_exchanger;
        parameter Modelica.SIunits.TemperatureDifference dT_int = 0
          "Temperature difference at the inlet of the evaporator" annotation (Dialog(tab="Addco"));
        parameter Boolean use_dT_int = false
          "true if temperature difference at the inlet of the evaporator is given"
                                                                                   annotation (Dialog(tab="Addco"));
        parameter Modelica.SIunits.SpecificEnthalpy hh_in_start
          "Inlet enthalpy evaporator hot side start value" annotation (Dialog(tab="Start"));
        Medium_c.SaturationProperties sat "Saturated state cold fluid";
        Medium_h.ThermodynamicState hot_sat
          "Hot fluid state corresponding to sat";
        Modelica.SIunits.SpecificEnthalpy hh_in_pre(start = hh_in_start)
          "Inlet enthalpy pre-heater hot side";
        Modelica.SIunits.SpecificEnthalpy hh_in_eva
          "Inlet enthalpy evaporator hot side ";
        Modelica.SIunits.HeatFlowRate qdot_pre "Heat rate pre-heater";
        Modelica.SIunits.HeatFlowRate qdot_eva "Heat rate evaporator";
        Modelica.SIunits.HeatFlowRate qdot_sup "Heat rate superheater";
        Modelica.SIunits.TemperatureDifference dT_eva
          "Temperature difference at the inlet of the evaporator";
      equation
        //Boundary equations
        if use_dT_int then
          dT_int         = dT_eva;
        end if;

        //Component equations
        sat              = Medium_c.setSat_p(node_c_in.p);
        hot_sat          = Medium_h.setState_phX(node_h_in.p, hh_in_pre);
        dT_eva           = hot_sat.T - sat.Tsat;
        qdot_pre         = node_c_in.m_flow*(sat.hl - node_c_in.h);
        qdot_eva         = node_c_in.m_flow*(sat.hv - sat.hl);
        qdot_sup         = node_c_in.m_flow*(node_c_out.h - sat.hv);
        qdot_pre         = node_h_in.m_flow*(hh_in_pre - node_h_out.h);
        qdot_eva         = node_h_in.m_flow*(hh_in_eva - hh_in_pre);

        //Check second principle of Thermodynamics
        assert(dT_eva > 0, "Second principle of Thermodynamics not respected");
        assert(dT_eva > 5, "Internal temperature difference lower than 5 K",
        AssertionLevel.warning);

        annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
                  -100},{100,100}}), graphics), Icon(coordinateSystem(
                preserveAspectRatio=false, extent={{-100,-100},{100,100}}), graphics));
      end Evaporator;

      model Condenser "Condenser model"
        extends CycleTempo.Components.HEX.Heat_exchanger;
        parameter Modelica.SIunits.TemperatureDifference dT_int = 0
          "Temperature difference at the inlet of the condensation region" annotation (Dialog(tab="Addco"));
        parameter Boolean use_dT_int = false
          "true if temperature difference at the inlet of the condensation region is given"
        annotation (Dialog(tab="Addco"));
        parameter Modelica.SIunits.SpecificEnthalpy hc_in_start
          "Inlet enthalpy condensation region cold side start value" annotation (Dialog(tab="Start"));
        Modelica.SIunits.SpecificEnthalpy hc_in_con( start = hc_in_start)
          "Inlet enthalpy condensation region cold side";
        Modelica.SIunits.HeatFlowRate qdot_con "Heat rate condensation region";
        Modelica.SIunits.TemperatureDifference dT_con
          "Temperature difference at the inlet of the condensation region";
        Medium_h.SaturationProperties sat "Saturated state hot fluid";
        Medium_c.ThermodynamicState cold_sat
          "Cold fluid state corresponing to sat";
      equation

        //Boundary equations
        if use_dT_int then
          dT_int         = dT_con;
        end if;

        //Component equations
        sat          = Medium_h.setSat_p(node_h_out.p);
        node_h_out.h = sat.hl;
        qdot_con     = node_c_in.m_flow*(hc_in_con - node_c_in.h);
        qdot_con     = node_h_in.m_flow*(sat.hv - sat.hl);
        cold_sat     = Medium_c.setState_phX(node_c_in.p, hc_in_con);
        dT_con       = sat.Tsat - cold_sat.T;

        //Check second principle of Thermodynamics
        assert(dT_con > 0, "Second principle of Thermodynamics not respected");
        assert(dT_con > 5, "Internal temperature difference lower than 5 K",
        AssertionLevel.warning);

        annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,
                  -100},{100,100}}),       graphics));
      end Condenser;

      model Heat_exchanger "Heat exchanger for single phase heat transfer"
        //Medium model
        replaceable package Medium_h = Test.Media.OneRandomOrganicFluid constrainedby
          Modelica.Media.Interfaces.PartialMedium "Medium model hot cells" annotation(choicesAllMatching = true);
        replaceable package Medium_c = Test.Media.OneRandomOrganicFluid constrainedby
          Modelica.Media.Interfaces.PartialMedium "Medium model hot cells" annotation(choicesAllMatching = true);
        //Input parameters
        parameter Modelica.SIunits.AbsolutePressure dp_h
          "Pressure drop hot side";
        parameter Modelica.SIunits.AbsolutePressure dp_c
          "Pressure drop cold side";
        parameter Modelica.SIunits.TemperatureDifference dT_in = 0
          "Inlet temperature difference" annotation (Dialog(tab="Addco"));
        parameter Modelica.SIunits.TemperatureDifference dT_out = 0
          "Outlet temperature difference" annotation (Dialog(tab="Addco"));
        parameter Modelica.SIunits.HeatFlowRate power = 0 "Heat rate" annotation (Dialog(tab="Addco"));
        parameter Boolean use_dT_in = false
          "true if inlet temperature difference is given" annotation (Dialog(tab="Addco"));
        parameter Boolean use_dT_out = false
          "true if outlet temperature difference is given" annotation (Dialog(tab="Addco"));
        parameter Boolean use_power = false "true if heat rate is given"  annotation (Dialog(tab="Addco"));
        Modelica.SIunits.HeatFlowRate qdot "Heat rate";
        Nodes.Node_in node_c_in(redeclare package Medium = Medium_c)
          "Inlet node cold side" annotation (Placement(transformation(extent={{90,-70},{110,-50}}),
              iconTransformation(
              extent={{-10,-10},{10,10}},
              rotation=-90,
              origin={-140,0})));
        Nodes.Node_out node_c_out(redeclare package Medium = Medium_c)
          "Outlet node cold side"  annotation (Placement(transformation(extent={{90,50},{110,70}}),
              iconTransformation(
              extent={{-10,-10},{10,10}},
              rotation=-90,
              origin={141,0})));
        Nodes.Node_in node_h_in(redeclare package Medium = Medium_h)
          "Inlet node hot side"   annotation (Placement(transformation(extent={{-10,-110},{10,-90}}),
              iconTransformation(extent={{-10,-10},{10,10}},
              rotation=-90,
              origin={0,142})));
        Nodes.Node_out node_h_out(redeclare package Medium = Medium_h)
          "Outlet node hot side"      annotation (Placement(transformation(extent={{-10,88},{10,108}}),
              iconTransformation(extent={{-10,-10},{10,10}},
              rotation=-90,
              origin={0,-140})));
      equation

        //Mass balance
        node_c_in.m_flow     = node_c_out.m_flow;
        node_h_in.m_flow     = node_h_out.m_flow;

        //Energy balance
        qdot                 = node_c_in.m_flow*(node_c_out.h - node_c_in.h);
        qdot                 = node_h_in.m_flow*(node_h_in.h - node_h_out.h);

        //Outlet pressures
        node_c_out.p         = node_c_in.p - dp_c;
        node_h_out.p         = node_h_in.p - dp_h;

        //Boundary equations
        if use_power then
          qdot  = power;
        end if;

        if use_dT_in then
          dT_in              = Medium_h.temperature(Medium_h.setState_phX(node_h_in.p,
          node_h_in.h)) - Medium_c.temperature(Medium_c.setState_phX(node_c_out.p,
          node_c_out.h));
        end if;

        if use_dT_out then
          dT_out             = Medium_h.temperature(Medium_h.setState_phX(node_h_out.p,
          node_h_out.h)) - Medium_c.temperature(Medium_c.setState_phX(node_c_in.p,
          node_c_in.h));
        end if;

        connect(node_c_in, node_c_in) annotation (Line(
            points={{100,-60},{100,-60}},
            color={0,127,255},
            smooth=Smooth.None));
        annotation(Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
                  -100},{100,100}}),
                           graphics), Icon(coordinateSystem(extent={{-100,
                  -100},{100,100}},                                                              preserveAspectRatio=false,  initialScale = 0.1, grid = {2, 2}), graphics={
                                                                                                  Text(lineColor=
                    {0,0,0},                                                                                                    extent={{
                    -100,15},{100,-15}},                                                                                                    textString=  "%name",
                origin={-128,134},
                rotation=0),
              Rectangle(extent={{-100,100},{100,-100}}, lineColor={0,0,0},
                fillColor={255,255,255},
                fillPattern=FillPattern.Solid),
              Line(points={{44,28},{44,86}},   color={255,0,0}),
              Line(points={{44,28},{64,64}},   color={255,0,0}),
              Line(points={{44,28},{24,64}},   color={255,0,0}),
              Line(points={{-48,0},{-74,-16}}, color={0,128,255}),
              Line(points={{-74,16},{-48,0}}, color={0,128,255}),
              Line(
                points={{-140,0},{-40,0},{0,-60},{0,60},{40,0},{140,0}},
                color={0,128,255}),
              Line(
                points={{0,-140},{0,-100}},
                color={255,0,0},
                smooth=Smooth.None),
              Line(
                points={{0,142},{0,100}},
                color={255,0,0},
                smooth=Smooth.None)}));
      end Heat_exchanger;
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

    package Flags "Components to impose or initialize variables"

      model START
        "Flag to initialize a condition related to a process variable"
        replaceable package Medium = Test.Media.OneRandomOrganicFluid constrainedby
          Modelica.Media.Interfaces.PartialMedium "Medium model hot cells" annotation(choicesAllMatching = true);
        parameter Modelica.SIunits.AbsolutePressure p = 1e5 "Pressure" annotation (Dialog(tab="Start"));
        parameter Modelica.SIunits.Temperature T = 273.15 "Temperature" annotation (Dialog(tab="Start"));
        parameter Modelica.SIunits.MassFlowRate m_flow = 1 "Mass flow" annotation (Dialog(tab="Start"));
        Nodes.Node_in node(redeclare package Medium = Medium, p(start = p),
          h(start = Medium.specificEnthalpy(Medium.setState_pT(p,T))), m_flow(start = m_flow))
          annotation (Placement(transformation(extent={{90,-10},{110,10}}),
              iconTransformation(extent={{90,-10},{110,10}})));

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
        replaceable package Medium = Test.Media.OneRandomOrganicFluid constrainedby
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
        Nodes.Node_in node(redeclare package Medium = Medium)
          annotation (Placement(transformation(extent={{91,-10},{111,10}}),
              iconTransformation(extent={{91,-10},{111,10}})));
      equation
        //Component equations
        state = Medium.setState_phX(node.p, node.h);

        //Boundary equations
        if use_T then
          node.h =  Medium.specificEnthalpy(Medium.setState_pTX(node.p, T));
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
                    -112,27.5},{112,-27.5}},                                                                                                textString=  "%name",
                origin={-162,150.5},
                rotation=0)}));
      end ADDCO;

      model ADDCOW
        "Flag to impose and visualize the power (work, eletric. etc.)"
        parameter Boolean use_W = false "True if power is given"  annotation (Dialog(tab="Addco"));
        parameter Modelica.SIunits.Power W = 0 "Power production" annotation (Dialog(tab="Addco"));
        parameter Boolean visible = true "if true display mode on";
        CycleTempo.Components.Nodes.terminal terminal
          annotation (Placement(transformation(extent={{90,-10},{110,10}})));
      equation

        //Boundary equations
        if use_W then
          W =  terminal.W;
        end if;

        annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
                  -100},{100,100}}), graphics), Icon(coordinateSystem(
                preserveAspectRatio=false, extent={{-100,-100},{100,100}}), graphics={
              Rectangle(extent={{120,50},{-120,120}}, lineColor={0,0,0},
                fillColor={255,255,255},
                fillPattern=FillPattern.Solid),
              Text(visible=visible,
              extent={{10,26},{101,131}},
              textString=DynamicSelect("kW", String("kW"))),
              Text(visible=visible,
              extent={{-98,26},{-7,131}},
              textString=DynamicSelect("0.0", String(1e-3*terminal.W, significantDigits=2, format=  ".2f"))),
              Rectangle(
              extent={{-40,-50},{60,50}},
              lineColor={0,0,0},
              fillColor={0,128,255},
              fillPattern=FillPattern.Solid), Line(
              points={{101,0},{60,0}},
              color={0,0,0},
              smooth=Smooth.None),
              Text(
                extent={{-30,36},{48,-43}},
                lineColor={0,0,0},
                textString="AW"),                                                                 Text(lineColor=
                    {0,0,0},                                                                                                    extent={{
                    -112,27.5},{112,-27.5}},
                origin={0,146.5},
                rotation=0,
                textString="%name")}));
      end ADDCOW;

      model CC "A component to close the cycle and avoid circular equalities"
        replaceable package Medium = Test.Media.OneRandomOrganicFluid constrainedby
          Modelica.Media.Interfaces.PartialMedium "Medium model hot cells" annotation(choicesAllMatching = true);

        Nodes.Node_in node_in(redeclare package Medium = Medium)
          annotation (Placement(transformation(extent={{-68,-10},{-48,10}}),
              iconTransformation(extent={{-68,-10},{-48,10}})));
        Nodes.Node_out node_out(redeclare package Medium = Medium)
          annotation (Placement(transformation(extent={{50,-10},{70,10}}),
              iconTransformation(extent={{50,-10},{70,10}})));
      equation
        //Component equations
        node_in.h = node_out.h;
        node_in.p = node_out.p;
        annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
                  -100},{100,100}}), graphics), Icon(coordinateSystem(
                preserveAspectRatio=false, extent={{-100,-100},{100,100}}), graphics={Ellipse(
                extent={{-60,60},{60,-60}},
                lineColor={255,0,0},
                fillColor={255,0,0},
                fillPattern=FillPattern.Sphere)}));
      end CC;
    end Flags;

    package Electrics "Package containing electric components"
      model Generator "Electric generator"
        parameter Real eta_el "Electric efficiency";
        Nodes.terminal terminal_in "Inlet terminal"
          annotation (Placement(transformation(extent={{-108,-10},{-88,10}})));
        Nodes.terminal terminal_out "Outlet terminal"
          annotation (Placement(transformation(extent={{90,-10},{110,10}})));
      equation
        terminal_out.W = terminal_in.W*eta_el;
        annotation(Documentation(info = "<HTML>
                                                                                              <p>This is a cycle model for pumps.
                                                                                              </html>"),  Icon(coordinateSystem(extent={{-100,
                  -100},{100,100}},                                                                                                    preserveAspectRatio=false,  initialScale = 0.1, grid = {2, 2}), graphics={                       Ellipse(
                                                                                                  extent={{
                    -60,60},{60,-60}},                                                                                                    endAngle=  360,
                lineColor={0,0,0},
                fillColor={255,255,255},
                fillPattern=FillPattern.Solid),                                                                                   Text(lineColor=
                    {0,0,0},                                                                                                    extent={{
                    -99,24},{99,-24}},
                origin={1,2},
                rotation=0,
                textString="G"),
              Line(
                points={{60,0},{100,0}},
                color={0,0,0},
                smooth=Smooth.None),
              Line(
                points={{-100,0},{-60,0}},
                color={0,0,0},
                smooth=Smooth.None)}),
          Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
                  -100},{100,100}}),
                          graphics));
      end Generator;

      model Motor "Electric motor"
        extends CycleTempo.Components.Electrics.Generator;
        annotation(Documentation(info = "<HTML>
                                                                                              <p>This is a cycle model for pumps.
                                                                                              </html>"),  Icon(coordinateSystem(extent={{-100,
                  -100},{100,100}},                                                                                                    preserveAspectRatio=false,  initialScale = 0.1, grid = {2, 2}), graphics={                       Ellipse(
                                                                                                  extent={{
                    -60,60},{60,-60}},                                                                                                    endAngle=  360,
                lineColor={0,0,0},
                fillColor={255,255,255},
                fillPattern=FillPattern.Solid),                                                                                   Text(lineColor=
                    {0,0,0},                                                                                                    extent={{
                    -99,24},{99,-24}},
                origin={1,2},
                rotation=0,
                textString="M"),
              Line(
                points={{60,0},{100,0}},
                color={0,0,0},
                smooth=Smooth.None),
              Line(
                points={{-100,0},{-60,0}},
                color={0,0,0},
                smooth=Smooth.None)}),
          Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,
                  100}}), graphics));
      end Motor;
    end Electrics;

    package Flow "Flow model"
      model Mixer "Mixer model"
        //Medium model
        replaceable package Medium = Test.Media.OneRandomOrganicFluid constrainedby
          Modelica.Media.Interfaces.PartialMedium "Medium model" annotation(choicesAllMatching = true);
        Nodes.Node_in node_in_1(redeclare package Medium = Medium)
          "Inlet node 1"                                                          annotation (Placement(transformation(extent={{50,-110},
                  {70,-90}}),
              iconTransformation(
              extent={{-10,-10},{10,10}},
              rotation=-90,
              origin={40,-100})));
        Nodes.Node_in node_in_2(redeclare package Medium = Medium)
          "Inlet node 2"                                                          annotation (Placement(transformation(extent={{-70,-110},
                  {-50,-90}}),
              iconTransformation(
              extent={{-10,-10},{10,10}},
              rotation=-90,
              origin={-40,-100})));
        Nodes.Node_out node_out(redeclare package Medium = Medium)
          "Outlet node"                                                          annotation (Placement(transformation(extent={{-10,-10},
                  {10,10}}), iconTransformation(
              extent={{-10,-10},{10,10}},
              rotation=-90,
              origin={0,20})));
      equation
        node_out.p      = node_in_1.p;
        node_out.m_flow = node_in_1.m_flow + node_in_2.m_flow;
        node_out.h      = (node_in_1.m_flow*node_in_1.h +
        node_in_2.m_flow*node_in_2.h)/node_out.m_flow;
        connect(node_out, node_out) annotation (Line(
            points={{0,0},{4,0},{4,0},{0,0}},
            color={0,0,0},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
                  -100},{100,100}}), graphics), Icon(coordinateSystem(
                preserveAspectRatio=false, extent={{-100,-100},{100,100}}), graphics={
              Line(
                points={{-50,-100},{-50,10}},
                color={0,0,0},
                smooth=Smooth.None),
              Line(
                points={{0,30},{-50,10}},
                color={0,0,0},
                smooth=Smooth.None),
              Line(
                points={{50,-100},{50,10}},
                color={0,0,0},
                smooth=Smooth.None),
              Line(
                points={{-30,-100},{-30,-10}},
                color={0,0,0},
                smooth=Smooth.None),
              Line(
                points={{0,10},{-30,-10}},
                color={0,0,0},
                smooth=Smooth.None),
              Line(
                points={{30,-100},{30,-10}},
                color={0,0,0},
                smooth=Smooth.None),
              Line(
                points={{0,10},{30,-10}},
                color={0,0,0},
                smooth=Smooth.None),
              Line(
                points={{0,30},{50,10}},
                color={0,0,0},
                smooth=Smooth.None),
              Line(points={{40,-40},{40,-76}}, color={0,0,0}),
              Line(points={{40,-40},{48,-54}}, color={0,0,0}),
              Line(points={{40,-40},{32,-54}}, color={0,0,0}),
              Line(points={{-40,-40},{-40,-76}},
                                               color={0,0,0}),
              Line(points={{-40,-40},{-32,-54}},
                                               color={0,0,0}),
              Line(points={{-40,-40},{-48,-54}},
                                               color={0,0,0})}));
      end Mixer;

      model Splitter "Splitter model"
        //Medium model
        replaceable package Medium = Test.Media.OneRandomOrganicFluid constrainedby
          Modelica.Media.Interfaces.PartialMedium "Medium model" annotation(choicesAllMatching = true);
        Nodes.Node_out node_out_1(redeclare package Medium = Medium)
          "Outlet node 1"
          annotation (Placement(transformation(extent={{50,-110},{70,-90}}),
              iconTransformation(
              extent={{-10,-10},{10,10}},
              rotation=-90,
              origin={40,-100})));
        Nodes.Node_out node_out_2(redeclare package Medium = Medium)
          "Outlet node 2"
          annotation (Placement(transformation(extent={{-70,-110},{-50,-90}}),
              iconTransformation(
              extent={{-10,-10},{10,10}},
              rotation=-90,
              origin={-40,-100})));
        Nodes.Node_in node_in(redeclare package Medium = Medium) "Inlet node"
          annotation (Placement(transformation(extent={{-10,-10},{10,10}}),
              iconTransformation(
              extent={{-10,-10},{10,10}},
              rotation=-90,
              origin={0,20})));
      equation
        node_in.m_flow = node_out_1.m_flow + node_out_2.m_flow;
        node_in.p      = node_out_1.p;
        node_in.p      = node_out_2.p;
        node_in.h      = node_out_1.h;
        node_in.h      = node_out_2.h;

        connect(node_in, node_in) annotation (Line(
            points={{0,0},{0,0}},
            color={0,0,0},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
                  -100},{100,100}}), graphics), Icon(coordinateSystem(
                preserveAspectRatio=false, extent={{-100,-100},{100,100}}), graphics={
              Line(
                points={{-50,-100},{-50,10}},
                color={0,0,0},
                smooth=Smooth.None),
              Line(
                points={{0,30},{-50,10}},
                color={0,0,0},
                smooth=Smooth.None),
              Line(
                points={{50,-100},{50,10}},
                color={0,0,0},
                smooth=Smooth.None),
              Line(
                points={{-30,-100},{-30,-10}},
                color={0,0,0},
                smooth=Smooth.None),
              Line(
                points={{0,10},{-30,-10}},
                color={0,0,0},
                smooth=Smooth.None),
              Line(
                points={{30,-100},{30,-10}},
                color={0,0,0},
                smooth=Smooth.None),
              Line(
                points={{0,10},{30,-10}},
                color={0,0,0},
                smooth=Smooth.None),
              Line(
                points={{0,30},{50,10}},
                color={0,0,0},
                smooth=Smooth.None),
              Line(points={{40,-40},{40,-76}}, color={0,0,0}),
              Line(points={{40,-76},{48,-62}}, color={0,0,0}),
              Line(points={{40,-76},{32,-62}}, color={0,0,0}),
              Line(points={{-40,-40},{-40,-76}},
                                               color={0,0,0}),
              Line(points={{-40,-76},{-32,-62}},
                                               color={0,0,0}),
              Line(points={{-40,-76},{-48,-62}},
                                               color={0,0,0})}));
      end Splitter;
    end Flow;

    package Mechanics "Package containing mechanic components"

      model GearBox "Gear box model"
        parameter Real eta_m "Mechanical efficiency";
        Nodes.terminal terminal_in "Inlet terminal"
          annotation (Placement(transformation(extent={{-108,-10},{-88,10}}),
              iconTransformation(extent={{-110,0},{-90,20}})));
        Nodes.terminal terminal_out "Outlet terminal"
          annotation (Placement(transformation(extent={{90,-10},{110,10}}),
              iconTransformation(extent={{110,0},{130,20}})));
      equation
        terminal_out.W = terminal_in.W*eta_m;
        annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},
                  {100,100}}), graphics={
            Rectangle(origin={30,10},
              lineColor={64,64,64},
              fillColor={192,192,192},
              fillPattern=FillPattern.HorizontalCylinder,
              extent={{-120.0,-10.0},{-80.0,10.0}}),
            Polygon(fillColor={192,192,192},
              fillPattern=FillPattern.HorizontalCylinder,
              points={{-50,20},{-50,30},{-30,50},{-30,-30},{-50,-10},{-50,20}}),
            Rectangle(lineColor={64,64,64},
              fillColor={255,255,255},
              fillPattern=FillPattern.HorizontalCylinder,
              extent={{-30,-50},{50,70}},
              radius=10.0),
            Rectangle(lineColor={64,64,64},
              fillPattern=FillPattern.None,
              extent={{-30,-50},{50,70}},
              radius=10.0),
            Polygon(fillColor={192,192,192},
              fillPattern=FillPattern.HorizontalCylinder,
              points={{70,30},{50,50},{50,-30},{70,-10},{70,30}}),
            Rectangle(origin={-10,10},
              lineColor={64,64,64},
              fillColor={192,192,192},
              fillPattern=FillPattern.HorizontalCylinder,
              extent={{80.0,-10.0},{120.0,10.0}}),
            Polygon(origin={10,20},
              fillColor={64,64,64},
              fillPattern=FillPattern.Solid,
              points={{-60.0,-90.0},{-50.0,-90.0},{-20.0,-30.0},{20.0,-30.0},{48.0,-90.0},{60.0,-90.0},{60.0,-100.0},{-60.0,-100.0},{-60.0,-90.0}}),
                                        Text(
                    extent={{-140,160},{160,120}},
                    lineColor={0,0,0},
                textString="%name")}));
      end GearBox;
    end Mechanics;
  end Components;
  annotation (uses(Modelica(version="3.2.1")));
end CycleTempo;
