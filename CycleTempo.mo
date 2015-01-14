within ;
package CycleTempo "Cycle Tempo 2_0"
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

      connector power "A node with a thermodynamic state"
        Modelica.SIunits.Power W "Power";
      end power;

      connector terminal "Electric terminal"
        extends power;
        annotation(Diagram(coordinateSystem(preserveAspectRatio = true, extent = {{-100, -100}, {100, 100}}), graphics={         Text(extent={{
                    -98,202},{102,108}},                                                                                                    textString = "%name", lineColor=
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
        Modelica.SIunits.Power W "Power consumption";
        Modelica.SIunits.SpecificEnthalpy h_is
          "Outlet isentropic specific enthalpy";
        Modelica.SIunits.SpecificEntropy s "Inlet specific entropy";
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
        s              = Medium.specificEntropy_ph(node_in.p, node_in.h);
        h_is           = Medium.specificEnthalpy_ps(node_out.p, s);
        node_out.h     = (h_is - node_in.h)/eta_is + node_in.h;

        annotation(Documentation(info = "<HTML>
                                                                                              <p>This is a cycle model for pumps.
                                                                                              </html>"),  Icon(coordinateSystem(extent={{-100,
                  -100},{100,100}},                                                                                                    preserveAspectRatio=false,  initialScale = 0.1, grid = {2, 2}), graphics={                       Ellipse(
                                                                                                  extent={{
                    -60,60},{60,-60}},                                                                                                    endAngle = 360,
                lineColor={0,0,0},
                fillColor={255,255,255},
                fillPattern=FillPattern.Solid),                                                                                                    Polygon(origin={
                    -20.716,0.7156},                                                                                                    rotation = 180, fillColor=
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
        Modelica.SIunits.Power W "Power production";
        Modelica.SIunits.SpecificEnthalpy h_is
          "Outlet isentropic specific enthalpy";
        Modelica.SIunits.SpecificEntropy s "Inlet specific entropy";

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
        s              = Medium.specificEntropy_ph(node_in.p, node_in.h);
        h_is           = Medium.specificEnthalpy_ps(node_out.p, s);
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
                    -100,15},{100,-15}},                                                                                                    textString = "%name",
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
    end Turbomachinery;

    package Sources "Soruces and sinks"
      model Source
        replaceable package Medium = Test.Media.OneRandomOrganicFluid constrainedby
          Modelica.Media.Interfaces.PartialMedium "Medium model hot cells" annotation(choicesAllMatching = true);

        parameter Modelica.SIunits.MassFlowRate mdot "Mass flow rate" annotation(Dialog(tab="Addco"));
        parameter Modelica.SIunits.Temperature T "Temperature" annotation(Dialog(tab="Addco"));
        parameter Modelica.SIunits.AbsolutePressure p "Pressure" annotation(Dialog(tab="Addco"));
        parameter Modelica.SIunits.SpecificEnthalpy h "Enthalpy" annotation(Dialog(tab="Addco"));
        parameter Boolean use_mdot = false "True if mass flow is given" annotation(Dialog(tab="Addco"));
        parameter Boolean use_T = false "True if temperature is given" annotation(Dialog(tab="Addco"));
        parameter Boolean use_p = false "True if pressure is given" annotation(Dialog(tab="Addco"));
        parameter Boolean use_h = false "True if enthalpy is given" annotation(Dialog(tab="Addco"));
        Nodes.Node_out node_out(redeclare package Medium = Medium)
          annotation (Placement(transformation(extent={{90,-10},{110,10}})));

      equation
        //Boundary equations
        if use_mdot then
          node_out.m_flow = mdot;
        end if;
        if use_T then
          T               = Medium.temperature(Medium.setState_ph(node_out.p,
          node_out.h));
        end if;
        if use_p then
          node_out.p      = p;
        end if;
        if use_h then
          node_out.h      = h;
        end if;
        annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},
                  {100,100}}), graphics={                                                                                                    Polygon(lineColor=  {255, 255, 255}, fillColor=  {255, 255, 255},
                  fillPattern=                                                                                                    FillPattern.Solid, points={{
                    20,-75},{50,-85},{20,-95},{20,-75}}),                                                              Rectangle(extent={{
                    20,60},{100,-60}},                                                                                                    lineColor=  {0, 0, 0},
                  fillPattern=                                                                                                    FillPattern.HorizontalCylinder, fillColor=  {192, 192, 192}), Rectangle(extent={{
                    38,40},{100,-40}},                                                                                                    lineColor=  {0, 0, 0},
                  fillPattern=                                                                                                    FillPattern.HorizontalCylinder, fillColor=  {0, 127, 255}), Ellipse(extent={{
                    -100,80},{60,-80}},                                                                                                    fillColor=  {255, 255, 255},
                  fillPattern=                                                                                                    FillPattern.Solid, lineColor=  {0, 0, 255}), Polygon(points={{
                    -60,70},{60,0},{-60,-68},{-60,70}},                                                                                                    lineColor=  {0, 0, 255}, fillColor=  {0, 0, 255},
                  fillPattern=                                                                                                    FillPattern.Solid), Text(extent={{
                    -54,32},{16,-30}},                                                                                                    lineColor=  {255, 0, 0}, textString=  "m"), Ellipse(extent={{
                    -26,30},{-18,22}},                                                                                                    lineColor=  {255, 0, 0}, fillColor=  {255, 0, 0},
                  fillPattern=                                                                                                    FillPattern.Solid)}));
      end Source;

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

      model Sink
        replaceable package Medium = Test.Media.OneRandomOrganicFluid constrainedby
          Modelica.Media.Interfaces.PartialMedium "Medium model hot cells" annotation(choicesAllMatching = true);

        parameter Modelica.SIunits.MassFlowRate mdot "Mass flow rate" annotation(Dialog(tab="Addco"));
        parameter Modelica.SIunits.Temperature T "Temperature" annotation(Dialog(tab="Addco"));
        parameter Modelica.SIunits.AbsolutePressure p "Pressure" annotation(Dialog(tab="Addco"));
        parameter Modelica.SIunits.SpecificEnthalpy h "Enthalpy" annotation(Dialog(tab="Addco"));
        parameter Boolean use_mdot = false "True if mass flow is given" annotation(Dialog(tab="Addco"));
        parameter Boolean use_T = false "True if temperature is given" annotation(Dialog(tab="Addco"));
        parameter Boolean use_p = false "True if pressure is given" annotation(Dialog(tab="Addco"));
        parameter Boolean use_h = false "True if enthalpy is given" annotation(Dialog(tab="Addco"));
        Nodes.Node_in node_in(redeclare package Medium = Medium)
          annotation (Placement(transformation(extent={{90,-10},{110,10}})));

      equation
        //Boundary equations
        if use_mdot then
          node_in.m_flow = mdot;
        end if;
        if use_T then
          T               = Medium.temperature(Medium.setState_ph(node_in.p,
          node_in.h));
        end if;
        if use_p then
          node_in.p      = p;
        end if;
        if use_h then
          node_in.h      = h;
        end if;

        annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},
                  {100,100}}), graphics={                                                                                                    Polygon(lineColor=  {255, 255, 255}, fillColor=  {255, 255, 255},
                  fillPattern=                                                                                                    FillPattern.Solid, points={{
                    20,-75},{50,-85},{20,-95},{20,-75}}),                                                              Rectangle(extent={{
                    20,60},{100,-60}},                                                                                                    lineColor=  {0, 0, 0},
                  fillPattern=                                                                                                    FillPattern.HorizontalCylinder, fillColor=  {192, 192, 192}), Rectangle(extent={{
                    38,40},{100,-40}},                                                                                                    lineColor=  {0, 0, 0},
                  fillPattern=                                                                                                    FillPattern.HorizontalCylinder, fillColor=  {0, 127, 255}), Ellipse(extent={{
                    -100,80},{60,-80}},                                                                                                    fillColor=  {255, 255, 255},
                  fillPattern=                                                                                                    FillPattern.Solid, lineColor=  {0, 0, 255}), Polygon(points={{
                    -60,70},{60,0},{-60,-68},{-60,70}},                                                                                                    lineColor=  {0, 0, 255}, fillColor=  {0, 0, 255},
                  fillPattern=                                                                                                    FillPattern.Solid), Text(extent={{
                    -54,32},{16,-30}},                                                                                                    lineColor=  {255, 0, 0}, textString=  "m"), Ellipse(extent={{
                    -26,30},{-18,22}},                                                                                                    lineColor=  {255, 0, 0}, fillColor=  {255, 0, 0},
                  fillPattern=                                                                                                    FillPattern.Solid)}));
      end Sink;
    end Sources;

    package HEX "Heat exchangers"
      model Evaporator
        "Evaporator with three zones preheater, evaporator, superheater"
        extends CycleTempo.Components.HEX.Heat_exchanger;
        parameter Modelica.SIunits.TemperatureDifference dT_int = 0
          "Temperature difference at the inlet of the evaporator" annotation (Dialog(tab="Addco"));
        parameter Boolean use_dT_int = false
          "true if temperature difference at the inlet of the evaporator is given";
        parameter Modelica.SIunits.SpecificEnthalpy hh_in_start=
        Medium_h.specificEnthalpy_pT(1e5, 500)
          "Inlet enthalpy evaporator hot side start value" annotation (Dialog(tab="Start"));
        Medium_c.SaturationProperties sat "Saturated state hot fluid";
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
        dT_eva           = Medium_h.temperature(Medium_h.setState_ph(node_h_in.p,
                           hh_in_pre)) - sat.Tsat;
        qdot_pre         = node_c_in.m_flow*(sat.hl - node_c_in.h);
        qdot_eva         = node_c_in.m_flow*(sat.hv - sat.hl);
        qdot_sup         = node_c_in.m_flow*(node_c_out.h - sat.hv);
        qdot_pre         = node_h_in.m_flow*(hh_in_pre - node_h_out.h);
        qdot_eva         = node_h_in.m_flow*(hh_in_eva - hh_in_pre);
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
        parameter Modelica.SIunits.SpecificEnthalpy hc_in_start=
        Medium_c.specificEnthalpy_pT(1e5, 298.15)
          "Inlet enthalpy condensation region cold side start value" annotation (Dialog(tab="Start"));
        Modelica.SIunits.SpecificEnthalpy hc_in_con( start = hc_in_start)
          "Inlet enthalpy condensation region cold side";
        Modelica.SIunits.HeatFlowRate qdot_con "Heat rate condensation region";
        Modelica.SIunits.TemperatureDifference dT_con
          "Temperature difference at the inlet of the condensation region";
        Medium_h.SaturationProperties sat "Saturated state hot fluid";
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
        dT_con       = sat.Tsat - Medium_c.temperature(Medium_c.setState_ph(node_c_in.p,
                           hc_in_con));
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
          dT_in              = Medium_h.temperature(Medium_h.setState_ph(node_h_in.p,
          node_h_in.h)) - Medium_c.temperature(Medium_c.setState_ph(node_c_out.p,
          node_c_out.h));
        end if;

        if use_dT_out then
          dT_out             = Medium_h.temperature(Medium_h.setState_ph(node_h_out.p,
          node_h_out.h)) - Medium_c.temperature(Medium_c.setState_ph(node_c_in.p,
          node_c_in.h));
        end if;

        connect(node_c_in, node_c_in) annotation (Line(
            points={{100,-60},{100,-60}},
            color={0,127,255},
            smooth=Smooth.None));
        annotation(Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},
                  {100,100}}),
                           graphics), Icon(coordinateSystem(extent={{-100,
                  -100},{100,100}},                                                              preserveAspectRatio=false,  initialScale = 0.1, grid = {2, 2}), graphics={
                                                                                                  Text(lineColor=
                    {0,0,0},                                                                                                    extent={{
                    -100,15},{100,-15}},                                                                                                    textString = "%name",
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
        parameter Modelica.SIunits.SpecificEnthalpy h = 1e5 "Enthalpy" annotation (Dialog(tab="Start"));
        parameter Modelica.SIunits.MassFlowRate m_flow = 1 "Mass flow" annotation (Dialog(tab="Start"));
        Nodes.Node_in node(redeclare package Medium = Medium, p(start = p), h(start = h), m_flow(start = m_flow))
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
                    -60,60},{60,-60}},                                                                                                    endAngle = 360,
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
                    -60,60},{60,-60}},                                                                                                    endAngle = 360,
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
  end Components;
  annotation (uses(Modelica(version="3.2.1")));
end CycleTempo;
