within VIP;
package CycleTempo "Cycle Tempo 2_0"
  package Components "Component library"
    package Nodes "Node package"
      extends Modelica.Icons.Package;

      connector hp_mdot "A node with a thermodynamic state"
        replaceable package Medium = VIP.Media.OneRandomOrganicFluid  constrainedby
          Modelica.Media.Interfaces.PartialMedium "Medium model hot cells" annotation(choicesAllMatching = true);
        Medium.MassFlowRate m_flow;
        Modelica.SIunits.SpecificEnthalpy h "Specific enthalpy";
        Modelica.SIunits.AbsolutePressure p "Pressure";
      end hp_mdot;

      connector Node_in "Fluid connector with filled icon"
        extends hp_mdot;
        annotation(Diagram(coordinateSystem(preserveAspectRatio = true, extent = {{-100, -100}, {100, 100}}), graphics={  Ellipse(extent = {{-100, 100}, {100, -100}}, lineColor = {0, 127, 255}, fillColor = {0, 127, 255},
                  fillPattern =                                                                                                    FillPattern.Solid), Ellipse(extent = {{-100, 100}, {100, -100}}, lineColor = {0, 0, 0}, fillColor = {0, 127, 255},
                  fillPattern =                                                                                                    FillPattern.Solid), Text(extent = {{-88, 206}, {112, 112}}, textString = "%name", lineColor = {0, 0, 255})}), Icon(coordinateSystem(preserveAspectRatio = true, extent = {{-100, -100}, {100, 100}}), graphics={  Ellipse(extent = {{-100, 100}, {100, -100}}, lineColor = {0, 127, 255}, fillColor = {0, 127, 255},
                  fillPattern =                                                                                                    FillPattern.Solid), Ellipse(extent = {{-100, 100}, {100, -100}}, lineColor = {0, 0, 0}, fillColor = {0, 127, 255},
                  fillPattern =                                                                                                    FillPattern.Solid)}), Documentation(info = "<html>Modelica.Media.Examples.Tests.Components.FluidPort_a
                                                                                                </html>"));
      end Node_in;

      connector Node_out "Fluid connector with outlined icon"
        extends hp_mdot;
        annotation(Diagram(coordinateSystem(preserveAspectRatio = true, extent = {{-100, -100}, {100, 100}}), graphics={  Ellipse(extent = {{-100, 100}, {100, -100}}, lineColor = {0, 127, 255}, fillColor = {0, 127, 255},
                  fillPattern =                                                                                                    FillPattern.Solid), Ellipse(extent = {{-100, 100}, {100, -100}}, lineColor = {0, 0, 0}, fillColor = {0, 127, 255},
                  fillPattern =                                                                                                    FillPattern.Solid), Ellipse(extent = {{-80, 80}, {80, -80}}, lineColor = {0, 127, 255}, fillColor = {255, 255, 255},
                  fillPattern =                                                                                                    FillPattern.Solid), Text(extent = {{-88, 192}, {112, 98}}, textString = "%name", lineColor = {0, 0, 255})}), Icon(coordinateSystem(preserveAspectRatio = true, extent = {{-100, -100}, {100, 100}}), graphics={  Ellipse(extent = {{-100, 100}, {100, -100}}, lineColor = {0, 127, 255}, fillColor = {0, 127, 255},
                  fillPattern =                                                                                                    FillPattern.Solid), Ellipse(extent = {{-100, 100}, {100, -100}}, lineColor = {0, 0, 0}, fillColor = {0, 127, 255},
                  fillPattern =                                                                                                    FillPattern.Solid), Ellipse(extent = {{-80, 80}, {80, -80}}, lineColor = {0, 127, 255}, fillColor = {255, 255, 255},
                  fillPattern =                                                                                                    FillPattern.Solid)}), Documentation(info = "<html>

                                                                                                </html>"));
      end Node_out;

    end Nodes;

    package Turbomachinery "Pump and turbine models"
      model Pump
        replaceable package Medium = VIP.Media.OneRandomOrganicFluid  constrainedby
          Modelica.Media.Interfaces.PartialMedium "Medium model hot cells" annotation(choicesAllMatching = true);
        parameter Real eta_is "Isentropic efficiency";
        parameter Real eta_m "Mechanical efficiency";
        parameter Real eta_el "Electric efficiency";
        parameter Modelica.SIunits.AbsolutePressure p_out = 1e5
          "Outlet pressure"
                           annotation (Dialog(tab="Addco"));
        parameter Modelica.SIunits.Temperature T_out_start = 298.15
          "Outlet temperature start value" annotation (Dialog(tab="Start"));
        parameter Modelica.SIunits.AbsolutePressure p_out_start = 1e5
          "Outlet pressure start value" annotation (Dialog(tab="Start"));
        parameter Medium.SpecificEnthalpy h_out_start = Medium.specificEnthalpy_pT(
        p_out_start, T_out_start) "Outlet enthalpy start value" annotation (Dialog(tab="Start"));
        parameter Boolean use_p_out = false "True if pressure is given"  annotation (Dialog(tab="Addco"));
        Modelica.SIunits.Power W "Power consumption";
        Modelica.SIunits.SpecificEnthalpy h_is
          "Outlet isentropic specific enthalpy";
        Modelica.SIunits.SpecificEntropy s "Inlet specific entropy";
        Nodes.Node_in node_in(redeclare package Medium = Medium) "Inlet node"
          annotation (Placement(transformation(extent={{-108,-10},{-88,10}}),
              iconTransformation(extent={{-108,-10},{-88,10}})));
        Nodes.Node_out node_out(redeclare package Medium = Medium, h(start=
        h_out_start), p(start = p_out_start)) "Outlet node"
          annotation (Placement(transformation(extent={{90,-10},{110,10}}),
              iconTransformation(extent={{90,-10},{110,10}})));
      equation

        //Mass balance
        node_in.m_flow = node_out.m_flow;

        //Energy balance
        W              = node_in.m_flow*(node_out.h - node_in.h)/(eta_m*eta_el);

        //Boundary equations
        if use_p_out then
          node_out.p   = p_out;
        end if;

        //Component equations
        s              = Medium.specificEntropy_ph(node_in.p, node_in.h);
        h_is           = Medium.specificEnthalpy_ps(node_out.p, s);
        node_out.h     = (h_is - node_in.h)/eta_is + node_in.h;

        annotation(Documentation(info = "<HTML>
                                                                                              <p>This is a cycle model for pumps.
                                                                                              </html>"),  Icon(coordinateSystem(extent={{-100,
                  -100},{100,100}},                                                                                                    preserveAspectRatio=false,  initialScale = 0.1, grid = {2, 2}), graphics={  Rectangle(fillColor = {0, 127, 255},
                  fillPattern =                                                                                                    FillPattern.HorizontalCylinder, extent = {{-100, 46}, {100, -46}}), Polygon(lineColor = {0, 0, 255}, pattern = LinePattern.None,
                  fillPattern =                                                                                                    FillPattern.VerticalCylinder, points = {{-48, -60}, {-72, -100}, {72, -100}, {48, -60}, {-48, -60}}), Ellipse(fillColor = {0, 100, 199},
                  fillPattern =                                                                                                    FillPattern.Sphere, extent = {{-80, 80}, {80, -80}}, endAngle = 360), Polygon(origin = {-0.715564, 0.715564}, rotation = 180, fillColor = {255, 255, 255}, pattern = LinePattern.None,
                  fillPattern =                                                                                                    FillPattern.HorizontalCylinder, points={{
                    28,30},{28,-30},{-50,-2},{28,30}})}),
          Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,
                  100}}), graphics));
      end Pump;

      model Turbine "Turbine model"

        replaceable package Medium = VIP.Media.OneRandomOrganicFluid  constrainedby
          Modelica.Media.Interfaces.PartialMedium "Medium model hot cells" annotation(choicesAllMatching = true);
        parameter Real eta_is "Isentropic efficiency";
        parameter Real eta_m "Mechanical efficiency";
        parameter Modelica.SIunits.Temperature T_out = 298.15
          "Outlet temperature"                                                            annotation (Dialog(tab="Addco"));
        parameter Modelica.SIunits.AbsolutePressure p_out = 1e5
          "Outlet pressure" annotation (Dialog(tab="Addco"));
        parameter Boolean use_T_out = false "True if temperature is given"  annotation (Dialog(tab="Addco"));
        parameter Boolean use_p_out = false "True if pressure is given"  annotation (Dialog(tab="Addco"));
        parameter Modelica.SIunits.Temperature T_out_start = 298.15
          "Outlet temperature start value" annotation (Dialog(tab="Start"));
        parameter Modelica.SIunits.AbsolutePressure p_out_start = 1e5
          "Outlet pressure start value" annotation (Dialog(tab="Start"));
        parameter Medium.SpecificEnthalpy h_out_start = Medium.specificEnthalpy_pT(
        p_out_start, T_out_start) "Outlet enthalpy start value" annotation (Dialog(tab="Start"));
        Modelica.SIunits.Power W "Power production";
        Modelica.SIunits.SpecificEnthalpy h_is
          "Outlet isentropic specific enthalpy";
        Modelica.SIunits.SpecificEntropy s "Inlet specific entropy";

        Nodes.Node_in node_in(redeclare package Medium = Medium) "Inlet node"
          annotation (Placement(transformation(extent={{-110,70},{-90,90}}),
              iconTransformation(extent={{-110,70},{-90,90}})));
        Nodes.Node_out node_out(redeclare package Medium = Medium, h(start=
        h_out_start), p(start = p_out_start)) "Outlet node"
          annotation (Placement(transformation(extent={{90,70},{110,90}}),
              iconTransformation(extent={{90,70},{110,90}})));
      equation

        //Mass balance
        node_in.m_flow = node_out.m_flow;

        //Energy balance
        W              = node_in.m_flow*(node_in.h - node_out.h)*eta_m;

        //Boundary equations
        if use_T_out then
          T_out        = Medium.temperature_ph(node_out.p, node_out.h);
        end if;
        if use_p_out then
          node_out.p   = p_out;
        end if;

        //Component equations
        s              = Medium.specificEntropy_ph(node_in.p, node_in.h);
        h_is           = Medium.specificEnthalpy_ps(node_out.p, s);
        node_out.h     = node_in.h - (node_in.h - h_is)*eta_is;
        annotation(Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},
                  {100,100}}),
                        graphics={  Polygon(points={{-28,76},{-28,28},{-22,28},{-22,82},
                    {-90,82},{-90,76},{-28,76}},                                                                                    lineColor = {0, 0, 0},
                  lineThickness =                                                                                                    0.5, fillColor = {0, 0, 255},
                  fillPattern =                                                                                                    FillPattern.Solid), Polygon(points={{
                    26,56},{32,56},{32,76},{90,76},{90,82},{26,82},{26,56}},                                                                                                    lineColor = {0, 0, 0},
                  lineThickness =                                                                                                    0.5, fillColor = {0, 0, 255},
                  fillPattern =                                                                                                    FillPattern.Solid), Rectangle(extent = {{-60, 8}, {60, -8}}, lineColor = {0, 0, 0},
                  fillPattern =                                                                                                    FillPattern.Sphere, fillColor = {160, 160, 164}), Polygon(points = {{-28, 28}, {-28, -26}, {32, -60}, {32, 60}, {-28, 28}}, lineColor = {0, 0, 0},
                  lineThickness =                                                                                                    0.5, fillColor = {0, 0, 255},
                  fillPattern =                                                                                                    FillPattern.Solid)}), Diagram(graphics), Documentation(info = "<html>
                                                                                                   <p>This is a cycle model for turbine
                                                                                                   </html>"));
      end Turbine;
    end Turbomachinery;

    package Sources "Soruces and sinks"
      model Source
        replaceable package Medium = VIP.Media.OneRandomOrganicFluid  constrainedby
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
                  {100,100}}), graphics={                                                                                                    Polygon(lineColor = {255, 255, 255}, fillColor = {255, 255, 255},
                  fillPattern =                                                                                                    FillPattern.Solid, points={{
                    20,-75},{50,-85},{20,-95},{20,-75}}),                                                              Rectangle(extent={{
                    20,60},{100,-60}},                                                                                                    lineColor = {0, 0, 0},
                  fillPattern =                                                                                                    FillPattern.HorizontalCylinder, fillColor = {192, 192, 192}), Rectangle(extent={{
                    38,40},{100,-40}},                                                                                                    lineColor = {0, 0, 0},
                  fillPattern =                                                                                                    FillPattern.HorizontalCylinder, fillColor = {0, 127, 255}), Ellipse(extent={{
                    -100,80},{60,-80}},                                                                                                    fillColor = {255, 255, 255},
                  fillPattern =                                                                                                    FillPattern.Solid, lineColor = {0, 0, 255}), Polygon(points={{
                    -60,70},{60,0},{-60,-68},{-60,70}},                                                                                                    lineColor = {0, 0, 255}, fillColor = {0, 0, 255},
                  fillPattern =                                                                                                    FillPattern.Solid), Text(extent={{
                    -54,32},{16,-30}},                                                                                                    lineColor = {255, 0, 0}, textString = "m"), Ellipse(extent={{
                    -26,30},{-18,22}},                                                                                                    lineColor = {255, 0, 0}, fillColor = {255, 0, 0},
                  fillPattern =                                                                                                    FillPattern.Solid)}));
      end Source;

      model CC "A component to close the cycle and avoid circular equalities"
        replaceable package Medium = VIP.Media.OneRandomOrganicFluid  constrainedby
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
        replaceable package Medium = VIP.Media.OneRandomOrganicFluid  constrainedby
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
                  {100,100}}), graphics={                                                                                                    Polygon(lineColor = {255, 255, 255}, fillColor = {255, 255, 255},
                  fillPattern =                                                                                                    FillPattern.Solid, points={{
                    20,-75},{50,-85},{20,-95},{20,-75}}),                                                              Rectangle(extent={{
                    20,60},{100,-60}},                                                                                                    lineColor = {0, 0, 0},
                  fillPattern =                                                                                                    FillPattern.HorizontalCylinder, fillColor = {192, 192, 192}), Rectangle(extent={{
                    38,40},{100,-40}},                                                                                                    lineColor = {0, 0, 0},
                  fillPattern =                                                                                                    FillPattern.HorizontalCylinder, fillColor = {0, 127, 255}), Ellipse(extent={{
                    -100,80},{60,-80}},                                                                                                    fillColor = {255, 255, 255},
                  fillPattern =                                                                                                    FillPattern.Solid, lineColor = {0, 0, 255}), Polygon(points={{
                    -60,70},{60,0},{-60,-68},{-60,70}},                                                                                                    lineColor = {0, 0, 255}, fillColor = {0, 0, 255},
                  fillPattern =                                                                                                    FillPattern.Solid), Text(extent={{
                    -54,32},{16,-30}},                                                                                                    lineColor = {255, 0, 0}, textString = "m"), Ellipse(extent={{
                    -26,30},{-18,22}},                                                                                                    lineColor = {255, 0, 0}, fillColor = {255, 0, 0},
                  fillPattern =                                                                                                    FillPattern.Solid)}));
      end Sink;
    end Sources;

    package HEX "Heat exchangers"
      model Evaporator
        "Evaporator with three zones preheater, evaporator, superheater"
        extends VIP.CycleTempo.Components.HEX.Heat_exchanger;
        parameter Modelica.SIunits.TemperatureDifference dT_int
          "Temperature difference at the inlet of the evaporator";
        parameter Boolean use_dT_int = false
          "true if temperature difference at the inlet of the evaporator is given";
        Medium_c.SaturationProperties sat "Saturated state hot fluid";
        Modelica.SIunits.HeatFlowRate qdot_pre "Heat rate pre-heater";
        Modelica.SIunits.HeatFlowRate qdot_eva "Heat rate evaporator";
        Modelica.SIunits.HeatFlowRate qdot_sup "Heat rate superheater";
        Nodes.Node_out node_c_in_eva(redeclare package Medium = Medium_c)
          "Inlet node evaporating section cold side"
          annotation (Placement(transformation(extent={{-66,-38},{-54,-26}}),
              iconTransformation(extent={{-66,-38},{-54,-26}})));
        Nodes.Node_out node_h_in_eva(redeclare package Medium = Medium_h)
          "Inlet node evaporating section hot side" annotation (Placement(
              transformation(extent={{20,-6},{32,6}}), iconTransformation(extent={{20,-6},
                  {32,6}})));
        Nodes.Node_out node_c_in_sup(redeclare package Medium = Medium_c)
          "Inlet node superheater section cold side" annotation (Placement(
              transformation(extent={{34,54},{46,66}}),     iconTransformation(extent=
                 {{54,-38},{66,-26}})));
        Nodes.Node_out node_h_in_pre(redeclare package Medium = Medium_h)
          "Inlet node superheater section hot side" annotation (Placement(
              transformation(extent={{-6,34},{6,46}}), iconTransformation(extent={{-34,
                  -6},{-22,6}})));
      equation

        //Boundary equations
        if use_dT_int then
          dT_int         = Medium_h.temperature(Medium_h.setState_ph(node_h_in_pre.p,
            node_h_in_pre.h)) - sat.Tsat;
        end if;

        //Component equations
        sat              = Medium_c.setSat_p(node_c_in.p);
        node_h_in.m_flow = node_h_in_pre.m_flow;
        node_h_in.m_flow = node_h_in_eva.m_flow;
        node_h_in.p      = node_h_in_pre.p;
        node_h_in.p      = node_h_in_eva.p;
        node_c_in.m_flow = node_c_in_sup.m_flow;
        node_c_in.m_flow = node_c_in_eva.m_flow;
        node_c_in.p      = node_c_in_sup.p;
        node_c_in.p      = node_c_in_eva.p;
        sat.hv           = node_c_in_sup.h;
        sat.hl           = node_c_in_eva.h;
        qdot_pre         = node_c_in.m_flow*(sat.hl - node_c_in.h);
        qdot_eva         = node_c_in.m_flow*(sat.hv - sat.hl);
        qdot_sup         = node_c_in.m_flow*(node_c_out.h - sat.hv);
        qdot_pre         = node_h_in.m_flow*(node_h_in_pre.h - node_h_out.h);
        qdot_eva         = node_h_in.m_flow*(node_h_in_eva.h - node_h_in_pre.h);
        annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
                  -100},{100,100}}), graphics), Icon(coordinateSystem(
                preserveAspectRatio=false, extent={{-100,-100},{100,100}}), graphics));
      end Evaporator;

      model Condenser "Condenser model"
        extends VIP.CycleTempo.Components.HEX.Heat_exchanger;
        Medium_h.SaturationProperties sat "Saturated state hot fluid";
      equation
        sat          = Medium_h.setSat_p(node_h_out.p);
        node_h_out.h = sat.hl;

      end Condenser;

      model Heat_exchanger "Heat exchanger for single phase heat transfer"
        //Medium model
        replaceable package Medium_h = VIP.Media.OneRandomOrganicFluid  constrainedby
          Modelica.Media.Interfaces.PartialMedium "Medium model hot cells" annotation(choicesAllMatching = true);
        replaceable package Medium_c = VIP.Media.OneRandomOrganicFluid  constrainedby
          Modelica.Media.Interfaces.PartialMedium "Medium model hot cells" annotation(choicesAllMatching = true);
        //Input parameters
        parameter Modelica.SIunits.Temperature T_h_out
          "Outlet temperature hot side"                                              annotation (Dialog(tab="Addco"));
        parameter Modelica.SIunits.AbsolutePressure dp_h
          "Pressure drop hot side";
        parameter Modelica.SIunits.Temperature T_c_out
          "Outlet temperature cold side"                                              annotation (Dialog(tab="Addco"));
        parameter Modelica.SIunits.AbsolutePressure dp_c
          "Pressure drop cold side";
        parameter Modelica.SIunits.TemperatureDifference dT_in
          "Inlet temperature difference" annotation (Dialog(tab="Addco"));
        parameter Modelica.SIunits.TemperatureDifference dT_out
          "Outlet temperature difference" annotation (Dialog(tab="Addco"));
        parameter Boolean use_T_h_out = false
          "true if outlet temperature hot side is given" annotation (Dialog(tab="Addco"));
        parameter Boolean use_T_c_out = false
          "true if outlet temperature cold side is given" annotation (Dialog(tab="Addco"));
        parameter Boolean use_dT_in = false
          "true if inlet temperature difference is given" annotation (Dialog(tab="Addco"));
        parameter Boolean use_dT_out = false
          "true if outlet temperature difference is given" annotation (Dialog(tab="Addco"));
        parameter Modelica.SIunits.Temperature T_h_in_start = 323.15
          "Inlet temperature hot side start value" annotation (Dialog(tab="Start"));
        parameter Modelica.SIunits.Temperature T_h_out_start = 303.15
          "Outlet temperature hot side start value" annotation (Dialog(tab="Start"));
        parameter Modelica.SIunits.Temperature T_c_in_start = 298.15
          "Inlet temperature cold side start value"  annotation (Dialog(tab="Start"));
        parameter Modelica.SIunits.Temperature T_c_out_start = 313.15
          "Outlet temperature cold side start value"  annotation (Dialog(tab="Start"));
        parameter Modelica.SIunits.AbsolutePressure p_h_in_start = 1e5
          "Inlet pressure hot side start value" annotation (Dialog(tab="Start"));
        parameter Modelica.SIunits.AbsolutePressure p_h_out_start = p_h_in_start
          "Outlet pressure hot side start value" annotation (Dialog(tab="Start"));
        parameter Modelica.SIunits.AbsolutePressure p_c_in_start = 1e5
          "Inlet pressure cold side start value"  annotation (Dialog(tab="Start"));
        parameter Modelica.SIunits.AbsolutePressure p_c_out_start = p_c_in_start
          "Outlet pressure cold side start value"  annotation (Dialog(tab="Start"));
        parameter Medium_h.SpecificEnthalpy h_h_in_start = Medium_h.specificEnthalpy_pT(
         p_h_in_start, T_h_in_start) "Inlet enthalpy start value hot side" annotation (Dialog(tab="Start"));
        parameter Medium_h.SpecificEnthalpy h_h_out_start = Medium_h.specificEnthalpy_pT(
         p_h_out_start, T_h_out_start) "Outlet enthalpy start value hot side" annotation (Dialog(tab="Start"));
        parameter Medium_c.SpecificEnthalpy h_c_in_start = Medium_c.specificEnthalpy_pT(
         p_c_in_start, T_c_in_start) "Inlet enthalpy start value cold side" annotation (Dialog(tab="Start"));
        parameter Medium_c.SpecificEnthalpy h_c_out_start = Medium_c.specificEnthalpy_pT(
         p_c_out_start, T_c_out_start) "Outlet enthalpy start value cold side" annotation (Dialog(tab="Start"));
        Modelica.SIunits.HeatFlowRate qdot "Heat rate";

        Nodes.Node_in node_c_in(redeclare package Medium = Medium_c,
        h(start = h_c_in_start), p(start = p_c_in_start))
          "Inlet node cold side"                                                  annotation (Placement(transformation(extent={{90,-70},{110,-50}}),
              iconTransformation(
              extent={{-10,-10},{10,10}},
              rotation=-90,
              origin={-60,-100})));
        Nodes.Node_out node_c_out(redeclare package Medium = Medium_c,
        h(start = h_c_out_start), p(start = p_c_out_start))
          "Outlet node cold side"                                                   annotation (Placement(transformation(extent={{90,50},{110,70}}),
              iconTransformation(
              extent={{-10,-10},{10,10}},
              rotation=-90,
              origin={60,-100})));
        Nodes.Node_in node_h_in(redeclare package Medium = Medium_h,
        h(start = h_h_in_start), p(start = p_h_in_start)) "Inlet node hot side"
                                                                                annotation (Placement(transformation(extent={{-10,-110},{10,-90}}),
              iconTransformation(extent={{-10,-10},{10,10}},
              rotation=-90,
              origin={100,0})));
        Nodes.Node_out node_h_out(redeclare package Medium = Medium_h,
        h(start = h_h_out_start), p(start = p_h_out_start))
          "Outlet node hot side"                                                   annotation (Placement(transformation(extent={{-10,88},{10,108}}),
              iconTransformation(extent={{-10,-10},{10,10}},
              rotation=-90,
              origin={-100,0})));
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
        if use_T_h_out then
          T_h_out  = Medium_h.temperature(Medium_h.setState_ph(node_h_out.p,
          node_h_out.h));
        end if;

        if use_T_c_out then
          T_c_out  = Medium_c.temperature(Medium_c.setState_ph(node_c_out.p,
          node_c_out.h));
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
          Medium_c.temperature_ph(node_c_in.p, node_c_in.h);
        end if;

        connect(node_c_in, node_c_in) annotation (Line(
            points={{100,-60},{100,-60}},
            color={0,127,255},
            smooth=Smooth.None));
        annotation(Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},
                  {100,100}}),
                           graphics), Icon(coordinateSystem(extent={{-100,-100},{100,100}},      preserveAspectRatio=false,  initialScale = 0.1, grid = {2, 2}), graphics={  Rectangle(lineColor = {0, 0, 255}, fillColor = {230, 230, 230},
                  fillPattern =                                                                                                    FillPattern.Solid, extent={{
                    -100,100},{100,-100}},
                origin={1.42109e-014,-1.42109e-014},
                rotation=-90),                                                                                                    Line(points={{
                    100,-61},{-58,-61},{-11,0},{-58,59},{102,59}},                                                                                                   color = {0, 0, 255}, thickness = 0.5,
                origin={1,2},
                rotation=-90),                                                                                                    Text(lineColor = {85, 170, 255}, extent={{
                    -100,15},{100,-15}},                                                                                                    textString = "%name",
                origin={0,112},
                rotation=0)}));
      end Heat_exchanger;
    end HEX;

    package Sensors "Sensors for the system"
      model State "Sensor of thermodynamic state"
        replaceable package Medium = VIP.Media.OneRandomOrganicFluid  constrainedby
          Modelica.Media.Interfaces.PartialMedium "Medium model hot cells" annotation(choicesAllMatching = true);
        Medium.ThermodynamicState state "Thermodynamic state";
        Nodes.Node_in node_in(redeclare package Medium = Medium)
          annotation (Placement(transformation(extent={{-70,-10},{-50,10}}),
              iconTransformation(extent={{-70,-10},{-50,10}})));
        Nodes.Node_out node_out(redeclare package Medium = Medium)
          annotation (Placement(transformation(extent={{50,-10},{70,10}}),
              iconTransformation(extent={{50,-10},{70,10}})));
      equation
        //Component equations
        state = Medium.setState_ph(node_in.p, node_in.h);
      //   node_in.p      = node_out.p;
      //   node_in.h      = node_out.h;
      //   node_in.m_flow = node_out.m_flow;
        connect(node_in, node_out) annotation (Line(
            points={{-60,0},{60,0}},
            color={0,127,255},
            smooth=Smooth.None));
        connect(node_in, node_in) annotation (Line(
            points={{-60,0},{-40,0},{-40,0},{-60,0}},
            color={0,127,255},
            smooth=Smooth.None));
        annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
                  -100},{100,100}}), graphics), Icon(coordinateSystem(
                preserveAspectRatio=false, extent={{-100,-100},{100,100}}), graphics={
              Ellipse(extent={{-50,20},{50,-20}}, lineColor={0,0,0},
                fillColor={85,170,255},
                fillPattern=FillPattern.Solid), Text(
                extent={{-26,10},{30,-8}},
                lineColor={0,0,0},
                textString="%name")}));
      end State;
    end Sensors;
  end Components;
end CycleTempo;
