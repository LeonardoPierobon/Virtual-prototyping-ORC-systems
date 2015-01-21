within ;
package Test "Here we make experiments on the VIP "
extends Modelica.Icons.ExamplesPackage;

  package Media "Fluid library"
    package OneRandomOrganicFluid
      "Change the name to something more appropriate"
        extends ExternalMedia.Media.FluidPropMedium(
        mediumName = "Name of the fluid for documentation purposes",
        libraryName = "FluidProp.RefProp",
        substanceNames = {"cyclopentane"},
        ThermoStates = Modelica.Media.Interfaces.PartialMedium.Choices.IndependentVariables.ph);
    end OneRandomOrganicFluid;

    package CoolProp "Media defined using CoolProp"
      package MM "Hexamethyldisiloxane (Coolprop)."
          extends ExternalMedia.Media.CoolPropMedium(
          mediumName = "MM",
          substanceNames = {"MM"},
          ThermoStates = Modelica.Media.Interfaces.PartialMedium.Choices.IndependentVariables.ph);
      end MM;

      package Pentane "CoolProp model of Pentane"
        extends ExternalMedia.Media.CoolPropMedium(
          mediumName = "Pentane",
          substanceNames = {"Pentane"},
          ThermoStates = Modelica.Media.Interfaces.PartialMedium.Choices.IndependentVariables.ph);
      end Pentane;

      package Cyclopentane "CoolProp model of Cyclopentane"
        extends ExternalMedia.Media.CoolPropMedium(
          mediumName = "Cyclopentane",
          substanceNames = {"Cyclopentane"},
          ThermoStates = Modelica.Media.Interfaces.PartialMedium.Choices.IndependentVariables.ph);
      end Cyclopentane;

      package Water "CoolProp model of Water"
        extends ExternalMedia.Media.CoolPropMedium(
          mediumName = "Water",
          substanceNames = {"Water"},
          ThermoStates = Modelica.Media.Interfaces.PartialMedium.Choices.IndependentVariables.ph);
      end Water;

      package Methanol "CoolProp model of Methanol"
        extends ExternalMedia.Media.CoolPropMedium(
          mediumName = "Methanol",
          substanceNames = {"Methanol"},
          ThermoStates = Modelica.Media.Interfaces.PartialMedium.Choices.IndependentVariables.ph);
      end Methanol;

      package Water_TTSE "CoolProp model of Water. TTSE is on. "
        extends ExternalMedia.Media.CoolPropMedium(
          mediumName = "Water",
          substanceNames = {"Water|enable_TTSE=1|enable_EXTTP=1"},
          ThermoStates = Modelica.Media.Interfaces.PartialMedium.Choices.IndependentVariables.ph);
      end Water_TTSE;

      package Mixture_propane_butane
        "Refprop model of a mixture of propane and butane"
        extends ExternalMedia.Media.CoolPropMedium(
          mediumName = "Propane-Butane",
          substanceNames = {"REFPROP-MIX:Butane[0.5]&Propane[0.5]"});
      end Mixture_propane_butane;

      package Methanol_TTSE "CoolProp model of Methanol. TTSE is on."
        extends ExternalMedia.Media.CoolPropMedium(
          mediumName = "Methanol",
          substanceNames = {"Methanol|enable_TTSE=1|enable_EXTTP=1"},
          ThermoStates = Modelica.Media.Interfaces.PartialMedium.Choices.IndependentVariables.ph);
      end Methanol_TTSE;

      package THERM66 "Therminol 66 properties from CoolProp"
        extends ExternalMedia.Media.IncompressibleCoolPropMedium(
          mediumName="T66",
          substanceNames={"T66|calc_transport=1"},
          ThermoStates=Modelica.Media.Interfaces.Choices.IndependentVariables.pT);
      end THERM66;

      package MM_TTSE "Hexamethyldisiloxane (Coolprop). Tables are on."
          extends ExternalMedia.Media.CoolPropMedium(
          mediumName = "MM",
          substanceNames = {"MM|enable_TTSE=1|enable_EXTTP=1"},
          ThermoStates = Modelica.Media.Interfaces.PartialMedium.Choices.IndependentVariables.ph);
      end MM_TTSE;

      package Air "CoolProp model of Air"
        extends ExternalMedia.Media.CoolPropMedium(
          mediumName = "Air",
          substanceNames = {"Air"},
          ThermoStates = Modelica.Media.Interfaces.PartialMedium.Choices.IndependentVariables.ph);
      end Air;

      package Air_TTSE "Air (Coolprop). Tables are on."
          extends ExternalMedia.Media.CoolPropMedium(
          mediumName = "Air",
          substanceNames = {"Air|enable_TTSE=1|enable_EXTTP=1"},
          ThermoStates = Modelica.Media.Interfaces.PartialMedium.Choices.IndependentVariables.ph);
      end Air_TTSE;

      package Toluene "Toluene (Coolprop)."
          extends ExternalMedia.Media.CoolPropMedium(
          mediumName = "Toluene",
          substanceNames = {"Toluene"},
          ThermoStates = Modelica.Media.Interfaces.PartialMedium.Choices.IndependentVariables.ph);
      end Toluene;

      package Toluene_TTSE "Toluene (Coolprop). TTSE is on."
          extends ExternalMedia.Media.CoolPropMedium(
          mediumName = "Toluene",
          substanceNames = {"Toluene|enable_TTSE=1|enable_EXTTP=1"},
          ThermoStates = Modelica.Media.Interfaces.PartialMedium.Choices.IndependentVariables.ph);
      end Toluene_TTSE;
    end CoolProp;

    package RefProp "Media defined using RefProp (via CoolProp)"
      package MM "Hexamethyldisiloxane (Refprop). "
          extends ExternalMedia.Media.CoolPropMedium(
          mediumName = "MM",
          substanceNames = {"REFPROP-MM"},
          ThermoStates = Modelica.Media.Interfaces.PartialMedium.Choices.IndependentVariables.ph);
      end MM;

      package Pentane "CoolProp model of Pentane"
        extends ExternalMedia.Media.CoolPropMedium(
          mediumName = "Pentane",
          substanceNames = {"REFPROP-Pentane"},
          ThermoStates = Modelica.Media.Interfaces.PartialMedium.Choices.IndependentVariables.ph);
      end Pentane;

      package Cyclopentane "CoolProp model of Cyclopentane"
        extends ExternalMedia.Media.CoolPropMedium(
          mediumName = "cyclopentane",
          substanceNames = {"REFPROP-cyclopentane"},
          ThermoStates = Modelica.Media.Interfaces.PartialMedium.Choices.IndependentVariables.ph);
      end Cyclopentane;

      package Water "water (Refprop)."
        extends ExternalMedia.Media.CoolPropMedium(
          mediumName = "water",
          substanceNames = {"REFPROP-water"},
          ThermoStates = Modelica.Media.Interfaces.PartialMedium.Choices.IndependentVariables.ph);
      end Water;

      package Methanol "methanol (Refprop). "
          extends ExternalMedia.Media.CoolPropMedium(
          mediumName = "methanol",
          substanceNames = {"REFPROP-methanol"},
          ThermoStates = Modelica.Media.Interfaces.PartialMedium.Choices.IndependentVariables.ph);
      end Methanol;
    end RefProp;

    package FluidProp "Media defined using FluidProp"
      package MM "Hexamethyldisiloxane "
          extends ExternalMedia.Media.FluidPropMedium(
          mediumName = "MM",
          substanceNames = {"MM"},
          libraryName = "FluidProp.Refprop");
      end MM;

      package Pentane "CoolProp model of Pentane"
        extends ExternalMedia.Media.FluidPropMedium(
          mediumName = "Pentane",
          libraryName = "FluidProp.StanMix3",
          substanceNames = {"Pentane"},
          ThermoStates = Modelica.Media.Interfaces.PartialMedium.Choices.IndependentVariables.ph);
      end Pentane;

      package Cyclopentane "CoolProp model of Cyclopentane"
        extends ExternalMedia.Media.FluidPropMedium(
          mediumName = "Cyclopentane",
          libraryName = "FluidProp.Refprop",
          substanceNames = {"Cyclopentane"},
          ThermoStates = Modelica.Media.Interfaces.PartialMedium.Choices.IndependentVariables.ph);
      end Cyclopentane;

      package Water "FluidProp model of Water"
        extends ExternalMedia.Media.FluidPropMedium(
          mediumName = "Water",
          libraryName = "FluidProp.Refprop",
          substanceNames = {"Water"},
          ThermoStates = Modelica.Media.Interfaces.PartialMedium.Choices.IndependentVariables.ph);
      end Water;

      package Methanol "CoolProp model of Methanol"
        extends ExternalMedia.Media.FluidPropMedium(
          mediumName = "methanol",
          libraryName = "FluidProp.StanMix3",
          substanceNames = {"methanol"},
          ThermoStates = Modelica.Media.Interfaces.PartialMedium.Choices.IndependentVariables.ph);
      end Methanol;

      package Mixture_pentane_ethane
        "Refprop model of a mixture of pentane and ethane"
        extends ExternalMedia.Media.FluidPropMedium(
          mediumName = "Pentane-Hexane",
          libraryName = "FluidProp.StanMix3",
          substanceNames = {"Pentane-Hexane"},
          ThermoStates = Modelica.Media.Interfaces.PartialMedium.Choices.IndependentVariables.ph);
      //     libraryName = "FluidProp.StanMix",
      end Mixture_pentane_ethane;

      package Mixture_propane_butane
        "Refprop model of a mixture of propane and butane"
        extends ExternalMedia.Media.FluidPropMedium(
          mediumName = "Propane-Butane",
          libraryName = "FluidProp.StanMix3",
          substanceNames = {"Propane-Butane"},
          ThermoStates = Modelica.Media.Interfaces.PartialMedium.Choices.IndependentVariables.ph);
      //     libraryName = "FluidProp.StanMix",
      end Mixture_propane_butane;
    end FluidProp;

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

  package Cycle "Tests for library CycleTempo 2.0"
    model ORCHID1 "Design of the ORCHID test rig. Therminol as heat source."
      parameter Medium.SaturationProperties  sat =  Medium.setSat_p(0.33e5)
        "Saturation properties";
      package Hot = Media.CoolProp.THERM66;
      package Medium = Media.CoolProp.MM_TTSE;
      package Coolant = Media.CoolProp.Water_TTSE;

      CycleTempo.Components.Turbomachinery.Pump Pump(
        redeclare package Medium = Medium,
        eta_is=0.9,
        eta_m=0.9) "Centrifugal pump"
        annotation (Placement(transformation(extent={{13,14},{-13,-14}},
            rotation=-90,
            origin={77,-48})));
      CycleTempo.Components.HEX.Evaporator evaporator(
        dp_h=0,
        redeclare package Medium_c = Medium,
        redeclare package Medium_h = Hot,
        dp_c=0,
        hh_in_start=Hot.specificEnthalpy(Hot.setState_pT(2e5, 500)))
                        annotation (Placement(transformation(
            extent={{-18,-15.5},{18,15.5}},
            rotation=0,
            origin={0,59.5})));
      CycleTempo.Components.Turbomachinery.Turbine turbine(
        redeclare package Medium = Medium,
        eta_m=0.9,
        eta_is=0.8) "Axial turbine"
        annotation (Placement(transformation(extent={{53,17},{92,50}})));
      CycleTempo.Components.HEX.Condenser condenser(
        dp_h=0,
        redeclare package Medium_h = Medium,
        redeclare package Medium_c = Coolant,
        dp_c=0) "condenser model"           annotation (Placement(transformation(
            extent={{-12,-11},{12,11}},
            rotation=0,
            origin={0,-72})));
      CycleTempo.Components.Flags.CC CC(redeclare package Medium = Medium)
        annotation (Placement(transformation(
            extent={{10,-10},{-10,10}},
            rotation=-90,
            origin={77,-78})));
      CycleTempo.Components.HEX.Heat_exchanger recuperator(
        redeclare package Medium_h = Medium,
        redeclare package Medium_c = Medium,
        dp_h=0,
        dp_c=0,
        dT_out=20,
        use_dT_out=true)
        annotation (Placement(transformation(extent={{11,-32},{-11,-12}})));
      CycleTempo.Components.Flags.ADDCO N3(redeclare package Medium = Medium)
        annotation (Placement(transformation(extent={{-98,8},{-78,30}})));
      CycleTempo.Components.Flags.ADDCO N20(
        redeclare package Medium = Coolant,
        use_T=true,
        use_p=true,
        T=278.15,
        p=100000)
        annotation (Placement(transformation(extent={{-97,-82},{-79,-63}})));
      CycleTempo.Components.Flags.ADDCO N21(
        redeclare package Medium = Coolant,
        use_T=true,
        T=313.15) annotation (Placement(transformation(extent={{45,-81},{23,-63}})));
      CycleTempo.Components.Flags.ADDCO N2(
        redeclare package Medium = Medium,
        use_p=true,
        p=1460000) annotation (Placement(transformation(
            extent={{10,-10},{-10,10}},
            rotation=0,
            origin={87,-22})));
      CycleTempo.Components.Flags.ADDCO N6(
        redeclare package Medium = Medium,
        use_T=false,
        T=343.15) annotation (Placement(transformation(
            extent={{9,9},{-9,-9}},
            rotation=180,
            origin={-88,-41})));
      CycleTempo.Components.Flags.ADDCO N10(
        redeclare package Medium = Hot,
        m_flow=0.5,
        use_T=true,
        use_p=true,
        use_m_flow=true,
        T=618.15,
        p=200000)         annotation (Placement(transformation(
            extent={{-11,-10},{11,10}},
            rotation=0,
            origin={-11,96})));
      CycleTempo.Components.Flags.ADDCO N11(
        redeclare package Medium = Hot,
        visible=true,
        use_T=true,
        T=513.15) annotation (Placement(transformation(
            extent={{-10,-9},{10,9}},
            rotation=0,
            origin={-10,20})));
      CycleTempo.Components.Flags.ADDCO N4(
        redeclare package Medium = Medium,
        use_T=true,
        T=593.15) annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=0,
            origin={44,96})));
      CycleTempo.Components.Flags.ADDCO N5(redeclare package Medium = Medium,
          visible=true) annotation (Placement(transformation(
            extent={{10.5,-9.5},{-10.5,9.5}},
            rotation=0,
            origin={134.5,3.5})));
      CycleTempo.Components.Flags.ADDCO N1(
        redeclare package Medium = Medium,
        use_T=false,
        use_p=true,
        p=33000) annotation (Placement(transformation(
            extent={{-9,-9},{9,9}},
            rotation=0,
            origin={-88,-105})));
      CycleTempo.Components.Flags.START START1(
        redeclare package Medium = Medium,
        m_flow=0.15,
        p=1460000,
        T=439.15)
        annotation (Placement(transformation(extent={{-96,32},{-78,50}})));
      CycleTempo.Components.Electrics.Generator generator(eta_el=0.98)
        annotation (Placement(transformation(extent={{102,18},{133,49}})));
      CycleTempo.Components.Flags.ADDCOW Pout(W=33e3, use_W=false)
        annotation (Placement(transformation(extent={{162,23},{142,44}})));
      CycleTempo.Components.Flags.ADDCOW Ppump(W=33e3, use_W=false)
        annotation (Placement(transformation(extent={{152,-58},{132,-37}})));
      CycleTempo.Components.Electrics.Motor motor(eta_el=0.9)
        annotation (Placement(transformation(extent={{97,-60},{122,-36}})));
      CycleTempo.Components.Flags.START START2(   redeclare package Medium = Medium,
        m_flow=0.15,
        p=33000,
        T=343.15)
        annotation (Placement(transformation(extent={{26,-50},{8,-32}})));
    equation
      connect(CC.node_in, N1.node) annotation (Line(
          points={{77,-83.8},{77,-105},{-78.91,-105}},
          color={0,0,0},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(N11.node, evaporator.node_h_out) annotation (Line(
          points={{0.1,20},{0.1,30},{0,30},{0,37.8}},
          color={0,0,0},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(N10.node, evaporator.node_h_in) annotation (Line(
          points={{0.11,96},{0.11,86},{0,86},{0,81.51}},
          color={0,0,0},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(N5.node, recuperator.node_h_in) annotation (Line(
          points={{123.895,3.5},{0,3.5},{0,-7.8}},
          color={0,0,0},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(condenser.node_h_out, N1.node) annotation (Line(
          points={{0,-87.4},{0,-105},{-78.91,-105}},
          color={0,0,0},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(CC.node_out, Pump.node_in) annotation (Line(
          points={{77,-72},{77,-60.74}},
          color={0,0,0},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(Pump.node_out, N2.node) annotation (Line(
          points={{77,-35},{77,-22},{76.9,-22}},
          color={0,0,0},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(recuperator.node_c_in, N2.node) annotation (Line(
          points={{15.4,-22},{76.9,-22}},
          color={0,0,0},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(recuperator.node_c_out, N3.node) annotation (Line(
          points={{-15.51,-22},{-77.9,-22},{-77.9,19}},
          color={0,0,0},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(N20.node, condenser.node_c_in) annotation (Line(
          points={{-78.91,-72.5},{-26,-72.5},{-26,-72},{-16.8,-72}},
          color={0,0,0},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(condenser.node_c_out, N21.node) annotation (Line(
          points={{16.92,-72},{22.89,-72}},
          color={0,0,0},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(N6.node, condenser.node_h_in) annotation (Line(
          points={{-78.91,-41},{0,-41},{0,-56.38}},
          color={0,0,0},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(N3.node, START1.node) annotation (Line(
          points={{-77.9,19},{-77.9,30.5},{-78,30.5},{-78,41}},
          color={0,0,0},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(START1.node, evaporator.node_c_in) annotation (Line(
          points={{-78,41},{-78,59.5},{-25.2,59.5}},
          color={0,0,0},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(evaporator.node_c_out, turbine.node_in) annotation (Line(
          points={{25.38,59.5},{64.7,59.5},{64.7,50}},
          color={0,0,0},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(turbine.node_out, N5.node) annotation (Line(
          points={{75.23,13.7},{75.23,10.75},{123.895,10.75},{123.895,3.5}},
          color={0,0,0},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(turbine.terminal, generator.terminal_in) annotation (Line(
          points={{91.9025,33.5},{102.31,33.5}},
          color={0,0,0},
          smooth=Smooth.None));
      connect(generator.terminal_out, Pout.terminal) annotation (Line(
          points={{133,33.5},{142,33.5}},
          color={0,0,0},
          smooth=Smooth.None));
      connect(Ppump.terminal, motor.terminal_out) annotation (Line(
          points={{132,-47.5},{122,-47.5},{122,-48}},
          color={0,0,0},
          smooth=Smooth.None));
      connect(Pump.terminal, motor.terminal_in) annotation (Line(
          points={{91,-48.065},{95.5,-48.065},{95.5,-48},{97.25,-48}},
          color={0,0,0},
          smooth=Smooth.None));
      connect(N4.node, turbine.node_in) annotation (Line(
          points={{54.1,96},{65,96},{65,60},{64.7,60},{64.7,50}},
          color={0,0,0},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(recuperator.node_h_out, condenser.node_h_in) annotation (Line(
          points={{0,-36},{0,-56.38}},
          color={0,0,0},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(START2.node, condenser.node_h_in) annotation (Line(
          points={{8,-41},{0,-41},{0,-56.38}},
          color={0,0,0},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-120,
                -120},{180,120}}),      graphics), Icon(coordinateSystem(extent={{-120,
                -120},{180,120}})),
        experiment(__Dymola_NumberOfIntervals=1, __Dymola_Algorithm="Dassl"),
        __Dymola_experimentSetupOutput);
    end ORCHID1;

    model ORCHID2 "Design of the ORCHID test rig. Air as heat source."
      parameter Medium.SaturationProperties  sat =  Medium.setSat_p(0.33e5)
        "Saturation properties";
      package Hot = Media.CoolProp.Air "Medium model hot cells";
      package Medium = Media.RefProp.MM "Medium model hot cells";
      package Coolant = Media.CoolProp.Water "Medium model hot cells";

      CycleTempo.Components.Turbomachinery.Pump Pump(
        redeclare package Medium = Medium,
        eta_is=0.9,
        eta_m=1) "Centrifugal pump"
        annotation (Placement(transformation(extent={{13,14},{-13,-14}},
            rotation=-90,
            origin={77,-48})));
      CycleTempo.Components.HEX.Evaporator evaporator(
        dp_h=0,
        redeclare package Medium_c = Medium,
        redeclare package Medium_h = Hot,
        dp_c=0,
        use_dT_int=true,
        dT_int=30,
        hh_in_start=7e5)
                        annotation (Placement(transformation(
            extent={{-18,-15.5},{18,15.5}},
            rotation=0,
            origin={0,59.5})));
      CycleTempo.Components.Turbomachinery.Turbine turbine(
        redeclare package Medium = Medium,
        eta_is=0.8,
        eta_m=1) "Axial turbine"
        annotation (Placement(transformation(extent={{53,17},{92,50}})));
      CycleTempo.Components.HEX.Condenser condenser(
        dp_h=0,
        redeclare package Medium_h = Medium,
        redeclare package Medium_c = Coolant,
        dp_c=0,
        dT_int=20,
        use_dT_int=true) "condenser model"  annotation (Placement(transformation(
            extent={{-12,-11},{12,11}},
            rotation=0,
            origin={0,-72})));
      CycleTempo.Components.Flags.CC CC(redeclare package Medium = Medium)
        annotation (Placement(transformation(
            extent={{10,-10},{-10,10}},
            rotation=-90,
            origin={77,-78})));
      CycleTempo.Components.HEX.Heat_exchanger recuperator(
        redeclare package Medium_h = Medium,
        redeclare package Medium_c = Medium,
        dp_h=0,
        dp_c=0,
        dT_out=20,
        use_dT_out=true)
        annotation (Placement(transformation(extent={{11,-32},{-11,-12}})));
      CycleTempo.Components.Flags.ADDCO N3(redeclare package Medium = Medium)
        annotation (Placement(transformation(extent={{-98,8},{-78,30}})));
      CycleTempo.Components.Flags.ADDCO N20(
        redeclare package Medium = Coolant,
        use_T=true,
        use_p=true,
        T=278.15,
        p=100000)
        annotation (Placement(transformation(extent={{-97,-82},{-79,-63}})));
      CycleTempo.Components.Flags.ADDCO N21(
        redeclare package Medium = Coolant,
        T=313.15,
        use_T=false)
                  annotation (Placement(transformation(extent={{45,-81},{23,-63}})));
      CycleTempo.Components.Flags.ADDCO N2(
        redeclare package Medium = Medium,
        use_p=true,
        p=1460000) annotation (Placement(transformation(
            extent={{10,-10},{-10,10}},
            rotation=0,
            origin={87,-22})));
      CycleTempo.Components.Flags.ADDCO N6(
        redeclare package Medium = Medium,
        use_T=false,
        T=343.15) annotation (Placement(transformation(
            extent={{9,9},{-9,-9}},
            rotation=180,
            origin={-88,-41})));
      CycleTempo.Components.Flags.ADDCO N10(
        redeclare package Medium = Hot,
        m_flow=0.5,
        use_T=true,
        use_p=true,
        use_m_flow=true,
        T=618.15,
        p=200000)         annotation (Placement(transformation(
            extent={{-11,-10},{11,10}},
            rotation=0,
            origin={-11,96})));
      CycleTempo.Components.Flags.ADDCO N11(
        redeclare package Medium = Hot,
        visible=true,
        use_T=true,
        T=513.15) annotation (Placement(transformation(
            extent={{-10,-9},{10,9}},
            rotation=0,
            origin={-10,20})));
      CycleTempo.Components.Flags.ADDCO N4(
        redeclare package Medium = Medium,
        use_T=false,
        T=593.15) annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=0,
            origin={44,96})));
      CycleTempo.Components.Flags.ADDCO N5(redeclare package Medium = Medium,
          visible=true) annotation (Placement(transformation(
            extent={{10.5,-9.5},{-10.5,9.5}},
            rotation=0,
            origin={134.5,3.5})));
      CycleTempo.Components.Flags.ADDCO N1(
        redeclare package Medium = Medium,
        use_T=false,
        use_p=true,
        p=33000) annotation (Placement(transformation(
            extent={{-9,-9},{9,9}},
            rotation=0,
            origin={-88,-105})));
      CycleTempo.Components.Flags.START START1(
        redeclare package Medium = Medium,
        m_flow=0.15,
        p=1460000,
        T=503.15)
        annotation (Placement(transformation(extent={{-96,32},{-78,50}})));
      CycleTempo.Components.Electrics.Generator generator(eta_el=0.98)
        annotation (Placement(transformation(extent={{102,18},{133,49}})));
      CycleTempo.Components.Flags.ADDCOW Pout(W=33e3, use_W=false)
        annotation (Placement(transformation(extent={{162,23},{142,44}})));
      CycleTempo.Components.Flags.ADDCOW Ppump(W=33e3, use_W=false)
        annotation (Placement(transformation(extent={{152,-58},{132,-37}})));
      CycleTempo.Components.Electrics.Motor motor(eta_el=0.9)
        annotation (Placement(transformation(extent={{97,-60},{122,-36}})));
      CycleTempo.Components.Flags.START START2(   redeclare package Medium = Medium,
        m_flow=0.15,
        p=33000,
        T=343.15)
        annotation (Placement(transformation(extent={{26,-50},{8,-32}})));
      CycleTempo.Components.Flags.START START3(
        redeclare package Medium = Medium,
        m_flow=0.15,
        p=1460000,
        T=593.15)
        annotation (Placement(transformation(extent={{88,101},{70,119}})));
    equation
      connect(CC.node_in, N1.node) annotation (Line(
          points={{77,-83.8},{77,-105},{-78.91,-105}},
          color={0,0,0},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(N11.node, evaporator.node_h_out) annotation (Line(
          points={{0.1,20},{0.1,30},{0,30},{0,37.8}},
          color={0,0,0},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(N10.node, evaporator.node_h_in) annotation (Line(
          points={{0.11,96},{0.11,86},{0,86},{0,81.51}},
          color={0,0,0},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(N5.node, recuperator.node_h_in) annotation (Line(
          points={{123.895,3.5},{0,3.5},{0,-7.8}},
          color={0,0,0},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(condenser.node_h_out, N1.node) annotation (Line(
          points={{0,-87.4},{0,-105},{-78.91,-105}},
          color={0,0,0},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(CC.node_out, Pump.node_in) annotation (Line(
          points={{77,-72},{77,-60.74}},
          color={0,0,0},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(Pump.node_out, N2.node) annotation (Line(
          points={{77,-35},{77,-22},{76.9,-22}},
          color={0,0,0},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(recuperator.node_c_in, N2.node) annotation (Line(
          points={{15.4,-22},{76.9,-22}},
          color={0,0,0},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(recuperator.node_c_out, N3.node) annotation (Line(
          points={{-15.51,-22},{-77.9,-22},{-77.9,19}},
          color={0,0,0},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(N20.node, condenser.node_c_in) annotation (Line(
          points={{-78.91,-72.5},{-26,-72.5},{-26,-72},{-16.8,-72}},
          color={0,0,0},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(condenser.node_c_out, N21.node) annotation (Line(
          points={{16.92,-72},{22.89,-72}},
          color={0,0,0},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(N6.node, condenser.node_h_in) annotation (Line(
          points={{-78.91,-41},{0,-41},{0,-56.38}},
          color={0,0,0},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(N3.node, START1.node) annotation (Line(
          points={{-77.9,19},{-77.9,30.5},{-78,30.5},{-78,41}},
          color={0,0,0},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(START1.node, evaporator.node_c_in) annotation (Line(
          points={{-78,41},{-78,59.5},{-25.2,59.5}},
          color={0,0,0},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(evaporator.node_c_out, turbine.node_in) annotation (Line(
          points={{25.38,59.5},{64.7,59.5},{64.7,50}},
          color={0,0,0},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(turbine.node_out, N5.node) annotation (Line(
          points={{75.23,13.7},{75.23,10.75},{123.895,10.75},{123.895,3.5}},
          color={0,0,0},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(turbine.terminal, generator.terminal_in) annotation (Line(
          points={{91.9025,33.5},{102.31,33.5}},
          color={0,0,0},
          smooth=Smooth.None));
      connect(generator.terminal_out, Pout.terminal) annotation (Line(
          points={{133,33.5},{142,33.5}},
          color={0,0,0},
          smooth=Smooth.None));
      connect(Ppump.terminal, motor.terminal_out) annotation (Line(
          points={{132,-47.5},{122,-47.5},{122,-48}},
          color={0,0,0},
          smooth=Smooth.None));
      connect(Pump.terminal, motor.terminal_in) annotation (Line(
          points={{91,-48.065},{95.5,-48.065},{95.5,-48},{97.25,-48}},
          color={0,0,0},
          smooth=Smooth.None));
      connect(N4.node, turbine.node_in) annotation (Line(
          points={{54.1,96},{65,96},{65,60},{64.7,60},{64.7,50}},
          color={0,0,0},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(recuperator.node_h_out, condenser.node_h_in) annotation (Line(
          points={{0,-36},{0,-56.38}},
          color={0,0,0},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(START2.node, condenser.node_h_in) annotation (Line(
          points={{8,-41},{0,-41},{0,-56.38}},
          color={0,0,0},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(START3.node, N4.node) annotation (Line(
          points={{70,110},{54.1,110},{54.1,96}},
          color={0,0,0},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-120,
                -120},{180,120}}),      graphics), Icon(coordinateSystem(extent={{-120,
                -120},{180,120}})));
    end ORCHID2;
  end Cycle;

  package Project "Tests for library Design 1.0"

    package Shell_and_tube "Tests for shell and tube heat exchangers"
      package Single_phase "Tests for one phase to one phase heat exchangers"
        model Coulson_Kern
          "Verification with the results given by Coulson et al. using the Kern method"
          package Medium_s = Media.RefProp.Methanol "Medium model";
          package Medium_t = Media.CoolProp.Water "Medium model";

          Design.Components.HEX.shell_and_tube ST(
            redeclare package Medium_s = Medium_s,
            redeclare package Medium_t = Medium_t,
            Dhyd=16e-3,
            l=4.83,
            pitch_f=1.25,
            layout=1,
            N_passes=2,
            N_baffles=27,
            redeclare Design.Materials.S_AISI_1040 Material_t,
            redeclare Design.Materials.S_AISI_1040 Material_s,
            redeclare function bundle_clearance =
                Design.Miscellanea.Shell_clearance.SRFH,
            thick_t=2e-3,
            thick_s=10e-3,
            N_tubes_start=600,
            t_s_in_start=368.15,
            t_s_out_start=313.15,
            t_t_in_start=298.15,
            t_t_out_start=313.15,
            p_s_start=320000,
            p_t_start=100000,
            m_s_start=27.8,
            m_t_start=68.9)
            annotation (Placement(transformation(extent={{-62,-55},{86,71}})));

          Design.Components.Flags.ADDCO HOT_OUT(
            use_T=true,
            redeclare package Medium = Medium_s,
            T=313.15) annotation (Placement(transformation(extent={{-109,-112},{-80,-84}})));
          Design.Components.Flags.ADDCO COLD_IN(
            m_flow=68.9,
            use_T=true,
            use_p=true,
            use_m_flow=true,
            redeclare package Medium = Medium_t,
            T=298.15,
            p=100000)
            annotation (Placement(transformation(extent={{-109,-51},{-80,-23}})));
          Design.Components.Flags.ADDCO COLD_OUT(
            redeclare package Medium = Medium_t,
            T=313.15,
            use_T=false)
            annotation (Placement(transformation(extent={{-109,40},{-80,68}})));
          Design.Components.Flags.ADDCO HOT_IN(
            use_T=true,
            use_p=true,
            use_m_flow=true,
            redeclare package Medium = Medium_s,
            m_flow=27.8,
            T=368.15,
            p=320000) annotation (Placement(transformation(extent={{11,69},{40,97}})));
        equation
          connect(COLD_IN.node, ST.tube_in) annotation (Line(
              points={{-79.855,-37},{-61.9275,-37},{-61.9275,-10.9},{-62,-10.9}},
              color={0,0,0},
              pattern=LinePattern.None,
              smooth=Smooth.None));
          connect(COLD_OUT.node, ST.tube_out) annotation (Line(
              points={{-79.855,54},{-62,54},{-62,27.215}},
              color={0,0,0},
              pattern=LinePattern.None,
              smooth=Smooth.None));
          connect(HOT_OUT.node, ST.shell_out) annotation (Line(
              points={{-79.855,-98},{-25,-98},{-25,-55}},
              color={0,0,0},
              pattern=LinePattern.None,
              smooth=Smooth.None));
          connect(HOT_IN.node, ST.shell_in) annotation (Line(
              points={{40.145,83},{63.8,83},{63.8,71}},
              color={0,0,0},
              pattern=LinePattern.None,
              smooth=Smooth.None));
         annotation (Placement(transformation(extent={{-108,-74},{88,66}})),
            experiment,
            __Dymola_experimentSetupOutput,
            Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-140,-120},{100,
                    100}}), graphics),
            Icon(coordinateSystem(extent={{-140,-120},{100,100}})));
        end Coulson_Kern;

        model Coulson_Bell_Delaware
          "Verification with the results given by Coulson using the Bell Delaware method."

          package Medium_s = Media.RefProp.Methanol "Medium model";
          package Medium_t = Media.CoolProp.Water_TTSE "Medium model";
          Design.Components.HEX.shell_and_tube ST(
            redeclare package Medium_s = Medium_s,
            redeclare package Medium_t = Medium_t,
            Dhyd=16e-3,
            l=4.83,
            pitch_f=1.25,
            layout=1,
            N_passes=2,
            N_baffles=27,
            redeclare Design.Materials.S_AISI_1040 Material_t,
            redeclare Design.Materials.S_AISI_1040 Material_s,
            redeclare function bundle_clearance =
                Design.Miscellanea.Shell_clearance.SRFH,
            redeclare Design.Heat_transfer.Shell.single_phase_Bell_Delaware
              hT_shell(
              pitch_f=ST.pitch_f,
              Dhyd_o=ST.bundle.Dhyd_o,
              N_tubes=ST.N_tubes,
              d_s=ST.d_s,
              d_b=ST.bundle.d_b,
              Dhyd=ST.d_s,
              layout=ST.layout,
              l_b=ST.l_b,
              N_ss=5,
              N_passes=ST.N_passes,
              ff=1.0,
              b_cut=0.25),
            redeclare Design.Pressure_drops.Shell.single_phase_Bell_Delaware
              dp_shell(
              pitch_f=ST.pitch_f,
              Dhyd_o=ST.bundle.Dhyd_o,
              N_tubes=ST.N_tubes,
              d_s=ST.d_s,
              d_b=ST.bundle.d_b,
              Dhyd=ST.d_s,
              layout=ST.layout,
              l_b=ST.l_b,
              N_ss=5,
              eta_wall=ST.shell.state[1].eta*ones(ST.Ncell),
              N_passes=ST.N_passes,
              N_baffles_d=ST.N_baffles_d,
              N_baffles=ST.N_baffles,
              ff=1.0,
              b_cut=0.25),
            thick_t=2e-3,
            thick_s=10e-3,
            N_tubes_start=600,
            m_s_start=27.8,
            m_t_start=68.9,
            t_s_in_start=368.15,
            t_s_out_start=313.15,
            t_t_in_start=298.15,
            t_t_out_start=313.15,
            p_s_start=320000,
            p_t_start=100000)
            annotation (Placement(transformation(extent={{-62,-55},{86,71}})));
          Design.Components.Flags.ADDCO HOT_OUT(
            use_T=true,
            redeclare package Medium = Test.Media.RefProp.Methanol,
            T=313.15) annotation (Placement(transformation(extent={{-109,-112},{-80,-84}})));
          Design.Components.Flags.ADDCO COLD_IN(
            m_flow=68.9,
            use_T=true,
            use_p=true,
            use_m_flow=true,
            redeclare package Medium = Test.Media.CoolProp.Water_TTSE,
            T=298.15,
            p=100000)
            annotation (Placement(transformation(extent={{-109,-51},{-80,-23}})));
          Design.Components.Flags.ADDCO COLD_OUT(
            redeclare package Medium = Test.Media.CoolProp.Water_TTSE,
            T=313.15,
            use_T=false)
            annotation (Placement(transformation(extent={{-109,40},{-80,68}})));
          Design.Components.Flags.ADDCO HOT_IN(
            use_T=true,
            use_p=true,
            use_m_flow=true,
            redeclare package Medium = Test.Media.RefProp.Methanol,
            m_flow=27.8,
            T=368.15,
            p=320000) annotation (Placement(transformation(extent={{11,69},{40,97}})));
        equation
          connect(COLD_IN.node,ST. tube_in) annotation (Line(
              points={{-79.855,-37},{-61.9275,-37},{-61.9275,-10.9},{-62,-10.9}},
              color={0,0,0},
              pattern=LinePattern.None,
              smooth=Smooth.None));
          connect(COLD_OUT.node,ST. tube_out) annotation (Line(
              points={{-79.855,54},{-62,54},{-62,27.215}},
              color={0,0,0},
              pattern=LinePattern.None,
              smooth=Smooth.None));
          connect(HOT_OUT.node,ST. shell_out) annotation (Line(
              points={{-79.855,-98},{-25,-98},{-25,-55}},
              color={0,0,0},
              pattern=LinePattern.None,
              smooth=Smooth.None));
          connect(HOT_IN.node,ST. shell_in) annotation (Line(
              points={{40.145,83},{63.8,83},{63.8,71}},
              color={0,0,0},
              pattern=LinePattern.None,
              smooth=Smooth.None));
         annotation (Placement(transformation(extent={{-108,-74},{88,66}})),
            experiment,
            __Dymola_experimentSetupOutput,
            Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-140,-120},{100,
                    100}}), graphics),
            Icon(coordinateSystem(extent={{-140,-120},{100,100}})));
        end Coulson_Bell_Delaware;

        model Aspen_Bell_Delaware
          "Verification with the results given by Aspen using the Bell Delaware method."

          package Medium_s = Media.RefProp.Methanol "Medium model";
          package Medium_t = Media.CoolProp.Water_TTSE "Medium model";

          Design.Components.HEX.shell_and_tube ST(
            redeclare package Medium_s = Medium_s,
            redeclare package Medium_t = Medium_t,
            Dhyd=16e-3,
            thick_t=2e-3,
            thick_s=10e-3,
            l=4.65,
            pitch_f=1.25,
            layout=1,
            N_passes=2,
            N_baffles=24,
            N_baffles_d=25,
            redeclare Design.Materials.S_AISI_1040 Material_t,
            redeclare Design.Materials.S_AISI_1040 Material_s,
            redeclare Design.Heat_transfer.Shell.single_phase_Bell_Delaware hT_shell(
              pitch_f=ST.pitch_f,
              Dhyd_o=ST.bundle.Dhyd_o,
              N_tubes=ST.N_tubes,
              d_s=ST.d_s,
              d_b=ST.bundle.d_b,
              Dhyd=ST.d_s,
              layout=ST.layout,
              l_b=ST.l_b,
              N_ss=5,
              N_passes=ST.N_passes,
              ff=1.0,
              b_cut=0.1588),
            redeclare Design.Pressure_drops.Shell.single_phase_Bell_Delaware dp_shell(
              pitch_f=ST.pitch_f,
              Dhyd_o=ST.bundle.Dhyd_o,
              N_tubes=ST.N_tubes,
              d_s=ST.d_s,
              d_b=ST.bundle.d_b,
              Dhyd=ST.d_s,
              layout=ST.layout,
              l_b=ST.l_b,
              N_ss=5,
              eta_wall=ST.shell.state[1].eta*ones(ST.Ncell),
              N_passes=ST.N_passes,
              N_baffles_d=ST.N_baffles_d,
              N_baffles=ST.N_baffles,
              ff=1.0,
              b_cut=0.1588),
            redeclare function bundle_clearance =
                Design.Miscellanea.Shell_clearance.Fixed_Utube,
            N_tubes_start=800,
            t_s_in_start=368.15,
            t_s_out_start=313.15,
            t_t_in_start=298.15,
            t_t_out_start=313.15,
            p_s_start=320000,
            p_t_start=100000,
            m_s_start=27.8,
            m_t_start=68.9)
            annotation (Placement(transformation(extent={{-42,-53},{106,73}})));
          Design.Components.Flags.ADDCO HOT_OUT(
            use_T=true,
            redeclare package Medium = Medium_s,
            T=313.15) annotation (Placement(transformation(extent={{-89,-110},{-60,-82}})));
          Design.Components.Flags.ADDCO COLD_IN(
            m_flow=68.9,
            use_T=true,
            use_p=true,
            use_m_flow=true,
            redeclare package Medium = Medium_t,
            T=298.15,
            p=100000)
            annotation (Placement(transformation(extent={{-89,-49},{-60,-21}})));
          Design.Components.Flags.ADDCO COLD_OUT(
            redeclare package Medium = Medium_t,
            T=313.15,
            use_T=false)
            annotation (Placement(transformation(extent={{-89,42},{-60,70}})));
          Design.Components.Flags.ADDCO HOT_IN(
            use_T=true,
            use_p=true,
            use_m_flow=true,
            redeclare package Medium = Medium_s,
            m_flow=27.8,
            T=368.15,
            p=3.2e5) annotation (Placement(transformation(extent={{31,71},{60,99}})));
        equation
          connect(COLD_IN.node, ST.tube_in) annotation (Line(
              points={{-59.855,-35},{-41.9275,-35},{-41.9275,-8.9},{-42,-8.9}},
              color={0,0,0},
              pattern=LinePattern.None,
              smooth=Smooth.None));
          connect(COLD_OUT.node, ST.tube_out) annotation (Line(
              points={{-59.855,56},{-42,56},{-42,29.215}},
              color={0,0,0},
              pattern=LinePattern.None,
              smooth=Smooth.None));
          connect(HOT_OUT.node, ST.shell_out) annotation (Line(
              points={{-59.855,-96},{-5,-96},{-5,-53}},
              color={0,0,0},
              pattern=LinePattern.None,
              smooth=Smooth.None));
          connect(HOT_IN.node, ST.shell_in) annotation (Line(
              points={{60.145,85},{83.8,85},{83.8,73}},
              color={0,0,0},
              pattern=LinePattern.None,
              smooth=Smooth.None));
         annotation (Placement(transformation(extent={{-108,-74},{88,66}})), Icon(
                coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{
                    100,100}}), graphics));
        end Aspen_Bell_Delaware;

        model Aspen_Bell_Delaware_inv
          "Verification with the results given by Aspen using the Bell Delaware method. Fluids are inverted."

          package Medium_t = Media.RefProp.Methanol "Medium model";
          package Medium_s = Media.CoolProp.Water_TTSE "Medium model";

          Design.Components.HEX.shell_and_tube ST(
            redeclare package Medium_s = Medium_s,
            redeclare package Medium_t = Medium_t,
            Dhyd=16e-3,
            thick_t=2e-3,
            thick_s=10e-3,
            l=4.65,
            pitch_f=1.25,
            layout=1,
            N_passes=4,
            N_baffles=24,
            redeclare Design.Materials.S_AISI_1040 Material_t,
            redeclare Design.Materials.S_AISI_1040 Material_s,
            redeclare Design.Heat_transfer.Shell.single_phase_Bell_Delaware hT_shell(
              pitch_f=ST.pitch_f,
              Dhyd_o=ST.bundle.Dhyd_o,
              N_tubes=ST.N_tubes,
              d_s=ST.d_s,
              d_b=ST.bundle.d_b,
              Dhyd=ST.d_s,
              layout=ST.layout,
              l_b=ST.l_b,
              N_ss=5,
              N_passes=ST.N_passes,
              ff=1.0,
              b_cut=0.1588),
            redeclare Design.Pressure_drops.Shell.single_phase_Bell_Delaware dp_shell(
              pitch_f=ST.pitch_f,
              Dhyd_o=ST.bundle.Dhyd_o,
              N_tubes=ST.N_tubes,
              d_s=ST.d_s,
              d_b=ST.bundle.d_b,
              Dhyd=ST.d_s,
              layout=ST.layout,
              l_b=ST.l_b,
              N_ss=5,
              eta_wall=ST.shell.state[1].eta*ones(ST.Ncell),
              N_passes=ST.N_passes,
              N_baffles_d=ST.N_baffles_d,
              N_baffles=ST.N_baffles,
              ff=1.0,
              b_cut=0.1588),
            redeclare function bundle_clearance =
                Design.Miscellanea.Shell_clearance.Fixed_Utube,
            N_tubes_start=800,
            pin_s=-1,
            pin_t=1,
            t_s_in_start=298.15,
            t_s_out_start=313.15,
            t_t_in_start=368.15,
            t_t_out_start=313.15,
            p_s_start=100000,
            p_t_start=320000,
            m_s_start=68.9,
            m_t_start=27.8)
            annotation (Placement(transformation(extent={{-42,-53},{106,73}})));
          Design.Components.Flags.ADDCO HOT_OUT(
            use_T=true,
            redeclare package Medium = Medium_t,
            T=313.15) annotation (Placement(transformation(extent={{-99,50},{-70,78}})));
          Design.Components.Flags.ADDCO COLD_IN(
            m_flow=68.9,
            use_T=true,
            use_p=true,
            use_m_flow=true,
            redeclare package Medium = Medium_s,
            T=298.15,
            p=100000)
            annotation (Placement(transformation(extent={{9,59},{38,87}})));
          Design.Components.Flags.ADDCO COLD_OUT(
            redeclare package Medium = Medium_s,
            T=313.15,
            use_T=false)
            annotation (Placement(transformation(extent={{-99,-97},{-70,-69}})));
          Design.Components.Flags.ADDCO HOT_IN(
            use_T=true,
            use_p=true,
            use_m_flow=true,
            redeclare package Medium = Medium_t,
            m_flow=27.8,
            T=368.15,
            p=3.2e5) annotation (Placement(transformation(extent={{-99,-23},{-70,5}})));
        equation
          connect(HOT_IN.node, ST.tube_in) annotation (Line(
              points={{-69.855,-9},{-56.9275,-9},{-56.9275,-8.9},{-42,-8.9}},
              color={0,0,0},
              pattern=LinePattern.None,
              smooth=Smooth.None));
          connect(HOT_OUT.node, ST.tube_out) annotation (Line(
              points={{-69.855,64},{-42,64},{-42,29.215}},
              color={0,0,0},
              pattern=LinePattern.None,
              smooth=Smooth.None));
          connect(COLD_OUT.node, ST.shell_out) annotation (Line(
              points={{-69.855,-83},{-5,-83},{-5,-53}},
              color={0,0,0},
              pattern=LinePattern.None,
              smooth=Smooth.None));
          connect(COLD_IN.node, ST.shell_in) annotation (Line(
              points={{38.145,73},{83.8,73}},
              color={0,0,0},
              pattern=LinePattern.None,
              smooth=Smooth.None));
         annotation (Placement(transformation(extent={{-108,-74},{88,66}})), Icon(
                coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{
                    100,100}}), graphics),
            Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,
                    100}}), graphics));
        end Aspen_Bell_Delaware_inv;
      end Single_phase;

      package Condensation "Tests for condensers"

        model Aspen_Kern
          "Verification with the results given by Aspen using the Kern method."
          package Medium_s = Media.RefProp.MM "Medium model";
          package Medium_t = Media.CoolProp.Water "Medium model";
          parameter Medium_s.SaturationProperties sat =  Medium_s.setSat_p(0.5e5)
            "Saturation properties";

           Design.Components.HEX.shell_and_tube ST(
            redeclare package Medium_s = Medium_s,
            redeclare package Medium_t = Medium_t,
            Dhyd=17.05e-3,
            l=1.8,
            pitch_f=1.25,
            layout=1,
            N_passes=4,
            N_baffles=5,
            N_baffles_d=6,
            redeclare Design.Heat_transfer.Shell.condensation_Kern hT_shell(
              Dhyd_o=ST.bundle.Dhyd_o,
              l=ST.l,
              N_tubes=ST.N_tubes,
              shear_vapor=true,
              pitch_f=ST.pitch_f,
              p_in = 0.5e5),
            redeclare Design.Pressure_drops.Shell.condensation_Kern dp_shell(
              l=ST.l/ST.N_baffles_d/ST.N_passes,
              d_s=ST.d_s,
              X=0.4,
              Dhyd=ST.bundle.Dhyd_o,
              l_b=ST.l_b,
              p_in = 0.5e5,
              eta_wall=ST.shell.state[1].eta*ones(ST.Ncell)),
            redeclare function bundle_clearance =
                Design.Miscellanea.Shell_clearance.Fixed_Utube,
            redeclare Design.Materials.S_AISI_1040 Material_t,
            redeclare Design.Materials.S_AISI_1040 Material_s,
            thick_t=2e-3,
            thick_s=10e-3,
            ht_t_f1=1e9,
            ht_s_f1=1e9,
            N_tubes_start=90,
            m_t_start=4.36,
            m_s_start=1.5365,
            t_s_in_start=373.15,
            t_s_out_start=sat.Tsat,
            t_t_in_start=313.15,
            t_t_out_start=333.15,
            p_s_start=50000,
            p_t_start=100000)
            "Shell and tube model of a condenser modeled with the Kern method."
            annotation (Placement(transformation(extent={{-68,-69},{80,57}})));

          Design.Components.Flags.ADDCO COLD_IN(
            m_flow=4.36,
            use_T=true,
            use_p=true,
            use_m_flow=true,
            redeclare package Medium = Medium_t,
            T=313.15,
            p=100000)
            annotation (Placement(transformation(extent={{-111,-39},{-82,-11}})));
          Design.Components.Flags.ADDCO HOT_OUT(
            use_h=true,
            redeclare package Medium = Medium_s,
            h= sat.hl) annotation (Placement(transformation(extent={{-111,-100},{-82,-72}})));
          Design.Components.Flags.ADDCO COLD_OUT(
            redeclare package Medium = Medium_t,
            T=333.15,
            use_T=false)
            annotation (Placement(transformation(extent={{-111,52},{-82,80}})));
          Design.Components.Flags.ADDCO HOT_IN(
            use_T=true,
            use_p=true,
            use_m_flow=true,
            redeclare package Medium = Medium_s,
            m_flow=1.5365,
            T= 99.1 + 273.15,
            p=0.5e5) annotation (Placement(transformation(extent={{9,81},{38,109}})));
        equation
          connect(COLD_IN.node,ST. tube_in) annotation (Line(
              points={{-81.855,-25},{-75.9275,-25},{-75.9275,-24.9},{-68,-24.9}},
              color={0,0,0},
              pattern=LinePattern.None,
              smooth=Smooth.None));
          connect(HOT_OUT.node,ST. shell_out) annotation (Line(
              points={{-81.855,-86},{-31,-86},{-31,-69}},
              color={0,0,0},
              pattern=LinePattern.None,
              smooth=Smooth.None));
          connect(HOT_IN.node,ST. shell_in) annotation (Line(
              points={{38.145,95},{57.8,95},{57.8,57}},
              color={0,0,0},
              pattern=LinePattern.None,
              smooth=Smooth.None));
          connect(COLD_OUT.node,ST. tube_out) annotation (Line(
              points={{-81.855,66},{-68,66},{-68,13.215}},
              color={0,0,0},
              pattern=LinePattern.None,
              smooth=Smooth.None));
         annotation (Placement(transformation(extent={{-108,-74},{88,66}})), Diagram(
                coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,100}}),
                graphics));
        end Aspen_Kern;

        model Aspen_Bell_Delaware
          "Verification with the results given by Aspen using the Bell Delaware method."

          package Medium_s = Media.RefProp.MM "Medium model";
          package Medium_t = Media.CoolProp.Water "Medium model";
          parameter Medium_s.SaturationProperties sat =  Medium_s.setSat_p(0.5e5)
            "Saturation properties";

           Design.Components.HEX.shell_and_tube ST(
            redeclare package Medium_s = Medium_s,
            redeclare package Medium_t = Medium_t,
            Dhyd=17.05e-3,
            l=1.8,
            pitch_f=1.25,
            layout=1,
            N_passes=4,
            N_baffles=5,
            N_baffles_d=6,
            m_t_start=4.36,
            m_s_start=1.5365,
            t_s_in_start=373.15,
            t_s_out_start=sat.Tsat,
            t_t_in_start=313.15,
            t_t_out_start=333.15,
            p_s_start=50000,
            p_t_start=100000,
            redeclare Design.Heat_transfer.Shell.condensation_Bell_Delaware
              hT_shell(
              pitch_f=ST.pitch_f,
              Dhyd_o=ST.bundle.Dhyd_o,
              N_tubes=ST.N_tubes,
              d_s=ST.d_s,
              d_b=ST.bundle.d_b,
              Dhyd=ST.d_s,
              l=ST.l,
              bts=3.18e-3,
              layout=ST.layout,
              l_b=ST.l_b,
              N_ss=2,
              p_in = 0.5e5,
              N_passes=ST.N_passes,
              ff=1.0,
              b_cut=0.3912),
            redeclare Design.Pressure_drops.Shell.condensation_Bell_Delaware
              dp_shell(
              pitch_f=ST.pitch_f,
              Dhyd_o=ST.bundle.Dhyd_o,
              N_tubes=ST.N_tubes,
              d_s=ST.d_s,
              d_b=ST.bundle.d_b,
              Dhyd=ST.d_s,
              layout=ST.layout,
              l_b=ST.l_b,
              N_ss=2,
              N_passes=ST.N_passes,
              N_baffles_d=ST.N_baffles_d,
              N_baffles=ST.N_baffles,
              ff=1.0,
              p_in = 0.5e5,
              b_cut=0.3912,
              bts=3.18e-3),
            redeclare function bundle_clearance =
                Design.Miscellanea.Shell_clearance.Fixed_Utube,
            redeclare Design.Materials.S_AISI_1040 Material_t,
            redeclare Design.Materials.S_AISI_1040 Material_s,
            thick_t=2e-3,
            thick_s=10e-3,
            N_tubes_start = 100,
            ht_t_f1=1e9,
            ht_s_f1=1e9,
            use_dp=true)
            "Shell and tube model of a condenser modeled with the Bell method."
            annotation (Placement(transformation(extent={{-80,-68},{68,58}})));

          Design.Components.Flags.ADDCO COLD_IN(
            m_flow=4.36,
            use_T=true,
            use_p=true,
            use_m_flow=true,
            redeclare package Medium = Medium_t,
            T=313.15,
            p=100000)
            annotation (Placement(transformation(extent={{-133,-38},{-104,-10}})));
          Design.Components.Flags.ADDCO HOT_OUT(
            use_h=true,
            redeclare package Medium = Medium_s,
            h= sat.hl) annotation (Placement(transformation(extent={{-133,-99},{-104,-71}})));
          Design.Components.Flags.ADDCO COLD_OUT(
            redeclare package Medium = Medium_t,
            T=333.15,
            use_T=true)
            annotation (Placement(transformation(extent={{-133,53},{-104,81}})));
          Design.Components.Flags.ADDCO HOT_IN(
            use_p=true,
            use_m_flow=true,
            redeclare package Medium = Medium_s,
            m_flow=1.5365,
            T= 99.1 + 273.15,
            p=50000,
            use_T=false)
                     annotation (Placement(transformation(extent={{-13,84},{16,112}})));
        equation
          connect(COLD_IN.node,ST. tube_in) annotation (Line(
              points={{-103.855,-24},{-92,-24},{-92,-23.9},{-80,-23.9}},
              color={0,0,0},
              pattern=LinePattern.None,
              smooth=Smooth.None));
          connect(HOT_OUT.node,ST. shell_out) annotation (Line(
              points={{-103.855,-85},{-43,-85},{-43,-68}},
              color={0,0,0},
              pattern=LinePattern.None,
              smooth=Smooth.None));
          connect(COLD_OUT.node,ST. tube_out) annotation (Line(
              points={{-103.855,67},{-80,67},{-80,14.215}},
              color={0,0,0},
              pattern=LinePattern.None,
              smooth=Smooth.None));
          connect(HOT_IN.node,ST. shell_in) annotation (Line(
              points={{16.145,98},{45.8,98},{45.8,58}},
              color={0,0,0},
              pattern=LinePattern.None,
              smooth=Smooth.None));
         annotation (Placement(transformation(extent={{-108,-74},{88,66}})),
            experiment,
            __Dymola_experimentSetupOutput,
            Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-180,
                    -120},{100,140}}),
                            graphics),
            Icon(coordinateSystem(extent={{-180,-120},{100,140}})));
        end Aspen_Bell_Delaware;

        model Aspen_Bell_Delaware_od
          "Verification with the results given by Aspen using the Bell Delaware method. Off-design mode is on."

          package Medium_s = Media.RefProp.MM "Medium model";
          package Medium_t = Media.CoolProp.Water "Medium model";
          parameter Medium_s.SaturationProperties sat =  Medium_s.setSat_p(0.5e5)
            "Saturation properties";

           Design.Components.HEX.shell_and_tube ST(
            redeclare package Medium_s = Medium_s,
            redeclare package Medium_t = Medium_t,
            Dhyd=17.05e-3,
            l=1.8,
            pitch_f=1.25,
            layout=1,
            N_passes=4,
            N_baffles=5,
            N_baffles_d=6,
            m_t_start=4.36,
            m_s_start=1.5365,
            t_s_in_start=373.15,
            t_s_out_start=sat.Tsat,
            t_t_in_start=313.15,
            t_t_out_start=333.15,
            p_s_start=50000,
            p_t_start=100000,
            redeclare Design.Heat_transfer.Shell.condensation_Bell_Delaware
              hT_shell(
              pitch_f=ST.pitch_f,
              Dhyd_o=ST.bundle.Dhyd_o,
              N_tubes=ST.N_tubes,
              d_s=ST.d_s,
              d_b=ST.bundle.d_b,
              Dhyd=ST.d_s,
              l=ST.l,
              bts=3.18e-3,
              layout=ST.layout,
              l_b=ST.l_b,
              N_ss=2,
              p_in = 0.5e5,
              N_passes=ST.N_passes,
              ff=1.0,
              b_cut=0.3912),
            redeclare Design.Pressure_drops.Shell.condensation_Bell_Delaware
              dp_shell(
              pitch_f=ST.pitch_f,
              Dhyd_o=ST.bundle.Dhyd_o,
              N_tubes=ST.N_tubes,
              d_s=ST.d_s,
              d_b=ST.bundle.d_b,
              Dhyd=ST.d_s,
              layout=ST.layout,
              l_b=ST.l_b,
              N_ss=2,
              N_passes=ST.N_passes,
              N_baffles_d=ST.N_baffles_d,
              N_baffles=ST.N_baffles,
              ff=1.0,
              p_in = 0.5e5,
              b_cut=0.3912,
              bts=3.18e-3),
            redeclare function bundle_clearance =
                Design.Miscellanea.Shell_clearance.Fixed_Utube,
            redeclare Design.Materials.S_AISI_1040 Material_t,
            redeclare Design.Materials.S_AISI_1040 Material_s,
            thick_t=2e-3,
            thick_s=10e-3,
            N_tubes_start = 100,
            ht_t_f1=1e9,
            ht_s_f1=1e9,
            offdesign=true,
            N_tubes_od=97)
            "Shell and tube model of a condenser modeled with the Bell method. Optimization mode is on."
            annotation (Placement(transformation(extent={{-78,-68},{70,58}})));

          Design.Components.Flags.ADDCO COLD_IN(
            m_flow=4.36,
            use_T=true,
            use_p=true,
            redeclare package Medium = Medium_t,
            use_m_flow=true,
            T=313.15,
            p=100000)
            annotation (Placement(transformation(extent={{-133,-38},{-104,-10}})));
          Design.Components.Flags.ADDCO HOT_OUT(
            redeclare package Medium = Medium_s,
            h= sat.hl,
            use_h=false)
                       annotation (Placement(transformation(extent={{-133,-99},{-104,-71}})));
          Design.Components.Flags.ADDCO COLD_OUT(
            redeclare package Medium = Medium_t,
            T=333.15,
            use_T=false)
            annotation (Placement(transformation(extent={{-135,53},{-106,81}})));
          Design.Components.Flags.ADDCO HOT_IN(
            use_p=true,
            use_m_flow=true,
            redeclare package Medium = Medium_s,
            use_T=true,
            T=372.25,
            p=50000,
            m_flow=0.9*1.5365)
                     annotation (Placement(transformation(extent={{-13,84},{16,112}})));
        equation
          connect(COLD_IN.node,ST. tube_in) annotation (Line(
              points={{-103.855,-24},{-92,-24},{-92,-23.9},{-78,-23.9}},
              color={0,0,0},
              pattern=LinePattern.None,
              smooth=Smooth.None));
          connect(HOT_OUT.node,ST. shell_out) annotation (Line(
              points={{-103.855,-85},{-41,-85},{-41,-68}},
              color={0,0,0},
              pattern=LinePattern.None,
              smooth=Smooth.None));
          connect(COLD_OUT.node,ST. tube_out) annotation (Line(
              points={{-105.855,67},{-78,67},{-78,14.215}},
              color={0,0,0},
              pattern=LinePattern.None,
              smooth=Smooth.None));
          connect(HOT_IN.node,ST. shell_in) annotation (Line(
              points={{16.145,98},{47.8,98},{47.8,58}},
              color={0,0,0},
              pattern=LinePattern.None,
              smooth=Smooth.None));
         annotation (Placement(transformation(extent={{-108,-74},{88,66}})),
            experiment,
            __Dymola_experimentSetupOutput,
            Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-180,-120},{100,
                    140}}), graphics),
            Icon(coordinateSystem(extent={{-180,-120},{100,140}})));
        end Aspen_Bell_Delaware_od;
      end Condensation;
    end Shell_and_tube;

    package Flat_plate "Tests for flat plate heat exchangers"

      package Single_phase "Tests for one phase to one phase heat exchangers"
        model Coulson "Verification with the results given by Coulson et al."

          package Medium_hot = Media.RefProp.Methanol "Medium model";
          package Medium_cold = Media.CoolProp.Water_TTSE "Medium model";

            Design.Components.HEX.Flat_plate FP(
            b=3e-3,
            w=0.5,
            X=9e-3,
            ht_hot_f1=1e4,
            ht_cold_f1=6e3,
            thick=0.75e-3,
            redeclare Design.Heat_transfer.Plates.single_phase_Coulson ht_cold,
            redeclare Design.Heat_transfer.Plates.single_phase_Coulson ht_hot,
            redeclare Design.Materials.T_R50250 material,
            d_pt=0.1,
            N_cell_pc=2,
            redeclare function cost = Design.Miscellanea.Cost.FP_Rafferty,
            N_ch_p=60,
            redeclare Design.Pressure_drops.Plates.single_phase_Coulson dp_hot,
            redeclare Design.Pressure_drops.Plates.single_phase_Coulson dp_cold,
            redeclare package Medium_hot = Test.Media.RefProp.Methanol,
            redeclare package Medium_cold = Test.Media.CoolProp.Water_TTSE,
            l_start=0.6,
            beta=0.87266462599716,
            t_hot_in_start=368.15,
            t_hot_out_start=313.15,
            t_cold_in_start=283.15,
            t_cold_out_start=313.15,
            p_hot_start=320000,
            p_cold_start=320000,
            mdot_hot_start=27.8,
            mdot_cold_start=68.9)
            "Model of a flat plate heat exchanger. Water cools down Methanol from saturated liquid."
            annotation (Placement(transformation(extent={{-92,-74},{74,48}})));

          Design.Components.Flags.ADDCO HOT_OUT(
            use_T=true,
            redeclare package Medium = Test.Media.RefProp.Methanol,
            T=313.15) annotation (Placement(transformation(extent={{-121,4},{-92,32}})));
          Design.Components.Flags.ADDCO COLD_IN(
            m_flow=68.9,
            use_T=true,
            use_p=true,
            use_m_flow=true,
            redeclare package Medium = Test.Media.CoolProp.Water_TTSE,
            T=298.15,
            p=100000)
            annotation (Placement(transformation(extent={{-121,-104},{-92,-76}})));
          Design.Components.Flags.ADDCO COLD_OUT(
            redeclare package Medium = Test.Media.CoolProp.Water_TTSE,
            T=313.15,
            use_T=false)
            annotation (Placement(transformation(extent={{23,-114},{52,-86}})));
          Design.Components.Flags.ADDCO HOT_IN(
            use_T=true,
            use_p=true,
            use_m_flow=true,
            redeclare package Medium = Test.Media.RefProp.Methanol,
            m_flow=27.8,
            T=368.15,
            p=320000) annotation (Placement(transformation(extent={{63,69},{92,97}})));
        equation
          connect(HOT_OUT.node,FP. node_h_out) annotation (Line(
              points={{-91.855,18},{-91.855,-49.6},{-92,-49.6}},
              color={0,0,0},
              pattern=LinePattern.None,
              smooth=Smooth.None));
          connect(COLD_IN.node,FP. node_c_in) annotation (Line(
              points={{-91.855,-90},{-91.855,-81},{-92,-81},{-92,-68.51}},
              color={0,0,0},
              pattern=LinePattern.None,
              smooth=Smooth.None));
          connect(FP.node_c_out, COLD_OUT.node) annotation (Line(
              points={{74,17.195},{86,17.195},{86,16},{100,16},{100,-100},{52.145,-100}},
              color={0,0,0},
              pattern=LinePattern.None,
              smooth=Smooth.None));

          connect(HOT_IN.node,FP. node_h_in) annotation (Line(
              points={{92.145,83},{98,83},{98,35.8},{74,35.8}},
              color={0,0,0},
              pattern=LinePattern.None,
              smooth=Smooth.None));

         annotation (Placement(transformation(extent={{-108,-74},{88,66}})),
            experiment(__Dymola_NumberOfIntervals=1, Tolerance=1e-006),
            __Dymola_experimentSetupOutput,
            Diagram(coordinateSystem(extent={{-160,-120},{100,100}},
                  preserveAspectRatio=false), graphics),
            Icon(coordinateSystem(extent={{-160,-120},{100,100}})));
        end Coulson;

        model Rossetto "Verification with the results given by Rossetto et al."

          package Medium_hot = Media.CoolProp.Water_TTSE "Medium model";
          package Medium_cold = Media.CoolProp.Water_TTSE "Medium model";

          Design.Components.HEX.Flat_plate FP(
            redeclare package Medium_hot = Medium_hot,
            redeclare package Medium_cold = Medium_cold,
            X=1.025,
            redeclare function cost = Design.Miscellanea.Cost.FP_Rafferty,
            redeclare Design.Materials.SS_AISI_410 material,
            thick=0.4e-3,
            b=3.4e-3,
            w=0.65,
            d_pt=0.25,
            ht_hot_f1=2e4,
            ht_cold_f1=1e12,
            N_cell_pc=2,
            N_ch_p=204,
            l_start=1,
            beta=0.78539816339745,
            t_hot_in_start=339.15,
            t_hot_out_start=306.15,
            t_cold_in_start=288.15,
            t_cold_out_start=300.15,
            p_hot_start=100000,
            p_cold_start=100000,
            mdot_hot_start=140,
            mdot_cold_start=140)
            "\"Model of a flat plate heat exchanger. Water cools down water.\""
            annotation (Placement(transformation(extent={{-82,-66},{86,54}})));

          Design.Components.Flags.ADDCO HOT_OUT(
            use_T=true,
            redeclare package Medium = Media.CoolProp.Water_TTSE,
            T=306.15) annotation (Placement(transformation(extent={{-111,14},{-82,42}})));
          Design.Components.Flags.ADDCO HOT_IN(
            use_T=true,
            use_p=true,
            use_m_flow=true,
            redeclare package Medium = Media.CoolProp.Water_TTSE,
            T=333.15,
            p=100000,
            m_flow=140)
                      annotation (Placement(transformation(extent={{73,79},{102,107}})));
          Design.Components.Flags.ADDCO COLD_IN(
            use_T=true,
            use_p=true,
            use_m_flow=true,
            redeclare package Medium = Test.Media.CoolProp.Water_TTSE,
            T=293.15,
            p=100000,
            m_flow=140)
            annotation (Placement(transformation(extent={{-111,-94},{-82,-66}})));
          Design.Components.Flags.ADDCO COLD_OUT(
            redeclare package Medium = Test.Media.CoolProp.Water_TTSE,
            use_T=false,
            T=320.15)
            annotation (Placement(transformation(extent={{33,-104},{62,-76}})));
        equation
          connect(COLD_IN.node,FP. node_c_in) annotation (Line(
              points={{-81.855,-80},{-82,-80},{-82,-60.6}},
              color={0,0,0},
              pattern=LinePattern.None,
              smooth=Smooth.None));
          connect(HOT_OUT.node,FP. node_h_out) annotation (Line(
              points={{-81.855,28},{-82,28},{-82,-42}},
              color={0,0,0},
              pattern=LinePattern.None,
              smooth=Smooth.None));
          connect(HOT_IN.node,FP. node_h_in) annotation (Line(
              points={{102.145,93},{122,93},{122,42},{86,42}},
              color={0,0,0},
              pattern=LinePattern.None,
              smooth=Smooth.None));
          connect(COLD_OUT.node,FP. node_c_out) annotation (Line(
              points={{62.145,-90},{120,-90},{120,23.7},{86,23.7}},
              color={0,0,0},
              pattern=LinePattern.None,
              smooth=Smooth.None));
         annotation (Placement(transformation(extent={{-108,-74},{88,66}})),
            experiment(__Dymola_NumberOfIntervals=1),
            __Dymola_experimentSetupOutput,
            Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-160,-120},{140,
                    120}}), graphics),
            Icon(coordinateSystem(extent={{-160,-120},{140,120}})));
        end Rossetto;

        model Aspen_rec
          "Verification with the results given by Aspen. We test recuperation in this example."
          package Medium_hot = Media.FluidProp.MM "Medium model";
          package Medium_cold = Media.FluidProp.MM "Medium model";
          parameter Medium_hot.ThermodynamicState hot_in = Medium_hot.setState_pT(0.33e5, 276.9 + 273.15)
            "Thermodynamic state at the inlet of the hot side";
          parameter Medium_hot.ThermodynamicState hot_out = Medium_hot.setState_pT(0.33e5, 71.4 + 273.15)
            "Thermodynamic state at the outlet of the hot side";
           parameter Medium_cold.ThermodynamicState cold_in = Medium_cold.setState_pT(14.6e5, 51.4 + 273.15)
            "Thermodynamic state at the inlet of the hot side";
           parameter Medium_cold.ThermodynamicState cold_out = Medium_cold.setState_pT(14.6e5, 222.9 + 273.15)
            "Thermodynamic state at the outlet of the hot side";

             Design.Components.HEX.Flat_plate fp(
            redeclare package Medium_hot = Medium_hot,
            redeclare package Medium_cold = Medium_cold,
            redeclare function cost = Design.Miscellanea.Cost.FP_Rafferty,
            ht_hot_f1=1e9,
            ht_cold_f1=1e9,
            X=1,
            redeclare Design.Materials.SS_AISI_304 material,
            N_ch_p=18,
            thick=0.9e-3,
            b=3.12e-3,
            w=325e-3,
            d_pt=100e-3,
            redeclare Design.Miscellanea.topology_PHE.two_pass_two_pass tpg_hot(N_ch_p=
                  fp.N_ch_p, stype=1),
            redeclare Design.Miscellanea.topology_PHE.two_pass_two_pass tpg_cold(N_ch_p=
                 fp.N_ch_p, stype=2),
            N_cell_pc=3,
            beta=1.0594148559606,
            t_hot_in_start=550.15,
            t_hot_out_start=343.15,
            t_cold_in_start=333.15,
            t_cold_out_start=473.15,
            p_hot_start=33000,
            p_cold_start=1460000,
            mdot_hot_start=0.15,
            mdot_cold_start=0.15)
            annotation (Placement(transformation(extent={{-82,-66},{86,54}})));

          Design.Components.Flags.ADDCO HOT_OUT(
            use_T=true,
            redeclare package Medium = Media.FluidProp.MM,
            T=344.55) annotation (Placement(transformation(extent={{-112,2},{-83,30}})));
          Design.Components.Flags.ADDCO COLD_IN(
            use_T=true,
            use_p=true,
            use_m_flow=true,
            redeclare package Medium = Test.Media.FluidProp.MM,
            T=324.55,
            p=1460000,
            m_flow=0.15)
            annotation (Placement(transformation(extent={{-112,-106},{-83,-78}})));
          Design.Components.Flags.ADDCO COLD_OUT(
            redeclare package Medium = Test.Media.FluidProp.MM,
            use_T=false,
            T=496.05)
            annotation (Placement(transformation(extent={{32,-116},{64,-86}})));
          Design.Components.Flags.ADDCO HOT_IN(
            use_T=true,
            use_p=true,
            use_m_flow=true,
            redeclare package Medium = Media.FluidProp.MM,
            T=550.05,
            p=33000,
            m_flow=0.15)
                      annotation (Placement(transformation(extent={{72,67},{101,95}})));
        equation
          connect(fp.node_c_in, COLD_IN.node) annotation (Line(
              points={{-82,-60.6},{-82,-92},{-82.855,-92}},
              color={0,0,0},
              pattern=LinePattern.None,
              smooth=Smooth.None));
          connect(HOT_OUT.node, fp.node_h_out) annotation (Line(
              points={{-82.855,16},{-82,16},{-82,-42}},
              color={0,0,0},
              pattern=LinePattern.None,
              smooth=Smooth.None));
          connect(COLD_OUT.node, fp.node_c_out) annotation (Line(
              points={{64.16,-101},{122,-101},{122,23.7},{86,23.7}},
              color={0,0,0},
              pattern=LinePattern.None,
              smooth=Smooth.None));
          connect(HOT_IN.node, fp.node_h_in) annotation (Line(
              points={{101.145,81},{114,81},{114,42},{86,42}},
              color={0,0,0},
              pattern=LinePattern.None,
              smooth=Smooth.None));
         annotation (Placement(transformation(extent={{-108,-74},{88,66}})),
            experiment(__Dymola_NumberOfIntervals=1),
            __Dymola_experimentSetupOutput,
            Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,
                    100}}), graphics));
        end Aspen_rec;
      end Single_phase;

      package Evaporation "Tests for evaporators"

        model Aspen_eva
          "Verification with the results given by Aspen. We test evaporation in this example."
          package Medium_hot =
               Media.CoolProp.THERM66 "Medium model";
          package Medium_cold = Media.RefProp.MM "Medium model";
          parameter Medium_hot.ThermodynamicState hot_in = Medium_hot.setState_pT(2e5, 345 + 273.15)
            "Thermodynamic state at the inlet of the hot side";
          parameter Medium_hot.ThermodynamicState hot_out = Medium_hot.setState_pT(2e5, 231.9 + 273.15)
            "Thermodynamic state at the outlet of the hot side";
           parameter Medium_cold.ThermodynamicState cold_in = Medium_cold.setState_pT(14.6e5, 216.9 + 273.15)
            "Thermodynamic state at the inlet of the hot side";
           parameter Medium_cold.ThermodynamicState cold_out = Medium_cold.setState_pT(14.6e5, 325 + 273.15)
            "Thermodynamic state at the outlet of the hot side";

            Design.Components.HEX.Flat_plate FP(
            redeclare package Medium_hot = Medium_hot,
            redeclare package Medium_cold = Medium_cold,
            redeclare function cost = Design.Miscellanea.Cost.FP_Rafferty,
            thick=0.8e-3,
            d_pt=32e-3,
            ht_hot_f1=1e9,
            ht_cold_f1=1e9,
            X=1,
            N_ch_p=14,
            w=0.1325,
            redeclare Design.Heat_transfer.Plates.evaporation_Martin ht_cold(
              p_in=FP.node_c_in.p,
              At=FP.At,
              beta=FP.beta,
              qdot=FP.plate.qdot_cold,
              qdot_tilde_start=FP.mdot_hot_start*(FP.h_hot_in_start - FP.h_hot_out_start)/(2*FP.l_start*FP.w*FP.N_ch_p)),
            redeclare Design.Miscellanea.topology_PHE.two_pass_two_pass tpg_hot(N_ch_p=
                  FP.N_ch_p, stype=1),
            redeclare Design.Miscellanea.topology_PHE.two_pass_two_pass tpg_cold(N_ch_p=
                 FP.N_ch_p, stype=2),
            redeclare Design.Materials.SS_AISI_304 material,
            b=2.3e-3,
            redeclare Design.Pressure_drops.Plates.evaporation_Martin dp_cold(
              p_in=FP.node_c_in.p,
              beta=FP.beta),
            redeclare Design.Miscellanea.check_velocity check_hot(T_sat=500),
            redeclare Design.Miscellanea.check_velocity check_cold(umin=min(FP.ht_cold.u)),
            h_hot_in_start=hot_in.h,
            h_hot_out_start=hot_out.h,
            l_start=0.5,
            N_cell_pc=10,
            offdesign=false,
            mdot_cold_start=0.15,
            beta=1.0471975511966,
            t_hot_in_start=618.15,
            t_hot_out_start=513.15,
            t_cold_in_start=500.15,
            t_cold_out_start=593.15,
            p_hot_start=200000,
            p_cold_start=1460000,
            mdot_hot_start=0.19)
            "Model of a flat plate heat exchanger. Therminol 66 evaporates MM."
            annotation (Placement(transformation(extent={{-90,-82},{86,60}})));

          Design.Components.Flags.ADDCO HOT_OUT(
            redeclare package Medium = Test.Media.CoolProp.THERM66,
            T=505.05,
            use_T=false)
            annotation (Placement(transformation(extent={{-133,-8},{-100,22}})));
          Design.Components.Flags.ADDCO COLD_IN(
            use_T=true,
            use_p=true,
            redeclare package Medium = Test.Media.RefProp.MM,
            T=490.05,
            p=1460000,
            m_flow=0.15,
            use_m_flow=true)
            annotation (Placement(transformation(extent={{-133,-116},{-104,-88}})));
          Design.Components.Flags.ADDCO COLD_OUT(
            redeclare package Medium = Test.Media.RefProp.MM,
            T=598.15,
            use_T=true)
            annotation (Placement(transformation(extent={{79,-116},{108,-88}})));
          Design.Components.Flags.ADDCO HOT_IN(
            use_T=true,
            use_p=true,
            use_m_flow=true,
            redeclare package Medium = Test.Media.CoolProp.THERM66,
            T=618.15,
            p=200000,
            m_flow=0.19)
                      annotation (Placement(transformation(extent={{49,100},{78,128}})));
        equation
          connect(HOT_OUT.node,FP. node_h_out) annotation (Line(
              points={{-99.835,7},{-99.835,-53.6},{-90,-53.6}},
              color={0,0,0},
              pattern=LinePattern.None,
              smooth=Smooth.None));
          connect(COLD_IN.node,FP. node_c_in) annotation (Line(
              points={{-103.855,-102},{-104,-102},{-104,-75.61},{-90,-75.61}},
              color={0,0,0},
              pattern=LinePattern.None,
              smooth=Smooth.None));
          connect(FP.node_c_out, COLD_OUT.node) annotation (Line(
              points={{86,24.145},{130,24.145},{130,-102},{108.145,-102}},
              color={0,0,0},
              pattern=LinePattern.None,
              smooth=Smooth.None));
          connect(HOT_IN.node,FP. node_h_in) annotation (Line(
              points={{78.145,114},{130,114},{130,45.8},{86,45.8}},
              color={0,0,0},
              pattern=LinePattern.None,
              smooth=Smooth.None));
         annotation (Placement(transformation(extent={{-108,-74},{88,66}})),
            experiment(__Dymola_NumberOfIntervals=1),
            __Dymola_experimentSetupOutput,
            Diagram(coordinateSystem(extent={{-160,-120},{140,140}},
                  preserveAspectRatio=false), graphics),
            Icon(coordinateSystem(extent={{-160,-120},{140,140}})));
        end Aspen_eva;

        model Aspen_eva_od
          "Verification with the results given by Aspen. We test evaporation in this example. Off-design mode on."
          package Medium_hot =
               Media.CoolProp.THERM66 "Medium model";
          package Medium_cold = Media.RefProp.MM "Medium model";
          parameter Medium_hot.ThermodynamicState hot_in = Medium_hot.setState_pT(2e5, 345 + 273.15)
            "Thermodynamic state at the inlet of the hot side";
          parameter Medium_hot.ThermodynamicState hot_out = Medium_hot.setState_pT(2e5, 231.9 + 273.15)
            "Thermodynamic state at the outlet of the hot side";

            Design.Components.HEX.Flat_plate FP(
            redeclare package Medium_hot = Medium_hot,
            redeclare package Medium_cold = Medium_cold,
            redeclare function cost = Design.Miscellanea.Cost.FP_Rafferty,
            thick=0.8e-3,
            d_pt=32e-3,
            ht_hot_f1=1e9,
            ht_cold_f1=1e9,
            X=1,
            N_ch_p=14,
            w=0.1325,
            redeclare Design.Heat_transfer.Plates.evaporation_Martin ht_cold(
              p_in=FP.node_c_in.p,
              At=FP.At,
              beta=FP.beta,
              qdot=FP.plate.qdot_cold,
              qdot_tilde_start=FP.mdot_hot_start*(FP.h_hot_in_start - FP.h_hot_out_start)/(2*FP.l_start*FP.w*FP.N_ch_p)),
            redeclare Design.Miscellanea.topology_PHE.two_pass_two_pass tpg_hot(N_ch_p=
                  FP.N_ch_p, stype=1),
            redeclare Design.Miscellanea.topology_PHE.two_pass_two_pass tpg_cold(N_ch_p=
                 FP.N_ch_p, stype=2),
            redeclare Design.Materials.SS_AISI_304 material,
            b=2.3e-3,
            redeclare Design.Pressure_drops.Plates.evaporation_Martin dp_cold(
              p_in=FP.node_c_in.p,
              beta=FP.beta),
            redeclare Design.Miscellanea.check_velocity check_hot(T_sat=500),
            redeclare Design.Miscellanea.check_velocity check_cold(umin=min(FP.ht_cold.u)),
            h_hot_in_start=hot_in.h,
            h_hot_out_start=hot_out.h,
            l_start=0.5,
            l_od=0.512939,
            offdesign=true,
            plate(mdot_cold(start=0.15/FP.N_ch_p)),
            N_cell_pc=7,
            mdot_hot_start=0.19,
            beta=1.0471975511966,
            t_hot_in_start=618.15,
            t_hot_out_start=513.15,
            t_cold_in_start=489.15,
            t_cold_out_start=598.15,
            p_hot_start=200000,
            p_cold_start=1460000,
            mdot_cold_start=0.14)
            "Model of a flat plate heat exchanger. Therminol 66 evaporates MM."
            annotation (Placement(transformation(extent={{-90,-82},{86,60}})));

          Design.Components.Flags.ADDCO HOT_OUT(
            redeclare package Medium = Test.Media.CoolProp.THERM66,
            T=505.05,
            use_T=false)
            annotation (Placement(transformation(extent={{-133,-8},{-100,22}})));
          Design.Components.Flags.ADDCO COLD_IN(
            use_T=true,
            use_p=true,
            redeclare package Medium = Test.Media.RefProp.MM,
            m_flow=0.14,
            use_m_flow=false,
            T=490.05,
            p=1460000)
            annotation (Placement(transformation(extent={{-133,-116},{-104,-88}})));
          Design.Components.Flags.ADDCO COLD_OUT(
            redeclare package Medium = Test.Media.RefProp.MM,
            use_T=true,
            T=598.15)
            annotation (Placement(transformation(extent={{77,-116},{106,-88}})));
          Design.Components.Flags.ADDCO HOT_IN(
            use_T=true,
            use_p=true,
            use_m_flow=true,
            redeclare package Medium = Test.Media.CoolProp.THERM66,
            T=618.15,
            p=200000,
            m_flow=0.19)
                      annotation (Placement(transformation(extent={{49,100},{78,128}})));
        equation
          connect(HOT_OUT.node,FP. node_h_out) annotation (Line(
              points={{-99.835,7},{-99.835,-53.6},{-90,-53.6}},
              color={0,0,0},
              pattern=LinePattern.None,
              smooth=Smooth.None));
          connect(COLD_IN.node,FP. node_c_in) annotation (Line(
              points={{-103.855,-102},{-104,-102},{-104,-75.61},{-90,-75.61}},
              color={0,0,0},
              pattern=LinePattern.None,
              smooth=Smooth.None));
          connect(FP.node_c_out, COLD_OUT.node) annotation (Line(
              points={{86,24.145},{130,24.145},{130,-102},{106.145,-102}},
              color={0,0,0},
              pattern=LinePattern.None,
              smooth=Smooth.None));
          connect(HOT_IN.node,FP. node_h_in) annotation (Line(
              points={{78.145,114},{130,114},{130,45.8},{86,45.8}},
              color={0,0,0},
              pattern=LinePattern.None,
              smooth=Smooth.None));
         annotation (Placement(transformation(extent={{-108,-74},{88,66}})),
            experiment(__Dymola_NumberOfIntervals=1),
            __Dymola_experimentSetupOutput,
            Diagram(coordinateSystem(extent={{-160,-120},{140,140}},
                  preserveAspectRatio=false), graphics),
            Icon(coordinateSystem(extent={{-160,-120},{140,140}})));
        end Aspen_eva_od;
      end Evaporation;

      package Condensation "Tests for condensers"

        model Aspen_con
          "Verification with the results given by Aspen. We test condensation in this example."
          replaceable package Medium = Media.CoolProp.Toluene_TTSE  constrainedby
            Modelica.Media.Interfaces.PartialMedium "Medium model hot cells" annotation(choicesAllMatching = true);
          replaceable package Coolant = Media.CoolProp.Water_TTSE  constrainedby
            Modelica.Media.Interfaces.PartialMedium "Medium model hot cells" annotation(choicesAllMatching = true);

            Design.Components.HEX.Flat_plate FP(
            redeclare package Medium_hot = Medium,
            redeclare package Medium_cold = Coolant,
            redeclare function cost = Design.Miscellanea.Cost.FP_Rafferty,
            ht_hot_f1=1e9,
            ht_cold_f1=1e9,
            w=325e-3,
            d_pt=100e-3,
            redeclare Design.Materials.SS_AISI_304 material,
            thick=0.9e-3,
            b=3.12e-3,
            redeclare Design.Miscellanea.check_velocity check_hot(T_sat=50+273.15),
            redeclare Design.Miscellanea.check_velocity check_cold(umin=min(FP.ht_cold.u)),
            l_start=0.5,
            offdesign=false,
            redeclare Design.Heat_transfer.Plates.condensation_Longo ht_hot(beta=FP.beta,
            p_in=FP.node_h_in.p),
            mdot_hot_start=0.22,
            mdot_cold_start=0.5,
            redeclare Design.Pressure_drops.Plates.evaporation_Martin dp_hot(
              p_in=FP.node_h_in.p,
              beta=FP.beta),
            N_ch_p=20,
            use_dp=false,
            N_cell_pc=12,
            X=1,
            beta=1.0471975511966,
            t_hot_in_start=343.15,
            t_hot_out_start=313.15,
            t_cold_in_start=278.15,
            t_cold_out_start=303.15,
            p_hot_start=10000,
            p_cold_start=100000)
            "Model of a flat plate heat exchanger. Therminol 66 evaporates MM."
            annotation (Placement(transformation(extent={{-92,-82},{84,60}})));

          Design.Components.Flags.ADDCO HOT_OUT(
            redeclare package Medium = Medium,
            use_T=true,
            use_p=true,
            T=313.15,
            p=10000)
            annotation (Placement(transformation(extent={{-133,-8},{-100,22}})));
          Design.Components.Flags.ADDCO COLD_IN(
            use_T=true,
            use_p=true,
            redeclare package Medium = Coolant,
            use_m_flow=true,
            T=278.15,
            p=100000,
            m_flow=0.8)
            annotation (Placement(transformation(extent={{-133,-116},{-104,-88}})));
          Design.Components.Flags.ADDCO COLD_OUT(
            redeclare package Medium = Coolant,
            T=598.15,
            use_T=false)
            annotation (Placement(transformation(extent={{79,-116},{108,-88}})));
          Design.Components.Flags.ADDCO HOT_IN(
            use_T=true,
            use_m_flow=true,
            redeclare package Medium = Medium,
            m_flow=0.22,
            T=343.15,
            p=10000,
            use_p=false)
                      annotation (Placement(transformation(extent={{49,100},{78,128}})));
        equation
          connect(HOT_OUT.node,FP. node_h_out) annotation (Line(
              points={{-99.835,7},{-99.835,-53.6},{-92,-53.6}},
              color={0,0,0},
              pattern=LinePattern.None,
              smooth=Smooth.None));
          connect(COLD_IN.node,FP. node_c_in) annotation (Line(
              points={{-103.855,-102},{-104,-102},{-104,-75.61},{-92,-75.61}},
              color={0,0,0},
              pattern=LinePattern.None,
              smooth=Smooth.None));
          connect(FP.node_c_out, COLD_OUT.node) annotation (Line(
              points={{84,24.145},{130,24.145},{130,-102},{108.145,-102}},
              color={0,0,0},
              pattern=LinePattern.None,
              smooth=Smooth.None));
          connect(HOT_IN.node,FP. node_h_in) annotation (Line(
              points={{78.145,114},{130,114},{130,45.8},{84,45.8}},
              color={0,0,0},
              pattern=LinePattern.None,
              smooth=Smooth.None));
         annotation (Placement(transformation(extent={{-108,-74},{88,66}})),
            experiment(__Dymola_NumberOfIntervals=1),
            __Dymola_experimentSetupOutput,
            Diagram(coordinateSystem(extent={{-160,-120},{140,140}},
                  preserveAspectRatio=false), graphics),
            Icon(coordinateSystem(extent={{-160,-120},{140,140}})));
        end Aspen_con;

      end Condensation;
    end Flat_plate;

    package Cell_method "Tests of the cell method for various configurations"
      model cell_method_VDI_example
        "Cell method verified with an example from VDI C1 pag. 48"
        parameter Integer N_passes = 2 "Number of cell passes";
        parameter Integer N_baffles = 3 "Number of cell baffles";
        parameter Integer Ncell = N_passes*N_baffles "Number of cell elements";
        parameter Modelica.SIunits.Temp_C t_s_in = 100
          "Inlet temperature of the hot source";
        parameter Modelica.SIunits.Temp_C t_t_in = 20
          "Inlet temperature of the cold source";
        parameter Modelica.SIunits.ThermalConductance kA_tot = 4749/Ncell
          "Thermal conductance of each cell";
        parameter Modelica.SIunits.HeatCapacity C_hot = 3.5e3
          "Capacity of the hot stream";
        parameter Modelica.SIunits.HeatCapacity C_cold = 3.5e3
          "Capacity of the cold stream";
        Modelica.SIunits.HeatFlowRate qdot_hot[Ncell] "Heat rate hot side";
        Modelica.SIunits.HeatFlowRate qdot_cold[Ncell] "Heat rate cold side";
        Modelica.SIunits.Temp_C t_hot[Ncell+1](start=linspace(t_s_in,40,Ncell+1))
          "Hot stream temperature array";
        Modelica.SIunits.Temp_C t_cold[Ncell+1](start=linspace(t_t_in,40,Ncell+1))
          "Cold stream temperature array";
        Modelica.SIunits.Temp_C t_av_hot[Ncell](start=linspace(t_s_in,40,Ncell))
          "Cell hot stream temperature array";
        Modelica.SIunits.Temp_C t_av_cold[Ncell](start=linspace(t_t_in,40,Ncell))
          "Cell cold stream temperature array";
        Real adim_hot[Ncell] "Adimensional temperatures of the hot stream";
        Real adim_cold[Ncell] "Adimensional temperatures of the cold stream";
        Integer k[Ncell] "Index for the shell to tube topology";

      equation
          for j in 1:Ncell loop
            qdot_hot[k[j]] = qdot_cold[j];
            qdot_cold[j]   = kA_tot*(t_av_hot[k[j]] - t_av_cold[j]);
            qdot_hot[j]    = C_hot*(t_hot[j] - t_hot[j+1]);
            t_av_hot[j]    = 0.5*(t_hot[j+1] + t_hot[j]);
            qdot_cold[j]   = C_cold*(t_cold[j+1] - t_cold[j]);
            t_av_cold[j]   = 0.5*(t_cold[j+1] + t_cold[j]);
            adim_hot[j]    = (t_hot[j] - t_t_in)/(t_s_in - t_t_in);
            adim_cold[j]   = (t_cold[j] - t_t_in)/(t_s_in - t_t_in);
          end for;

          //Things get tricky here. We need to get the index of the shell cell seen by the tube cell
          for j in 1:N_passes loop
               for i in 1+(j-1)*N_baffles:2:j*N_baffles loop
                  if mod(j,2) <> 0 then
                    if mod(N_baffles,2) <> 0 then
                      k[i]   = Ncell - N_passes*(i-(1+(j-1)*N_baffles)) - (j-1);
                    else
                      k[i]   = Ncell - N_passes*(i-(1+(j-1)*N_baffles)) - (N_passes
                      - j);
                    end if;
                  else
                      k[i]   = (N_passes + 1 - j) + N_passes*(i-(1+(j-1)*N_baffles));
                  end if;
               end for;
           end for;

           for j in 1:N_passes loop
               for i in 2+(j-1)*N_baffles:2:j*N_baffles loop
                  if mod(j,2) <> 0 then
                    if mod(N_baffles,2) <> 0 then
                      k[i]   = Ncell - (2*N_passes - 1) + (j-1) - N_passes*(i-(2+(j-1)
                      *N_baffles));
                    else
                      k[i]   = Ncell - N_passes - (j-1) - N_passes*(i-(2+(j-1)
                      *N_baffles));
                    end if;
                  else
                      k[i]   = N_passes + j + N_passes*(i-(2+(j-1)*N_baffles));
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

      model cell_method_flat_plate "Cell method for flat plate heat exchangers"
        parameter Integer layout = 2 "Flow path, 1 = parallel, 2 = series";
        parameter Integer N_passes_d = 10 "Number of passes hot and cold side";
        parameter Integer N_plates_d = 2*N_passes_d - 1
          "Discretized number of active plates";
        parameter Integer N_cell_pc = 5 "Number of cells per channels";
        parameter Integer Ncell = N_plates_d*N_cell_pc
          "Number of cell elements";
        parameter Modelica.SIunits.Temp_C t_hot_in = 90
          "Inlet temperature of the hot source";
        parameter Modelica.SIunits.Temp_C t_cold_in = 18
          "Inlet temperature of the cold source";
        parameter Modelica.SIunits.ThermalConductance kA_tot = 4185.5*4.9940*0.05/Ncell
          "Thermal conductance of each cell";
        parameter Real C_hot = 4185.5*0.05 "Capacity of the hot stream";
        parameter Real C_cold = 4185.5*0.05 "Capacity of the cold stream";
        final parameter Modelica.SIunits.Temp_C t_hot_start[N_cell_pc + 1]=
          linspace(t_hot_in, 40, N_cell_pc + 1);
        final parameter Modelica.SIunits.Temp_C t_cold_start[N_cell_pc + 1]=
          linspace(t_cold_in, 40, N_cell_pc + 1);
        Modelica.SIunits.Temp_C t_hot_out
          "Outlet temperature of the hot source";
        Modelica.SIunits.Temp_C t_cold_out
          "Outlet temperature of the cold source";
        Modelica.SIunits.HeatFlowRate qdot_tot "Total heat rate";
        Modelica.SIunits.HeatFlowRate qdot_hot[N_passes_d, N_cell_pc]
          "Heat rate hot side";
        Modelica.SIunits.HeatFlowRate qdot_cold[N_passes_d, N_cell_pc]
          "Heat rate cold side";
        Modelica.SIunits.HeatFlowRate qdot_wall[2, N_passes_d, N_cell_pc]
          "Heat rate metal wall";
        Modelica.SIunits.Temp_C t_hot[N_passes_d, N_cell_pc + 1](
         start=fill(t_hot_start, N_passes_d)) "Hot stream temperature matrix";
        Modelica.SIunits.Temp_C t_cold[N_passes_d, N_cell_pc + 1](
         start=fill(t_cold_start, N_passes_d)) "Cold stream temperature matrix";
        Modelica.SIunits.Temp_C t_av_hot[N_passes_d, N_cell_pc]
          "Cell hot stream temperature matrix";
        Modelica.SIunits.Temp_C t_av_cold[N_passes_d, N_cell_pc]
          "Cell cold stream temperature matrix";
        Modelica.SIunits.TemperatureDifference  DTML
          "Logarithmic mean temperature difference";
        Modelica.SIunits.TemperatureDifference  DTML_tilde
          "Logarithmic mean temperature difference corrected";
        Integer k1[N_passes_d] "Index for the flat plate topology columns";
        Integer k2[N_cell_pc] "Index for the flat plate topology rows";

      equation
            for j in 1:N_passes_d loop
              k1[j] = N_passes_d + 1 - j;
              for i in 1:N_cell_pc loop
                //Average temperature in the volumes
                t_av_hot[j, i]   = 0.5*(t_hot[j, i] + t_hot[j, i + 1]);
                t_av_cold[j, i]  = 0.5*(t_cold[j, i] + t_cold[j, i + 1]);

                if (layout == 1) then
                  //Heat flow rate cold and hot side
                  qdot_hot[j, i]  = C_hot*(t_hot[j, i] - t_hot[j, i + 1]);
                  qdot_cold[j, i] = C_cold*(t_cold[j, i + 1] - t_cold[j, i]);
                else
                  if (mod(j,2) <> 0) then
                  qdot_hot[j, i]  = C_hot*(t_hot[j, i] - t_hot[j, i + 1]);
                  qdot_cold[j, i] = C_cold*(t_cold[j, i + 1] - t_cold[j, i]);
                  else
                  qdot_hot[j, i]  = C_hot*(t_hot[j, i + 1] - t_hot[j, i]);
                  qdot_cold[j, i] = C_cold*(t_cold[j, i] - t_cold[j, i + 1]);
                  end if;
                end if;
                //Boundary conditions
                if (i == 1) then
                  if (layout == 1) then
                    t_hot[j, i]  = t_hot_in;
                    t_cold[j, i] = t_cold_in;
                  else
                    if (j == 1) then
                      t_hot_in   = t_hot[1, 1];
                      t_cold_in  = t_cold[1, 1];
                      else
                      if (mod(j,2) <> 0) then
                        t_hot[j, 1]              = t_hot[j - 1, 1];
                        t_cold[j, 1]             = t_cold[j - 1, 1];
                      else
                        t_hot[j, N_cell_pc + 1]  = t_hot[j - 1, N_cell_pc + 1];
                        t_cold[j, N_cell_pc + 1] = t_cold[j - 1, N_cell_pc + 1];
                      end if;
                    end if;
                  end if;
                end if;

                //The heat rate of the metal wall follows the hot side convention
                qdot_wall[1, j, i] = kA_tot*(t_av_hot[j, i] - t_av_cold[k1[j], k2[i]]);
                qdot_hot[j,i] = qdot_wall[1, j, i] + qdot_wall[2, j, i];

                //First row hot side and last row cold side have only one heat flux
                if (j == 1) then
                  k2[i] = N_cell_pc + 1 - i;
                  qdot_wall[1, j, i] = qdot_cold[k1[j], N_cell_pc + 1 - i];
                  qdot_wall[2, j, i] = kA_tot*(t_av_hot[j, i] - t_av_cold[k1[j] - 1,
                  k2[i]]);
                elseif (j == N_passes_d) then
                  qdot_wall[2, j, i] = 0;
                  qdot_wall[1, j, i] + qdot_wall[2, j - 1, i] = qdot_cold[k1[j],
                    k2[i]];
                else
                  qdot_wall[2, j, i] = kA_tot*(t_av_hot[j, i] - t_av_cold[k1[j] - 1,
                  k2[i]]);
                  qdot_wall[1, j, i] + qdot_wall[2, j - 1, i] = qdot_cold[k1[j],
                    k2[i]];
                end if;
              end for;
            end for;

            qdot_tot   = C_hot*(t_hot_in - t_hot_out);
            qdot_tot   = C_cold*(t_cold_out - t_cold_in);
            qdot_tot   = kA_tot*N_plates_d*N_cell_pc*DTML_tilde;
        DTML = Design.Miscellanea.log_mean_delta_T(
                t_hot_in,
                t_hot_out,
                t_cold_in,
                t_cold_out);

            if (layout == 1) then
              t_hot_out  = sum(t_hot[:, N_cell_pc + 1])/N_passes_d;
            else
              t_hot_out  = t_hot[N_passes_d, N_cell_pc + 1];
            end if;
      //       t_cold_out = sum(t_cold[:, N_cell_pc + 1])/N_passes_d;

        annotation (experiment(Tolerance=1e-006), __Dymola_experimentSetupOutput);
      end cell_method_flat_plate;
    end Cell_method;
  end Project;

  package Merge "Merge cycle and component design"

    model ORCHID "Design of the ORCHID test rig. No design of the components."
      parameter Medium.SaturationProperties  sat =  Medium.setSat_p(0.33e5)
        "Saturation properties";
      package Hot = Media.CoolProp.THERM66 "Medium model hot cells";
      package Medium = Media.CoolProp.Toluene "Medium model hot cells";
      package Coolant = Media.CoolProp.Water_TTSE "Medium model hot cells";

      CycleTempo.Components.Turbomachinery.Pump Pump(
        redeclare package Medium = Medium,
        eta_is=0.9,
        eta_m=1) "Centrifugal pump"
        annotation (Placement(transformation(extent={{13,14},{-13,-14}},
            rotation=-90,
            origin={77,-48})));
      CycleTempo.Components.HEX.Evaporator evaporator(
        dp_h=0,
        redeclare package Medium_c = Medium,
        redeclare package Medium_h = Hot,
        dp_c=0,
        dT_int=30,
        use_dT_in=false,
        use_dT_int=false,
        hh_in_start=1e5)
                        annotation (Placement(transformation(
            extent={{-18,-15.5},{18,15.5}},
            rotation=0,
            origin={0,59.5})));
      CycleTempo.Components.Turbomachinery.Turbine turbine(
        redeclare package Medium = Medium,
        eta_is=0.8,
        eta_m=1) "Axial turbine"
        annotation (Placement(transformation(extent={{53,17},{92,50}})));
      CycleTempo.Components.HEX.Condenser condenser(
        dp_h=0,
        redeclare package Medium_h = Medium,
        redeclare package Medium_c = Media.CoolProp.Water,
        dp_c=0,
        dT_int=20,
        use_dT_int=true) "condenser model"  annotation (Placement(transformation(
            extent={{-12,-11},{12,11}},
            rotation=0,
            origin={0,-72})));
      CycleTempo.Components.Flags.CC CC(redeclare package Medium = Medium)
        annotation (Placement(transformation(
            extent={{10,-10},{-10,10}},
            rotation=-90,
            origin={77,-78})));
      CycleTempo.Components.HEX.Heat_exchanger recuperator(
        redeclare package Medium_h = Medium,
        redeclare package Medium_c = Medium,
        dp_h=0,
        dp_c=0,
        dT_out=20,
        use_dT_out=false)
        annotation (Placement(transformation(extent={{11,-32},{-11,-12}})));
      CycleTempo.Components.Flags.ADDCO N3(redeclare package Medium = Medium)
        annotation (Placement(transformation(extent={{-170,4},{-150,26}})));
      CycleTempo.Components.Flags.ADDCO N20(
        redeclare package Medium = Coolant,
        use_T=true,
        use_p=true,
        T=278.15,
        p=100000)
        annotation (Placement(transformation(extent={{-97,-82},{-79,-63}})));
      CycleTempo.Components.Flags.ADDCO N21(
        redeclare package Medium = Coolant,
        T=313.15,
        use_T=false)
                  annotation (Placement(transformation(extent={{45,-81},{23,-63}})));
      CycleTempo.Components.Flags.ADDCO N2(
        redeclare package Medium = Medium,
        use_p=true,
        p=1500000) annotation (Placement(transformation(
            extent={{10,-10},{-10,10}},
            rotation=0,
            origin={87,-22})));
      CycleTempo.Components.Flags.ADDCO N6(
        redeclare package Medium = Medium,
        use_T=true,
        T=343.15) annotation (Placement(transformation(
            extent={{9,9},{-9,-9}},
            rotation=180,
            origin={-158,-41})));
      CycleTempo.Components.Flags.ADDCO N10(
        redeclare package Medium = Hot,
        m_flow=0.5,
        use_T=true,
        use_p=true,
        use_m_flow=true,
        T=618.15,
        p=200000)         annotation (Placement(transformation(
            extent={{-11,-10},{11,10}},
            rotation=0,
            origin={-11,96})));
      CycleTempo.Components.Flags.ADDCO N11(
        redeclare package Medium = Hot,
        visible=true,
        use_T=true,
        T=513.15) annotation (Placement(transformation(
            extent={{-10,-9},{10,9}},
            rotation=0,
            origin={-10,20})));
      CycleTempo.Components.Flags.ADDCO N4(
        redeclare package Medium = Medium,
        T=593.15,
        use_T=true)
                  annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=0,
            origin={44,96})));
      CycleTempo.Components.Flags.ADDCO N5(redeclare package Medium = Medium,
          visible=true) annotation (Placement(transformation(
            extent={{10.5,-9.5},{-10.5,9.5}},
            rotation=0,
            origin={134.5,3.5})));
      CycleTempo.Components.Flags.ADDCO N1(
        redeclare package Medium = Medium,
        T=323.15,
        p=33000,
        use_T=true,
        use_p=false)
                 annotation (Placement(transformation(
            extent={{-9,-9},{9,9}},
            rotation=0,
            origin={-88,-105})));
      CycleTempo.Components.Flags.START START1(
        redeclare package Medium = Medium,
        m_flow=0.15,
        p=1460000,
        T=453.15)
        annotation (Placement(transformation(extent={{-120,36},{-102,54}})));
      CycleTempo.Components.Electrics.Generator generator(eta_el=0.98)
        annotation (Placement(transformation(extent={{102,18},{133,49}})));
      CycleTempo.Components.Flags.ADDCOW Pout(W=33e3, use_W=false)
        annotation (Placement(transformation(extent={{162,23},{142,44}})));
      CycleTempo.Components.Flags.ADDCOW Ppump(W=33e3, use_W=false)
        annotation (Placement(transformation(extent={{152,-58},{132,-37}})));
      CycleTempo.Components.Electrics.Motor motor(eta_el=0.9)
        annotation (Placement(transformation(extent={{97,-60},{122,-36}})));
      CycleTempo.Components.Flags.START START2(   redeclare package Medium = Medium,
        m_flow=0.15,
        p=33000,
        T=343.15)
        annotation (Placement(transformation(extent={{26,-50},{8,-32}})));
      CycleTempo.Components.Flags.START START3(   redeclare package Medium = Medium,
        m_flow=0.15,
        p=33000,
        T=593.15)
        annotation (Placement(transformation(extent={{34,41},{52,59}})));

        parameter Medium.ThermodynamicState hot_in = Medium.setState_pT(0.33e5, 276.9 + 273.15)
        "Thermodynamic state at the inlet of the hot side";
        parameter Medium.ThermodynamicState hot_out = Medium.setState_pT(0.33e5, 71.4 + 273.15)
        "Thermodynamic state at the outlet of the hot side";
        parameter Medium.ThermodynamicState cold_in = Medium.setState_pT(14.6e5, 51.4 + 273.15)
        "Thermodynamic state at the inlet of the hot side";
        parameter Medium.ThermodynamicState cold_out = Medium.setState_pT(14.6e5, 222.9 + 273.15)
        "Thermodynamic state at the outlet of the hot side";
      CycleTempo.Components.Flags.START START4(   redeclare package Medium = Medium,
        m_flow=0.15,
        p=33000)
        annotation (Placement(transformation(extent={{20,19},{38,37}})));
      CycleTempo.Components.Flags.START START5(   redeclare package Medium = Medium,
        m_flow=0.15,
        p=33000,
        T=313.15)
        annotation (Placement(transformation(extent={{118,-78},{100,-60}})));
    equation
      connect(CC.node_in, N1.node) annotation (Line(
          points={{77,-83.8},{77,-105},{-78.91,-105}},
          color={0,0,0},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(N11.node, evaporator.node_h_out) annotation (Line(
          points={{0.1,20},{0.1,30},{0,30},{0,37.8}},
          color={0,0,0},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(N10.node, evaporator.node_h_in) annotation (Line(
          points={{0.11,96},{0.11,86},{0,86},{0,81.51}},
          color={0,0,0},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(N5.node, recuperator.node_h_in) annotation (Line(
          points={{123.895,3.5},{0,3.5},{0,-7.8}},
          color={0,0,0},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(condenser.node_h_out, N1.node) annotation (Line(
          points={{0,-87.4},{0,-105},{-78.91,-105}},
          color={0,0,0},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(CC.node_out, Pump.node_in) annotation (Line(
          points={{77,-72},{77,-60.74}},
          color={0,0,0},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(Pump.node_out, N2.node) annotation (Line(
          points={{77,-35},{77,-22},{76.9,-22}},
          color={0,0,0},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(recuperator.node_c_in, N2.node) annotation (Line(
          points={{15.4,-22},{76.9,-22}},
          color={0,0,0},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(recuperator.node_c_out, N3.node) annotation (Line(
          points={{-15.51,-22},{-149.9,-22},{-149.9,15}},
          color={0,0,0},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(N20.node, condenser.node_c_in) annotation (Line(
          points={{-78.91,-72.5},{-26,-72.5},{-26,-72},{-16.8,-72}},
          color={0,0,0},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(condenser.node_c_out, N21.node) annotation (Line(
          points={{16.92,-72},{22.89,-72}},
          color={0,0,0},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(N6.node, condenser.node_h_in) annotation (Line(
          points={{-148.91,-41},{0,-41},{0,-56.38}},
          color={0,0,0},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(N3.node, START1.node) annotation (Line(
          points={{-149.9,15},{-149.9,14.5},{-102,14.5},{-102,45}},
          color={0,0,0},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(START1.node, evaporator.node_c_in) annotation (Line(
          points={{-102,45},{-102,59.5},{-25.2,59.5}},
          color={0,0,0},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(evaporator.node_c_out, turbine.node_in) annotation (Line(
          points={{25.38,59.5},{64.7,59.5},{64.7,50}},
          color={0,0,0},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(turbine.node_out, N5.node) annotation (Line(
          points={{75.23,13.7},{75.23,10.75},{123.895,10.75},{123.895,3.5}},
          color={0,0,0},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(turbine.terminal, generator.terminal_in) annotation (Line(
          points={{91.9025,33.5},{102.31,33.5}},
          color={0,0,0},
          smooth=Smooth.None));
      connect(generator.terminal_out, Pout.terminal) annotation (Line(
          points={{133,33.5},{142,33.5}},
          color={0,0,0},
          smooth=Smooth.None));
      connect(Ppump.terminal, motor.terminal_out) annotation (Line(
          points={{132,-47.5},{122,-47.5},{122,-48}},
          color={0,0,0},
          smooth=Smooth.None));
      connect(Pump.terminal, motor.terminal_in) annotation (Line(
          points={{91,-48.065},{95.5,-48.065},{95.5,-48},{97.25,-48}},
          color={0,0,0},
          smooth=Smooth.None));
      connect(N4.node, turbine.node_in) annotation (Line(
          points={{54.1,96},{65,96},{65,60},{64.7,60},{64.7,50}},
          color={0,0,0},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(recuperator.node_h_out, condenser.node_h_in) annotation (Line(
          points={{0,-36},{0,-56.38}},
          color={0,0,0},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(START2.node, condenser.node_h_in) annotation (Line(
          points={{8,-41},{0,-41},{0,-56.38}},
          color={0,0,0},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(START3.node, turbine.node_in) annotation (Line(
          points={{52,50},{64.7,50}},
          color={0,0,0},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(START4.node, turbine.node_out) annotation (Line(
          points={{38,28},{46,28},{46,16},{75.23,16},{75.23,13.7}},
          color={0,0,0},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(Pump.node_in, START5.node) annotation (Line(
          points={{77,-60.74},{90,-60.74},{90,-69},{100,-69}},
          color={0,0,0},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-180,
                -120},{180,120}}),      graphics), Icon(coordinateSystem(extent={{-180,
                -120},{180,120}})));
    end ORCHID;

    model ORCHID_design
      "Design of the ORCHID test rig. Design of the components."

      package Hot = Media.CoolProp.THERM66;
      package Medium = Media.CoolProp.Toluene_TTSE;
      package Coolant = Media.CoolProp.Water_TTSE;

      parameter Hot.ThermodynamicState hot_in = Hot.setState_pT(2e5, 345 + 273.15)
        "Thermodynamic state at the inlet of the hot side";
      parameter Hot.ThermodynamicState hot_out = Hot.setState_pT(2e5, 240 + 273.15)
        "Thermodynamic state at the outlet of the hot side";
      parameter Medium.SaturationProperties sat = Medium.setSat_T(50 + 273.15);

        Design.Components.HEX.Flat_plate evaporator(
        redeclare package Medium_hot = Hot,
        redeclare package Medium_cold = Medium,
        redeclare function cost = Design.Miscellanea.Cost.FP_Rafferty,
        thick=0.8e-3,
        b=2.3e-3,
        ht_hot_f1=1e9,
        ht_cold_f1=1e9,
        X=1,
        N_ch_p=14,
        w=0.1325,
        d_pt=32e-3,
        redeclare Design.Heat_transfer.Plates.evaporation_Martin ht_cold(
          p_in=evaporator.node_c_in.p,
          At=evaporator.At,
          beta=evaporator.beta,
          qdot=evaporator.plate.qdot_cold,
          qdot_tilde_start=evaporator.mdot_hot_start*(evaporator.h_hot_in_start - evaporator.h_hot_out_start)/(2*evaporator.l_start*evaporator.w*evaporator.N_ch_p)),
        redeclare Design.Materials.SS_AISI_304 material,
        redeclare Design.Pressure_drops.Plates.evaporation_Martin dp_cold(p_in=evaporator.node_c_in.p,
            beta=evaporator.beta),
        redeclare Design.Miscellanea.check_velocity check_hot(T_sat=500),
        redeclare Design.Miscellanea.check_velocity check_cold(umin=max(evaporator.ht_cold.u)),
        N_cell_pc=8,
        redeclare Design.Miscellanea.topology_PHE.parallel tpg_hot,
        redeclare Design.Miscellanea.topology_PHE.parallel tpg_cold,
        h_hot_in_start=hot_in.h,
        h_hot_out_start=hot_out.h,
        l_od=0.6,
        offdesign=false,
        l_start=1,
        use_dp=true,
        beta=1.0471975511966,
        t_hot_in_start=618.15,
        t_hot_out_start=513.15,
        t_cold_in_start=439.15,
        t_cold_out_start=593.15,
        p_hot_start=200000,
        p_cold_start=1500000,
        mdot_hot_start=0.5,
        mdot_cold_start=0.22)
        "Model of a flat plate heat exchanger. Therminol 66 evaporates the working fluid."
        annotation (Placement(transformation(extent={{-84,-86},{12,-20}})));

      Design.Components.Flags.ADDCO N_11(
        redeclare package Medium = Hot,
        T=513.15,
        use_T=true)
        annotation (Placement(transformation(extent={{-145,-28},{-112,2}})));
      Design.Components.Flags.ADDCO N_3(
        redeclare package Medium = Medium,
        m_flow=0.22,
        use_p=false,
        use_T=false,
        T=440.15,
        p=1500000,
        use_m_flow=false)
        annotation (Placement(transformation(extent={{-159,-97},{-130,-69}})));
      Design.Components.Flags.ADDCO N_4(
        redeclare package Medium = Medium,
        use_T=true,
        T=593.15)
        annotation (Placement(transformation(extent={{51,-104},{80,-76}})));
      Design.Components.Flags.ADDCO N_10(
        use_T=true,
        use_p=true,
        use_m_flow=true,
        redeclare package Medium = Hot,
        m_flow=0.5,
        T=618.15,
        p=200000) annotation (Placement(transformation(extent={{55,-8},{84,20}})));
      CycleTempo.Components.Turbomachinery.Turbine turbine(
        redeclare package Medium = Medium,
        eta_is=0.8,
        eta_m=1) "Axial turbine"
        annotation (Placement(transformation(extent={{143,-163},{212,-106}})));
      CycleTempo.Components.Electrics.Generator generator(eta_el=0.98)
        annotation (Placement(transformation(extent={{234,-161},{285,-108}})));
      CycleTempo.Components.Flags.ADDCOW Pout(W=33e3, use_W=false)
        annotation (Placement(transformation(extent={{346,-156},{302,-114}})));
      CycleTempo.Components.Flags.ADDCO N5(redeclare package Medium = Medium,
          visible=true,
        use_p=false,
        p=10000)        annotation (Placement(transformation(
            extent={{16.5,-16},{-16.5,16}},
            rotation=0,
            origin={274.5,-202})));
        Design.Components.HEX.Flat_plate recuperator(
        redeclare package Medium_hot = Medium,
        redeclare package Medium_cold = Medium,
        redeclare function cost = Design.Miscellanea.Cost.FP_Rafferty,
        w=325e-3,
        d_pt=100e-3,
        ht_hot_f1=1e9,
        ht_cold_f1=1e9,
        X=1,
        redeclare Design.Miscellanea.topology_PHE.two_pass_two_pass tpg_hot(N_ch_p=
              recuperator.N_ch_p, stype=1),
        redeclare Design.Miscellanea.topology_PHE.two_pass_two_pass tpg_cold(N_ch_p=
             recuperator.N_ch_p, stype=2),
        redeclare Design.Materials.SS_AISI_304 material,
        redeclare Design.Miscellanea.check_velocity check_hot(T_sat=500),
        redeclare Design.Miscellanea.check_velocity check_cold(umin=max(evaporator.ht_cold.u)),
        l_od=0.6,
        offdesign=false,
        l_start=0.6,
        N_cell_pc=1,
        thick=0.9,
        b=3.2e-3,
        N_ch_p=36,
        use_dp=true,
        beta=1.0471975511966,
        t_hot_in_start=493.15,
        t_hot_out_start=343.15,
        t_cold_in_start=333.15,
        t_cold_out_start=439.15,
        p_hot_start=10000,
        p_cold_start=1500000,
        mdot_hot_start=0.22,
        mdot_cold_start=0.22)
        "Model of a flat plate heat exchanger. Recuperator of the ORC unit."
        annotation (Placement(transformation(extent={{53,-239},{149,-173}})));
      CycleTempo.Components.Flags.ADDCO N6(
        redeclare package Medium = Medium,
        use_T=true,
        T=343.15,
        p=10000,
        use_p=false)
                  annotation (Placement(transformation(
            extent={{15,14},{-15,-14}},
            rotation=180,
            origin={23,-152})));
      CycleTempo.Components.Flags.ADDCO N2(
        redeclare package Medium = Medium,
        p=1500000,
        use_p=true)
                   annotation (Placement(transformation(
            extent={{15.5,-15},{-15.5,15}},
            rotation=0,
            origin={114.5,-265})));
        Design.Components.HEX.Flat_plate condenser(
        redeclare package Medium_hot = Medium,
        redeclare package Medium_cold = Coolant,
        redeclare function cost = Design.Miscellanea.Cost.FP_Rafferty,
        ht_hot_f1=1e9,
        ht_cold_f1=1e9,
        w=325e-3,
        d_pt=100e-3,
        redeclare Design.Materials.SS_AISI_304 material,
        thick=0.9e-3,
        b=3.12e-3,
        redeclare Design.Miscellanea.check_velocity check_hot(T_sat=50 + 273.15),
        redeclare Design.Miscellanea.check_velocity check_cold(umin=min(condenser.ht_cold.u)),
        l_start=0.5,
        offdesign=false,
        redeclare Design.Heat_transfer.Plates.condensation_Longo ht_hot(beta=condenser.beta,
            p_in=condenser.node_h_in.p),
        redeclare Design.Pressure_drops.Plates.evaporation_Martin dp_hot(p_in=condenser.node_h_in.p,
            beta=condenser.beta),
        N_ch_p=20,
        X=1,
        N_cell_pc=7,
        use_dp=true,
        beta=1.0471975511966,
        t_hot_in_start=343.15,
        t_hot_out_start=313.15,
        t_cold_in_start=278.15,
        t_cold_out_start=303.15,
        p_hot_start=10000,
        p_cold_start=100000,
        mdot_hot_start=0.22,
        mdot_cold_start=0.8)
        "Model of a flat plate condenser. Water condenses Toluene."
        annotation (Placement(transformation(extent={{-118,-239},{-22,-173}})));

      Design.Components.Flags.ADDCO COLD_IN(
        use_T=true,
        use_p=true,
        redeclare package Medium = Coolant,
        use_m_flow=true,
        T=278.15,
        p=100000,
        m_flow=0.8)
        annotation (Placement(transformation(extent={{-193,-294},{-164,-266}})));
      Design.Components.Flags.ADDCO COLD_OUT(
        redeclare package Medium = Coolant,
        T=598.15,
        use_T=false)
        annotation (Placement(transformation(extent={{7,-272},{-22,-244}})));
      Design.Components.Flags.ADDCO HOT_OUT(
        redeclare package Medium = Medium,
        h=sat.hl,
        use_h=true,
        T=313.15,
        p=sat.psat,
        use_T=false,
        use_p=true)
        annotation (Placement(transformation(extent={{-205,-186},{-172,-156}})));
      CycleTempo.Components.Flags.CC CC(redeclare package Medium = Medium)
        annotation (Placement(transformation(
            extent={{21,-22},{-21,22}},
            rotation=180,
            origin={-178,-361})));
      CycleTempo.Components.Turbomachinery.Pump Pump(
        redeclare package Medium = Medium,
        eta_is=0.9,
        eta_m=1) "Centrifugal pump"
        annotation (Placement(transformation(extent={{-21,22},{21,-22}},
            rotation=0,
            origin={-67,-361})));
      CycleTempo.Components.Flags.START sTART(redeclare package Medium = Medium,
        m_flow=0.22,
        p=1500000,
        T=439.15)
        annotation (Placement(transformation(extent={{-242,-58},{-202,-22}})));
      CycleTempo.Components.Flags.ADDCOW Ppump(W=33e3, use_W=false)
        annotation (Placement(transformation(extent={{70,-358},{32,-321}})));
      CycleTempo.Components.Electrics.Motor motor(eta_el=0.9)
        annotation (Placement(transformation(extent={{-31,-356},{3,-322}})));
    equation

      connect(N_11.node, evaporator.node_h_out) annotation (Line(
          points={{-111.835,-13},{-111.835,-72.8},{-84,-72.8}},
          color={0,0,0},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(N_3.node, evaporator.node_c_in) annotation (Line(
          points={{-129.855,-83},{-104,-83},{-104,-83.03},{-84,-83.03}},
          color={0,0,0},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(evaporator.node_c_out, N_4.node) annotation (Line(
          points={{12,-36.665},{130,-36.665},{130,-90},{80.145,-90}},
          color={0,0,0},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(N_10.node, evaporator.node_h_in) annotation (Line(
          points={{84.145,6},{130,6},{130,-26.6},{12,-26.6}},
          color={0,0,0},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(turbine.node_in, N_4.node) annotation (Line(
          points={{163.7,-106},{80.145,-106},{80.145,-90}},
          color={0,0,0},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(turbine.terminal, generator.terminal_in) annotation (Line(
          points={{211.827,-134.5},{234.51,-134.5}},
          color={0,0,0},
          smooth=Smooth.None));

      connect(generator.terminal_out, Pout.terminal) annotation (Line(
          points={{285,-134.5},{294.5,-134.5},{294.5,-135},{302,-135}},
          color={0,0,0},
          smooth=Smooth.None));
      connect(turbine.node_out, N5.node) annotation (Line(
          points={{182.33,-168.7},{218,-168.7},{218,-202},{257.835,-202}},
          color={0,0,0},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(recuperator.node_h_in, turbine.node_out) annotation (Line(
          points={{149,-179.6},{182.33,-179.6},{182.33,-168.7}},
          color={0,0,0},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(N6.node, recuperator.node_h_out) annotation (Line(
          points={{38.15,-152},{53,-152},{53,-225.8}},
          color={0,0,0},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(condenser.node_h_in, recuperator.node_h_out) annotation (Line(
          points={{-22,-179.6},{18,-179.6},{18,-225.8},{53,-225.8}},
          color={0,0,0},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(recuperator.node_c_out, evaporator.node_c_in) annotation (Line(
          points={{149,-189.665},{200,-189.665},{200,-380},{-278,-380},{-278,
              -120},{-84,-120},{-84,-83.03}},
          color={0,0,0},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(HOT_OUT.node, condenser.node_h_out) annotation (Line(
          points={{-171.835,-171},{-118,-171},{-118,-225.8}},
          color={0,0,0},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(condenser.node_h_out, CC.node_in) annotation (Line(
          points={{-118,-225.8},{-260,-225.8},{-260,-361},{-190.18,-361}},
          color={0,0,0},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(CC.node_out, Pump.node_in) annotation (Line(
          points={{-165.4,-361},{-87.58,-361}},
          color={0,0,0},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(COLD_IN.node, condenser.node_c_in) annotation (Line(
          points={{-163.855,-280},{-118,-280},{-118,-236.03}},
          color={0,0,0},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(N2.node, recuperator.node_c_in) annotation (Line(
          points={{98.845,-265},{98.845,-264},{53,-264},{53,-236.03}},
          color={0,0,0},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(sTART.node, N_3.node) annotation (Line(
          points={{-202,-40},{-129.855,-40},{-129.855,-83}},
          color={0,0,0},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(motor.terminal_out, Ppump.terminal) annotation (Line(
          points={{3,-339},{3,-339.5},{32,-339.5}},
          color={0,0,0},
          smooth=Smooth.None));
      connect(condenser.node_c_out, COLD_OUT.node) annotation (Line(
          points={{-22,-189.665},{-22,-258},{-22.145,-258}},
          color={0,0,0},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(Pump.terminal, motor.terminal_in) annotation (Line(
          points={{-67.105,-339},{-48.5525,-339},{-48.5525,-339},{-30.66,-339}},
          color={0,0,0},
          smooth=Smooth.None));

      connect(Pump.node_out, N2.node) annotation (Line(
          points={{-46,-361},{98.845,-361},{98.845,-265}},
          color={0,0,0},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-280,
                -380},{320,40}}),       graphics), Icon(coordinateSystem(extent={{-280,
                -380},{320,40}})),
        experiment(__Dymola_NumberOfIntervals=1, __Dymola_Algorithm="Dassl"),
        __Dymola_experimentSetupOutput);
    end ORCHID_design;

    model ORCHID_od "Design of the ORCHID test rig. Part-load."

      package Hot = Media.CoolProp.THERM66;
      package Medium = Media.CoolProp.Toluene_TTSE;
      package Coolant = Media.CoolProp.Water_TTSE;

      parameter Hot.ThermodynamicState hot_in = Hot.setState_pT(2e5, 345 + 273.15)
        "Thermodynamic state at the inlet of the hot side";
      parameter Hot.ThermodynamicState hot_out = Hot.setState_pT(2e5, 240 + 273.15)
        "Thermodynamic state at the outlet of the hot side";
      parameter Medium.SaturationProperties sat = Medium.setSat_T(50 + 273.15);
      parameter Modelica.SIunits.MassFlowRate mdot_Medium = 0.15
        "Start value organic medium";
      parameter Modelica.SIunits.MassFlowRate mdot_Hot = 0.2
        "Start value organic hot fluid";
      parameter Modelica.SIunits.MassFlowRate mdot_Coolant = 0.5
        "Start value organic coolant";

        Design.Components.HEX.Flat_plate evaporator(
        redeclare package Medium_hot = Hot,
        redeclare package Medium_cold = Medium,
        redeclare function cost = Design.Miscellanea.Cost.FP_Rafferty,
        thick=0.8e-3,
        b=2.3e-3,
        ht_hot_f1=1e9,
        ht_cold_f1=1e9,
        X=1,
        N_ch_p=14,
        w=0.1325,
        d_pt=32e-3,
        redeclare Design.Heat_transfer.Plates.evaporation_Martin ht_cold(
          p_in=evaporator.node_c_in.p,
          At=evaporator.At,
          beta=evaporator.beta,
          qdot=evaporator.plate.qdot_cold,
          qdot_tilde_start=evaporator.mdot_hot_start*(evaporator.h_hot_in_start - evaporator.h_hot_out_start)/(2*evaporator.l_start*evaporator.w*evaporator.N_ch_p)),
        redeclare Design.Materials.SS_AISI_304 material,
        redeclare Design.Pressure_drops.Plates.evaporation_Martin dp_cold(p_in=evaporator.node_c_in.p,
            beta=evaporator.beta),
        redeclare Design.Miscellanea.check_velocity check_hot(T_sat=500),
        redeclare Design.Miscellanea.check_velocity check_cold(umin=max(evaporator.ht_cold.u)),
        N_cell_pc=8,
        redeclare Design.Miscellanea.topology_PHE.parallel tpg_hot,
        redeclare Design.Miscellanea.topology_PHE.parallel tpg_cold,
        mdot_hot_start=mdot_Hot,
        h_hot_in_start=hot_in.h,
        h_hot_out_start=hot_out.h,
        mdot_cold_start=mdot_Medium,
        offdesign=true,
        l_od=1.16127,
        use_dp=true,
        l_start=1.16127,
        beta=1.0471975511966,
        t_hot_in_start=618.15,
        t_hot_out_start=513.15,
        t_cold_in_start=439.15,
        t_cold_out_start=593.15,
        p_hot_start=200000,
        p_cold_start=1500000,
        dp_hot_start=5656.17,
        dp_cold_start=7119.53)
        "Model of a flat plate heat exchanger. Therminol 66 evaporates the working fluid."
        annotation (Placement(transformation(extent={{-84,-86},{12,-20}})));

      Design.Components.Flags.ADDCO N_11(
        redeclare package Medium = Hot,
        T=513.15,
        use_T=false)
        annotation (Placement(transformation(extent={{-145,-28},{-112,2}})));
      Design.Components.Flags.ADDCO N_3(
        redeclare package Medium = Medium,
        m_flow=mdot_Medium,
        use_p=false,
        use_T=false,
        T=440.15,
        p=1500000,
        use_m_flow=false)
        annotation (Placement(transformation(extent={{-159,-97},{-130,-69}})));
      Design.Components.Flags.ADDCO N_4(
        redeclare package Medium = Medium,
        T=593.15,
        use_T=true)
        annotation (Placement(transformation(extent={{51,-104},{80,-76}})));
      Design.Components.Flags.ADDCO N_10(
        use_T=true,
        use_p=true,
        use_m_flow=true,
        redeclare package Medium = Hot,
        T=618.15,
        p=200000,
        m_flow=0.21)
                  annotation (Placement(transformation(extent={{55,-8},{84,20}})));
      CycleTempo.Components.Turbomachinery.Turbine turbine(
        redeclare package Medium = Medium,
        eta_is=0.8,
        eta_m=1,
        off_design(CT_d=3.94737e-06),
        activated=true) "Axial turbine"
        annotation (Placement(transformation(extent={{143,-163},{212,-106}})));
      CycleTempo.Components.Electrics.Generator generator(eta_el=0.98)
        annotation (Placement(transformation(extent={{234,-161},{285,-108}})));
      CycleTempo.Components.Flags.ADDCOW Pout(W=33e3, use_W=false)
        annotation (Placement(transformation(extent={{346,-155},{302,-113}})));
      CycleTempo.Components.Flags.ADDCO N5(redeclare package Medium = Medium,
          visible=true,
        use_p=false,
        p=10000)        annotation (Placement(transformation(
            extent={{16.5,-16},{-16.5,16}},
            rotation=0,
            origin={274.5,-202})));
        Design.Components.HEX.Flat_plate recuperator(
        redeclare package Medium_hot = Medium,
        redeclare package Medium_cold = Medium,
        redeclare function cost = Design.Miscellanea.Cost.FP_Rafferty,
        w=325e-3,
        d_pt=100e-3,
        ht_hot_f1=1e9,
        ht_cold_f1=1e9,
        X=1,
        redeclare Design.Miscellanea.topology_PHE.two_pass_two_pass tpg_hot(N_ch_p=
              recuperator.N_ch_p, stype=1),
        redeclare Design.Miscellanea.topology_PHE.two_pass_two_pass tpg_cold(N_ch_p=
             recuperator.N_ch_p, stype=2),
        redeclare Design.Materials.SS_AISI_304 material,
        redeclare Design.Miscellanea.check_velocity check_hot(T_sat=500),
        redeclare Design.Miscellanea.check_velocity check_cold(umin=max(evaporator.ht_cold.u)),
        mdot_cold_start=mdot_Medium,
        mdot_hot_start=mdot_Medium,
        N_cell_pc=1,
        thick=0.9,
        b=3.2e-3,
        N_ch_p=36,
        offdesign=true,
        l_od=2.21,
        use_dp=true,
        l_start=2.21,
        beta=1.0471975511966,
        t_hot_in_start=453.15,
        t_hot_out_start=343.15,
        t_cold_in_start=318.15,
        t_cold_out_start=439.15,
        p_hot_start=10000,
        p_cold_start=1500000,
        dp_hot_start=29144.4,
        dp_cold_start=449.385)
        "Model of a flat plate heat exchanger. Recuperator of the ORC unit."
        annotation (Placement(transformation(extent={{53,-239},{149,-173}})));
      CycleTempo.Components.Flags.ADDCO N6(
        redeclare package Medium = Medium,
        use_p=false,
        T=343.15,
        p=10000,
        use_T=false)
                  annotation (Placement(transformation(
            extent={{15,14},{-15,-14}},
            rotation=180,
            origin={23,-152})));
      CycleTempo.Components.Flags.ADDCO N2(
        redeclare package Medium = Medium,
        p=1500000,
        use_p=false)
                   annotation (Placement(transformation(
            extent={{15.5,-15},{-15.5,15}},
            rotation=0,
            origin={114.5,-265})));
        Design.Components.HEX.Flat_plate condenser(
        redeclare package Medium_hot = Medium,
        redeclare package Medium_cold = Coolant,
        redeclare function cost = Design.Miscellanea.Cost.FP_Rafferty,
        ht_hot_f1=1e9,
        ht_cold_f1=1e9,
        w=325e-3,
        d_pt=100e-3,
        redeclare Design.Materials.SS_AISI_304 material,
        thick=0.9e-3,
        b=3.12e-3,
        redeclare Design.Miscellanea.check_velocity check_hot(T_sat=50 + 273.15),
        redeclare Design.Miscellanea.check_velocity check_cold(umin=min(condenser.ht_cold.u)),
        redeclare Design.Heat_transfer.Plates.condensation_Longo ht_hot(beta=condenser.beta,
            p_in=condenser.node_h_in.p),
        mdot_hot_start=mdot_Medium,
        redeclare Design.Pressure_drops.Plates.evaporation_Martin dp_hot(p_in=condenser.node_h_in.p,
            beta=condenser.beta),
        N_ch_p=20,
        X=1,
        offdesign=true,
        l_od=0.31601,
        mdot_cold_start=mdot_Coolant,
        l_start=0.31601,
        p_hot_start=sat.psat,
        N_cell_pc=10,
        use_dp=true,
        beta=1.0471975511966,
        t_hot_in_start=343.15,
        t_hot_out_start=308.15,
        t_cold_in_start=278.15,
        t_cold_out_start=303.15,
        p_cold_start=100000,
        dp_hot_start=3619.37,
        dp_cold_start=162.125)
        "Model of a flat plate condenser. Water condenses Toluene."
        annotation (Placement(transformation(extent={{-120,-239},{-24,-173}})));

      Design.Components.Flags.ADDCO COLD_IN(
        use_T=true,
        use_p=true,
        redeclare package Medium = Coolant,
        m_flow=mdot_Coolant,
        use_m_flow=false,
        T=278.15,
        p=100000)
        annotation (Placement(transformation(extent={{-193,-294},{-164,-266}})));
      Design.Components.Flags.ADDCO COLD_OUT(
        redeclare package Medium = Coolant,
        T=598.15,
        use_T=false)
        annotation (Placement(transformation(extent={{7,-272},{-22,-244}})));
      Design.Components.Flags.ADDCO HOT_OUT(
        redeclare package Medium = Medium,
        h=sat.hl,
        use_h=true,
        p=sat.psat,
        use_T=false,
        use_p=true,
        T=313.15)
        annotation (Placement(transformation(extent={{-205,-186},{-172,-156}})));
      CycleTempo.Components.Flags.CC CC(redeclare package Medium = Medium)
        annotation (Placement(transformation(
            extent={{21,-22},{-21,22}},
            rotation=180,
            origin={-178,-361})));
      CycleTempo.Components.Turbomachinery.Pump Pump(
        redeclare package Medium = Medium,
        eta_is=0.9,
        eta_m=1,
        head_d=200.861,
        n_d=3000,
        V_flow_d=0.000288295,
        constant_n=false,
        activated=true) "Centrifugal pump"
        annotation (Placement(transformation(extent={{-21,22},{21,-22}},
            rotation=0,
            origin={-67,-361})));
      CycleTempo.Components.Flags.START sTART(redeclare package Medium = Medium,
        m_flow=mdot_Medium,
        p=1500000,
        T=453.15)
        annotation (Placement(transformation(extent={{-242,-58},{-202,-22}})));
      CycleTempo.Components.Flags.ADDCOW Ppump(W=33e3, use_W=false)
        annotation (Placement(transformation(extent={{70,-358},{32,-321}})));
      CycleTempo.Components.Electrics.Motor motor(eta_el=0.9)
        annotation (Placement(transformation(extent={{-31,-356},{3,-322}})));
      CycleTempo.Components.Flags.START sTART1(
                                              redeclare package Medium = Medium,
        m_flow=mdot_Medium,
        p=1500000,
        T=333.15)
        annotation (Placement(transformation(extent={{162,-360},{122,-324}})));
      CycleTempo.Components.Flags.START sTART4(
                                              redeclare package Medium = Coolant,
        m_flow=mdot_Coolant,
        p=80000,
        T=303.15)
        annotation (Placement(transformation(extent={{22,-318},{-18,-282}})));
    equation

      connect(N_11.node, evaporator.node_h_out) annotation (Line(
          points={{-111.835,-13},{-111.835,-72.8},{-84,-72.8}},
          color={0,0,0},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(N_3.node, evaporator.node_c_in) annotation (Line(
          points={{-129.855,-83},{-104,-83},{-104,-83.03},{-84,-83.03}},
          color={0,0,0},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(evaporator.node_c_out, N_4.node) annotation (Line(
          points={{12,-36.665},{130,-36.665},{130,-90},{80.145,-90}},
          color={0,0,0},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(N_10.node, evaporator.node_h_in) annotation (Line(
          points={{84.145,6},{130,6},{130,-26.6},{12,-26.6}},
          color={0,0,0},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(turbine.node_in, N_4.node) annotation (Line(
          points={{163.7,-106},{80.145,-106},{80.145,-90}},
          color={0,0,0},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(turbine.terminal, generator.terminal_in) annotation (Line(
          points={{211.827,-134.5},{234.51,-134.5}},
          color={0,0,0},
          smooth=Smooth.None));

      connect(generator.terminal_out, Pout.terminal) annotation (Line(
          points={{285,-134.5},{294.5,-134.5},{294.5,-134},{302,-134}},
          color={0,0,0},
          smooth=Smooth.None));
      connect(turbine.node_out, N5.node) annotation (Line(
          points={{182.33,-168.7},{218,-168.7},{218,-202},{257.835,-202}},
          color={0,0,0},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(recuperator.node_h_in, turbine.node_out) annotation (Line(
          points={{149,-179.6},{182.33,-179.6},{182.33,-168.7}},
          color={0,0,0},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(N6.node, recuperator.node_h_out) annotation (Line(
          points={{38.15,-152},{53,-152},{53,-225.8}},
          color={0,0,0},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(condenser.node_h_in, recuperator.node_h_out) annotation (Line(
          points={{-24,-179.6},{18,-179.6},{18,-225.8},{53,-225.8}},
          color={0,0,0},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(recuperator.node_c_out, evaporator.node_c_in) annotation (Line(
          points={{149,-189.665},{200,-189.665},{200,-380},{-278,-380},{-278,-120},{
              -84,-120},{-84,-83.03}},
          color={0,0,0},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(HOT_OUT.node, condenser.node_h_out) annotation (Line(
          points={{-171.835,-171},{-120,-171},{-120,-225.8}},
          color={0,0,0},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(condenser.node_h_out, CC.node_in) annotation (Line(
          points={{-120,-225.8},{-260,-225.8},{-260,-361},{-190.18,-361}},
          color={0,0,0},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(CC.node_out, Pump.node_in) annotation (Line(
          points={{-165.4,-361},{-87.58,-361}},
          color={0,0,0},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(COLD_IN.node, condenser.node_c_in) annotation (Line(
          points={{-163.855,-280},{-120,-280},{-120,-236.03}},
          color={0,0,0},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(N2.node, recuperator.node_c_in) annotation (Line(
          points={{98.845,-265},{98.845,-264},{53,-264},{53,-236.03}},
          color={0,0,0},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(sTART.node, N_3.node) annotation (Line(
          points={{-202,-40},{-129.855,-40},{-129.855,-83}},
          color={0,0,0},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(motor.terminal_out, Ppump.terminal) annotation (Line(
          points={{3,-339},{3,-339.5},{32,-339.5}},
          color={0,0,0},
          smooth=Smooth.None));
      connect(condenser.node_c_out, COLD_OUT.node) annotation (Line(
          points={{-24,-189.665},{-24,-258},{-22.145,-258}},
          color={0,0,0},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(Pump.terminal, motor.terminal_in) annotation (Line(
          points={{-67.105,-339},{-30.66,-339}},
          color={0,0,0},
          smooth=Smooth.None));
      connect(Pump.node_out, N2.node) annotation (Line(
          points={{-46,-361},{98.845,-361},{98.845,-265}},
          color={0,0,0},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(Pump.node_out, sTART1.node) annotation (Line(
          points={{-46,-361},{-46,-372},{122,-372},{122,-342}},
          color={0,0,0},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(sTART4.node, COLD_OUT.node) annotation (Line(
          points={{-18,-300},{-22,-300},{-22,-258},{-22.145,-258}},
          color={0,0,0},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-280,
                -380},{320,40}}),       graphics), Icon(coordinateSystem(extent={{-280,
                -380},{320,40}})),
        experiment(__Dymola_NumberOfIntervals=1, __Dymola_Algorithm="Dassl"),
        __Dymola_experimentSetupOutput);
    end ORCHID_od;
  end Merge;

  package Part_load "Tests with part-load models"
    model Pump "Test on the part-load model of the pump"

      CycleTempo.Components.Turbomachinery.Pump pump(
        redeclare package Medium = Test.Media.CoolProp.Water,
        eta_m=0.9,
        eta_is=0.7,
        head_d=131.082,
        activated=true,
        n_d=3000,
        V_flow_d=1e-3,
        constant_n=true)
        annotation (Placement(transformation(extent={{-28,-18},{8,16}})));
      CycleTempo.Components.Flags.ADDCO IN(
        use_T=true,
        use_p=true,
        use_m_flow=true,
        redeclare package Medium = ExternalMedia.Examples.WaterCoolProp,
        T=278.15,
        p=100000,
        m_flow=1)
        annotation (Placement(transformation(extent={{-70,-11},{-50,9}})));
      CycleTempo.Components.Flags.START sTART(
        m_flow=1,
        redeclare package Medium = Test.Media.CoolProp.Water,
        p=1000000,
        T=333.15)
        annotation (Placement(transformation(extent={{-98,24},{-62,56}})));
      CycleTempo.Components.Flags.ADDCO OUT(redeclare package Medium =
            ExternalMedia.Examples.WaterCoolProp)
        annotation (Placement(transformation(extent={{-82,68},{-62,88}})));
      CycleTempo.Components.Flags.ADDCOW Ppump
        annotation (Placement(transformation(extent={{-50,-50},{-30,-30}})));
    equation
      connect(IN.node, pump.node_in) annotation (Line(
          points={{-49.9,-1},{-27.64,-1}},
          color={0,0,0},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(sTART.node, pump.node_out) annotation (Line(
          points={{-62,40},{8,40},{8,-1}},
          color={0,0,0},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(OUT.node, sTART.node) annotation (Line(
          points={{-61.9,78},{-62,78},{-62,40}},
          color={0,0,0},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(Ppump.terminal, pump.terminal) annotation (Line(
          points={{-30,-40},{-10.09,-40},{-10.09,-18}},
          color={0,0,0},
          smooth=Smooth.None));
      annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
                -100},{100,100}}),      graphics));
    end Pump;

    model Turbine "Test on the part-load model of the turbine"

      CycleTempo.Components.Turbomachinery.Turbine turbine(
        redeclare package Medium = Test.Media.CoolProp.Cyclopentane,
        eta_m=0.9,
        eta_is=0.7,
        activated=true,
        redeclare
          CycleTempo.Components.Turbomachinery.Part_load.Turbine.de_Laval
          off_design(A_d=0.00486741, p_is_t_start=1000000))
        annotation (Placement(transformation(extent={{-28,-18},{8,16}})));
      CycleTempo.Components.Flags.ADDCO IN(
        use_T=true,
        use_m_flow=true,
        redeclare package Medium = Test.Media.CoolProp.Cyclopentane,
        use_p=false,
        T=513.15,
        p=3000000,
        m_flow=10)
        annotation (Placement(transformation(extent={{-68,31},{-48,51}})));
      CycleTempo.Components.Flags.ADDCO OUT(
        redeclare package Medium = Test.Media.CoolProp.Cyclopentane,
        use_p=true,
        p=100000)
        annotation (Placement(transformation(extent={{-68,-78},{-48,-58}})));
      CycleTempo.Components.Flags.ADDCOW Pt
        annotation (Placement(transformation(extent={{54,-11},{34,9}})));
      CycleTempo.Components.Flags.START START(
        redeclare package Medium = Test.Media.CoolProp.Cyclopentane,
        m_flow=20,
        p=3000000,
        T=513.15)
                 annotation (Placement(transformation(extent={{-92,60},{-56,92}})));
    equation
      connect(IN.node, turbine.node_in) annotation (Line(
          points={{-47.9,41},{-17,41},{-17,16},{-17.2,16}},
          color={0,0,0},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(Pt.terminal, turbine.terminal) annotation (Line(
          points={{34,-1},{7.91,-1}},
          color={0,0,0},
          smooth=Smooth.None));
      connect(OUT.node, turbine.node_out) annotation (Line(
          points={{-47.9,-68},{-7.48,-68},{-7.48,-21.4}},
          color={0,0,0},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(START.node, IN.node) annotation (Line(
          points={{-56,76},{-47.9,76},{-47.9,41}},
          color={0,0,0},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
                -100},{100,100}}),      graphics));
    end Turbine;
  end Part_load;
  annotation (uses(Modelica(version="3.2.1")));
end Test;
