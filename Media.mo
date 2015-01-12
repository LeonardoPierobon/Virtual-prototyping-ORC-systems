within VIP;
package Media "I am the packag containing the definition of the fluids"
  package OneRandomOrganicFluid "Change the name to something more appropriate"
      extends ExternalMedia.Media.FluidPropMedium(
      mediumName = "Name of the fluid for documentation purposes",
      libraryName = "FluidProp.RefProp",
      substanceNames = {"cyclopentane"},
      ThermoStates = Modelica.Media.Interfaces.PartialMedium.Choices.IndependentVariables.ph);
  end OneRandomOrganicFluid;

  package CoolProp "Media defined using CoolProp"
    package MM "Hexamethyldisiloxane "
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

    package Water_TTSE "CoolProp model of Water using tables"
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
  end CoolProp;

  package RefProp "Media defined using RefProp (via CoolProp)"
    package MM "Hexamethyldisiloxane "
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
        mediumName = "Cyclopentane",
        substanceNames = {"REFPROP-Cyclopentane"},
        ThermoStates = Modelica.Media.Interfaces.PartialMedium.Choices.IndependentVariables.ph);
    end Cyclopentane;

    package Water "CoolProp model of Water"
      extends ExternalMedia.Media.CoolPropMedium(
        mediumName = "Water",
        substanceNames = {"REFPROP-Water"},
        ThermoStates = Modelica.Media.Interfaces.PartialMedium.Choices.IndependentVariables.ph);
    end Water;

    package Methanol "methanol"
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
        libraryName = "FluidProp.StanMix3",
        substanceNames = {"Cyclopentane"},
        ThermoStates = Modelica.Media.Interfaces.PartialMedium.Choices.IndependentVariables.ph);
    end Cyclopentane;

    package Water "CoolProp model of Water"
      extends ExternalMedia.Media.FluidPropMedium(
        mediumName = "Water",
        libraryName = "FluidProp.StanMix3",
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
