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
          reference_X={0},
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

    package Gas "Ideal gases package"
      package FlueGas "flue gas"
        extends Test.Media.Gas.MixtureGasNasa(
          mediumName="FlueGas",
          data={Modelica.Media.IdealGases.Common.SingleGasesData.O2,
                Modelica.Media.IdealGases.Common.SingleGasesData.Ar,
                Modelica.Media.IdealGases.Common.SingleGasesData.H2O,
                Modelica.Media.IdealGases.Common.SingleGasesData.CO2,
                Modelica.Media.IdealGases.Common.SingleGasesData.N2},
          fluidConstants={
                Modelica.Media.IdealGases.Common.FluidData.O2,
                Modelica.Media.IdealGases.Common.FluidData.Ar,
                Modelica.Media.IdealGases.Common.FluidData.H2O,
                Modelica.Media.IdealGases.Common.FluidData.CO2,
                Modelica.Media.IdealGases.Common.FluidData.N2},
          substanceNames={"Oxygen","Argon","Water","Carbondioxide","Nitrogen"},
          reference_X={0.23,0.02,0.01,0.04,0.7});

      end FlueGas;

    partial package MixtureGasNasa
        "Medium model of a mixture of ideal gases based on NASA source"

      import Modelica.Math;
      import Modelica.Media.Interfaces.Choices.ReferenceEnthalpy;

      extends Modelica.Media.Interfaces.PartialMixtureMedium(
         ThermoStates=Modelica.Media.Interfaces.Choices.IndependentVariables.pTX,
         substanceNames=data[:].name,
         reducedX = false,
         singleState=false,
         reference_X=fill(1/nX,nX),
         SpecificEnthalpy(start=if referenceChoice==ReferenceEnthalpy.ZeroAt0K then 3e5 else
            if referenceChoice==ReferenceEnthalpy.UserDefined then h_offset else 0, nominal=1.0e5),
         Density(start=10, nominal=10),
         AbsolutePressure(start=10e5, nominal=10e5),
         Temperature(min=200, max=6000, start=500, nominal=500));

        redeclare record ThermodynamicState "Thermodynamic state variables"
          AbsolutePressure p "pressure";
          Temperature T "temperature";
          SpecificEntropy s "specific entropy";
          MassFraction X[nX]
            "Mass fractions (= (component mass)/total mass  m_i/m)";
        //   VelocityOfSound a "velocity of sound";
        //   Modelica.SIunits.CubicExpansionCoefficient beta
        //     "isobaric expansion coefficient";
          SpecificHeatCapacity cp "specific heat capacity cp";
          SpecificHeatCapacity cv "specific heat capacity cv";
          Density d "density";
        //   DerDensityByEnthalpy ddhp
        //     "derivative of density wrt enthalpy at constant pressure";
        //   DerDensityByPressure ddph
        //     "derivative of density wrt pressure at constant enthalpy";
           DynamicViscosity eta "dynamic viscosity";
           SpecificEnthalpy h "specific enthalpy";
        //   Modelica.SIunits.Compressibility kappa "compressibility";
           ThermalConductivity lambda "thermal conductivity";

        end ThermodynamicState;

    //   redeclare record extends FluidConstants "Fluid constants"
    //   end FluidConstants;

      constant Modelica.Media.IdealGases.Common.DataRecord[:] data
          "Data records of ideal gas substances";
        // ={Common.SingleGasesData.N2,Common.SingleGasesData.O2}

      constant Boolean excludeEnthalpyOfFormation=true
          "If true, enthalpy of formation Hf is not included in specific enthalpy h";
      constant ReferenceEnthalpy referenceChoice=ReferenceEnthalpy.ZeroAt0K
          "Choice of reference enthalpy";
      constant SpecificEnthalpy h_offset=0.0
          "User defined offset for reference enthalpy, if referenceChoice = UserDefined";

    //   constant FluidConstants[nX] fluidConstants
    //     "Additional data needed for transport properties";
      constant MolarMass[nX] MMX=data[:].MM "Molar masses of components";
      constant Integer methodForThermalConductivity(min=1,max=2)=1;
      redeclare replaceable model extends BaseProperties(
        T(stateSelect=if preferredMediumStates then StateSelect.prefer else StateSelect.default),
        p(stateSelect=if preferredMediumStates then StateSelect.prefer else StateSelect.default),
        Xi(each stateSelect=if preferredMediumStates then StateSelect.prefer else StateSelect.default),
        final standardOrderComponents=true)
          "Base properties (p, d, T, h, u, R, MM, X, and Xi of NASA mixture gas"
      equation
        assert(T >= 200 and T <= 6000, "
Temperature T (="     + String(T) + " K = 200 K) is not in the allowed range
200 K <= T <= 6000 K
required from medium model \""     + mediumName + "\".");

        MM = molarMass(state);
        h = h_TX(T, X);
        R = data.R*X;
        u = h - R*T;
        d = p/(R*T);
        // connect state with BaseProperties
        state.T = T;
        state.p = p;
        state.X = if fixedX then reference_X else X;
      end BaseProperties;

        redeclare function setState_pTX
          "Return thermodynamic state as function of p, T and composition X"
          extends Modelica.Icons.Function;
          input AbsolutePressure p "Pressure";
          input Temperature T "Temperature";
          input MassFraction X[:]=reference_X "Mass fractions";
          output ThermodynamicState state;
        algorithm
        //   state := if size(X,1) == 0 then ThermodynamicState(p=p,T=T,X=reference_X, a=  Test.Media.Gas.MixtureGasNasa.velocityOfSound(state),
        //   beta=  Test.Media.Gas.MixtureGasNasa.isobaricExpansionCoefficient(state),cp=Test.Media.Gas.MixtureGasNasa.heatCapacity_cp(state), cv=Test.Media.Gas.MixtureGasNasa.heatCapacity_cv(state),
        //   d=Test.Media.Gas.MixtureGasNasa.density_pTX(p,state.T,reference_X),eta=Test.Media.Gas.MixtureGasNasa.dynamicViscosity(state), h=Test.Media.Gas.MixtureGasNasa.specificEnthalpy(state), kappa=Test.Media.Gas.MixtureGasNasa.kappa(state),
        //   lambda=Test.Media.Gas.MixtureGasNasa.thermalConductivity(state), s=  Test.Media.Gas.MixtureGasNasa.specificEntropy(state))
        //  else
        //   if size(X,1) == nX then ThermodynamicState(p=p,T=T,X=X,a=  Test.Media.Gas.MixtureGasNasa.velocityOfSound(state),beta=  Test.Media.Gas.MixtureGasNasa.isobaricExpansionCoefficient(state),
        //   cp=Test.Media.Gas.MixtureGasNasa.heatCapacity_cp(state), cv=Test.Media.Gas.MixtureGasNasa.heatCapacity_cv(state),d=Test.Media.Gas.MixtureGasNasa.density_pTX(p,state.T,reference_X),
        //   eta=Test.Media.Gas.MixtureGasNasa.dynamicViscosity(state), h=Test.Media.Gas.MixtureGasNasa.specificEnthalpy(state), kappa=Test.Media.Gas.MixtureGasNasa.kappa(state), lambda=Test.Media.Gas.MixtureGasNasa.thermalConductivity(state),
        //   s=  Test.Media.Gas.MixtureGasNasa.specificEntropy(state)) else
        //   ThermodynamicState(p=p,T=T, X=cat(1,X,{1-sum(X)}), a=  Test.Media.Gas.MixtureGasNasa.velocityOfSound(state),beta=  Test.Media.Gas.MixtureGasNasa.isobaricExpansionCoefficient(state),
        //   cp=Test.Media.Gas.MixtureGasNasa.heatCapacity_cp(state), cv=Test.Media.Gas.MixtureGasNasa.heatCapacity_cv(state),d=Test.Media.Gas.MixtureGasNasa.density_pTX(p,state.T,reference_X),
        //   eta=Test.Media.Gas.MixtureGasNasa.dynamicViscosity(state), h=Test.Media.Gas.MixtureGasNasa.specificEnthalpy(state), kappa=Test.Media.Gas.MixtureGasNasa.kappa(state), lambda=Test.Media.Gas.MixtureGasNasa.thermalConductivity(state),
        //   s=  Test.Media.Gas.MixtureGasNasa.specificEntropy(state));
        //   state :=if size(X, 1) == 0 then ThermodynamicState(
        //     p=p,
        //     T=T,
        //     X=reference_X,
        //     cp=specificHeatCapacityCp_TX(T, reference_X),
        //     cv=specificHeatCapacityCv_TX(T, reference_X),
        //     d=Test.Media.Gas.MixtureGasNasa.density_pTX2(T,p,reference_X))
        //     else if size(X, 1) == nX then ThermodynamicState(
        //     p=p,
        //     T=T,
        //     X=X,
        //     cp=specificHeatCapacityCp_TX(T, X),
        //     cv=specificHeatCapacityCv_TX(T, X),
        //     d=Test.Media.Gas.MixtureGasNasa.density_pTX2(T,p,X))
        //     else ThermodynamicState(
        //     p=p,
        //     T=T,
        //     X=cat(
        //       1,
        //       X,
        //       {1 - sum(X)}),
        //     cp=specificHeatCapacityCp_TX(T, X),
        //     cv=specificHeatCapacityCv_TX(T, X),
        //     d=Test.Media.Gas.MixtureGasNasa.density_pTX2(T,p,X));
          state.T       := T;
          state.p       := p;
          state.h       := h_TX(T, X);
          state.X       := X;
          state.cp      := specificHeatCapacityCp_TX(T, X);
          state.cv      := specificHeatCapacityCv_TX(T, X);
          state.d       := p/(X*data.R*T);
          state.lambda  := thermalConductivity_TX(T, X);
          state.s       := s_TX(T, X);
          state.eta     := dynamicViscosity_pTX(T, X);

        annotation(Inline=true,smoothOrder=2);
        end setState_pTX;

        redeclare function setState_phX
          "Return thermodynamic state as function of p, h and composition X"
          extends Modelica.Icons.Function;
          input AbsolutePressure p "Pressure";
          input SpecificEnthalpy h "Specific enthalpy";
          input MassFraction X[:]=reference_X "Mass fractions";
          output ThermodynamicState state;
        algorithm
        //   state := if size(X,1) == 0 then ThermodynamicState(p=p,T=T_hX(h,reference_X),X=reference_X, a=  Test.Media.Gas.MixtureGasNasa.velocityOfSound(state),
        //   beta=  Test.Media.Gas.MixtureGasNasa.isobaricExpansionCoefficient(state),cp=Test.Media.Gas.MixtureGasNasa.heatCapacity_cp(state), cv=Test.Media.Gas.MixtureGasNasa.heatCapacity_cv(state),
        //   d=Test.Media.Gas.MixtureGasNasa.density(state),eta=Test.Media.Gas.MixtureGasNasa.dynamicViscosity(state), h=h, kappa=Test.Media.Gas.MixtureGasNasa.kappa(state),
        //   lambda=Test.Media.Gas.MixtureGasNasa.thermalConductivity(state), s=  Test.Media.Gas.MixtureGasNasa.specificEntropy(state))
        //  else
        //   if size(X,1) == nX then ThermodynamicState(p=p,T=T_hX(h,X),X=X,a=  Test.Media.Gas.MixtureGasNasa.velocityOfSound(state),beta=  Test.Media.Gas.MixtureGasNasa.isobaricExpansionCoefficient(state),
        //   cp=Test.Media.Gas.MixtureGasNasa.specificHeatCapacityCp_TX(T_hX(h,X),X), cv=Test.Media.Gas.MixtureGasNasa.specificHeatCapacityCv_TX(T_hX(h,X),X),d=Test.Media.Gas.MixtureGasNasa.density(state),
        //   eta=Test.Media.Gas.MixtureGasNasa.dynamicViscosity(state), h=h, kappa=Test.Media.Gas.MixtureGasNasa.kappa(state), lambda=Test.Media.Gas.MixtureGasNasa.thermalConductivity(state),
        //   s=  Test.Media.Gas.MixtureGasNasa.specificEntropy(state)) else
        //   ThermodynamicState(p=p,T=T_hX(h,X), X=cat(1,X,{1-sum(X)}), a=  Test.Media.Gas.MixtureGasNasa.velocityOfSound(state),beta=  Test.Media.Gas.MixtureGasNasa.isobaricExpansionCoefficient(state),
        //   cp=Test.Media.Gas.MixtureGasNasa.heatCapacity_cp(state), cv=Test.Media.Gas.MixtureGasNasa.heatCapacity_cv(state),d=Test.Media.Gas.MixtureGasNasa.density(state),
        //   eta=Test.Media.Gas.MixtureGasNasa.dynamicViscosity(state), h=h, kappa=Test.Media.Gas.MixtureGasNasa.kappa(state), lambda=Test.Media.Gas.MixtureGasNasa.thermalConductivity(state),
        //   s=  Test.Media.Gas.MixtureGasNasa.specificEntropy(state));

          state.T      :=T_hX(h, X);
          state.p      :=p;
          state.h      :=h;
          state.X      :=X;
          state.cp     :=specificHeatCapacityCp_TX(state.T, X);
          state.cv     :=specificHeatCapacityCv_TX(state.T, X);
          state.d      :=p/(X*data.R*state.T);
          state.lambda := thermalConductivity_TX(state.T, X);
          state.s      :=s_TX(state.T, X);
          state.eta    := dynamicViscosity_pTX(state.T, X);

        //   state :=if size(X, 1) == 0 then ThermodynamicState(
        //     p=p,
        //     T=T_hX(h, reference_X),
        //     X=reference_X,
        //     cp=specificHeatCapacityCp_TX(T_hX(h, reference_X), reference_X),
        //     cv=specificHeatCapacityCv_TX(T_hX(h, reference_X), reference_X),
        //     d=density_pTX2(T_hX(h, reference_X),p,reference_X))
        //     else if size(X, 1) == nX then ThermodynamicState(
        //     p=p,
        //     T=T_hX(h, X),
        //     X=X,
        //     cp=specificHeatCapacityCp_TX(T_hX(h, X), X),
        //     cv=specificHeatCapacityCv_TX(T_hX(h, X), X),
        //     d=Test.Media.Gas.MixtureGasNasa.density_pTX2(T_hX(h, X),p,X))
        //     else ThermodynamicState(
        //     p=p,
        //     T=T_hX(h, X),
        //     X=cat(
        //       1,
        //       X,
        //       {1 - sum(X)}),
        //     cp=specificHeatCapacityCp_TX(T_hX(h, X), X),
        //     cv=specificHeatCapacityCv_TX(T_hX(h, X), X),
        //     d=Test.Media.Gas.MixtureGasNasa.density_pTX2(T_hX(h, X),p,X));
          annotation(Inline=true,smoothOrder=2);
        end setState_phX;

        redeclare function setState_psX
          "Return thermodynamic state as function of p, s and composition X"
          extends Modelica.Icons.Function;
          input AbsolutePressure p "Pressure";
          input SpecificEntropy s "Specific entropy";
          input MassFraction X[:]=reference_X "Mass fractions";
          output ThermodynamicState state;
        algorithm
          state := if size(X,1) == 0 then ThermodynamicState(p=p,T=T_psX(p,s,reference_X),X=reference_X) else if size(X,1) == nX then ThermodynamicState(p=p,T=T_psX(p,s,X),X=X) else
                 ThermodynamicState(p=p,T=T_psX(p,s,X), X=cat(1,X,{1-sum(X)}));
          annotation(Inline=true,smoothOrder=2);
        end setState_psX;

        redeclare function setState_dTX
          "Return thermodynamic state as function of d, T and composition X"
          extends Modelica.Icons.Function;
          input Density d "Density";
          input Temperature T "Temperature";
          input MassFraction X[:]=reference_X "Mass fractions";
          output ThermodynamicState state;
        algorithm
          state := if size(X,1) == 0 then ThermodynamicState(p=d*(data.R*reference_X)*T,T=T,X=reference_X) else if size(X,1) == nX then ThermodynamicState(p=d*(data.R*X)*T,T=T,X=X) else
                 ThermodynamicState(p=d*(data.R*cat(1,X,{1-sum(X)}))*T,T=T, X=cat(1,X,{1-sum(X)}));
          annotation(Inline=true,smoothOrder=2);
        end setState_dTX;

          redeclare function extends setSmoothState
          "Return thermodynamic state so that it smoothly approximates: if x > 0 then state_a else state_b"
          algorithm
          state := ThermodynamicState(
                      p=Modelica.Media.Common.smoothStep(
                        x,
                        state_a.p,
                        state_b.p,
                        x_small),
                      T=Modelica.Media.Common.smoothStep(
                        x,
                        state_a.T,
                        state_b.T,
                        x_small),
                      X=Modelica.Media.Common.smoothStep(
                        x,
                        state_a.X,
                        state_b.X,
                        x_small));
            annotation(Inline=true,smoothOrder=2);
          end setSmoothState;

        redeclare function extends pressure "Return pressure of ideal gas"
        algorithm
          p := state.p;
          annotation(Inline=true,smoothOrder=2);
        end pressure;

        redeclare function extends temperature
          "Return temperature of ideal gas"
        algorithm
          T := state.T;
          annotation(Inline=true,smoothOrder=2);
        end temperature;

        redeclare function extends density "Return density of ideal gas"
        algorithm
          d := state.p/((state.X*data.R)*state.T);
          annotation(Inline = true, smoothOrder = 3);
        end density;

      redeclare function extends specificEnthalpy "Return specific enthalpy"
        extends Modelica.Icons.Function;
      algorithm
        h := h_TX(state.T,state.X);
        annotation(Inline=true,smoothOrder=2);
      end specificEnthalpy;

      redeclare function extends specificInternalEnergy
          "Return specific internal energy"
        extends Modelica.Icons.Function;
      algorithm
        u := h_TX(state.T,state.X) - gasConstant(state)*state.T;
        annotation(Inline=true,smoothOrder=2);
      end specificInternalEnergy;

      redeclare function extends specificEntropy "Return specific entropy"
        protected
        Real[nX] Y(unit="mol/mol")=massToMoleFractions(state.X, data.MM)
            "Molar fractions";
      algorithm
      s :=  s_TX(state.T, state.X) - sum(state.X[i]*Modelica.Constants.R/MMX[i]*
          (if state.X[i]<Modelica.Constants.eps then Y[i] else
          Modelica.Math.log(Y[i]*state.p/reference_p)) for i in 1:nX);
        annotation(Inline=true,smoothOrder=2);
      end specificEntropy;

      redeclare function extends specificGibbsEnergy
          "Return specific Gibbs energy"
        extends Modelica.Icons.Function;
      algorithm
        g := h_TX(state.T,state.X) - state.T*specificEntropy(state);
        annotation(Inline=true,smoothOrder=2);
      end specificGibbsEnergy;

      redeclare function extends specificHelmholtzEnergy
          "Return specific Helmholtz energy"
        extends Modelica.Icons.Function;
      algorithm
        f := h_TX(state.T,state.X) - gasConstant(state)*state.T - state.T*specificEntropy(state);
        annotation(Inline=true,smoothOrder=2);
      end specificHelmholtzEnergy;

      function h_TX "Return specific enthalpy"
        import Modelica.Media.Interfaces.Choices;
         extends Modelica.Icons.Function;
         input Modelica.SIunits.Temperature T "Temperature";
         input MassFraction X[:]=reference_X
            "Independent Mass fractions of gas mixture";
         input Boolean exclEnthForm=excludeEnthalpyOfFormation
            "If true, enthalpy of formation Hf is not included in specific enthalpy h";
         input Modelica.Media.Interfaces.Choices.ReferenceEnthalpy
                                         refChoice=referenceChoice
            "Choice of reference enthalpy";
         input Modelica.SIunits.SpecificEnthalpy h_off=h_offset
            "User defined offset for reference enthalpy, if referenceChoice = UserDefined";
         output Modelica.SIunits.SpecificEnthalpy h
            "Specific enthalpy at temperature T";
      algorithm
        h :=(if fixedX then reference_X else X)*
             {Modelica.Media.IdealGases.Common.Functions.h_T(
                                data[i], T, exclEnthForm, refChoice, h_off) for i in 1:nX};
        annotation(Inline=false,smoothOrder=2);
      end h_TX;

      function h_TX_der "Return specific enthalpy derivative"
        import Modelica.Media.Interfaces.Choices;
         extends Modelica.Icons.Function;
         input Modelica.SIunits.Temperature T "Temperature";
         input MassFraction X[nX] "Independent Mass fractions of gas mixture";
         input Boolean exclEnthForm=excludeEnthalpyOfFormation
            "If true, enthalpy of formation Hf is not included in specific enthalpy h";
         input Modelica.Media.Interfaces.Choices.ReferenceEnthalpy
                                         refChoice=referenceChoice
            "Choice of reference enthalpy";
         input Modelica.SIunits.SpecificEnthalpy h_off=h_offset
            "User defined offset for reference enthalpy, if referenceChoice = UserDefined";
        input Real dT "Temperature derivative";
        input Real dX[nX] "Independent mass fraction derivative";
        output Real h_der "Specific enthalpy at temperature T";
      algorithm
        h_der := if fixedX then
          dT*sum((Modelica.Media.IdealGases.Common.Functions.cp_T(
                                     data[i], T)*reference_X[i]) for i in 1:nX) else
          dT*sum((Modelica.Media.IdealGases.Common.Functions.cp_T(
                                     data[i], T)*X[i]) for i in 1:nX)+
          sum((Modelica.Media.IdealGases.Common.Functions.h_T(
                                 data[i], T)*dX[i]) for i in 1:nX);
        annotation (Inline = false, smoothOrder=1);
      end h_TX_der;

      redeclare function extends gasConstant "Return gasConstant"
      algorithm
        R := data.R*state.X;
        annotation(Inline = true, smoothOrder = 3);
      end gasConstant;

      redeclare function extends specificHeatCapacityCp
          "Return specific heat capacity at constant pressure"
      algorithm
        cp := {Modelica.Media.IdealGases.Common.Functions.cp_T(
                                  data[i], state.T) for i in 1:nX}*state.X;
        annotation(Inline=true,smoothOrder=1);
      end specificHeatCapacityCp;

      redeclare function extends specificHeatCapacityCv
          "Return specific heat capacity at constant volume from temperature and gas data"
      algorithm
        cv := {Modelica.Media.IdealGases.Common.Functions.cp_T(
                                  data[i], state.T) for i in 1:nX}*state.X -data.R*state.X;
        annotation(Inline=true, smoothOrder = 1);
      end specificHeatCapacityCv;

      function MixEntropy "Return mixing entropy of ideal gases / R"
        extends Modelica.Icons.Function;
        input Modelica.SIunits.MoleFraction x[:] "Mole fraction of mixture";
        output Real smix "Mixing entropy contribution, divided by gas constant";
      algorithm
        smix := sum(if x[i] > Modelica.Constants.eps then -x[i]*Modelica.Math.log(x[i]) else
                         x[i] for i in 1:size(x,1));
        annotation(Inline=true,smoothOrder=2);
      end MixEntropy;

      function s_TX
          "Return temperature dependent part of the entropy, expects full entropy vector"
        extends Modelica.Icons.Function;
        input Temperature T "Temperature";
        input MassFraction[nX] X "Mass fraction";
        output SpecificEntropy s "Specific entropy";
      algorithm
        s := sum(Modelica.Media.IdealGases.Common.Functions.s0_T(
                                    data[i], T)*X[i] for i in 1:size(X,1));
        annotation(Inline=true,smoothOrder=2);
      end s_TX;

      redeclare function extends isentropicExponent
          "Return isentropic exponent"
      algorithm
        gamma := specificHeatCapacityCp(state)/specificHeatCapacityCv(state);
        annotation(Inline=true,smoothOrder=2);
      end isentropicExponent;

      redeclare function extends velocityOfSound "Return velocity of sound"
        extends Modelica.Icons.Function;
        input ThermodynamicState state "Properties at upstream location";
      algorithm
        a := sqrt(max(0,gasConstant(state)*state.T*specificHeatCapacityCp(state)/specificHeatCapacityCv(state)));
        annotation(Inline=true,smoothOrder=2);
      end velocityOfSound;

      function isentropicEnthalpyApproximation
          "Approximate method of calculating h_is from upstream properties and downstream pressure"
        extends Modelica.Icons.Function;
        input AbsolutePressure p2 "Downstream pressure";
        input ThermodynamicState state
            "Thermodynamic state at upstream location";
        output SpecificEnthalpy h_is "Isentropic enthalpy";
        protected
        SpecificEnthalpy h "Specific enthalpy at upstream location";
        SpecificEnthalpy h_component[nX]
            "Specific enthalpy at upstream location";
        IsentropicExponent gamma =  isentropicExponent(state)
            "Isentropic exponent";
        protected
        MassFraction[nX] X "Complete X-vector";
      algorithm
        X := if reducedX then cat(1,state.X,{1-sum(state.X)}) else state.X;
        h_component :={Modelica.Media.IdealGases.Common.Functions.h_T(
                                         data[i], state.T, excludeEnthalpyOfFormation,
          referenceChoice, h_offset) for i in 1:nX};
        h :=h_component*X;
        h_is := h + gamma/(gamma - 1.0)*(state.T*gasConstant(state))*
          ((p2/state.p)^((gamma - 1)/gamma) - 1.0);
        annotation(smoothOrder=2);
      end isentropicEnthalpyApproximation;

      redeclare function extends isentropicEnthalpy
          "Return isentropic enthalpy"
        input Boolean exact = false
            "Flag whether exact or approximate version should be used";
      algorithm
        h_is := if exact then specificEnthalpy_psX(p_downstream,specificEntropy(refState),refState.X) else
               isentropicEnthalpyApproximation(p_downstream,refState);
        annotation(Inline=true,smoothOrder=2);
      end isentropicEnthalpy;

    function gasMixtureViscosity
          "Return viscosities of gas mixtures at low pressures (Wilke method)"
      extends Modelica.Icons.Function;
      input MoleFraction[:] yi "Mole fractions";
      input MolarMass[:] M "Mole masses";
      input DynamicViscosity[:] eta "Pure component viscosities";
      output DynamicViscosity etam "Viscosity of the mixture";
        protected
      Real fi[size(yi,1),size(yi,1)];
    algorithm
      for i in 1:size(eta,1) loop
        assert(fluidConstants[i].hasDipoleMoment,"Dipole moment for " + fluidConstants[i].chemicalFormula +
           " not known. Can not compute viscosity.");
        assert(fluidConstants[i].hasCriticalData, "Critical data for "+ fluidConstants[i].chemicalFormula +
           " not known. Can not compute viscosity.");
        for j in 1:size(eta,1) loop
          if i==1 then
            fi[i,j] := (1 + (eta[i]/eta[j])^(1/2)*(M[j]/M[i])^(1/4))^2/(8*(1 + M[i]/M[j]))^(1/2);
          elseif j<i then
              fi[i,j] := eta[i]/eta[j]*M[j]/M[i]*fi[j,i];
            else
              fi[i,j] := (1 + (eta[i]/eta[j])^(1/2)*(M[j]/M[i])^(1/4))^2/(8*(1 + M[i]/M[j]))^(1/2);
          end if;
        end for;
      end for;
      etam := sum(yi[i]*eta[i]/sum(yi[j]*fi[i,j] for j in 1:size(eta,1)) for i in 1:size(eta,1));

      annotation (smoothOrder=2,
                 Documentation(info="<html>

<p>
Simplification of the kinetic theory (Chapman and Enskog theory)
approach neglecting the second-order effects.<br>
<br>
This equation has been extensively tested (Amdur and Mason, 1958;
Bromley and Wilke, 1951; Cheung, 1958; Dahler, 1959; Gandhi and Saxena,
1964; Ranz and Brodowsky, 1962; Saxena and Gambhir, 1963a; Strunk, et
al., 1964; Vanderslice, et al. 1962; Wright and Gray, 1962). In most
cases, only nonpolar mixtures were compared, and very good results
obtained. For some systems containing hydrogen as one component, less
satisfactory agreement was noted. Wilke's method predicted mixture
viscosities that were larger than experimental for the H2-N2 system,
but for H2-NH3, it underestimated the viscosities. <br>
Gururaja, et al. (1967) found that this method also overpredicted in
the H2-O2 case but was quite accurate for the H2-CO2 system. <br>
Wilke's approximation has proved reliable even for polar-polar gas
mixtures of aliphatic alcohols (Reid and Belenyessy, 1960). The
principal reservation appears to lie in those cases where Mi&gt;&gt;Mj
and etai&gt;&gt;etaj.<br>
</p>

</html>"));
    end gasMixtureViscosity;

        function dynamicViscosity_pTX "Return mixture dynamic viscosity"
          input Temperature T "Temperature";
          input MassFraction X[nX] "Mass fractions";
          output DynamicViscosity eta "dynamic viscosity";
        protected
          DynamicViscosity[nX] etaX "Component dynamic viscosities";
        algorithm
          for i in 1:nX loop
        etaX[i] := Modelica.Media.IdealGases.Common.Functions.dynamicViscosityLowPressure(
                                                             T,
                           fluidConstants[i].criticalTemperature,
                           fluidConstants[i].molarMass,
                           fluidConstants[i].criticalMolarVolume,
                           fluidConstants[i].acentricFactor,
                           fluidConstants[i].dipoleMoment);
          end for;
          eta := gasMixtureViscosity(massToMoleFractions(X,
                                 fluidConstants[:].molarMass),
                     fluidConstants[:].molarMass,
                     etaX);
          annotation (smoothOrder=2);
        end dynamicViscosity_pTX;

      function mixtureViscosityChung
          "Return the viscosity of gas mixtures without access to component viscosities (Chung, et. al. rules)"
      extends Modelica.Icons.Function;

        input Temperature T "Temperature";
        input Temperature[:] Tc "Critical temperatures";
        input MolarVolume[:] Vcrit "Critical volumes (m3/mol)";
        input Real[:] w "Acentric factors";
        input Real[:] mu "Dipole moments (debyes)";
        input MolarMass[:] MolecularWeights "Molecular weights (kg/mol)";
        input MoleFraction[:] y "Molar Fractions";
        input Real[:] kappa =  zeros(nX) "Association Factors";
        output DynamicViscosity etaMixture "Mixture viscosity (Pa.s)";
        protected
      constant Real[size(y,1)] Vc =  Vcrit*1000000 "Critical volumes (cm3/mol)";
      constant Real[size(y,1)] M =  MolecularWeights*1000
            "Molecular weights (g/mol)";
      Integer n = size(y,1) "Number of mixed elements";
      Real sigmam3 "Mixture sigma3 in Angstrom";
      Real sigma[size(y,1),size(y,1)];
      Real edivkm;
      Real edivk[size(y,1),size(y,1)];
      Real Mm;
      Real Mij[size(y,1),size(y,1)];
      Real wm "Accentric factor";
      Real wij[size(y,1),size(y,1)];
      Real kappam
            "Correlation for highly polar substances such as alcohols and acids";
      Real kappaij[size(y,1),size(y,1)];
      Real mum;
      Real Vcm;
      Real Tcm;
      Real murm "Dimensionless dipole moment of the mixture";
      Real Fcm "Factor to correct for shape and polarity";
      Real omegav;
      Real Tmstar;
      Real etam "Mixture viscosity in microP";
      algorithm
      //combining rules
      for i in 1:n loop
        for j in 1:n loop
          Mij[i,j] := 2*M[i]*M[j]/(M[i]+M[j]);
          if i==j then
            sigma[i,j] := 0.809*Vc[i]^(1/3);
            edivk[i,j] := Tc[i]/1.2593;
            wij[i,j] := w[i];
            kappaij[i,j] := kappa[i];
          else
            sigma[i,j] := (0.809*Vc[i]^(1/3)*0.809*Vc[j]^(1/3))^(1/2);
            edivk[i,j] := (Tc[i]/1.2593*Tc[j]/1.2593)^(1/2);
            wij[i,j] := (w[i] + w[j])/2;
            kappaij[i,j] := (kappa[i]*kappa[j])^(1/2);
          end if;
        end for;
      end for;
      //mixing rules
      sigmam3 := (sum(sum(y[i]*y[j]*sigma[i,j]^3 for j in 1:n) for i in 1:n));
      //(epsilon/k)m
      edivkm := (sum(sum(y[i]*y[j]*edivk[i,j]*sigma[i,j]^3 for j in 1:n) for i in 1:n))/sigmam3;
      Mm := ((sum(sum(y[i]*y[j]*edivk[i,j]*sigma[i,j]^2*Mij[i,j]^(1/2) for j in 1:n) for i in 1:n))/(edivkm*sigmam3^(2/3)))^2;
      wm := (sum(sum(y[i]*y[j]*wij[i,j]*sigma[i,j]^3 for j in 1:n) for i in 1:n))/sigmam3;
      mum := (sigmam3*(sum(sum(y[i]*y[j]*mu[i]^2*mu[j]^2/sigma[i,j]^3 for j in 1:n) for i in 1:n)))^(1/4);
      Vcm := sigmam3/(0.809)^3;
      Tcm := 1.2593*edivkm;
      murm := 131.3*mum/(Vcm*Tcm)^(1/2);
      kappam := (sigmam3*(sum(sum(y[i]*y[j]*kappaij[i,j] for j in 1:n) for i in 1:n)));
      Fcm := 1 - 0.275*wm + 0.059035*murm^4 + kappam;
      Tmstar := T/edivkm;
      omegav := 1.16145*(Tmstar)^(-0.14874) + 0.52487*Math.exp(-0.77320*Tmstar) + 2.16178*Math.exp(-2.43787*Tmstar);
      etam := 26.69*Fcm*(Mm*T)^(1/2)/(sigmam3^(2/3)*omegav);
      etaMixture := etam*1e7;

        annotation (smoothOrder=2,
                  Documentation(info="<html>

<p>
Equation to estimate the viscosity of gas mixtures at low pressures.<br>
It is a simplification of an extension of the rigorous kinetic theory
of Chapman and Enskog to determine the viscosity of multicomponent
mixtures, at low pressures and with a factor to correct for molecule
shape and polarity.
</p>

<p>
The input argument Kappa is a special correction for highly polar substances such as
alcohols and acids.<br>
Values of kappa for a few such materials:
</p>

<table style=\"text-align: left; width: 302px; height: 200px;\" border=\"1\"
cellspacing=\"0\" cellpadding=\"2\">
<tbody>
<tr>
<td style=\"vertical-align: top;\">Compound <br>
</td>
<td style=\"vertical-align: top; text-align: center;\">Kappa<br>
</td>
<td style=\"vertical-align: top;\">Compound<br>
</td>
<td style=\"vertical-align: top;\">Kappa<br>
</td>
</tr>
<tr>
<td style=\"vertical-align: top;\">Methanol<br>
</td>
<td style=\"vertical-align: top;\">0.215<br>
</td>
<td style=\"vertical-align: top;\">n-Pentanol<br>
</td>
<td style=\"vertical-align: top;\">0.122<br>
</td>
</tr>
<tr>
<td style=\"vertical-align: top;\">Ethanol<br>
</td>
<td style=\"vertical-align: top;\">0.175<br>
</td>
<td style=\"vertical-align: top;\">n-Hexanol<br>
</td>
<td style=\"vertical-align: top;\">0.114<br>
</td>
</tr>
<tr>
<td style=\"vertical-align: top;\">n-Propanol<br>
</td>
<td style=\"vertical-align: top;\">0.143<br>
</td>
<td style=\"vertical-align: top;\">n-Heptanol<br>
</td>
<td style=\"vertical-align: top;\">0.109<br>
</td>
</tr>
<tr>
<td style=\"vertical-align: top;\">i-Propanol<br>
</td>
<td style=\"vertical-align: top;\">0.143<br>
</td>
<td style=\"vertical-align: top;\">Acetic Acid<br>
</td>
<td style=\"vertical-align: top;\">0.0916<br>
</td>
</tr>
<tr>
<td style=\"vertical-align: top;\">n-Butanol<br>
</td>
<td style=\"vertical-align: top;\">0.132<br>
</td>
<td style=\"vertical-align: top;\">Water<br>
</td>
<td style=\"vertical-align: top;\">0.076<br>
</td>
</tr>
<tr>
<td style=\"vertical-align: top;\">i-Butanol<br>
</td>
<td style=\"vertical-align: top;\">0.132</td>
<td style=\"vertical-align: top;\"><br>
</td>
<td style=\"vertical-align: top;\"><br>
</td>
</tr>
</tbody>
</table>
<p>
Chung, et al. (1984) suggest that for other alcohols not shown in the
table:<br>
&nbsp;&nbsp;&nbsp;&nbsp; <br>
&nbsp;&nbsp;&nbsp; kappa = 0.0682 + 4.704*[(number of -OH
groups)]/[molecular weight]<br>
<br>
<span style=\"font-weight: normal;\">S.I. units relation for the
debyes:&nbsp;</span><br>
&nbsp;&nbsp;&nbsp; &nbsp;&nbsp; &nbsp;&nbsp; &nbsp;&nbsp; &nbsp;&nbsp;
&nbsp;&nbsp; &nbsp;&nbsp; &nbsp;&nbsp; &nbsp;&nbsp; &nbsp;&nbsp;
&nbsp;&nbsp; &nbsp;&nbsp; &nbsp;&nbsp; &nbsp;&nbsp; &nbsp;&nbsp;
&nbsp;&nbsp; &nbsp;&nbsp; &nbsp;&nbsp; 1 debye = 3.162e-25 (J.m^3)^(1/2)<br>
</p>
<h4>References</h4>
<p>
[1] THE PROPERTIES OF GASES AND LIQUIDS, Fifth Edition,<br>
&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp; &nbsp; Bruce E. Poling, John M.
Prausnitz, John P. O'Connell.<br>
[2] Chung, T.-H., M. Ajlan, L. L. Lee, and K. E. Starling: Ind. Eng.
Chem. Res., 27: 671 (1988).<br>
[3] Chung, T.-H., L. L. Lee, and K. E. Starling; Ing. Eng. Chem.
Fundam., 23: 3 ()1984).<br>
</p>
</html>"));
      end mixtureViscosityChung;

    function lowPressureThermalConductivity
          "Return thermal conductivities of low-pressure gas mixtures (Mason and Saxena Modification)"
      extends Modelica.Icons.Function;
      input MoleFraction[:] y
            "Mole fraction of the components in the gas mixture";
      input Temperature T "Temperature";
      input Temperature[:] Tc "Critical temperatures";
      input AbsolutePressure[:] Pc "Critical pressures";
      input MolarMass[:] M "Molecular weights";
      input ThermalConductivity[:] lambda
            "Thermal conductivities of the pure gases";
      output ThermalConductivity lambdam
            "Thermal conductivity of the gas mixture";
        protected
      MolarMass[size(y,1)] gamma;
      Real[size(y,1)] Tr "Reduced temperature";
      Real[size(y,1),size(y,1)] A "Mason and Saxena Modification";
      constant Real epsilon =  1.0 "Numerical constant near unity";
    algorithm
      for i in 1:size(y,1) loop
        gamma[i] := 210*(Tc[i]*M[i]^3/Pc[i]^4)^(1/6);
        Tr[i] := T/Tc[i];
      end for;
      for i in 1:size(y,1) loop
        for j in 1:size(y,1) loop
          A[i,j] := epsilon*(1 + (gamma[j]*(Math.exp(0.0464*Tr[i]) - Math.exp(-0.2412*Tr[i]))/
          (gamma[i]*(Math.exp(0.0464*Tr[j]) - Math.exp(-0.2412*Tr[j]))))^(1/2)*(M[i]/M[j])^(1/4))^2/
          (8*(1 + M[i]/M[j]))^(1/2);
        end for;
      end for;
      lambdam := sum(y[i]*lambda[i]/(sum(y[j]*A[i,j] for j in 1:size(y,1))) for i in 1:size(y,1));

      annotation (smoothOrder=2,
                  Documentation(info="<html>

<p>
This function applies the Masson and Saxena modification of the
Wassiljewa Equation for the thermal conductivity for gas mixtures of
n elements at low pressure.
</p>

<p>
For nonpolar gas mixtures errors will generally be less than 3 to 4%.
For mixtures of nonpolar-polar and polar-polar gases, errors greater
than 5 to 8% may be expected. For mixtures in which the sizes and
polarities of the constituent molecules are not greatly different, the
thermal conductivity can be estimated satisfactorily by a mole fraction
average of the pure component conductivities.
</p>

</html>"));
    end lowPressureThermalConductivity;

        redeclare replaceable function extends thermalConductivity
          "Return thermal conductivity for low pressure gas mixtures"
          input Integer method=methodForThermalConductivity
            "Method to compute single component thermal conductivity";
        protected
          ThermalConductivity[nX] lambdaX "Component thermal conductivities";
          DynamicViscosity[nX] eta "Component thermal dynamic viscosities";
          SpecificHeatCapacity[nX] cp "Component heat capacity";
        algorithm
          for i in 1:nX loop
        assert(fluidConstants[i].hasCriticalData, "Critical data for "+ fluidConstants[i].chemicalFormula +
           " not known. Can not compute thermal conductivity.");
        eta[i] := Modelica.Media.IdealGases.Common.Functions.dynamicViscosityLowPressure(
                                                            state.T,
                           fluidConstants[i].criticalTemperature,
                           fluidConstants[i].molarMass,
                           fluidConstants[i].criticalMolarVolume,
                           fluidConstants[i].acentricFactor,
                           fluidConstants[i].dipoleMoment);
        cp[i] := Modelica.Media.IdealGases.Common.Functions.cp_T(
                                    data[i],state.T);
        lambdaX[i] :=Modelica.Media.IdealGases.Common.Functions.thermalConductivityEstimate(
                                                               Cp=cp[i], eta=
              eta[i], method=method,data=data[i]);
          end for;
          lambda := lowPressureThermalConductivity(massToMoleFractions(state.X,
                                       fluidConstants[:].molarMass),
                               state.T,
                               fluidConstants[:].criticalTemperature,
                               fluidConstants[:].criticalPressure,
                               fluidConstants[:].molarMass,
                               lambdaX);
          annotation (smoothOrder=2);
        end thermalConductivity;

      redeclare function extends isobaricExpansionCoefficient
          "Return isobaric expansion coefficient beta"
      algorithm
        beta := 1/state.T;
        annotation(Inline=true,smoothOrder=2);
      end isobaricExpansionCoefficient;

      redeclare function extends isothermalCompressibility
          "Return isothermal compressibility factor"
      algorithm
        kappa := 1.0/state.p;
        annotation(Inline=true,smoothOrder=2);
      end isothermalCompressibility;

      redeclare function extends density_derp_T
          "Return density derivative by pressure at constant temperature"
      algorithm
        ddpT := 1/(state.T*gasConstant(state));
        annotation(Inline=true,smoothOrder=2);
      end density_derp_T;

      redeclare function extends density_derT_p
          "Return density derivative by temperature at constant pressure"
      algorithm
        ddTp := -state.p/(state.T*state.T*gasConstant(state));
        annotation(Inline=true,smoothOrder=2);
      end density_derT_p;

      redeclare function density_derX
          "Return density derivative by mass fraction"
        extends Modelica.Icons.Function;
        input ThermodynamicState state "Thermodynamic state record";
        output Density[nX] dddX "Derivative of density w.r.t. mass fraction";
      algorithm
        dddX := {-state.p/(state.T*gasConstant(state))*molarMass(state)/data[
          i].MM for i in 1:nX};
        annotation(Inline=true,smoothOrder=2);
      end density_derX;

      redeclare function extends molarMass "Return molar mass of mixture"
      algorithm
        MM := 1/sum(state.X[j]/data[j].MM for j in 1:size(state.X, 1));
        annotation(Inline=true,smoothOrder=2);
      end molarMass;

      function T_hX
          "Return temperature from specific enthalpy and mass fraction"
        extends Modelica.Icons.Function;
        input SpecificEnthalpy h "Specific enthalpy";
        input MassFraction[:] X "Mass fractions of composition";
         input Boolean exclEnthForm=excludeEnthalpyOfFormation
            "If true, enthalpy of formation Hf is not included in specific enthalpy h";
         input Modelica.Media.Interfaces.Choices.ReferenceEnthalpy
                                         refChoice=referenceChoice
            "Choice of reference enthalpy";
         input Modelica.SIunits.SpecificEnthalpy h_off=h_offset
            "User defined offset for reference enthalpy, if referenceChoice = UserDefined";
        output Temperature T "Temperature";
        protected
        MassFraction[nX] Xfull = if size(X,1) == nX then X else cat(1,X,{1-sum(X)});
      package Internal
            "Solve h(data,T) for T with given h (use only indirectly via temperature_phX)"
        extends Modelica.Media.Common.OneNonLinearEquation;
        redeclare record extends f_nonlinear_Data
              "Data to be passed to non-linear function"
          extends Modelica.Media.IdealGases.Common.DataRecord;
        end f_nonlinear_Data;

        redeclare function extends f_nonlinear
        algorithm
            y := h_TX(x,X);
        end f_nonlinear;

        // Dummy definition has to be added for current Dymola
        redeclare function extends solve
        end solve;
      end Internal;

      algorithm
        T := Internal.solve(h, 200, 6000, 1.0e5, Xfull, data[1]);
        annotation(inverse(h = h_TX(T,X,exclEnthForm,refChoice,h_off)));
      end T_hX;

      function T_psX
          "Return temperature from pressure, specific entropy and mass fraction"
        extends Modelica.Icons.Function;
        input AbsolutePressure p "Pressure";
        input SpecificEntropy s "Specific entropy";
        input MassFraction[:] X "Mass fractions of composition";
        output Temperature T "Temperature";
        protected
        MassFraction[nX] Xfull = if size(X,1) == nX then X else cat(1,X,{1-sum(X)});
      package Internal
            "Solve h(data,T) for T with given h (use only indirectly via temperature_phX)"
        extends Modelica.Media.Common.OneNonLinearEquation;
        redeclare record extends f_nonlinear_Data
              "Data to be passed to non-linear function"
          extends Modelica.Media.IdealGases.Common.DataRecord;
        end f_nonlinear_Data;

        redeclare function extends f_nonlinear
              "Note that this function always sees the complete mass fraction vector"
            protected
        MassFraction[nX] Xfull = if size(X,1) == nX then X else cat(1,X,{1-sum(X)});
        Real[nX] Y(unit="mol/mol")=massToMoleFractions(if size(X,1) == nX then X else cat(1,X,{1-sum(X)}), data.MM)
                "Molar fractions";
        algorithm
          y := s_TX(x,Xfull) - sum(Xfull[i]*Modelica.Constants.R/MMX[i]*
          (if Xfull[i]<Modelica.Constants.eps then Y[i] else
          Modelica.Math.log(Y[i]*p/reference_p)) for i in 1:nX);
            // s_TX(x,X)- data[:].R*X*(Modelica.Math.log(p/reference_p)
            //       + MixEntropy(massToMoleFractions(X,data[:].MM)));
        end f_nonlinear;

        // Dummy definition has to be added for current Dymola
        redeclare function extends solve
        end solve;
      end Internal;

      algorithm
        T := Internal.solve(s, 200, 6000, p, Xfull, data[1]);
      end T_psX;

    //   redeclare function extends specificEnthalpy_psX
    //   protected
    //     Temperature T "Temperature";
    //   algorithm
    //     T := temperature_psX(p,s,X);
    //     h := specificEnthalpy_pTX(p,T,X);
    //   end extends;

    //   redeclare function extends density_phX
    //     "Compute density from pressure, specific enthalpy and mass fraction"
    //     protected
    //     Temperature T "Temperature";
    //     SpecificHeatCapacity R "Gas constant";
    //   algorithm
    //     T := temperature_phX(p,h,X);
    //     R := if (not reducedX) then
    //       sum(data[i].R*X[i] for i in 1:size(substanceNames, 1)) else
    //       sum(data[i].R*X[i] for i in 1:size(substanceNames, 1)-1) + data[end].R*(1-sum(X[i]));
    //     d := p/(R*T);
    //   end density_phX;

      function velocityOfSound_pTX "Return velocity of sound"
        extends Modelica.Icons.Function;
        input ThermodynamicState state "Properties at upstream location";
      algorithm
        a := sqrt(max(0,gasConstant(state)*state.T*specificHeatCapacityCp(state)/specificHeatCapacityCv(state)));
        annotation(Inline=true,smoothOrder=2);
      end velocityOfSound_pTX;

      function specificHeatCapacityCp_TX
          "Return specific heat capacity at constant pressure"
        input Temperature T "Temperature";
        input MassFraction X[nX] "Mass fractions";
        output Modelica.SIunits.SpecificHeatCapacityAtConstantPressure cp
            "specific heat capacity cp";
      algorithm
        cp := {Modelica.Media.IdealGases.Common.Functions.cp_T(data[i], T)
        for i in 1:nX}*X;
        annotation(Inline=true,smoothOrder=1);
      end specificHeatCapacityCp_TX;

      function specificHeatCapacityCv_TX
          "Return specific heat capacity at constant volume from temperature and gas data"
        input Temperature T "Temperature";
        input MassFraction X[nX] "Mass fractions";
        output Modelica.SIunits.SpecificHeatCapacityAtConstantVolume cv
            "specific heat capacity cv";
      algorithm
        cv := {Modelica.Media.IdealGases.Common.Functions.cp_T(data[i], T)
        for i in 1:nX}*X -data.R*X;
        annotation(Inline=true, smoothOrder = 1);
      end specificHeatCapacityCv_TX;

        function thermalConductivity_TX
          "Return thermal conductivity for low pressure gas mixtures"
          input Temperature T "Temperature";
          input MassFraction X[nX] "Mass fractions";
          output ThermalConductivity lambda "thermal conductivity";
        protected
          ThermalConductivity[nX] lambdaX "Component thermal conductivities";
          DynamicViscosity[nX] eta "Component thermal dynamic viscosities";
          SpecificHeatCapacity[nX] cp "Component heat capacity";
        algorithm
          for i in 1:nX loop
        assert(fluidConstants[i].hasCriticalData, "Critical data for "+ fluidConstants[i].chemicalFormula +
           " not known. Can not compute thermal conductivity.");
        eta[i] := Modelica.Media.IdealGases.Common.Functions.dynamicViscosityLowPressure(
                                                            T,
                           fluidConstants[i].criticalTemperature,
                           fluidConstants[i].molarMass,
                           fluidConstants[i].criticalMolarVolume,
                           fluidConstants[i].acentricFactor,
                           fluidConstants[i].dipoleMoment);
        cp[i] := Modelica.Media.IdealGases.Common.Functions.cp_T(
                                    data[i],T);
        lambdaX[i] :=Modelica.Media.IdealGases.Common.Functions.thermalConductivityEstimate(
                                                               Cp=cp[i], eta=
              eta[i], method=methodForThermalConductivity,data=data[i]);
          end for;
          lambda := lowPressureThermalConductivity(massToMoleFractions(X,
                                       fluidConstants[:].molarMass),
                               T,
                               fluidConstants[:].criticalTemperature,
                               fluidConstants[:].criticalPressure,
                               fluidConstants[:].molarMass,
                               lambdaX);
          annotation (smoothOrder=2);
        end thermalConductivity_TX;

        function density_pTX2 "Return density of ideal gas"
          input Temperature T "Temperature";
          input AbsolutePressure p "pressure";
          input MassFraction X[nX] "Mass fractions";
          output Modelica.SIunits.Density d "density";
        algorithm
          d := p/(X*data.R*T);
          annotation(Inline = true, smoothOrder = 3);
        end density_pTX2;
      annotation (Documentation(info="<HTML>
<p>
This model calculates the medium properties for single component ideal gases.
</p>
<p>
<b>Sources for model and literature:</b><br>
Original Data: Computer program for calculation of complex chemical
equilibrium compositions and applications. Part 1: Analysis
Document ID: 19950013764 N (95N20180) File Series: NASA Technical Reports
Report Number: NASA-RP-1311  E-8017  NAS 1.61:1311
Authors: Gordon, Sanford (NASA Lewis Research Center)
 Mcbride, Bonnie J. (NASA Lewis Research Center)
Published: Oct 01, 1994.
</p>
<p><b>Known limits of validity:</b></br>
The data is valid for
temperatures between 200 K and 6000 K.  A few of the data sets for
monatomic gases have a discontinuous 1st derivative at 1000 K, but
this never caused problems so far.
</p>
<p>
This model has been copied from the ThermoFluid library.
It has been developed by Hubertus Tummescheit.
</p>
</HTML>"));
    end MixtureGasNasa;
    end Gas;
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
        hh_in_start=1e3)
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
        dp_c=0,
        hc_in_start=1e3) "condenser model"  annotation (Placement(transformation(
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
        m_flow=0.5,
        use_T=true,
        use_p=true,
        use_m_flow=true,
        redeclare package Medium = Hot,
        T=618.15,
        p=200000) annotation (Placement(transformation(
            extent={{-11,-10},{11,10}},
            rotation=0,
            origin={-11,98})));
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
          points={{0.11,98},{0.11,86},{0,86},{0,81.51}},
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
      package Hot = Media.Gas.FlueGas "Medium model hot cells";
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
        hh_in_start=4.8e5)
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
        use_dT_int=true,
        hc_in_start=1e3) "condenser model"  annotation (Placement(transformation(
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
          package Medium_hot = Media.CoolProp.THERM66 "Medium model";
          package Medium_cold = Media.RefProp.MM "Medium model";

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
            N_cell_pc=7,
            mdot_cold_start=0.15,
            mdot_hot_start=0.19,
            l_start=0.6,
            beta=1.0471975511966,
            t_hot_in_start=618.15,
            t_hot_out_start=505.05,
            t_cold_in_start=490.05,
            t_cold_out_start=598.15,
            p_hot_start=200000,
            p_cold_start=1460000)
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
            t_hot_out_start=505.05,
            t_cold_in_start=490.05,
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

    model ORCHID "Design of the ORCHID test rig. No design of the components."
      parameter Modelica.SIunits.Temperature TIT = 320 + 273.15
        "Temperature difference at the inlet of the evaporator";
      parameter Modelica.SIunits.TemperatureDifference dT_int_rec = 20
        "Temperature difference at the inlet of the recuperator";
      parameter Modelica.SIunits.AbsolutePressure pmax = 15e5
        "Evaporation pressure";
      parameter Modelica.SIunits.AbsolutePressure pmin = 0.33e5
        "Condensation pressure";
      parameter Medium.SaturationProperties  sat =  Medium.setSat_p(pmin)
        "Saturation properties";
      replaceable package Hot = Media.CoolProp.THERM66 "Medium model hot cells";
      replaceable package Medium = Media.CoolProp.Toluene
        "Medium model hot cells";
      replaceable package Coolant = Media.CoolProp.Water_TTSE
        "Medium model hot cells";

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
        use_dT_int=true,
        hc_in_start=1E3) "condenser model"  annotation (Placement(transformation(
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
        dT_out=dT_int_rec,
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
        p=pmax) annotation (Placement(transformation(
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
        use_T=true,
        T=593.15,
        use_h=false)
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
        p=pmin,
        use_T=true,
        use_p=false)
                 annotation (Placement(transformation(
            extent={{-9,-9},{9,9}},
            rotation=0,
            origin={-88,-105})));
      CycleTempo.Components.Flags.START START1(
        redeclare package Medium = Medium,
        m_flow=0.15,
        p=pmax,
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
        p=pmin,
        T=343.15)
        annotation (Placement(transformation(extent={{26,-50},{8,-32}})));
      CycleTempo.Components.Flags.START START3(   redeclare package Medium = Medium,
        m_flow=0.15,
        p=pmin,
        T=593.15)
        annotation (Placement(transformation(extent={{34,41},{52,59}})));

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

    model ORCHID_opti
      "Design of the ORCHID test rig. No design of the components. Optimization on."
      replaceable package Hot = Media.CoolProp.THERM66 "Medium model hot cells";
      replaceable package Medium = Media.CoolProp.Toluene_TTSE
        "Medium model hot cells";
      replaceable package Coolant = Media.CoolProp.Water_TTSE
        "Medium model hot cells";
      input Modelica.SIunits.Temperature TIT
        "Temperature difference at the inlet of the evaporator";
      input Modelica.SIunits.TemperatureDifference dT_int_rec
        "Temperature difference at the inlet of the recuperator";
      input Modelica.SIunits.AbsolutePressure pmax "Evaporation pressure";
      input Modelica.SIunits.AbsolutePressure pmin "Condensation pressure";
      Medium.SaturationProperties  sat "Saturation properties";

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
        use_dT_int=true,
        hc_in_start=1E3) "condenser model"  annotation (Placement(transformation(
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
        dp_c=0)
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
        redeclare package Medium = Medium)
                annotation (Placement(transformation(
            extent={{10,-10},{-10,10}},
            rotation=0,
            origin={87,-22})));
      CycleTempo.Components.Flags.ADDCO N6(
        redeclare package Medium = Medium)
                  annotation (Placement(transformation(
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
        use_T=false)
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
        redeclare package Medium = Medium)
                 annotation (Placement(transformation(
            extent={{-9,-9},{9,9}},
            rotation=0,
            origin={-88,-105})));
      CycleTempo.Components.Flags.START START1(
        redeclare package Medium = Medium,
        m_flow=0.15,
        p=15e5,
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
        p=1e4,
        T=343.15)
        annotation (Placement(transformation(extent={{26,-50},{8,-32}})));
      CycleTempo.Components.Flags.START START3(   redeclare package Medium = Medium,
        m_flow=0.15,
        p=1e4,
        T=593.15)
        annotation (Placement(transformation(extent={{34,41},{52,59}})));

      CycleTempo.Components.Flags.START START4(   redeclare package Medium = Medium,
        m_flow=0.15,
        p=1e4)
        annotation (Placement(transformation(extent={{20,19},{38,37}})));
      CycleTempo.Components.Flags.START START5(   redeclare package Medium = Medium,
        m_flow=0.15,
        p=1e4,
        T=323)
        annotation (Placement(transformation(extent={{118,-78},{100,-60}})));
    equation
      sat        =  Medium.setSat_p(pmin);
      pmax       = N2.node.p;
      N4.node.h  = Medium.specificEnthalpy_pT(N4.node.p, TIT);
      dT_int_rec = Medium.temperature_ph(N6.node.p, N6.node.h) - Medium.temperature_ph(N2.node.p, N2.node.h);
      pmin       = N1.node.p;

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
    end ORCHID_opti;
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

  package Optimization "Package for testing GenOpt"
    model Optimize_ORCHID "Optimization of the ORCHID plant"
      parameter String parameterFileName = "modelicaParameters.txt"
        "File on which data is present";
      parameter String inputFileName = "modelicaSchedule.txt"
        "File on which data is present";
      parameter String resultFileName = "result.txt"
        "File on which data is present";
      parameter String header = "Objective function value"
        "Header for result file";
      Modelica.Blocks.Sources.CombiTimeTable X(
        tableOnFile=true,
        columns[:]=2:5,
        fileName=inputFileName,
        tableName="tab1") "Table with control input";
      Merge.ORCHID_opti ORCHID(TIT=Modelica.SIunits.Conversions.from_degC(X.y[2]), dT_int_rec=X.y[3], pmax=Modelica.SIunits.Conversions.from_bar(X.y[1]), pmin=Modelica.SIunits.Conversions.from_bar(X.y[4]))
        "ORCHID ORC power system"
        annotation (Placement(transformation(extent={{-34,-30},{36,30}})));
    initial algorithm
     if (resultFileName <> "") then
        Modelica.Utilities.Files.removeFile(resultFileName);
      end if;
      Modelica.Utilities.Streams.print(fileName=resultFileName, string=header);
    equation
      when terminal() then
      Modelica.Utilities.Streams.print("f(x) = " +
      realString(number=-ORCHID.generator.terminal_out.W, minimumWidth=1, precision=16), resultFileName);
      end when;
    end Optimize_ORCHID;
  end Optimization;

  package WHR_TRUCK "Waste heat recovery from truck engines"

    model Cycle "Waste Heat
Recovery From a Heavy-Duty
Truck Engine by Means of an
ORC Turbogenerator"
      package Hot = Media.CoolProp.Air_TTSE;
      package Medium = Media.CoolProp.Toluene_TTSE;
      package Coolant = Media.CoolProp.Water_TTSE;

      CycleTempo.Components.HEX.Evaporator EVA_EGR(
        dp_h=0,
        redeclare package Medium_c = Medium,
        redeclare package Medium_h = Hot,
        dT_int=30,
        use_dT_int=false,
        hh_in_start=5.5e5,
        dp_c=7000) "Evaporator exhaust gas ricirculation"        annotation (
          Placement(transformation(
            extent={{-18,-15.5},{18,15.5}},
            rotation=0,
            origin={22,71.5})));
      CycleTempo.Components.Flags.ADDCO N10(
        redeclare package Medium = Hot,
        use_T=true,
        use_p=true,
        use_m_flow=true,
        m_flow=0.066,
        T=673.15,
        p=120000)         annotation (Placement(transformation(
            extent={{-11,-10},{11,10}},
            rotation=0,
            origin={11,122})));
      CycleTempo.Components.Flags.ADDCO N11(
        redeclare package Medium = Hot,
        visible=true,
        use_T=true,
        T=473.15) annotation (Placement(transformation(
            extent={{-10,-9},{10,9}},
            rotation=0,
            origin={12,22})));
      CycleTempo.Components.Flow.Splitter splitter1(redeclare package Medium =
            Test.Media.CoolProp.Toluene_TTSE)
        annotation (Placement(transformation(extent={{14,-15},{-14,15}},
            rotation=90,
            origin={-74,139})));
      CycleTempo.Components.Flags.ADDCO N1(
        redeclare package Medium = Hot,
        visible=true,
        use_T=true,
        T=473.15) annotation (Placement(transformation(
            extent={{-10,-9},{10,9}},
            rotation=0,
            origin={10,154})));
      CycleTempo.Components.HEX.Evaporator EVA_exhaust(
        redeclare package Medium_c = Medium,
        redeclare package Medium_h = Hot,
        dT_int=30,
        use_dT_int=false,
        hh_in_start=5.5e5,
        dp_h=0,
        dp_c=7000) "Evaporator fed by the exhaust gases"        annotation (
          Placement(transformation(
            extent={{-18,-15.5},{18,15.5}},
            rotation=0,
            origin={20,203.5})));
      CycleTempo.Components.Flags.ADDCO N2(
        redeclare package Medium = Hot,
        use_T=true,
        use_p=true,
        use_m_flow=true,
        m_flow=0.131,
        T=587.15,
        p=120000)         annotation (Placement(transformation(
            extent={{-11,-10},{11,10}},
            rotation=0,
            origin={9,254})));
      CycleTempo.Components.Flags.ADDCO N4(
        redeclare package Medium = Medium,
        use_T=true,
        T=553.15) annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=0,
            origin={130,226})));
      CycleTempo.Components.Flags.ADDCO N5(
        redeclare package Medium = Medium,
        use_p=true,
        use_T=false,
        T=433.15,
        p=700000) annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=0,
            origin={-114,139})));
      CycleTempo.Components.Flow.Mixer mixer(redeclare package Medium = Medium)
        annotation (Placement(transformation(
            extent={{-17,17},{17,-17}},
            rotation=90,
            origin={157,140})));
      CycleTempo.Components.Flags.ADDCO N3(
        redeclare package Medium = Medium,
        use_T=true,
        T=553.15) annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=0,
            origin={130,48})));
      CycleTempo.Components.Turbomachinery.Turbine turbine(
        redeclare package Medium = Medium,
        eta_m=1,
        eta_is=0.7) "Axial turbine"
        annotation (Placement(transformation(extent={{166,6},{230,60}})));
      CycleTempo.Components.Flags.ADDCO N6(
        redeclare package Medium = Medium,
        use_p=false,
        T=358.15,
        p=10000,
        use_T=true)
                  annotation (Placement(transformation(
            extent={{10,-10},{-10,10}},
            rotation=0,
            origin={282,-115})));
      CycleTempo.Components.Flags.ADDCOW Pout(W=33e3, use_W=false)
        annotation (Placement(transformation(extent={{334,11},{288,55}})));
      CycleTempo.Components.HEX.Heat_exchanger recuperator(
        redeclare package Medium_h = Medium,
        redeclare package Medium_c = Medium,
        dT_out=20,
        use_dT_out=false,
        dp_h=500,
        dp_c=7000)
        annotation (Placement(transformation(extent={{35,-47},{6,-19}})));
      CycleTempo.Components.HEX.Condenser condenser(
        redeclare package Medium_h = Medium,
        redeclare package Medium_c = Coolant,
        dp_c=0,
        use_dT_int=true,
        hc_in_start=1e3,
        dp_h=500,
        dT_int=10) "condenser model"        annotation (Placement(transformation(
            extent={{15,-14},{-15,14}},
            rotation=0,
            origin={20,-98})));
      CycleTempo.Components.Turbomachinery.Pump Pump(
        redeclare package Medium = Medium,
        eta_m=1,
        eta_is=0.65) "Centrifugal pump"
        annotation (Placement(transformation(extent={{19.5,20.5},{-19.5,-20.5}},
            rotation=-90,
            origin={200.5,-71.5})));
      CycleTempo.Components.Flags.CC CC(redeclare package Medium = Medium)
        annotation (Placement(transformation(
            extent={{13,-12.5},{-13,12.5}},
            rotation=-90,
            origin={200.5,-117})));
      CycleTempo.Components.Flags.ADDCO N20(
        redeclare package Medium = Coolant,
        use_T=true,
        use_p=true,
        T=343.15,
        p=120000)
        annotation (Placement(transformation(extent={{113,-108},{95,-89}})));
      CycleTempo.Components.Flags.ADDCO N21(
        redeclare package Medium = Coolant,
        T=333.15,
        use_T=false)
                  annotation (Placement(transformation(extent={{-77,-107},{-55,-89}})));
      CycleTempo.Components.Flags.ADDCO N7(
        redeclare package Medium = Medium,
        use_T=true,
        T=373.15) annotation (Placement(transformation(
            extent={{10.5,10},{-10.5,-10}},
            rotation=180,
            origin={-66.5,-55})));
      CycleTempo.Components.Flags.START START1(
        redeclare package Medium = Medium,
        m_flow=0.03,
        p=1460000,
        T=439.15)
        annotation (Placement(transformation(extent={{-104,195},{-70,225}})));
      CycleTempo.Components.Flags.START START2(
        redeclare package Medium = Medium,
        p=1460000,
        T=439.15,
        m_flow=0.06)
        annotation (Placement(transformation(extent={{-143,173},{-109,203}})));
      CycleTempo.Components.Electrics.Motor motor(eta_el=0.9)
        annotation (Placement(transformation(extent={{230,-87},{260,-56}})));
      CycleTempo.Components.Flags.ADDCOW Ppump(W=33e3, use_W=false)
        annotation (Placement(transformation(extent={{298,-87},{270,-56}})));
      CycleTempo.Components.Flags.START START3(
        redeclare package Medium = Medium,
        m_flow=0.06,
        p=30000,
        T=373.15)
        annotation (Placement(transformation(extent={{143,-105},{177,-75}})));
      CycleTempo.Components.Mechanics.GearBox gearBox(eta_m=0.98)
        annotation (Placement(transformation(extent={{242,16},{274,47}})));
      CycleTempo.Components.Flags.ADDCO N8(
        redeclare package Medium = Medium,
        T=553.15,
        use_T=false)
                  annotation (Placement(transformation(
            extent={{10,-10},{-10,10}},
            rotation=0,
            origin={280,0})));
      CycleTempo.Components.Flags.ADDCO N9(
        redeclare package Medium = Medium,
        use_T=false,
        T=553.15) annotation (Placement(transformation(
            extent={{10,-10},{-10,10}},
            rotation=0,
            origin={234,-38})));
    equation
      connect(N11.node, EVA_EGR.node_h_out) annotation (Line(
          points={{22.1,22},{22.1,42},{22,42},{22,49.8}},
          color={0,0,0},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(N10.node, EVA_EGR.node_h_in) annotation (Line(
          points={{22.11,122},{22.11,98},{22,98},{22,93.51}},
          color={0,0,0},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(N1.node, EVA_exhaust.node_h_out) annotation (Line(
          points={{20.1,154},{20.1,168},{20,168},{20,181.8}},
          color={0,0,0},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(N2.node, EVA_exhaust.node_h_in) annotation (Line(
          points={{20.11,254},{20,254},{20,225.51}},
          color={0,0,0},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(N5.node, splitter1.node_in) annotation (Line(
          points={{-103.9,139},{-77,139}},
          color={0,0,0},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(splitter1.node_out_1, EVA_EGR.node_c_in) annotation (Line(
          points={{-59,133.4},{-59,71.5},{-3.2,71.5}},
          color={0,0,0},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(splitter1.node_out_2, EVA_exhaust.node_c_in) annotation (Line(
          points={{-59,144.6},{-59,203.5},{-5.2,203.5}},
          color={0,0,0},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(N4.node, mixer.node_in_1) annotation (Line(
          points={{140.1,226},{140,226},{140,146.8}},
          color={0,0,0},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(EVA_exhaust.node_c_out, mixer.node_in_1) annotation (Line(
          points={{45.38,203.5},{100,203.5},{100,146.8},{140,146.8}},
          color={0,0,0},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(EVA_EGR.node_c_out, mixer.node_in_2) annotation (Line(
          points={{47.38,71.5},{100,71.5},{100,133.2},{140,133.2}},
          color={0,0,0},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(N3.node, mixer.node_in_2) annotation (Line(
          points={{140.1,48},{140,48},{140,133.2}},
          color={0,0,0},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(mixer.node_out, turbine.node_in) annotation (Line(
          points={{160.4,140},{185.2,140},{185.2,60}},
          color={0,0,0},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(recuperator.node_h_in, turbine.node_out) annotation (Line(
          points={{20.5,-13.12},{20.5,0.6},{202.48,0.6}},
          color={0,0,0},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(recuperator.node_h_out, condenser.node_h_in) annotation (Line(
          points={{20.5,-52.6},{20.5,-68.8},{20,-68.8},{20,-78.12}},
          color={0,0,0},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(recuperator.node_c_out, N5.node) annotation (Line(
          points={{0.055,-33},{-103.9,-33},{-103.9,139}},
          color={0,0,0},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(recuperator.node_c_in, Pump.node_out) annotation (Line(
          points={{40.8,-33},{200.5,-33},{200.5,-52}},
          color={0,0,0},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(Pump.node_in, CC.node_out) annotation (Line(
          points={{200.5,-90.61},{200.5,-109.2}},
          color={0,0,0},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(CC.node_in, condenser.node_h_out) annotation (Line(
          points={{200.5,-124.54},{200.5,-138},{20,-138},{20,-117.6}},
          color={0,0,0},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(Pump.node_in, N6.node) annotation (Line(
          points={{200.5,-90.61},{241.25,-90.61},{241.25,-115},{271.9,-115}},
          color={0,0,0},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(condenser.node_c_in, N20.node) annotation (Line(
          points={{41,-98},{68,-98},{68,-98.5},{94.91,-98.5}},
          color={0,0,0},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(N21.node, condenser.node_c_out) annotation (Line(
          points={{-54.89,-98},{-1.15,-98}},
          color={0,0,0},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(N7.node, recuperator.node_h_out) annotation (Line(
          points={{-55.895,-55},{-18,-55},{-18,-52.6},{20.5,-52.6}},
          color={0,0,0},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(START1.node, EVA_exhaust.node_c_in) annotation (Line(
          points={{-70,210},{-5.2,210},{-5.2,203.5}},
          color={0,0,0},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(START2.node, N5.node) annotation (Line(
          points={{-109,188},{-103.9,188},{-103.9,139}},
          color={0,0,0},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(Pump.terminal, motor.terminal_in) annotation (Line(
          points={{221,-71.5975},{225.5,-71.5975},{225.5,-71.5},{230.3,-71.5}},
          color={0,0,0},
          smooth=Smooth.None));
      connect(motor.terminal_out, Ppump.terminal) annotation (Line(
          points={{260,-71.5},{270,-71.5}},
          color={0,0,0},
          smooth=Smooth.None));
      connect(START3.node, Pump.node_in) annotation (Line(
          points={{177,-90},{187.5,-90},{187.5,-90.61},{200.5,-90.61}},
          color={0,0,0},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(turbine.terminal, gearBox.terminal_in) annotation (Line(
          points={{229.84,33},{235.92,33},{235.92,33.05},{242,33.05}},
          color={0,0,0},
          smooth=Smooth.None));
      connect(gearBox.terminal_out, Pout.terminal) annotation (Line(
          points={{277.2,33.05},{281.5,33.05},{281.5,33},{288,33}},
          color={0,0,0},
          smooth=Smooth.None));
      connect(N8.node, turbine.node_out) annotation (Line(
          points={{269.9,0},{236,0},{236,0.6},{202.48,0.6}},
          color={0,0,0},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(N9.node, Pump.node_out) annotation (Line(
          points={{223.9,-38},{212,-38},{212,-52},{200.5,-52}},
          color={0,0,0},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-160,
                -140},{320,280}}),      graphics), Icon(coordinateSystem(extent={{-160,
                -140},{320,280}})),
        experiment(__Dymola_NumberOfIntervals=1, __Dymola_Algorithm="Dassl"),
        __Dymola_experimentSetupOutput);
    end Cycle;

    model EVA_EGR "Design of the evaporator for EGR"
      package Hot = Test.Media.Gas.FlueGas;
      package Cold = Media.CoolProp.Toluene_TTSE;

        Design.Components.HEX.Flat_plate EVA_EGR(
        redeclare package Medium_hot = Hot,
        redeclare package Medium_cold = Cold,
        redeclare function cost = Design.Miscellanea.Cost.FP_Rafferty,
        ht_hot_f1=1e9,
        ht_cold_f1=1e9,
        X=1,
        w=0.325,
        d_pt=100e-3,
        redeclare Design.Heat_transfer.Plates.evaporation_Martin ht_cold(
          p_in=EVA_EGR.node_c_in.p,
          At=EVA_EGR.At,
          beta=EVA_EGR.beta,
          qdot=EVA_EGR.plate.qdot_cold,
          qdot_tilde_start=EVA_EGR.mdot_hot_start*(EVA_EGR.h_hot_in_start - EVA_EGR.h_hot_out_start)
              /(EVA_EGR.l_start*EVA_EGR.w*EVA_EGR.N_ch_p)),
        redeclare Design.Materials.SS_AISI_304 material,
        thick=0.8e-3,
        b=2.3e-3,
        redeclare Design.Pressure_drops.Plates.evaporation_Martin dp_cold(p_in=
              EVA_EGR.node_c_in.p, beta=EVA_EGR.beta),
        redeclare Design.Miscellanea.check_velocity check_hot(T_sat=500),
        redeclare Design.Miscellanea.check_velocity check_cold(umin=min(EVA_EGR.ht_cold.u)),
        mdot_cold_start=0.0286474,
        mdot_hot_start=0.066,
        N_ch_p=8,
        l_start=0.5,
        N_cell_pc=7,
        use_dp=true,
        beta=1.0471975511966,
        t_hot_in_start=673.15,
        t_hot_out_start=473.15,
        t_cold_in_start=460.55,
        t_cold_out_start=553.15,
        p_hot_start=120000,
        p_cold_start=700000)
        "Model of a flat plate heat exchanger. Evaporator for the ricirculated exhaust gases."
        annotation (Placement(transformation(extent={{-88,-82},{88,60}})));

    //      redeclare Design.Miscellanea.topology_PHE.two_pass_two_pass tpg_hot(N_ch_p=
    //            EVA_EGR.N_ch_p, stype=1),
    //      redeclare Design.Miscellanea.topology_PHE.two_pass_two_pass tpg_cold(N_ch_p=
    //           EVA_EGR.N_ch_p, stype=2),

      Design.Components.Flags.ADDCO HOT_OUT(
        redeclare package Medium = Hot,
        T=473.15,
        use_T=false)
        annotation (Placement(transformation(extent={{-133,-8},{-100,22}})));
      Design.Components.Flags.ADDCO COLD_IN(
        use_T=true,
        use_p=true,
        redeclare package Medium = Cold,
        m_flow=0.0286474,
        use_m_flow=true,
        T=460.55,
        p=700000)
        annotation (Placement(transformation(extent={{-133,-116},{-104,-88}})));
      Design.Components.Flags.ADDCO COLD_OUT(
        redeclare package Medium = Cold,
        T=553.15,
        use_T=true)
        annotation (Placement(transformation(extent={{79,-116},{108,-88}})));
      Design.Components.Flags.ADDCO HOT_IN(
        use_T=true,
        use_p=true,
        redeclare package Medium = Hot,
        m_flow=0.066,
        T=673.15,
        p=120000,
        use_m_flow=true)
                  annotation (Placement(transformation(extent={{49,100},{78,128}})));
    equation
      connect(HOT_OUT.node, EVA_EGR.node_h_out) annotation (Line(
          points={{-99.835,7},{-99.835,-53.6},{-88,-53.6}},
          color={0,0,0},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(COLD_IN.node, EVA_EGR.node_c_in) annotation (Line(
          points={{-103.855,-102},{-104,-102},{-104,-75.61},{-88,-75.61}},
          color={0,0,0},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(EVA_EGR.node_c_out, COLD_OUT.node) annotation (Line(
          points={{88,24.145},{130,24.145},{130,-102},{108.145,-102}},
          color={0,0,0},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(HOT_IN.node, EVA_EGR.node_h_in) annotation (Line(
          points={{78.145,114},{130,114},{130,45.8},{88,45.8}},
          color={0,0,0},
          pattern=LinePattern.None,
          smooth=Smooth.None));
     annotation (Placement(transformation(extent={{-108,-74},{88,66}})),
        experiment(__Dymola_NumberOfIntervals=1),
        __Dymola_experimentSetupOutput,
        Diagram(coordinateSystem(extent={{-160,-120},{140,140}},
              preserveAspectRatio=false), graphics),
        Icon(coordinateSystem(extent={{-160,-120},{140,140}})));
    end EVA_EGR;

    model EVA_Exh "Design of the evaporator for the exhaust gases"
      package Hot = Test.Media.Gas.FlueGas;
      package Cold = Media.CoolProp.Toluene_TTSE;

        Design.Components.HEX.Flat_plate EVA_EXH(
        redeclare package Medium_hot = Hot,
        redeclare package Medium_cold = Cold,
        redeclare function cost = Design.Miscellanea.Cost.FP_Rafferty,
        ht_hot_f1=1e9,
        ht_cold_f1=1e9,
        X=1,
        w=0.325,
        d_pt=100e-3,
        redeclare Design.Heat_transfer.Plates.evaporation_Martin ht_cold(
          p_in=EVA_EXH.node_c_in.p,
          At=EVA_EXH.At,
          beta=EVA_EXH.beta,
          qdot=EVA_EXH.plate.qdot_cold,
          qdot_tilde_start=EVA_EXH.mdot_hot_start*(EVA_EXH.h_hot_in_start - EVA_EXH.h_hot_out_start)
              /(2*EVA_EXH.l_start*EVA_EXH.w*EVA_EXH.N_ch_p)),
        redeclare Design.Materials.SS_AISI_304 material,
        thick=0.8e-3,
        b=2.3e-3,
        redeclare Design.Pressure_drops.Plates.evaporation_Martin dp_cold(p_in=
              EVA_EXH.node_c_in.p, beta=EVA_EXH.beta),
        redeclare Design.Miscellanea.check_velocity check_hot(T_sat=500),
        redeclare Design.Miscellanea.check_velocity check_cold(umin=min(EVA_EXH.ht_cold.u)),
        N_ch_p=16,
        mdot_hot_start=0.131,
        mdot_cold_start=0.032,
        use_dp=true,
        l_start=0.5,
        N_cell_pc=7,
        beta=1.0471975511966,
        t_hot_in_start=587.15,
        t_hot_out_start=473.15,
        t_cold_in_start=460.55,
        t_cold_out_start=553.15,
        p_hot_start=120000,
        p_cold_start=700000,
        dp_hot_start=20000)
        "Model of a flat plate heat exchanger. Evaporator for the exhaust gases."
        annotation (Placement(transformation(extent={{-88,-82},{88,60}})));

      Design.Components.Flags.ADDCO HOT_OUT(
        redeclare package Medium = Hot,
        use_T=false,
        T=473.15,
        p=100000,
        use_p=true)
        annotation (Placement(transformation(extent={{-133,-8},{-100,22}})));
      Design.Components.Flags.ADDCO COLD_IN(
        use_T=true,
        use_p=true,
        redeclare package Medium = Cold,
        use_m_flow=true,
        m_flow=0.0321154,
        T=460.55,
        p=700000)
        annotation (Placement(transformation(extent={{-133,-116},{-104,-88}})));
      Design.Components.Flags.ADDCO COLD_OUT(
        redeclare package Medium = Cold,
        T=553.15,
        use_T=true)
        annotation (Placement(transformation(extent={{79,-116},{108,-88}})));
      Design.Components.Flags.ADDCO HOT_IN(
        use_T=true,
        redeclare package Medium = Hot,
        use_m_flow=true,
        m_flow=0.131,
        T=587.15,
        p=120000,
        use_p=false)
                  annotation (Placement(transformation(extent={{49,100},{78,128}})));
    equation
      connect(HOT_OUT.node,EVA_EXH. node_h_out) annotation (Line(
          points={{-99.835,7},{-99.835,-53.6},{-88,-53.6}},
          color={0,0,0},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(COLD_IN.node,EVA_EXH. node_c_in) annotation (Line(
          points={{-103.855,-102},{-104,-102},{-104,-75.61},{-88,-75.61}},
          color={0,0,0},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(EVA_EXH.node_c_out, COLD_OUT.node) annotation (Line(
          points={{88,24.145},{130,24.145},{130,-102},{108.145,-102}},
          color={0,0,0},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(HOT_IN.node,EVA_EXH. node_h_in) annotation (Line(
          points={{78.145,114},{130,114},{130,45.8},{88,45.8}},
          color={0,0,0},
          pattern=LinePattern.None,
          smooth=Smooth.None));
     annotation (Placement(transformation(extent={{-108,-74},{88,66}})),
        experiment(__Dymola_NumberOfIntervals=1),
        __Dymola_experimentSetupOutput,
        Diagram(coordinateSystem(extent={{-160,-120},{140,140}},
              preserveAspectRatio=false), graphics),
        Icon(coordinateSystem(extent={{-160,-120},{140,140}})));
    end EVA_Exh;

    model Recuperator "Design of the recuperator"
      package Hot = Media.CoolProp.Toluene_TTSE;
      package Cold = Media.CoolProp.Toluene_TTSE;

        Design.Components.HEX.Flat_plate REC(
        redeclare package Medium_hot = Hot,
        redeclare package Medium_cold = Cold,
        redeclare function cost = Design.Miscellanea.Cost.FP_Rafferty,
        ht_hot_f1=1e9,
        ht_cold_f1=1e9,
        X=1,
        w=0.325,
        d_pt=100e-3,
        redeclare Design.Materials.SS_AISI_304 material,
        thick=0.8e-3,
        b=2.3e-3,
        redeclare Design.Miscellanea.check_velocity check_hot(T_sat=500),
        redeclare Design.Miscellanea.check_velocity check_cold(umin=min(REC.ht_cold.u)),
        N_ch_p=8,
        l_start=0.5,
        N_cell_pc=3,
        use_dp=true,
        mdot_hot_start=0.0607628,
        mdot_cold_start=0.0607628,
        beta=1.0471975511966,
        t_hot_in_start=503.15,
        t_hot_out_start=373.15,
        t_cold_in_start=358.15,
        t_cold_out_start=460.15,
        p_hot_start=47000,
        p_cold_start=700000)
        "Model of a flat plate heat exchanger. Recuperator."
        annotation (Placement(transformation(extent={{-88,-82},{88,60}})));

      Design.Components.Flags.ADDCO HOT_OUT(
        redeclare package Medium = Hot,
        T=373.15,
        use_T=true)
        annotation (Placement(transformation(extent={{-133,-8},{-100,22}})));
      Design.Components.Flags.ADDCO COLD_IN(
        use_T=true,
        use_p=true,
        redeclare package Medium = Cold,
        use_m_flow=true,
        T=358.5704,
        p=700000,
        m_flow=0.0607628)
        annotation (Placement(transformation(extent={{-133,-116},{-104,-88}})));
      Design.Components.Flags.ADDCO COLD_OUT(
        redeclare package Medium = Cold,
        T=553.15,
        use_T=false)
        annotation (Placement(transformation(extent={{79,-116},{108,-88}})));
      Design.Components.Flags.ADDCO HOT_IN(
        use_T=true,
        use_p=true,
        redeclare package Medium = Hot,
        use_m_flow=true,
        m_flow=0.0607628,
        T=502.732,
        p=47000)  annotation (Placement(transformation(extent={{49,100},{78,128}})));
    equation
      connect(HOT_OUT.node, REC.node_h_out) annotation (Line(
          points={{-99.835,7},{-99.835,-53.6},{-88,-53.6}},
          color={0,0,0},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(COLD_IN.node, REC.node_c_in) annotation (Line(
          points={{-103.855,-102},{-104,-102},{-104,-75.61},{-88,-75.61}},
          color={0,0,0},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(REC.node_c_out, COLD_OUT.node) annotation (Line(
          points={{88,24.145},{130,24.145},{130,-102},{108.145,-102}},
          color={0,0,0},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(HOT_IN.node, REC.node_h_in) annotation (Line(
          points={{78.145,114},{130,114},{130,45.8},{88,45.8}},
          color={0,0,0},
          pattern=LinePattern.None,
          smooth=Smooth.None));
     annotation (Placement(transformation(extent={{-108,-74},{88,66}})),
        experiment(__Dymola_NumberOfIntervals=1),
        __Dymola_experimentSetupOutput,
        Diagram(coordinateSystem(extent={{-160,-120},{140,140}},
              preserveAspectRatio=false), graphics),
        Icon(coordinateSystem(extent={{-160,-120},{140,140}})));
    end Recuperator;

    model Condenser "Design of the recuperator"
      package Hot = Media.CoolProp.Toluene_TTSE;
      package Cold = Media.CoolProp.Water_TTSE;

        Design.Components.HEX.Flat_plate COND(
        redeclare package Medium_hot = Hot,
        redeclare package Medium_cold = Cold,
        redeclare function cost = Design.Miscellanea.Cost.FP_Rafferty,
        ht_hot_f1=1e9,
        ht_cold_f1=1e9,
        X=1,
        w=0.325,
        d_pt=100e-3,
        redeclare Design.Heat_transfer.Plates.condensation_Longo ht_hot(beta=COND.beta,
        p_in=COND.node_h_in.p),
        redeclare Design.Pressure_drops.Plates.evaporation_Martin dp_hot(
          p_in=COND.node_h_in.p,
          beta=COND.beta),
        redeclare Design.Materials.SS_AISI_304 material,
        thick=0.8e-3,
        b=2.3e-3,
        redeclare Design.Miscellanea.check_velocity check_hot(T_sat=500),
        redeclare Design.Miscellanea.check_velocity check_cold(umin=min(COND.ht_cold.u)),
        l_start=0.5,
        N_cell_pc=10,
        mdot_hot_start=0.0607628,
        mdot_cold_start=1.09,
        use_dp=false,
        N_ch_p=6,
        beta=1.0471975511966,
        t_hot_in_start=373.15,
        t_hot_out_start=358.15,
        t_cold_in_start=343.15,
        t_cold_out_start=353.15,
        p_hot_start=47000,
        p_cold_start=120000) "Model of a flat plate heat exchanger. Condenser."
        annotation (Placement(transformation(extent={{-88,-82},{88,60}})));

      Design.Components.Flags.ADDCO HOT_OUT(
        redeclare package Medium = Hot,
        T=353.15,
        h=-50078.8,
        use_T=false,
        use_h=true)
        annotation (Placement(transformation(extent={{-133,-8},{-100,22}})));
      Design.Components.Flags.ADDCO COLD_IN(
        use_T=true,
        use_p=true,
        redeclare package Medium = Cold,
        use_m_flow=true,
        T=343.15,
        p=120000,
        m_flow=1.09)
        annotation (Placement(transformation(extent={{-133,-116},{-104,-88}})));
      Design.Components.Flags.ADDCO COLD_OUT(
        redeclare package Medium = Cold,
        T=553.15,
        use_T=false)
        annotation (Placement(transformation(extent={{79,-116},{108,-88}})));
      Design.Components.Flags.ADDCO HOT_IN(
        use_T=true,
        use_p=true,
        redeclare package Medium = Hot,
        use_m_flow=true,
        m_flow=0.0607628,
        T=373.15,
        p=47000)  annotation (Placement(transformation(extent={{49,100},{78,128}})));
    equation
      connect(HOT_OUT.node, COND.node_h_out) annotation (Line(
          points={{-99.835,7},{-99.835,-53.6},{-88,-53.6}},
          color={0,0,0},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(COLD_IN.node, COND.node_c_in) annotation (Line(
          points={{-103.855,-102},{-104,-102},{-104,-75.61},{-88,-75.61}},
          color={0,0,0},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(COND.node_c_out, COLD_OUT.node) annotation (Line(
          points={{88,24.145},{130,24.145},{130,-102},{108.145,-102}},
          color={0,0,0},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(HOT_IN.node, COND.node_h_in) annotation (Line(
          points={{78.145,114},{130,114},{130,45.8},{88,45.8}},
          color={0,0,0},
          pattern=LinePattern.None,
          smooth=Smooth.None));
     annotation (Placement(transformation(extent={{-108,-74},{88,66}})),
        experiment(__Dymola_NumberOfIntervals=1),
        __Dymola_experimentSetupOutput,
        Diagram(coordinateSystem(extent={{-160,-120},{140,140}},
              preserveAspectRatio=false), graphics),
        Icon(coordinateSystem(extent={{-160,-120},{140,140}})));
    end Condenser;
  end WHR_TRUCK;
  annotation (uses(Modelica(version="3.2.1")));
end Test;
