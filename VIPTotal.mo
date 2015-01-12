package ModelicaServices
  "(version = 3.2.1, target = \"Dymola\") Models and functions used in the Modelica Standard Library requiring a tool specific implementation"

package Machine

  final constant Real eps=1.e-15 "Biggest number such that 1.0 + eps = 1.0";

  final constant Real small=1.e-60
  "Smallest number such that small and -small are representable on the machine";
  annotation (Documentation(info="<html>
<p>
Package in which processor specific constants are defined that are needed
by numerical algorithms. Typically these constants are not directly used,
but indirectly via the alias definition in
<a href=\"modelica://Modelica.Constants\">Modelica.Constants</a>.
</p>
</html>"));
end Machine;
annotation (
  Protection(access=Access.hide),
  preferredView="info",
  version="3.2.1",
  versionDate="2013-01-17",
  versionBuild=1,
  uses(Modelica(version="3.2.1")),
  conversion(
    noneFromVersion="1.0",
    noneFromVersion="1.1",
    noneFromVersion="1.2"),
  Documentation(info="<html>
<p>
This package contains a set of functions and models to be used in the
Modelica Standard Library that requires a tool specific implementation.
These are:
</p>

<ul>
<li> <a href=\"modelica://ModelicaServices.Animation.Shape\">Shape</a>
     provides a 3-dim. visualization of elementary
     mechanical objects. It is used in
<a href=\"modelica://Modelica.Mechanics.MultiBody.Visualizers.Advanced.Shape\">Modelica.Mechanics.MultiBody.Visualizers.Advanced.Shape</a>
     via inheritance.</li>

<li> <a href=\"modelica://ModelicaServices.Animation.Surface\">Surface</a>
     provides a 3-dim. visualization of
     moveable parameterized surface. It is used in
<a href=\"modelica://Modelica.Mechanics.MultiBody.Visualizers.Advanced.Surface\">Modelica.Mechanics.MultiBody.Visualizers.Advanced.Surface</a>
     via inheritance.</li>

<li> <a href=\"modelica://ModelicaServices.ExternalReferences.loadResource\">loadResource</a>
     provides a function to return the absolute path name of an URI or a local file name. It is used in
<a href=\"modelica://Modelica.Utilities.Files.loadResource\">Modelica.Utilities.Files.loadResource</a>
     via inheritance.</li>

<li> <a href=\"modelica://ModelicaServices.Machine\">ModelicaServices.Machine</a>
     provides a package of machine constants. It is used in
<a href=\"modelica://Modelica.Constants\">Modelica.Constants</a>.</li>

<li> <a href=\"modelica://ModelicaServices.Types.SolverMethod\">Types.SolverMethod</a>
     provides a string defining the integration method to solve differential equations in
     a clocked discretized continuous-time partition (see Modelica 3.3 language specification).
     It is not yet used in the Modelica Standard Library, but in the Modelica_Synchronous library
     that provides convenience blocks for the clock operators of Modelica version &ge; 3.3.</li>
</ul>

<p>
This implementation is targeted for Dymola.
</p>

<p>
<b>Licensed by DLR and Dassault Syst&egrave;mes AB under the Modelica License 2</b><br>
Copyright &copy; 2009-2013, DLR and Dassault Syst&egrave;mes AB.
</p>

<p>
<i>This Modelica package is <u>free</u> software and the use is completely at <u>your own risk</u>; it can be redistributed and/or modified under the terms of the Modelica License 2. For license conditions (including the disclaimer of warranty) see <a href=\"modelica://Modelica.UsersGuide.ModelicaLicense2\">Modelica.UsersGuide.ModelicaLicense2</a> or visit <a href=\"http://www.modelica.org/licenses/ModelicaLicense2\"> http://www.modelica.org/licenses/ModelicaLicense2</a>.</i>
</p>

</html>"));
end ModelicaServices;

package ExternalMedia
  extends Modelica.Icons.Package;
  import SI = Modelica.SIunits;

  package Common "Package with common definitions"
    extends Modelica.Icons.Package;

    type InputChoice = enumeration(
      dT "(d,T) as inputs",
      hs "(h,s) as inputs",
      ph "(p,h) as inputs",
      ps "(p,s) as inputs",
      pT "(p,T) as inputs");

    type InputChoiceIncompressible = enumeration(
      ph "(p,h) as inputs",
      pT "(p,T) as inputs",
      phX "(p,h,X) as inputs",
      pTX "(p,T,X) as inputs");

    function XtoName "A function to convert concentration to substance name"
      extends Modelica.Icons.Function;
      input String substanceName = "";
      input Real[:] composition = {0.0};
      input String delimiter = "|";
      input Boolean debug = false;
      output String result;
  protected
      Integer nextIndex;
      Integer inLength;
      String name;
      String rest;
    algorithm
      if noEvent(size(composition,1) <= 0) then
        assert(not debug, "You are passing an empty composition vector, returning name only: "+substanceName, level=  AssertionLevel.warning);
        result :=substanceName;
      else
        assert(noEvent(size(composition,1)==1), "Your mixture has more than two components, ignoring all but the first element.", level=  AssertionLevel.warning);
        inLength  := Modelica.Utilities.Strings.length(substanceName);
        nextIndex := Modelica.Utilities.Strings.find(substanceName, delimiter);
        if noEvent(nextIndex<2) then
          // Assuming there are no special options
          name   := substanceName;
          rest   := "";
        else
          name   := Modelica.Utilities.Strings.substring(substanceName, 1, nextIndex-1);
          rest   := Modelica.Utilities.Strings.substring(substanceName, nextIndex, inLength);
        end if;
        if noEvent(noEvent(composition[1]<=0) or noEvent(composition[1]>=1)) then
          result := substanceName;
        else
          result := name + "-" + String(composition[1]) + rest;
        end if;
      end if;
      if noEvent(debug) then
        Modelica.Utilities.Streams.print(result+" --- "+substanceName);
      end if;
    end XtoName;

    function CheckCoolPropOptions
    "A function to extract and check the options passed to CoolProp"
      extends Modelica.Icons.Function;
      input String substance = "";
      input Boolean debug = false;
      output String result;

  protected
      Integer nextIndex;
      Integer intVal;
      Integer length;
      String name;
      String rest;
      // used to process the option
      String option;
      String value;
      // gather all valid options
      String options;
      // accept these inputs and set the default parameters
      String[:] allowedOptions = {
        "calc_transport",
        "enable_TTSE",
        "enable_BICUBIC",
        "enable_EXTTP",
        "twophase_derivsmoothing_xend",
        "rho_smoothing_xend",
        "debug"};
      String[:] defaultOptions = {
        "1",
        "0",
        "0",
        "1",
        "0.0",
        "0.0",
        "0"};
      // predefined delimiters
      String delimiter1 = "|";
      String delimiter2 = "=";

    algorithm
      if noEvent(debug) then
        Modelica.Utilities.Streams.print("input  = " + substance);
      end if;

      name := substance;

      for i in 1:size(allowedOptions,1) loop
        nextIndex := Modelica.Utilities.Strings.find(name, allowedOptions[i]);     // 0 if not found
        if nextIndex==0 then // not found
          name := name+delimiter1+allowedOptions[i]+delimiter2+defaultOptions[i];
        end if;
      end for;

      nextIndex := Modelica.Utilities.Strings.find(name, delimiter1);     // 0 if not found
      if nextIndex > 0 then
        // separate fluid name and options
        length  := Modelica.Utilities.Strings.length(name);
        rest    := Modelica.Utilities.Strings.substring(name, nextIndex+1, length);
        name    := Modelica.Utilities.Strings.substring(name, 1, nextIndex-1);
        options := "";

        while (nextIndex > 0) loop
          nextIndex := Modelica.Utilities.Strings.find(rest, delimiter1);     // 0 if not found
          if nextIndex > 0 then
            option  := Modelica.Utilities.Strings.substring(rest, 1, nextIndex-1);
            length  := Modelica.Utilities.Strings.length(rest);
            rest    := Modelica.Utilities.Strings.substring(rest, nextIndex+1, length);
          else
            option  := rest;
          end if;
          // now option contains enable_TTSE=1 or enable_TTSE
          intVal    := Modelica.Utilities.Strings.find(option, delimiter2);     // 0 if not found
          if intVal > 0 then // found "="
            length  := Modelica.Utilities.Strings.length(option);
            value   := Modelica.Utilities.Strings.substring(option, intVal+1, length);
            option  := Modelica.Utilities.Strings.substring(option, 1, intVal-1);
          else  // enable option by default
            value   := "1";
          end if;
          // now option contains only enable_TTSE
          intVal :=1;
          for i in 1:size(allowedOptions,1) loop
            if Modelica.Utilities.Strings.compare(option,allowedOptions[i])==Modelica.Utilities.Types.Compare.Equal then
              intVal := intVal - 1;
            end if;
          end for;
          if intVal <> 0 then
            assert(false, "Your option (" + option + ") is unknown.");
          else
            options := options+delimiter1+option+delimiter2+value;
          end if;
        end while;
      else
        // Assuming there are no special options
        name   := substance;
        options:= "";
      end if;

      result := name+options;
      if noEvent(debug) then
        Modelica.Utilities.Streams.print("output = " + result);
      end if;
    end CheckCoolPropOptions;
  end Common;

  package Media "Medium packages compatible with Modelica.Media"
    extends Modelica.Icons.Package;

    package FluidPropMedium "Medium package accessing the FluidProp solver"
      extends BaseClasses.ExternalTwoPhaseMedium;
      redeclare replaceable function setBubbleState
      "Set the thermodynamic state on the bubble line"
        extends Modelica.Icons.Function;
        input SaturationProperties sat "saturation point";
        input FixedPhase phase = 0 "phase flag";
        output ThermodynamicState state "complete thermodynamic state info";
        // Standard definition
        external "C" TwoPhaseMedium_setBubbleState_C_impl(sat, phase, state, mediumName, libraryName, substanceName)
          annotation(Include="#include \"externalmedialib.h\"", Library="ExternalMediaLib", IncludeDirectory="modelica://ExternalMedia/Resources/Include", LibraryDirectory="modelica://ExternalMedia/Resources/Library");
        annotation(Inline = true);
      end setBubbleState;

      redeclare replaceable function setDewState
      "Set the thermodynamic state on the dew line"
        extends Modelica.Icons.Function;
        input SaturationProperties sat "saturation point";
        input FixedPhase phase = 0 "phase flag";
        output ThermodynamicState state "complete thermodynamic state info";
        // Standard definition
        external "C" TwoPhaseMedium_setDewState_C_impl(sat, phase, state, mediumName, libraryName, substanceName)
          annotation(Include="#include \"externalmedialib.h\"", Library="ExternalMediaLib", IncludeDirectory="modelica://ExternalMedia/Resources/Include", LibraryDirectory="modelica://ExternalMedia/Resources/Library");
        annotation(Inline = true);
      end setDewState;

      redeclare function bubbleEntropy "Return bubble point specific entropy"
        input SaturationProperties sat "saturation property record";
        output SI.SpecificEntropy sl "boiling curve specific entropy";
      algorithm
        sl := specificEntropy(setBubbleState(sat));
      end bubbleEntropy;

      redeclare function dewEntropy "Return dew point specific entropy"
        input SaturationProperties sat "saturation property record";
        output SI.SpecificEntropy sv "dew curve specific entropy";
      algorithm
        sv := specificEntropy(setDewState(sat));
      end dewEntropy;

      redeclare function surfaceTension
        extends Modelica.Icons.Function;
        input SaturationProperties sat "saturation property record";
        output SurfaceTension sigma
        "Surface tension sigma in the two phase region";
      algorithm
        assert(false, "The FluidProp solver does not provide surface tension");
      end surfaceTension;
    end FluidPropMedium;

    package CoolPropMedium "Medium package accessing the CoolProp solver"
      extends BaseClasses.ExternalTwoPhaseMedium(
        final libraryName = "CoolProp",
        final substanceName = ExternalMedia.Common.CheckCoolPropOptions(substanceNames[1],debug=false));

      redeclare replaceable function isentropicEnthalpy
        input AbsolutePressure p_downstream "downstream pressure";
        input ThermodynamicState refState "reference state for entropy";
        output SpecificEnthalpy h_is "Isentropic enthalpy";
    protected
        SpecificEntropy s_ideal;
        ThermodynamicState state_ideal;
      algorithm
        s_ideal := specificEntropy(refState);
        state_ideal := setState_psX(p_downstream, s_ideal);
        h_is := specificEnthalpy(state_ideal);
      end isentropicEnthalpy;

      redeclare replaceable function setBubbleState
      "Set the thermodynamic state on the bubble line"
        extends Modelica.Icons.Function;
        input SaturationProperties sat "saturation point";
        input FixedPhase phase=0 "phase flag";
        output ThermodynamicState state "complete thermodynamic state info";
        // Standard definition
      external"C" TwoPhaseMedium_setBubbleState_C_impl(
            sat,
            phase,
            state,
            mediumName,
            libraryName,
            substanceName)
          annotation(Include="#include \"externalmedialib.h\"", Library="ExternalMediaLib", IncludeDirectory="modelica://ExternalMedia/Resources/Include", LibraryDirectory="modelica://ExternalMedia/Resources/Library");
        annotation (Inline=true);
      end setBubbleState;

      redeclare replaceable function setDewState
      "Set the thermodynamic state on the dew line"
        extends Modelica.Icons.Function;
        input SaturationProperties sat "saturation point";
        input FixedPhase phase=0 "phase flag";
        output ThermodynamicState state "complete thermodynamic state info";
        // Standard definition
      external"C" TwoPhaseMedium_setDewState_C_impl(
            sat,
            phase,
            state,
            mediumName,
            libraryName,
            substanceName)
          annotation(Include="#include \"externalmedialib.h\"", Library="ExternalMediaLib", IncludeDirectory="modelica://ExternalMedia/Resources/Include", LibraryDirectory="modelica://ExternalMedia/Resources/Library");
        annotation (Inline=true);
      end setDewState;

      redeclare function bubbleEntropy "Return bubble point specific entropy"
        input SaturationProperties sat "saturation property record";
        output SI.SpecificEntropy sl "boiling curve specific entropy";
      algorithm
        sl := specificEntropy(setBubbleState(sat));
      end bubbleEntropy;

      redeclare function dewEntropy "Return dew point specific entropy"
        input SaturationProperties sat "saturation property record";
        output SI.SpecificEntropy sv "dew curve specific entropy";
      algorithm
        sv := specificEntropy(setDewState(sat));
      end dewEntropy;

      redeclare function surfaceTension
        extends Modelica.Icons.Function;
        input SaturationProperties sat "saturation property record";
        output SurfaceTension sigma
        "Surface tension sigma in the two phase region";
      algorithm
        assert(false, "The CoolProp solver does not provide surface tension");
      end surfaceTension;

    end CoolPropMedium;

    partial package IncompressibleCoolPropMedium
    "External incompressible medium with up to two components using CoolProp"
      extends Modelica.Media.Interfaces.PartialMedium(
        mediumName =  "ExternalMedium",
        singleState = true,
        reducedX =    true);
      import ExternalMedia.Common.InputChoiceIncompressible;
      constant String libraryName = "CoolProp"
      "Name of the external fluid property computation library";
      constant String substanceName = ExternalMedia.Common.CheckCoolPropOptions(substanceNames[1],debug=false)
      "Only one substance can be specified, predefined mixture in CoolProp";
      redeclare record extends FluidConstants "external fluid constants"
        MolarMass molarMass "molecular mass";
        Temperature criticalTemperature "critical temperature";
        AbsolutePressure criticalPressure "critical pressure";
        MolarVolume criticalMolarVolume "critical molar Volume";
      end FluidConstants;
      constant InputChoiceIncompressible inputChoice=InputChoiceIncompressible.pTX
      "Default choice of input variables for property computations, incompressibles are in p,T";
      redeclare replaceable record ThermodynamicState =
      ExternalMedia.Media.BaseClasses.ExternalTwoPhaseMedium.ThermodynamicState;

      redeclare replaceable model extends BaseProperties(
        p(stateSelect = if preferredMediumStates and
                           (basePropertiesInputChoice == InputChoiceIncompressible.phX or
                            basePropertiesInputChoice == InputChoiceIncompressible.pTX or
                            basePropertiesInputChoice == InputChoiceIncompressible.ph or
                            basePropertiesInputChoice == InputChoiceIncompressible.pT) then
                                StateSelect.prefer else StateSelect.default),
        T(stateSelect = if preferredMediumStates and
                           (basePropertiesInputChoice == InputChoiceIncompressible.pTX or
                            basePropertiesInputChoice == InputChoiceIncompressible.pT) then
                             StateSelect.prefer else StateSelect.default),
        h(stateSelect = if preferredMediumStates and
                           (basePropertiesInputChoice == InputChoiceIncompressible.phX or
                            basePropertiesInputChoice == InputChoiceIncompressible.ph) then
                             StateSelect.prefer else StateSelect.default))
        import ExternalMedia.Common.InputChoiceIncompressible;
        parameter InputChoiceIncompressible basePropertiesInputChoice=inputChoice
        "Choice of input variables for property computations";
        Integer phaseInput
        "Phase input for property computation functions, 2 for two-phase, 1 for one-phase, 0 if not known";
        Integer phaseOutput
        "Phase output for medium, 2 for two-phase, 1 for one-phase";
        SpecificEntropy s "Specific entropy";
        //SaturationProperties sat "saturation property record";
      equation
        phaseInput = 1 "Force one-phase property computation";
        R  = Modelica.Constants.small "Gas constant (of mixture if applicable)";
        MM = 0.001 "Molar mass (of mixture or single fluid)";
        if (basePropertiesInputChoice == InputChoiceIncompressible.phX or
            basePropertiesInputChoice == InputChoiceIncompressible.ph) then
          state = setState_phX(p, h, Xi, phaseInput);
          d = density_phX(p, h, Xi, phaseInput);
          s = specificEntropy_phX(p, h, Xi, phaseInput);
          T = temperature_phX(p, h, Xi, phaseInput);
        elseif (basePropertiesInputChoice == InputChoiceIncompressible.pTX or
                basePropertiesInputChoice == InputChoiceIncompressible.pT) then
          state = setState_pTX(p, T, Xi, phaseInput);
          d = density_pTX(p, T, Xi, phaseInput);
          h = specificEnthalpy_pTX(p, T, Xi, phaseInput);
          s = specificEntropy_pTX(p, T, Xi, phaseInput);
        end if;
        // Compute the internal energy
        u = h - p/d;
        // Compute the saturation properties record
        //sat = setSat_p_state(state);
        // No phase boundary crossing
        phaseOutput = 1;
      end BaseProperties;

      replaceable function setState_ph
      "Return thermodynamic state record from p and h"
        extends Modelica.Icons.Function;
        input AbsolutePressure p "pressure";
        input SpecificEnthalpy h "specific enthalpy";
        input Integer phase = 1
        "2 for two-phase, 1 for one-phase, 0 if not known";
        output ThermodynamicState state;
    protected
        String name;
      algorithm
      state := setState_ph_library(p, h, phase, substanceName);
      end setState_ph;

      redeclare replaceable function setState_phX
      "Return thermodynamic state record from p and h"
        extends Modelica.Icons.Function;
        input AbsolutePressure p "pressure";
        input SpecificEnthalpy h "specific enthalpy";
        input MassFraction X[nX] "Mass fractions";
        input Integer phase = 1
        "2 for two-phase, 1 for one-phase, 0 if not known";
        output ThermodynamicState state;
    protected
        String name;
      algorithm
      name := ExternalMedia.Common.XtoName(substanceName,X);
      state := setState_ph_library(p, h, phase, name);
      end setState_phX;

      function setState_ph_library
      "Return thermodynamic state record from p and h"
        extends Modelica.Icons.Function;
        input AbsolutePressure p "pressure";
        input SpecificEnthalpy h "specific enthalpy";
        input Integer phase = 1
        "2 for two-phase, 1 for one-phase, 0 if not known";
        input String name "name and mass fractions";
        output ThermodynamicState state;
      external "C" TwoPhaseMedium_setState_ph_C_impl(p, h, phase, state, mediumName, libraryName, name)
        annotation(Include="#include \"externalmedialib.h\"", Library="ExternalMediaLib", IncludeDirectory="modelica://ExternalMedia/Resources/Include", LibraryDirectory="modelica://ExternalMedia/Resources/Library");
      end setState_ph_library;

      replaceable function setState_pT
      "Return thermodynamic state record from p and T"
        extends Modelica.Icons.Function;
        input AbsolutePressure p "pressure";
        input Temperature T "temperature";
        input Integer phase = 1
        "2 for two-phase, 1 for one-phase, 0 if not known";
        output ThermodynamicState state;
    protected
        String name;
      algorithm
      state := setState_pT_library(p, T, phase, substanceName);
      end setState_pT;

      redeclare replaceable function setState_pTX
      "Return thermodynamic state record from p and T"
        extends Modelica.Icons.Function;
        input AbsolutePressure p "pressure";
        input Temperature T "temperature";
        input MassFraction X[:] "Mass fractions";
        input Integer phase = 1
        "2 for two-phase, 1 for one-phase, 0 if not known";
        output ThermodynamicState state;
    protected
        String name;
      algorithm
      name := ExternalMedia.Common.XtoName(substanceName,X);
      state := setState_pT_library(p, T, phase, name);
      end setState_pTX;

      function setState_pT_library
      "Return thermodynamic state record from p and T"
        extends Modelica.Icons.Function;
        input AbsolutePressure p "pressure";
        input Temperature T "temperature";
        input Integer phase = 1
        "2 for two-phase, 1 for one-phase, 0 if not known";
        input String name "name and mass fractions";
        output ThermodynamicState state;
      external "C" TwoPhaseMedium_setState_pT_C_impl(p, T, state, mediumName, libraryName, name)
        annotation(Include="#include \"externalmedialib.h\"", Library="ExternalMediaLib", IncludeDirectory="modelica://ExternalMedia/Resources/Include", LibraryDirectory="modelica://ExternalMedia/Resources/Library");
      end setState_pT_library;

      redeclare replaceable function setState_psX
      "Return thermodynamic state record from p and s"
        extends Modelica.Icons.Function;
        input AbsolutePressure p "pressure";
        input SpecificEntropy s "specific entropy";
        input MassFraction X[nX] "Mass fractions";
        input Integer phase = 1
        "2 for two-phase, 1 for one-phase, 0 if not known";
        output ThermodynamicState state;
    protected
        String in1 = ExternalMedia.Common.XtoName(substanceName,X);
        //assert(false, "Incompressibles only support pT and ph as inputs!", level=AssertionLevel.error);
      external "C" TwoPhaseMedium_setState_ps_C_impl(p, s, phase, state, mediumName, libraryName, in1)
        annotation(Include="#include \"externalmedialib.h\"", Library="ExternalMediaLib", IncludeDirectory="modelica://ExternalMedia/Resources/Include", LibraryDirectory="modelica://ExternalMedia/Resources/Library");
      end setState_psX;

      redeclare function density_phX "returns density for given p and h"
        extends Modelica.Icons.Function;
        input AbsolutePressure p "Pressure";
        input SpecificEnthalpy h "Enthalpy";
        input MassFraction X[nX] "Mass fractions";
        input Integer phase=1
        "2 for two-phase, 1 for one-phase, 0 if not known";
      //input ThermodynamicState state;
        output Density d "density";
      algorithm
        d := density_phX_state(p=p, h=h, X=X, state=setState_phX(p=p, h=h, X=X, phase=phase));
      annotation (
        Inline=true);
      end density_phX;

      function density_phX_state "returns density for given p and h"
        extends Modelica.Icons.Function;
        input AbsolutePressure p "Pressure";
        input SpecificEnthalpy h "Enthalpy";
        input MassFraction X[nX] "Mass fractions";
        input ThermodynamicState state;
        output Density d "density";
      algorithm
        d := density(state);
      annotation (
        Inline=false,
        LateInline=true,
        derivative(noDerivative=state,noDerivative=X)=density_phX_der);
      end density_phX_state;

      replaceable function density_phX_der "Total derivative of density_ph"
        extends Modelica.Icons.Function;
        input AbsolutePressure p "Pressure";
        input SpecificEnthalpy h "Specific enthalpy";
        input MassFraction X[nX] "Mass fractions";
        input ThermodynamicState state;
        input Real p_der "time derivative of pressure";
        input Real h_der "time derivative of specific enthalpy";
        output Real d_der "time derivative of density";
      algorithm
        d_der := p_der*density_derp_h(state=state)
               + h_der*density_derh_p(state=state);
      annotation (Inline=true);
      end density_phX_der;

      redeclare replaceable function extends density_derh_p
      "Return derivative of density wrt enthalpy at constant pressure from state"
        // Standard definition
      algorithm
        ddhp := state.ddhp;
        annotation(Inline = true);
      end density_derh_p;

      redeclare replaceable function extends density_derp_h
      "Return derivative of density wrt pressure at constant enthalpy from state"
        // Standard definition
      algorithm
        ddph := state.ddph;
        annotation(Inline = true);
      end density_derp_h;

      redeclare function temperature_phX
      "returns temperature for given p and h"
        extends Modelica.Icons.Function;
        input AbsolutePressure p "Pressure";
        input SpecificEnthalpy h "Enthalpy";
        input MassFraction X[nX] "Mass fractions";
        input Integer phase=1
        "2 for two-phase, 1 for one-phase, 0 if not known";
        output Temperature T "Temperature";
      algorithm
        T := temperature_phX_state(p=p, h=h, X=X, state=setState_phX(p=p, h=h, X=X, phase=phase));
      annotation (
        Inline=true,
        inverse(h=specificEnthalpy_pTX(p=p, T=T, X=X, phase=phase)));
      end temperature_phX;

      function temperature_phX_state "returns temperature for given p and h"
        extends Modelica.Icons.Function;
        input AbsolutePressure p "Pressure";
        input SpecificEnthalpy h "Enthalpy";
        input MassFraction X[nX] "Mass fractions";
        input ThermodynamicState state;
        output Temperature T "Temperature";
      algorithm
        T := temperature(state);
      annotation (
        Inline=false,
        LateInline=true,
        inverse(h=specificEnthalpy_pTX_state(p=p, T=T, X=X, state=state)));
      end temperature_phX_state;

        function specificEntropy_phX
      "returns specific entropy for a given p and h"
        extends Modelica.Icons.Function;
        input AbsolutePressure p "Pressure";
        input SpecificEnthalpy h "Specific Enthalpy";
        input MassFraction X[nX] "Mass fractions";
        input Integer phase=1
        "2 for two-phase, 1 for one-phase, 0 if not known";
        output SpecificEntropy s "Specific Entropy";
        algorithm
        s := specificEntropy_phX_state(p=p, h=h, X=X, state=setState_phX(p=p, h=h, X=X, phase=phase));
        annotation (
        Inline=true);
        end specificEntropy_phX;

      function specificEntropy_phX_state
      "returns specific entropy for a given p and h"
        extends Modelica.Icons.Function;
        input AbsolutePressure p "Pressure";
        input SpecificEnthalpy h "Specific Enthalpy";
        input MassFraction X[nX] "Mass fractions";
        input ThermodynamicState state;
        output SpecificEntropy s "Specific Entropy";
      algorithm
        s := specificEntropy(state);
      annotation (
        Inline=false,
        LateInline=true,
        derivative(noDerivative=state,noDerivative=X)=specificEntropy_phX_der);
      end specificEntropy_phX_state;

      function specificEntropy_phX_der "time derivative of specificEntropy_phX"
        extends Modelica.Icons.Function;
        input AbsolutePressure p;
        input SpecificEnthalpy h;
        input MassFraction X[nX] "Mass fractions";
        input ThermodynamicState state;
        input Real p_der "time derivative of pressure";
        input Real h_der "time derivative of specific enthalpy";
        output Real s_der "time derivative of specific entropy";
      algorithm
        s_der := p_der*(-1.0/(state.d*state.T))
               + h_der*( 1.0/state.T);
      annotation (
        Inline=true);
      end specificEntropy_phX_der;

      redeclare function density_pTX "Return density from p and T"
        extends Modelica.Icons.Function;
        input AbsolutePressure p "Pressure";
        input Temperature T "Temperature";
        input MassFraction X[nX] "Mass fractions";
        input Integer phase=1
        "2 for two-phase, 1 for one-phase, 0 if not known";
        output Density d "Density";
      algorithm
        d := density_pTX_state(p=p, T=T, X=X, state=setState_pTX(p=p, T=T, X=X, phase=phase));
      annotation (
        Inline=true);
      end density_pTX;

      function density_pTX_state
        extends Modelica.Icons.Function;
        input AbsolutePressure p "Pressure";
        input Temperature T "Temperature";
        input MassFraction X[nX] "Mass fractions";
        input ThermodynamicState state;
        output Density d "Density";
      algorithm
        d := density(state);
      annotation (
        Inline=false,
        LateInline=true);
      end density_pTX_state;

      redeclare function specificEnthalpy_pTX
      "returns specific enthalpy for given p and T"
        extends Modelica.Icons.Function;
        input AbsolutePressure p "Pressure";
        input Temperature T "Temperature";
        input MassFraction X[nX] "Mass fractions";
        input Integer phase=1
        "2 for two-phase, 1 for one-phase, 0 if not known";
        output SpecificEnthalpy h "specific enthalpy";
      algorithm
        h := specificEnthalpy_pTX_state(p=p, T=T, X=X, state=setState_pTX(p=p, T=T, X=X, phase=phase));
      annotation (
        Inline=true,
        inverse(T=temperature_phX(p=p, h=h, X=X, phase=phase)));
      end specificEnthalpy_pTX;

      function specificEnthalpy_pTX_state
      "returns specific enthalpy for given p and T"
        extends Modelica.Icons.Function;
        input AbsolutePressure p "Pressure";
        input Temperature T "Temperature";
        input MassFraction X[nX] "Mass fractions";
        input ThermodynamicState state;
        output SpecificEnthalpy h "specific enthalpy";
      algorithm
        h := specificEnthalpy(state);
      annotation (
        Inline=false,
        LateInline=true,
        inverse(T=temperature_phX_state(p=p, h=h, X=X, state=state)));
      end specificEnthalpy_pTX_state;

      redeclare function specificEntropy_pTX
      "returns specific entropy for a given p and T"
      extends Modelica.Icons.Function;
      input AbsolutePressure p "Pressure";
      input Temperature T "Temperature";
      input MassFraction X[nX] "Mass fractions";
      input Integer phase=1 "2 for two-phase, 1 for one-phase, 0 if not known";
      output SpecificEntropy s "Specific Entropy";
      algorithm
      s := specificEntropy_pTX_state(p=p, T=T, X=X, state=setState_pTX(p=p, T=T, X=X, phase=phase));
        annotation (
          Inline=true);
      end specificEntropy_pTX;

      function specificEntropy_pTX_state
      "returns specific entropy for a given p and T"
        extends Modelica.Icons.Function;
        input AbsolutePressure p "Pressure";
        input Temperature T "Temperature";
        input MassFraction X[nX] "Mass fractions";
        input ThermodynamicState state;
        output SpecificEntropy s "Specific Entropy";
      algorithm
        s := specificEntropy(state);
      annotation (
        Inline=false,
        LateInline=true);
      end specificEntropy_pTX_state;

      redeclare replaceable function extends density
      "Return density from state"
        // Standard definition
      algorithm
        d := state.d;
        annotation(Inline = true);
      end density;

      redeclare replaceable function extends pressure
      "Return pressure from state"
        // Standard definition
      algorithm
        p := state.p;
        /*  // If special definition in "C"
  external "C" p=  TwoPhaseMedium_pressure_(state, mediumName, libraryName, substanceName)
    annotation(Include="#include \"externalmedialib.h\"", Library="ExternalMediaLib", IncludeDirectory="modelica://ExternalMedia/Resources/Include", LibraryDirectory="modelica://ExternalMedia/Resources/Library");
*/
        annotation(Inline = true);
      end pressure;

      redeclare replaceable function extends specificEnthalpy
      "Return specific enthalpy from state"
        // Standard definition
      algorithm
        h := state.h;
        annotation(Inline = true);
      end specificEnthalpy;

      redeclare replaceable function extends specificEntropy
      "Return specific entropy from state"
        // Standard definition
      algorithm
        s := state.s;
        annotation(Inline = true);
      end specificEntropy;

      redeclare replaceable function extends temperature
      "Return temperature from state"
        // Standard definition
      algorithm
        T := state.T;
        annotation(Inline = true);
      end temperature;

      redeclare function extends prandtlNumber
        /*  // If special definition in "C"
  external "C" T=  TwoPhaseMedium_prandtlNumber_(state, mediumName, libraryName, substanceName)
    annotation(Include="#include \"externalmedialib.h\"", Library="ExternalMediaLib", IncludeDirectory="modelica://ExternalMedia/Resources/Include", LibraryDirectory="modelica://ExternalMedia/Resources/Library");
*/
        annotation(Inline = true);
      end prandtlNumber;

      redeclare replaceable function extends velocityOfSound
      "Return velocity of sound from state"
        // Standard definition
      algorithm
        a := state.a;
        /*  // If special definition in "C"
  external "C" a=  TwoPhaseMedium_velocityOfSound_(state, mediumName, libraryName, substanceName)
    annotation(Include="#include \"externalmedialib.h\"", Library="ExternalMediaLib", IncludeDirectory="modelica://ExternalMedia/Resources/Include", LibraryDirectory="modelica://ExternalMedia/Resources/Library");
*/
        annotation(Inline = true);
      end velocityOfSound;

      redeclare replaceable function extends specificHeatCapacityCp
      "Return specific heat capacity cp from state"
        // Standard definition
      algorithm
        cp := state.cp;
        /*  // If special definition in "C"
  external "C" cp=  TwoPhaseMedium_specificHeatCapacityCp_(state, mediumName, libraryName, substanceName)
    annotation(Include="#include \"externalmedialib.h\"", Library="ExternalMediaLib", IncludeDirectory="modelica://ExternalMedia/Resources/Include", LibraryDirectory="modelica://ExternalMedia/Resources/Library");
*/
        annotation(Inline = true);
      end specificHeatCapacityCp;

      redeclare replaceable function extends specificHeatCapacityCv
      "Return specific heat capacity cv from state"
        // Standard definition
      algorithm
        cv := state.cv;
        /*  // If special definition in "C"
  external "C" cv=  TwoPhaseMedium_specificHeatCapacityCv_(state, mediumName, libraryName, substanceName)
    annotation(Include="#include \"externalmedialib.h\"", Library="ExternalMediaLib", IncludeDirectory="modelica://ExternalMedia/Resources/Include", LibraryDirectory="modelica://ExternalMedia/Resources/Library");
*/
        annotation(Inline = true);
      end specificHeatCapacityCv;

      redeclare replaceable function extends dynamicViscosity
      "Return dynamic viscosity from state"
        // Standard definition
      algorithm
        eta := state.eta;
        /*  // If special definition in "C"
  external "C" eta=  TwoPhaseMedium_dynamicViscosity_(state, mediumName, libraryName, substanceName)
    annotation(Include="#include \"externalmedialib.h\"", Library="ExternalMediaLib", IncludeDirectory="modelica://ExternalMedia/Resources/Include", LibraryDirectory="modelica://ExternalMedia/Resources/Library");
*/
        annotation(Inline = true);
      end dynamicViscosity;

      redeclare replaceable function extends thermalConductivity
      "Return thermal conductivity from state"
        // Standard definition
      algorithm
        lambda := state.lambda;
        /*  // If special definition in "C"
  external "C" lambda=  TwoPhaseMedium_thermalConductivity_(state, mediumName, libraryName, substanceName)
    annotation(Include="#include \"externalmedialib.h\"", Library="ExternalMediaLib", IncludeDirectory="modelica://ExternalMedia/Resources/Include", LibraryDirectory="modelica://ExternalMedia/Resources/Library");
*/
        annotation(Inline = true);
      end thermalConductivity;
    end IncompressibleCoolPropMedium;

    package BaseClasses "Base classes for external media packages"
      extends Modelica.Icons.BasesPackage;

      package ExternalTwoPhaseMedium
      "Generic external two phase medium package"
        extends Modelica.Media.Interfaces.PartialTwoPhaseMedium(
          singleState = false,
          onePhase = false,
          smoothModel = false,
          fluidConstants = {externalFluidConstants});
        import ExternalMedia.Common.InputChoice;
        // mediumName is declared here instead of in the extends clause
        // to break a circular dependency in redeclaration that OpenModelica
        // cannot yet handle
        constant String mediumName="unusablePartialMedium" "Name of the medium";
        constant String libraryName = "UnusableExternalMedium"
        "Name of the external fluid property computation library";
        constant String substanceName = substanceNames[1]
        "Only one substance can be specified";
        constant FluidConstants externalFluidConstants = FluidConstants(
          iupacName=  "unknown",
          casRegistryNumber=  "unknown",
          chemicalFormula=  "unknown",
          structureFormula=  "unknown",
          molarMass=  getMolarMass(),
          criticalTemperature=  getCriticalTemperature(),
          criticalPressure=  getCriticalPressure(),
          criticalMolarVolume=  getCriticalMolarVolume(),
          acentricFactor=  0,
          triplePointTemperature=  280.0,
          triplePointPressure=  500.0,
          meltingPoint=  280,
          normalBoilingPoint=  380.0,
          dipoleMoment=  2.0);

        constant InputChoice inputChoice=InputChoice.ph
        "Default choice of input variables for property computations";
        redeclare replaceable record ThermodynamicState
          // Fields in ASCII lexicographical order to work in Dymola
          Temperature T "temperature";
          VelocityOfSound a "velocity of sound";
          Modelica.SIunits.CubicExpansionCoefficient beta
          "isobaric expansion coefficient";
          SpecificHeatCapacity cp "specific heat capacity cp";
          SpecificHeatCapacity cv "specific heat capacity cv";
          Density d "density";
          DerDensityByEnthalpy ddhp
          "derivative of density wrt enthalpy at constant pressure";
          DerDensityByPressure ddph
          "derivative of density wrt pressure at constant enthalpy";
          DynamicViscosity eta "dynamic viscosity";
          SpecificEnthalpy h "specific enthalpy";
          Modelica.SIunits.Compressibility kappa "compressibility";
          ThermalConductivity lambda "thermal conductivity";
          AbsolutePressure p "pressure";
          FixedPhase phase(min=0, max=2)
          "phase flag: 2 for two-phase, 1 for one-phase";
          SpecificEntropy s "specific entropy";
        end ThermodynamicState;

        redeclare record SaturationProperties
          // Fields in ASCII lexicographical order to work in Dymola
          Temperature Tsat "saturation temperature";
          Real dTp "derivative of Ts wrt pressure";
          DerDensityByPressure ddldp "derivative of dls wrt pressure";
          DerDensityByPressure ddvdp "derivative of dvs wrt pressure";
          DerEnthalpyByPressure dhldp "derivative of hls wrt pressure";
          DerEnthalpyByPressure dhvdp "derivative of hvs wrt pressure";
          Density dl "density at bubble line (for pressure ps)";
          Density dv "density at dew line (for pressure ps)";
          SpecificEnthalpy hl
          "specific enthalpy at bubble line (for pressure ps)";
          SpecificEnthalpy hv "specific enthalpy at dew line (for pressure ps)";
          AbsolutePressure psat "saturation pressure";
          SurfaceTension sigma "surface tension";
          SpecificEntropy sl
          "specific entropy at bubble line (for pressure ps)";
          SpecificEntropy sv "specific entropy at dew line (for pressure ps)";
        end SaturationProperties;

        redeclare replaceable model extends BaseProperties(
          p(stateSelect = if preferredMediumStates and
                             (basePropertiesInputChoice == InputChoice.ph or
                              basePropertiesInputChoice == InputChoice.pT or
                              basePropertiesInputChoice == InputChoice.ps) then
                                  StateSelect.prefer else StateSelect.default),
          T(stateSelect = if preferredMediumStates and
                             (basePropertiesInputChoice == InputChoice.pT or
                              basePropertiesInputChoice == InputChoice.dT) then
                               StateSelect.prefer else StateSelect.default),
          h(stateSelect = if preferredMediumStates and
                             (basePropertiesInputChoice == InputChoice.hs or
                              basePropertiesInputChoice == InputChoice.ph) then
                               StateSelect.prefer else StateSelect.default),
          d(stateSelect = if preferredMediumStates and
                             basePropertiesInputChoice == InputChoice.dT then
                               StateSelect.prefer else StateSelect.default))
          import ExternalMedia.Common.InputChoice;
          parameter InputChoice basePropertiesInputChoice=inputChoice
          "Choice of input variables for property computations";
          FixedPhase phaseInput
          "Phase input for property computation functions, 2 for two-phase, 1 for one-phase, 0 if not known";
          Integer phaseOutput
          "Phase output for medium, 2 for two-phase, 1 for one-phase";
          SpecificEntropy s(
            stateSelect = if (basePropertiesInputChoice == InputChoice.hs or
                              basePropertiesInputChoice == InputChoice.ps) then
                             StateSelect.prefer else StateSelect.default)
          "Specific entropy";
          SaturationProperties sat "saturation property record";
        equation
          MM = externalFluidConstants.molarMass;
          R = Modelica.Constants.R/MM;
          if (onePhase or (basePropertiesInputChoice == InputChoice.pT)) then
            phaseInput = 1 "Force one-phase property computation";
          else
            phaseInput = 0 "Unknown phase";
          end if;
          if (basePropertiesInputChoice == InputChoice.ph) then
            // Compute the state record (including the unique ID)
            state = setState_ph(p, h, phaseInput);
            // Compute the remaining variables.
            // It is not possible to use the standard functions like
            // d = density(state), because differentiation for index
            // reduction and change of state variables would not be supported
            // density_ph(), which has an appropriate derivative annotation,
            // is used instead. The implementation of density_ph() uses
            // setState with the same inputs, so there's no actual overhead
            d = density_ph(p, h, phaseInput);
            s = specificEntropy_ph(p, h, phaseInput);
            T = temperature_ph(p, h, phaseInput);
          elseif (basePropertiesInputChoice == InputChoice.dT) then
            state = setState_dT(d, T, phaseInput);
            h = specificEnthalpy(state);
            p = pressure(state);
            s = specificEntropy(state);
          elseif (basePropertiesInputChoice == InputChoice.pT) then
            state = setState_pT(p, T, phaseInput);
            d = density(state);
            h = specificEnthalpy(state);
            s = specificEntropy(state);
          elseif (basePropertiesInputChoice == InputChoice.ps) then
            state = setState_ps(p, s, phaseInput);
            d = density(state);
            h = specificEnthalpy(state);
            T = temperature(state);
          elseif (basePropertiesInputChoice == InputChoice.hs) then
            state = setState_hs(h, s, phaseInput);
            d = density(state);
            p = pressure(state);
            T = temperature(state);
          end if;
          // Compute the internal energy
          u = h - p/d;
          // Compute the saturation properties record only if below critical point
          //sat = setSat_p(min(p,fluidConstants[1].criticalPressure));
          sat = setSat_p_state(state);
          // Event generation for phase boundary crossing
          if smoothModel then
            // No event generation
            phaseOutput = state.phase;
          else
            // Event generation at phase boundary crossing
            if basePropertiesInputChoice == InputChoice.ph then
              phaseOutput = if ((h > bubbleEnthalpy(sat) and h < dewEnthalpy(sat)) and
                                 p < fluidConstants[1].criticalPressure) then 2 else 1;
            elseif basePropertiesInputChoice == InputChoice.dT then
              phaseOutput = if  ((d < bubbleDensity(sat) and d > dewDensity(sat)) and
                                  T < fluidConstants[1].criticalTemperature) then 2 else 1;
            elseif basePropertiesInputChoice == InputChoice.ps then
              phaseOutput = if ((s > bubbleEntropy(sat) and s < dewEntropy(sat)) and
                                 p < fluidConstants[1].criticalPressure) then 2 else 1;
            elseif basePropertiesInputChoice == InputChoice.hs then
              phaseOutput = if ((s > bubbleEntropy(sat)  and s < dewEntropy(sat)) and
                                (h > bubbleEnthalpy(sat) and h < dewEnthalpy(sat))) then 2 else 1;
            elseif basePropertiesInputChoice == InputChoice.pT then
              phaseOutput = 1;
            else
              assert(false, "You are using an unsupported pair of inputs.");
            end if;
          end if;
        end BaseProperties;

        redeclare function molarMass "Return the molar mass of the medium"
            input ThermodynamicState state;
            output MolarMass MM "Mixture molar mass";
        algorithm
          MM := fluidConstants[1].molarMass;
        end molarMass;

        replaceable function getMolarMass
          output MolarMass MM "molar mass";
          external "C" MM = TwoPhaseMedium_getMolarMass_C_impl(mediumName, libraryName, substanceName)
            annotation(Include="#include \"externalmedialib.h\"", Library="ExternalMediaLib", IncludeDirectory="modelica://ExternalMedia/Resources/Include", LibraryDirectory="modelica://ExternalMedia/Resources/Library");
        end getMolarMass;

        replaceable function getCriticalTemperature
          output Temperature Tc "Critical temperature";
          external "C" Tc = TwoPhaseMedium_getCriticalTemperature_C_impl(mediumName, libraryName, substanceName)
            annotation(Include="#include \"externalmedialib.h\"", Library="ExternalMediaLib", IncludeDirectory="modelica://ExternalMedia/Resources/Include", LibraryDirectory="modelica://ExternalMedia/Resources/Library");
        end getCriticalTemperature;

        replaceable function getCriticalPressure
          output AbsolutePressure pc "Critical temperature";
          external "C" pc = TwoPhaseMedium_getCriticalPressure_C_impl(mediumName, libraryName, substanceName)
            annotation(Include="#include \"externalmedialib.h\"", Library="ExternalMediaLib", IncludeDirectory="modelica://ExternalMedia/Resources/Include", LibraryDirectory="modelica://ExternalMedia/Resources/Library");
        end getCriticalPressure;

        replaceable function getCriticalMolarVolume
          output MolarVolume vc "Critical molar volume";
          external "C" vc = TwoPhaseMedium_getCriticalMolarVolume_C_impl(mediumName, libraryName, substanceName)
            annotation(Include="#include \"externalmedialib.h\"", Library="ExternalMediaLib", IncludeDirectory="modelica://ExternalMedia/Resources/Include", LibraryDirectory="modelica://ExternalMedia/Resources/Library");
        end getCriticalMolarVolume;

        redeclare replaceable function setState_ph
        "Return thermodynamic state record from p and h"
          extends Modelica.Icons.Function;
          input AbsolutePressure p "pressure";
          input SpecificEnthalpy h "specific enthalpy";
          input FixedPhase phase = 0
          "2 for two-phase, 1 for one-phase, 0 if not known";
          output ThermodynamicState state;
        external "C" TwoPhaseMedium_setState_ph_C_impl(p, h, phase, state, mediumName, libraryName, substanceName)
          annotation(Include="#include \"externalmedialib.h\"", Library="ExternalMediaLib", IncludeDirectory="modelica://ExternalMedia/Resources/Include", LibraryDirectory="modelica://ExternalMedia/Resources/Library");
        end setState_ph;

        redeclare replaceable function setState_pT
        "Return thermodynamic state record from p and T"
          extends Modelica.Icons.Function;
          input AbsolutePressure p "pressure";
          input Temperature T "temperature";
          input FixedPhase phase = 0
          "2 for two-phase, 1 for one-phase, 0 if not known";
          output ThermodynamicState state;
        external "C" TwoPhaseMedium_setState_pT_C_impl(p, T, state, mediumName, libraryName, substanceName)
          annotation(Include="#include \"externalmedialib.h\"", Library="ExternalMediaLib", IncludeDirectory="modelica://ExternalMedia/Resources/Include", LibraryDirectory="modelica://ExternalMedia/Resources/Library");
        end setState_pT;

        redeclare replaceable function setState_dT
        "Return thermodynamic state record from d and T"
          extends Modelica.Icons.Function;
          input Density d "density";
          input Temperature T "temperature";
          input FixedPhase phase = 0
          "2 for two-phase, 1 for one-phase, 0 if not known";
          output ThermodynamicState state;
        external "C" TwoPhaseMedium_setState_dT_C_impl(d, T, phase, state, mediumName, libraryName, substanceName)
          annotation(Include="#include \"externalmedialib.h\"", Library="ExternalMediaLib", IncludeDirectory="modelica://ExternalMedia/Resources/Include", LibraryDirectory="modelica://ExternalMedia/Resources/Library");
        end setState_dT;

        redeclare replaceable function setState_ps
        "Return thermodynamic state record from p and s"
          extends Modelica.Icons.Function;
          input AbsolutePressure p "pressure";
          input SpecificEntropy s "specific entropy";
          input FixedPhase phase = 0
          "2 for two-phase, 1 for one-phase, 0 if not known";
          output ThermodynamicState state;
        external "C" TwoPhaseMedium_setState_ps_C_impl(p, s, phase, state, mediumName, libraryName, substanceName)
          annotation(Include="#include \"externalmedialib.h\"", Library="ExternalMediaLib", IncludeDirectory="modelica://ExternalMedia/Resources/Include", LibraryDirectory="modelica://ExternalMedia/Resources/Library");
        end setState_ps;

        replaceable function setState_hs
        "Return thermodynamic state record from h and s"
          extends Modelica.Icons.Function;
          input SpecificEnthalpy h "specific enthalpy";
          input SpecificEntropy s "specific entropy";
          input FixedPhase phase = 0
          "2 for two-phase, 1 for one-phase, 0 if not known";
          output ThermodynamicState state;
        external "C" TwoPhaseMedium_setState_hs_C_impl(h, s, phase, state, mediumName, libraryName, substanceName)
          annotation(Include="#include \"externalmedialib.h\"", Library="ExternalMediaLib", IncludeDirectory="modelica://ExternalMedia/Resources/Include", LibraryDirectory="modelica://ExternalMedia/Resources/Library");
        end setState_hs;

        redeclare function extends setState_phX
        algorithm
          // The composition is an empty vector
          state :=setState_ph(p, h, phase);
        end setState_phX;

        redeclare function extends setState_pTX
        algorithm
          // The composition is an empty vector
          state :=setState_pT(p, T, phase);
        end setState_pTX;

        redeclare function extends setState_dTX
        algorithm
          // The composition is an empty vector
          state :=setState_dT(d, T, phase);
        end setState_dTX;

        redeclare function extends setState_psX
        algorithm
          // The composition is an empty vector
          state :=setState_ps(p, s, phase);
        end setState_psX;

        replaceable function setState_hsX
                                          extends Modelica.Icons.Function;
          input SpecificEnthalpy h "specific enthalpy";
          input SpecificEntropy s "specific entropy";
          input FixedPhase phase = 0
          "2 for two-phase, 1 for one-phase, 0 if not known";
          output ThermodynamicState state;
        algorithm
          // The composition is an empty vector
          state :=setState_hs(h, s, phase);
        end setState_hsX;

        redeclare replaceable function density_ph "Return density from p and h"
          extends Modelica.Icons.Function;
          input AbsolutePressure p "Pressure";
          input SpecificEnthalpy h "Specific enthalpy";
          input FixedPhase phase = 0
          "2 for two-phase, 1 for one-phase, 0 if not known";
          output Density d "Density";
        algorithm
          d := density_ph_state(p=p, h=h, state=setState_ph(p=p, h=h, phase=phase));
        annotation (Inline = true);
        end density_ph;

        function density_ph_state "returns density for given p and h"
          extends Modelica.Icons.Function;
          input AbsolutePressure p "Pressure";
          input SpecificEnthalpy h "Enthalpy";
          input ThermodynamicState state;
          output Density d "density";
        algorithm
          d := density(state);
        annotation (
          Inline=false,
          LateInline=true,
          derivative(noDerivative=state)=density_ph_der);
        end density_ph_state;

        replaceable function density_ph_der "Total derivative of density_ph"
          extends Modelica.Icons.Function;
          input AbsolutePressure p "Pressure";
          input SpecificEnthalpy h "Specific enthalpy";
          input ThermodynamicState state;
          input Real p_der "time derivative of pressure";
          input Real h_der "time derivative of specific enthalpy";
          output Real d_der "time derivative of density";
        algorithm
          d_der := p_der*density_derp_h(state=state)
                 + h_der*density_derh_p(state=state);
        annotation (Inline=true);
        end density_ph_der;

        redeclare replaceable function temperature_ph
        "Return temperature from p and h"
          extends Modelica.Icons.Function;
          input AbsolutePressure p "Pressure";
          input SpecificEnthalpy h "Specific enthalpy";
          input FixedPhase phase = 0
          "2 for two-phase, 1 for one-phase, 0 if not known";
          output Temperature T "Temperature";
        algorithm
          T := temperature_ph_state(p=p, h=h, state=setState_ph(p=p, h=h, phase=phase));
        annotation (
          Inline=true,
          inverse(h=specificEnthalpy_pT(p=p, T=T, phase=phase)));
        end temperature_ph;

        function temperature_ph_state "returns temperature for given p and h"
          extends Modelica.Icons.Function;
          input AbsolutePressure p "Pressure";
          input SpecificEnthalpy h "Enthalpy";
          input ThermodynamicState state;
          output Temperature T "Temperature";
        algorithm
          T := temperature(state);
        annotation (
          Inline=false,
          LateInline=true,
          inverse(h=specificEnthalpy_pT_state(p=p, T=T, state=state)));
        end temperature_ph_state;

        replaceable function specificEntropy_ph
        "Return specific entropy from p and h"
          extends Modelica.Icons.Function;
          input AbsolutePressure p "Pressure";
          input SpecificEnthalpy h "Specific enthalpy";
          input FixedPhase phase = 0
          "2 for two-phase, 1 for one-phase, 0 if not known";
          output SpecificEntropy s "specific entropy";
        algorithm
          s := specificEntropy_ph_state(p=p, h=h, state=setState_ph(p=p, h=h, phase=phase));
          annotation (
          Inline=true,
          inverse(h=specificEnthalpy_ps(p=p, s=s, phase=phase)));
        end specificEntropy_ph;

        function specificEntropy_ph_state
        "returns specific entropy for a given p and h"
          extends Modelica.Icons.Function;
          input AbsolutePressure p "Pressure";
          input SpecificEnthalpy h "Specific Enthalpy";
          input ThermodynamicState state;
          output SpecificEntropy s "Specific Entropy";
        algorithm
          s := specificEntropy(state);
        annotation (
          Inline=false,
          LateInline=true,
          derivative(noDerivative=state)=specificEntropy_ph_der);
        end specificEntropy_ph_state;

        function specificEntropy_ph_der "time derivative of specificEntropy_ph"
          extends Modelica.Icons.Function;
          input AbsolutePressure p;
          input SpecificEnthalpy h;
          input ThermodynamicState state;
          input Real p_der "time derivative of pressure";
          input Real h_der "time derivative of specific enthalpy";
          output Real s_der "time derivative of specific entropy";
        algorithm
          s_der := p_der*(-1.0/(state.d*state.T))
                 + h_der*( 1.0/state.T);
        annotation (Inline = true);
        end specificEntropy_ph_der;

        redeclare replaceable function density_pT "Return density from p and T"
          extends Modelica.Icons.Function;
          input AbsolutePressure p "Pressure";
          input Temperature T "Temperature";
          input FixedPhase phase = 0
          "2 for two-phase, 1 for one-phase, 0 if not known";
          output Density d "Density";
        algorithm
          d := density_pT_state(p=p, T=T, state=setState_pT(p=p, T=T, phase=phase));
        annotation (
          Inline=true,
          inverse(p=pressure_dT(d=d, T=T, phase=phase)));
        end density_pT;

        function density_pT_state
          extends Modelica.Icons.Function;
          input AbsolutePressure p "Pressure";
          input Temperature T "Temperature";
          input ThermodynamicState state;
          output Density d "Density";
        algorithm
          d := density(state);
        annotation (
          Inline=false,
          LateInline=true,
          inverse(p=pressure_dT_state(d=d, T=T, state=state)));
        end density_pT_state;

        replaceable function density_pT_der "Total derivative of density_pT"
          extends Modelica.Icons.Function;
          input AbsolutePressure p "Pressure";
          input Temperature T "Temperature";
          input FixedPhase phase
          "2 for two-phase, 1 for one-phase, 0 if not known";
          input Real p_der;
          input Real T_der;
          output Real d_der;
        algorithm
          d_der:=density_derp_T(setState_pT(p, T))*p_der +
                 density_derT_p(setState_pT(p, T))*T_der;
          /*  // If special definition in "C"
    external "C" d_der=  TwoPhaseMedium_density_pT_der_C_impl(state, mediumName, libraryName, substanceName)
    annotation(Include="#include \"externalmedialib.h\"", Library="ExternalMediaLib", IncludeDirectory="modelica://ExternalMedia/Resources/Include", LibraryDirectory="modelica://ExternalMedia/Resources/Library");
    */
          annotation(Inline = true);
        end density_pT_der;

        redeclare replaceable function specificEnthalpy_pT
        "Return specific enthalpy from p and T"
          extends Modelica.Icons.Function;
          input AbsolutePressure p "Pressure";
          input Temperature T "Temperature";
          input FixedPhase phase = 0
          "2 for two-phase, 1 for one-phase, 0 if not known";
          output SpecificEnthalpy h "specific enthalpy";
        algorithm
          h := specificEnthalpy_pT_state(p=p, T=T, state=setState_pT(p=p, T=T, phase=phase));
        annotation (
          Inline=true,
          inverse(T=temperature_ph(p=p, h=h, phase=phase)));
        end specificEnthalpy_pT;

        function specificEnthalpy_pT_state
        "returns specific enthalpy for given p and T"
          extends Modelica.Icons.Function;
          input AbsolutePressure p "Pressure";
          input Temperature T "Temperature";
          input ThermodynamicState state;
          output SpecificEnthalpy h "specific enthalpy";
        algorithm
          h := specificEnthalpy(state);
        annotation (
          Inline=false,
          LateInline=true,
          inverse(T=temperature_ph_state(p=p, h=h, state=state)));
        end specificEnthalpy_pT_state;

        function specificEntropy_pT
        "returns specific entropy for a given p and T"
          extends Modelica.Icons.Function;
          input AbsolutePressure p "Pressure";
          input Temperature T "Temperature";
          input FixedPhase phase=0
          "2 for two-phase, 1 for one-phase, 0 if not known";
          output SpecificEntropy s "Specific Entropy";
        algorithm
          s := specificEntropy_pT_state(p=p, T=T, state=setState_pT(p=p, T=T, phase=phase));
        annotation (
          Inline=true,
          inverse(T=temperature_ps(p=p, s=s, phase=phase)));
        end specificEntropy_pT;

        function specificEntropy_pT_state
        "returns specific entropy for a given p and T"
          extends Modelica.Icons.Function;
          input AbsolutePressure p "Pressure";
          input Temperature T "Temperature";
          input ThermodynamicState state;
          output SpecificEntropy s "Specific Entropy";
        algorithm
          s := specificEntropy(state);
        annotation (
          Inline=false,
          LateInline=true);
        end specificEntropy_pT_state;

        redeclare replaceable function pressure_dT
        "Return pressure from d and T"
          extends Modelica.Icons.Function;
          input Density d "Density";
          input Temperature T "Temperature";
          input FixedPhase phase = 0
          "2 for two-phase, 1 for one-phase, 0 if not known";
          output AbsolutePressure p "Pressure";
        algorithm
          p := pressure_dT_state(d=d, T=T, state=setState_dT(d=d, T=T, phase=phase));
          annotation (
          Inline=true,
          inverse(d=density_pT(p=p, T=T, phase=phase)));
        end pressure_dT;

        function pressure_dT_state
          extends Modelica.Icons.Function;
          input Density d "Density";
          input Temperature T "Temperature";
          input ThermodynamicState state;
          output AbsolutePressure p "pressure";
        algorithm
          p := pressure(state);
        annotation (
          Inline=false,
          LateInline=true,
          inverse(d=density_pT_state(p=p, T=T, state=state)));
        end pressure_dT_state;

        redeclare replaceable function specificEnthalpy_dT
        "Return specific enthalpy from d and T"
          extends Modelica.Icons.Function;
          input Density d "Density";
          input Temperature T "Temperature";
          input FixedPhase phase = 0
          "2 for two-phase, 1 for one-phase, 0 if not known";
          output SpecificEnthalpy h "specific enthalpy";
        algorithm
          h := specificEnthalpy_dT_state(d=d, T=T, state=setState_dT(d=d, T=T, phase=phase));
        annotation (
          Inline=true);
        end specificEnthalpy_dT;

        function specificEnthalpy_dT_state
          extends Modelica.Icons.Function;
          input Density d "Density";
          input Temperature T "Temperature";
          input ThermodynamicState state;
          output SpecificEnthalpy h "SpecificEnthalpy";
        algorithm
          h := specificEnthalpy(state);
        annotation (
          Inline=false,
          LateInline=true);
        end specificEnthalpy_dT_state;

        function specificEntropy_dT
        "returns specific entropy for a given d and T"
          extends Modelica.Icons.Function;
          input Density d "Density";
          input Temperature T "Temperature";
          input FixedPhase phase=0
          "2 for two-phase, 1 for one-phase, 0 if not known";
        output SpecificEntropy s "Specific Entropy";
        algorithm
          s := specificEntropy_dT_state(d=d, T=T, state=setState_dT(d=d, T=T, phase=phase));
        annotation (Inline=true);
        end specificEntropy_dT;

        function specificEntropy_dT_state
        "returns specific entropy for a given d and T"
          extends Modelica.Icons.Function;
          input Density d "Density";
          input Temperature T "Temperature";
          input ThermodynamicState state;
          output SpecificEntropy s "Specific Entropy";
        algorithm
          s := specificEntropy(state);
        annotation (
          Inline=false,
          LateInline=true);
        end specificEntropy_dT_state;

        redeclare replaceable function density_ps "Return density from p and s"
          extends Modelica.Icons.Function;
          input AbsolutePressure p "Pressure";
          input SpecificEntropy s "Specific entropy";
          input FixedPhase phase = 0
          "2 for two-phase, 1 for one-phase, 0 if not known";
          output Density d "Density";
        algorithm
          d := density_ps_state(p=p, s=s, state=setState_ps(p=p, s=s, phase=phase));
        annotation (
          Inline=true);
        end density_ps;

        function density_ps_state "Return density from p and s"
          extends Modelica.Icons.Function;
          input AbsolutePressure p "Pressure";
          input SpecificEntropy s "Specific entropy";
          input ThermodynamicState state;
          output Density d "Density";
        algorithm
          d := density(state);
        annotation (
          Inline=false,
          LateInline=true,
          derivative(noDerivative=state) = density_ps_der);
        end density_ps_state;

        replaceable partial function density_ps_der
        "Total derivative of density_ps"
          extends Modelica.Icons.Function;
          input AbsolutePressure p "Pressure";
          input SpecificEntropy s "Specific entropy";
          input ThermodynamicState state;
          input Real p_der;
          input Real h_der;
          output Real d_der;
          // To be implemented
          annotation(Inline = true);
        end density_ps_der;

        redeclare replaceable function temperature_ps
        "Return temperature from p and s"
          extends Modelica.Icons.Function;
          input AbsolutePressure p "Pressure";
          input SpecificEntropy s "Specific entropy";
          input FixedPhase phase = 0
          "2 for two-phase, 1 for one-phase, 0 if not known";
          output Temperature T "Temperature";
        algorithm
          T := temperature_ps_state(p=p, s=s, state=setState_ps(p=p, s=s, phase=phase));
        annotation (
          Inline=true,
          inverse(s=specificEntropy_pT(p=p, T=T, phase=phase)));
        end temperature_ps;

        function temperature_ps_state "returns temperature for given p and s"
          extends Modelica.Icons.Function;
          input AbsolutePressure p "Pressure";
          input SpecificEntropy s "Specific entropy";
          input ThermodynamicState state;
          output Temperature T "Temperature";
        algorithm
          T := temperature(state);
        annotation (
          Inline=false,
          LateInline=true,
          inverse(s=specificEntropy_pT_state(p=p, T=T, state=state)));
        end temperature_ps_state;

        redeclare replaceable function specificEnthalpy_ps
        "Return specific enthalpy from p and s"
          extends Modelica.Icons.Function;
          input AbsolutePressure p "Pressure";
          input SpecificEntropy s "Specific entropy";
          input FixedPhase phase = 0
          "2 for two-phase, 1 for one-phase, 0 if not known";
          output SpecificEnthalpy h "specific enthalpy";
        algorithm
          h := specificEnthalpy_ps_state(p=p, s=s, state=setState_ps(p=p, s=s, phase=phase));
          annotation (
          Inline = true,
          inverse(s=specificEntropy_ph(p=p, h=h, phase=phase)));
        end specificEnthalpy_ps;

        function specificEnthalpy_ps_state "Return enthalpy from p and s"
          extends Modelica.Icons.Function;
          input AbsolutePressure p "Pressure";
          input SpecificEntropy s "Specific entropy";
          input ThermodynamicState state;
          output SpecificEnthalpy h "Enthalpy";
        algorithm
          h := specificEnthalpy(state);
        annotation (
          Inline=false,
          LateInline=true,
          inverse(s=specificEntropy_ph_state(p=p, h=h, state=state)));
        end specificEnthalpy_ps_state;

        function density_hs "Return density for given h and s"
          extends Modelica.Icons.Function;
          input SpecificEnthalpy h "Enthalpy";
          input SpecificEntropy s "Specific entropy";
          input FixedPhase phase=0
          "2 for two-phase, 1 for one-phase, 0 if not known";
          output Density d "density";
        algorithm
          d := density_hs_state(h=h, s=s, state=setState_hs(h=h, s=s, phase=phase));
        annotation (
          Inline=true);
        end density_hs;

        function density_hs_state "Return density for given h and s"
          extends Modelica.Icons.Function;
          input SpecificEnthalpy h "Enthalpy";
          input SpecificEntropy s "Specific entropy";
          input ThermodynamicState state;
          output Density d "density";
        algorithm
          d := density(state);
        annotation (
          Inline=false,
          LateInline=true);
        end density_hs_state;

        replaceable function pressure_hs "Return pressure from h and s"
          extends Modelica.Icons.Function;
          input SpecificEnthalpy h "Specific enthalpy";
          input SpecificEntropy s "Specific entropy";
          input FixedPhase phase = 0
          "2 for two-phase, 1 for one-phase, 0 if not known";
          output AbsolutePressure p "Pressure";
        algorithm
          p := pressure_hs_state(h=h, s=s, state=setState_hs(h=h, s=s, phase=phase));
          annotation (
            Inline = true,
            inverse(
              h=specificEnthalpy_ps(p=p, s=s, phase=phase),
              s=specificEntropy_ph(p=p, h=h, phase=phase)));
        end pressure_hs;

        function pressure_hs_state "Return pressure for given h and s"
          extends Modelica.Icons.Function;
          input SpecificEnthalpy h "Enthalpy";
          input SpecificEntropy s "Specific entropy";
          input ThermodynamicState state;
          output AbsolutePressure p "Pressure";
        algorithm
          p := pressure(state);
        annotation (
          Inline=false,
          LateInline=true);
        end pressure_hs_state;

        replaceable function temperature_hs "Return temperature from h and s"
          extends Modelica.Icons.Function;
          input SpecificEnthalpy h "Specific enthalpy";
          input SpecificEntropy s "Specific entropy";
          input FixedPhase phase = 0
          "2 for two-phase, 1 for one-phase, 0 if not known";
          output Temperature T "Temperature";
        algorithm
          T := temperature_hs_state(h=h, s=s, state=setState_hs(h=h, s=s, phase=phase));
          annotation (
            Inline = true);
        end temperature_hs;

        function temperature_hs_state "Return temperature for given h and s"
          extends Modelica.Icons.Function;
          input SpecificEnthalpy h "Enthalpy";
          input SpecificEntropy s "Specific entropy";
          input ThermodynamicState state;
          output Temperature T "Temperature";
        algorithm
          T := temperature(state);
        annotation (
          Inline=false,
          LateInline=true);
        end temperature_hs_state;

        redeclare function extends prandtlNumber "Returns Prandtl number"
          /*  // If special definition in "C"
  external "C" T=  TwoPhaseMedium_prandtlNumber_C_impl(state, mediumName, libraryName, substanceName)
    annotation(Include="#include \"externalmedialib.h\"", Library="ExternalMediaLib", IncludeDirectory="modelica://ExternalMedia/Resources/Include", LibraryDirectory="modelica://ExternalMedia/Resources/Library");
*/
          annotation(Inline = true);
        end prandtlNumber;

        redeclare replaceable function extends temperature
        "Return temperature from state"
          // Standard definition
        algorithm
          T := state.T;
          /*  // If special definition in "C"
  external "C" T=  TwoPhaseMedium_temperature_C_impl(state, mediumName, libraryName, substanceName)
    annotation(Include="#include \"externalmedialib.h\"", Library="ExternalMediaLib", IncludeDirectory="modelica://ExternalMedia/Resources/Include", LibraryDirectory="modelica://ExternalMedia/Resources/Library");
*/
          annotation(Inline = true);
        end temperature;

        redeclare replaceable function extends velocityOfSound
        "Return velocity of sound from state"
          // Standard definition
        algorithm
          a := state.a;
          /*  // If special definition in "C"
  external "C" a=  TwoPhaseMedium_velocityOfSound_C_impl(state, mediumName, libraryName, substanceName)
    annotation(Include="#include \"externalmedialib.h\"", Library="ExternalMediaLib", IncludeDirectory="modelica://ExternalMedia/Resources/Include", LibraryDirectory="modelica://ExternalMedia/Resources/Library");
*/
          annotation(Inline = true);
        end velocityOfSound;

        redeclare replaceable function extends isobaricExpansionCoefficient
        "Return isobaric expansion coefficient from state"
          // Standard definition
        algorithm
          beta := state.beta;
          /*  // If special definition in "C"
  external "C" beta=  TwoPhaseMedium_isobaricExpansionCoefficient_C_impl(state, mediumName, libraryName, substanceName)
    annotation(Include="#include \"externalmedialib.h\"", Library="ExternalMediaLib", IncludeDirectory="modelica://ExternalMedia/Resources/Include", LibraryDirectory="modelica://ExternalMedia/Resources/Library");
*/
          annotation(Inline = true);
        end isobaricExpansionCoefficient;

        redeclare replaceable function extends isentropicExponent
        "Return isentropic exponent"
          extends Modelica.Icons.Function;
        algorithm
          gamma := density(state) / pressure(state) * velocityOfSound(state) * velocityOfSound(state);
        end isentropicExponent;

        redeclare replaceable function extends specificHeatCapacityCp
        "Return specific heat capacity cp from state"
          // Standard definition
        algorithm
          cp := state.cp;
          /*  // If special definition in "C"
  external "C" cp=  TwoPhaseMedium_specificHeatCapacityCp_C_impl(state, mediumName, libraryName, substanceName)
    annotation(Include="#include \"externalmedialib.h\"", Library="ExternalMediaLib", IncludeDirectory="modelica://ExternalMedia/Resources/Include", LibraryDirectory="modelica://ExternalMedia/Resources/Library");
*/
          annotation(Inline = true);
        end specificHeatCapacityCp;

        redeclare replaceable function extends specificHeatCapacityCv
        "Return specific heat capacity cv from state"
          // Standard definition
        algorithm
          cv := state.cv;
          /*  // If special definition in "C"
  external "C" cv=  TwoPhaseMedium_specificHeatCapacityCv_C_impl(state, mediumName, libraryName, substanceName)
    annotation(Include="#include \"externalmedialib.h\"", Library="ExternalMediaLib", IncludeDirectory="modelica://ExternalMedia/Resources/Include", LibraryDirectory="modelica://ExternalMedia/Resources/Library");
*/
          annotation(Inline = true);
        end specificHeatCapacityCv;

        redeclare replaceable function extends density
        "Return density from state"
          // Standard definition
        algorithm
          d := state.d;
          /*  // If special definition in "C"
  external "C" d=  TwoPhaseMedium_density_C_impl(state, mediumName, libraryName, substanceName)
    annotation(Include="#include \"externalmedialib.h\"", Library="ExternalMediaLib", IncludeDirectory="modelica://ExternalMedia/Resources/Include", LibraryDirectory="modelica://ExternalMedia/Resources/Library");
*/
          annotation(Inline = true);
        end density;

        redeclare replaceable function extends density_derh_p
        "Return derivative of density wrt enthalpy at constant pressure from state"
          // Standard definition
        algorithm
          ddhp := state.ddhp;
          /*  // If special definition in "C"
  external "C" ddhp=  TwoPhaseMedium_density_derh_p_C_impl(state, mediumName, libraryName, substanceName)
    annotation(Include="#include \"externalmedialib.h\"", Library="ExternalMediaLib", IncludeDirectory="modelica://ExternalMedia/Resources/Include", LibraryDirectory="modelica://ExternalMedia/Resources/Library");
*/
          annotation(Inline = true);
        end density_derh_p;

        redeclare replaceable function extends density_derp_h
        "Return derivative of density wrt pressure at constant enthalpy from state"
          // Standard definition
        algorithm
          ddph := state.ddph;
          /*  // If special definition in "C"
  external "C" ddph=  TwoPhaseMedium_density_derp_h_C_impl(state, mediumName, libraryName, substanceName)
    annotation(Include="#include \"externalmedialib.h\"", Library="ExternalMediaLib", IncludeDirectory="modelica://ExternalMedia/Resources/Include", LibraryDirectory="modelica://ExternalMedia/Resources/Library");
*/
          annotation(Inline = true);
        end density_derp_h;

        redeclare replaceable function extends density_derp_T
        algorithm
          ddpT := state.kappa*state.d;
        end density_derp_T;

        redeclare replaceable function extends density_derT_p
        algorithm
          ddTp :=-state.beta*state.d;
        end density_derT_p;

        redeclare replaceable function extends dynamicViscosity
        "Return dynamic viscosity from state"
          // Standard definition
        algorithm
          eta := state.eta;
          /*  // If special definition in "C"
  external "C" eta=  TwoPhaseMedium_dynamicViscosity_C_impl(state, mediumName, libraryName, substanceName)
    annotation(Include="#include \"externalmedialib.h\"", Library="ExternalMediaLib", IncludeDirectory="modelica://ExternalMedia/Resources/Include", LibraryDirectory="modelica://ExternalMedia/Resources/Library");
*/
          annotation(Inline = true);
        end dynamicViscosity;

        redeclare replaceable function extends specificEnthalpy
        "Return specific enthalpy from state"
          // Standard definition
        algorithm
          h := state.h;
          /*  // If special definition in "C"
  external "C" h=  TwoPhaseMedium_specificEnthalpy_C_impl(state, mediumName, libraryName, substanceName)
    annotation(Include="#include \"externalmedialib.h\"", Library="ExternalMediaLib", IncludeDirectory="modelica://ExternalMedia/Resources/Include", LibraryDirectory="modelica://ExternalMedia/Resources/Library");
*/
          annotation(Inline = true);
        end specificEnthalpy;

        redeclare replaceable function extends specificInternalEnergy
        "Returns specific internal energy"
          extends Modelica.Icons.Function;
        algorithm
          u := specificEnthalpy(state) - pressure(state)/density(state);
        end specificInternalEnergy;

        redeclare replaceable function extends isothermalCompressibility
        "Return isothermal compressibility from state"
          // Standard definition
        algorithm
          kappa := state.kappa;
          /*  // If special definition in "C"
  external "C" kappa=  TwoPhaseMedium_isothermalCompressibility_C_impl(state, mediumName, libraryName, substanceName)
    annotation(Include="#include \"externalmedialib.h\"", Library="ExternalMediaLib", IncludeDirectory="modelica://ExternalMedia/Resources/Include", LibraryDirectory="modelica://ExternalMedia/Resources/Library");
*/
          annotation(Inline = true);
        end isothermalCompressibility;

        redeclare replaceable function extends thermalConductivity
        "Return thermal conductivity from state"
          // Standard definition
        algorithm
          lambda := state.lambda;
          /*  // If special definition in "C"
  external "C" lambda=  TwoPhaseMedium_thermalConductivity_C_impl(state, mediumName, libraryName, substanceName)
    annotation(Include="#include \"externalmedialib.h\"", Library="ExternalMediaLib", IncludeDirectory="modelica://ExternalMedia/Resources/Include", LibraryDirectory="modelica://ExternalMedia/Resources/Library");
*/
          annotation(Inline = true);
        end thermalConductivity;

        redeclare replaceable function extends pressure
        "Return pressure from state"
          // Standard definition
        algorithm
          p := state.p;
          /*  // If special definition in "C"
  external "C" p=  TwoPhaseMedium_pressure_C_impl(state, mediumName, libraryName, substanceName)
    annotation(Include="#include \"externalmedialib.h\"", Library="ExternalMediaLib", IncludeDirectory="modelica://ExternalMedia/Resources/Include", LibraryDirectory="modelica://ExternalMedia/Resources/Library");
*/
          annotation(Inline = true);
        end pressure;

        redeclare replaceable function extends specificEntropy
        "Return specific entropy from state"
          // Standard definition
        algorithm
          s := state.s;
          /*  // If special definition in "C"
    external "C" s=  TwoPhaseMedium_specificEntropy_C_impl(state, mediumName, libraryName, substanceName)
    annotation(Include="#include \"externalmedialib.h\"", Library="ExternalMediaLib", IncludeDirectory="modelica://ExternalMedia/Resources/Include", LibraryDirectory="modelica://ExternalMedia/Resources/Library");
*/
          annotation(Inline = true);
        end specificEntropy;

        redeclare replaceable function extends isentropicEnthalpy
        external "C" h_is = TwoPhaseMedium_isentropicEnthalpy_C_impl(p_downstream, refState,
         mediumName, libraryName, substanceName)
          annotation(Include="#include \"externalmedialib.h\"", Library="ExternalMediaLib", IncludeDirectory="modelica://ExternalMedia/Resources/Include", LibraryDirectory="modelica://ExternalMedia/Resources/Library");
        end isentropicEnthalpy;

        redeclare replaceable function setSat_p
        "Return saturation properties from p"
          extends Modelica.Icons.Function;
          input AbsolutePressure p "pressure";
          output SaturationProperties sat "saturation property record";
        external "C" TwoPhaseMedium_setSat_p_C_impl(p, sat, mediumName, libraryName, substanceName)
          annotation(Include="#include \"externalmedialib.h\"", Library="ExternalMediaLib", IncludeDirectory="modelica://ExternalMedia/Resources/Include", LibraryDirectory="modelica://ExternalMedia/Resources/Library");
        end setSat_p;

        replaceable function setSat_p_state
        "Return saturation properties from the state"
          extends Modelica.Icons.Function;
          input ThermodynamicState state;
          output SaturationProperties sat "saturation property record";
          // Standard definition
        algorithm
          sat:=setSat_p(state.p);
          //Redeclare this function for more efficient implementations avoiding the repeated computation of saturation properties
        /*  // If special definition in "C"
  external "C" TwoPhaseMedium_setSat_p_state_C_impl(state, sat)
    annotation(Include="#include \"externalmedialib.h\"", Library="ExternalMediaLib", IncludeDirectory="modelica://ExternalMedia/Resources/Include", LibraryDirectory="modelica://ExternalMedia/Resources/Library");
*/
          annotation(Inline = true);
        end setSat_p_state;

        redeclare replaceable function setSat_T
        "Return saturation properties from p"
          extends Modelica.Icons.Function;
          input Temperature T "temperature";
          output SaturationProperties sat "saturation property record";
        external "C" TwoPhaseMedium_setSat_T_C_impl(T, sat, mediumName, libraryName, substanceName)
          annotation(Include="#include \"externalmedialib.h\"", Library="ExternalMediaLib", IncludeDirectory="modelica://ExternalMedia/Resources/Include", LibraryDirectory="modelica://ExternalMedia/Resources/Library");
        end setSat_T;

        replaceable function setSat_T_state
        "Return saturation properties from the state"
          extends Modelica.Icons.Function;
          input ThermodynamicState state;
          output SaturationProperties sat "saturation property record";
          // Standard definition
        algorithm
          sat:=setSat_T(state.T);
          //Redeclare this function for more efficient implementations avoiding the repeated computation of saturation properties
        /*  // If special definition in "C"
  external "C" TwoPhaseMedium_setSat_T_state_C_impl(state, sat)
    annotation(Include="#include \"externalmedialib.h\"", Library="ExternalMediaLib", IncludeDirectory="modelica://ExternalMedia/Resources/Include", LibraryDirectory="modelica://ExternalMedia/Resources/Library");
*/
          annotation(Inline = true);
        end setSat_T_state;

        redeclare replaceable function extends setBubbleState
        "set the thermodynamic state on the bubble line"
          extends Modelica.Icons.Function;
          input SaturationProperties sat "saturation point";
          input FixedPhase phase(min = 1, max = 2) =  1
          "phase: default is one phase";
          output ThermodynamicState state "complete thermodynamic state info";
          // Standard definition
        algorithm
          state :=setState_ph(sat.psat, sat.hl, phase);
          /*  // If special definition in "C"
  external "C" TwoPhaseMedium_setBubbleState_C_impl(sat, phase, state, mediumName, libraryName, substanceName)
    annotation(Include="#include \"externalmedialib.h\"", Library="ExternalMediaLib", IncludeDirectory="modelica://ExternalMedia/Resources/Include", LibraryDirectory="modelica://ExternalMedia/Resources/Library");
*/
          annotation(Inline = true);
        end setBubbleState;

        redeclare replaceable function extends setDewState
        "set the thermodynamic state on the dew line"
          extends Modelica.Icons.Function;
          input SaturationProperties sat "saturation point";
          input FixedPhase phase(min = 1, max = 2) = 1
          "phase: default is one phase";
          output ThermodynamicState state "complete thermodynamic state info";
          // Standard definition
        algorithm
          state :=setState_ph(sat.psat, sat.hv, phase);
          /*  // If special definition in "C"
  external "C" TwoPhaseMedium_setDewState_C_impl(sat, phase, state, mediumName, libraryName, substanceName)
    annotation(Include="#include \"externalmedialib.h\"", Library="ExternalMediaLib", IncludeDirectory="modelica://ExternalMedia/Resources/Include", LibraryDirectory="modelica://ExternalMedia/Resources/Library");
*/
          annotation(Inline = true);
        end setDewState;

        redeclare replaceable function extends saturationTemperature
          // Standard definition
        algorithm
          T :=saturationTemperature_sat(setSat_p(p));
          /*  // If special definition in "C"
  external "C" T=  TwoPhaseMedium_saturationTemperature_C_impl(p, mediumName, libraryName, substanceName)
    annotation(Include="#include \"externalmedialib.h\"", Library="ExternalMediaLib", IncludeDirectory="modelica://ExternalMedia/Resources/Include", LibraryDirectory="modelica://ExternalMedia/Resources/Library");
*/
          annotation(Inline = true);
        end saturationTemperature;

        redeclare function extends saturationTemperature_sat

          annotation(Inline = true);
        end saturationTemperature_sat;

        redeclare replaceable function extends saturationTemperature_derp "Returns derivative of saturation temperature w.r.t.. pressureBeing this function inefficient, it is strongly recommended to use saturationTemperature_derp_sat
     and never use saturationTemperature_derp directly"
        external "C" dTp = TwoPhaseMedium_saturationTemperature_derp_C_impl(p, mediumName, libraryName, substanceName)
          annotation(Include="#include \"externalmedialib.h\"", Library="ExternalMediaLib", IncludeDirectory="modelica://ExternalMedia/Resources/Include", LibraryDirectory="modelica://ExternalMedia/Resources/Library");
        end saturationTemperature_derp;

        redeclare replaceable function saturationTemperature_derp_sat
        "Returns derivative of saturation temperature w.r.t.. pressure"
          extends Modelica.Icons.Function;
          input SaturationProperties sat "saturation property record";
          output Real dTp
          "derivative of saturation temperature w.r.t. pressure";
          // Standard definition
        algorithm
          dTp := sat.dTp;
          /*  // If special definition in "C"
  external "C" dTp=  TwoPhaseMedium_saturationTemperature_derp_sat_C_impl(sat.psat, sat.Tsat, sat.uniqueID, mediumName, libraryName, substanceName)
    annotation(Include="#include \"externalmedialib.h\"", Library="ExternalMediaLib", IncludeDirectory="modelica://ExternalMedia/Resources/Include", LibraryDirectory="modelica://ExternalMedia/Resources/Library");
*/
          annotation(Inline = true);
        end saturationTemperature_derp_sat;

        redeclare replaceable function extends dBubbleDensity_dPressure
        "Returns bubble point density derivative"
          // Standard definition
        algorithm
          ddldp := sat.ddldp;
          /*  // If special definition in "C"
  external "C" ddldp=  TwoPhaseMedium_dBubbleDensity_dPressure_C_impl(sat, mediumName, libraryName, substanceName)
    annotation(Include="#include \"externalmedialib.h\"", Library="ExternalMediaLib", IncludeDirectory="modelica://ExternalMedia/Resources/Include", LibraryDirectory="modelica://ExternalMedia/Resources/Library");
*/
          annotation(Inline = true);
        end dBubbleDensity_dPressure;

        redeclare replaceable function extends dDewDensity_dPressure
        "Returns dew point density derivative"
          // Standard definition
        algorithm
          ddvdp := sat.ddvdp;
          /*  // If special definition in "C"
  external "C" ddvdp=  TwoPhaseMedium_dDewDensity_dPressure_C_impl(sat, mediumName, libraryName, substanceName)
    annotation(Include="#include \"externalmedialib.h\"", Library="ExternalMediaLib", IncludeDirectory="modelica://ExternalMedia/Resources/Include", LibraryDirectory="modelica://ExternalMedia/Resources/Library");
*/
          annotation(Inline = true);
        end dDewDensity_dPressure;

        redeclare replaceable function extends dBubbleEnthalpy_dPressure
        "Returns bubble point specific enthalpy derivative"
          // Standard definition
        algorithm
          dhldp := sat.dhldp;
          /*  // If special definition in "C"
  external "C" dhldp=  TwoPhaseMedium_dBubbleEnthalpy_dPressure_C_impl(sat, mediumName, libraryName, substanceName)
    annotation(Include="#include \"externalmedialib.h\"", Library="ExternalMediaLib", IncludeDirectory="modelica://ExternalMedia/Resources/Include", LibraryDirectory="modelica://ExternalMedia/Resources/Library");
*/
          annotation(Inline = true);
        end dBubbleEnthalpy_dPressure;

        redeclare replaceable function extends dDewEnthalpy_dPressure
        "Returns dew point specific enthalpy derivative"
          // Standard definition
        algorithm
          dhvdp := sat.dhvdp;
          /*  // If special definition in "C"
  external "C" dhvdp=  TwoPhaseMedium_dDewEnthalpy_dPressure_C_impl(sat, mediumName, libraryName, substanceName)
    annotation(Include="#include \"externalmedialib.h\"", Library="ExternalMediaLib", IncludeDirectory="modelica://ExternalMedia/Resources/Include", LibraryDirectory="modelica://ExternalMedia/Resources/Library");
*/
          annotation(Inline = true);
        end dDewEnthalpy_dPressure;

        redeclare replaceable function extends bubbleDensity
        "Returns bubble point density"
          // Standard definition
        algorithm
          dl := sat.dl;
          /*  // If special definition in "C"
  external "C" dl=  TwoPhaseMedium_bubbleDensity_C_impl(sat, mediumName, libraryName, substanceName)
    annotation(Include="#include \"externalmedialib.h\"", Library="ExternalMediaLib", IncludeDirectory="modelica://ExternalMedia/Resources/Include", LibraryDirectory="modelica://ExternalMedia/Resources/Library");
*/
          annotation(Inline = true);
        end bubbleDensity;

        redeclare replaceable function extends dewDensity
        "Returns dew point density"
          // Standard definition
        algorithm
          dv := sat.dv;
          /*  // If special definition in "C"
  external "C" dv=  TwoPhaseMedium_dewDensity_C_impl(sat, mediumName, libraryName, substanceName)
    annotation(Include="#include \"externalmedialib.h\"", Library="ExternalMediaLib", IncludeDirectory="modelica://ExternalMedia/Resources/Include", LibraryDirectory="modelica://ExternalMedia/Resources/Library");
*/
          annotation(Inline = true);
        end dewDensity;

        redeclare replaceable function extends bubbleEnthalpy
        "Returns bubble point specific enthalpy"
          // Standard definition
        algorithm
          hl := sat.hl;
          /*  // If special definition in "C"
  external "C" hl=  TwoPhaseMedium_bubbleEnthalpy_C_impl(sat, mediumName, libraryName, substanceName)
    annotation(Include="#include \"externalmedialib.h\"", Library="ExternalMediaLib", IncludeDirectory="modelica://ExternalMedia/Resources/Include", LibraryDirectory="modelica://ExternalMedia/Resources/Library");
*/
          annotation(Inline = true);
        end bubbleEnthalpy;

        redeclare replaceable function extends dewEnthalpy
        "Returns dew point specific enthalpy"
          // Standard definition
        algorithm
          hv := sat.hv;
          /*  // If special definition in "C"
  external "C" hv=  TwoPhaseMedium_dewEnthalpy_C_impl(sat, mediumName, libraryName, substanceName)
    annotation(Include="#include \"externalmedialib.h\"", Library="ExternalMediaLib", IncludeDirectory="modelica://ExternalMedia/Resources/Include", LibraryDirectory="modelica://ExternalMedia/Resources/Library");
*/
          annotation(Inline = true);
        end dewEnthalpy;

        redeclare replaceable function extends saturationPressure
          // Standard definition
        algorithm
          p :=saturationPressure_sat(setSat_T(T));
          /*  // If special definition in "C"
  external "C" p=  TwoPhaseMedium_saturationPressure_C_impl(T, mediumName, libraryName, substanceName)
    annotation(Include="#include \"externalmedialib.h\"", Library="ExternalMediaLib", IncludeDirectory="modelica://ExternalMedia/Resources/Include", LibraryDirectory="modelica://ExternalMedia/Resources/Library");
*/
          annotation(Inline = false,
                     LateInline = true,
                     derivative = saturationPressure_der);
        end saturationPressure;

        function saturationPressure_der
        "Return saturation pressure time derivative"
          extends Modelica.Icons.Function;
          input Temperature T "temperature";
          input Real T_der "Temperature derivative";
          output Real p_der "saturation pressure derivative";
          // Standard definition
        algorithm
          p_der :=T_der/saturationTemperature_derp_sat(setSat_T(T));
          annotation(Inline = true);
        end saturationPressure_der;

        redeclare function extends saturationPressure_sat

          annotation(Inline = true);
        end saturationPressure_sat;

        redeclare replaceable function extends surfaceTension
        "Returns surface tension sigma in the two phase region"
          //Standard definition
        algorithm
          sigma := sat.sigma;
          /*  //If special definition in "C"
  external "C" sigma=  TwoPhaseMedium_surfaceTension_C_impl(sat, mediumName, libraryName, substanceName)
    annotation(Include="#include \"externalmedialib.h\"", Library="ExternalMediaLib", IncludeDirectory="modelica://ExternalMedia/Resources/Include", LibraryDirectory="modelica://ExternalMedia/Resources/Library");
*/
          annotation(Inline = true);
        end surfaceTension;

        redeclare replaceable function extends bubbleEntropy
        "Returns bubble point specific entropy"
          //Standard definition
        algorithm
          sl := specificEntropy(setBubbleState(sat));
          /*  //If special definition in "C"
  external "C" sl=  TwoPhaseMedium_bubbleEntropy_C_impl(sat, mediumName, libraryName, substanceName)
    annotation(Include="#include \"externalmedialib.h\"", Library="ExternalMediaLib", IncludeDirectory="modelica://ExternalMedia/Resources/Include", LibraryDirectory="modelica://ExternalMedia/Resources/Library");
*/
          annotation(Inline = true);
        end bubbleEntropy;

        redeclare replaceable function extends dewEntropy
        "Returns dew point specific entropy"
          //Standard definition
        algorithm
          sv := specificEntropy(setDewState(sat));
          /*  //If special definition in "C"
  external "C" sv=  TwoPhaseMedium_dewEntropy_C_impl(sat, mediumName, libraryName, substanceName)
    annotation(Include="#include \"externalmedialib.h\"", Library="ExternalMediaLib", IncludeDirectory="modelica://ExternalMedia/Resources/Include", LibraryDirectory="modelica://ExternalMedia/Resources/Library");
*/
          annotation(Inline = true);
        end dewEntropy;
      end ExternalTwoPhaseMedium;
    end BaseClasses;
  end Media;
  annotation(uses(Modelica(version="3.2.1")),
  Documentation(info="<html>
<p>The <b>ExternalMedia</b> library provides a framework for interfacing external codes computing fluid properties to Modelica.Media-compatible component models. The library has been designed with two main goals: maximizing the efficiency of the code, while minimizing the amount of extra code required to interface existing external codes to the library.</p>
<p>The library covers pure fluids models, possibly two-phase, compliant with the <a href=\"modelica://Modelica.Media.Interfaces.PartialTwoPhaseMedium\">Modelica.Media.Interfaces.PartialTwoPhaseMedium</a> interface. </p>
<p>Two external softwares for fluid property computation are currently suppored by the ExternalMedia library:</p>
<ul>
<li><a href=\"http://www.fluidprop.com\">FluidProp</a>, formerly developed at TU Delft and currently devloped and maintained by Asimptote</li>
<li><a href=\"http://coolprop.org\">CoolProp</a>, developed at the University of Liege and at the Technical University of Denmark (DTU)</li>
</ul>
<p>The library has been tested with the Dymola and OpenModelica tools under the Windows operating system. If you are interested in the support of other tools, operating systems, and external fluid property computation codes, please contact the developers.</p>
<p>Main contributors: Francesco Casella, Christoph Richter, Roberto Bonifetto, Ian Bell.</p>
<p><b>The code is licensed under the Modelica License 2. </b>For license conditions (including the disclaimer of warranty) visit <a href=\"https://www.modelica.org/licenses/ModelicaLicense2\">https://www.modelica.org/licenses/ModelicaLicense2</a>. </p>
<p>Copyright &copy; 2006-2014, Politecnico di Milano, TU Braunschweig, Politecnico di Torino, Universit&eacute; de Liege.</p>
</html>"));
end ExternalMedia;

package Modelica "Modelica Standard Library - Version 3.2.1 (Build 2)"
extends Modelica.Icons.Package;

  package Media "Library of media property models"
  extends Modelica.Icons.Package;
  import SI = Modelica.SIunits;
  import Cv = Modelica.SIunits.Conversions;

  package Interfaces "Interfaces for media models"
    extends Modelica.Icons.InterfacesPackage;

    partial package PartialMedium
    "Partial medium properties (base package of all media packages)"
      extends Modelica.Media.Interfaces.Types;
      extends Modelica.Icons.MaterialPropertiesPackage;

      // Constants to be set in Medium
      constant Modelica.Media.Interfaces.Choices.IndependentVariables
        ThermoStates "Enumeration type for independent variables";
      constant String mediumName="unusablePartialMedium" "Name of the medium";
      constant String substanceNames[:]={mediumName}
      "Names of the mixture substances. Set substanceNames={mediumName} if only one substance.";
      constant String extraPropertiesNames[:]=fill("", 0)
      "Names of the additional (extra) transported properties. Set extraPropertiesNames=fill(\"\",0) if unused";
      constant Boolean singleState
      "= true, if u and d are not a function of pressure";
      constant Boolean reducedX=true
      "= true if medium contains the equation sum(X) = 1.0; set reducedX=true if only one substance (see docu for details)";
      constant Boolean fixedX=false
      "= true if medium contains the equation X = reference_X";
      constant AbsolutePressure reference_p=101325
      "Reference pressure of Medium: default 1 atmosphere";
      constant Temperature reference_T=298.15
      "Reference temperature of Medium: default 25 deg Celsius";
      constant MassFraction reference_X[nX]=fill(1/nX, nX)
      "Default mass fractions of medium";
      constant AbsolutePressure p_default=101325
      "Default value for pressure of medium (for initialization)";
      constant Temperature T_default=Modelica.SIunits.Conversions.from_degC(20)
      "Default value for temperature of medium (for initialization)";
      constant SpecificEnthalpy h_default=specificEnthalpy_pTX(
              p_default,
              T_default,
              X_default)
      "Default value for specific enthalpy of medium (for initialization)";
      constant MassFraction X_default[nX]=reference_X
      "Default value for mass fractions of medium (for initialization)";

      final constant Integer nS=size(substanceNames, 1) "Number of substances"
        annotation (Evaluate=true);
      constant Integer nX=nS "Number of mass fractions" annotation (Evaluate=true);
      constant Integer nXi=if fixedX then 0 else if reducedX then nS - 1 else nS
      "Number of structurally independent mass fractions (see docu for details)"
        annotation (Evaluate=true);

      final constant Integer nC=size(extraPropertiesNames, 1)
      "Number of extra (outside of standard mass-balance) transported properties"
        annotation (Evaluate=true);
      constant Real C_nominal[nC](min=fill(Modelica.Constants.eps, nC)) = 1.0e-6*
        ones(nC) "Default for the nominal values for the extra properties";
      replaceable record FluidConstants =
          Modelica.Media.Interfaces.Types.Basic.FluidConstants
      "Critical, triple, molecular and other standard data of fluid";

      replaceable record ThermodynamicState
      "Minimal variable set that is available as input argument to every medium function"
        extends Modelica.Icons.Record;
      end ThermodynamicState;

      replaceable partial model BaseProperties
      "Base properties (p, d, T, h, u, R, MM and, if applicable, X and Xi) of a medium"
        InputAbsolutePressure p "Absolute pressure of medium";
        InputMassFraction[nXi] Xi(start=reference_X[1:nXi])
        "Structurally independent mass fractions";
        InputSpecificEnthalpy h "Specific enthalpy of medium";
        Density d "Density of medium";
        Temperature T "Temperature of medium";
        MassFraction[nX] X(start=reference_X)
        "Mass fractions (= (component mass)/total mass  m_i/m)";
        SpecificInternalEnergy u "Specific internal energy of medium";
        SpecificHeatCapacity R "Gas constant (of mixture if applicable)";
        MolarMass MM "Molar mass (of mixture or single fluid)";
        ThermodynamicState state
        "Thermodynamic state record for optional functions";
        parameter Boolean preferredMediumStates=false
        "= true if StateSelect.prefer shall be used for the independent property variables of the medium"
          annotation (Evaluate=true, Dialog(tab="Advanced"));
        parameter Boolean standardOrderComponents=true
        "If true, and reducedX = true, the last element of X will be computed from the other ones";
        SI.Conversions.NonSIunits.Temperature_degC T_degC=
            Modelica.SIunits.Conversions.to_degC(T)
        "Temperature of medium in [degC]";
        SI.Conversions.NonSIunits.Pressure_bar p_bar=
            Modelica.SIunits.Conversions.to_bar(p)
        "Absolute pressure of medium in [bar]";

        // Local connector definition, used for equation balancing check
        connector InputAbsolutePressure = input SI.AbsolutePressure
        "Pressure as input signal connector";
        connector InputSpecificEnthalpy = input SI.SpecificEnthalpy
        "Specific enthalpy as input signal connector";
        connector InputMassFraction = input SI.MassFraction
        "Mass fraction as input signal connector";

      equation
        if standardOrderComponents then
          Xi = X[1:nXi];

          if fixedX then
            X = reference_X;
          end if;
          if reducedX and not fixedX then
            X[nX] = 1 - sum(Xi);
          end if;
          for i in 1:nX loop
            assert(X[i] >= -1.e-5 and X[i] <= 1 + 1.e-5, "Mass fraction X[" +
              String(i) + "] = " + String(X[i]) + "of substance " +
              substanceNames[i] + "\nof medium " + mediumName +
              " is not in the range 0..1");
          end for;

        end if;

        assert(p >= 0.0, "Pressure (= " + String(p) + " Pa) of medium \"" +
          mediumName + "\" is negative\n(Temperature = " + String(T) + " K)");
        annotation (Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,
                  -100},{100,100}}), graphics={Rectangle(
                extent={{-100,100},{100,-100}},
                fillColor={255,255,255},
                fillPattern=FillPattern.Solid,
                lineColor={0,0,255}), Text(
                extent={{-152,164},{152,102}},
                textString="%name",
                lineColor={0,0,255})}), Documentation(info="<html>
<p>
Model <b>BaseProperties</b> is a model within package <b>PartialMedium</b>
and contains the <b>declarations</b> of the minimum number of
variables that every medium model is supposed to support.
A specific medium inherits from model <b>BaseProperties</b> and provides
the equations for the basic properties.</p>
<p>
The BaseProperties model contains the following <b>7+nXi variables</b>
(nXi is the number of independent mass fractions defined in package
PartialMedium):
</p>
<table border=1 cellspacing=0 cellpadding=2>
  <tr><td valign=\"top\"><b>Variable</b></td>
      <td valign=\"top\"><b>Unit</b></td>
      <td valign=\"top\"><b>Description</b></td></tr>
  <tr><td valign=\"top\">T</td>
      <td valign=\"top\">K</td>
      <td valign=\"top\">temperature</td></tr>
  <tr><td valign=\"top\">p</td>
      <td valign=\"top\">Pa</td>
      <td valign=\"top\">absolute pressure</td></tr>
  <tr><td valign=\"top\">d</td>
      <td valign=\"top\">kg/m3</td>
      <td valign=\"top\">density</td></tr>
  <tr><td valign=\"top\">h</td>
      <td valign=\"top\">J/kg</td>
      <td valign=\"top\">specific enthalpy</td></tr>
  <tr><td valign=\"top\">u</td>
      <td valign=\"top\">J/kg</td>
      <td valign=\"top\">specific internal energy</td></tr>
  <tr><td valign=\"top\">Xi[nXi]</td>
      <td valign=\"top\">kg/kg</td>
      <td valign=\"top\">independent mass fractions m_i/m</td></tr>
  <tr><td valign=\"top\">R</td>
      <td valign=\"top\">J/kg.K</td>
      <td valign=\"top\">gas constant</td></tr>
  <tr><td valign=\"top\">M</td>
      <td valign=\"top\">kg/mol</td>
      <td valign=\"top\">molar mass</td></tr>
</table>
<p>
In order to implement an actual medium model, one can extend from this
base model and add <b>5 equations</b> that provide relations among
these variables. Equations will also have to be added in order to
set all the variables within the ThermodynamicState record state.</p>
<p>
If standardOrderComponents=true, the full composition vector X[nX]
is determined by the equations contained in this base class, depending
on the independent mass fraction vector Xi[nXi].</p>
<p>Additional <b>2 + nXi</b> equations will have to be provided
when using the BaseProperties model, in order to fully specify the
thermodynamic conditions. The input connector qualifier applied to
p, h, and nXi indirectly declares the number of missing equations,
permitting advanced equation balance checking by Modelica tools.
Please note that this doesn't mean that the additional equations
should be connection equations, nor that exactly those variables
should be supplied, in order to complete the model.
For further information, see the Modelica.Media User's guide, and
Section 4.7 (Balanced Models) of the Modelica 3.0 specification.</p>
</html>"));
      end BaseProperties;

      replaceable partial function setState_pTX
      "Return thermodynamic state as function of p, T and composition X or Xi"
        extends Modelica.Icons.Function;
        input AbsolutePressure p "Pressure";
        input Temperature T "Temperature";
        input MassFraction X[:]=reference_X "Mass fractions";
        output ThermodynamicState state "Thermodynamic state record";
      end setState_pTX;

      replaceable partial function setState_phX
      "Return thermodynamic state as function of p, h and composition X or Xi"
        extends Modelica.Icons.Function;
        input AbsolutePressure p "Pressure";
        input SpecificEnthalpy h "Specific enthalpy";
        input MassFraction X[:]=reference_X "Mass fractions";
        output ThermodynamicState state "Thermodynamic state record";
      end setState_phX;

      replaceable partial function setState_psX
      "Return thermodynamic state as function of p, s and composition X or Xi"
        extends Modelica.Icons.Function;
        input AbsolutePressure p "Pressure";
        input SpecificEntropy s "Specific entropy";
        input MassFraction X[:]=reference_X "Mass fractions";
        output ThermodynamicState state "Thermodynamic state record";
      end setState_psX;

      replaceable partial function setState_dTX
      "Return thermodynamic state as function of d, T and composition X or Xi"
        extends Modelica.Icons.Function;
        input Density d "Density";
        input Temperature T "Temperature";
        input MassFraction X[:]=reference_X "Mass fractions";
        output ThermodynamicState state "Thermodynamic state record";
      end setState_dTX;

      replaceable partial function setSmoothState
      "Return thermodynamic state so that it smoothly approximates: if x > 0 then state_a else state_b"
        extends Modelica.Icons.Function;
        input Real x "m_flow or dp";
        input ThermodynamicState state_a "Thermodynamic state if x > 0";
        input ThermodynamicState state_b "Thermodynamic state if x < 0";
        input Real x_small(min=0)
        "Smooth transition in the region -x_small < x < x_small";
        output ThermodynamicState state
        "Smooth thermodynamic state for all x (continuous and differentiable)";
        annotation (Documentation(info="<html>
<p>
This function is used to approximate the equation
</p>
<pre>
    state = <b>if</b> x &gt; 0 <b>then</b> state_a <b>else</b> state_b;
</pre>

<p>
by a smooth characteristic, so that the expression is continuous and differentiable:
</p>

<pre>
   state := <b>smooth</b>(1, <b>if</b> x &gt;  x_small <b>then</b> state_a <b>else</b>
                      <b>if</b> x &lt; -x_small <b>then</b> state_b <b>else</b> f(state_a, state_b));
</pre>

<p>
This is performed by applying function <b>Media.Common.smoothStep</b>(..)
on every element of the thermodynamic state record.
</p>

<p>
If <b>mass fractions</b> X[:] are approximated with this function then this can be performed
for all <b>nX</b> mass fractions, instead of applying it for nX-1 mass fractions and computing
the last one by the mass fraction constraint sum(X)=1. The reason is that the approximating function has the
property that sum(state.X) = 1, provided sum(state_a.X) = sum(state_b.X) = 1.
This can be shown by evaluating the approximating function in the abs(x) &lt; x_small
region (otherwise state.X is either state_a.X or state_b.X):
</p>

<pre>
    X[1]  = smoothStep(x, X_a[1] , X_b[1] , x_small);
    X[2]  = smoothStep(x, X_a[2] , X_b[2] , x_small);
       ...
    X[nX] = smoothStep(x, X_a[nX], X_b[nX], x_small);
</pre>

<p>
or
</p>

<pre>
    X[1]  = c*(X_a[1]  - X_b[1])  + (X_a[1]  + X_b[1])/2
    X[2]  = c*(X_a[2]  - X_b[2])  + (X_a[2]  + X_b[2])/2;
       ...
    X[nX] = c*(X_a[nX] - X_b[nX]) + (X_a[nX] + X_b[nX])/2;
    c     = (x/x_small)*((x/x_small)^2 - 3)/4
</pre>

<p>
Summing all mass fractions together results in
</p>

<pre>
    sum(X) = c*(sum(X_a) - sum(X_b)) + (sum(X_a) + sum(X_b))/2
           = c*(1 - 1) + (1 + 1)/2
           = 1
</pre>

</html>"));
      end setSmoothState;

      replaceable partial function dynamicViscosity "Return dynamic viscosity"
        extends Modelica.Icons.Function;
        input ThermodynamicState state "Thermodynamic state record";
        output DynamicViscosity eta "Dynamic viscosity";
      end dynamicViscosity;

      replaceable partial function thermalConductivity
      "Return thermal conductivity"
        extends Modelica.Icons.Function;
        input ThermodynamicState state "Thermodynamic state record";
        output ThermalConductivity lambda "Thermal conductivity";
      end thermalConductivity;

      replaceable function prandtlNumber "Return the Prandtl number"
        extends Modelica.Icons.Function;
        input ThermodynamicState state "Thermodynamic state record";
        output PrandtlNumber Pr "Prandtl number";
      algorithm
        Pr := dynamicViscosity(state)*specificHeatCapacityCp(state)/
          thermalConductivity(state);
      end prandtlNumber;

      replaceable partial function pressure "Return pressure"
        extends Modelica.Icons.Function;
        input ThermodynamicState state "Thermodynamic state record";
        output AbsolutePressure p "Pressure";
      end pressure;

      replaceable partial function temperature "Return temperature"
        extends Modelica.Icons.Function;
        input ThermodynamicState state "Thermodynamic state record";
        output Temperature T "Temperature";
      end temperature;

      replaceable partial function density "Return density"
        extends Modelica.Icons.Function;
        input ThermodynamicState state "Thermodynamic state record";
        output Density d "Density";
      end density;

      replaceable partial function specificEnthalpy "Return specific enthalpy"
        extends Modelica.Icons.Function;
        input ThermodynamicState state "Thermodynamic state record";
        output SpecificEnthalpy h "Specific enthalpy";
      end specificEnthalpy;

      replaceable partial function specificInternalEnergy
      "Return specific internal energy"
        extends Modelica.Icons.Function;
        input ThermodynamicState state "Thermodynamic state record";
        output SpecificEnergy u "Specific internal energy";
      end specificInternalEnergy;

      replaceable partial function specificEntropy "Return specific entropy"
        extends Modelica.Icons.Function;
        input ThermodynamicState state "Thermodynamic state record";
        output SpecificEntropy s "Specific entropy";
      end specificEntropy;

      replaceable partial function specificGibbsEnergy
      "Return specific Gibbs energy"
        extends Modelica.Icons.Function;
        input ThermodynamicState state "Thermodynamic state record";
        output SpecificEnergy g "Specific Gibbs energy";
      end specificGibbsEnergy;

      replaceable partial function specificHelmholtzEnergy
      "Return specific Helmholtz energy"
        extends Modelica.Icons.Function;
        input ThermodynamicState state "Thermodynamic state record";
        output SpecificEnergy f "Specific Helmholtz energy";
      end specificHelmholtzEnergy;

      replaceable partial function specificHeatCapacityCp
      "Return specific heat capacity at constant pressure"
        extends Modelica.Icons.Function;
        input ThermodynamicState state "Thermodynamic state record";
        output SpecificHeatCapacity cp
        "Specific heat capacity at constant pressure";
      end specificHeatCapacityCp;

      function heatCapacity_cp = specificHeatCapacityCp
      "Alias for deprecated name";

      replaceable partial function specificHeatCapacityCv
      "Return specific heat capacity at constant volume"
        extends Modelica.Icons.Function;
        input ThermodynamicState state "Thermodynamic state record";
        output SpecificHeatCapacity cv
        "Specific heat capacity at constant volume";
      end specificHeatCapacityCv;

      function heatCapacity_cv = specificHeatCapacityCv
      "Alias for deprecated name";

      replaceable partial function isentropicExponent
      "Return isentropic exponent"
        extends Modelica.Icons.Function;
        input ThermodynamicState state "Thermodynamic state record";
        output IsentropicExponent gamma "Isentropic exponent";
      end isentropicExponent;

      replaceable partial function isentropicEnthalpy
      "Return isentropic enthalpy"
        extends Modelica.Icons.Function;
        input AbsolutePressure p_downstream "Downstream pressure";
        input ThermodynamicState refState "Reference state for entropy";
        output SpecificEnthalpy h_is "Isentropic enthalpy";
        annotation (Documentation(info="<html>
<p>
This function computes an isentropic state transformation:
</p>
<ol>
<li> A medium is in a particular state, refState.</li>
<li> The enthalpy at another state (h_is) shall be computed
     under the assumption that the state transformation from refState to h_is
     is performed with a change of specific entropy ds = 0 and the pressure of state h_is
     is p_downstream and the composition X upstream and downstream is assumed to be the same.</li>
</ol>

</html>"));
      end isentropicEnthalpy;

      replaceable partial function velocityOfSound "Return velocity of sound"
        extends Modelica.Icons.Function;
        input ThermodynamicState state "Thermodynamic state record";
        output VelocityOfSound a "Velocity of sound";
      end velocityOfSound;

      replaceable partial function isobaricExpansionCoefficient
      "Return overall the isobaric expansion coefficient beta"
        extends Modelica.Icons.Function;
        input ThermodynamicState state "Thermodynamic state record";
        output IsobaricExpansionCoefficient beta
        "Isobaric expansion coefficient";
        annotation (Documentation(info="<html>
<pre>
beta is defined as  1/v * der(v,T), with v = 1/d, at constant pressure p.
</pre>
</html>"));
      end isobaricExpansionCoefficient;

      function beta = isobaricExpansionCoefficient
      "Alias for isobaricExpansionCoefficient for user convenience";

      replaceable partial function isothermalCompressibility
      "Return overall the isothermal compressibility factor"
        extends Modelica.Icons.Function;
        input ThermodynamicState state "Thermodynamic state record";
        output SI.IsothermalCompressibility kappa "Isothermal compressibility";
        annotation (Documentation(info="<html>
<pre>

kappa is defined as - 1/v * der(v,p), with v = 1/d at constant temperature T.

</pre>
</html>"));
      end isothermalCompressibility;

      function kappa = isothermalCompressibility
      "Alias of isothermalCompressibility for user convenience";

      // explicit derivative functions for finite element models
      replaceable partial function density_derp_h
      "Return density derivative w.r.t. pressure at const specific enthalpy"
        extends Modelica.Icons.Function;
        input ThermodynamicState state "Thermodynamic state record";
        output DerDensityByPressure ddph "Density derivative w.r.t. pressure";
      end density_derp_h;

      replaceable partial function density_derh_p
      "Return density derivative w.r.t. specific enthalpy at constant pressure"
        extends Modelica.Icons.Function;
        input ThermodynamicState state "Thermodynamic state record";
        output DerDensityByEnthalpy ddhp
        "Density derivative w.r.t. specific enthalpy";
      end density_derh_p;

      replaceable partial function density_derp_T
      "Return density derivative w.r.t. pressure at const temperature"
        extends Modelica.Icons.Function;
        input ThermodynamicState state "Thermodynamic state record";
        output DerDensityByPressure ddpT "Density derivative w.r.t. pressure";
      end density_derp_T;

      replaceable partial function density_derT_p
      "Return density derivative w.r.t. temperature at constant pressure"
        extends Modelica.Icons.Function;
        input ThermodynamicState state "Thermodynamic state record";
        output DerDensityByTemperature ddTp
        "Density derivative w.r.t. temperature";
      end density_derT_p;

      replaceable partial function density_derX
      "Return density derivative w.r.t. mass fraction"
        extends Modelica.Icons.Function;
        input ThermodynamicState state "Thermodynamic state record";
        output Density[nX] dddX "Derivative of density w.r.t. mass fraction";
      end density_derX;

      replaceable partial function molarMass
      "Return the molar mass of the medium"
        extends Modelica.Icons.Function;
        input ThermodynamicState state "Thermodynamic state record";
        output MolarMass MM "Mixture molar mass";
      end molarMass;

      replaceable function specificEnthalpy_pTX
      "Return specific enthalpy from p, T, and X or Xi"
        extends Modelica.Icons.Function;
        input AbsolutePressure p "Pressure";
        input Temperature T "Temperature";
        input MassFraction X[:]=reference_X "Mass fractions";
        output SpecificEnthalpy h "Specific enthalpy";
      algorithm
        h := specificEnthalpy(setState_pTX(
                p,
                T,
                X));
        annotation (inverse(T=temperature_phX(
                      p,
                      h,
                      X)));
      end specificEnthalpy_pTX;

      replaceable function specificEntropy_pTX
      "Return specific enthalpy from p, T, and X or Xi"
        extends Modelica.Icons.Function;
        input AbsolutePressure p "Pressure";
        input Temperature T "Temperature";
        input MassFraction X[:]=reference_X "Mass fractions";
        output SpecificEntropy s "Specific entropy";
      algorithm
        s := specificEntropy(setState_pTX(
                p,
                T,
                X));

        annotation (inverse(T=temperature_psX(
                      p,
                      s,
                      X)));
      end specificEntropy_pTX;

      replaceable function density_pTX "Return density from p, T, and X or Xi"
        extends Modelica.Icons.Function;
        input AbsolutePressure p "Pressure";
        input Temperature T "Temperature";
        input MassFraction X[:] "Mass fractions";
        output Density d "Density";
      algorithm
        d := density(setState_pTX(
                p,
                T,
                X));
      end density_pTX;

      replaceable function temperature_phX
      "Return temperature from p, h, and X or Xi"
        extends Modelica.Icons.Function;
        input AbsolutePressure p "Pressure";
        input SpecificEnthalpy h "Specific enthalpy";
        input MassFraction X[:]=reference_X "Mass fractions";
        output Temperature T "Temperature";
      algorithm
        T := temperature(setState_phX(
                p,
                h,
                X));
      end temperature_phX;

      replaceable function density_phX "Return density from p, h, and X or Xi"
        extends Modelica.Icons.Function;
        input AbsolutePressure p "Pressure";
        input SpecificEnthalpy h "Specific enthalpy";
        input MassFraction X[:]=reference_X "Mass fractions";
        output Density d "Density";
      algorithm
        d := density(setState_phX(
                p,
                h,
                X));
      end density_phX;

      replaceable function temperature_psX
      "Return temperature from p,s, and X or Xi"
        extends Modelica.Icons.Function;
        input AbsolutePressure p "Pressure";
        input SpecificEntropy s "Specific entropy";
        input MassFraction X[:]=reference_X "Mass fractions";
        output Temperature T "Temperature";
      algorithm
        T := temperature(setState_psX(
                p,
                s,
                X));
        annotation (inverse(s=specificEntropy_pTX(
                      p,
                      T,
                      X)));
      end temperature_psX;

      replaceable function density_psX "Return density from p, s, and X or Xi"
        extends Modelica.Icons.Function;
        input AbsolutePressure p "Pressure";
        input SpecificEntropy s "Specific entropy";
        input MassFraction X[:]=reference_X "Mass fractions";
        output Density d "Density";
      algorithm
        d := density(setState_psX(
                p,
                s,
                X));
      end density_psX;

      replaceable function specificEnthalpy_psX
      "Return specific enthalpy from p, s, and X or Xi"
        extends Modelica.Icons.Function;
        input AbsolutePressure p "Pressure";
        input SpecificEntropy s "Specific entropy";
        input MassFraction X[:]=reference_X "Mass fractions";
        output SpecificEnthalpy h "Specific enthalpy";
      algorithm
        h := specificEnthalpy(setState_psX(
                p,
                s,
                X));
      end specificEnthalpy_psX;

      type MassFlowRate = SI.MassFlowRate (
          quantity="MassFlowRate." + mediumName,
          min=-1.0e5,
          max=1.e5) "Type for mass flow rate with medium specific attributes";

      // Only for backwards compatibility to version 3.2 (
      // (do not use these definitions in new models, but use Modelica.Media.Interfaces.Choices instead)
      package Choices = Modelica.Media.Interfaces.Choices annotation (obsolete=
            "Use Modelica.Media.Interfaces.Choices");

      annotation (Documentation(info="<html>
<p>
<b>PartialMedium</b> is a package and contains all <b>declarations</b> for
a medium. This means that constants, models, and functions
are defined that every medium is supposed to support
(some of them are optional). A medium package
inherits from <b>PartialMedium</b> and provides the
equations for the medium. The details of this package
are described in
<a href=\"modelica://Modelica.Media.UsersGuide\">Modelica.Media.UsersGuide</a>.
</p>
</html>",   revisions="<html>

</html>"));
    end PartialMedium;

    partial package PartialPureSubstance
    "Base class for pure substances of one chemical substance"
      extends PartialMedium(final reducedX=true, final fixedX=true);

      replaceable function setState_pT
      "Return thermodynamic state from p and T"
        extends Modelica.Icons.Function;
        input AbsolutePressure p "Pressure";
        input Temperature T "Temperature";
        output ThermodynamicState state "Thermodynamic state record";
      algorithm
        state := setState_pTX(
                p,
                T,
                fill(0, 0));
      end setState_pT;

      replaceable function setState_ph
      "Return thermodynamic state from p and h"
        extends Modelica.Icons.Function;
        input AbsolutePressure p "Pressure";
        input SpecificEnthalpy h "Specific enthalpy";
        output ThermodynamicState state "Thermodynamic state record";
      algorithm
        state := setState_phX(
                p,
                h,
                fill(0, 0));
      end setState_ph;

      replaceable function setState_ps
      "Return thermodynamic state from p and s"
        extends Modelica.Icons.Function;
        input AbsolutePressure p "Pressure";
        input SpecificEntropy s "Specific entropy";
        output ThermodynamicState state "Thermodynamic state record";
      algorithm
        state := setState_psX(
                p,
                s,
                fill(0, 0));
      end setState_ps;

      replaceable function setState_dT
      "Return thermodynamic state from d and T"
        extends Modelica.Icons.Function;
        input Density d "Density";
        input Temperature T "Temperature";
        output ThermodynamicState state "Thermodynamic state record";
      algorithm
        state := setState_dTX(
                d,
                T,
                fill(0, 0));
      end setState_dT;

      replaceable function density_ph "Return density from p and h"
        extends Modelica.Icons.Function;
        input AbsolutePressure p "Pressure";
        input SpecificEnthalpy h "Specific enthalpy";
        output Density d "Density";
      algorithm
        d := density_phX(
                p,
                h,
                fill(0, 0));
      end density_ph;

      replaceable function temperature_ph "Return temperature from p and h"
        extends Modelica.Icons.Function;
        input AbsolutePressure p "Pressure";
        input SpecificEnthalpy h "Specific enthalpy";
        output Temperature T "Temperature";
      algorithm
        T := temperature_phX(
                p,
                h,
                fill(0, 0));
      end temperature_ph;

      replaceable function pressure_dT "Return pressure from d and T"
        extends Modelica.Icons.Function;
        input Density d "Density";
        input Temperature T "Temperature";
        output AbsolutePressure p "Pressure";
      algorithm
        p := pressure(setState_dTX(
                d,
                T,
                fill(0, 0)));
      end pressure_dT;

      replaceable function specificEnthalpy_dT
      "Return specific enthalpy from d and T"
        extends Modelica.Icons.Function;
        input Density d "Density";
        input Temperature T "Temperature";
        output SpecificEnthalpy h "Specific enthalpy";
      algorithm
        h := specificEnthalpy(setState_dTX(
                d,
                T,
                fill(0, 0)));
      end specificEnthalpy_dT;

      replaceable function specificEnthalpy_ps
      "Return specific enthalpy from p and s"
        extends Modelica.Icons.Function;
        input AbsolutePressure p "Pressure";
        input SpecificEntropy s "Specific entropy";
        output SpecificEnthalpy h "Specific enthalpy";
      algorithm
        h := specificEnthalpy_psX(
                p,
                s,
                fill(0, 0));
      end specificEnthalpy_ps;

      replaceable function temperature_ps "Return temperature from p and s"
        extends Modelica.Icons.Function;
        input AbsolutePressure p "Pressure";
        input SpecificEntropy s "Specific entropy";
        output Temperature T "Temperature";
      algorithm
        T := temperature_psX(
                p,
                s,
                fill(0, 0));
      end temperature_ps;

      replaceable function density_ps "Return density from p and s"
        extends Modelica.Icons.Function;
        input AbsolutePressure p "Pressure";
        input SpecificEntropy s "Specific entropy";
        output Density d "Density";
      algorithm
        d := density_psX(
                p,
                s,
                fill(0, 0));
      end density_ps;

      replaceable function specificEnthalpy_pT
      "Return specific enthalpy from p and T"
        extends Modelica.Icons.Function;
        input AbsolutePressure p "Pressure";
        input Temperature T "Temperature";
        output SpecificEnthalpy h "Specific enthalpy";
      algorithm
        h := specificEnthalpy_pTX(
                p,
                T,
                fill(0, 0));
      end specificEnthalpy_pT;

      replaceable function density_pT "Return density from p and T"
        extends Modelica.Icons.Function;
        input AbsolutePressure p "Pressure";
        input Temperature T "Temperature";
        output Density d "Density";
      algorithm
        d := density(setState_pTX(
                p,
                T,
                fill(0, 0)));
      end density_pT;

      redeclare replaceable partial model extends BaseProperties(final
          standardOrderComponents=true)
      end BaseProperties;
    end PartialPureSubstance;

    partial package PartialTwoPhaseMedium
    "Base class for two phase medium of one substance"
      extends PartialPureSubstance(redeclare record FluidConstants =
            Modelica.Media.Interfaces.Types.TwoPhase.FluidConstants);
      constant Boolean smoothModel=false
      "True if the (derived) model should not generate state events";
      constant Boolean onePhase=false
      "True if the (derived) model should never be called with two-phase inputs";

      constant FluidConstants[nS] fluidConstants "Constant data for the fluid";

      redeclare replaceable record extends ThermodynamicState
      "Thermodynamic state of two phase medium"
        FixedPhase phase(min=0, max=2)
        "Phase of the fluid: 1 for 1-phase, 2 for two-phase, 0 for not known, e.g., interactive use";
      end ThermodynamicState;

      redeclare replaceable partial model extends BaseProperties
      "Base properties (p, d, T, h, u, R, MM, sat) of two phase medium"
        SaturationProperties sat "Saturation properties at the medium pressure";
      end BaseProperties;

      replaceable partial function setDewState
      "Return the thermodynamic state on the dew line"
        extends Modelica.Icons.Function;
        input SaturationProperties sat "Saturation point";
        input FixedPhase phase(
          min=1,
          max=2) = 1 "Phase: default is one phase";
        output ThermodynamicState state "Complete thermodynamic state info";
      end setDewState;

      replaceable partial function setBubbleState
      "Return the thermodynamic state on the bubble line"
        extends Modelica.Icons.Function;
        input SaturationProperties sat "Saturation point";
        input FixedPhase phase(
          min=1,
          max=2) = 1 "Phase: default is one phase";
        output ThermodynamicState state "Complete thermodynamic state info";
      end setBubbleState;

      redeclare replaceable partial function extends setState_dTX
      "Return thermodynamic state as function of d, T and composition X or Xi"
        input FixedPhase phase=0
        "2 for two-phase, 1 for one-phase, 0 if not known";
      end setState_dTX;

      redeclare replaceable partial function extends setState_phX
      "Return thermodynamic state as function of p, h and composition X or Xi"
        input FixedPhase phase=0
        "2 for two-phase, 1 for one-phase, 0 if not known";
      end setState_phX;

      redeclare replaceable partial function extends setState_psX
      "Return thermodynamic state as function of p, s and composition X or Xi"
        input FixedPhase phase=0
        "2 for two-phase, 1 for one-phase, 0 if not known";
      end setState_psX;

      redeclare replaceable partial function extends setState_pTX
      "Return thermodynamic state as function of p, T and composition X or Xi"
        input FixedPhase phase=0
        "2 for two-phase, 1 for one-phase, 0 if not known";
      end setState_pTX;

      replaceable function setSat_T
      "Return saturation property record from temperature"
        extends Modelica.Icons.Function;
        input Temperature T "Temperature";
        output SaturationProperties sat "Saturation property record";
      algorithm
        sat.Tsat := T;
        sat.psat := saturationPressure(T);
      end setSat_T;

      replaceable function setSat_p
      "Return saturation property record from pressure"
        extends Modelica.Icons.Function;
        input AbsolutePressure p "Pressure";
        output SaturationProperties sat "Saturation property record";
      algorithm
        sat.psat := p;
        sat.Tsat := saturationTemperature(p);
      end setSat_p;

      replaceable partial function bubbleEnthalpy
      "Return bubble point specific enthalpy"
        extends Modelica.Icons.Function;
        input SaturationProperties sat "Saturation property record";
        output SI.SpecificEnthalpy hl "Boiling curve specific enthalpy";
      end bubbleEnthalpy;

      replaceable partial function dewEnthalpy
      "Return dew point specific enthalpy"
        extends Modelica.Icons.Function;
        input SaturationProperties sat "Saturation property record";
        output SI.SpecificEnthalpy hv "Dew curve specific enthalpy";
      end dewEnthalpy;

      replaceable partial function bubbleEntropy
      "Return bubble point specific entropy"
        extends Modelica.Icons.Function;
        input SaturationProperties sat "Saturation property record";
        output SI.SpecificEntropy sl "Boiling curve specific entropy";
      end bubbleEntropy;

      replaceable partial function dewEntropy
      "Return dew point specific entropy"
        extends Modelica.Icons.Function;
        input SaturationProperties sat "Saturation property record";
        output SI.SpecificEntropy sv "Dew curve specific entropy";
      end dewEntropy;

      replaceable partial function bubbleDensity "Return bubble point density"
        extends Modelica.Icons.Function;
        input SaturationProperties sat "Saturation property record";
        output Density dl "Boiling curve density";
      end bubbleDensity;

      replaceable partial function dewDensity "Return dew point density"
        extends Modelica.Icons.Function;
        input SaturationProperties sat "Saturation property record";
        output Density dv "Dew curve density";
      end dewDensity;

      replaceable partial function saturationPressure
      "Return saturation pressure"
        extends Modelica.Icons.Function;
        input Temperature T "Temperature";
        output AbsolutePressure p "Saturation pressure";
      end saturationPressure;

      replaceable partial function saturationTemperature
      "Return saturation temperature"
        extends Modelica.Icons.Function;
        input AbsolutePressure p "Pressure";
        output Temperature T "Saturation temperature";
      end saturationTemperature;

      replaceable function saturationPressure_sat
      "Return saturation temperature"
        extends Modelica.Icons.Function;
        input SaturationProperties sat "Saturation property record";
        output AbsolutePressure p "Saturation pressure";
      algorithm
        p := sat.psat;
      end saturationPressure_sat;

      replaceable function saturationTemperature_sat
      "Return saturation temperature"
        extends Modelica.Icons.Function;
        input SaturationProperties sat "Saturation property record";
        output Temperature T "Saturation temperature";
      algorithm
        T := sat.Tsat;
      end saturationTemperature_sat;

      replaceable partial function saturationTemperature_derp
      "Return derivative of saturation temperature w.r.t. pressure"
        extends Modelica.Icons.Function;
        input AbsolutePressure p "Pressure";
        output DerTemperatureByPressure dTp
        "Derivative of saturation temperature w.r.t. pressure";
      end saturationTemperature_derp;

      replaceable function saturationTemperature_derp_sat
      "Return derivative of saturation temperature w.r.t. pressure"
        extends Modelica.Icons.Function;
        input SaturationProperties sat "Saturation property record";
        output DerTemperatureByPressure dTp
        "Derivative of saturation temperature w.r.t. pressure";
      algorithm
        dTp := saturationTemperature_derp(sat.psat);
      end saturationTemperature_derp_sat;

      replaceable partial function surfaceTension
      "Return surface tension sigma in the two phase region"
        extends Modelica.Icons.Function;
        input SaturationProperties sat "Saturation property record";
        output SurfaceTension sigma
        "Surface tension sigma in the two phase region";
      end surfaceTension;

      redeclare replaceable function extends molarMass
      "Return the molar mass of the medium"
      algorithm
        MM := fluidConstants[1].molarMass;
      end molarMass;

      replaceable partial function dBubbleDensity_dPressure
      "Return bubble point density derivative"
        extends Modelica.Icons.Function;
        input SaturationProperties sat "Saturation property record";
        output DerDensityByPressure ddldp "Boiling curve density derivative";
      end dBubbleDensity_dPressure;

      replaceable partial function dDewDensity_dPressure
      "Return dew point density derivative"
        extends Modelica.Icons.Function;
        input SaturationProperties sat "Saturation property record";
        output DerDensityByPressure ddvdp "Saturated steam density derivative";
      end dDewDensity_dPressure;

      replaceable partial function dBubbleEnthalpy_dPressure
      "Return bubble point specific enthalpy derivative"
        extends Modelica.Icons.Function;
        input SaturationProperties sat "Saturation property record";
        output DerEnthalpyByPressure dhldp
        "Boiling curve specific enthalpy derivative";
      end dBubbleEnthalpy_dPressure;

      replaceable partial function dDewEnthalpy_dPressure
      "Return dew point specific enthalpy derivative"
        extends Modelica.Icons.Function;

        input SaturationProperties sat "Saturation property record";
        output DerEnthalpyByPressure dhvdp
        "Saturated steam specific enthalpy derivative";
      end dDewEnthalpy_dPressure;

      redeclare replaceable function specificEnthalpy_pTX
      "Return specific enthalpy from pressure, temperature and mass fraction"
        extends Modelica.Icons.Function;
        input AbsolutePressure p "Pressure";
        input Temperature T "Temperature";
        input MassFraction X[:] "Mass fractions";
        input FixedPhase phase=0
        "2 for two-phase, 1 for one-phase, 0 if not known";
        output SpecificEnthalpy h "Specific enthalpy at p, T, X";
      algorithm
        h := specificEnthalpy(setState_pTX(
                p,
                T,
                X,
                phase));
      end specificEnthalpy_pTX;

      redeclare replaceable function temperature_phX
      "Return temperature from p, h, and X or Xi"
        extends Modelica.Icons.Function;
        input AbsolutePressure p "Pressure";
        input SpecificEnthalpy h "Specific enthalpy";
        input MassFraction X[:] "Mass fractions";
        input FixedPhase phase=0
        "2 for two-phase, 1 for one-phase, 0 if not known";
        output Temperature T "Temperature";
      algorithm
        T := temperature(setState_phX(
                p,
                h,
                X,
                phase));
      end temperature_phX;

      redeclare replaceable function density_phX
      "Return density from p, h, and X or Xi"
        extends Modelica.Icons.Function;
        input AbsolutePressure p "Pressure";
        input SpecificEnthalpy h "Specific enthalpy";
        input MassFraction X[:] "Mass fractions";
        input FixedPhase phase=0
        "2 for two-phase, 1 for one-phase, 0 if not known";
        output Density d "Density";
      algorithm
        d := density(setState_phX(
                p,
                h,
                X,
                phase));
      end density_phX;

      redeclare replaceable function temperature_psX
      "Return temperature from p, s, and X or Xi"
        extends Modelica.Icons.Function;
        input AbsolutePressure p "Pressure";
        input SpecificEntropy s "Specific entropy";
        input MassFraction X[:] "Mass fractions";
        input FixedPhase phase=0
        "2 for two-phase, 1 for one-phase, 0 if not known";
        output Temperature T "Temperature";
      algorithm
        T := temperature(setState_psX(
                p,
                s,
                X,
                phase));
      end temperature_psX;

      redeclare replaceable function density_psX
      "Return density from p, s, and X or Xi"
        extends Modelica.Icons.Function;
        input AbsolutePressure p "Pressure";
        input SpecificEntropy s "Specific entropy";
        input MassFraction X[:] "Mass fractions";
        input FixedPhase phase=0
        "2 for two-phase, 1 for one-phase, 0 if not known";
        output Density d "Density";
      algorithm
        d := density(setState_psX(
                p,
                s,
                X,
                phase));
      end density_psX;

      redeclare replaceable function specificEnthalpy_psX
      "Return specific enthalpy from p, s, and X or Xi"
        extends Modelica.Icons.Function;
        input AbsolutePressure p "Pressure";
        input SpecificEntropy s "Specific entropy";
        input MassFraction X[:] "Mass fractions";
        input FixedPhase phase=0
        "2 for two-phase, 1 for one-phase, 0 if not known";
        output SpecificEnthalpy h "Specific enthalpy";
      algorithm
        h := specificEnthalpy(setState_psX(
                p,
                s,
                X,
                phase));
      end specificEnthalpy_psX;

      redeclare replaceable function setState_pT
      "Return thermodynamic state from p and T"
        extends Modelica.Icons.Function;
        input AbsolutePressure p "Pressure";
        input Temperature T "Temperature";
        input FixedPhase phase=0
        "2 for two-phase, 1 for one-phase, 0 if not known";
        output ThermodynamicState state "Thermodynamic state record";
      algorithm
        state := setState_pTX(
                p,
                T,
                fill(0, 0),
                phase);
      end setState_pT;

      redeclare replaceable function setState_ph
      "Return thermodynamic state from p and h"
        extends Modelica.Icons.Function;
        input AbsolutePressure p "Pressure";
        input SpecificEnthalpy h "Specific enthalpy";
        input FixedPhase phase=0
        "2 for two-phase, 1 for one-phase, 0 if not known";
        output ThermodynamicState state "Thermodynamic state record";
      algorithm
        state := setState_phX(
                p,
                h,
                fill(0, 0),
                phase);
      end setState_ph;

      redeclare replaceable function setState_ps
      "Return thermodynamic state from p and s"
        extends Modelica.Icons.Function;
        input AbsolutePressure p "Pressure";
        input SpecificEntropy s "Specific entropy";
        input FixedPhase phase=0
        "2 for two-phase, 1 for one-phase, 0 if not known";
        output ThermodynamicState state "Thermodynamic state record";
      algorithm
        state := setState_psX(
                p,
                s,
                fill(0, 0),
                phase);
      end setState_ps;

      redeclare replaceable function setState_dT
      "Return thermodynamic state from d and T"
        extends Modelica.Icons.Function;
        input Density d "Density";
        input Temperature T "Temperature";
        input FixedPhase phase=0
        "2 for two-phase, 1 for one-phase, 0 if not known";
        output ThermodynamicState state "Thermodynamic state record";
      algorithm
        state := setState_dTX(
                d,
                T,
                fill(0, 0),
                phase);
      end setState_dT;

      replaceable function setState_px
      "Return thermodynamic state from pressure and vapour quality"
        extends Modelica.Icons.Function;
        input AbsolutePressure p "Pressure";
        input MassFraction x "Vapour quality";
        output ThermodynamicState state "Thermodynamic state record";
      algorithm
        state := setState_ph(
                p,
                (1 - x)*bubbleEnthalpy(setSat_p(p)) + x*dewEnthalpy(setSat_p(p)),
                2);
      end setState_px;

      replaceable function setState_Tx
      "Return thermodynamic state from temperature and vapour quality"
        extends Modelica.Icons.Function;
        input Temperature T "Temperature";
        input MassFraction x "Vapour quality";
        output ThermodynamicState state "Thermodynamic state record";
      algorithm
        state := setState_ph(
                saturationPressure_sat(setSat_T(T)),
                (1 - x)*bubbleEnthalpy(setSat_T(T)) + x*dewEnthalpy(setSat_T(T)),
                2);
      end setState_Tx;

      replaceable function vapourQuality "Return vapour quality"
        extends Modelica.Icons.Function;
        input ThermodynamicState state "Thermodynamic state record";
        output MassFraction x "Vapour quality";
    protected
        constant SpecificEnthalpy eps=1e-8;
      algorithm
        x := min(max((specificEnthalpy(state) - bubbleEnthalpy(setSat_p(pressure(
          state))))/(dewEnthalpy(setSat_p(pressure(state))) - bubbleEnthalpy(
          setSat_p(pressure(state))) + eps), 0), 1);
      end vapourQuality;

      redeclare replaceable function density_ph "Return density from p and h"
        extends Modelica.Icons.Function;
        input AbsolutePressure p "Pressure";
        input SpecificEnthalpy h "Specific enthalpy";
        input FixedPhase phase=0
        "2 for two-phase, 1 for one-phase, 0 if not known";
        output Density d "Density";
      algorithm
        d := density_phX(
                p,
                h,
                fill(0, 0),
                phase);
      end density_ph;

      redeclare replaceable function temperature_ph
      "Return temperature from p and h"
        extends Modelica.Icons.Function;
        input AbsolutePressure p "Pressure";
        input SpecificEnthalpy h "Specific enthalpy";
        input FixedPhase phase=0
        "2 for two-phase, 1 for one-phase, 0 if not known";
        output Temperature T "Temperature";
      algorithm
        T := temperature_phX(
                p,
                h,
                fill(0, 0),
                phase);
      end temperature_ph;

      redeclare replaceable function pressure_dT "Return pressure from d and T"
        extends Modelica.Icons.Function;
        input Density d "Density";
        input Temperature T "Temperature";
        input FixedPhase phase=0
        "2 for two-phase, 1 for one-phase, 0 if not known";
        output AbsolutePressure p "Pressure";
      algorithm
        p := pressure(setState_dTX(
                d,
                T,
                fill(0, 0),
                phase));
      end pressure_dT;

      redeclare replaceable function specificEnthalpy_dT
      "Return specific enthalpy from d and T"
        extends Modelica.Icons.Function;
        input Density d "Density";
        input Temperature T "Temperature";
        input FixedPhase phase=0
        "2 for two-phase, 1 for one-phase, 0 if not known";
        output SpecificEnthalpy h "Specific enthalpy";
      algorithm
        h := specificEnthalpy(setState_dTX(
                d,
                T,
                fill(0, 0),
                phase));
      end specificEnthalpy_dT;

      redeclare replaceable function specificEnthalpy_ps
      "Return specific enthalpy from p and s"
        extends Modelica.Icons.Function;
        input AbsolutePressure p "Pressure";
        input SpecificEntropy s "Specific entropy";
        input FixedPhase phase=0
        "2 for two-phase, 1 for one-phase, 0 if not known";
        output SpecificEnthalpy h "Specific enthalpy";
      algorithm
        h := specificEnthalpy_psX(
                p,
                s,
                fill(0, 0));
      end specificEnthalpy_ps;

      redeclare replaceable function temperature_ps
      "Return temperature from p and s"
        extends Modelica.Icons.Function;
        input AbsolutePressure p "Pressure";
        input SpecificEntropy s "Specific entropy";
        input FixedPhase phase=0
        "2 for two-phase, 1 for one-phase, 0 if not known";
        output Temperature T "Temperature";
      algorithm
        T := temperature_psX(
                p,
                s,
                fill(0, 0),
                phase);
      end temperature_ps;

      redeclare replaceable function density_ps "Return density from p and s"
        extends Modelica.Icons.Function;
        input AbsolutePressure p "Pressure";
        input SpecificEntropy s "Specific entropy";
        input FixedPhase phase=0
        "2 for two-phase, 1 for one-phase, 0 if not known";
        output Density d "Density";
      algorithm
        d := density_psX(
                p,
                s,
                fill(0, 0),
                phase);
      end density_ps;

      redeclare replaceable function specificEnthalpy_pT
      "Return specific enthalpy from p and T"
        extends Modelica.Icons.Function;
        input AbsolutePressure p "Pressure";
        input Temperature T "Temperature";
        input FixedPhase phase=0
        "2 for two-phase, 1 for one-phase, 0 if not known";
        output SpecificEnthalpy h "Specific enthalpy";
      algorithm
        h := specificEnthalpy_pTX(
                p,
                T,
                fill(0, 0),
                phase);
      end specificEnthalpy_pT;

      redeclare replaceable function density_pT "Return density from p and T"
        extends Modelica.Icons.Function;
        input AbsolutePressure p "Pressure";
        input Temperature T "Temperature";
        input FixedPhase phase=0
        "2 for two-phase, 1 for one-phase, 0 if not known";
        output Density d "Density";
      algorithm
        d := density(setState_pTX(
                p,
                T,
                fill(0, 0),
                phase));
      end density_pT;
    end PartialTwoPhaseMedium;

    package Choices "Types, constants to define menu choices"
      extends Modelica.Icons.Package;

      type IndependentVariables = enumeration(
        T "Temperature",
        pT "Pressure, Temperature",
        ph "Pressure, Specific Enthalpy",
        phX "Pressure, Specific Enthalpy, Mass Fraction",
        pTX "Pressure, Temperature, Mass Fractions",
        dTX "Density, Temperature, Mass Fractions")
      "Enumeration defining the independent variables of a medium";
      annotation (Documentation(info="<html>
<p>
Enumerations and data types for all types of fluids
</p>

<p>
Note: Reference enthalpy might have to be extended with enthalpy of formation.
</p>
</html>"));
    end Choices;

    package Types "Types to be used in fluid models"
      extends Modelica.Icons.Package;

      type AbsolutePressure = SI.AbsolutePressure (
          min=0,
          max=1.e8,
          nominal=1.e5,
          start=1.e5)
      "Type for absolute pressure with medium specific attributes";

      type Density = SI.Density (
          min=0,
          max=1.e5,
          nominal=1,
          start=1) "Type for density with medium specific attributes";

      type DynamicViscosity = SI.DynamicViscosity (
          min=0,
          max=1.e8,
          nominal=1.e-3,
          start=1.e-3)
      "Type for dynamic viscosity with medium specific attributes";

      type MassFraction = Real (
          quantity="MassFraction",
          final unit="kg/kg",
          min=0,
          max=1,
          nominal=0.1) "Type for mass fraction with medium specific attributes";

      type MolarMass = SI.MolarMass (
          min=0.001,
          max=0.25,
          nominal=0.032) "Type for molar mass with medium specific attributes";

      type MolarVolume = SI.MolarVolume (
          min=1e-6,
          max=1.0e6,
          nominal=1.0) "Type for molar volume with medium specific attributes";

      type IsentropicExponent = SI.RatioOfSpecificHeatCapacities (
          min=1,
          max=500000,
          nominal=1.2,
          start=1.2)
      "Type for isentropic exponent with medium specific attributes";

      type SpecificEnergy = SI.SpecificEnergy (
          min=-1.0e8,
          max=1.e8,
          nominal=1.e6)
      "Type for specific energy with medium specific attributes";

      type SpecificInternalEnergy = SpecificEnergy
      "Type for specific internal energy with medium specific attributes";

      type SpecificEnthalpy = SI.SpecificEnthalpy (
          min=-1.0e10,
          max=1.e10,
          nominal=1.e6)
      "Type for specific enthalpy with medium specific attributes";

      type SpecificEntropy = SI.SpecificEntropy (
          min=-1.e7,
          max=1.e7,
          nominal=1.e3)
      "Type for specific entropy with medium specific attributes";

      type SpecificHeatCapacity = SI.SpecificHeatCapacity (
          min=0,
          max=1.e7,
          nominal=1.e3,
          start=1.e3)
      "Type for specific heat capacity with medium specific attributes";

      type SurfaceTension = SI.SurfaceTension
      "Type for surface tension with medium specific attributes";

      type Temperature = SI.Temperature (
          min=1,
          max=1.e4,
          nominal=300,
          start=300) "Type for temperature with medium specific attributes";

      type ThermalConductivity = SI.ThermalConductivity (
          min=0,
          max=500,
          nominal=1,
          start=1)
      "Type for thermal conductivity with medium specific attributes";

      type PrandtlNumber = SI.PrandtlNumber (
          min=1e-3,
          max=1e5,
          nominal=1.0)
      "Type for Prandtl number with medium specific attributes";

      type VelocityOfSound = SI.Velocity (
          min=0,
          max=1.e5,
          nominal=1000,
          start=1000)
      "Type for velocity of sound with medium specific attributes";

      type IsobaricExpansionCoefficient = Real (
          min=0,
          max=1.0e8,
          unit="1/K")
      "Type for isobaric expansion coefficient with medium specific attributes";

      type DipoleMoment = Real (
          min=0.0,
          max=2.0,
          unit="debye",
          quantity="ElectricDipoleMoment")
      "Type for dipole moment with medium specific attributes";

      type DerDensityByPressure = SI.DerDensityByPressure
      "Type for partial derivative of density with respect to pressure with medium specific attributes";

      type DerDensityByEnthalpy = SI.DerDensityByEnthalpy
      "Type for partial derivative of density with respect to enthalpy with medium specific attributes";

      type DerEnthalpyByPressure = SI.DerEnthalpyByPressure
      "Type for partial derivative of enthalpy with respect to pressure with medium specific attributes";

      type DerDensityByTemperature = SI.DerDensityByTemperature
      "Type for partial derivative of density with respect to temperature with medium specific attributes";

      type DerTemperatureByPressure = Real (final unit="K/Pa")
      "Type for partial derivative of temperature with respect to pressure with medium specific attributes";

      replaceable record SaturationProperties
      "Saturation properties of two phase medium"
        extends Modelica.Icons.Record;
        AbsolutePressure psat "Saturation pressure";
        Temperature Tsat "Saturation temperature";
      end SaturationProperties;

      type FixedPhase = Integer (min=0, max=2)
      "Phase of the fluid: 1 for 1-phase, 2 for two-phase, 0 for not known, e.g., interactive use";

      package Basic
      "The most basic version of a record used in several degrees of detail"
        extends Icons.Package;

        record FluidConstants
        "Critical, triple, molecular and other standard data of fluid"
          extends Modelica.Icons.Record;
          String iupacName
          "Complete IUPAC name (or common name, if non-existent)";
          String casRegistryNumber
          "Chemical abstracts sequencing number (if it exists)";
          String chemicalFormula
          "Chemical formula, (brutto, nomenclature according to Hill";
          String structureFormula "Chemical structure formula";
          MolarMass molarMass "Molar mass";
        end FluidConstants;
      end Basic;

      package TwoPhase
      "The two phase fluid version of a record used in several degrees of detail"
        extends Icons.Package;

        record FluidConstants "Extended fluid constants"
          extends Modelica.Media.Interfaces.Types.Basic.FluidConstants;
          Temperature criticalTemperature "Critical temperature";
          AbsolutePressure criticalPressure "Critical pressure";
          MolarVolume criticalMolarVolume "Critical molar Volume";
          Real acentricFactor "Pitzer acentric factor";
          Temperature triplePointTemperature "Triple point temperature";
          AbsolutePressure triplePointPressure "Triple point pressure";
          Temperature meltingPoint "Melting point at 101325 Pa";
          Temperature normalBoilingPoint "Normal boiling point (at 101325 Pa)";
          DipoleMoment dipoleMoment
          "Dipole moment of molecule in Debye (1 debye = 3.33564e10-30 C.m)";
          Boolean hasIdealGasHeatCapacity=false
          "True if ideal gas heat capacity is available";
          Boolean hasCriticalData=false "True if critical data are known";
          Boolean hasDipoleMoment=false "True if a dipole moment known";
          Boolean hasFundamentalEquation=false "True if a fundamental equation";
          Boolean hasLiquidHeatCapacity=false
          "True if liquid heat capacity is available";
          Boolean hasSolidHeatCapacity=false
          "True if solid heat capacity is available";
          Boolean hasAccurateViscosityData=false
          "True if accurate data for a viscosity function is available";
          Boolean hasAccurateConductivityData=false
          "True if accurate data for thermal conductivity is available";
          Boolean hasVapourPressureCurve=false
          "True if vapour pressure data, e.g., Antoine coefficents are known";
          Boolean hasAcentricFactor=false
          "True if Pitzer accentric factor is known";
          SpecificEnthalpy HCRIT0=0.0
          "Critical specific enthalpy of the fundamental equation";
          SpecificEntropy SCRIT0=0.0
          "Critical specific entropy of the fundamental equation";
          SpecificEnthalpy deltah=0.0
          "Difference between specific enthalpy model (h_m) and f.eq. (h_f) (h_m - h_f)";
          SpecificEntropy deltas=0.0
          "Difference between specific enthalpy model (s_m) and f.eq. (s_f) (s_m - s_f)";
        end FluidConstants;
      end TwoPhase;
    end Types;
    annotation (Documentation(info="<HTML>
<p>
This package provides basic interfaces definitions of media models for different
kind of media.
</p>
</HTML>"));
  end Interfaces;
  annotation (preferredView="info",Documentation(info="<HTML>
<p>
This library contains <a href=\"modelica://Modelica.Media.Interfaces\">interface</a>
definitions for media and the following <b>property</b> models for
single and multiple substance fluids with one and multiple phases:
</p>
<ul>
<li> <a href=\"modelica://Modelica.Media.IdealGases\">Ideal gases:</a><br>
     1241 high precision gas models based on the
     NASA Glenn coefficients, plus ideal gas mixture models based
     on the same data.</li>
<li> <a href=\"modelica://Modelica.Media.Water\">Water models:</a><br>
     ConstantPropertyLiquidWater, WaterIF97 (high precision
     water model according to the IAPWS/IF97 standard)</li>
<li> <a href=\"modelica://Modelica.Media.Air\">Air models:</a><br>
     SimpleAir, DryAirNasa, ReferenceAir, MoistAir, ReferenceMoistAir.</li>
<li> <a href=\"modelica://Modelica.Media.Incompressible\">
     Incompressible media:</a><br>
     TableBased incompressible fluid models (properties are defined by tables rho(T),
     HeatCapacity_cp(T), etc.)</li>
<li> <a href=\"modelica://Modelica.Media.CompressibleLiquids\">
     Compressible liquids:</a><br>
     Simple liquid models with linear compressibility</li>
<li> <a href=\"modelica://Modelica.Media.R134a\">Refrigerant Tetrafluoroethane (R134a)</a>.</li>
</ul>
<p>
The following parts are useful, when newly starting with this library:
<ul>
<li> <a href=\"modelica://Modelica.Media.UsersGuide\">Modelica.Media.UsersGuide</a>.</li>
<li> <a href=\"modelica://Modelica.Media.UsersGuide.MediumUsage\">Modelica.Media.UsersGuide.MediumUsage</a>
     describes how to use a medium model in a component model.</li>
<li> <a href=\"modelica://Modelica.Media.UsersGuide.MediumDefinition\">
     Modelica.Media.UsersGuide.MediumDefinition</a>
     describes how a new fluid medium model has to be implemented.</li>
<li> <a href=\"modelica://Modelica.Media.UsersGuide.ReleaseNotes\">Modelica.Media.UsersGuide.ReleaseNotes</a>
     summarizes the changes of the library releases.</li>
<li> <a href=\"modelica://Modelica.Media.Examples\">Modelica.Media.Examples</a>
     contains examples that demonstrate the usage of this library.</li>
</ul>
<p>
Copyright &copy; 1998-2013, Modelica Association.
</p>
<p>
<i>This Modelica package is <u>free</u> software and the use is completely at <u>your own risk</u>; it can be redistributed and/or modified under the terms of the Modelica License 2. For license conditions (including the disclaimer of warranty) see <a href=\"modelica://Modelica.UsersGuide.ModelicaLicense2\">Modelica.UsersGuide.ModelicaLicense2</a> or visit <a href=\"https://www.modelica.org/licenses/ModelicaLicense2\"> https://www.modelica.org/licenses/ModelicaLicense2</a>.</i>
</p>
</HTML>",   revisions="<html>
<ul>
<li><i>May 16, 2013</i> by Stefan Wischhusen (XRG Simulation):<br/>
    Added new media models Air.ReferenceMoistAir, Air.ReferenceAir, R134a.</li>
<li><i>May 25, 2011</i> by Francesco Casella:<br/>Added min/max attributes to Water, TableBased, MixtureGasNasa, SimpleAir and MoistAir local types.</li>
<li><i>May 25, 2011</i> by Stefan Wischhusen:<br/>Added individual settings for polynomial fittings of properties.</li>
</ul>
</html>"),
      Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,100}}),
          graphics={
          Line(
            points = {{-76,-80},{-62,-30},{-32,40},{4,66},{48,66},{73,45},{62,-8},{48,-50},{38,-80}},
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

  package Math
  "Library of mathematical functions (e.g., sin, cos) and of functions operating on vectors and matrices"
  import SI = Modelica.SIunits;
  extends Modelica.Icons.Package;

  package Icons "Icons for Math"
    extends Modelica.Icons.IconsPackage;

    partial function AxisLeft
    "Basic icon for mathematical function with y-axis on left side"

      annotation (
        Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},{100,
                100}}), graphics={
            Rectangle(
              extent={{-100,100},{100,-100}},
              lineColor={0,0,0},
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid),
            Line(points={{-80,-80},{-80,68}}, color={192,192,192}),
            Polygon(
              points={{-80,90},{-88,68},{-72,68},{-80,90}},
              lineColor={192,192,192},
              fillColor={192,192,192},
              fillPattern=FillPattern.Solid),
            Text(
              extent={{-150,150},{150,110}},
              textString="%name",
              lineColor={0,0,255})}),
        Diagram(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},{
                100,100}}), graphics={Line(points={{-80,80},{-88,80}}, color={95,
              95,95}),Line(points={{-80,-80},{-88,-80}}, color={95,95,95}),Line(
              points={{-80,-90},{-80,84}}, color={95,95,95}),Text(
                  extent={{-75,104},{-55,84}},
                  lineColor={95,95,95},
                  textString="y"),Polygon(
                  points={{-80,98},{-86,82},{-74,82},{-80,98}},
                  lineColor={95,95,95},
                  fillColor={95,95,95},
                  fillPattern=FillPattern.Solid)}),
        Documentation(info="<html>
<p>
Icon for a mathematical function, consisting of an y-axis on the left side.
It is expected, that an x-axis is added and a plot of the function.
</p>
</html>"));
    end AxisLeft;

    partial function AxisCenter
    "Basic icon for mathematical function with y-axis in the center"

      annotation (
        Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},{100,
                100}}), graphics={
            Rectangle(
              extent={{-100,100},{100,-100}},
              lineColor={0,0,0},
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid),
            Line(points={{0,-80},{0,68}}, color={192,192,192}),
            Polygon(
              points={{0,90},{-8,68},{8,68},{0,90}},
              lineColor={192,192,192},
              fillColor={192,192,192},
              fillPattern=FillPattern.Solid),
            Text(
              extent={{-150,150},{150,110}},
              textString="%name",
              lineColor={0,0,255})}),
        Diagram(graphics={Line(points={{0,80},{-8,80}}, color={95,95,95}),Line(
              points={{0,-80},{-8,-80}}, color={95,95,95}),Line(points={{0,-90},{
              0,84}}, color={95,95,95}),Text(
                  extent={{5,104},{25,84}},
                  lineColor={95,95,95},
                  textString="y"),Polygon(
                  points={{0,98},{-6,82},{6,82},{0,98}},
                  lineColor={95,95,95},
                  fillColor={95,95,95},
                  fillPattern=FillPattern.Solid)}),
        Documentation(info="<html>
<p>
Icon for a mathematical function, consisting of an y-axis in the middle.
It is expected, that an x-axis is added and a plot of the function.
</p>
</html>"));
    end AxisCenter;
  end Icons;

  function asin "Inverse sine (-1 <= u <= 1)"
    extends Modelica.Math.Icons.AxisCenter;
    input Real u;
    output SI.Angle y;

  external "builtin" y=  asin(u);
    annotation (
      Icon(coordinateSystem(
          preserveAspectRatio=true,
          extent={{-100,-100},{100,100}}), graphics={
          Line(points={{-90,0},{68,0}}, color={192,192,192}),
          Polygon(
            points={{90,0},{68,8},{68,-8},{90,0}},
            lineColor={192,192,192},
            fillColor={192,192,192},
            fillPattern=FillPattern.Solid),
          Line(points={{-80,-80},{-79.2,-72.8},{-77.6,-67.5},{-73.6,-59.4},{-66.3,
                -49.8},{-53.5,-37.3},{-30.2,-19.7},{37.4,24.8},{57.5,40.8},{68.7,
                52.7},{75.2,62.2},{77.6,67.5},{80,80}}, color={0,0,0}),
          Text(
            extent={{-88,78},{-16,30}},
            lineColor={192,192,192},
            textString="asin")}),
      Diagram(coordinateSystem(
          preserveAspectRatio=true,
          extent={{-100,-100},{100,100}}), graphics={Text(
              extent={{-40,-72},{-15,-88}},
              textString="-pi/2",
              lineColor={0,0,255}),Text(
              extent={{-38,88},{-13,72}},
              textString=" pi/2",
              lineColor={0,0,255}),Text(
              extent={{68,-9},{88,-29}},
              textString="+1",
              lineColor={0,0,255}),Text(
              extent={{-90,21},{-70,1}},
              textString="-1",
              lineColor={0,0,255}),Line(points={{-100,0},{84,0}}, color={95,95,95}),
            Polygon(
              points={{98,0},{82,6},{82,-6},{98,0}},
              lineColor={95,95,95},
              fillColor={95,95,95},
              fillPattern=FillPattern.Solid),Line(
              points={{-80,-80},{-79.2,-72.8},{-77.6,-67.5},{-73.6,-59.4},{-66.3,
              -49.8},{-53.5,-37.3},{-30.2,-19.7},{37.4,24.8},{57.5,40.8},{68.7,
              52.7},{75.2,62.2},{77.6,67.5},{80,80}},
              color={0,0,255},
              thickness=0.5),Text(
              extent={{82,24},{102,4}},
              lineColor={95,95,95},
              textString="u"),Line(
              points={{0,80},{86,80}},
              color={175,175,175},
              smooth=Smooth.None),Line(
              points={{80,86},{80,-10}},
              color={175,175,175},
              smooth=Smooth.None)}),
      Documentation(info="<html>
<p>
This function returns y = asin(u), with -1 &le; u &le; +1:
</p>

<p>
<img src=\"modelica://Modelica/Resources/Images/Math/asin.png\">
</p>
</html>"));
  end asin;

  function acos "Inverse cosine (-1 <= u <= 1)"
    extends Modelica.Math.Icons.AxisCenter;
    input Real u;
    output SI.Angle y;

  external "builtin" y=  acos(u);
    annotation (
      Icon(coordinateSystem(
          preserveAspectRatio=true,
          extent={{-100,-100},{100,100}}), graphics={
          Line(points={{-90,-80},{68,-80}}, color={192,192,192}),
          Polygon(
            points={{90,-80},{68,-72},{68,-88},{90,-80}},
            lineColor={192,192,192},
            fillColor={192,192,192},
            fillPattern=FillPattern.Solid),
          Line(points={{-80,80},{-79.2,72.8},{-77.6,67.5},{-73.6,59.4},{-66.3,
                49.8},{-53.5,37.3},{-30.2,19.7},{37.4,-24.8},{57.5,-40.8},{68.7,-52.7},
                {75.2,-62.2},{77.6,-67.5},{80,-80}}, color={0,0,0}),
          Text(
            extent={{-86,-14},{-14,-62}},
            lineColor={192,192,192},
            textString="acos")}),
      Diagram(coordinateSystem(
          preserveAspectRatio=true,
          extent={{-100,-100},{100,100}}), graphics={Line(points={{-100,-80},{84,-80}}, color={95,95,
            95}),Polygon(
              points={{98,-80},{82,-74},{82,-86},{98,-80}},
              lineColor={95,95,95},
              fillColor={95,95,95},
              fillPattern=FillPattern.Solid),Line(
              points={{-80,80},{-79.2,72.8},{-77.6,67.5},{-73.6,59.4},{-66.3,49.8},
              {-53.5,37.3},{-30.2,19.7},{37.4,-24.8},{57.5,-40.8},{68.7,-52.7},{
              75.2,-62.2},{77.6,-67.5},{80,-80}},
              color={0,0,255},
              thickness=0.5),Text(
              extent={{-30,88},{-5,72}},
              textString=" pi",
              lineColor={0,0,255}),Text(
              extent={{-94,-57},{-74,-77}},
              textString="-1",
              lineColor={0,0,255}),Text(
              extent={{60,-81},{80,-101}},
              textString="+1",
              lineColor={0,0,255}),Text(
              extent={{82,-56},{102,-76}},
              lineColor={95,95,95},
              textString="u"),Line(
              points={{-2,80},{84,80}},
              color={175,175,175},
              smooth=Smooth.None),Line(
              points={{80,82},{80,-86}},
              color={175,175,175},
              smooth=Smooth.None)}),
      Documentation(info="<html>
<p>
This function returns y = acos(u), with -1 &le; u &le; +1:
</p>

<p>
<img src=\"modelica://Modelica/Resources/Images/Math/acos.png\">
</p>
</html>"));
  end acos;

  function log "Natural (base e) logarithm (u shall be > 0)"
    extends Modelica.Math.Icons.AxisLeft;
    input Real u;
    output Real y;

  external "builtin" y=  log(u);
    annotation (
      Icon(coordinateSystem(
          preserveAspectRatio=true,
          extent={{-100,-100},{100,100}}), graphics={
          Line(points={{-90,0},{68,0}}, color={192,192,192}),
          Polygon(
            points={{90,0},{68,8},{68,-8},{90,0}},
            lineColor={192,192,192},
            fillColor={192,192,192},
            fillPattern=FillPattern.Solid),
          Line(points={{-80,-80},{-79.2,-50.6},{-78.4,-37},{-77.6,-28},{-76.8,-21.3},
                {-75.2,-11.4},{-72.8,-1.31},{-69.5,8.08},{-64.7,17.9},{-57.5,28},
                {-47,38.1},{-31.8,48.1},{-10.1,58},{22.1,68},{68.7,78.1},{80,80}},
              color={0,0,0}),
          Text(
            extent={{-6,-24},{66,-72}},
            lineColor={192,192,192},
            textString="log")}),
      Diagram(coordinateSystem(
          preserveAspectRatio=true,
          extent={{-100,-100},{100,100}}), graphics={Line(points={{-100,0},{84,0}}, color={95,95,95}),
            Polygon(
              points={{100,0},{84,6},{84,-6},{100,0}},
              lineColor={95,95,95},
              fillColor={95,95,95},
              fillPattern=FillPattern.Solid),Line(
              points={{-78,-80},{-77.2,-50.6},{-76.4,-37},{-75.6,-28},{-74.8,-21.3},
              {-73.2,-11.4},{-70.8,-1.31},{-67.5,8.08},{-62.7,17.9},{-55.5,28},{-45,
              38.1},{-29.8,48.1},{-8.1,58},{24.1,68},{70.7,78.1},{82,80}},
              color={0,0,255},
              thickness=0.5),Text(
              extent={{-105,72},{-85,88}},
              textString="3",
              lineColor={0,0,255}),Text(
              extent={{60,-3},{80,-23}},
              textString="20",
              lineColor={0,0,255}),Text(
              extent={{-78,-7},{-58,-27}},
              textString="1",
              lineColor={0,0,255}),Text(
              extent={{84,26},{104,6}},
              lineColor={95,95,95},
              textString="u"),Text(
              extent={{-100,9},{-80,-11}},
              textString="0",
              lineColor={0,0,255}),Line(
              points={{-80,80},{84,80}},
              color={175,175,175},
              smooth=Smooth.None),Line(
              points={{82,82},{82,-6}},
              color={175,175,175},
              smooth=Smooth.None)}),
      Documentation(info="<html>
<p>
This function returns y = log(10) (the natural logarithm of u),
with u &gt; 0:
</p>

<p>
<img src=\"modelica://Modelica/Resources/Images/Math/log.png\">
</p>
</html>"));
  end log;

  function log10 "Base 10 logarithm (u shall be > 0)"
    extends Modelica.Math.Icons.AxisLeft;
    input Real u;
    output Real y;

  external "builtin" y=  log10(u);
    annotation (
      Icon(coordinateSystem(
          preserveAspectRatio=true,
          extent={{-100,-100},{100,100}}), graphics={
          Line(points={{-90,0},{68,0}}, color={192,192,192}),
          Polygon(
            points={{90,0},{68,8},{68,-8},{90,0}},
            lineColor={192,192,192},
            fillColor={192,192,192},
            fillPattern=FillPattern.Solid),
          Line(points={{-79.8,-80},{-79.2,-50.6},{-78.4,-37},{-77.6,-28},{-76.8,-21.3},
                {-75.2,-11.4},{-72.8,-1.31},{-69.5,8.08},{-64.7,17.9},{-57.5,28},
                {-47,38.1},{-31.8,48.1},{-10.1,58},{22.1,68},{68.7,78.1},{80,80}},
              color={0,0,0}),
          Text(
            extent={{-30,-22},{60,-70}},
            lineColor={192,192,192},
            textString="log10")}),
      Diagram(coordinateSystem(
          preserveAspectRatio=true,
          extent={{-100,-100},{100,100}}), graphics={Line(points={{-100,0},{84,0}}, color={95,95,95}),
            Polygon(
              points={{98,0},{82,6},{82,-6},{98,0}},
              lineColor={95,95,95},
              fillColor={95,95,95},
              fillPattern=FillPattern.Solid),Line(
              points={{-77.8,-80},{-77.2,-50.6},{-76.4,-37},{-75.6,-28},{-74.8,-21.3},
              {-73.2,-11.4},{-70.8,-1.31},{-67.5,8.08},{-62.7,17.9},{-55.5,28},{-45,
              38.1},{-29.8,48.1},{-8.1,58},{24.1,68},{70.7,78.1},{82,80}},
              color={0,0,255},
              thickness=0.5),Text(
              extent={{66,-13},{86,-33}},
              textString="20",
              lineColor={0,0,255}),Text(
              extent={{-78,-1},{-58,-21}},
              textString="1",
              lineColor={0,0,255}),Text(
              extent={{-83,62},{-63,78}},
              textString=" 1.3",
              lineColor={0,0,255}),Text(
              extent={{80,24},{100,4}},
              lineColor={95,95,95},
              textString="u"),Text(
              extent={{-100,9},{-80,-11}},
              textString="0",
              lineColor={0,0,255}),Line(
              points={{-80,80},{86,80}},
              color={175,175,175},
              smooth=Smooth.None),Line(
              points={{80,92},{80,-12}},
              color={175,175,175},
              smooth=Smooth.None)}),
      Documentation(info="<html>
<p>
This function returns y = log10(u),
with u &gt; 0:
</p>

<p>
<img src=\"modelica://Modelica/Resources/Images/Math/log10.png\">
</p>
</html>"));
  end log10;
  annotation (Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},
            {100,100}}), graphics={Line(points={{-80,0},{-68.7,34.2},{-61.5,53.1},
              {-55.1,66.4},{-49.4,74.6},{-43.8,79.1},{-38.2,79.8},{-32.6,76.6},{
              -26.9,69.7},{-21.3,59.4},{-14.9,44.1},{-6.83,21.2},{10.1,-30.8},{17.3,
              -50.2},{23.7,-64.2},{29.3,-73.1},{35,-78.4},{40.6,-80},{46.2,-77.6},
              {51.9,-71.5},{57.5,-61.9},{63.9,-47.2},{72,-24.8},{80,0}}, color={
              0,0,0}, smooth=Smooth.Bezier)}), Documentation(info="<HTML>
<p>
This package contains <b>basic mathematical functions</b> (such as sin(..)),
as well as functions operating on
<a href=\"modelica://Modelica.Math.Vectors\">vectors</a>,
<a href=\"modelica://Modelica.Math.Matrices\">matrices</a>,
<a href=\"modelica://Modelica.Math.Nonlinear\">nonlinear functions</a>, and
<a href=\"modelica://Modelica.Math.BooleanVectors\">Boolean vectors</a>.
</p>

<dl>
<dt><b>Main Authors:</b>
<dd><a href=\"http://www.robotic.dlr.de/Martin.Otter/\">Martin Otter</a> and
    Marcus Baur<br>
    Deutsches Zentrum f&uuml;r Luft und Raumfahrt e.V. (DLR)<br>
    Institut f&uuml;r Robotik und Mechatronik<br>
    Postfach 1116<br>
    D-82230 Wessling<br>
    Germany<br>
    email: <A HREF=\"mailto:Martin.Otter@dlr.de\">Martin.Otter@dlr.de</A><br>
</dl>

<p>
Copyright &copy; 1998-2013, Modelica Association and DLR.
</p>
<p>
<i>This Modelica package is <u>free</u> software and the use is completely at <u>your own risk</u>; it can be redistributed and/or modified under the terms of the Modelica License 2. For license conditions (including the disclaimer of warranty) see <a href=\"modelica://Modelica.UsersGuide.ModelicaLicense2\">Modelica.UsersGuide.ModelicaLicense2</a> or visit <a href=\"https://www.modelica.org/licenses/ModelicaLicense2\"> https://www.modelica.org/licenses/ModelicaLicense2</a>.</i>
</p>
</html>",   revisions="<html>
<ul>
<li><i>October 21, 2002</i>
       by <a href=\"http://www.robotic.dlr.de/Martin.Otter/\">Martin Otter</a>
       and <a href=\"http://www.robotic.dlr.de/Christian.Schweiger/\">Christian Schweiger</a>:<br>
       Function tempInterpol2 added.</li>
<li><i>Oct. 24, 1999</i>
       by <a href=\"http://www.robotic.dlr.de/Martin.Otter/\">Martin Otter</a>:<br>
       Icons for icon and diagram level introduced.</li>
<li><i>June 30, 1999</i>
       by <a href=\"http://www.robotic.dlr.de/Martin.Otter/\">Martin Otter</a>:<br>
       Realized.</li>
</ul>

</html>"));
  end Math;

  package Utilities
  "Library of utility functions dedicated to scripting (operating on files, streams, strings, system)"
    extends Modelica.Icons.Package;

    package Streams "Read from files and write to files"
      extends Modelica.Icons.Package;

      function print "Print string to terminal or file"
        extends Modelica.Icons.Function;
        input String string="" "String to be printed";
        input String fileName=""
        "File where to print (empty string is the terminal)"
                     annotation(Dialog(saveSelector(filter="Text files (*.txt)",
                            caption="Text file to store the output of print(..)")));
      external "C" ModelicaInternal_print(string, fileName) annotation(Library="ModelicaExternalC");
        annotation (Documentation(info="<HTML>
<h4>Syntax</h4>
<blockquote><pre>
Streams.<b>print</b>(string);
Streams.<b>print</b>(string,fileName);
</pre></blockquote>
<h4>Description</h4>
<p>
Function <b>print</b>(..) opens automatically the given file, if
it is not yet open. If the file does not exist, it is created.
If the file does exist, the given string is appended to the file.
If this is not desired, call \"Files.remove(fileName)\" before calling print
(\"remove(..)\" is silent, if the file does not exist).
The Modelica environment may close the file whenever appropriate.
This can be enforced by calling <b>Streams.close</b>(fileName).
After every call of \"print(..)\" a \"new line\" is printed automatically.
</p>
<h4>Example</h4>
<blockquote><pre>
  Streams.print(\"x = \" + String(x));
  Streams.print(\"y = \" + String(y));
  Streams.print(\"x = \" + String(y), \"mytestfile.txt\");
</pre></blockquote>
<h4>See also</h4>
<p>
<a href=\"modelica://Modelica.Utilities.Streams\">Streams</a>,
<a href=\"modelica://Modelica.Utilities.Streams.error\">Streams.error</a>,
<a href=\"modelica://ModelicaReference.Operators.'String()'\">ModelicaReference.Operators.'String()'</a>
</p>
</HTML>"));
      end print;
      annotation (
        Documentation(info="<HTML>
<h4>Library content</h4>
<p>
Package <b>Streams</b> contains functions to input and output strings
to a message window or on files. Note that a string is interpreted
and displayed as html text (e.g., with print(..) or error(..))
if it is enclosed with the Modelica html quotation, e.g.,
</p>
<center>
string = \"&lt;html&gt; first line &lt;br&gt; second line &lt;/html&gt;\".
</center>
<p>
It is a quality of implementation, whether (a) all tags of html are supported
or only a subset, (b) how html tags are interpreted if the output device
does not allow to display formatted text.
</p>
<p>
In the table below an example call to every function is given:
</p>
<table border=1 cellspacing=0 cellpadding=2>
  <tr><th><b><i>Function/type</i></b></th><th><b><i>Description</i></b></th></tr>
  <tr><td valign=\"top\"><a href=\"modelica://Modelica.Utilities.Streams.print\">print</a>(string)<br>
          <a href=\"modelica://Modelica.Utilities.Streams.print\">print</a>(string,fileName)</td>
      <td valign=\"top\"> Print string \"string\" or vector of strings to message window or on
           file \"fileName\".</td>
  </tr>
  <tr><td valign=\"top\">stringVector =
         <a href=\"modelica://Modelica.Utilities.Streams.readFile\">readFile</a>(fileName)</td>
      <td valign=\"top\"> Read complete text file and return it as a vector of strings.</td>
  </tr>
  <tr><td valign=\"top\">(string, endOfFile) =
         <a href=\"modelica://Modelica.Utilities.Streams.readLine\">readLine</a>(fileName, lineNumber)</td>
      <td valign=\"top\">Returns from the file the content of line lineNumber.</td>
  </tr>
  <tr><td valign=\"top\">lines =
         <a href=\"modelica://Modelica.Utilities.Streams.countLines\">countLines</a>(fileName)</td>
      <td valign=\"top\">Returns the number of lines in a file.</td>
  </tr>
  <tr><td valign=\"top\"><a href=\"modelica://Modelica.Utilities.Streams.error\">error</a>(string)</td>
      <td valign=\"top\"> Print error message \"string\" to message window
           and cancel all actions</td>
  </tr>
  <tr><td valign=\"top\"><a href=\"modelica://Modelica.Utilities.Streams.close\">close</a>(fileName)</td>
      <td valign=\"top\"> Close file if it is still open. Ignore call if
           file is already closed or does not exist. </td>
  </tr>
</table>
<p>
Use functions <b>scanXXX</b> from package
<a href=\"modelica://Modelica.Utilities.Strings\">Strings</a>
to parse a string.
</p>
<p>
If Real, Integer or Boolean values shall be printed
or used in an error message, they have to be first converted
to strings with the builtin operator
<a href=\"modelica://ModelicaReference.Operators.'String()'\">ModelicaReference.Operators.'String()'</a>(...).
Example:
</p>
<pre>
  <b>if</b> x &lt; 0 <b>or</b> x &gt; 1 <b>then</b>
     Streams.error(\"x (= \" + String(x) + \") has to be in the range 0 .. 1\");
  <b>end if</b>;
</pre>
</html>"));
    end Streams;

    package Strings "Operations on strings"
      extends Modelica.Icons.Package;

      function length "Returns length of string"
        extends Modelica.Icons.Function;
        input String string;
        output Integer result "Number of characters of string";
      external "C" result = ModelicaStrings_length(string) annotation(Library="ModelicaExternalC");
        annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
Strings.<b>length</b>(string);
</pre></blockquote>
<h4>Description</h4>
<p>
Returns the number of characters of \"string\".
</p>
</html>"));
      end length;

      function substring "Returns a substring defined by start and end index"

        extends Modelica.Icons.Function;
        input String string "String from which a substring is inquired";
        input Integer startIndex(min=1)
        "Character position of substring begin (index=1 is first character in string)";
        input Integer endIndex(min=1) "Character position of substring end";
        output String result
        "String containing substring string[startIndex:endIndex]";
      external "C" result =
                          ModelicaStrings_substring(string,startIndex,endIndex) annotation(Library="ModelicaExternalC");
        annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
string2 = Strings.<b>substring</b>(string, startIndex, endIndex);
</pre></blockquote>
<h4>Description</h4>
<p>
This function returns
the substring from position startIndex
up to and including position endIndex of \"string\" .
</p>
<p>
If index, startIndex, or endIndex are not correct, e.g.,
if endIndex &gt; length(string), an assert is triggered.
</p>
<h4>Example</h4>
<blockquote><pre>
  string1 := \"This is line 111\";
  string2 := Strings.substring(string1,9,12); // string2 = \"line\"
</pre></blockquote>
</html>"));
      end substring;

      function compare "Compare two strings lexicographically"
        extends Modelica.Icons.Function;
        input String string1;
        input String string2;
        input Boolean caseSensitive=true
        "= false, if case of letters is ignored";
        output Modelica.Utilities.Types.Compare result "Result of comparison";
      external "C" result = ModelicaStrings_compare(string1, string2, caseSensitive) annotation(Library="ModelicaExternalC");
        annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
result = Strings.<b>compare</b>(string1, string2);
result = Strings.<b>compare</b>(string1, string2, caseSensitive=true);
</pre></blockquote>
<h4>Description</h4>
<p>
Compares two strings. If the optional argument caseSensitive=false,
upper case letters are treated as if they would be lower case letters.
The result of the comparison is returned as:
</p>
<pre>
  result = Modelica.Utilities.Types.Compare.Less     // string1 &lt; string2
         = Modelica.Utilities.Types.Compare.Equal    // string1 = string2
         = Modelica.Utilities.Types.Compare.Greater  // string1 &gt; string2
</pre>
<p>
Comparison is with regards to lexicographical order,
e.g., \"a\" &lt; \"b\";
</p>
</html>"));
      end compare;

      function isEqual "Determine whether two strings are identical"
        extends Modelica.Icons.Function;
        input String string1;
        input String string2;
        input Boolean caseSensitive=true
        "= false, if lower and upper case are ignored for the comparison";
        output Boolean identical "True, if string1 is identical to string2";
      algorithm
        identical :=compare(string1, string2, caseSensitive) == Types.Compare.Equal;
        annotation (
      Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
Strings.<b>isEqual</b>(string1, string2);
Strings.<b>isEqual</b>(string1, string2, caseSensitive=true);
</pre></blockquote>
<h4>Description</h4>
<p>
Compare whether two strings are identical,
optionally ignoring case.
</p>
</html>"));
      end isEqual;

      function find "Find first occurrence of a string within another string"
        extends Modelica.Icons.Function;
        input String string "String that is analyzed";
        input String searchString "String that is searched for in string";
        input Integer startIndex(min=1)=1 "Start search at index startIndex";
        input Boolean caseSensitive=true
        "= false, if lower and upper case are ignored for the search";
         output Integer index
        "Index of the beginning of the first occurrence of 'searchString' within 'string', or zero if not present";
    protected
        Integer lengthSearchString = length(searchString);
        Integer len = lengthSearchString-1;
        Integer i = startIndex;
        Integer i_max = length(string) - lengthSearchString + 1;
      algorithm
        index := 0;
        while i <= i_max loop
           if isEqual(substring(string,i,i+len),
                      searchString, caseSensitive) then
              index := i;
              i := i_max + 1;
           else
              i := i+1;
           end if;
        end while;
        annotation (
      Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
index = Strings.<b>find</b>(string, searchString);
index = Strings.<b>find</b>(string, searchString, startIndex=1,
                     caseSensitive=true);
</pre></blockquote>
<h4>Description</h4>
<p>
Finds first occurrence of \"searchString\" within \"string\"
and return the corresponding index.
Start search at index \"startIndex\" (default = 1).
If the optional argument \"caseSensitive\" is false, lower
and upper case are ignored for the search.
If \"searchString\" is not found, a value of \"0\" is returned.
</p>
</html>"));
      end find;
      annotation (
        Documentation(info="<HTML>
<h4>Library content</h4>
<p>
Package <b>Strings</b> contains functions to manipulate strings.
</p>
<p>
In the table below an example
call to every function is given using the <b>default</b> options.
</p>
<table border=1 cellspacing=0 cellpadding=2>
  <tr><th><b><i>Function</i></b></th><th><b><i>Description</i></b></th></tr>
  <tr><td valign=\"top\">len = <a href=\"modelica://Modelica.Utilities.Strings.length\">length</a>(string)</td>
      <td valign=\"top\">Returns length of string</td></tr>
  <tr><td valign=\"top\">string2 = <a href=\"modelica://Modelica.Utilities.Strings.substring\">substring</a>(string1,startIndex,endIndex)
       </td>
      <td valign=\"top\">Returns a substring defined by start and end index</td></tr>
  <tr><td valign=\"top\">result = <a href=\"modelica://Modelica.Utilities.Strings.repeat\">repeat</a>(n)<br>
 result = <a href=\"modelica://Modelica.Utilities.Strings.repeat\">repeat</a>(n,string)</td>
      <td valign=\"top\">Repeat a blank or a string n times.</td></tr>
  <tr><td valign=\"top\">result = <a href=\"modelica://Modelica.Utilities.Strings.compare\">compare</a>(string1, string2)</td>
      <td valign=\"top\">Compares two substrings with regards to alphabetical order</td></tr>
  <tr><td valign=\"top\">identical =
<a href=\"modelica://Modelica.Utilities.Strings.isEqual\">isEqual</a>(string1,string2)</td>
      <td valign=\"top\">Determine whether two strings are identical</td></tr>
  <tr><td valign=\"top\">result = <a href=\"modelica://Modelica.Utilities.Strings.count\">count</a>(string,searchString)</td>
      <td valign=\"top\">Count the number of occurrences of a string</td></tr>
  <tr>
<td valign=\"top\">index = <a href=\"modelica://Modelica.Utilities.Strings.find\">find</a>(string,searchString)</td>
      <td valign=\"top\">Find first occurrence of a string in another string</td></tr>
<tr>
<td valign=\"top\">index = <a href=\"modelica://Modelica.Utilities.Strings.findLast\">findLast</a>(string,searchString)</td>
      <td valign=\"top\">Find last occurrence of a string in another string</td></tr>
  <tr><td valign=\"top\">string2 = <a href=\"modelica://Modelica.Utilities.Strings.replace\">replace</a>(string,searchString,replaceString)</td>
      <td valign=\"top\">Replace one or all occurrences of a string</td></tr>
  <tr><td valign=\"top\">stringVector2 = <a href=\"modelica://Modelica.Utilities.Strings.sort\">sort</a>(stringVector1)</td>
      <td valign=\"top\">Sort vector of strings in alphabetic order</td></tr>
  <tr><td valign=\"top\">(token, index) = <a href=\"modelica://Modelica.Utilities.Strings.scanToken\">scanToken</a>(string,startIndex)</td>
      <td valign=\"top\">Scan for a token (Real/Integer/Boolean/String/Identifier/Delimiter/NoToken)</td></tr>
  <tr><td valign=\"top\">(number, index) = <a href=\"modelica://Modelica.Utilities.Strings.scanReal\">scanReal</a>(string,startIndex)</td>
      <td valign=\"top\">Scan for a Real constant</td></tr>
  <tr><td valign=\"top\">(number, index) = <a href=\"modelica://Modelica.Utilities.Strings.scanInteger\">scanInteger</a>(string,startIndex)</td>
      <td valign=\"top\">Scan for an Integer constant</td></tr>
  <tr><td valign=\"top\">(boolean, index) = <a href=\"modelica://Modelica.Utilities.Strings.scanBoolean\">scanBoolean</a>(string,startIndex)</td>
      <td valign=\"top\">Scan for a Boolean constant</td></tr>
  <tr><td valign=\"top\">(string2, index) = <a href=\"modelica://Modelica.Utilities.Strings.scanString\">scanString</a>(string,startIndex)</td>
      <td valign=\"top\">Scan for a String constant</td></tr>
  <tr><td valign=\"top\">(identifier, index) = <a href=\"modelica://Modelica.Utilities.Strings.scanIdentifier\">scanIdentifier</a>(string,startIndex)</td>
      <td valign=\"top\">Scan for an identifier</td></tr>
  <tr><td valign=\"top\">(delimiter, index) = <a href=\"modelica://Modelica.Utilities.Strings.scanDelimiter\">scanDelimiter</a>(string,startIndex)</td>
      <td valign=\"top\">Scan for delimiters</td></tr>
  <tr><td valign=\"top\"><a href=\"modelica://Modelica.Utilities.Strings.scanNoToken\">scanNoToken</a>(string,startIndex)</td>
      <td valign=\"top\">Check that remaining part of string consists solely of <br>
          white space or line comments (\"// ...\\n\").</td></tr>
  <tr><td valign=\"top\"><a href=\"modelica://Modelica.Utilities.Strings.syntaxError\">syntaxError</a>(string,index,message)</td>
      <td valign=\"top\"> Print a \"syntax error message\" as well as a string and the <br>
           index at which scanning detected an error</td></tr>
</table>
<p>
The functions \"compare\", \"isEqual\", \"count\", \"find\", \"findLast\", \"replace\", \"sort\"
have the optional
input argument <b>caseSensitive</b> with default <b>true</b>.
If <b>false</b>, the operation is carried out without taking
into account whether a character is upper or lower case.
</p>
</HTML>"));
    end Strings;

    package Types "Type definitions used in package Modelica.Utilities"
      extends Modelica.Icons.TypesPackage;

      type Compare = enumeration(
        Less "String 1 is lexicographically less than string 2",
        Equal "String 1 is identical to string 2",
        Greater "String 1 is lexicographically greater than string 2")
      "Enumeration defining comparison of two strings";
      annotation (Documentation(info="<html>
<p>
This package contains type definitions used in Modelica.Utilities.
</p>

</html>"));
    end Types;
      annotation (
  Icon(coordinateSystem(extent={{-100.0,-100.0},{100.0,100.0}}), graphics={
      Polygon(
        origin={1.3835,-4.1418},
        rotation=45.0,
        fillColor={64,64,64},
        pattern=LinePattern.None,
        fillPattern=FillPattern.Solid,
        points={{-15.0,93.333},{-15.0,68.333},{0.0,58.333},{15.0,68.333},{15.0,93.333},{20.0,93.333},{25.0,83.333},{25.0,58.333},{10.0,43.333},{10.0,-41.667},{25.0,-56.667},{25.0,-76.667},{10.0,-91.667},{0.0,-91.667},{0.0,-81.667},{5.0,-81.667},{15.0,-71.667},{15.0,-61.667},{5.0,-51.667},{-5.0,-51.667},{-15.0,-61.667},{-15.0,-71.667},{-5.0,-81.667},{0.0,-81.667},{0.0,-91.667},{-10.0,-91.667},{-25.0,-76.667},{-25.0,-56.667},{-10.0,-41.667},{-10.0,43.333},{-25.0,58.333},{-25.0,83.333},{-20.0,93.333}}),
      Polygon(
        origin={10.1018,5.218},
        rotation=-45.0,
        fillColor={255,255,255},
        fillPattern=FillPattern.Solid,
        points={{-15.0,87.273},{15.0,87.273},{20.0,82.273},{20.0,27.273},{10.0,17.273},{10.0,7.273},{20.0,2.273},{20.0,-2.727},{5.0,-2.727},{5.0,-77.727},{10.0,-87.727},{5.0,-112.727},{-5.0,-112.727},{-10.0,-87.727},{-5.0,-77.727},{-5.0,-2.727},{-20.0,-2.727},{-20.0,2.273},{-10.0,7.273},{-10.0,17.273},{-20.0,27.273},{-20.0,82.273}})}),
  Documentation(info="<html>
<p>
This package contains Modelica <b>functions</b> that are
especially suited for <b>scripting</b>. The functions might
be used to work with strings, read data from file, write data
to file or copy, move and remove files.
</p>
<p>
For an introduction, have especially a look at:
</p>
<ul>
<li> <a href=\"modelica://Modelica.Utilities.UsersGuide\">Modelica.Utilities.User's Guide</a>
     discusses the most important aspects of this library.</li>
<li> <a href=\"modelica://Modelica.Utilities.Examples\">Modelica.Utilities.Examples</a>
     contains examples that demonstrate the usage of this library.</li>
</ul>
<p>
The following main sublibraries are available:
</p>
<ul>
<li> <a href=\"modelica://Modelica.Utilities.Files\">Files</a>
     provides functions to operate on files and directories, e.g.,
     to copy, move, remove files.</li>
<li> <a href=\"modelica://Modelica.Utilities.Streams\">Streams</a>
     provides functions to read from files and write to files.</li>
<li> <a href=\"modelica://Modelica.Utilities.Strings\">Strings</a>
     provides functions to operate on strings. E.g.
     substring, find, replace, sort, scanToken.</li>
<li> <a href=\"modelica://Modelica.Utilities.System\">System</a>
     provides functions to interact with the environment.
     E.g., get or set the working directory or environment
     variables and to send a command to the default shell.</li>
</ul>

<p>
Copyright &copy; 1998-2013, Modelica Association, DLR, and Dassault Syst&egrave;mes AB.
</p>

<p>
<i>This Modelica package is <u>free</u> software and the use is completely at <u>your own risk</u>; it can be redistributed and/or modified under the terms of the Modelica License 2. For license conditions (including the disclaimer of warranty) see <a href=\"modelica://Modelica.UsersGuide.ModelicaLicense2\">Modelica.UsersGuide.ModelicaLicense2</a> or visit <a href=\"https://www.modelica.org/licenses/ModelicaLicense2\"> https://www.modelica.org/licenses/ModelicaLicense2</a>.</i>
</p>

</html>"));
  end Utilities;

  package Constants
  "Library of mathematical constants and constants of nature (e.g., pi, eps, R, sigma)"
    import SI = Modelica.SIunits;
    import NonSI = Modelica.SIunits.Conversions.NonSIunits;
    extends Modelica.Icons.Package;

    final constant Real pi=2*Modelica.Math.asin(1.0);

    final constant Real eps=ModelicaServices.Machine.eps
    "Biggest number such that 1.0 + eps = 1.0";

    final constant Real small=ModelicaServices.Machine.small
    "Smallest number such that small and -small are representable on the machine";

    final constant SI.Acceleration g_n=9.80665
    "Standard acceleration of gravity on earth";

    final constant Real R(final unit="J/(mol.K)") = 8.314472
    "Molar gas constant";

    final constant NonSI.Temperature_degC T_zero=-273.15
    "Absolute zero temperature";
    annotation (
      Documentation(info="<html>
<p>
This package provides often needed constants from mathematics, machine
dependent constants and constants from nature. The latter constants
(name, value, description) are from the following source:
</p>

<dl>
<dt>Peter J. Mohr and Barry N. Taylor (1999):</dt>
<dd><b>CODATA Recommended Values of the Fundamental Physical Constants: 1998</b>.
    Journal of Physical and Chemical Reference Data, Vol. 28, No. 6, 1999 and
    Reviews of Modern Physics, Vol. 72, No. 2, 2000. See also <a href=
\"http://physics.nist.gov/cuu/Constants/\">http://physics.nist.gov/cuu/Constants/</a></dd>
</dl>

<p>CODATA is the Committee on Data for Science and Technology.</p>

<dl>
<dt><b>Main Author:</b></dt>
<dd><a href=\"http://www.robotic.dlr.de/Martin.Otter/\">Martin Otter</a><br>
    Deutsches Zentrum f&uuml;r Luft und Raumfahrt e. V. (DLR)<br>
    Oberpfaffenhofen<br>
    Postfach 11 16<br>
    D-82230 We&szlig;ling<br>
    email: <a href=\"mailto:Martin.Otter@dlr.de\">Martin.Otter@dlr.de</a></dd>
</dl>

<p>
Copyright &copy; 1998-2013, Modelica Association and DLR.
</p>
<p>
<i>This Modelica package is <u>free</u> software and the use is completely at <u>your own risk</u>; it can be redistributed and/or modified under the terms of the Modelica License 2. For license conditions (including the disclaimer of warranty) see <a href=\"modelica://Modelica.UsersGuide.ModelicaLicense2\">Modelica.UsersGuide.ModelicaLicense2</a> or visit <a href=\"https://www.modelica.org/licenses/ModelicaLicense2\"> https://www.modelica.org/licenses/ModelicaLicense2</a>.</i>
</p>
</html>",   revisions="<html>
<ul>
<li><i>Nov 8, 2004</i>
       by <a href=\"http://www.robotic.dlr.de/Christian.Schweiger/\">Christian Schweiger</a>:<br>
       Constants updated according to 2002 CODATA values.</li>
<li><i>Dec 9, 1999</i>
       by <a href=\"http://www.robotic.dlr.de/Martin.Otter/\">Martin Otter</a>:<br>
       Constants updated according to 1998 CODATA values. Using names, values
       and description text from this source. Included magnetic and
       electric constant.</li>
<li><i>Sep 18, 1999</i>
       by <a href=\"http://www.robotic.dlr.de/Martin.Otter/\">Martin Otter</a>:<br>
       Constants eps, inf, small introduced.</li>
<li><i>Nov 15, 1997</i>
       by <a href=\"http://www.robotic.dlr.de/Martin.Otter/\">Martin Otter</a>:<br>
       Realized.</li>
</ul>
</html>"),
      Icon(coordinateSystem(extent={{-100.0,-100.0},{100.0,100.0}}), graphics={
        Polygon(
          origin={-9.2597,25.6673},
          fillColor={102,102,102},
          pattern=LinePattern.None,
          fillPattern=FillPattern.Solid,
          points={{48.017,11.336},{48.017,11.336},{10.766,11.336},{-25.684,10.95},{-34.944,-15.111},{-34.944,-15.111},{-32.298,-15.244},{-32.298,-15.244},{-22.112,0.168},{11.292,0.234},{48.267,-0.097},{48.267,-0.097}},
          smooth=Smooth.Bezier),
        Polygon(
          origin={-19.9923,-8.3993},
          fillColor={102,102,102},
          pattern=LinePattern.None,
          fillPattern=FillPattern.Solid,
          points={{3.239,37.343},{3.305,37.343},{-0.399,2.683},{-16.936,-20.071},{-7.808,-28.604},{6.811,-22.519},{9.986,37.145},{9.986,37.145}},
          smooth=Smooth.Bezier),
        Polygon(
          origin={23.753,-11.5422},
          fillColor={102,102,102},
          pattern=LinePattern.None,
          fillPattern=FillPattern.Solid,
          points={{-10.873,41.478},{-10.873,41.478},{-14.048,-4.162},{-9.352,-24.8},{7.912,-24.469},{16.247,0.27},{16.247,0.27},{13.336,0.071},{13.336,0.071},{7.515,-9.983},{-3.134,-7.271},{-2.671,41.214},{-2.671,41.214}},
          smooth=Smooth.Bezier)}));
  end Constants;

  package Icons "Library of icons"
    extends Icons.Package;

    partial package ExamplesPackage
    "Icon for packages containing runnable examples"
      extends Modelica.Icons.Package;
      annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,
                -100},{100,100}}), graphics={
            Polygon(
              origin={8.0,14.0},
              lineColor={78,138,73},
              fillColor={78,138,73},
              pattern=LinePattern.None,
              fillPattern=FillPattern.Solid,
              points={{-58.0,46.0},{42.0,-14.0},{-58.0,-74.0},{-58.0,46.0}})}), Documentation(info="<html>
<p>This icon indicates a package that contains executable examples.</p>
</html>"));
    end ExamplesPackage;

    partial package Package "Icon for standard packages"

      annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,100}}), graphics={
            Rectangle(
              lineColor={200,200,200},
              fillColor={248,248,248},
              fillPattern=FillPattern.HorizontalCylinder,
              extent={{-100.0,-100.0},{100.0,100.0}},
              radius=25.0),
            Rectangle(
              lineColor={128,128,128},
              fillPattern=FillPattern.None,
              extent={{-100.0,-100.0},{100.0,100.0}},
              radius=25.0)}),   Documentation(info="<html>
<p>Standard package icon.</p>
</html>"));
    end Package;

    partial package BasesPackage "Icon for packages containing base classes"
      extends Modelica.Icons.Package;
      annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,
                -100},{100,100}}), graphics={
            Ellipse(
              extent={{-30.0,-30.0},{30.0,30.0}},
              lineColor={128,128,128},
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid)}),
                                Documentation(info="<html>
<p>This icon shall be used for a package/library that contains base models and classes, respectively.</p>
</html>"));
    end BasesPackage;

    partial package InterfacesPackage "Icon for packages containing interfaces"
      extends Modelica.Icons.Package;
      annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,
                -100},{100,100}}), graphics={
            Polygon(origin={20.0,0.0},
              lineColor={64,64,64},
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid,
              points={{-10.0,70.0},{10.0,70.0},{40.0,20.0},{80.0,20.0},{80.0,-20.0},{40.0,-20.0},{10.0,-70.0},{-10.0,-70.0}}),
            Polygon(fillColor={102,102,102},
              pattern=LinePattern.None,
              fillPattern=FillPattern.Solid,
              points={{-100.0,20.0},{-60.0,20.0},{-30.0,70.0},{-10.0,70.0},{-10.0,-70.0},{-30.0,-70.0},{-60.0,-20.0},{-100.0,-20.0}})}),
                                Documentation(info="<html>
<p>This icon indicates packages containing interfaces.</p>
</html>"));
    end InterfacesPackage;

    partial package TypesPackage
    "Icon for packages containing type definitions"
      extends Modelica.Icons.Package;
      annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,
                -100},{100,100}}), graphics={Polygon(
              origin={-12.167,-23},
              fillColor={128,128,128},
              pattern=LinePattern.None,
              fillPattern=FillPattern.Solid,
              points={{12.167,65},{14.167,93},{36.167,89},{24.167,20},{4.167,-30},
                  {14.167,-30},{24.167,-30},{24.167,-40},{-5.833,-50},{-15.833,
                  -30},{4.167,20},{12.167,65}},
              smooth=Smooth.Bezier,
              lineColor={0,0,0}), Polygon(
              origin={2.7403,1.6673},
              fillColor={128,128,128},
              pattern=LinePattern.None,
              fillPattern=FillPattern.Solid,
              points={{49.2597,22.3327},{31.2597,24.3327},{7.2597,18.3327},{-26.7403,
                10.3327},{-46.7403,14.3327},{-48.7403,6.3327},{-32.7403,0.3327},{-6.7403,
                4.3327},{33.2597,14.3327},{49.2597,14.3327},{49.2597,22.3327}},
              smooth=Smooth.Bezier)}));
    end TypesPackage;

    partial package IconsPackage "Icon for packages containing icons"
      extends Modelica.Icons.Package;
      annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,
                -100},{100,100}}), graphics={Polygon(
              origin={-8.167,-17},
              fillColor={128,128,128},
              pattern=LinePattern.None,
              fillPattern=FillPattern.Solid,
              points={{-15.833,20.0},{-15.833,30.0},{14.167,40.0},{24.167,20.0},{
                  4.167,-30.0},{14.167,-30.0},{24.167,-30.0},{24.167,-40.0},{-5.833,
                  -50.0},{-15.833,-30.0},{4.167,20.0},{-5.833,20.0}},
              smooth=Smooth.Bezier,
              lineColor={0,0,0}), Ellipse(
              origin={-0.5,56.5},
              fillColor={128,128,128},
              pattern=LinePattern.None,
              fillPattern=FillPattern.Solid,
              extent={{-12.5,-12.5},{12.5,12.5}},
              lineColor={0,0,0})}));
    end IconsPackage;

    partial package MaterialPropertiesPackage
    "Icon for package containing property classes"
      extends Modelica.Icons.Package;
      annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,
                -100},{100,100}}), graphics={
            Ellipse(
              lineColor={102,102,102},
              fillColor={204,204,204},
              pattern=LinePattern.None,
              fillPattern=FillPattern.Sphere,
              extent={{-60.0,-60.0},{60.0,60.0}})}),
                                Documentation(info="<html>
<p>This icon indicates a package that contains properties</p>
</html>"));
    end MaterialPropertiesPackage;

    partial class MaterialProperty "Icon for property classes"

      annotation (Icon(coordinateSystem(preserveAspectRatio=true,  extent={{-100,-100},{100,100}}), graphics={
            Ellipse(lineColor={102,102,102},
              fillColor={204,204,204},
              pattern=LinePattern.None,
              fillPattern=FillPattern.Sphere,
              extent={{-100.0,-100.0},{100.0,100.0}})}),
                                Documentation(info="<html>
<p>This icon indicates a property class.</p>
</html>"));
    end MaterialProperty;

    partial function Function "Icon for functions"

      annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,100}}), graphics={
            Text(
              lineColor={0,0,255},
              extent={{-150,105},{150,145}},
              textString="%name"),
            Ellipse(
              lineColor = {108,88,49},
              fillColor = {255,215,136},
              fillPattern = FillPattern.Solid,
              extent = {{-100,-100},{100,100}}),
            Text(
              lineColor={108,88,49},
              extent={{-90.0,-90.0},{90.0,90.0}},
              textString="f")}),
    Documentation(info="<html>
<p>This icon indicates Modelica functions.</p>
</html>"));
    end Function;

    partial record Record "Icon for records"

      annotation (Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},{100,100}}), graphics={
            Text(
              lineColor={0,0,255},
              extent={{-150,60},{150,100}},
              textString="%name"),
            Rectangle(
              origin={0.0,-25.0},
              lineColor={64,64,64},
              fillColor={255,215,136},
              fillPattern=FillPattern.Solid,
              extent={{-100.0,-75.0},{100.0,75.0}},
              radius=25.0),
            Line(
              points={{-100.0,0.0},{100.0,0.0}},
              color={64,64,64}),
            Line(
              origin={0.0,-50.0},
              points={{-100.0,0.0},{100.0,0.0}},
              color={64,64,64}),
            Line(
              origin={0.0,-25.0},
              points={{0.0,75.0},{0.0,-75.0}},
              color={64,64,64})}),                        Documentation(info="<html>
<p>
This icon is indicates a record.
</p>
</html>"));
    end Record;
    annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,
                -100},{100,100}}), graphics={Polygon(
              origin={-8.167,-17},
              fillColor={128,128,128},
              pattern=LinePattern.None,
              fillPattern=FillPattern.Solid,
              points={{-15.833,20.0},{-15.833,30.0},{14.167,40.0},{24.167,20.0},{
                  4.167,-30.0},{14.167,-30.0},{24.167,-30.0},{24.167,-40.0},{-5.833,
                  -50.0},{-15.833,-30.0},{4.167,20.0},{-5.833,20.0}},
              smooth=Smooth.Bezier,
              lineColor={0,0,0}), Ellipse(
              origin={-0.5,56.5},
              fillColor={128,128,128},
              pattern=LinePattern.None,
              fillPattern=FillPattern.Solid,
              extent={{-12.5,-12.5},{12.5,12.5}},
              lineColor={0,0,0})}), Documentation(info="<html>
<p>This package contains definitions for the graphical layout of components which may be used in different libraries. The icons can be utilized by inheriting them in the desired class using &quot;extends&quot; or by directly copying the &quot;icon&quot; layer. </p>

<h4>Main Authors:</h4>

<dl>
<dt><a href=\"http://www.robotic.dlr.de/Martin.Otter/\">Martin Otter</a></dt>
    <dd>Deutsches Zentrum fuer Luft und Raumfahrt e.V. (DLR)</dd>
    <dd>Oberpfaffenhofen</dd>
    <dd>Postfach 1116</dd>
    <dd>D-82230 Wessling</dd>
    <dd>email: <a href=\"mailto:Martin.Otter@dlr.de\">Martin.Otter@dlr.de</a></dd>
<dt>Christian Kral</dt>
    <dd><a href=\"http://www.ait.ac.at/\">Austrian Institute of Technology, AIT</a></dd>
    <dd>Mobility Department</dd><dd>Giefinggasse 2</dd>
    <dd>1210 Vienna, Austria</dd>
    <dd>email: <a href=\"mailto:dr.christian.kral@gmail.com\">dr.christian.kral@gmail.com</a></dd>
<dt>Johan Andreasson</dt>
    <dd><a href=\"http://www.modelon.se/\">Modelon AB</a></dd>
    <dd>Ideon Science Park</dd>
    <dd>22370 Lund, Sweden</dd>
    <dd>email: <a href=\"mailto:johan.andreasson@modelon.se\">johan.andreasson@modelon.se</a></dd>
</dl>

<p>Copyright &copy; 1998-2013, Modelica Association, DLR, AIT, and Modelon AB. </p>
<p><i>This Modelica package is <b>free</b> software; it can be redistributed and/or modified under the terms of the <b>Modelica license</b>, see the license conditions and the accompanying <b>disclaimer</b> in <a href=\"modelica://Modelica.UsersGuide.ModelicaLicense2\">Modelica.UsersGuide.ModelicaLicense2</a>.</i> </p>
</html>"));
  end Icons;

  package SIunits
  "Library of type and unit definitions based on SI units according to ISO 31-1992"
    extends Modelica.Icons.Package;

    package Icons "Icons for SIunits"
      extends Modelica.Icons.IconsPackage;

      partial function Conversion "Base icon for conversion functions"

        annotation (Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,
                  -100},{100,100}}), graphics={
              Rectangle(
                extent={{-100,100},{100,-100}},
                lineColor={191,0,0},
                fillColor={255,255,255},
                fillPattern=FillPattern.Solid),
              Line(points={{-90,0},{30,0}}, color={191,0,0}),
              Polygon(
                points={{90,0},{30,20},{30,-20},{90,0}},
                lineColor={191,0,0},
                fillColor={191,0,0},
                fillPattern=FillPattern.Solid),
              Text(
                extent={{-115,155},{115,105}},
                textString="%name",
                lineColor={0,0,255})}));
      end Conversion;
    end Icons;

    package Conversions
    "Conversion functions to/from non SI units and type definitions of non SI units"
      extends Modelica.Icons.Package;

      package NonSIunits "Type definitions of non SI units"
        extends Modelica.Icons.Package;

        type Temperature_degC = Real (final quantity="ThermodynamicTemperature",
              final unit="degC")
        "Absolute temperature in degree Celsius (for relative temperature use SIunits.TemperatureDifference)"
                                                                                                            annotation(absoluteValue=true);

        type Pressure_bar = Real (final quantity="Pressure", final unit="bar")
        "Absolute pressure in bar";
        annotation (Documentation(info="<HTML>
<p>
This package provides predefined types, such as <b>Angle_deg</b> (angle in
degree), <b>AngularVelocity_rpm</b> (angular velocity in revolutions per
minute) or <b>Temperature_degF</b> (temperature in degree Fahrenheit),
which are in common use but are not part of the international standard on
units according to ISO 31-1992 \"General principles concerning quantities,
units and symbols\" and ISO 1000-1992 \"SI units and recommendations for
the use of their multiples and of certain other units\".</p>
<p>If possible, the types in this package should not be used. Use instead
types of package Modelica.SIunits. For more information on units, see also
the book of Francois Cardarelli <b>Scientific Unit Conversion - A
Practical Guide to Metrication</b> (Springer 1997).</p>
<p>Some units, such as <b>Temperature_degC/Temp_C</b> are both defined in
Modelica.SIunits and in Modelica.Conversions.NonSIunits. The reason is that these
definitions have been placed erroneously in Modelica.SIunits although they
are not SIunits. For backward compatibility, these type definitions are
still kept in Modelica.SIunits.</p>
</html>"),   Icon(coordinateSystem(extent={{-100,-100},{100,100}}), graphics={
        Text(
          origin={15.0,51.8518},
          extent={{-105.0,-86.8518},{75.0,-16.8518}},
          lineColor={0,0,0},
          textString="[km/h]")}));
      end NonSIunits;

      function to_degC "Convert from Kelvin to degCelsius"
        extends Modelica.SIunits.Icons.Conversion;
        input Temperature Kelvin "Kelvin value";
        output NonSIunits.Temperature_degC Celsius "Celsius value";
      algorithm
        Celsius := Kelvin + Modelica.Constants.T_zero;
        annotation (Inline=true,Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,
                  -100},{100,100}}), graphics={Text(
                extent={{-20,100},{-100,20}},
                lineColor={0,0,0},
                textString="K"), Text(
                extent={{100,-20},{20,-100}},
                lineColor={0,0,0},
                textString="degC")}));
      end to_degC;

      function from_degC "Convert from degCelsius to Kelvin"
        extends Modelica.SIunits.Icons.Conversion;
        input NonSIunits.Temperature_degC Celsius "Celsius value";
        output Temperature Kelvin "Kelvin value";
      algorithm
        Kelvin := Celsius - Modelica.Constants.T_zero;
        annotation (Inline=true,Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,
                  -100},{100,100}}), graphics={Text(
                extent={{-20,100},{-100,20}},
                lineColor={0,0,0},
                textString="degC"),  Text(
                extent={{100,-20},{20,-100}},
                lineColor={0,0,0},
                textString="K")}));
      end from_degC;

      function to_bar "Convert from Pascal to bar"
        extends Modelica.SIunits.Icons.Conversion;
        input Pressure Pa "Pascal value";
        output NonSIunits.Pressure_bar bar "bar value";
      algorithm
        bar := Pa/1e5;
        annotation (Inline=true,Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,
                  -100},{100,100}}), graphics={Text(
                extent={{-12,100},{-100,56}},
                lineColor={0,0,0},
                textString="Pa"),     Text(
                extent={{98,-52},{-4,-100}},
                lineColor={0,0,0},
                textString="bar")}));
      end to_bar;
      annotation (                              Documentation(info="<HTML>
<p>This package provides conversion functions from the non SI Units
defined in package Modelica.SIunits.Conversions.NonSIunits to the
corresponding SI Units defined in package Modelica.SIunits and vice
versa. It is recommended to use these functions in the following
way (note, that all functions have one Real input and one Real output
argument):</p>
<pre>
  <b>import</b> SI = Modelica.SIunits;
  <b>import</b> Modelica.SIunits.Conversions.*;
     ...
  <b>parameter</b> SI.Temperature     T   = from_degC(25);   // convert 25 degree Celsius to Kelvin
  <b>parameter</b> SI.Angle           phi = from_deg(180);   // convert 180 degree to radian
  <b>parameter</b> SI.AngularVelocity w   = from_rpm(3600);  // convert 3600 revolutions per minutes
                                                      // to radian per seconds
</pre>

</html>"));
    end Conversions;

    type Angle = Real (
        final quantity="Angle",
        final unit="rad",
        displayUnit="deg");

    type Length = Real (final quantity="Length", final unit="m");

    type Area = Real (final quantity="Area", final unit="m2");

    type Velocity = Real (final quantity="Velocity", final unit="m/s");

    type Acceleration = Real (final quantity="Acceleration", final unit="m/s2");

    type Mass = Real (
        quantity="Mass",
        final unit="kg",
        min=0);

    type Density = Real (
        final quantity="Density",
        final unit="kg/m3",
        displayUnit="g/cm3",
        min=0.0);

    type Pressure = Real (
        final quantity="Pressure",
        final unit="Pa",
        displayUnit="bar");

    type AbsolutePressure = Pressure (min=0.0, nominal = 1e5);

    type DynamicViscosity = Real (
        final quantity="DynamicViscosity",
        final unit="Pa.s",
        min=0);

    type SurfaceTension = Real (final quantity="SurfaceTension", final unit="N/m");

    type Power = Real (final quantity="Power", final unit="W");

    type MassFlowRate = Real (quantity="MassFlowRate", final unit="kg/s");

    type ThermodynamicTemperature = Real (
        final quantity="ThermodynamicTemperature",
        final unit="K",
        min = 0.0,
        start = 288.15,
        nominal = 300,
        displayUnit="degC")
    "Absolute temperature (use type TemperatureDifference for relative temperatures)"
                                                                                                        annotation(absoluteValue=true);

    type Temperature = ThermodynamicTemperature;

    type TemperatureDifference = Real (
        final quantity="ThermodynamicTemperature",
        final unit="K") annotation(absoluteValue=false);

    type Temp_C = SIunits.Conversions.NonSIunits.Temperature_degC;

    type CubicExpansionCoefficient = Real (final quantity=
            "CubicExpansionCoefficient", final unit="1/K");

    type Compressibility = Real (final quantity="Compressibility", final unit=
            "1/Pa");

    type IsothermalCompressibility = Compressibility;

    type HeatFlowRate = Real (final quantity="Power", final unit="W");

    type HeatFlux = Real (final quantity="HeatFlux", final unit="W/m2");

    type ThermalConductivity = Real (final quantity="ThermalConductivity", final unit=
               "W/(m.K)");

    type CoefficientOfHeatTransfer = Real (final quantity=
            "CoefficientOfHeatTransfer", final unit="W/(m2.K)");

    type ThermalConductance = Real (final quantity="ThermalConductance", final unit=
               "W/K");

    type HeatCapacity = Real (final quantity="HeatCapacity", final unit="J/K");

    type SpecificHeatCapacity = Real (final quantity="SpecificHeatCapacity",
          final unit="J/(kg.K)");

    type RatioOfSpecificHeatCapacities = Real (final quantity=
            "RatioOfSpecificHeatCapacities", final unit="1");

    type SpecificEntropy = Real (final quantity="SpecificEntropy",
                                 final unit="J/(kg.K)");

    type SpecificEnergy = Real (final quantity="SpecificEnergy",
                                final unit="J/kg");

    type SpecificEnthalpy = SpecificEnergy;

    type DerDensityByEnthalpy = Real (final unit="kg.s2/m5");

    type DerDensityByPressure = Real (final unit="s2/m2");

    type DerDensityByTemperature = Real (final unit="kg/(m3.K)");

    type DerEnthalpyByPressure = Real (final unit="J.m.s2/kg2");

    type MolarMass = Real (final quantity="MolarMass", final unit="kg/mol",min=0);

    type MolarVolume = Real (final quantity="MolarVolume", final unit="m3/mol", min=0);

    type MassFraction = Real (final quantity="MassFraction", final unit="1",
                              min=0, max=1);

    type ReynoldsNumber = Real (final quantity="ReynoldsNumber", final unit="1");

    type NusseltNumber = Real (final quantity="NusseltNumber", final unit="1");

    type PrandtlNumber = Real (final quantity="PrandtlNumber", final unit="1");
    annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,
              -100},{100,100}}), graphics={
          Line(
            points={{-66,78},{-66,-40}},
            color={64,64,64},
            smooth=Smooth.None),
          Ellipse(
            extent={{12,36},{68,-38}},
            lineColor={64,64,64},
            fillColor={175,175,175},
            fillPattern=FillPattern.Solid),
          Rectangle(
            extent={{-74,78},{-66,-40}},
            lineColor={64,64,64},
            fillColor={175,175,175},
            fillPattern=FillPattern.Solid),
          Polygon(
            points={{-66,-4},{-66,6},{-16,56},{-16,46},{-66,-4}},
            lineColor={64,64,64},
            smooth=Smooth.None,
            fillColor={175,175,175},
            fillPattern=FillPattern.Solid),
          Polygon(
            points={{-46,16},{-40,22},{-2,-40},{-10,-40},{-46,16}},
            lineColor={64,64,64},
            smooth=Smooth.None,
            fillColor={175,175,175},
            fillPattern=FillPattern.Solid),
          Ellipse(
            extent={{22,26},{58,-28}},
            lineColor={64,64,64},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid),
          Polygon(
            points={{68,2},{68,-46},{64,-60},{58,-68},{48,-72},{18,-72},{18,-64},
                {46,-64},{54,-60},{58,-54},{60,-46},{60,-26},{64,-20},{68,-6},{68,
                2}},
            lineColor={64,64,64},
            smooth=Smooth.Bezier,
            fillColor={175,175,175},
            fillPattern=FillPattern.Solid)}), Documentation(info="<html>
<p>This package provides predefined types, such as <i>Mass</i>,
<i>Angle</i>, <i>Time</i>, based on the international standard
on units, e.g.,
</p>

<pre>   <b>type</b> Angle = Real(<b>final</b> quantity = \"Angle\",
                     <b>final</b> unit     = \"rad\",
                     displayUnit    = \"deg\");
</pre>

<p>
as well as conversion functions from non SI-units to SI-units
and vice versa in subpackage
<a href=\"modelica://Modelica.SIunits.Conversions\">Conversions</a>.
</p>

<p>
For an introduction how units are used in the Modelica standard library
with package SIunits, have a look at:
<a href=\"modelica://Modelica.SIunits.UsersGuide.HowToUseSIunits\">How to use SIunits</a>.
</p>

<p>
Copyright &copy; 1998-2013, Modelica Association and DLR.
</p>
<p>
<i>This Modelica package is <u>free</u> software and the use is completely at <u>your own risk</u>; it can be redistributed and/or modified under the terms of the Modelica License 2. For license conditions (including the disclaimer of warranty) see <a href=\"modelica://Modelica.UsersGuide.ModelicaLicense2\">Modelica.UsersGuide.ModelicaLicense2</a> or visit <a href=\"https://www.modelica.org/licenses/ModelicaLicense2\"> https://www.modelica.org/licenses/ModelicaLicense2</a>.</i>
</p>
</html>",   revisions="<html>
<ul>
<li><i>May 25, 2011</i> by Stefan Wischhusen:<br/>Added molar units for energy and enthalpy.</li>
<li><i>Jan. 27, 2010</i> by Christian Kral:<br/>Added complex units.</li>
<li><i>Dec. 14, 2005</i> by <a href=\"http://www.robotic.dlr.de/Martin.Otter/\">Martin Otter</a>:<br/>Add User&#39;;s Guide and removed &quot;min&quot; values for Resistance and Conductance.</li>
<li><i>October 21, 2002</i> by <a href=\"http://www.robotic.dlr.de/Martin.Otter/\">Martin Otter</a> and <a href=\"http://www.robotic.dlr.de/Christian.Schweiger/\">Christian Schweiger</a>:<br/>Added new package <b>Conversions</b>. Corrected typo <i>Wavelenght</i>.</li>
<li><i>June 6, 2000</i> by <a href=\"http://www.robotic.dlr.de/Martin.Otter/\">Martin Otter</a>:<br/>Introduced the following new types<br/>type Temperature = ThermodynamicTemperature;<br/>types DerDensityByEnthalpy, DerDensityByPressure, DerDensityByTemperature, DerEnthalpyByPressure, DerEnergyByDensity, DerEnergyByPressure<br/>Attribute &quot;final&quot; removed from min and max values in order that these values can still be changed to narrow the allowed range of values.<br/>Quantity=&quot;Stress&quot; removed from type &quot;Stress&quot;, in order that a type &quot;Stress&quot; can be connected to a type &quot;Pressure&quot;.</li>
<li><i>Oct. 27, 1999</i> by <a href=\"http://www.robotic.dlr.de/Martin.Otter/\">Martin Otter</a>:<br/>New types due to electrical library: Transconductance, InversePotential, Damping.</li>
<li><i>Sept. 18, 1999</i> by <a href=\"http://www.robotic.dlr.de/Martin.Otter/\">Martin Otter</a>:<br/>Renamed from SIunit to SIunits. Subpackages expanded, i.e., the SIunits package, does no longer contain subpackages.</li>
<li><i>Aug 12, 1999</i> by <a href=\"http://www.robotic.dlr.de/Martin.Otter/\">Martin Otter</a>:<br/>Type &quot;Pressure&quot; renamed to &quot;AbsolutePressure&quot; and introduced a new type &quot;Pressure&quot; which does not contain a minimum of zero in order to allow convenient handling of relative pressure. Redefined BulkModulus as an alias to AbsolutePressure instead of Stress, since needed in hydraulics.</li>
<li><i>June 29, 1999</i> by <a href=\"http://www.robotic.dlr.de/Martin.Otter/\">Martin Otter</a>:<br/>Bug-fix: Double definition of &quot;Compressibility&quot; removed and appropriate &quot;extends Heat&quot; clause introduced in package SolidStatePhysics to incorporate ThermodynamicTemperature.</li>
<li><i>April 8, 1998</i> by <a href=\"http://www.robotic.dlr.de/Martin.Otter/\">Martin Otter</a> and Astrid Jaschinski:<br/>Complete ISO 31 chapters realized.</li>
<li><i>Nov. 15, 1997</i> by <a href=\"http://www.robotic.dlr.de/Martin.Otter/\">Martin Otter</a> and <a href=\"http://www.control.lth.se/~hubertus/\">Hubertus Tummescheit</a>:<br/>Some chapters realized.</li>
</ul>
</html>"));
  end SIunits;
annotation (
preferredView="info",
version="3.2.1",
versionBuild=2,
versionDate="2013-08-14",
dateModified = "2013-08-14 08:44:41Z",
revisionId="$Id:: package.mo 6947 2013-08-23 07:41:37Z #$",
uses(Complex(version="3.2.1"), ModelicaServices(version="3.2.1")),
conversion(
 noneFromVersion="3.2",
 noneFromVersion="3.1",
 noneFromVersion="3.0.1",
 noneFromVersion="3.0",
 from(version="2.1", script="modelica://Modelica/Resources/Scripts/Dymola/ConvertModelica_from_2.2.2_to_3.0.mos"),
 from(version="2.2", script="modelica://Modelica/Resources/Scripts/Dymola/ConvertModelica_from_2.2.2_to_3.0.mos"),
 from(version="2.2.1", script="modelica://Modelica/Resources/Scripts/Dymola/ConvertModelica_from_2.2.2_to_3.0.mos"),
 from(version="2.2.2", script="modelica://Modelica/Resources/Scripts/Dymola/ConvertModelica_from_2.2.2_to_3.0.mos")),
Icon(coordinateSystem(extent={{-100.0,-100.0},{100.0,100.0}}), graphics={
  Polygon(
    origin={-6.9888,20.048},
    fillColor={0,0,0},
    pattern=LinePattern.None,
    fillPattern=FillPattern.Solid,
    points={{-93.0112,10.3188},{-93.0112,10.3188},{-73.011,24.6},{-63.011,31.221},{-51.219,36.777},{-39.842,38.629},{-31.376,36.248},{-25.819,29.369},{-24.232,22.49},{-23.703,17.463},{-15.501,25.135},{-6.24,32.015},{3.02,36.777},{15.191,39.423},{27.097,37.306},{32.653,29.633},{35.035,20.108},{43.501,28.046},{54.085,35.19},{65.991,39.952},{77.897,39.688},{87.422,33.338},{91.126,21.696},{90.068,9.525},{86.099,-1.058},{79.749,-10.054},{71.283,-21.431},{62.816,-33.337},{60.964,-32.808},{70.489,-16.14},{77.368,-2.381},{81.072,10.054},{79.749,19.05},{72.605,24.342},{61.758,23.019},{49.587,14.817},{39.003,4.763},{29.214,-6.085},{21.012,-16.669},{13.339,-26.458},{5.401,-36.777},{-1.213,-46.037},{-6.24,-53.446},{-8.092,-52.387},{-0.684,-40.746},{5.401,-30.692},{12.81,-17.198},{19.424,-3.969},{23.658,7.938},{22.335,18.785},{16.514,23.283},{8.047,23.019},{-1.478,19.05},{-11.267,11.113},{-19.734,2.381},{-29.259,-8.202},{-38.519,-19.579},{-48.044,-31.221},{-56.511,-43.392},{-64.449,-55.298},{-72.386,-66.939},{-77.678,-74.612},{-79.53,-74.083},{-71.857,-61.383},{-62.861,-46.037},{-52.278,-28.046},{-44.869,-15.346},{-38.784,-2.117},{-35.344,8.731},{-36.403,19.844},{-42.488,23.813},{-52.013,22.49},{-60.744,16.933},{-68.947,10.054},{-76.884,2.646},{-93.0112,-12.1707},{-93.0112,-12.1707}},
    smooth=Smooth.Bezier),
  Ellipse(
    origin={40.8208,-37.7602},
    fillColor={161,0,4},
    pattern=LinePattern.None,
    fillPattern=FillPattern.Solid,
    extent={{-17.8562,-17.8563},{17.8563,17.8562}})}),
Documentation(info="<HTML>
<p>
Package <b>Modelica&reg;</b> is a <b>standardized</b> and <b>free</b> package
that is developed together with the Modelica&reg; language from the
Modelica Association, see
<a href=\"https://www.Modelica.org\">https://www.Modelica.org</a>.
It is also called <b>Modelica Standard Library</b>.
It provides model components in many domains that are based on
standardized interface definitions. Some typical examples are shown
in the next figure:
</p>

<p>
<img src=\"modelica://Modelica/Resources/Images/UsersGuide/ModelicaLibraries.png\">
</p>

<p>
For an introduction, have especially a look at:
</p>
<ul>
<li> <a href=\"modelica://Modelica.UsersGuide.Overview\">Overview</a>
  provides an overview of the Modelica Standard Library
  inside the <a href=\"modelica://Modelica.UsersGuide\">User's Guide</a>.</li>
<li><a href=\"modelica://Modelica.UsersGuide.ReleaseNotes\">Release Notes</a>
 summarizes the changes of new versions of this package.</li>
<li> <a href=\"modelica://Modelica.UsersGuide.Contact\">Contact</a>
  lists the contributors of the Modelica Standard Library.</li>
<li> The <b>Examples</b> packages in the various libraries, demonstrate
  how to use the components of the corresponding sublibrary.</li>
</ul>

<p>
This version of the Modelica Standard Library consists of
</p>
<ul>
<li><b>1360</b> models and blocks, and</li>
<li><b>1280</b> functions</li>
</ul>
<p>
that are directly usable (= number of public, non-partial classes). It is fully compliant
to <a href=\"https://www.modelica.org/documents/ModelicaSpec32Revision2.pdf\">Modelica Specification Version 3.2 Revision 2</a>
and it has been tested with Modelica tools from different vendors.
</p>

<p>
<b>Licensed by the Modelica Association under the Modelica License 2</b><br>
Copyright &copy; 1998-2013, ABB, AIT, T.&nbsp;B&ouml;drich, DLR, Dassault Syst&egrave;mes AB, Fraunhofer, A.Haumer, ITI, Modelon,
TU Hamburg-Harburg, Politecnico di Milano, XRG Simulation.
</p>

<p>
<i>This Modelica package is <u>free</u> software and the use is completely at <u>your own risk</u>; it can be redistributed and/or modified under the terms of the Modelica License 2. For license conditions (including the disclaimer of warranty) see <a href=\"modelica://Modelica.UsersGuide.ModelicaLicense2\">Modelica.UsersGuide.ModelicaLicense2</a> or visit <a href=\"https://www.modelica.org/licenses/ModelicaLicense2\"> https://www.modelica.org/licenses/ModelicaLicense2</a>.</i>
</p>

<p>
<b>Modelica&reg;</b> is a registered trademark of the Modelica Association.
</p>
</html>"));
end Modelica;

package VIP "I am a package for the Virtual Prototyping Environment"
  import Modelica.Math.*;
  import Modelica.SIunits.*;
  import Modelica.Constants.*;

  package Design "Package for the component design"
    package Tests "I test the components here"
    extends Modelica.Icons.ExamplesPackage;

      package Shell_and_tube "Tests for shell and tube heat exchangers"
        package Single_phase "Tests for one phase to one phase heat exchangers"
          model Coulson_Kern
            "Verification with the results given by Coulson et al. using the Kern method"
            import VIP;

            replaceable package Medium_s = VIP.Media.RefProp.Methanol
              "Medium model";
            replaceable package Medium_t = VIP.Media.CoolProp.Water_TTSE
              "Medium model";
            parameter Modelica.SIunits.SpecificEnthalpy h_s_in= Medium_s.specificEnthalpy_pT(3.2e5,368.15)
              "Inlet specific enthalpy shell side";
            parameter Modelica.SIunits.SpecificEnthalpy h_s_out= Medium_s.specificEnthalpy_pT(3.2e5,313.15)
              "Outlet specific enthalpy shell side";
            parameter Modelica.SIunits.SpecificEnthalpy h_t_in = Medium_t.specificEnthalpy_pT(1e5,298.15)
              "Inlet specific enthalpy tube side";
            parameter Modelica.SIunits.SpecificEnthalpy h_t_out = Medium_t.specificEnthalpy_pT(1e5,313.15)
              "Outlet specific enthalpy tube side";

            Components.HEX.shell_and_tube             shell_and_tube(
              redeclare package Medium_s = Medium_s,
              redeclare package Medium_t = Medium_t,
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
              redeclare VIP.Design.Materials.S_AISI_1040 Material_t,
              redeclare VIP.Design.Materials.S_AISI_1040 Material_s,
              redeclare function bundle_clearance =
                  VIP.Design.Miscellanea.Shell_clearance.SRFH,
              pin_s=1,
              pin_t=-1,
              thick_t=2e-3,
              thick_s=10e-3,
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

            replaceable package Medium_s = VIP.Media.RefProp.Methanol
              "Medium model";
            replaceable package Medium_t = VIP.Media.CoolProp.Water
              "Medium model";
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
               redeclare Heat_transfer.Shell.single_phase_Bell_Delaware hT_shell(pitch_f = st.pitch_f, Dhyd_o = st.bundle.Dhyd_o,
               N_tubes = st.N_tubes, d_s = st.d_s, d_b = st.bundle.d_b, Dhyd = st.d_s,
               layout = st.layout, l_b = st.l_b, N_ss = 5,
               N_passes = st.N_passes, ff =  1.0, b_cut = 0.25),
               redeclare Pressure_drops.Shell.single_phase_Bell_Delaware dp_shell(
               pitch_f=st.pitch_f,
               Dhyd_o=st.bundle.Dhyd_o,
               N_tubes=st.N_tubes,
               d_s=st.d_s,
               d_b=st.bundle.d_b,
               Dhyd=st.d_s,
               layout=st.layout,
               l_b=st.l_b,
               N_ss=5,
               eta_wall=st.shell.state[1].eta*ones(st.Ncell),
               N_passes=st.N_passes,
               N_baffles_d=st.N_baffles_d,
               N_baffles=st.N_baffles,
               ff=1.0,
               b_cut=0.25),
              redeclare VIP.Design.Materials.S_AISI_1040 Material_t,
              redeclare VIP.Design.Materials.S_AISI_1040 Material_s,
               redeclare function bundle_clearance =
               VIP.Design.Miscellanea.Shell_clearance.SRFH,
               thick_t=2e-3,
               thick_s=10e-3,
               p_s_in=301435,
               p_s_out=301435,
               p_t_in=100000,
               p_t_out=100000)
              annotation (Placement(transformation(extent={{-80,-68},{68,58}})));

           annotation (Placement(transformation(extent={{-108,-74},{88,66}})));
          end Coulson_Bell_Delaware;

          model Aspen_Bell_Delaware_upsidedown
            "Verification with the results given by Aspen using the Bell Delaware method. We switch the fluid allocation here."

            replaceable package Medium_t = VIP.Media.RefProp.Methanol
              "Medium model";
            replaceable package Medium_s = VIP.Media.CoolProp.Water
              "Medium model";
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
               redeclare Heat_transfer.Shell.single_phase_Bell_Delaware hT_shell(pitch_f = st.pitch_f, Dhyd_o = st.bundle.Dhyd_o,
               N_tubes = st.N_tubes, d_s = st.d_s, d_b = st.bundle.d_b, Dhyd = st.d_s,
               layout = st.layout, l_b = st.l_b, N_ss = 5,
               N_passes = st.N_passes, ff =  1.0, b_cut = 0.1588),
               redeclare Pressure_drops.Shell.single_phase_Bell_Delaware dp_shell(
               pitch_f=st.pitch_f,
               Dhyd_o=st.bundle.Dhyd_o,
               N_tubes=st.N_tubes,
               d_s=st.d_s,
               d_b=st.bundle.d_b,
               Dhyd=st.d_s,
               layout=st.layout,
               l_b=st.l_b,
               N_ss=5,
               eta_wall=st.shell.state[1].eta*ones(st.Ncell),
               N_passes=st.N_passes,
               N_baffles_d=st.N_baffles_d,
               N_baffles=st.N_baffles,
               ff=1.0,
               b_cut=0.1588),
              redeclare VIP.Design.Materials.S_AISI_1040 Material_t,
              redeclare VIP.Design.Materials.S_AISI_1040 Material_s,
               redeclare function bundle_clearance =
               VIP.Design.Miscellanea.Shell_clearance.Fixed_Utube,
               thick_t=2e-3,
               thick_s=10e-3,
               p_s_in=100000,
               p_s_out=100000,
               p_t_in=301435,
               p_t_out=301435)
              annotation (Placement(transformation(extent={{-80,-68},{68,58}})));

           annotation (Placement(transformation(extent={{-108,-74},{88,66}})));
          end Aspen_Bell_Delaware_upsidedown;

          model Aspen_Bell_Delaware
            "Verification with the results given by Aspen using the Bell Delaware method."

            replaceable package Medium_s = VIP.Media.RefProp.Methanol
              "Medium model";
            replaceable package Medium_t = VIP.Media.CoolProp.Water
              "Medium model";
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
               redeclare Heat_transfer.Shell.single_phase_Bell_Delaware hT_shell(pitch_f = st.pitch_f, Dhyd_o = st.bundle.Dhyd_o,
               N_tubes = st.N_tubes, d_s = st.d_s, d_b = st.bundle.d_b, Dhyd = st.d_s,
               layout = st.layout, l_b = st.l_b, N_ss = 5,
               N_passes = st.N_passes, ff =  1.0, b_cut = 0.1588),
               redeclare Pressure_drops.Shell.single_phase_Bell_Delaware dp_shell(
               pitch_f=st.pitch_f,
               Dhyd_o=st.bundle.Dhyd_o,
               N_tubes=st.N_tubes,
               d_s=st.d_s,
               d_b=st.bundle.d_b,
               Dhyd=st.d_s,
               layout=st.layout,
               l_b=st.l_b,
               N_ss=5,
               eta_wall=st.shell.state[1].eta*ones(st.Ncell),
               N_passes=st.N_passes,
               N_baffles_d=st.N_baffles_d,
               N_baffles=st.N_baffles,
               ff=1.0,
               b_cut=0.1588),
              redeclare VIP.Design.Materials.S_AISI_1040 Material_t,
              redeclare VIP.Design.Materials.S_AISI_1040 Material_s,
               redeclare function bundle_clearance =
               VIP.Design.Miscellanea.Shell_clearance.Fixed_Utube,
               p_s_in=301435,
               p_s_out=301435,
               p_t_in=100000,
               p_t_out=100000,
               redeclare function cost = VIP.Design.Miscellanea.Cost.CS_Hall)
              annotation (Placement(transformation(extent={{-80,-68},{68,58}})));

           annotation (Placement(transformation(extent={{-108,-74},{88,66}})), Icon(
                  coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{
                      100,100}}), graphics));
          end Aspen_Bell_Delaware;

          model DYNDES_Bell_Delaware
            "Verification with the results given by the DYNDES simulation tool using the Bell Delaware method."

            replaceable package Medium_s = VIP.Media.CoolProp.Cyclopentane
              "Medium model";
            replaceable package Medium_t = VIP.Media.CoolProp.Cyclopentane
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
              redeclare Heat_transfer.Shell.single_phase_Bell_Delaware   hT_shell(pitch_f = st.pitch_f, Dhyd_o = st.bundle.Dhyd_o,
              N_tubes = st.N_tubes, d_s = st.d_s, d_b = st.bundle.d_b, Dhyd = st.d_s,
              layout = st.layout, l_b = st.l_b, N_ss = 5,
              N_passes = st.N_passes, ff =  1.0, b_cut = 0.25),
              redeclare Pressure_drops.Shell.single_phase_Bell_Delaware dp_shell(
              pitch_f=st.pitch_f,
              Dhyd_o=st.bundle.Dhyd_o,
              N_tubes=st.N_tubes,
              d_s=st.d_s,
              d_b=st.bundle.d_b,
              Dhyd=st.d_s,
              layout=st.layout,
              l_b=st.l_b,
              N_ss=5,
              eta_wall=st.shell.state[1].eta*ones(st.Ncell),
              N_baffles=st.N_baffles,
              N_baffles_d=st.N_baffles_d,
              N_passes=st.N_passes,
              ff=1.0,
              b_cut=0.25),
              redeclare function bundle_clearance =
              VIP.Design.Miscellanea.Shell_clearance.SRFH,
              redeclare VIP.Design.Materials.S_AISI_1040 Material_t,
              redeclare VIP.Design.Materials.S_AISI_1040 Material_s,
              thick_t=2e-3,
              thick_s=10e-3,
              p_s_in=110000,
              p_s_out=110000,
              p_t_in=3000000,
              p_t_out=3000000)
              annotation (Placement(transformation(extent={{-80,-68},{68,58}})));

           annotation (Placement(transformation(extent={{-108,-74},{88,66}})));
          end DYNDES_Bell_Delaware;
        end Single_phase;

        package Condensation "Tests for condensers"
          model Coulson_Kern
            "Verification with the results given by Coulson using the Kern method."
            import VIP;

            replaceable package Medium_s =
                VIP.Media.FluidProp.Mixture_pentane_ethane "Medium model";
            replaceable package Medium_t = VIP.Media.CoolProp.Water
              "Medium model";
            parameter Modelica.SIunits.SpecificEnthalpy h_s_in= Medium_s.specificEnthalpy_pT(1e6,333.15)
              "Inlet specific enthalpy shell side";
            parameter Modelica.SIunits.SpecificEnthalpy h_s_out= Medium_s.specificEnthalpy_pT(1e6,318.15)
              "Outlet specific enthalpy shell side";
            parameter Modelica.SIunits.SpecificEnthalpy h_t_in = Medium_t.specificEnthalpy_pT(1.2e5,303.15)
              "Inlet specific enthalpy tube side";
            parameter Modelica.SIunits.SpecificEnthalpy h_t_out = Medium_t.specificEnthalpy_pT(1.2e5,313.15)
              "Outlet specific enthalpy tube side";

            Components.HEX.shell_and_tube  st(
              redeclare package Medium_s = Medium_s,
              redeclare package Medium_t = Medium_t,
              Dhyd=16e-3,
              l=4.88,
              pitch_f=1.25,
              layout=1,
              N_passes=2,
              N_baffles=4,
              N_baffles_d = 9,
              h_s_in=h_s_in,
              h_s_out=h_s_out,
              h_t_in=h_t_in,
              h_t_out=h_t_out,
              m_s = 12.5,
              m_t = st.m_s*(h_s_in - h_s_out)/(h_t_out - h_t_in),
              DTML = 20,
              redeclare function bundle_clearance =
              VIP.Design.Miscellanea.Shell_clearance.SRFH,
              redeclare VIP.Design.Materials.S_AISI_1040 Material_t,
              redeclare VIP.Design.Materials.S_AISI_1040 Material_s,
              pin_s=1,
              pin_t=-1,
              thick_t=2e-3,
              thick_s=10e-3,
              p_s_in=1e6,
              p_s_out=1e6,
              p_t_in=1.2e5,
              p_t_out=1.2e5)
              annotation (Placement(transformation(extent={{-80,-68},{68,58}})));

           annotation (Placement(transformation(extent={{-108,-74},{88,66}})));
          end Coulson_Kern;

          model Coulson_Kern_2
            "Verification with the results given by Coulson using the Kern method."
            import VIP;

            replaceable package Medium_s =
                VIP.Media.CoolProp.Mixture_propane_butane "Medium model";
            replaceable package Medium_t = VIP.Media.CoolProp.Water
              "Medium model";
            parameter Modelica.SIunits.MassFlowRate m_s = 12.5
              "Mass flow shell side";
            parameter Modelica.SIunits.SpecificEnthalpy h_s_in= Medium_s.specificEnthalpy_pT(1e5,333.15)
              "Inlet specific enthalpy shell side";
            parameter Modelica.SIunits.SpecificEnthalpy h_s_out= Medium_s.specificEnthalpy_pT(1e5,318.15)
              "Outlet specific enthalpy shell side";
            parameter Modelica.SIunits.AbsolutePressure p_s_in = 10e5
              "Inlet pressure shell side";
            parameter Modelica.SIunits.AbsolutePressure p_s_out = 10e5
              "Outlet pressure shell side";
            parameter Modelica.SIunits.SpecificEnthalpy h_t_in = Medium_t.specificEnthalpy_pT(p_t_in,303.15)
              "Inlet specific enthalpy tube side";
            parameter Modelica.SIunits.SpecificEnthalpy h_t_out = Medium_t.specificEnthalpy_pT(p_t_out,313.15)
              "Outlet specific enthalpy tube side";
            parameter Modelica.SIunits.AbsolutePressure p_t_in = 2e5
              "Inlet pressure tube side";
            parameter Modelica.SIunits.AbsolutePressure p_t_out = 2e5
              "Outlet pressure tube side";
          //   parameter Medium_s.SaturationProperties sat =  Medium_s.setSat_p(p_s_in)
          //     "Saturation properties";
          //   parameter Medium_s.ThermodynamicState state_l = Medium_s.setState_ph(p_s_in, sat.hl)
          //     "Thermodynamic state in saturated liquid";
          //   parameter Medium_s.ThermodynamicState state_v = Medium_s.setState_ph(p_s_in, sat.hv)
          //     "Thermodynamic state in saturated vapour";

          //   Components.HEX.shell_and_tube  st(
          //     redeclare package Medium_s = Medium_s,
          //     redeclare package Medium_t = Medium_t,
          //     Dhyd=16e-3,
          //     l=4.88,
          //     pitch_f=1.25,
          //     layout=1,
          //     N_passes=2,
          //     N_baffles=6,
          //     h_s_in=h_s_in,
          //     h_s_out=h_s_out,
          //     h_t_in=h_t_in,
          //     h_t_out=h_t_out,
          //     m_s = m_s,
          //     m_t = st.m_s*(h_s_in - h_s_out)/(h_t_out - h_t_in),
          //     redeclare function bundle_clearance =
          //       VIP.Miscellanea.Shell_clearance.SRFH,
          //     redeclare VIP.Materials.S_AISI_1040 Material_t,
          //     redeclare VIP.Materials.S_AISI_1040 Material_s,
          //     pin_s=1,
          //     pin_t=-1,
          //     thick_t=2e-3,
          //     thick_s=10e-3,
          //     N_heads=1,
          //     p_s_in=p_s_in,
          //     p_s_out=p_s_out,
          //     p_t_in=p_t_in,
          //     p_t_out=p_t_out,
          //     ht_t_f1=6e3,
          //     ht_s_f1=6e3,
          //     N_tubes(start=500, fixed=true))
           //   annotation (Placement(transformation(extent={{-80,-68},{68,58}})));
          //     redeclare Heat_transfer.Shell.condensation_Kern hT_shell(Dhyd_o=st.bundle.Dhyd_o, l = st.l,
          //        N_tubes=st.N_tubes, pitch_f = st.pitch_f, d_b = st.bundle.d_b, sat = sat, state_l = state_l,
          //        q_tilde=st.shell.qdot/st.shell.At),

           annotation (Placement(transformation(extent={{-108,-74},{88,66}})));
          end Coulson_Kern_2;

          model Aspen_Kern
            "Verification with the results given by Aspen using the Kern method."
            import VIP;

            replaceable package Medium_s = VIP.Media.RefProp.MM "Medium model";
            replaceable package Medium_t = VIP.Media.CoolProp.Water
              "Medium model";
            parameter Modelica.SIunits.MassFlowRate m_s = 1.5365
              "Mass flow shell side";
            parameter Modelica.SIunits.SpecificEnthalpy h_s_in= Medium_s.specificEnthalpy_pT(p_s_in,99.1 +273.15)
              "Inlet specific enthalpy shell side";
            parameter Modelica.SIunits.SpecificEnthalpy h_s_out= sat.hl
              "Outlet specific enthalpy shell side";
            parameter Modelica.SIunits.AbsolutePressure p_s_in = 0.5e5
              "Inlet pressure shell side";
            parameter Modelica.SIunits.AbsolutePressure p_s_out = 0.5e5
              "Outlet pressure shell side";
            parameter Modelica.SIunits.SpecificEnthalpy h_t_in = Medium_t.specificEnthalpy_pT(p_t_in,313.15)
              "Inlet specific enthalpy tube side";
            parameter Modelica.SIunits.SpecificEnthalpy h_t_out = Medium_t.specificEnthalpy_pT(p_t_out,333.15)
              "Outlet specific enthalpy tube side";
            parameter Modelica.SIunits.AbsolutePressure p_t_in = 2e5
              "Inlet pressure tube side";
            parameter Modelica.SIunits.AbsolutePressure p_t_out = 2e5
              "Outlet pressure tube side";
            parameter Medium_s.SaturationProperties sat =  Medium_s.setSat_p(p_s_in)
              "Saturation properties";
            parameter Medium_s.ThermodynamicState state_l = Medium_s.setState_ph(p_s_in, sat.hl)
              "Thermodynamic state in saturated liquid";
            parameter Medium_s.ThermodynamicState state_v = Medium_s.setState_ph(p_s_in, sat.hv)
              "Thermodynamic state in saturated vapor";

             Components.HEX.shell_and_tube  st(
               redeclare package Medium_s = Medium_s,
               redeclare package Medium_t = Medium_t,
               Dhyd=17.05e-3,
               l=1.6887,
               pitch_f=1.25,
               layout=1,
               N_passes=4,
               N_baffles=5,
               N_baffles_d=6,
               h_s_in=h_s_in,
               h_s_out=h_s_out,
               h_t_in=h_t_in,
               h_t_out=h_t_out,
               m_s = m_s,
               m_t = st.m_s*(h_s_in - h_s_out)/(h_t_out - h_t_in),
               redeclare Heat_transfer.Shell.condensation_Kern hT_shell(
               Dhyd_o = st.bundle.Dhyd_o, l = st.l, N_tubes=st.N_tubes, shear_vapor = true,
               pitch_f = st.pitch_f, state_l = state_l, state_v = state_v),
               redeclare Pressure_drops.Shell.condensation_Kern dp_shell(
               l = st.l/st.N_baffles_d/st.N_passes, d_s = st.d_s, X = 0.4,
               Dhyd = st.bundle.Dhyd_o, l_b = st.l_b, state_l = state_l,
               state_v = state_v, eta_wall = st.shell.state[1].eta*ones(st.Ncell)),
               redeclare function bundle_clearance =
               VIP.Design.Miscellanea.Shell_clearance.Fixed_Utube,
              redeclare VIP.Design.Materials.S_AISI_1040 Material_t,
              redeclare VIP.Design.Materials.S_AISI_1040 Material_s,
               pin_s=1,
               pin_t=-1,
               thick_t=2e-3,
               thick_s=10e-3,
               p_s_in=p_s_in,
               p_s_out=p_s_out,
               p_t_in=p_t_in,
               p_t_out=p_t_out,
               ht_t_f1=1e12,
               ht_s_f1=1e12,
               N_tubes(start=300, fixed=true))
              annotation (Placement(transformation(extent={{-80,-68},{68,58}})));

           annotation (Placement(transformation(extent={{-108,-74},{88,66}})));
          end Aspen_Kern;

          model Aspen_Bell_Delaware
            "Verification with the results given by Aspen using the Bell Delaware method."
            import VIP;

            replaceable package Medium_s = VIP.Media.RefProp.MM "Medium model";
            replaceable package Medium_t = VIP.Media.CoolProp.Water
              "Medium model";
            parameter Modelica.SIunits.MassFlowRate m_s = 1.5365
              "Mass flow shell side";
            parameter Modelica.SIunits.SpecificEnthalpy h_s_in= Medium_s.specificEnthalpy_pT(p_s_in,99.1 +273.15)
              "Inlet specific enthalpy shell side";
            parameter Modelica.SIunits.SpecificEnthalpy h_s_out= sat.hl
              "Outlet specific enthalpy shell side";
            parameter Modelica.SIunits.AbsolutePressure p_s_in = 0.5e5
              "Inlet pressure shell side";
            parameter Modelica.SIunits.AbsolutePressure p_s_out = 0.5e5
              "Outlet pressure shell side";
            parameter Modelica.SIunits.SpecificEnthalpy h_t_in = Medium_t.specificEnthalpy_pT(p_t_in,313.15)
              "Inlet specific enthalpy tube side";
            parameter Modelica.SIunits.SpecificEnthalpy h_t_out = Medium_t.specificEnthalpy_pT(p_t_out,333.15)
              "Outlet specific enthalpy tube side";
            parameter Modelica.SIunits.AbsolutePressure p_t_in = 2e5
              "Inlet pressure tube side";
            parameter Modelica.SIunits.AbsolutePressure p_t_out = 2e5
              "Outlet pressure tube side";
            parameter Medium_s.SaturationProperties sat =  Medium_s.setSat_p(p_s_in)
              "Saturation properties";
            parameter Medium_s.ThermodynamicState state_l = Medium_s.setState_ph(p_s_in, sat.hl)
              "Thermodynamic state in saturated liquid";
            parameter Medium_s.ThermodynamicState state_v = Medium_s.setState_ph(p_s_in, sat.hv)
              "Thermodynamic state in saturated vapor";

             Components.HEX.shell_and_tube  st(
               redeclare package Medium_s = Medium_s,
               redeclare package Medium_t = Medium_t,
               Dhyd=17.05e-3,
               l=1.8,
               pitch_f=1.25,
               layout=1,
               N_passes=4,
               N_baffles=5,
               N_baffles_d=6,
               h_s_in=h_s_in,
               h_s_out=h_s_out,
               h_t_in=h_t_in,
               h_t_out=h_t_out,
               m_s = m_s,
               m_t = st.m_s*(h_s_in - h_s_out)/(h_t_out - h_t_in),
               redeclare Heat_transfer.Shell.condensation_Bell_Delaware hT_shell(
               pitch_f = st.pitch_f, Dhyd_o = st.bundle.Dhyd_o, N_tubes = st.N_tubes,
               d_s = st.d_s, d_b = st.bundle.d_b, Dhyd = st.d_s, l = st.l, bts = 3.18e-3,
               layout = st.layout, l_b = st.l_b, N_ss = 2, state_l = state_l,
               state_v = state_v, N_passes = st.N_passes, ff =  1.0, b_cut = 0.3912),
               redeclare Pressure_drops.Shell.condensation_Bell_Delaware dp_shell(
               pitch_f = st.pitch_f, Dhyd_o = st.bundle.Dhyd_o, N_tubes = st.N_tubes,
               d_s = st.d_s, d_b = st.bundle.d_b, Dhyd = st.d_s, layout = st.layout,
               l_b = st.l_b, N_ss = 2, N_passes = st.N_passes,
               N_baffles_d = st.N_baffles_d, N_baffles = st.N_baffles, ff = 1.0,
               b_cut = 0.3912, state_l = state_l, state_v = state_v, bts = 3.18e-3),
               redeclare function bundle_clearance =
               VIP.Design.Miscellanea.Shell_clearance.Fixed_Utube,
              redeclare VIP.Design.Materials.S_AISI_1040 Material_t,
              redeclare VIP.Design.Materials.S_AISI_1040 Material_s,
               pin_s=1,
               pin_t=-1,
               thick_t=2e-3,
               thick_s=10e-3,
               p_s_in=p_s_in,
               p_s_out=p_s_out,
               p_t_in=p_t_in,
               p_t_out=p_t_out,
               ht_t_f1=1e12,
               ht_s_f1=1e12,
               N_tubes(start=100, fixed=true))
              annotation (Placement(transformation(extent={{-80,-68},{68,58}})));

           annotation (Placement(transformation(extent={{-108,-74},{88,66}})),
              experiment,
              __Dymola_experimentSetupOutput);
          end Aspen_Bell_Delaware;
        end Condensation;
      end Shell_and_tube;

      package Flat_plate "Tests for flat plate heat exchangers"
        model Coulson "Verification with the results given by Coulson et al."
          import VIP;

          replaceable package Medium_hot = VIP.Media.RefProp.Methanol
            "Medium model";
          replaceable package Medium_cold = VIP.Media.CoolProp.Water_TTSE
            "Medium model";
          parameter Medium_hot.ThermodynamicState hot_in = Medium_hot.setState_pT(3.2e5, 95 + 273.15)
            "Thermodynamic state at the inlet of the hot side";
          parameter Medium_hot.ThermodynamicState hot_out = Medium_hot.setState_pT(3.2e5, 40 + 273.15)
            "Thermodynamic state at the outlet of the hot side";
           parameter Medium_cold.ThermodynamicState cold_in = Medium_cold.setState_pT(1e5, 25 + 273.15)
            "Thermodynamic state at the inlet of the hot side";
           parameter Medium_cold.ThermodynamicState cold_out = Medium_cold.setState_pT(1e5, 40 + 273.15)
            "Thermodynamic state at the outlet of the hot side";

        Components.HEX.Flat_plate fp(
          redeclare package Medium_hot = Medium_hot,
          redeclare package Medium_cold = Medium_cold,
          h_hot_in=hot_in.h,
          h_hot_out=hot_out.h,
          h_cold_in=cold_in.h,
          h_cold_out=cold_out.h,
          mdot_cold=fp.mdot_hot*(hot_in.h - hot_out.h)/(cold_out.h - cold_in.h),
          b=3e-3,
          w=0.5,
          X=9e-3,
          ht_hot_f1=1e4,
          ht_cold_f1=6e3,
          thick=0.75e-3,
          redeclare Heat_transfer.Plates.single_phase_Coulson ht_cold,
          redeclare Heat_transfer.Plates.single_phase_Coulson ht_hot,
            redeclare VIP.Design.Materials.T_R50250 material,
          d_pt=0.1,
          N_cell_pc=2,
          p_hot_in=hot_in.p,
          p_hot_out=hot_out.p,
          p_cold_in=cold_in.p,
          p_cold_out=cold_out.p,
          redeclare function cost = VIP.Design.Miscellanea.Cost.FP_Rafferty,
            mdot_hot=100/3.6,
            U_guess=600,
            N_ch_p=60,
            redeclare VIP.Design.Pressure_drops.Plates.single_phase_Coulson
              dp_hot,
            redeclare VIP.Design.Pressure_drops.Plates.single_phase_Coulson
              dp_cold,
            beta=0.87266462599716)
          annotation (Placement(transformation(extent={{-82,-66},{86,38}})));

         annotation (Placement(transformation(extent={{-108,-74},{88,66}})),
            experiment(__Dymola_NumberOfIntervals=1, Tolerance=1e-006),
            __Dymola_experimentSetupOutput);
        end Coulson;

        model Rossetto "Verification with the results given by Rossetto et al."
          import VIP;

          replaceable package Medium_hot = VIP.Media.CoolProp.Water_TTSE
            "Medium model";
          replaceable package Medium_cold = VIP.Media.CoolProp.Water_TTSE
            "Medium model";
          parameter Medium_hot.ThermodynamicState hot_in = Medium_hot.setState_pT(1e5, 60 + 273.15)
            "Thermodynamic state at the inlet of the hot side";
          parameter Medium_hot.ThermodynamicState hot_out = Medium_hot.setState_pT(1e5, 33 + 273.15)
            "Thermodynamic state at the outlet of the hot side";
           parameter Medium_cold.ThermodynamicState cold_in = Medium_cold.setState_pT(1e5, 20 + 273.15)
            "Thermodynamic state at the inlet of the hot side";
           parameter Medium_cold.ThermodynamicState cold_out = Medium_cold.setState_pT(1e5, 47 + 273.15)
            "Thermodynamic state at the outlet of the hot side";

          Components.HEX.Flat_plate fp(
            redeclare package Medium_hot = Medium_hot,
            redeclare package Medium_cold = Medium_cold,
            h_hot_in=hot_in.h,
            h_hot_out=hot_out.h,
            h_cold_in=cold_in.h,
            h_cold_out=cold_out.h,
            mdot_cold=fp.mdot_hot*(hot_in.h - hot_out.h)/(cold_out.h - cold_in.h),
            X=1.025,
            redeclare function cost = VIP.Design.Miscellanea.Cost.FP_Rafferty,
            redeclare VIP.Design.Materials.SS_AISI_410 material,
            thick=0.4e-3,
            b=3.4e-3,
            w=0.65,
            d_pt=0.25,
            mdot_hot=140,
            ht_hot_f1=2e4,
            ht_cold_f1=1e12,
            N_cell_pc=2,
            p_hot_in=hot_in.p,
            p_hot_out=hot_out.p,
            p_cold_in=cold_in.p,
            p_cold_out=cold_out.p,
            U_guess=600,
            DTML=13,
            N_ch_p=204,
            beta=0.78539816339745)
            annotation (Placement(transformation(extent={{-82,-66},{86,38}})));

         annotation (Placement(transformation(extent={{-108,-74},{88,66}})),
            experiment(__Dymola_NumberOfIntervals=1),
            __Dymola_experimentSetupOutput);
        end Rossetto;

        model Aspen_eva
          "Verification with the results given by Aspen. We test evaporation in this example."
          import VIP;

           replaceable package Medium_hot =
               VIP.Media.CoolProp.THERM66 "Medium model";
          replaceable package Medium_cold = VIP.Media.RefProp.MM "Medium model";
          parameter Medium_cold.SaturationProperties sat =  Medium_cold.setSat_p(14.6e5)
            "Saturation properties";
          parameter Medium_hot.ThermodynamicState hot_in = Medium_hot.setState_pT(1e5, 345 + 273.15)
            "Thermodynamic state at the inlet of the hot side";
          parameter Medium_hot.ThermodynamicState hot_out = Medium_hot.setState_pT(1e5, 231.9 + 273.15)
            "Thermodynamic state at the outlet of the hot side";
           parameter Medium_cold.ThermodynamicState cold_in = Medium_cold.setState_pT(14.6e5, 216.9 + 273.15)
            "Thermodynamic state at the inlet of the hot side";
           parameter Medium_cold.ThermodynamicState cold_out = Medium_cold.setState_pT(14.6e5, 325 + 273.15)
            "Thermodynamic state at the outlet of the hot side";
          parameter Medium_cold.ThermodynamicState state_l = Medium_cold.setState_ph(14.6e5, sat.hl)
            "Thermodynamic state in saturated liquid";
          parameter Medium_cold.ThermodynamicState state_v = Medium_cold.setState_ph(14.6e5, sat.hv)
            "Thermodynamic state in saturated vapor";

            Components.HEX.Flat_plate fp(
              redeclare package Medium_hot = Medium_hot,
              redeclare package Medium_cold = Medium_cold,
              h_hot_in=hot_in.h,
              h_hot_out=hot_out.h,
              h_cold_in=cold_in.h,
              h_cold_out=cold_out.h,
              mdot_cold=fp.mdot_hot*(hot_in.h - hot_out.h)/(cold_out.h - cold_in.h),
              redeclare function cost = Miscellanea.Cost.FP_Rafferty,
              thick=0.8e-3,
              d_pt=32e-3,
              mdot_hot=0.19,
              ht_hot_f1=1e9,
              ht_cold_f1=1e9,
              p_hot_in=hot_in.p,
              p_hot_out=hot_out.p,
              p_cold_in=cold_in.p,
              p_cold_out=cold_out.p,
              t_hot_in=hot_in.T,
              t_hot_out=hot_out.T,
              X = 1,
              N_ch_p=14,
              U_guess=600,
              w=0.1325,
              redeclare Heat_transfer.Plates.evaporation_Martin ht_cold(p_in=
              fp.p_cold_in, state_l = state_l, state_v = state_v, At = fp.At, l = fp.l,
              beta = fp.beta, qdot = fp.plate.qdot_cold,
              qdot_tilde_start = fp.qdot/(2*fp.w*fp.l_start*fp.N_ch*fp.N_cell_pc)),
              redeclare Miscellanea.topology_PHE.two_pass_two_pass tpg_hot(N_ch_p=
              fp.N_ch_p, stype = 1),
              redeclare Miscellanea.topology_PHE.two_pass_two_pass tpg_cold(N_ch_p=
              fp.N_ch_p, stype= 2),
            redeclare VIP.Design.Materials.SS_AISI_304 material,
              N_cell_pc=8,
              b = 2.3e-3,
            redeclare VIP.Design.Pressure_drops.Plates.evaporation_Martin dp_cold(
              state_l=state_l,
              state_v=state_v,
              beta=fp.beta),
              redeclare Miscellanea.check_velocity check_hot(T_sat = 500),
              redeclare Miscellanea.check_velocity check_cold(umin = max(fp.ht_cold.u)),
            beta=1.0471975511966)
            annotation (Placement(transformation(extent={{-82,-66},{86,38}})));

         annotation (Placement(transformation(extent={{-108,-74},{88,66}})),
            experiment(__Dymola_NumberOfIntervals=1),
            __Dymola_experimentSetupOutput);
        end Aspen_eva;

        model Rossetto_series
          "Verification with the results given by Rossetto et al. Series configuration."
          import VIP;

          replaceable package Medium_hot = VIP.Media.CoolProp.Water_TTSE
            "Medium model";
          replaceable package Medium_cold = VIP.Media.CoolProp.Water_TTSE
            "Medium model";
          parameter Medium_hot.ThermodynamicState hot_in = Medium_hot.setState_pT(1e5, 60 + 273.15)
            "Thermodynamic state at the inlet of the hot side";
          parameter Medium_hot.ThermodynamicState hot_out = Medium_hot.setState_pT(1e5, 33 + 273.15)
            "Thermodynamic state at the outlet of the hot side";
           parameter Medium_cold.ThermodynamicState cold_in = Medium_cold.setState_pT(1e5, 20 + 273.15)
            "Thermodynamic state at the inlet of the hot side";
           parameter Medium_cold.ThermodynamicState cold_out = Medium_cold.setState_pT(1e5, 47 + 273.15)
            "Thermodynamic state at the outlet of the hot side";

          Components.HEX.Flat_plate fp(
            redeclare package Medium_hot = Medium_hot,
            redeclare package Medium_cold = Medium_cold,
            h_hot_in=hot_in.h,
            h_hot_out=hot_out.h,
            h_cold_in=cold_in.h,
            h_cold_out=cold_out.h,
            mdot_cold=fp.mdot_hot*(hot_in.h - hot_out.h)/(cold_out.h - cold_in.h),
            X=1.025,
            redeclare function cost = VIP.Design.Miscellanea.Cost.FP_Rafferty,
            redeclare VIP.Design.Materials.SS_AISI_410 material,
            thick=0.4e-3,
            b=3.4e-3,
            w=0.65,
            d_pt=0.25,
            mdot_hot=14,
            ht_hot_f1=2e4,
            ht_cold_f1=1e12,
            N_cell_pc=2,
            p_hot_in=hot_in.p,
            p_hot_out=hot_out.p,
            p_cold_in=cold_in.p,
            p_cold_out=cold_out.p,
            t_hot_in=hot_in.T,
            t_hot_out=hot_out.T,
            t_cold_in=cold_in.T,
            t_cold_out=cold_out.T,
            U_guess=600,
            DTML=13,
            N_ch_p=20,
            beta=0.78539816339745,
            redeclare VIP.Design.Miscellanea.topology_PHE.series tpg_hot,
            redeclare VIP.Design.Miscellanea.topology_PHE.series tpg_cold)
            annotation (Placement(transformation(extent={{-82,-66},{86,38}})));

         annotation (Placement(transformation(extent={{-108,-74},{88,66}})),
            experiment(__Dymola_NumberOfIntervals=1),
            __Dymola_experimentSetupOutput);
        end Rossetto_series;

        model Rossetto_parallel
          "Verification with the results given by Rossetto et al. 2_2 configuration."
          import VIP;

          replaceable package Medium_hot = VIP.Media.CoolProp.Water_TTSE
            "Medium model";
          replaceable package Medium_cold = VIP.Media.CoolProp.Water_TTSE
            "Medium model";
          parameter Medium_hot.ThermodynamicState hot_in = Medium_hot.setState_pT(1e5, 60 + 273.15)
            "Thermodynamic state at the inlet of the hot side";
          parameter Medium_hot.ThermodynamicState hot_out = Medium_hot.setState_pT(1e5, 33 + 273.15)
            "Thermodynamic state at the outlet of the hot side";
           parameter Medium_cold.ThermodynamicState cold_in = Medium_cold.setState_pT(1e5, 20 + 273.15)
            "Thermodynamic state at the inlet of the hot side";
           parameter Medium_cold.ThermodynamicState cold_out = Medium_cold.setState_pT(1e5, 47 + 273.15)
            "Thermodynamic state at the outlet of the hot side";

          Components.HEX.Flat_plate fp(
            redeclare package Medium_hot = Medium_hot,
            redeclare package Medium_cold = Medium_cold,
            h_hot_in=hot_in.h,
            h_hot_out=hot_out.h,
            h_cold_in=cold_in.h,
            h_cold_out=cold_out.h,
            mdot_cold=fp.mdot_hot*(hot_in.h - hot_out.h)/(cold_out.h - cold_in.h),
            X=1.025,
            redeclare function cost = VIP.Design.Miscellanea.Cost.FP_Rafferty,
            redeclare VIP.Design.Materials.SS_AISI_410 material,
            thick=0.4e-3,
            b=3.4e-3,
            w=0.65,
            d_pt=0.25,
            ht_hot_f1=2e4,
            ht_cold_f1=1e12,
            p_hot_in=hot_in.p,
            p_hot_out=hot_out.p,
            p_cold_in=cold_in.p,
            p_cold_out=cold_out.p,
            t_hot_in=hot_in.T,
            t_hot_out=hot_out.T,
            t_cold_in=cold_in.T,
            t_cold_out=cold_out.T,
            U_guess=600,
            DTML=13,
            N_ch_p=3,
            mdot_hot=1.4,
            redeclare VIP.Design.Miscellanea.topology_PHE.two_pass_two_pass
              tpg_hot(N_ch_p=fp.N_ch_p, stype=1),
            redeclare VIP.Design.Miscellanea.topology_PHE.two_pass_two_pass
              tpg_cold(N_ch_p=fp.N_ch_p, stype=2),
            N_cell_pc=2,
            beta=0.78539816339745)
            annotation (Placement(transformation(extent={{-82,-66},{86,38}})));

         annotation (Placement(transformation(extent={{-108,-74},{88,66}})),
            experiment(__Dymola_NumberOfIntervals=1),
            __Dymola_experimentSetupOutput);
        end Rossetto_parallel;

        model Aspen_rec
          "Verification with the results given by Aspen. We test recuperation in this example."
          import VIP;

          package Medium_hot = VIP.Media.FluidProp.MM "Medium model";
          package Medium_cold = VIP.Media.FluidProp.MM "Medium model";
          parameter Medium_hot.ThermodynamicState hot_in = Medium_hot.setState_pT(0.33e5, 276.9 + 273.15)
            "Thermodynamic state at the inlet of the hot side";
          parameter Medium_hot.ThermodynamicState hot_out = Medium_hot.setState_pT(0.33e5, 71.4 + 273.15)
            "Thermodynamic state at the outlet of the hot side";
           parameter Medium_cold.ThermodynamicState cold_in = Medium_cold.setState_pT(14.6e5, 51.4 + 273.15)
            "Thermodynamic state at the inlet of the hot side";
           parameter Medium_cold.ThermodynamicState cold_out = Medium_cold.setState_pT(14.6e5, 222.9 + 273.15)
            "Thermodynamic state at the outlet of the hot side";
          parameter Medium_cold.SaturationProperties sat =  Medium_cold.setSat_p(14.6e5)
            "Saturation properties";
          parameter Medium_cold.ThermodynamicState state_l = Medium_cold.setState_ph(14.6e5, sat.hl)
            "Thermodynamic state in saturated liquid";
          parameter Medium_cold.ThermodynamicState state_v = Medium_cold.setState_ph(14.6e5, sat.hv)
            "Thermodynamic state in saturated vapor";

             Components.HEX.Flat_plate fp(
               redeclare package Medium_hot = Medium_hot,
               redeclare package Medium_cold = Medium_cold,
               h_hot_in=hot_in.h,
               h_hot_out=hot_out.h,
               h_cold_in=cold_in.h,
               h_cold_out=cold_in.h + hot_in.h - hot_out.h,
               mdot_cold=0.15,
               redeclare function cost = Miscellanea.Cost.FP_Rafferty,
               ht_hot_f1=1e9,
               ht_cold_f1=1e9,
               p_hot_in=hot_in.p,
               p_hot_out=hot_out.p,
               p_cold_in=cold_in.p,
               p_cold_out=cold_out.p,
               X = 1,
               U_guess=100,
            redeclare VIP.Design.Materials.SS_AISI_304 material,
               N_ch_p=18,
               thick=0.9e-3,
               b=3.12e-3,
               w=325e-3,
               d_pt=100e-3,
               mdot_hot=0.15,
               redeclare Miscellanea.topology_PHE.two_pass_two_pass tpg_hot(N_ch_p=
               fp.N_ch_p, stype = 1),
               redeclare Miscellanea.topology_PHE.two_pass_two_pass tpg_cold(N_ch_p=
               fp.N_ch_p, stype= 2),
               N_cell_pc=3,
               beta=1.0594148559606)
              annotation (Placement(transformation(extent={{-82,-66},{86,38}})));

         annotation (Placement(transformation(extent={{-108,-74},{88,66}})),
            experiment(__Dymola_NumberOfIntervals=1),
            __Dymola_experimentSetupOutput);
        end Aspen_rec;

        model rec_variant
          "We test a variant for the recuperation in this example."
          import VIP;

          package Medium_hot = VIP.Media.FluidProp.MM "Medium model";
          package Medium_cold = VIP.Media.FluidProp.MM "Medium model";
          parameter Medium_hot.ThermodynamicState hot_in = Medium_hot.setState_pT(0.33e5, 276.9 + 273.15)
            "Thermodynamic state at the inlet of the hot side";
          parameter Medium_hot.ThermodynamicState hot_out = Medium_hot.setState_pT(0.33e5, 71.4 + 273.15)
            "Thermodynamic state at the outlet of the hot side";
           parameter Medium_cold.ThermodynamicState cold_in = Medium_cold.setState_pT(14.6e5, 51.4 + 273.15)
            "Thermodynamic state at the inlet of the hot side";
           parameter Medium_cold.ThermodynamicState cold_out = Medium_cold.setState_pT(14.6e5, 222.9 + 273.15)
            "Thermodynamic state at the outlet of the hot side";
          parameter Medium_cold.SaturationProperties sat =  Medium_cold.setSat_p(14.6e5)
            "Saturation properties";
          parameter Medium_cold.ThermodynamicState state_l = Medium_cold.setState_ph(14.6e5, sat.hl)
            "Thermodynamic state in saturated liquid";
          parameter Medium_cold.ThermodynamicState state_v = Medium_cold.setState_ph(14.6e5, sat.hv)
            "Thermodynamic state in saturated vapor";

             Components.HEX.Flat_plate fp(
               redeclare package Medium_hot = Medium_hot,
               redeclare package Medium_cold = Medium_cold,
               h_hot_in=hot_in.h,
               h_hot_out=hot_out.h,
               h_cold_in=cold_in.h,
               h_cold_out=cold_in.h + hot_in.h - hot_out.h,
               mdot_cold=0.15,
               redeclare function cost = Miscellanea.Cost.FP_Rafferty,
               ht_hot_f1=1e9,
               ht_cold_f1=1e9,
               p_hot_in=hot_in.p,
               p_hot_out=hot_out.p,
               p_cold_in=cold_in.p,
               p_cold_out=cold_out.p,
               X = 1,
               U_guess=100,
            redeclare VIP.Design.Materials.SS_AISI_304 material,
               N_ch_p=18,
               thick=0.9e-3,
               b=2.3e-3,
               w=325e-3,
               d_pt=100e-3,
               mdot_hot=0.15,
               N_cell_pc=3,
               beta=1.0594148559606)
              annotation (Placement(transformation(extent={{-82,-66},{86,38}})));

         annotation (Placement(transformation(extent={{-108,-74},{88,66}})),
            experiment(__Dymola_NumberOfIntervals=1),
            __Dymola_experimentSetupOutput);
        end rec_variant;

      end Flat_plate;

      package Cell_method "Tests of the cell method for various configurations"
        model cell_method_VDI_example
          "Cell method verified with an example from VDI C1 pag. 48"
          parameter Integer N_passes = 2 "Number of cell passes";
          parameter Integer N_baffles = 3 "Number of cell baffles";
          parameter Integer Ncell = N_passes*N_baffles
            "Number of cell elements";
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

        model cell_method_flat_plate
          "Cell method for flat plate heat exchangers"
          parameter Integer layout = 2 "Flow path, 1 = parallel, 2 = series";
          parameter Integer N_passes_d = 10
            "Number of passes hot and cold side";
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
           start=fill(t_cold_start, N_passes_d))
            "Cold stream temperature matrix";
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
              DTML       = Miscellanea.log_mean_delta_T(t_hot_in, t_hot_out, t_cold_in,
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
    end Tests;

    package Components "Library with the design of the components "
      package HEX "Heat exchangers"
        extends VIP.Design.Icons.HEX;
        model shell_and_tube
          "Shell and tube heat exchanger where the hot fluid flows on the shell and enters from the top. The cold fluid enters at the bottom."
          import VIP;
          extends Icons.shell_tube;

          //THE WORKING FLUIDS
          replaceable package Medium_s = VIP.Media.OneRandomOrganicFluid
                                                                     constrainedby
            Modelica.Media.Interfaces.PartialMedium "Medium model shell" annotation(choicesAllMatching = true);
          replaceable package Medium_t = VIP.Media.OneRandomOrganicFluid
                                                                     constrainedby
            Modelica.Media.Interfaces.PartialMedium "Medium model tubes" annotation(choicesAllMatching = true);
          replaceable Materials.void_material Material_t "Material model shell"
                                                                                      annotation(choicesAllMatching = true);
          replaceable Materials.void_material Material_s "Material model shell"
                                                                                       annotation(choicesAllMatching = true);

          //GEOMETRY OF THE HEAT EXCHANGER
           parameter Modelica.SIunits.Length Dhyd "Hydraulic diameter";
           parameter Modelica.SIunits.Length thick_t "Tube thickness";
           parameter Modelica.SIunits.Length thick_s "Shell thickness";
           parameter Modelica.SIunits.Length l "Tube lenght";
           parameter Real pitch_f "Tube pitch";
           parameter Integer layout "Tube layout, 1 = triangular, 2 = squared";
           parameter Integer N_passes = 2 "Number of tube passes";
           parameter Integer N_baffles = 4 "Number of baffles";
           parameter Real pin_s
            "Pin for the heat flow in the shell -1 cold fluid else hot";
           parameter Real pin_t
            "Pin for the heat flow in the tubes -1 cold fluid else hot";
           parameter Modelica.SIunits.CoefficientOfHeatTransfer U_guess = 600
            "Guess value for the heat transfer coefficient";
           parameter Integer N_baffles_d = N_baffles + 1
            "The number of discretization volumes";
           final parameter Integer Ncell = N_baffles_d*N_passes
            "Number of cell elements";
           final parameter Modelica.SIunits.HeatFlowRate qdot = m_s*abs(h_s_in - h_s_out)
            "Heat flow rate";
           parameter Modelica.SIunits.Temp_C  DTML = Miscellanea.log_mean_delta_T(t_s_in,
             t_s_out, t_t_in, t_t_out)
            "Logarithmic mean temperature difference";
           //Boundary conditions at inlet and outlet
           parameter Modelica.SIunits.MassFlowRate m_s "Shell mass flow";
           parameter Modelica.SIunits.SpecificEnthalpy h_s_in
            "Inlet specific enthalpy shell side";
           parameter Modelica.SIunits.SpecificEnthalpy h_s_out
            "Outlet specific enthalpy shell side";
           parameter Modelica.SIunits.AbsolutePressure p_s_in
            "Inlet pressure shell side";
           parameter Modelica.SIunits.AbsolutePressure p_s_out
            "Outlet pressure shell side";
           parameter Modelica.SIunits.MassFlowRate m_t "Tube mass flow";
           parameter Modelica.SIunits.SpecificEnthalpy h_t_in
            "Inlet specific enthalpy tube side";
           parameter Modelica.SIunits.SpecificEnthalpy h_t_out
            "Outlet specific enthalpy tube side";
           parameter Modelica.SIunits.AbsolutePressure p_t_in
            "Inlet pressure tube side";
           parameter Modelica.SIunits.AbsolutePressure p_t_out
            "Outlet pressure tube side";
           parameter Modelica.SIunits.Temperature  t_s_in=
            Medium_s.temperature_ph(p_s_in, h_s_in)
            "Inlet temperature shell side";
           parameter Modelica.SIunits.Temperature  t_s_out=
            Medium_s.temperature_ph(p_s_out, h_s_out)
            "Outlet temperature shell side";
           parameter Modelica.SIunits.Temperature  t_t_in=
            Medium_t.temperature_ph(p_t_in, h_t_in)
            "Inlet temperature tube side";
           parameter Modelica.SIunits.Temperature  t_t_out=
            Medium_t.temperature_ph(p_t_out, h_t_out)
            "Outlet temperature tube side";
           parameter Modelica.SIunits.CoefficientOfHeatTransfer ht_t_f1 = 3e3
            "Tube fouling heat transfer coefficient";
           parameter Modelica.SIunits.CoefficientOfHeatTransfer ht_s_f1 = 5e3
            "Shell fouling heat transfer coefficient";
           Modelica.SIunits.Length d_s "Shell diameter";
           Modelica.SIunits.Length l_b "Baffle lenght";
           Real N_tubes(start = qdot/(DTML*U_guess)/(pi*l*(Dhyd + 2*thick_t)))
            "Number of tubes in the bundle";
           Real N_t_p_p(start = qdot/(DTML*U_guess)/(pi*l*(Dhyd + 2*thick_t)/N_passes))
            "Number of tubes per pass";
           Real bs_f "Baffle spacing in percent";
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
          replaceable Heat_transfer.Tubes.Sieder_Tate hT_tube(Dhyd=Dhyd,
          eta_wall=bundle.state[1].eta*ones(Ncell))
            constrainedby Heat_transfer.Tubes.Base_classes.base_ht(Medium = Medium_t,
            Ncell = Ncell, state = bundle.state, mdot = m_t/N_t_p_p,
            Aflow = bundle.Aflow) annotation(choicesAllMatching = true);

          replaceable Pressure_drops.Tubes.Frank dp_tube(Dhyd = Dhyd, l = l/N_baffles_d,
          heads = 2.5/Ncell)
            constrainedby Pressure_drops.Tubes.Base_classes.base_dp(Medium = Medium_t,
            Ncell = Ncell, state = bundle.state, mdot = m_t/N_t_p_p,
            Aflow = bundle.Aflow) annotation(choicesAllMatching = true);

          replaceable Heat_transfer.Shell.single_phase_Kern hT_shell(
          Dhyd_o = Dhyd + 2*thick_t, layout = layout, pitch_f = pitch_f)
            constrainedby Heat_transfer.Shell.Base_classes.base_ht(Medium = Medium_s,
            Ncell = Ncell, state = shell.state, mdot = m_s, Aflow = shell.Aflow) annotation(choicesAllMatching = true);

          replaceable Pressure_drops.Shell.single_phase_Johnston dp_shell(
            l = l/N_baffles_d/N_passes,
            d_s = d_s,
            Dhyd = hT_shell.d_s_eq,
            l_b = l_b,
            eta_wall = shell.state[1].eta*ones(Ncell)) constrainedby
            Pressure_drops.Shell.Base_classes.base_dp(
            Medium = Medium_s,
            Ncell = Ncell,
            state = shell.state,
            mdot = m_s,
            Aflow = shell.Aflow) annotation (choicesAllMatching=true);

          //Defining the model for the bundle clearance
          replaceable function bundle_clearance =
              Miscellanea.Shell_clearance.base_clearance
                                                annotation(choicesAllMatching = true);

          //Defining the model for the cost
          replaceable function cost =
              Miscellanea.Cost.base_cost        annotation(choicesAllMatching = true);

          //Definiing the tubes and the shell
          Objects.tube_bundle
                      bundle(redeclare package Medium = Medium_t,
                      h_in = h_t_in,
                      h_out = h_t_out,
                      Ncell = Ncell, Aflow = 0.25*pi*Dhyd^2, mdot = m_t,
                      Dhyd = Dhyd, thick = thick_t, lambda = Material_t.lambda,
                      rho = Material_t.rho, N_tubes = N_tubes,
                      mdot_pt = m_t/N_t_p_p,N_passes = N_passes, layout = layout,
                      pin = pin_t, pitch_f = pitch_f);

          Objects.tube shell(
            redeclare package Medium = Medium_s,
            h_in=h_s_in,
            h_out=h_s_out,
            Ncell=Ncell,
            Aflow=(1 - 1/pitch_f)*d_s*l_b,
            mdot=m_s,
            mdot_pt=m_s,
            Dhyd=d_s,
            thick=thick_s,
            lambda=Material_s.lambda,
            rho=Material_s.rho,
            pin=pin_s);

          Miscellanea.check_velocity check_shell(redeclare package Medium = Medium_s,
                      T = shell.state[1].T, umin = min(hT_shell.u),
                      umax = max(hT_shell.u), geometry = "shell", op_p = p_s_in);

          Miscellanea.check_velocity check_tube(redeclare package Medium = Medium_t,
                      T = bundle.state[1].T, umin = min(hT_tube.u),
                      umax = max(hT_tube.u), geometry = "tube", op_p = p_t_in);

        protected
           parameter Modelica.SIunits.CoefficientOfHeatTransfer ht_t_f[Ncell]=
            ones(Ncell)*ht_t_f1
            "Tube fouling heat transfer coefficient (array)";
           parameter Modelica.SIunits.CoefficientOfHeatTransfer ht_s_f[Ncell]=
           ones(Ncell)*ht_s_f1
            "Shell fouling heat transfer coefficient (array)";
           Integer k[Ncell] "Index for the shell to tube topology";
           Modelica.SIunits.ThermalConductance kA_tot[Ncell]
            "Overall heat transfer coefficient (array)";
           Modelica.SIunits.ThermalConductance htA_shell[Ncell]
            "Shell heat transfer coefficient (array)";
           Modelica.SIunits.ThermalConductance htA_tubes[Ncell]
            "Tube heat transfer coefficient (array)";
           Modelica.SIunits.ThermalConductance G_wall[Ncell]
            "Wall heat transfer coefficient (array)";

        equation
            //Set boundary conditions at the inlet and outelt
            bundle.p_in        = p_t_in;
            shell.p_in         = p_s_in;
            bundle.h[1]        = h_t_in;
            shell.h[1]         = h_s_in;
            shell.h[Ncell + 1] = h_s_out;

            //Things get tricky here. We need to get the index of the shell cell seen by the tube cell
            for j in 1:N_passes loop
                 for i in 1+(j-1)*N_baffles_d:2:j*N_baffles_d loop
                    if mod(j,2) <> 0 then
                      if mod(N_baffles_d,2) <> 0 then
                        k[i]   = Ncell - N_passes*(i-(1+(j-1)*N_baffles_d)) - (j-1);
                      else
                        k[i]   = Ncell - N_passes*(i-(1+(j-1)*N_baffles_d)) - (N_passes
                        - j);
                      end if;
                    else
                        k[i]   = (N_passes + 1 - j) + N_passes*(i-(1+(j-1)*N_baffles_d));
                    end if;
                 end for;
             end for;

             for j in 1:N_passes loop
                 for i in 2+(j-1)*N_baffles_d:2:j*N_baffles_d loop
                    if mod(j,2) <> 0 then
                      if mod(N_baffles_d,2) <> 0 then
                        k[i]   = Ncell - (2*N_passes - 1) + (j-1) - N_passes*(i-(2+(j-1)
                        *N_baffles_d));
                      else
                        k[i]   = Ncell - N_passes - (j-1) - N_passes*(i-(2+(j-1)
                        *N_baffles_d));
                      end if;
                    else
                        k[i]   = N_passes + j + N_passes*(i-(2+(j-1)*N_baffles_d));
                    end if;
                 end for;
             end for;

            //A sperate loop for the heat transfer and energy balances
            for i in 1:Ncell loop
              G_wall[i]          = bundle.At*bundle.G_wall/bundle.Dhyd;
              htA_tubes[i]       = bundle.At/(1/hT_tube.ht[i] + 1/ht_t_f[i]);
              htA_shell[k[i]]    = shell.At/(1/hT_shell.ht[k[i]] + 1/ht_s_f[k[i]]);
              kA_tot[i]          = 1/(1/htA_tubes[i] + 1/G_wall[i] + 1/htA_shell[k[i]]);
              bundle.qdot[i]     = shell.qdot[k[i]];
              bundle.qdot[i]     = kA_tot[i]*(pin_s*shell.state[k[i]].T + pin_t*
              bundle.state[i].T);

              //Fluid weight calculation
              bundle.W_fluids[i] = 0.25*pi*Dhyd^2*l/N_baffles_d/N_passes*N_tubes*
              bundle.state[i].d;
              shell.W_fluids[i]  = 0.25*pi*(shell.Dhyd^2 - N_tubes*bundle.Dhyd_o^2)*l/
              N_baffles_d/N_passes*shell.state[i].d;

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
            qdot           = DTML_tilde*U*A;

            //Weight calculation
            bundle.W_dry = 0.25*pi*(bundle.Dhyd_o^2 - bundle.Dhyd^2)*l*N_tubes
                             *bundle.rho;
            shell.W_dry  = 0.25*pi*(shell.Dhyd_o^2 - shell.Dhyd^2)*l*shell.rho;
            W_dry          = 1.25*bundle.W_dry + 1.3*shell.W_dry;
            W_fluids       = sum(bundle.W_fluids) + sum(shell.W_fluids);
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

        model Flat_plate "Flat plate heat exchanger"
          import VIP;
          extends VIP.Design.Icons.flat_plate;
          replaceable package Medium_hot = VIP.Media.OneRandomOrganicFluid
                                                                       constrainedby
            Modelica.Media.Interfaces.PartialMedium "Medium model hot cells" annotation(choicesAllMatching = true);
          replaceable package Medium_cold = VIP.Media.OneRandomOrganicFluid
                                                                        constrainedby
            Modelica.Media.Interfaces.PartialMedium "Medium model cold cells" annotation(choicesAllMatching = true);
          replaceable Materials.void_material material
            "Material model for the plate"
            annotation(choicesAllMatching = true);
          final parameter Integer N_ch = tpg_hot.N_passes*N_ch_p
            "Total number of channels";
          parameter Integer N_ch_p = 3 "Number of channels per pass";
          parameter Integer N_plates = 2*N_ch + 1 "Number of plates";
          parameter Integer N_cell_pc = 3 "Number of cells per channel";
          parameter Modelica.SIunits.Length thick "Plate thickness";
          parameter Modelica.SIunits.Length b "Flow thickness";
          parameter Modelica.SIunits.Length w "Width";
          parameter Real X "Wave number";
          parameter Modelica.SIunits.Angle beta "Plate inclination angle";
          parameter Modelica.SIunits.Length d_pt "Port diameter";
          parameter Modelica.SIunits.SpecificEnthalpy h_hot_in
            "Inlet specific enthalpy hot cells";
          parameter Modelica.SIunits.SpecificEnthalpy h_hot_out
            "Outlet specific enthalpy hot cells";
          parameter Modelica.SIunits.AbsolutePressure p_hot_in
            "Inlet pressure hot cells";
          parameter Modelica.SIunits.AbsolutePressure p_hot_out
            "Outlet pressure hot cells";
          parameter Modelica.SIunits.SpecificEnthalpy h_cold_in
            "Inlet specific enthalpy cold cells";
          parameter Modelica.SIunits.SpecificEnthalpy h_cold_out
            "Outlet specific enthalpy cold cells";
          parameter Modelica.SIunits.AbsolutePressure p_cold_in
            "Inlet pressure cold cells";
          parameter Modelica.SIunits.AbsolutePressure p_cold_out
            "Outlet pressure cold cells";
          parameter Modelica.SIunits.Temperature  t_hot_in=
          Medium_hot.temperature_ph(p_hot_in, h_hot_in)
            "Inlet temperature hot cells";
          parameter Modelica.SIunits.Temperature  t_hot_out=
          Medium_hot.temperature_ph(p_hot_out, h_hot_out)
            "Outlet temperature hot cells";
          parameter Modelica.SIunits.Temperature  t_cold_in=
          Medium_cold.temperature_ph(p_cold_in, h_cold_in)
            "Inlet temperature cold cells";
          parameter Modelica.SIunits.Temperature  t_cold_out=
          Medium_cold.temperature_ph(p_cold_out, h_cold_out)
            "Outlet temperature cold cells";
          parameter Modelica.SIunits.MassFlowRate mdot_hot
            "Mass flow rate hot cells";
          parameter Modelica.SIunits.MassFlowRate mdot_cold
            "Mass flow rate cold cells";
          final parameter Modelica.SIunits.HeatFlowRate qdot = mdot_hot*(h_hot_in -
          h_hot_out) "Heat rate";
          parameter Modelica.SIunits.Temp_C DTML = Miscellanea.log_mean_delta_T(t_hot_in,
            t_hot_out, t_cold_in, t_cold_out)
            "Logarithmic mean temperature difference";
          parameter Modelica.SIunits.CoefficientOfHeatTransfer ht_hot_f1
            "Fouling heat transfer coefficient hot cells";
          parameter Modelica.SIunits.CoefficientOfHeatTransfer ht_cold_f1
            "Fouling heat transfer coefficient cold cells";
          parameter Modelica.SIunits.CoefficientOfHeatTransfer U_guess
            "Guess value for the global heat transfer coefficient";
          parameter Modelica.SIunits.Length l_start = qdot/(2*w*DTML*U_guess*N_plates) "Length";
          Modelica.SIunits.Length l(start = l_start, fixed = true) "Length";
          Modelica.SIunits.CoefficientOfHeatTransfer U
            "Global heat transfer coefficient";
          Modelica.SIunits.Area A "Heat transfer area";
          Modelica.SIunits.Area At( start = l_start*w/N_cell_pc)
            "Area of one cell";
          Modelica.SIunits.Temp_C  DTML_tilde
            "Logarithmic mean temperature difference corrected";
          Modelica.SIunits.Mass W_dry "Dry weight of the heat exchanger";
          Modelica.SIunits.Mass W_fluids "Weight of the fluids";
          Modelica.SIunits.Mass W_wet "Wet weight of the heat exchanger";
          Real PEC "Purchased equipment cost";
          Real NTU "Number of transport units";

          //Plate model
          Objects.plate plate(redeclare package Medium_hot = Medium_hot,
          redeclare package Medium_cold = Medium_cold,
          h_hot_in = h_hot_in, h_hot_out = h_hot_out, h_cold_in = h_cold_in,
          h_cold_out = h_cold_out, b = b, N_cell_pc = N_cell_pc, X = X,
          N_ch = N_ch, mdot_hot = tpg_hot.mdot_p, pin_hot = tpg_hot.pin,
          mdot_cold = tpg_cold.mdot_p, p_hot_in = p_hot_in, p_cold_in = p_cold_in,
          pin_cold = tpg_cold.pin, h_hot_start = tpg_hot.h_start,
          h_cold_start = tpg_cold.h_start);

          //Topologies
          replaceable Miscellanea.topology_PHE.parallel tpg_hot
          constrainedby Miscellanea.topology_PHE.Base_classes.base_top(
          redeclare package Medium = Medium_hot,
            N_cell_pc = N_cell_pc,
            N_ch = N_ch,
            mdot = mdot_hot,
            h = plate.h_hot,
            h_in = h_hot_in,
            h_out = h_hot_out,
            p_in = p_hot_in,
            boundary=false,
            p_out = p_hot_out) annotation(choicesAllMatching = true);

          //Heat transfer model for the hot side
          replaceable Heat_transfer.Plates.single_phase_Martin ht_hot(beta=beta)
          constrainedby Heat_transfer.Plates.Base_classes.base_ht(
          redeclare package Medium = Medium_hot,
          N_cell_pc=N_cell_pc,
          N_ch=N_ch,
          mdot=tpg_hot.mdot_p,
          state=plate.state_hot,
          Aflow=w*b,
          Dhyd=plate.Dhyd) annotation (choicesAllMatching=true);

          //Pressure drop model for the hot side
          replaceable Pressure_drops.Plates.single_phase_Martin dp_hot(beta=beta)
            constrainedby Pressure_drops.Plates.Base_classes.base_dp(
            redeclare package Medium = Medium_hot,
            N_cell_pc=N_cell_pc,
            N_ch=N_ch,
            N_ch_p=N_ch_p,
            N_passes=tpg_hot.N_passes,
            state_p=tpg_hot.state_p,
            mdot=tpg_hot.mdot_p,
            state=plate.state_hot,
            Aflow=w*b,
            Dhyd=plate.Dhyd,
            l=l/N_cell_pc,
            A_pt=0.25*Modelica.Constants.pi
                        *d_pt^2) annotation (choicesAllMatching=true);

          replaceable Miscellanea.topology_PHE.parallel tpg_cold
          constrainedby Miscellanea.topology_PHE.Base_classes.base_top(
          redeclare package Medium = Medium_cold,
            N_cell_pc = N_cell_pc,
            N_ch = N_ch,
            mdot = mdot_cold,
            h = plate.h_cold,
            h_in = h_cold_in,
            h_out = h_cold_out,
            p_in = p_cold_in,
            boundary=true,
            p_out = p_cold_out) annotation(choicesAllMatching = true);

          //Heat transfer model for the cold side
          replaceable Heat_transfer.Plates.single_phase_Martin ht_cold(beta=beta)
          constrainedby Heat_transfer.Plates.Base_classes.base_ht(
          redeclare package Medium = Medium_cold,
          N_cell_pc=N_cell_pc,
          N_ch=N_ch,
          mdot=tpg_cold.mdot_p,
          state=plate.state_cold,
          Aflow=w*b,
          Dhyd=plate.Dhyd) annotation (choicesAllMatching=true);

          //Pressure drop model for the cold side
          replaceable Pressure_drops.Plates.single_phase_Martin dp_cold(beta=beta)
            constrainedby Pressure_drops.Plates.Base_classes.base_dp(
            redeclare package Medium = Medium_cold,
            N_cell_pc=N_cell_pc,
            N_ch=N_ch,
            N_ch_p=N_ch_p,
            N_passes=tpg_cold.N_passes,
            state_p=tpg_cold.state_p,
            mdot=tpg_cold.mdot_p,
            state=plate.state_cold,
            Aflow=w*b,
            Dhyd=plate.Dhyd,
            l=l/N_cell_pc,
            A_pt=0.25*Modelica.Constants.pi
                        *d_pt^2) annotation (choicesAllMatching=true);

          //Defining the model for the cost
          replaceable function cost = Miscellanea.Cost.base_cost annotation(choicesAllMatching = true);

          replaceable Miscellanea.check_velocity check_hot(redeclare package
              Medium =
          Medium_hot, T = t_hot_in, umin = min(ht_hot.u), umax = max(ht_hot.u),
          geometry = "tube", op_p = p_hot_in);

          replaceable Miscellanea.check_velocity check_cold(redeclare package
              Medium =
          Medium_cold, T = t_cold_in, umin = min(ht_cold.u), umax = max(ht_cold.u),
          geometry = "tube", op_p = p_cold_in);

        protected
          parameter Modelica.SIunits.CoefficientOfHeatTransfer ht_hot_f[N_ch, N_cell_pc]=
            fill(ht_hot_f1, N_ch, N_cell_pc)
            "Fouling heat transfer coefficient (array) hot side";
          parameter Modelica.SIunits.CoefficientOfHeatTransfer ht_cold_f[N_ch, N_cell_pc]=
            fill(ht_cold_f1, N_ch, N_cell_pc)
            "Fouling heat transfer coefficient (array) cold side";
           Integer k1[N_ch] "Index for the flat plate topology columns";
           Integer k2[N_cell_pc] "Index for the flat plate topology rows";
           Modelica.SIunits.CoefficientOfHeatTransfer h_wall
            "Thermal conductance of each wall cell";
           Modelica.SIunits.CoefficientOfHeatTransfer G_tot[2, N_ch, N_cell_pc]
            "Global thermal conductance of each cell";
        equation

              for j in 1:N_ch loop
                k1[j] = N_ch + 1 - j;
                for i in 1:N_cell_pc loop

                  //Heat transfer coefficient
                  G_tot[1, j, i]           = 1/(1/ht_hot.ht[j, i] + 1/ht_hot_f[j, i] +
                  1/ht_cold.ht[k1[j], k2[i]] + 1/ht_cold_f[k1[j], k2[i]] + 1/h_wall);

                  //The heat rate of the metal wall follows the hot side convention
                  plate.qdot_wall[1, j, i] = At*G_tot[1, j, i]*(plate.state_hot[j, i].T
                  - plate.state_cold[k1[j], k2[i]].T);
                  plate.qdot_hot[j, i]     = plate.qdot_wall[1, j, i] +
                  plate.qdot_wall[2, j, i];

                  //Weight of the fluids
                  plate.W_fluids[j, i]     = At*b*(plate.state_hot[j, i].d +
                  plate.state_cold[j, i].d);

                  //First row hot side and last row cold side have only one heat flux
                  if (j == 1) then
                    k2[i]                         = N_cell_pc + 1 - i;
                    plate.qdot_wall[1, j, i]      = plate.qdot_cold[k1[j], k2[i]];
                    plate.qdot_wall[2, j, i]      = At*G_tot[2, j, i]* (
                    plate.state_hot[j, i].T - plate.state_cold[k1[j] - 1, k2[i]].T);
                    G_tot[2, j, i]                = 1/(1/ht_hot.ht[j, i] +
                    1/ht_hot_f[j, i] + 1/ht_cold.ht[k1[j] - 1, k2[i]] +
                    1/ht_cold_f[k1[j] - 1, k2[i]] + 1/h_wall);
                    //Check second principle of Thermodynamics
                    assert(plate.state_hot[j, i].T > plate.state_cold[k1[j] - 1,
                    k2[i]].T, "Second principle of Thermodynamics not respected",
                    AssertionLevel.warning);
                  elseif (j == N_ch) then
                    G_tot[2, j, i]                = 0;
                    plate.qdot_wall[2, j, i]      = 0;
                    plate.qdot_cold[k1[j], k2[i]] = plate.qdot_wall[1, j, i] +
                    plate.qdot_wall[2, j - 1, i];
                  else
                    plate.qdot_wall[2, j, i]      =  At*G_tot[2, j, i]*(
                    plate.state_hot[j, i].T - plate.state_cold[k1[j] - 1, k2[i]].T);
                    plate.qdot_cold[k1[j], k2[i]] = plate.qdot_wall[1, j, i] +
                    plate.qdot_wall[2, j - 1, i];
                    G_tot[2, j, i]                = 1/(1/ht_hot.ht[j, i] +
                    1/ht_hot_f[j, i] + 1/ht_cold.ht[k1[j], k2[i]] +
                    1/ht_cold_f[k1[j], k2[i]] + 1/h_wall);
                    //Check second principle of Thermodynamics
                    assert(plate.state_hot[j, i].T > plate.state_cold[k1[j] - 1,
                    k2[i]].T, "Second principle of Thermodynamics not respected",
                    AssertionLevel.warning);
                  end if;

                  //Check second principle of Thermodynamics
                   assert(plate.state_hot[j, i].T > plate.state_cold[k1[j], k2[i]].T,
                   "Second principle of Thermodynamics not respected",
                   AssertionLevel.warning);
                end for;
              end for;

              //Area and thermal conductance of the metal wall
              At             = l*w/N_cell_pc;
              h_wall         = material.lambda/thick;

              //Weight calculation
              plate.W_dry    = l*w*thick*material.rho;
              W_dry          = N_plates*plate.W_dry;
              W_fluids       = sum(plate.W_fluids);
              W_wet          = W_dry + W_fluids;

              //Area, global heat transfer coefficient and corrected DMTL
              if (sum(plate.state_cold[:,:].cp)*mdot_cold > sum(plate.state_hot[:,:].cp)*mdot_hot) then
                NTU            = U*A/(sum(plate.state_hot[:,:].cp)*mdot_hot)*(N_ch*N_cell_pc);
              else
                NTU            = U*A/(sum(plate.state_cold[:,:].cp)*mdot_cold)*(N_ch*N_cell_pc);
              end if;
              A              = 2*l*w*N_ch;
              U              = At*sum(G_tot)/A;
              qdot           = DTML_tilde*U*A;

              //Cost calculation
              PEC            = cost(A);

          annotation (experiment(Tolerance=1e-006), __Dymola_experimentSetupOutput);
        end Flat_plate;
      end HEX;
    end Components;

    package Materials
      "Package containing the properties of different materials"
      extends Modelica.Icons.MaterialPropertiesPackage;
      class S_AISI_1010 "S_AISI_1010"
          extends VIP.Design.Materials.void_material(
          materialName="S_AISI_1010",
          lambda=65.2,
          rho=7850);
      end S_AISI_1010;

      class SS_AISI_410 "SS_AISI_410"
          extends VIP.Design.Materials.void_material(
          materialName="SS_AISI_410",
          lambda=24.9,
          rho=7740);
      end SS_AISI_410;

      class SS_AISI_304 "SS_AISI_304"
          extends VIP.Design.Materials.void_material(
          materialName="SS_AISI_304",
          lambda=18.3866,
          rho=8030);
      end SS_AISI_304;

      class S_AISI_1040 "S_AISI_1040"
          extends VIP.Design.Materials.void_material(
          materialName="S_AISI_1040",
          lambda=50.7,
          rho=7845);
      end S_AISI_1040;

      class void_material "Void material"
          extends Modelica.Icons.MaterialProperty;
          parameter String materialName = "Random";
          parameter Modelica.SIunits.ThermalConductivity lambda = 50;
          parameter Modelica.SIunits.Density rho = 8000;
      end void_material;

      class T_R50250 "Titanium alloy R50250"
          extends VIP.Design.Materials.void_material(
          materialName="T_R50250",
          lambda=20.8,
          rho=4510);
      end T_R50250;
    end Materials;

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
        parameter Modelica.SIunits.Temperature T_sat = Medium.saturationTemperature(op_p)
          "Saturation temperature";
        input Modelica.SIunits.Temperature T "Fluid temperature";

      equation
        if (T < T_sat) then
          if (geometry == "tube") then
            assert(umax < 4*2, "Velocity twice higher than 4 m/s", AssertionLevel.warning);
            assert(umin > 1*0.5,"Velocity twice lower than 0.8 m/s", AssertionLevel.warning);
          elseif (geometry == "shell") then
            assert(umax < 1.2*2, "Velocity twice higher than 1.2 m/s", AssertionLevel.warning);
            assert(umin > 0.3*0.5,"Velocity twice lower than 0.3 m/s", AssertionLevel.warning);
          end if;
        elseif (T > T_sat) then
          if (op_p > 1e6) then
            assert(umax < 70*2, "Velocity twice higher than 70 m/s", AssertionLevel.warning);
            assert(umin > 50*0.5,"Velocity twice lower than 50 m/s", AssertionLevel.warning);
          elseif (op_p < 1e4) then
            assert(umax < 10*2, "Velocity twice higher than 10 m/s", AssertionLevel.warning);
            assert(umin > 5*0.5,"Velocity twice lower than 5 m/s", AssertionLevel.warning);
          else
            assert(umax < 30*2, "Velocity twice higher than 30 m/s", AssertionLevel.warning);
            assert(umin > 10*0.5,"Velocity twice lower than 10 m/s", AssertionLevel.warning);
          end if;
        end if;
      end check_velocity;

      package Cost "Functions to evaluate the cost of heat exchangers "
        extends VIP.Design.Icons.cost;
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

        function FP_Rafferty
          "Stainless steel flat plate heat exchanger Rafferty et al."
           extends base_cost(a = {1378.1781, 10.8962, 1});
        end FP_Rafferty;
        annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
                  -100},{100,100}}), graphics), Icon(graphics));
      end Cost;

      package Shell_clearance
        "Functions to evaluate the shell clearance based on the diameter of the tube bundle "
        extends VIP.Design.Icons.clearance;
        function Fixed_Utube
          "Fixed_U_tube shape clearance shell-bundle diameter"
           extends base_clearance(a = {10, 10, 0.2});
        end Fixed_Utube;

        function OPH "Outside packed head clearance shell-bundle diameter"
           extends base_clearance(a = {38, 0, 0.2});
        end OPH;

        function SRFH
          "Split ring floating head clearance shell-bundle diameter"
           extends base_clearance(a = {50, 28, 0.2});
        end SRFH;

        function PTFH
          "Pull through floating head clearance shell-bundle diameter"
           extends base_clearance(a = {88, 11, 0.2});
        // algorithm
        //   clearance :=1e-3*(88 + 11*(d_b - 0.2));
        end PTFH;

        function base_clearance
          "Base class for clearance shell-bundle diameter"
          input Modelica.SIunits.Length d_b "bundle diameter";
          output Modelica.SIunits.Length clearance
            "Shell to bundle diameter clearance";
        protected
          parameter Real a[3] = {10, 0, 0};
        algorithm
            clearance :=1e-3*(a[1] + a[2]*(d_b - a[3]));
        end base_clearance;
      end Shell_clearance;

      package numbers "Package containing adimensional numbers"
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
      end numbers;

      package topology_PHE "Topologies for the flat plate heat exchanger"
        class series "Series"
          extends Base_classes.base_top;
          final parameter Integer N_passes = 1 "Number of passes";
          parameter Modelica.SIunits.MassFlowRate mdot_p = mdot
            "Mass flow rate plate channel";
          parameter Modelica.SIunits.SpecificEnthalpy h_start[N_ch, N_cell_pc + 1]=
          fill(linspace(h_in, h_out, N_cell_pc + 1), N_ch)
            "Start value for the hot stream enthalpy matrix";
          Real pin[N_ch, N_cell_pc] "Pin for the heat flow";
          Medium.ThermodynamicState state_p[1]
            "Thermodynamic states at the inlet port";

        equation
          for j in 1:N_ch loop
            for i in 1:N_cell_pc loop
              if (mod(j,2) <> 0) then
                pin[j, i]  = N_ch;
              else
                pin[j, i]  = -N_ch;
              end if;
              //Boundary conditions
              if (i == 1) then
                if (j == 1) then
                  h_in                 = h[1, 1];
                else
                  if (mod(j,2) <> 0) then
                    h[j, 1]              = h[j - 1, 1];
                  else
                    h[j, N_cell_pc + 1]  = h[j - 1, N_cell_pc + 1];
                  end if;
                end if;
              end if;
            end for;
          end for;

          if boundary then
            h_out      = h[N_ch, N_cell_pc + 1];
          end if;

          state_p[1] = Medium.setState_ph(p_in, h_in);

        end series;

        class parallel "Parallel"
          extends Base_classes.base_top;
          final parameter Integer N_passes = 1 "Number of passes";
          parameter Modelica.SIunits.MassFlowRate mdot_p = mdot/N_ch
            "Mass flow rate plate channel";
          parameter Modelica.SIunits.SpecificEnthalpy h_start[N_ch, N_cell_pc + 1]=
          fill(linspace(h_in, h_out, N_cell_pc + 1), N_ch)
            "Start value for the hot stream enthalpy matrix";
          Real pin[N_ch, N_cell_pc] "Pin for the heat flow";
          Medium.ThermodynamicState state_p[1]
            "Thermodynamic states at the inlet port";
        equation
          for j in 1:N_ch loop
            for i in 1:N_cell_pc loop
              pin[j, i]  = 1;
              //Boundary conditions
              if (i == 1) then
                h[j, i]       = h_in;
              end if;
            end for;
          end for;

          if boundary then
            h_out      = sum(h[:, N_cell_pc + 1])/N_ch;
          end if;

          state_p[1] = Medium.setState_ph(p_in, h_in);

        end parallel;

        class two_pass_two_pass "Two pass two pass configuration"
          extends Base_classes.base_top;
          parameter Modelica.SIunits.MassFlowRate mdot_p = mdot/N_ch_p
            "Mass flow rate plate channel";
          final parameter Integer N_passes = 2 "Number of passes";
          parameter Integer N_ch_p( start = 3) "Number of channels per pass";
          parameter Integer stype "Stream type";
          parameter Modelica.SIunits.SpecificEnthalpy
          h_start[N_ch, N_cell_pc + 1] = (if stype == 1 then
          [fill(linspace(h_out, 0.5*(h_out + h_in), N_cell_pc + 1), N_ch_p);
          fill(linspace(h_in, 0.5*(h_out + h_in), N_cell_pc + 1), N_ch_p)] else
          [fill(linspace(0.5*(h_out + h_in), h_out, N_cell_pc + 1), N_ch_p);
          fill(linspace(0.5*(h_out + h_in), h_in, N_cell_pc + 1), N_ch_p)])
            "Start value for the hot stream enthalpy matrix";
          Real pin[N_ch, N_cell_pc] "Pin for the heat flow";
          Medium.ThermodynamicState state_p[2]
            "Thermodynamic states at the ports";
        equation
          if (stype == 1) then
            for j in 1:N_ch loop
              for i in 1:N_cell_pc loop
                if (j <= N_ch_p) then
                  pin[j, i]             = -1;
                  if (i == N_cell_pc) then
                    h[j, i + 1] = sum(h[N_ch_p + 1:N_ch, N_cell_pc + 1])/N_ch_p;
                  end if;
                else
                  pin[j, i]             = 1;
                  if (i == 1) then
                    h[j, i]             = h_in;
                  end if;
                end if;
              end for;
            end for;
            if boundary then
              h_out = sum(h[1:N_ch_p, 1])/N_ch_p;
            end if;
            state_p[1] = Medium.setState_ph(p_in, sum(h[N_ch_p + 1:N_ch, N_cell_pc + 1]
            /N_ch_p));
            state_p[2] = Medium.setState_ph(p_in, h_in);
          else
            for j in 1:N_ch loop
              for i in 1:N_cell_pc loop
                if (j <= N_ch_p) then
                  pin[j, i]             = 1;
                  if (i == 1) then
                    h[j, i] = sum(h[N_ch_p + 1:N_ch, 1])/N_ch_p;
                  end if;
                else
                  pin[j, i]             = -1;
                  if (i == N_cell_pc) then
                    h[j, i + 1]         = h_in;
                  end if;
                end if;
              end for;
            end for;
            if boundary then
              h_out = sum(h[1:N_ch_p, N_cell_pc + 1])/N_ch_p;
            end if;
            state_p[1] = Medium.setState_ph(p_in, sum(h[N_ch_p + 1:N_ch, 1])/N_ch_p);
            state_p[2] = Medium.setState_ph(p_in, h_in);
          end if;

        end two_pass_two_pass;

        package Base_classes "Base classes for the topology of the PHE"
          class base_top "Base class topology PHE"
            import VIP;
            replaceable package Medium = VIP.Media.OneRandomOrganicFluid
              "Medium model";
            parameter Integer N_ch(start = 6) "Total number of channels";
            parameter Integer N_cell_pc = 3 "Number of cells per channel";
            parameter Modelica.SIunits.SpecificEnthalpy h_in
              "Inlet specific enthalpy";
            parameter Modelica.SIunits.SpecificEnthalpy h_out
              "Outlet specific enthalpy";
            parameter Modelica.SIunits.AbsolutePressure p_in "Inlet pressure";
            parameter Modelica.SIunits.AbsolutePressure p_out "Outlet pressure";
            parameter Boolean boundary "Boundary condition for the enthalpy";
            parameter Modelica.SIunits.MassFlowRate mdot "Mass flow rate";
            input Modelica.SIunits.SpecificEnthalpy h[N_ch, N_cell_pc + 1]
              "Stream enthalpy matrix";
          end base_top;
        end Base_classes;
      end topology_PHE;
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
            Ellipse(
              origin={10,10},
              fillColor={76,76,76},
              pattern=LinePattern.None,
              fillPattern=FillPattern.Solid,
              extent={{-80.0,-80.0},{-20.0,-20.0}}),
            Ellipse(
              origin={10,10},
              pattern=LinePattern.None,
              fillPattern=FillPattern.Solid,
              extent={{0.0,-80.0},{60.0,-20.0}}),
            Ellipse(
              origin={10,10},
              fillColor={128,128,128},
              pattern=LinePattern.None,
              fillPattern=FillPattern.Solid,
              extent={{0.0,0.0},{60.0,60.0}}),
            Ellipse(
              origin={10,10},
              lineColor={128,128,128},
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid,
              extent={{-80.0,0.0},{-20.0,60.0}})}));
    end Miscellanea;

    package Objects "Package containing all the objects of the VIP"

      class tube "I am a tube and I contain all the my relevant informations"
          replaceable package Medium = VIP.Media.OneRandomOrganicFluid
          "Medium model";
          parameter Integer Ncell(start=3) "Number of cell elements";
          parameter Modelica.SIunits.Length thick "Thickness";
          parameter Modelica.SIunits.ThermalConductivity lambda
          "Thermal conductivity of the wall";
          parameter Modelica.SIunits.Density rho "Density of the wall";
          parameter Modelica.SIunits.SpecificEnthalpy h_in
          "Inlet specific enthalpy";
          parameter Modelica.SIunits.SpecificEnthalpy h_out
          "Outlet specific enthalpy";
          parameter Modelica.SIunits.MassFlowRate mdot "Mass flow rate";
          parameter Real pin "This pin identifies if the fluid is hot or cold";
          input Modelica.SIunits.Length Dhyd "Hydraulic diameter";
          input Modelica.SIunits.MassFlowRate mdot_pt "Mass flow rate per tube";
          input Modelica.SIunits.Area Aflow "Flow area";
          Modelica.SIunits.Area At "Heat transfer area";
          Modelica.SIunits.Length Dhyd_o "Outer hydraulic diameter";
          Modelica.SIunits.ThermalConductivity G_wall
          "Thermal conductivity of the wall";
          Modelica.SIunits.AbsolutePressure p_in "Inlet pressure";
          Medium.ThermodynamicState state[Ncell]
          "Thermodynamic states of the cells";
          Modelica.SIunits.SpecificEnthalpy h[Ncell + 1](start=linspace(h_in, h_out, Ncell + 1))
          "Specific enthalpies";
          Modelica.SIunits.HeatFlowRate qdot[Ncell] "Heat rate";
          Modelica.SIunits.Mass W_dry "Dry weight of one element";
          Modelica.SIunits.Mass W_fluids[Ncell] "Fluids weight of one element";

      equation
          //Trivial calculations
          Dhyd_o = Dhyd + 2*thick;
          G_wall = 2*lambda/log(Dhyd_o/Dhyd);

          for i in 1:Ncell loop
            qdot[i]           = pin*mdot*(h[i] - h[i + 1]);
            state[i]          = Medium.setState_ph(p_in, 0.5*(h[i] + h[i + 1]));
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
      end tube;

      class tube_bundle
        "I am a tube bundle and I contain all the my relevant information"
          extends VIP.Design.Objects.tube;
          parameter Integer N_passes "Number of tube passes";
          parameter Integer layout "Tube layout, 1 = triangular, 2 = squared";
          parameter Real pitch_f "Tube pitch";
          input Real N_tubes "Number of tubes in the bundle";
          Modelica.SIunits.Length d_b "Bundle diameter";
      equation
          d_b  =Design.Miscellanea.bundle_diameter(
              N_tubes,
              N_passes,
              Dhyd_o,
              layout) "calculate the tube bundle";
      end tube_bundle;

      class plate "I am a plate and I contain all the my relevant informations"
        import VIP;
          extends VIP.Design.Icons.plate;
          replaceable package Medium_hot = VIP.Media.OneRandomOrganicFluid
          "Medium model hot cells";
          replaceable package Medium_cold = VIP.Media.OneRandomOrganicFluid
          "Medium model cold cells";
          parameter Integer N_ch(start = 6) "Total number of channels";
          parameter Integer N_cell_pc(start = 3) "Number of cells per channels";
          parameter Modelica.SIunits.Length b "Plate thickness";
          parameter Real X "Wavenumber";
          parameter Modelica.SIunits.SpecificEnthalpy h_hot_in
          "Inlet specific enthalpy hot side";
          parameter Modelica.SIunits.SpecificEnthalpy h_hot_out
          "Outlet specific enthalpy hot side";
          parameter Modelica.SIunits.SpecificEnthalpy h_cold_in
          "Inlet specific enthalpy cold side";
          parameter Modelica.SIunits.SpecificEnthalpy h_cold_out
          "Outlet specific enthalpy cold side";
          final parameter Modelica.SIunits.Length Dhyd = 2*b/phi
          "Hydraulic diameter";
          parameter Modelica.SIunits.AbsolutePressure p_hot_in
          "Inlet pressure hot cells";
          parameter Modelica.SIunits.AbsolutePressure p_cold_in
          "Inlet pressure cold cells";
          parameter Modelica.SIunits.SpecificEnthalpy h_hot_start[N_ch, N_cell_pc + 1]
          "Hot stream temperature matrix";
          parameter Modelica.SIunits.SpecificEnthalpy h_cold_start[N_ch, N_cell_pc + 1]
          "Cold stream temperature matrix";
          input Modelica.SIunits.MassFlowRate mdot_hot
          "Mass flow rate hot cells";
          input Modelica.SIunits.MassFlowRate mdot_cold
          "Mass flow rate cold cells";
          input Real pin_hot[N_ch, N_cell_pc] "Pin hot cells";
          input Real pin_cold[N_ch, N_cell_pc] "Pin cold cells";
          Medium_hot.ThermodynamicState state_hot[N_ch, N_cell_pc]
          "Thermodynamic states hot cells";
          Medium_cold.ThermodynamicState state_cold[N_ch, N_cell_pc]
          "Thermodynamic states cold cells";
          Modelica.SIunits.HeatFlowRate qdot_hot[N_ch, N_cell_pc]
          "Heat rate hot cells";
          Modelica.SIunits.HeatFlowRate qdot_cold[N_ch, N_cell_pc]
          "Heat rate cold cells";
          Modelica.SIunits.HeatFlowRate qdot_wall[2, N_ch, N_cell_pc]
          "Heat rate metal wall";
          Modelica.SIunits.SpecificEnthalpy h_hot[N_ch, N_cell_pc + 1](
           start = h_hot_start) "Hot stream enthalpy matrix";
          Modelica.SIunits.SpecificEnthalpy h_cold[N_ch, N_cell_pc + 1](
           start = h_cold_start) "Cold stream enthalpy matrix";
          Modelica.SIunits.Mass W_dry "Dry weight of the plate";
          Modelica.SIunits.Mass W_fluids[N_ch, N_cell_pc] "Fluid weights";
      protected
          final parameter Real phi = (1 + sqrt(1 + X^2) + 4*sqrt(1 + 0.5*X^2))/6
          "Corrugation pitch";

      equation
            for j in 1:N_ch loop
              for i in 1:N_cell_pc loop
                qdot_hot[j, i]     = pin_hot[j, i]*mdot_hot*(h_hot[j, i] -
                h_hot[j, i + 1]);
                qdot_cold[j, i]    = pin_cold[j, i]*mdot_cold*(h_cold[j, i + 1] -
                h_cold[j, i]);
                state_hot[j, i]    = Medium_hot.setState_ph(p_hot_in,
                  0.5*(h_hot[j, i] + h_hot[j, i + 1]));
                state_cold[j, i]   = Medium_cold.setState_ph(p_cold_in,
                  0.5*(h_cold[j, i] + h_cold[j, i + 1]));
              end for;
          end for;
        annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},
                  {100,100}}), graphics));
      end plate;
    end Objects;

    package Pressure_drops "A package containing pressure drops correlations"
      package Tubes "heat transfer correlations in tubes"
          extends VIP.Design.Icons.tube;
        class Frank "Frank correlation for tubes"
            extends VIP.Design.Pressure_drops.Tubes.Base_classes.base_dp;
            output Modelica.SIunits.AbsolutePressure dp[Ncell]
            "Pressure drops cells";
            output Modelica.SIunits.AbsolutePressure dp_tot
            "Pressure drops tubes";
            input Modelica.SIunits.Length l "Lenght";
            input Modelica.SIunits.Length Dhyd "Hydraulic diameter";
            input Modelica.SIunits.AbsolutePressure heads
            "Number of velocity heads";
            Modelica.SIunits.ReynoldsNumber Re[Ncell] "Reynolds number";
            Modelica.SIunits.Velocity u[Ncell](start=ones(Ncell)) "Velocity";
            Real csi[Ncell] "Friction factor";

        equation
          for i in 1:Ncell loop
            u[i]     = mdot/state[i].d/Aflow;
            Re[i]    =Design.Miscellanea.numbers.Reynolds(
                    u[i],
                    state[i].d,
                    state[i].eta,
                    Dhyd);
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

        package Base_classes "Base classes for the pressure drops"
          class base_dp "Basic pressure drop correlation"
              replaceable package Medium = VIP.Media.OneRandomOrganicFluid
              "Medium model";
              parameter Integer Ncell(start=3) "Number of cell elements";
              input Modelica.SIunits.Area Aflow "Cross-sectional area";
              input Modelica.SIunits.MassFlowRate mdot "Mass flow rate";
              input Medium.ThermodynamicState state[Ncell];
          end base_dp;
        end Base_classes;
      end Tubes;

      package Shell "heat transfer correlations in shells"
        extends VIP.Design.Icons.shell;

        class single_phase_Kern
          "Kern correlation for shell single phase to single phase"
            extends VIP.Design.Pressure_drops.Shell.Base_classes.base_Kern;
            Modelica.SIunits.AbsolutePressure dp[Ncell] "Pressure drops cells";
            Modelica.SIunits.AbsolutePressure dp_tot "Pressure drops tubes";
            input Modelica.SIunits.DynamicViscosity eta_wall[Ncell]
            "Viscosity of the fluid at the wall";
            Modelica.SIunits.ReynoldsNumber Re[Ncell](start=10e6*ones(Ncell))
            "Reynolds number";
            Modelica.SIunits.Velocity u[Ncell](start=ones(Ncell)) "Velocity";
            Real csi[Ncell] "Friction factor";

        equation
          for i in 1:Ncell loop
            u[i]     = mdot/state[i].d/Aflow;
            Re[i]    =Design.Miscellanea.numbers.Reynolds(
                    u[i],
                    state[i].d,
                    state[i].eta,
                    Dhyd);
            if (Re[i] < 3e2) then
              csi[i] = 18.4810 * Re[i]^(-0.9344);
            else
              csi[i] = 0.2581 * Re[i]^(-0.1759);
            end if;
            dp[i]    =  4*csi[i]*d_s/Dhyd*l/l_b*state[i].d*u[i]^2;
          end for;

          dp_tot = sum(dp);

        end single_phase_Kern;

        class condensation_Kern "Kern correlation for shell condensation"
            extends VIP.Design.Pressure_drops.Shell.Base_classes.base_Kern;
            parameter Medium.ThermodynamicState state_l
            "Thermodynamic state in saturated liquid";
            parameter Medium.ThermodynamicState state_v
            "Thermodynamic state in saturated vapor";
            Modelica.SIunits.AbsolutePressure dp[Ncell] "Pressure drops cells";
            Modelica.SIunits.AbsolutePressure dp_tot "Pressure drops tubes";
            parameter Real X "Factor for pressure drop in condensation";
            input Modelica.SIunits.DynamicViscosity eta_wall[Ncell]
            "exponent for the viscosity correction";
            Modelica.SIunits.ReynoldsNumber Re[Ncell](start=10e6*ones(Ncell))
            "Reynolds number";
            Modelica.SIunits.Velocity u[Ncell](start=ones(Ncell)) "Velocity";
            Real csi[Ncell] "Friction factor";

        equation
          for i in 1:Ncell loop
            u[i]     = mdot/state[i].d/Aflow "tube velocity";
            Re[i]    =Design.Miscellanea.numbers.Reynolds(
                    u[i],
                    state[i].d,
                    state[i].eta,
                    Dhyd);
            if (Re[i] < 3e2) then
              csi[i] = 18.4810 * Re[i]^(-0.9344);
            else
              csi[i] = 0.2581 * Re[i]^(-0.1759);
            end if;
            if (state[i].h < state_l.h or state[i].h > state_v.h) then
              dp[i]  =  4*csi[i]*d_s/Dhyd*l/l_b*state[i].d*u[i]^2;
            else
              dp[i]  =  4*X*csi[i]*d_s/Dhyd*l/l_b*state[i].d*u[i]^2;
            end if;
          end for;

          dp_tot = sum(dp);

        end condensation_Kern;

        class single_phase_Johnston "Wills and Johnston correlation for shell"
            extends VIP.Design.Pressure_drops.Shell.Base_classes.base_Kern;
            output Modelica.SIunits.AbsolutePressure dp[Ncell]
            "Pressure drops tubes";
            output Modelica.SIunits.AbsolutePressure dp_tot
            "Pressure drops tubes";
            input Modelica.SIunits.DynamicViscosity eta_wall[Ncell]
            "exponent for the viscosity correction";
            Modelica.SIunits.ReynoldsNumber Re[Ncell](start=10e6*ones(Ncell))
            "Reynolds number tubes (average)";
            Modelica.SIunits.Velocity u[Ncell](start=ones(Ncell))
            "Velocity (average)";
            Real csi[Ncell] "Friction factor";

        equation
          for i in 1:Ncell loop
            u[i]   = mdot/state[i].d/Aflow;
            Re[i]  =Design.Miscellanea.numbers.Reynolds(
                    u[i],
                    state[i].d,
                    state[i].eta,
                    Dhyd);
            assert(Re[i] > 4e2, "Reynolds number is lower than 1e2 to use Wills_Johnston", AssertionLevel.warning);
            csi[i] = 1.7789*Re[i]^(-0.195868);
            dp[i]  = 0.5*csi[i]*d_s/Dhyd*l/l_b*(eta_wall[i]/state[i].eta)^0.14*state[i].d*u[i]^2;
          end for;

          dp_tot = sum(dp);

        end single_phase_Johnston;

        class single_phase_Bell_Delaware
          "Bell Delaware correlation for shell side for single phase to single phase"
            extends
            VIP.Design.Pressure_drops.Shell.Base_classes.base_Bell_Delaware;
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

        protected
             Real b1[Ncell];
             Real b2[Ncell];
             Real b[Ncell];
        equation

          for i in 1:Ncell loop

            u[i]    = mdot/state[i].d/S_m;
            u_w[i]  = mdot/state[i].d/S_w;
            Re[i]   =Design.Miscellanea.numbers.Reynolds(
                    u[i],
                    state[i].d,
                    state[i].eta,
                    Dhyd_o);

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
                b1[i]  = 45.1;
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

            csi[i] = b1[i]*(1.33/pitch_f)^b[i]*Re[i]^b2[i];

            if (i > Ncell - N_passes) then
              dp_w[i]  = 0;
            else
              dp_w[i]  = (2 + 0.6*N_tcw)*mdot^2/(2*S_w*S_m*state[i].d)*(N_baffles
              /(N_passes*N_baffles_d));
            end if;

            if (i <= N_passes or i > Ncell - N_passes) then
              dp_c[i] = 0;
              dp_e[i] = 2*(N_tcc + N_tcw)*csi[i]*(state[i].d*u[i]^2/2)*(N_baffles
              /(N_passes*N_baffles_d));
            else
              dp_e[i] = 0;
              dp_c[i] = ff*2*N_tcc*csi[i]*(state[i].d*u[i]^2/2)*(N_baffles
              /(N_passes*N_baffles_d));
            end if;

          end for;

          dp_w_tot = sum(dp_w)*R_l;
          dp_c_tot = sum(dp_c)*R_b*R_l;
          dp_e_tot = sum(dp_e)*R_b;

          dp_tot = dp_w_tot + dp_c_tot + dp_e_tot;

        end single_phase_Bell_Delaware;

        class condensation_Bell_Delaware
          "Bell Delaware correlation for shell side for condensation"
            extends
            VIP.Design.Pressure_drops.Shell.Base_classes.base_Bell_Delaware;
            parameter Medium.ThermodynamicState state_l
            "Thermodynamic state in saturated liquid";
            parameter Medium.ThermodynamicState state_v
            "Thermodynamic state in saturated vapor";
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
            Modelica.SIunits.ReynoldsNumber Re[Ncell] "Reynolds number";
            Modelica.SIunits.ReynoldsNumber Re_l(start=1e5)
            "Reynolds number if all liquid phase";
            Modelica.SIunits.Velocity u[Ncell](start=ones(Ncell))
            "Velocities inside the shell";
            Modelica.SIunits.Velocity u_w[Ncell](start=ones(Ncell))
            "Velocities inside the shell";
            Modelica.SIunits.Velocity u_l(start=1)
            "Velocity if all liquid phase";
            Real csi[Ncell] "Friction factor";
            Real csi_l "Friction factor if all liquid phase";
            Real xq[Ncell] "Vapor quality";
            Real phi_lo[Ncell] "Two phase multiplier cross flow";
            Real phi_lo_w[Ncell]
            "Two phase multiplier baffle window, horizontal cut";
            Real Y "Parameter cross flow";
            Real Y_w "Parameter baffle window, horizontal cut";
        protected
             Real b1[Ncell];
             Real b2[Ncell];
             Real b[Ncell];
             Real b1_l;
             Real b2_l;
             Real b_l;
        equation

          u_l       = mdot/state_l.d/S_m;
          Y         = sqrt(state_l.d/state_v.d)*(state_v.eta/state_l.eta)^0.23;
          Y_w       = sqrt(state_l.d/state_v.d);
          Re_l      =Design.Miscellanea.numbers.Reynolds(
                  u_l,
                  state_l.d,
                  state_l.eta,
                  Dhyd_o);

            if (Re_l <= 1e1) then
              if (layout == 1) then
                b1_l  = 48;
                b2_l  = -1;
              else
                b1_l  = 35;
                b2_l  = -1;
              end if;
            elseif (Re_l > 1e1 and Re_l <= 1e2) then
              if (layout == 1) then
                b1_l  = 45;
                b2_l  = -0.973;
              else
                b1_l  = 32.1;
                b2_l  = -0.963;
              end if;
            elseif (Re_l > 1e2 and Re_l <= 1e3) then
              if (layout == 1) then
                b1_l  = 4.57;
                b2_l  = -0.476;
              else
                b1_l  = 6.09;
                b2_l  = -0.602;
              end if;
            elseif (Re_l > 1e3 and Re_l < 1e4) then
              if (layout == 1) then
                b1_l  = 0.486;
                b2_l  = -0.152;
              else
                b1_l  = 0.0815;
                b2_l  = 0.022;
              end if;
            else
              if (layout == 1) then
                b1_l  = 0.372;
                b2_l  = -0.123;
              else
                b1_l  = 0.391;
                b2_l  = -0.148;
              end if;
            end if;

            if (layout == 1) then
              b_l   = 7/(1 + 0.14*Re_l^0.5);
            else
              b_l   = 6.3/(1 + 0.14*Re_l^0.378);
            end if;

            csi_l = b1_l*(1.33/pitch_f)^b_l*Re_l^b2_l;

          for i in 1:Ncell loop

            u[i]    = mdot/state[i].d/S_m;
            u_w[i]  = mdot/state[i].d/S_w;
            Re[i]   =Design.Miscellanea.numbers.Reynolds(
                    u[i],
                    state[i].d,
                    state[i].eta,
                    Dhyd_o);

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

            csi[i] = b1[i]*(1.33/pitch_f)^b[i]*Re[i]^b2[i];

            if (state[i].h < state_l.h or state[i].h > state_v.h) then

              if (i > Ncell - N_passes) then
                dp_w[i]    = 0;
              else
                dp_w[i]    = (2 + 0.6*N_tcw)*mdot^2/(2*S_w*S_m*state[i].d)*(N_baffles/
                (N_passes*N_baffles_d));
              end if;

              if (i <= N_passes or i > Ncell - N_passes) then
                dp_c[i]   = 0;
                dp_e[i]   = 2*(N_tcc + N_tcw)*csi[i]*(state[i].d*u[i]^2/2)*(N_baffles/
                (N_passes*N_baffles_d));
              else
                dp_e[i]   = 0;
                dp_c[i]   = ff*2*N_tcc*csi[i]*(state[i].d*u[i]^2/2)*(N_baffles/
                (N_passes*N_baffles_d));
              end if;

              xq[i]       = 1;
              phi_lo[i]   = 1;
              phi_lo_w[i] = 1;

            else

              xq[i]       = (state[i].h - state_l.h)/(state_v.h - state_l.h);
              phi_lo[i]   = 1 + (Y^2 - 1)*(0.75*(xq[i]*(1 - xq[i]))^0.77 + xq[i]^1.54);
              phi_lo_w[i] = 1 + (Y_w^2 - 1)*(2*xq[i]/(Y_w + 1));

              if (i > Ncell - N_passes) then
                dp_w[i]   = 0;
              else
                dp_w[i]   = phi_lo_w[i]*(2 + 0.6*N_tcw)*mdot^2/(2*S_w*S_m*state_l.d)*
                (N_baffles/(N_passes*N_baffles_d));
              end if;

              if (i <= N_passes or i > Ncell - N_passes) then
                dp_c[i]   = 0;
                dp_e[i]   = phi_lo[i]*(N_tcc + N_tcw)*csi_l*(state_l.d*u_l^2)*
                (N_baffles/(N_passes*N_baffles_d));
              else
                dp_e[i]   = 0;
                dp_c[i]   = phi_lo[i]*ff*N_tcc*csi_l*(state_l.d*u_l^2)*(N_baffles/
                (N_passes*N_baffles_d));
              end if;
            end if;
          end for;

          dp_w_tot = sum(dp_w)*R_l;
          dp_c_tot = sum(dp_c)*R_b*R_l;
          dp_e_tot = sum(dp_e)*R_b;

          dp_tot = dp_w_tot + dp_c_tot + dp_e_tot;

        end condensation_Bell_Delaware;

        package Base_classes "Base classes for the pressure drops"
          class base_dp "Basic pressure drop correlation"
              replaceable package Medium = VIP.Media.OneRandomOrganicFluid
              "Medium model";
              parameter Integer Ncell(start=3) "Number of cell elements";
              input Modelica.SIunits.Area Aflow "Cross-sectional area";
              input Modelica.SIunits.MassFlowRate mdot "Mass flow rate";
              input Medium.ThermodynamicState state[Ncell];
          end base_dp;

          class base_Bell_Delaware
            "Bell Delaware correlation for shell side base class"
              extends VIP.Design.Pressure_drops.Shell.Base_classes.base_dp;
              input Modelica.SIunits.Length Dhyd "Hydraulic diameter";
              input Modelica.SIunits.Length Dhyd_o "Outer hydraulic diameter";
              parameter Real b_cut "Baffle cut";
              parameter Real ttb = 0.8e-3 "tube to baffle clearance";
              parameter Real bts = 4e-3 "tube to baffle clearance";
              parameter Real pitch_f "Tube pitch";
              parameter Real N_ss
              "The number of sealing strips (pairs) in one baffle spacing";
              parameter Integer layout
              "Tube layout 1 = triangular, 2 = squared";
              parameter Integer N_passes "Number of tube passes";
              parameter Real ff "Fouling factor (=1 if clean)";
              parameter Integer N_baffles "Number of baffles";
              parameter Integer N_baffles_d
              "The number of discretization volumes is always N_baffles + 1";
              input Modelica.SIunits.Length d_s "Shell diameter";
              input Modelica.SIunits.Length d_b "Bundle diameter";
              input Real N_tubes "Number of tubes in the bundle";
              input Modelica.SIunits.Length l_b "Baffle lenght";
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
              Modelica.SIunits.Angle teta_ctl
              "The upper centriangle of baffle cut";

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

          end base_Bell_Delaware;

          class base_Kern "Kern correlation for shell base class"
              extends VIP.Design.Pressure_drops.Shell.Base_classes.base_dp;
              parameter Modelica.SIunits.Length l "Lenght";
              input Modelica.SIunits.Length Dhyd "Hydraulic diameter";
              input Modelica.SIunits.Length d_s "Shell diameter";
              input Modelica.SIunits.Length l_b "Baffle length";
          end base_Kern;
        end Base_classes;
      end Shell;

      package Plates
        extends Icons.plate;
        class single_phase_Martin
          "Martin pressure drop correlation for single phase"
            extends Base_classes.base_dp;
            parameter Modelica.SIunits.Angle beta "Plate inclination angle";
            Modelica.SIunits.ReynoldsNumber Re[N_ch, N_cell_pc]
            "Reynolds number";
            Modelica.SIunits.Velocity u[N_ch, N_cell_pc] "Velocity";
            Modelica.SIunits.Velocity u_pt[N_passes] "Port velocity";
            Modelica.SIunits.AbsolutePressure dp[N_ch, N_cell_pc]
            "Pressure drops cells";
            Modelica.SIunits.AbsolutePressure dp_pt[N_passes]
            "Pressure drops ports";
            Modelica.SIunits.AbsolutePressure dp_plates[N_ch]
            "Pressure drops plates";
            Modelica.SIunits.AbsolutePressure dp_tot "Total pressure drops";

        protected
            Real csi[N_ch, N_cell_pc] "Friction factor";
            Real csi0[N_ch, N_cell_pc] "Friction factor 0";
            Real csi1[N_ch, N_cell_pc] "Friction factor 1";
        equation
            for j in 1:N_ch loop
              dp_plates[j]   = sum(dp[j,:]);
              for i in 1:N_cell_pc loop
                u[j, i]      = mdot/state[j, i].d/Aflow;
                Re[j, i]     =Design.Miscellanea.numbers.Reynolds(
                      u[j, i],
                      state[j, i].d,
                      state[j, i].eta,
                      Dhyd);
                if (Re[j, i] < 2e3) then
                  csi0[j, i] = 16/Re[j, i];
                  csi1[j, i] = 149.25/Re[j, i] + 0.9625;
                else
                  csi0[j, i] = 1/(1.56*log(Re[j, i]) - 3)^2;
                  csi1[j, i] = 9.75/Re[j, i]^0.289;
                end if;
                csi[j, i]    = 1/((1 - cos(beta))/sqrt(3.8*csi1[j, i]) + cos(beta)/
                sqrt(0.045*tan(beta) + 0.09*sin(beta) + csi0[j, i]/cos(beta)))^2;
                dp[j, i]     = 2*csi[j, i]*(l/Dhyd)*state[j, i].d*u[j, i]^2;
              end for;
            end for;

            for j in 1:N_passes loop
              u_pt[j]     = N_ch_p*mdot/(state_p[j].d*A_pt);
              dp_pt[j]    = 0.5*1.3*state_p[j].d*u_pt[j]^2;
            end for;

            dp_tot        = sum(dp_plates)/N_ch_p + sum(dp_pt);

        end single_phase_Martin;

        class evaporation_Martin
          "Lockhart-Martinelli pressure drop correlation for evaporation heat transfer. Martin for single phase."
            extends Base_classes.base_dp;
            parameter Modelica.SIunits.Angle beta "Plate inclination angle";
            parameter Medium.ThermodynamicState state_l
            "Thermodynamic state in saturated liquid";
            parameter Medium.ThermodynamicState state_v
            "Thermodynamic state in saturated vapor";
            Modelica.SIunits.ReynoldsNumber Re[N_ch, N_cell_pc]
            "Reynolds number";
            Modelica.SIunits.ReynoldsNumber Re_l[N_ch, N_cell_pc]
            "Reynolds number liquid phase";
            Modelica.SIunits.ReynoldsNumber Re_v[N_ch, N_cell_pc]
            "Reynolds number vapor phase";
            Modelica.SIunits.Velocity u[N_ch, N_cell_pc] "Velocity";
            Modelica.SIunits.Velocity u_pt[N_passes] "Port velocity";
            Real xq[N_ch, N_cell_pc] "Vapor quality";
            Modelica.SIunits.AbsolutePressure dp[N_ch, N_cell_pc]
            "Pressure drops cells";
            Modelica.SIunits.AbsolutePressure dp_pt[N_passes]
            "Pressure drops ports";
            Modelica.SIunits.AbsolutePressure dp_plates[N_ch]
            "Pressure drops plates";
            Modelica.SIunits.AbsolutePressure dp_tot "Total pressure drops";

        protected
            Modelica.SIunits.AbsolutePressure dp_l[N_ch, N_cell_pc]
            "Pressure drops liquid phase";
            Modelica.SIunits.AbsolutePressure dp_v[N_ch, N_cell_pc]
            "Pressure drops vapor phase";
            Real X_tt[N_ch, N_cell_pc] "Martinelli parameter";
            Real phi_l[N_ch, N_cell_pc] "Multiplier liquid phase";
            Real phi_v[N_ch, N_cell_pc] "Multiplier vapor phase";
            Modelica.SIunits.Velocity u_l[N_ch, N_cell_pc]
            "Velocity liquid phase";
            Modelica.SIunits.Velocity u_v[N_ch, N_cell_pc]
            "Velocity vapor phase";
            Real csi[N_ch, N_cell_pc] "Friction factor";
            Real csi_l[N_ch, N_cell_pc] "Friction factor liquid phase";
            Real csi_v[N_ch, N_cell_pc] "Friction factor vapor phase";
            Real csi0[N_ch, N_cell_pc] "Friction factor 0";
            Real csi1[N_ch, N_cell_pc] "Friction factor 1";
            Real C[N_ch, N_cell_pc] "Constant for the two phase multiplier";
        equation

            for j in 1:N_ch loop
              dp_plates[j]   = sum(dp[j,:]);
              for i in 1:N_cell_pc loop
                u[j, i]      = mdot/state[j, i].d/Aflow;
                Re[j, i]     =Design.Miscellanea.numbers.Reynolds(
                      u[j, i],
                      state[j, i].d,
                      state[j, i].eta,
                      Dhyd);
                if (Re[j, i] < 2e3) then
                  csi0[j, i] = 16/Re[j, i];
                  csi1[j, i] = 149.25/Re[j, i] + 0.9625;
                else
                  csi0[j, i] = 1/(1.56*log(Re[j, i]) - 3)^2;
                  csi1[j, i] = 9.75/Re[j, i]^0.289;
                end if;
                csi[j, i]    = 1/((1 - cos(beta))/sqrt(3.8*csi1[j, i]) + cos(beta)/
                sqrt(0.045*tan(beta) + 0.09*sin(beta) + csi0[j, i]/cos(beta)))^2;
                if (state[j, i].h < state_l.h or state[j, i].h > state_v.h) then
                  xq[j, i]     = 0;
                  u_l[j, i]    = 0;
                  u_v[j, i]    = 0;
                  Re_l[j, i]   = 0;
                  Re_v[j, i]   = 0;
                  X_tt[j, i]   = 0;
                  C[j, i]      = 0;
                  phi_l[j, i]  = 0;
                  phi_v[j, i]  = 0;
                  dp_l[j, i]   = 0;
                  dp_v[j, i]   = 0;
                  csi_l[j, i]  = 0;
                  csi_v[j, i]  = 0;
                  dp[j, i]     = 2*csi[j, i]*(l/Dhyd)*state[j, i].d*u[j, i]^2;
                else
                  xq[j, i]     = (state[j, i].h - state_l.h)/(state_v.h - state_l.h);
                  u_l[j, i]    = mdot*(1 - xq[j, i])/state_l.d/Aflow;
                  u_v[j, i]    = mdot*xq[j, i]/state_v.d/Aflow;
                  Re_l[j, i]   =Design.Miscellanea.numbers.Reynolds(
                        u_l[j, i],
                        state_l.d,
                        state_l.eta,
                        Dhyd);
                  Re_v[j, i]   =Design.Miscellanea.numbers.Reynolds(
                        u_v[j, i],
                        state_v.d,
                        state_v.eta,
                        Dhyd);
                  if (Re_l[j, i] <= 1e3) then
                    csi_l[j, i] = 16/Re_l[j, i];
                  elseif (Re_l[j, i] > 1e3 and Re_l[j, i] <= 2e3) then
                    csi_l[j, i] = (1e-3*Re_l[j, i] - 1)*(16/Re_l[j, i]) + (2 - 1e-3*
                    Re_l[j, i])*4.6e-2/Re_l[j, i]^0.2;
                  else
                    csi_l[j, i] = 4.6e-2/Re_l[j, i]^0.2;
                  end if;
                  if (Re_v[j, i] <= 1e3) then
                    csi_v[j, i] = 16/Re_v[j, i];
                  elseif (Re_v[j, i] > 1e3 and Re_v[j, i] <= 2e3) then
                    csi_v[j, i] = (1e-3*Re_v[j, i] - 1)*(16/Re_v[j, i]) + (2 - 1e-3*
                    Re_v[j, i])*4.6e-2/Re_v[j, i]^0.2;
                  else
                    csi_v[j, i] = 4.6e-2/Re_v[j, i]^0.2;
                  end if;
                  dp_l[j, i]   = 2*csi_l[j, i]*(l/Dhyd)*state_l.d*u_l[j, i]^2;
                  dp_v[j, i]   = 2*csi_v[j, i]*(l/Dhyd)*state_v.d*u_v[j, i]^2;
                  X_tt[j, i]   = sqrt(dp_l[j, i]/dp_v[j, i]);
                  if (Re_l[j, i] > 1.5e3 and Re_v[j, i] > 1.5e3) then
                    C[j, i]    = 20;
                  elseif (Re_l[j, i] <= 1.5e3 and Re_v[j, i] > 1.5e3) then
                    C[j, i]    = 12;
                  elseif (Re_l[j, i] > 1.5e3 and Re_v[j, i] <= 1.5e3) then
                    C[j, i]    = 10;
                  else
                    C[j, i]    = 5;
                  end if;
                  phi_l[j, i]  = 1 + C[j, i]/X_tt[j, i] + 1/X_tt[j, i]^2;
                  phi_v[j, i]  = 1 + C[j, i]*X_tt[j, i] + X_tt[j, i]^2;
                  if (dp_l[j, i]*phi_l[j, i] > dp_v[j, i]*phi_v[j, i]) then
                    dp[j, i] = dp_v[j, i]*phi_v[j, i];
                  else
                    dp[j, i] = dp_l[j, i]*phi_l[j, i];
                  end if;
                end if;
              end for;
            end for;

            for j in 1:N_passes loop
              u_pt[j]     = N_ch_p*mdot/(state_p[j].d*A_pt);
              dp_pt[j]    = 0.5*1.3*state_p[j].d*u_pt[j]^2;
            end for;

            dp_tot        = sum(dp_plates)/N_ch_p + sum(dp_pt);

        end evaporation_Martin;

        class single_phase_Coulson
          "Coulson pressure drop correlation for single phase"
            extends Base_classes.base_dp;
            Modelica.SIunits.ReynoldsNumber Re[N_ch, N_cell_pc]
            "Reynolds number";
            Modelica.SIunits.Velocity u[N_ch, N_cell_pc] "Velocity";
            Modelica.SIunits.Velocity u_pt[N_passes] "Port velocity";
            Modelica.SIunits.AbsolutePressure dp[N_ch, N_cell_pc]
            "Pressure drops cells";
            Modelica.SIunits.AbsolutePressure dp_pt[N_passes]
            "Pressure drops ports";
            Modelica.SIunits.AbsolutePressure dp_plates[N_ch]
            "Pressure drops plates";
            Modelica.SIunits.AbsolutePressure dp_tot "Total pressure drops";
        protected
            Real csi[N_ch, N_cell_pc] "Friction factor";
        equation

            for j in 1:N_ch loop
              dp_plates[j]   = sum(dp[j,:]);
              for i in 1:N_cell_pc loop
                u[j, i]   = mdot/state[j, i].d/Aflow;
                Re[j, i]  =Design.Miscellanea.numbers.Reynolds(
                      u[j, i],
                      state[j, i].d,
                      state[j, i].eta,
                      Dhyd);
                csi[j, i] = 0.6/Re[j, i]^0.3;
                dp[j, i]  = 4*csi[j, i]*(l/Dhyd)*state[j, i].d*u[j, i]^2;
              end for;
            end for;

            for j in 1:N_passes loop
              u_pt[j]     = N_ch_p*mdot/(state_p[j].d*A_pt);
              dp_pt[j]    = 0.5*1.3*state_p[j].d*u_pt[j]^2;
            end for;

            dp_tot        = sum(dp_plates)/N_ch_p + sum(dp_pt);

        end single_phase_Coulson;

        package Base_classes "Base classes for the heat transfer"
          class base_dp "Basic pressure drop correlation"
            import VIP;
              replaceable package Medium = VIP.Media.OneRandomOrganicFluid
              "Medium model";
              parameter Integer N_cell_pc( start = 3)
              "Number of cells per channel";
              parameter Integer N_ch(start = 3) "Total number of channels";
              parameter Modelica.SIunits.Area A_pt "Port area";
              parameter Modelica.SIunits.Length Dhyd "Hydraulic diameter";
              parameter Integer N_ch_p(start = 3) "Number of channels per pass";
              parameter Modelica.SIunits.MassFlowRate mdot "Mass flow rate";
              parameter Integer N_passes(start = 1) "Number of passes";
              input Medium.ThermodynamicState state_p[N_passes] "Index matrix";
              input Modelica.SIunits.Area Aflow "Cross-sectional area";
              input Modelica.SIunits.Length l "Length";
              input Medium.ThermodynamicState state[N_ch, N_cell_pc]
              "Thermodynamic states";
          end base_dp;
        end Base_classes;
      end Plates;
    annotation (Icon(graphics={
          Rectangle(
            extent={{70,88},{90,-88}},
            lineColor={95,95,95},
            fillColor={215,215,215},
            fillPattern=FillPattern.Solid),
          Polygon(
            points={{-20,-2},{-20,20},{20,0},{-20,-20},{-20,-2}},
            lineColor={0,128,255},
            smooth=Smooth.None,
            fillColor={0,128,255},
            fillPattern=FillPattern.Solid,
            origin={33,52},
            rotation=90),
          Rectangle(
            extent={{-50,6},{50,-6}},
            lineColor={0,128,255},
            fillColor={0,128,255},
            fillPattern=FillPattern.Solid,
            origin={33,-18},
            rotation=90),
          Polygon(
            points={{-20,-2},{-20,20},{20,0},{-20,-20},{-20,-2}},
            lineColor={0,128,255},
            smooth=Smooth.None,
            fillColor={0,128,255},
            fillPattern=FillPattern.Solid,
            origin={-41,52},
            rotation=90),
          Rectangle(
            extent={{-50,6},{50,-6}},
            lineColor={0,128,255},
            fillColor={0,128,255},
            fillPattern=FillPattern.Solid,
            origin={-41,-18},
            rotation=90)}));
    end Pressure_drops;

    package Heat_transfer "A package containing heat transfer correlations"
      package Tubes "heat transfer correlations in tubes"
        extends VIP.Design.Icons.tube;
        class Dittus_Boelter "Dittus Boelter correlation for tubes"
            extends Base_classes.base_ht;
            Modelica.SIunits.NusseltNumber Nu[Ncell] "Nusselt number";
            Modelica.SIunits.CoefficientOfHeatTransfer ht[Ncell]
            "Heat transfer coefficient";
            parameter Modelica.SIunits.Length Dhyd "Hydraulic diameter";
            parameter Real alfa "exponent for the Prandtl number";
            Modelica.SIunits.ReynoldsNumber Re[Ncell](start=1e5*ones(Ncell))
            "Reynolds number";
            Modelica.SIunits.PrandtlNumber Pr[Ncell](start=4*ones(Ncell))
            "Prandtl number";
            Modelica.SIunits.Velocity u[Ncell](start=ones(Ncell)) "Velocity";
        equation
              for i in 1:Ncell loop
                u[i]   = mdot/state[i].d/Aflow;
                Re[i]  =Design.Miscellanea.numbers.Reynolds(
                    u[i],
                    state[i].d,
                    state[i].eta,
                    Dhyd);
                Pr[i]  =Design.Miscellanea.numbers.Prandtl(
                    state[i].cp,
                    state[i].eta,
                    state[i].lambda);
                assert(Re[i] > 1e4, "Reynolds number is lower than 1e4 to use Dittus and Boelter",
                  AssertionLevel.warning);
                Nu[i]  =  2.3e-2*Re[i]^0.8*Pr[i]^alfa;
                ht[i]  =Design.Miscellanea.numbers.Nusselt(
                    Nu[i],
                    state[i].lambda,
                    Dhyd);
              end for;

        end Dittus_Boelter;

        class Sieder_Tate "Sieder Tate correlation for tubes"
            extends VIP.Design.Heat_transfer.Tubes.Base_classes.base_ht;
            Modelica.SIunits.NusseltNumber Nu[Ncell] "Nusselt number";
            Modelica.SIunits.CoefficientOfHeatTransfer ht[Ncell](start=1e3*ones(Ncell))
            "Heat transfer coefficient";
            parameter Modelica.SIunits.Length Dhyd "Hydraulic diameter";
            input Modelica.SIunits.DynamicViscosity eta_wall[Ncell]
            "exponent for the viscosity correction";
            Modelica.SIunits.ReynoldsNumber Re[Ncell](start=1e5*ones(Ncell))
            "Reynolds number";
            Modelica.SIunits.PrandtlNumber Pr[Ncell](start=4*ones(Ncell))
            "Prandtl number";
            Modelica.SIunits.Velocity u[Ncell](start=ones(Ncell)) "Velocity";
        equation
              for i in 1:Ncell loop
                u[i]  = mdot/state[i].d/Aflow;
                Re[i]  =Design.Miscellanea.numbers.Reynolds(
                    u[i],
                    state[i].d,
                    state[i].eta,
                    Dhyd);
                Pr[i]  =Design.Miscellanea.numbers.Prandtl(
                    state[i].cp,
                    state[i].eta,
                    state[i].lambda);
                assert(Re[i] > 1e4,"Reynolds number is lower than 1e4 to use Sieder and Tate",
                  AssertionLevel.warning);
                assert(Pr[i] > 0.6, "Prandtl number is lower than 0.6 to be use Sieder and Tate",
                  AssertionLevel.warning);
                Nu[i] = 2.3e-2*Re[i]^0.8*Pr[i]^(1.0/3)*(eta_wall[i]
                /state[i].eta)^0.14;
                ht[i]  =Design.Miscellanea.numbers.Nusselt(
                    Nu[i],
                    state[i].lambda,
                    Dhyd);
              end for;

        end Sieder_Tate;

        class Gnielinski "Gnielinski correlation for tubes"
            extends VIP.Design.Heat_transfer.Tubes.Base_classes.base_ht;
            parameter Modelica.SIunits.Length l "Lenght";
            parameter Modelica.SIunits.Length Dhyd "Hydraulic diameter";
            Modelica.SIunits.ReynoldsNumber Re[Ncell](start=1e5*ones(Ncell))
            "Reynolds number";
            Modelica.SIunits.PrandtlNumber Pr[Ncell](start=4*ones(Ncell))
            "Prandtl number";
            Modelica.SIunits.Velocity u[Ncell](start=ones(Ncell)) "Velocity";
            Modelica.SIunits.NusseltNumber Nu[Ncell] "Nusselt number";
            Modelica.SIunits.CoefficientOfHeatTransfer ht[Ncell]
            "Heat transfer coefficient";
        protected
            Real csi[Ncell] "Friction factor";
        equation
              for i in 1:Ncell loop
                u[i]   = mdot/state[i].d/Aflow;
                Re[i]  =Design.Miscellanea.numbers.Reynolds(
                    u[i],
                    state[i].d,
                    state[i].eta,
                    Dhyd);
                Pr[i]  =Design.Miscellanea.numbers.Prandtl(
                    state[i].cp,
                    state[i].eta,
                    state[i].lambda);
                csi[i] = 1/(0.78*log(Re[i]) - 1.5)^2;
                Nu[i]  = ((csi[i]/8)*Re[i]*Pr[i])/(1 + 12.7*sqrt(csi[i]/8)*(Pr[i]^
                (2.0/3) - 1))*(1 +  (Dhyd/l)^(2.0/3));
                ht[i]  =Design.Miscellanea.numbers.Nusselt(
                    Nu[i],
                    state[i].lambda,
                    Dhyd);
              end for;
        end Gnielinski;

        class EagleFerguson
          "Eagle-Ferguson heat transfer correlation (only for liquid water)"
            extends VIP.Design.Heat_transfer.Tubes.Base_classes.base_ht;
            Modelica.SIunits.NusseltNumber Nu[Ncell] "Nusselt number tubes ";
            parameter Modelica.SIunits.Length Dhyd "Hydraulic diameter";
            Modelica.SIunits.Velocity u[Ncell](start=ones(Ncell)) "Velocity";
            Modelica.SIunits.CoefficientOfHeatTransfer ht[Ncell]
            "Heat transfer coefficient tubes";

        equation
              for i in 1:Ncell loop
                u[i]  = mdot/state[i].d/Aflow;
                Nu[i] =  4.2e3*Dhyd/state[i].lambda*(1.35 + 2e-2*(state[i].T - 273.15))
                *u[i]^0.8/(1e3*Dhyd)^0.2;
                ht[i]  =Design.Miscellanea.numbers.Nusselt(
                    Nu[i],
                    state[i].lambda,
                    Dhyd);
              end for;
        end EagleFerguson;

        package Base_classes "Base classes for the heat transfer"
          class base_ht "Basic heat transfer correlation"
            import VIP;
              replaceable package Medium = VIP.Media.OneRandomOrganicFluid
              "Medium model";
              parameter Integer Ncell(start=3) "Number of cell elements";
              input Modelica.SIunits.Area Aflow
              "Cross-sectional area (single tube)";
              input Modelica.SIunits.MassFlowRate mdot "Mass flow rate";
              input Medium.ThermodynamicState state[Ncell]
              "Thermodynamic state";
          end base_ht;
        end Base_classes;
      end Tubes;

      package Shell "heat transfer correlations in shells"
        extends VIP.Design.Icons.shell;

        class single_phase_Kern
          "Kern method for shell side only single phase to single phase"
            extends VIP.Design.Heat_transfer.Shell.Base_classes.base_Kern;
            Modelica.SIunits.NusseltNumber Nu[Ncell] "Nusselt number";
            Modelica.SIunits.CoefficientOfHeatTransfer ht[Ncell](start=1e3*ones(Ncell))
            "Heat transfer coefficient";
            Modelica.SIunits.ReynoldsNumber Re[Ncell](start=1e5*ones(Ncell))
            "Reynolds number";
            Modelica.SIunits.PrandtlNumber Pr[Ncell](start=4*ones(Ncell))
            "Prandtl number ";
            Modelica.SIunits.Length d_s_eq "Equivalent shell diameter";
            Modelica.SIunits.Velocity u[Ncell](start=ones(Ncell)) "Velocity";

        equation
            for i in 1:Ncell loop
              u[i]  = mdot/state[i].d/Aflow;
              Re[i] =Design.Miscellanea.numbers.Reynolds(
                    u[i],
                    state[i].d,
                    state[i].eta,
                    d_s_eq);
              Pr[i] =Design.Miscellanea.numbers.Prandtl(
                    state[i].cp,
                    state[i].eta,
                    state[i].lambda);
              Nu[i] = 0.6246*Re[i]^(-0.4989)*Re[i]*Pr[i]^(1/3);
              ht[i] =Design.Miscellanea.numbers.Nusselt(
                    Nu[i],
                    state[i].lambda,
                    d_s_eq);
            end for;
        end single_phase_Kern;

        class condensation_Kern
          "Kern correlation for shell side for condensation outside tube bundles"
            extends VIP.Design.Heat_transfer.Shell.Base_classes.base_Kern;
            parameter Medium.ThermodynamicState state_l
            "Thermodynamic state in saturated liquid";
            parameter Medium.ThermodynamicState state_v
            "Thermodynamic state in saturated vapor";
            parameter Boolean shear_vapor = true
            "= true, if shear vapor effect is considered";
            input Modelica.SIunits.Length l "Lenght";
            input Real N_tubes "Number of tubes in the bundle";
            Modelica.SIunits.CoefficientOfHeatTransfer ht[Ncell]
            "Heat transfer coefficient";
            Modelica.SIunits.CoefficientOfHeatTransfer ht_l
            "Heat transfer coefficient if all liquid phase";
            Modelica.SIunits.NusseltNumber ht_g[Ncell]
            "Nusselt number gravity-controlled";
            Modelica.SIunits.NusseltNumber ht_sh[Ncell]
            "Nusselt number shear-controlled condensation";
            Modelica.SIunits.Length d_s_eq "Equivalent shell diameter";
            Modelica.SIunits.ReynoldsNumber Re[Ncell](start=1e5*ones(Ncell))
            "Reynolds number";
            Modelica.SIunits.PrandtlNumber Pr[Ncell](start=4*ones(Ncell))
            "Prandtl number";
            Modelica.SIunits.ReynoldsNumber Re_l(start=1e5)
            "Reynolds number if all liquid phase";
            Modelica.SIunits.PrandtlNumber Pr_l(start=4)
            "Prandtl number if all liquid phase";
            Modelica.SIunits.Velocity u[Ncell](start=ones(Ncell)) "Velocity";
            Modelica.SIunits.Velocity u_l(start=1)
            "Velocity if all liquid phase";
            Real X_tt[Ncell] "Martinelli parameter";
            Real xq[Ncell] "Vapor quality";

        equation
            u_l   = mdot/state_l.d/Aflow;
            Re_l  =Design.Miscellanea.numbers.Reynolds(
                  u_l,
                  state_l.d,
                  state_l.eta,
                  d_s_eq);
            Pr_l  =Design.Miscellanea.numbers.Prandtl(
                  state_l.cp,
                  state_l.eta,
                  state_l.lambda);
            ht_l  = 0.6246*Re_l^(-0.4989)*Re_l*Pr_l^(1/3)*state_l.lambda/d_s_eq;

            for i in 1:Ncell loop

              u[i]  = mdot/state[i].d/Aflow "tube velocity";
              Re[i] =Design.Miscellanea.numbers.Reynolds(
                    u[i],
                    state[i].d,
                    state[i].eta,
                    d_s_eq);
              Pr[i] =Design.Miscellanea.numbers.Prandtl(
                    state[i].cp,
                    state[i].eta,
                    state[i].lambda);

              if (state[i].h < state_l.h or state[i].h > state_v.h) then

                ht_g[i]  = 0.6246*Re[i]^(-0.4989)*Re[i]*Pr[i]^(1/3)*state[i].lambda/d_s_eq;
                xq[i]    = 1;
                X_tt[i]  = 0;
                ht_sh[i] = ht_g[i];

              else

                xq[i]    = (state[i].h - state_l.h)/(state_v.h - state_l.h);
                X_tt[i]  = (1/xq[i] - 1)^0.9*(state_v.d/state_l.d)^0.5*(state_l.eta
                /state_v.eta)^0.1;
                ht_sh[i] = 1.26*(1/X_tt[i])^0.78*ht_l;
                ht_g[i]  = 0.95*state_l.lambda*(g_n*l*N_tubes^(2/3)*state_l.d*(state_l.d
                - state_v.d)/(state_l.eta*mdot))^(1/3);

              end if;

              if shear_vapor then
                ht[i]    = sqrt(ht_g[i]^2 + ht_sh[i]^2);
              else
                ht[i]    = ht_g[i];
              end if;

            end for;

        end condensation_Kern;

        class single_phase_Bell_Delaware
          "Bell method for shell side only single phase to single phase"
            extends
            VIP.Design.Heat_transfer.Shell.Base_classes.base_Bell_Delaware;
            Modelica.SIunits.CoefficientOfHeatTransfer ht[Ncell]
            "Heat transfer coefficient";
            Modelica.SIunits.ReynoldsNumber Re[Ncell](start=1e5*ones(Ncell))
            "Reynolds number";
            Modelica.SIunits.PrandtlNumber Pr[Ncell](start=4*ones(Ncell))
            "Prandtl number";
            Modelica.SIunits.Velocity u[Ncell](start=ones(Ncell)) "Velocity";
        protected
            Real csi[Ncell] "Friction factor";
            Real a1[Ncell];
            Real a2[Ncell];
            Real a[Ncell];
        equation

          for i in 1:Ncell loop

            u[i]  = mdot/state[i].d/S_m;
            Re[i] =Design.Miscellanea.numbers.Reynolds(
                    u[i],
                    state[i].d,
                    state[i].eta,
                    Dhyd_o);
            Pr[i] =Design.Miscellanea.numbers.Prandtl(
                    state[i].cp,
                    state[i].eta,
                    state[i].lambda);

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

            csi[i] = a1[i]*(1.33/pitch_f)^a[i]*Re[i]^a2[i];
            ht[i]  = J_l*J_b*J_c*csi[i]*state[i].cp*mdot/(S_m*Pr[i]^(2/3));

          end for;

        end single_phase_Bell_Delaware;

        class condensation_Bell_Delaware
          "Bell method for shell side for condensation outside tube bundles"
            extends
            VIP.Design.Heat_transfer.Shell.Base_classes.base_Bell_Delaware;
            parameter Medium.ThermodynamicState state_l
            "Thermodynamic state in saturated liquid";
            parameter Medium.ThermodynamicState state_v
            "Thermodynamic state in saturated vapor";
            parameter Boolean shear_vapor = true
            "= true, if shear vapor effect is considered";
            input Modelica.SIunits.Length l "Lenght";
            input Real N_tubes "Number of tubes in the bundle";
            Modelica.SIunits.CoefficientOfHeatTransfer ht[Ncell]
            "Heat transfer coefficient";
            Modelica.SIunits.CoefficientOfHeatTransfer ht_l
            "Heat transfer coefficient if all liquid phase";
            Modelica.SIunits.NusseltNumber ht_g[Ncell]
            "Nusselt number gravity-controlled";
            Modelica.SIunits.NusseltNumber ht_sh[Ncell]
            "Nusselt number shear-controlled condensation";
            Modelica.SIunits.ReynoldsNumber Re[Ncell](start=1e5*ones(Ncell))
            "Reynolds number";
            Modelica.SIunits.PrandtlNumber Pr[Ncell](start=4*ones(Ncell))
            "Prandtl number";
            Modelica.SIunits.ReynoldsNumber Re_l(start=1e5)
            "Reynolds number if all liquid phase";
            Modelica.SIunits.PrandtlNumber Pr_l(start=4)
            "Prandtl number if all liquid phase";
            Modelica.SIunits.Velocity u[Ncell](start=ones(Ncell)) "Velocity";
            Modelica.SIunits.Velocity u_l(start=1)
            "Velocity if all liquid phase";
            Real X_tt[Ncell] "Martinelli parameter";
            Real xq[Ncell] "Vapor quality";

        protected
            Real csi[Ncell] "Friction factor";
            Real a1[Ncell];
            Real a2[Ncell];
            Real a[Ncell];
            Real csi_l "Friction factor if all liquid phase";
            Real a1_l;
            Real a2_l;
            Real a_l;
        equation

          u_l   = mdot/state_l.d/S_m;
          Re_l  =Design.Miscellanea.numbers.Reynolds(
                  u_l,
                  state_l.d,
                  state_l.eta,
                  Dhyd_o);
          Pr_l  =Design.Miscellanea.numbers.Prandtl(
                  state_l.cp,
                  state_l.eta,
                  state_l.lambda);

          if (Re_l <= 1e1) then
            if (layout == 1) then
              a1_l  = 1.4;
              a2_l  = -0.667;
            else
              a1_l  = 0.97;
              a2_l  = -0.667;
            end if;
          elseif (Re_l > 1e1 and Re_l <= 1e2) then
            if (layout == 1) then
              a1_l  = 1.36;
              a2_l  = -0.657;
            else
              a1_l  = 0.9;
              a2_l  = -0.631;
            end if;
          elseif (Re_l > 1e2 and Re_l <= 1e3) then
            if (layout == 1) then
              a1_l  = 0.593;
              a2_l  = -0.477;
            else
              a1_l  = 0.408;
              a2_l  = -0.460;
            end if;
          elseif (Re_l > 1e3 and Re_l < 1e4) then
            if (layout == 1) then
              a1_l  = 0.321;
              a2_l  = -0.388;
            else
              a1_l  = 0.107;
              a2_l  = -0.266;
            end if;
          else
            if (layout == 1) then
              a1_l  = 0.321;
              a2_l  = -0.388;
            else
              a1_l  = 0.370;
              a2_l  = -0.395;
            end if;
          end if;

          if (layout == 1) then
            a_l   = 1.450/(1 + 0.14*Re_l^0.519);
          else
            a_l   = 1.187/(1 + 0.14*Re_l^0.370);
          end if;

          csi_l = a1_l*(1.33/pitch_f)^a_l*Re_l^a2_l;
          ht_l  = J_l*J_b*J_c*csi_l*state_l.cp*mdot/(S_m*Pr_l^(2/3));

          for i in 1:Ncell loop

            u[i]  = mdot/state[i].d/S_m;
            Re[i] =Design.Miscellanea.numbers.Reynolds(
                    u[i],
                    state[i].d,
                    state[i].eta,
                    Dhyd_o);
            Pr[i] =Design.Miscellanea.numbers.Prandtl(
                    state[i].cp,
                    state[i].eta,
                    state[i].lambda);

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

            csi[i] = a1[i]*(1.33/pitch_f)^a[i]*Re[i]^a2[i];

            if (state[i].h < state_l.h or state[i].h > state_v.h) then

              ht_g[i]  = J_l*J_b*J_c*csi[i]*state[i].cp*mdot/(S_m*Pr[i]^(2/3));
              xq[i]    = 1;
              X_tt[i]  = 0;
              ht_sh[i] = ht_g[i];

            else

              xq[i]    = (state[i].h - state_l.h)/(state_v.h - state_l.h);
              X_tt[i]  = (1/xq[i] - 1)^0.9*(state_v.d/state_l.d)^0.5*(state_l.eta
              /state_v.eta)^0.1;
              ht_sh[i] = 1.26*(1/X_tt[i])^0.78*ht_l;
              ht_g[i]  = 0.95*state_l.lambda*(g_n*l*N_tubes^(2/3)*state_l.d*(state_l.d
                - state_v.d)/(state_l.eta*mdot))^(1/3);

            end if;

            if shear_vapor then
              ht[i]    = sqrt(ht_g[i]^2 + ht_sh[i]^2);
            else
              ht[i]    = ht_g[i];
            end if;

          end for;

        end condensation_Bell_Delaware;

        package Base_classes "Base classes for the heat transfer"
          class base_Kern "Kern correlation for shell side base class"
              extends VIP.Design.Heat_transfer.Shell.Base_classes.base_ht;
              input Modelica.SIunits.Length Dhyd_o "Outer hydraulic diameter";
              parameter Real pitch_f "Tube pitch";
              parameter Integer layout = 1
              "Tube layout 1 = triangular, 2 = squared";
              Modelica.SIunits.Length d_s_eq "Equivalent shell diameter";

          equation
              if layout == 1 then
                d_s_eq  = 1.1*Dhyd_o*(pitch_f^2 - 0.917);
              elseif layout == 2 then
                d_s_eq  = 1.27*Dhyd_o*(pitch_f^2 - 0.785);
              else
                d_s_eq  = 1;
              end if;
          end base_Kern;

          class base_Bell_Delaware "Bell Delaware correlation for shell side"
              extends VIP.Design.Heat_transfer.Shell.Base_classes.base_ht;
              input Modelica.SIunits.Length Dhyd "Hydraulic diameter";
              input Modelica.SIunits.Length Dhyd_o "Outer hydraulic diameter";
              parameter Real b_cut "Baffle cut";
              parameter Real ttb = 8e-4 "Tube to baffle clearance";
              parameter Real bts = 4.8e-3 "Baffle to shell clearance";
              parameter Real pitch_f "Tube pitch";
              parameter Real N_ss
              "The number of sealing strips (pairs) in one baffle spacing";
              parameter Integer layout
              "Tube layout, 1 = triangular, 2 = squared";
              parameter Integer N_passes "Number of tube passes";
              parameter Real ff = 1 "Fouling factor (=1 if clean)";
              input Modelica.SIunits.Length d_s "Shell diameter";
              input Modelica.SIunits.Length d_b "Bundle diameter";
              input Real N_tubes "Number of tubes in the bundle";
              input Modelica.SIunits.Length l_b "Baffle length";
              Modelica.SIunits.Length tpv "Vertical tube pitch";
              Real N_tcc "Number of tube rows between the baffle tips";
              Real N_tcw
              "The effective number of tube rows crossed in the baffle window";
              Real F_sbp
              "The ratio of the bypass area to the overall crossflow area";
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
              Modelica.SIunits.Angle teta_ctl
              "The upper centriangle of baffle cut";

          equation
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

          end base_Bell_Delaware;

          class base_ht "Basic heat transfer correlation"
              replaceable package Medium = VIP.Media.OneRandomOrganicFluid
              "Medium model";
              parameter Integer Ncell(start=3) "Number of cell elements";
              input Modelica.SIunits.Area Aflow
              "Cross-sectional area (single tube)";
              input Modelica.SIunits.MassFlowRate mdot "Mass flow rate";
              input Medium.ThermodynamicState state[Ncell];
          end base_ht;
        end Base_classes;
      end Shell;

      package Plates
        extends Icons.plate;
        class single_phase_Martin
          "Martin heat transfer correlation for single phase heat transfer"
            extends Base_classes.base_ht;
            parameter Modelica.SIunits.Angle beta "Plate inclination angle";
            Modelica.SIunits.ReynoldsNumber Re[N_ch, N_cell_pc]
            "Reynolds number";
            Modelica.SIunits.PrandtlNumber Pr[N_ch, N_cell_pc] "Prandtl number";
            Modelica.SIunits.Velocity u[N_ch, N_cell_pc] "Velocity";
            Modelica.SIunits.NusseltNumber Nu[N_ch, N_cell_pc] "Nusselt number";
            Modelica.SIunits.CoefficientOfHeatTransfer ht[N_ch, N_cell_pc]
            "Heat transfer coefficient";
        protected
            Real csi[N_ch, N_cell_pc] "Friction factor";
            Real csi0[N_ch, N_cell_pc] "Friction factor 0";
            Real csi1[N_ch, N_cell_pc] "Friction factor 1";
        equation
            for j in 1:N_ch loop
              for i in 1:N_cell_pc loop
                u[j, i]   = mdot/state[j, i].d/Aflow;
                Re[j, i]  =Design.Miscellanea.numbers.Reynolds(
                      u[j, i],
                      state[j, i].d,
                      state[j, i].eta,
                      Dhyd);
                Pr[j, i]  =Design.Miscellanea.numbers.Prandtl(
                      state[j, i].cp,
                      state[j, i].eta,
                      state[j, i].lambda);
                if (Re[j, i] < 2e3) then
                  csi0[j, i] = 16/Re[j, i];
                  csi1[j, i] = 149.25/Re[j, i] + 0.9625;
                else
                  csi0[j, i] = 1/(1.56*log(Re[j, i]) - 3)^2;
                  csi1[j, i] = 9.75/Re[j, i]^0.289;
                end if;
                csi[j, i]       = 1/((1 - cos(beta))/sqrt(3.8*csi1[j, i]) + cos(beta)/
                sqrt(0.045*tan(beta) + 0.09*sin(beta) + csi0[j, i]/cos(beta)))^2;
                Nu[j, i]  = 0.205*Pr[j, i]^(1/3)*(csi[j, i]*Re[j, i]^2*sin(2*beta))
                ^0.374;
                ht[j, i]  =Design.Miscellanea.numbers.Nusselt(
                      Nu[j, i],
                      state[j, i].lambda,
                      Dhyd);
              end for;
            end for;
        end single_phase_Martin;

        class evaporation_Coulson
          "Cooper heat transfer correlation for evaporation heat transfer. Coulson for single phase."
            extends Base_classes.base_ht;
            parameter Medium.ThermodynamicState state_l
            "Thermodynamic state in saturated liquid";
            parameter Medium.ThermodynamicState state_v
            "Thermodynamic state in saturated vapor";
            parameter Modelica.SIunits.Length R_p = 1e-6
            "Relative roughness of the surface";
            input Modelica.SIunits.Area At "Width";
            parameter Modelica.SIunits.MolarMass M = Medium.fluidConstants[1].molarMass
            "Molar mass of the fluid";
            parameter Real pc = Medium.fluidConstants[1].criticalPressure
            "Critical pressure";
            parameter Modelica.SIunits.AbsolutePressure p_in "Inlet pressure";
            parameter Modelica.SIunits.HeatFlux qdot_tilde_start
            "Heat flux start";
            input Modelica.SIunits.Length l "Length";
            input Modelica.SIunits.HeatFlowRate qdot[N_ch, N_cell_pc]
            "Heat rate";
            Modelica.SIunits.ReynoldsNumber Re[N_ch, N_cell_pc]
            "Reynolds number";
            Modelica.SIunits.PrandtlNumber Pr[N_ch, N_cell_pc] "Prandtl number";
            Modelica.SIunits.Velocity u[N_ch, N_cell_pc] "Velocity";
            Modelica.SIunits.CoefficientOfHeatTransfer ht[N_ch, N_cell_pc]
            "Heat transfer coefficient";
            Modelica.SIunits.HeatFlux qdot_tilde[N_ch, N_cell_pc](
             start = fill(qdot_tilde_start, N_ch, N_cell_pc)) "Heat flux";

        equation
            for j in 1:N_ch loop
              for i in 1:N_cell_pc loop
                qdot_tilde[j, i] = qdot[j, i]/(2*At);
                u[j, i]          = mdot/state[j, i].d/Aflow;
                Re[j, i]         =Design.Miscellanea.numbers.Reynolds(
                      u[j, i],
                      state[j, i].d,
                      state[j, i].eta,
                      Dhyd);
                Pr[j, i]         =Design.Miscellanea.numbers.Prandtl(
                      state[j, i].cp,
                      state[j, i].eta,
                      state[j, i].lambda);
                if (state[j, i].h < state_l.h or state[j, i].h > state_v.h) then
                  ht[j, i]       = 0.26*Pr[j, i]^0.4*Re[j, i]^0.65*state[j, i].lambda/Dhyd;
                else
                  ht[j, i]      = 55*(p_in/pc)^(0.12 - 0.2*log10(R_p))/(-log10(
                  p_in/pc))^0.55*abs(qdot_tilde[j, i])^0.67/sqrt(1e3*M);
                end if;
              end for;
            end for;
        end evaporation_Coulson;

        class evaporation_Martin
          "Cooper heat transfer correlation for evaporation heat transfer. Martin for single phase."
            extends Base_classes.base_ht;
            parameter Modelica.SIunits.Angle beta "Plate inclination angle";
            parameter Medium.ThermodynamicState state_l
            "Thermodynamic state in saturated liquid";
            parameter Medium.ThermodynamicState state_v
            "Thermodynamic state in saturated vapor";
            parameter Modelica.SIunits.Length R_p = 1e-6
            "Relative roughness of the surface";
            input Modelica.SIunits.Area At "Width";
            parameter Modelica.SIunits.MolarMass M = Medium.fluidConstants[1].molarMass
            "Molar mass of the fluid";
            parameter Real pc = Medium.fluidConstants[1].criticalPressure
            "Critical pressure";
            parameter Modelica.SIunits.AbsolutePressure p_in "Inlet pressure";
            parameter Modelica.SIunits.HeatFlux qdot_tilde_start
            "Heat flux start";
            input Modelica.SIunits.Length l "Length";
            input Modelica.SIunits.HeatFlowRate qdot[N_ch, N_cell_pc]
            "Heat rate";
            Modelica.SIunits.ReynoldsNumber Re[N_ch, N_cell_pc]
            "Reynolds number";
            Modelica.SIunits.PrandtlNumber Pr[N_ch, N_cell_pc] "Prandtl number";
            Modelica.SIunits.Velocity u[N_ch, N_cell_pc] "Velocity";
            Modelica.SIunits.CoefficientOfHeatTransfer ht[N_ch, N_cell_pc]
            "Heat transfer coefficient";
            Modelica.SIunits.HeatFlux qdot_tilde[N_ch, N_cell_pc](
             start = fill(qdot_tilde_start, N_ch, N_cell_pc)) "Heat flux";
        protected
             Real csi[N_ch, N_cell_pc] "Friction factor";
             Real csi0[N_ch, N_cell_pc] "Friction factor 0";
             Real csi1[N_ch, N_cell_pc] "Friction factor 1";
        equation

            for j in 1:N_ch loop
              for i in 1:N_cell_pc loop
                qdot_tilde[j, i] = qdot[j, i]/(2*At);
                u[j, i]          = mdot/state[j, i].d/Aflow;
                Re[j, i]         =Design.Miscellanea.numbers.Reynolds(
                      u[j, i],
                      state[j, i].d,
                      state[j, i].eta,
                      Dhyd);
                Pr[j, i]         =Design.Miscellanea.numbers.Prandtl(
                      state[j, i].cp,
                      state[j, i].eta,
                      state[j, i].lambda);
                if (Re[j, i] < 2e3) then
                  csi0[j, i] = 16/Re[j, i];
                  csi1[j, i] = 149.25/Re[j, i] + 0.9625;
                else
                  csi0[j, i] = 1/(1.56*log(Re[j, i]) - 3)^2;
                  csi1[j, i] = 9.75/Re[j, i]^0.289;
                end if;
                csi[j, i]       = 1/((1 - cos(beta))/sqrt(3.8*csi1[j, i]) + cos(beta)/
                sqrt(0.045*tan(beta) + 0.09*sin(beta) + csi0[j, i]/cos(beta)))^2;
                if (state[j, i].h < state_l.h or state[j, i].h > state_v.h) then
                  ht[j, i]      = 0.205*Pr[j, i]^(1/3)*(csi[j, i]*Re[j, i]^2*sin(2*beta))
                  ^0.374*state[j, i].lambda/Dhyd;
                else
                   ht[j, i]      = 55*(p_in/pc)^(0.12 - 0.2*log10(R_p))/(-log10(
                   p_in/pc))^0.55*abs(qdot_tilde[j, i])^0.67/sqrt(1e3*M);
                end if;
              end for;
            end for;
        end evaporation_Martin;

        class single_phase_Coulson
          "Coulson heat transfer correlation for single phase heat transfer"
            extends Base_classes.base_ht;
            Modelica.SIunits.ReynoldsNumber Re[N_ch, N_cell_pc]
            "Reynolds number";
            Modelica.SIunits.PrandtlNumber Pr[N_ch, N_cell_pc] "Prandtl number";
            Modelica.SIunits.Velocity u[N_ch, N_cell_pc] "Velocity";
            Modelica.SIunits.NusseltNumber Nu[N_ch, N_cell_pc] "Nusselt number";
            Modelica.SIunits.CoefficientOfHeatTransfer ht[N_ch, N_cell_pc]
            "Heat transfer coefficient";

        equation
            for j in 1:N_ch loop
              for i in 1:N_cell_pc loop
                u[j, i]   = mdot/state[j, i].d/Aflow;
                Re[j, i]  =Design.Miscellanea.numbers.Reynolds(
                      u[j, i],
                      state[j, i].d,
                      state[j, i].eta,
                      Dhyd);
                Pr[j, i]  =Design.Miscellanea.numbers.Prandtl(
                      state[j, i].cp,
                      state[j, i].eta,
                      state[j, i].lambda);
                Nu[j, i]  = 0.26*Pr[j, i]^0.4*Re[j, i]^0.65;
                ht[j, i]  =Design.Miscellanea.numbers.Nusselt(
                      Nu[j, i],
                      state[j, i].lambda,
                      Dhyd);
              end for;
            end for;
        end single_phase_Coulson;

        package Base_classes "Base classes for the heat transfer"
          class base_ht "Basic heat transfer correlation"
            import VIP;
              replaceable package Medium = VIP.Media.OneRandomOrganicFluid
              "Medium model";
              parameter Integer N_ch(start = 3) "Total number of channels";
              parameter Integer N_cell_pc( start = 3)
              "Number of cells per channel";
              parameter Modelica.SIunits.Length Dhyd "Hydraulic diameter";
              input Modelica.SIunits.Area Aflow "Cross-sectional area";
              parameter Modelica.SIunits.MassFlowRate mdot "Mass flow rate";
              input Medium.ThermodynamicState state[N_ch, N_cell_pc]
              "Thermodynamic states";
          end base_ht;
        end Base_classes;
      end Plates;
    annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,
              -100},{100,100}}), graphics={
          Rectangle(
            extent={{60,88},{80,-88}},
            lineColor={95,95,95},
            fillColor={215,215,215},
            fillPattern=FillPattern.Solid),
          Polygon(
            points={{6,58},{6,80},{46,60},{6,40},{6,58}},
            lineColor={255,0,0},
            smooth=Smooth.None,
            fillColor={255,0,0},
            fillPattern=FillPattern.Solid),
          Rectangle(
            extent={{-94,66},{6,54}},
            lineColor={255,0,0},
            fillColor={255,0,0},
            fillPattern=FillPattern.Solid),
          Polygon(
            points={{6,-42},{6,-20},{46,-40},{6,-60},{6,-42}},
            lineColor={255,0,0},
            smooth=Smooth.None,
            fillColor={255,0,0},
            fillPattern=FillPattern.Solid),
          Rectangle(
            extent={{-94,-34},{6,-46}},
            lineColor={255,0,0},
            fillColor={255,0,0},
            fillPattern=FillPattern.Solid)}));
    end Heat_transfer;

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

      partial model flat_plate "Flat plate heat exchanger icon"

        annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{
                  -100,-100},{100,100}}),
                         graphics={
              Line(
                points={{-60,80},{90,80}},
                color={255,0,0},
                smooth=Smooth.None),
              Line(
                points={{66,52}},
                color={0,0,0},
                smooth=Smooth.None),
              Line(
                points={{-60,80},{-60,-60}},
                color={255,0,0},
                smooth=Smooth.None),
              Line(
                points={{20,80},{20,-60}},
                color={255,0,0},
                smooth=Smooth.None),
              Line(
                points={{-90,-60},{20,-60}},
                color={255,0,0},
                smooth=Smooth.None),
              Line(
                points={{-20,49},{90,49}},
                color={0,128,255},
                smooth=Smooth.None),
              Line(
                points={{-90,-91},{60,-91}},
                color={0,128,255},
                smooth=Smooth.None),
              Line(
                points={{-66,6},{-60,0},{-54,6}},
                color={255,0,0},
                smooth=Smooth.None),
              Line(
                points={{14,6},{20,0},{26,6}},
                color={255,0,0},
                smooth=Smooth.None),
              Line(
                points={{-20,49},{-20,-91}},
                color={0,128,255},
                smooth=Smooth.None),
              Line(
                points={{6,3},{0,-3},{-6,3}},
                color={0,128,255},
                smooth=Smooth.None,
                origin={-20,-3},
                rotation=180),
              Line(
                points={{60,49},{60,-91}},
                color={0,128,255},
                smooth=Smooth.None),
              Line(
                points={{6,3},{0,-3},{-6,3}},
                color={0,128,255},
                smooth=Smooth.None,
                origin={60,-3},
                rotation=180),
              Rectangle(
                extent={{-42,46},{-38,-54}},
                lineColor={0,0,0},
                fillColor={215,215,215},
                fillPattern=FillPattern.Forward),
              Rectangle(
                extent={{-2,46},{2,-54}},
                lineColor={0,0,0},
                fillColor={215,215,215},
                fillPattern=FillPattern.Forward),
              Rectangle(
                extent={{38,46},{42,-54}},
                lineColor={0,0,0},
                fillColor={215,215,215},
                fillPattern=FillPattern.Forward)}));
      end flat_plate;

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

      partial package plate "Flat plate icon"

        annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{
                  -100,-100},{100,100}}),
                         graphics={
              Line(
                points={{66,52}},
                color={0,0,0},
                smooth=Smooth.None),     Polygon(
                points={{-30,-70},{30,-50},{30,90},{-30,70},{-30,-70}},
                lineColor={0,0,0},
                smooth=Smooth.None,
                fillColor={215,215,215},
                fillPattern=FillPattern.Forward)}));
      end plate;
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
              radius=25.0),                    Polygon(
                origin={-8.167,-17},
                fillColor={128,128,128},
                pattern=LinePattern.None,
                fillPattern=FillPattern.Solid,
                points={{-15.833,20.0},{-15.833,30.0},{14.167,40.0},{24.167,20.0},{
                    4.167,-30.0},{14.167,-30.0},{24.167,-30.0},{24.167,-40.0},{-5.833,
                    -50.0},{-15.833,-30.0},{4.167,20.0},{-5.833,20.0}},
                smooth=Smooth.Bezier,
                lineColor={0,0,0}), Ellipse(
                origin={-0.5,56.5},
                fillColor={128,128,128},
                pattern=LinePattern.None,
                fillPattern=FillPattern.Solid,
                extent={{-12.5,-12.5},{12.5,12.5}},
                lineColor={0,0,0})}));
    end Icons;
  end Design;

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
                                                                                              </html>"),    Icon(coordinateSystem(extent={{-100,
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
           p_h_out_start, T_h_out_start) "Outlet enthalpy start value hot side"
                                                                                annotation (Dialog(tab="Start"));
          parameter Medium_c.SpecificEnthalpy h_c_in_start = Medium_c.specificEnthalpy_pT(
           p_c_in_start, T_c_in_start) "Inlet enthalpy start value cold side" annotation (Dialog(tab="Start"));
          parameter Medium_c.SpecificEnthalpy h_c_out_start = Medium_c.specificEnthalpy_pT(
           p_c_out_start, T_c_out_start)
            "Outlet enthalpy start value cold side"                              annotation (Dialog(tab="Start"));
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
          h(start = h_h_in_start), p(start = p_h_in_start))
            "Inlet node hot side"                                                 annotation (Placement(transformation(extent={{-10,-110},{10,-90}}),
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

  package Media "I am the packag containing the definition of the fluids"
    package OneRandomOrganicFluid
      "Change the name to something more appropriate"
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
  annotation (uses(Modelica(version="3.2.1")),
                                             Icon(coordinateSystem(
          preserveAspectRatio=false, extent={{-100,-100},{100,100}}), graphics));
end VIP;
