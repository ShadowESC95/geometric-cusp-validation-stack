// COMSOL Java API template. Create a 1D component and solve
// -D*Sxx + m_eff(x)^2*S = sigmaIR(x).
// Run from the COMSOL Java API environment; API names can differ slightly by release.
import com.comsol.model.*;
import com.comsol.model.util.*;

public class GeometricBoundaryPDE {
  public static Model run() {
    Model model = ModelUtil.create("Model");
    model.modelPath(".");
    model.component().create("comp1", true);
    model.component("comp1").geom().create("geom1", 1);
    model.component("comp1").geom("geom1").lengthUnit("1");
    model.component("comp1").geom("geom1").feature().create("i1", "Interval");
    model.component("comp1").geom("geom1").feature("i1").set("size", 12.0);
    model.component("comp1").geom("geom1").run();
    model.component("comp1").physics().create("gpb", "GeneralFormPDE", "geom1");
    model.component("comp1").physics("gpb").field("dimensionless").component(new String[]{"S"});
    model.component("comp1").physics("gpb").feature("gfeq1").set("f", "exp(-x^2)");
    model.study().create("std1");
    model.study("std1").feature().create("stat", "Stationary");
    model.study("std1").run();
    return model;
  }
}
