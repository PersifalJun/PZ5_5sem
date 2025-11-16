package lab_NewtonCotes_Gauss.integration;

import lab_NewtonCotes_Gauss.point.Point;

import java.util.ArrayList;
import java.util.List;

public abstract class IntegrationScheme {

    protected List<Point> points = new ArrayList<>();

    protected List<Double> weight = new ArrayList<>();

    public enum IntegrationSchemeType {
        Gauss2,
        Simpson
    }
}