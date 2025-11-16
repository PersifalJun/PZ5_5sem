package lab_NewtonCotes_Gauss.integration;

import lab_NewtonCotes_Gauss.point.Point;

@FunctionalInterface
public interface Integrand {
    double apply(Point p);
}
