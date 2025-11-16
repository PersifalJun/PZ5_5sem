package lab_NewtonCotes_Gauss.integration;

import lab_NewtonCotes_Gauss.point.Point;

import java.util.Arrays;


public class IntegrationSchemeInterval extends IntegrationScheme {

    public IntegrationSchemeInterval(IntegrationSchemeType type) {
        switch (type) {
            case Gauss2: {
                double invSqrt3 = 1.0 / Math.sqrt(3.0);
                weight = Arrays.asList(1.0, 1.0);
                points = Arrays.asList(
                        new Point(-invSqrt3, 0.0, 0.0),
                        new Point( invSqrt3, 0.0, 0.0)
                );
                break;
            }
            case Simpson: {
                weight = Arrays.asList(1.0 / 3.0, 4.0 / 3.0, 1.0 / 3.0);
                points = Arrays.asList(
                        new Point(-1.0, 0.0, 0.0),
                        new Point( 0.0, 0.0, 0.0),
                        new Point( 1.0, 0.0, 0.0)
                );
                break;
            }
        }
    }

    public double calculateIntegral(
            Point begin,
            Point end,
            int numberSegments,
            Integrand func) {

        double result = 0.0;
        double a = begin.x();
        double b = end.x();

        double h = (b - a) / numberSegments;

        for (int i = 0; i < numberSegments; i++) {
            double x0 = a + i * h;

            for (int k = 0; k < points.size(); k++) {
                double t = points.get(k).x();
                double x = x0 + (1.0 + t) * h / 2.0;
                Point p = new Point(x, 0.0, 0.0);
                result += weight.get(k) * func.apply(p);
            }
        }

        return result * (h / 2.0);
    }
}
