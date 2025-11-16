package lab_NewtonCotes_Gauss;

import lab_NewtonCotes_Gauss.integration.Integrand;
import lab_NewtonCotes_Gauss.integration.IntegrationScheme;
import lab_NewtonCotes_Gauss.integration.IntegrationSchemeInterval;
import lab_NewtonCotes_Gauss.point.Point;

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Locale;

public class Main {
    public static void main(String[] args) throws IOException {

        double a = 0.03;
        double b = 0.30;

        Integrand phi = p -> Math.log(1.0 + p.x());

        double Istar = exactIntegral(a, b);

        int nBase = 10;

        int[] segments = { nBase, 2 * nBase, 4 * nBase, 8 * nBase };

        IntegrationScheme.IntegrationSchemeType[] types = {
                IntegrationScheme.IntegrationSchemeType.Gauss2,
                IntegrationScheme.IntegrationSchemeType.Simpson
        };

        String[] typeNames = { "gauss2", "simpson" };

        Point begin = new Point(a, 0.0, 0.0);
        Point end   = new Point(b, 0.0, 0.0);

        for (int t = 0; t < types.length; t++) {
            IntegrationScheme.IntegrationSchemeType type = types[t];
            IntegrationSchemeInterval scheme = new IntegrationSchemeInterval(type);

            double[] I = new double[segments.length];

            for (int i = 0; i < segments.length; i++) {
                int n = segments[i];
                double h = (b - a) / n;

                double value = scheme.calculateIntegral(begin, end, n, phi);
                I[i] = value;

                String fileName = String.format("%s_N%d.txt", typeNames[t], n);

                try (PrintWriter out =
                             new PrintWriter(new FileWriter(fileName))) {

                    out.printf(Locale.US, "Scheme: %s%n", type);
                    out.printf(Locale.US, "a = %.8f, b = %.8f%n", a, b);
                    out.printf(Locale.US, "N = %d%n", n);
                    out.printf(Locale.US, "h = %.12e%n", h);
                    out.printf(Locale.US, "Integral ≈ %.15e%n", value);
                }
            }

            double pRaw = Math.log(
                    Math.abs(I[0] - I[1]) / Math.abs(I[1] - I[2])
            ) / Math.log(2.0);

            long pRoundedInt = Math.round(pRaw);

            System.out.println("Scheme " + type + ": ");

            int k = 4;

            String richFileName = typeNames[t] + "_richardson.txt";
            try (PrintWriter out = new PrintWriter(new FileWriter(richFileName))) {

                out.printf("%10s %15s %25s %15s %25s %15s %15s%n",
                        "h",
                        "I*-I^h",
                        "(I*-I^h)/(I*-I^{h/2})",
                        "p(h)",
                        "(I^{h/2}-I^h)/(2^k-1)",
                        "I^R",
                        "I*-I^R");

                int m = segments.length;
                double[] hArr = new double[m];
                double[] EArr = new double[m];

                for (int i = 0; i < m; i++) {
                    hArr[i] = (b - a) / segments[i];
                    EArr[i] = Istar - I[i];
                }

                for (int i = 0; i < m - 1; i++) {
                    double hVal   = hArr[i];
                    double Eh     = EArr[i];
                    double EhNext = EArr[i + 1];

                    double ratio = Eh / EhNext;
                    double pLocal = Math.log(Math.abs(ratio)) / Math.log(2.0);

                    double richTerm = (I[i + 1] - I[i]) /
                            (Math.pow(2.0, k) - 1.0);
                    double IR   = I[i + 1] + richTerm;
                    double ErrR = Istar - IR;
                    System.out.println("h = " + hVal + ": p(h)" + pLocal + " ≈ " + pRoundedInt);

                    out.printf(Locale.US,
                            "%10.3e %15.6e %25.6e %15.6e %25.6e %15.6e %15.6e%n",
                            hVal, Eh, ratio, pLocal, richTerm, IR, ErrR);
                }
            }

            System.out.println();
        }
    }

    private static double exactIntegral(double a, double b) {
        return F(b) - F(a);
    }

    private static double F(double x) {
        return (1.0 + x) * Math.log(1.0 + x) - x;
    }
}
