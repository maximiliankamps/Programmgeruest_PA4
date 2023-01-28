package ode;

/**
 * Der klassische Runge-Kutta der Ordnung 4
 *
 * @author braeckle
 *
 */
public class RungeKutta4 implements Einschrittverfahren {

    @Override
    /**
     * {@inheritDoc}
     * Bei der Umsetzung koennen die Methoden addVectors und multScalar benutzt werden.
     */
    public double[] nextStep(double[] y_k, double t, double delta_t, ODE ode) {
        double[] T1 = multScalar(ode.auswerten(t, y_k), delta_t);
        double[] T2 = multScalar(ode.auswerten(t + (delta_t/2), addVectors(y_k, multScalar(T1, (delta_t / 2)))), delta_t);
        double[] T3 = multScalar(ode.auswerten(t + (delta_t/2), addVectors(y_k, multScalar(T2, (delta_t / 2)))), delta_t);
        double[] T4 = multScalar(ode.auswerten(t + delta_t, addVectors(y_k, multScalar(T3, (delta_t)))), delta_t);

        return addVectors(y_k, multScalar(addVector4(T1, T2, T3, T4), 1./6.));
    }

    /**
     * addiert die zwei Vektoren a und b
     */
    private double[] addVectors(double[] a, double[] b) {
        double[] erg = new double[a.length];
        for (int i = 0; i < a.length; i++) {
            erg[i] = a[i] + b[i];
        }
        return erg;
    }

    private double[] addVector4(double[] a, double[] b, double[] c, double[] d) {
        return addVectors(addVectors(addVectors(a, b), c), d);
    }

    /**
     * multipliziert den Skalar scalar auf den Vektor a
     */
    private double[] multScalar(double[] a, double scalar) {
        double[] erg = new double[a.length];
        for (int i = 0; i < a.length; i++) {
            erg[i] = scalar * a[i];
        }
        return erg;
    }

}

