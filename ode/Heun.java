package ode;

/**
 * Das Einschrittverfahren von Heun
 *
 * @author braeckle
 *
 */
public class Heun implements Einschrittverfahren {

    @Override
    /**
     * {@inheritDoc}
     * Nutzen Sie dabei geschickt den Expliziten Euler.
     */

    public double[] nextStep(double[] y_k, double t, double delta_t, ODE ode) {
        ExpliziterEuler euler = new ExpliziterEuler();

        double delta_t_half = delta_t / 2.0;
        double[] y_k_f = ode.auswerten(t, y_k);

        double[] euler_step = euler.nextStep(y_k, t, delta_t, ode);
        double[] y_k_f_f = ode.auswerten(delta_t + t, euler_step);

        for(int i = 0; i < y_k.length; i++) {
            y_k[i] = y_k[i] + (delta_t_half * (y_k_f[i] + y_k_f_f[i]));
        }
        return y_k;
    }

}

