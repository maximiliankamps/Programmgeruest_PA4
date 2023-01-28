package ode;

/**
 * Das Einschrittverfahren "Expliziter Euler"
 *
 * @author braeckle
 *
 */
public class ExpliziterEuler implements Einschrittverfahren {

    public double[] nextStep(double[] y_k, double t, double delta_t, ODE ode) { //Works
        double[] y_k_cpy = ode.auswerten(t, y_k);
        for(int i = 0; i < y_k.length; i++) {
            y_k[i] = y_k[i] + delta_t * y_k_cpy[i];
        }
        return y_k;
    }

}

