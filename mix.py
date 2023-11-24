

# IMPORTING FILES
import iwata
import enderling


def g_GLQ(t_prev, x_prev, a_G, b_G, D_LQ, alpha_LQ, beta_LQ):
    # COMBINES THE GOMPERTZ GROWTH FUNCTION FROM IWATA WITH THE
    # LIN-QUAD GROWTH FUNCTION FROM ENDERLING.
    #
    # INPUT:
    #   t_prev,...,beta_LQ :  PARAMETERS   USED  BY  IWATA  AND
    #                        ENDERLING     GROWTH    FUNCTIONS.
    #
    # OUTPUT:
    #   growth :  NUMPY  ARRAY.  GROWTH  COMPUTED  FROM  MIXING
    #             IWATA AND ENDERLING GROWTH RATES.

    # COMPUTING GROWTH
    growth = iwata.g_G(x_prev, a_G, b_G) + enderling.g_LQ(x_prev, t_prev, D_LQ, alpha_LQ, beta_LQ)

    return growth

def g_x_GLQ(t, x, a_G, b_G, D_LQ, alpha_LQ, beta_LQ):
    # COMBINES   GROWTH  DERIVATIVE  FUNCTION  FROM  IWATA  AND
    # ENDERLING.
    #
    # INPUT:
    #   t,...,beta_LQ :  PARAMETERS USED BY IWATA AND ENDERLING
    #   GROWTH DERIVATIVE FUNCTIONS.
    #
    # OUTPUT:
    #   growth_der : NUMPY ARRAY.  COMBINED  GROWTH DERIVATIVE.

    # COMPUTING COMBINED GROWTH DERIVATIVE
    growth_der = iwata.g_x_G(x, a_G, b_G) + enderling.g_x_LQ(t, D_LQ, alpha_LQ, beta_LQ)

    return growth_der

def x_max_GLQ(x_init, t, a_G, b_G, D_LQ, alpha_LQ, beta_LQ):
    # COMPUTES MAXIMUM TUMOR SIZE WHEN THE MAXIMUM TUMOR SIZES
    # FROM IWATA AND ENDERLING ARE COMBINED.
    #
    # INPUT:
    #   x_init,...,beta_LQ : PARAMETERS  USED  BY MAXIMUM TUMOR
    #                        SIZE   IN   IWATA  AND  ENDERLING.
    #
    # OUTPUT:
    #   maximum_x : NUMPY ARRAY. MAXIMUM TUMOR SIZE.

    # COMPUTING MAXIMUM TUMOR SIZE
    maximum_x = iwata.x_max_G(x_init, t, a_G, b_G) + enderling.x_max_LQ(x_init, t, D_LQ, alpha_LQ, beta_LQ)

    return maximum_x


