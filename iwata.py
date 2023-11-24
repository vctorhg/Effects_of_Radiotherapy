

# IMPORTING LIBRARIES
import numpy as np


# RHO'S INITIAL CONDITION.
def rho_at_x_t0(x, t):
    # INITIAL CONDITION FOR RHO AS SEEN IN IWATA.
    #
    # INPUT:
    #   x : NUMPY ARRAY. CONTAINS PRIMARY TUMOR SIZES AS FLOATS.
    #   t : FLOAT. DISCRETE TIME-STEP.
    #
    # OUTPUT:
    #   rho : NUMPY ARRAY. CONTAINS DENSITY VALUES.

    rho = np.zeros(len(x))

    return rho

# GOMPERTZ GROWTH FUNCTION
def g_G(x, a_G, b_G):
    # COMPUTES  GROWTH  USING  GOMPERTZ  GROWTH RATE AS SEEN IN
    # IWATA.
    #
    # INPUT:
    #   x : NUMPY ARRAY. CONTAINS PRIMARY TUMOR SIZES AS FLOATS
    #   a_G : FLOAT. GROWTH RATE CONSTANT.
    #   b_G : FLOAT. TUMOR SIZE AT SATURATED LEVEL.
    #
    # OUTPUT:
    #   growth : NUMPY ARRAY. COMPUTED GROWTH VALUES.

    growth = a_G * x * np.log(b_G/x)

    return growth

# GOMPERTZ GROWTH DERIVATIVE
def g_x_G(x, a_G, b_G):
    # COMPUTES GOMPERTZ GROWTH DERIVATIVE.
    #
    # INPUT:
    #   x : NUMPY ARRAY. CONTAINS PRIMARY TUMOR SIZES AS FLOATS
    #   a_G : FLOAT. GROWTH RATE CONSTANT.
    #   b_G : FLOAT. TUMOR SIZE AT SATURATED LEVEL.
    #
    # OUTPUT:
    #   growth_der : NUMPY  ARRAY.  COMPUTED  GROWTH DERIVATIVE
    #                VALUES.

    growth_der = (a_G*np.log(b_G/x))-a_G

    return growth_der

# ANALYTICAL SOLUTION OF THE GOMPERTZ GROWTH FUNCTION FOR COMPUTING MAX X.
def x_max_G(x_init, t, a_G, b_G):
    # COMPUTES x_M FOR THE GOMPERTZ GROWTH FUNCTION.
    #
    # INPUT:
    #   x_init : NUMPY ARRAY.  CONTAINS  PRIMARY TUMOR SIZES AS
    #            FLOATS.
    #   t : FLOAT. DISCRETE TIME-STEP.
    #   a_G : FLOAT. GROWTH RATE CONSTANT.
    #   b_G : FLOAT. TUMOR SIZE AT SATURATED LEVEL.
    #
    # OUTPUT:
    #   maximum_x : NUMPY ARRAY.  CONTAINS  MAXIMUM TUMOR SIZES
    #               AT TIME 't'.

    maximum_x = (b_G**(1-np.exp(-a_G*t)))*(x_init**(np.exp(-a_G*t)))

    return maximum_x

# NOTE: BETA FUNCTION FROM IWATA IS ALPHA FUNCTION FROM ANGULO.
def beta(x, m_G, alpha_G):
    # COMPUTES COLONIZATION RATE AS SEEN IN IWATA.
    #
    # INPUT:
    #   x : NUMPY ARRAY. CONTAINS PRIMARY TUMOR SIZES AS FLOATS.
    #   m_G : FLOAT. COLONIZATION COEFFICIENT.
    #   alpha_G : FLOAT. FRACTAL DIMENSION OF BLOOD VESSELS.
    #
    # OUTPUT:
    #   vals : NUMPY ARRAY. COMPUTED SEEDING VALUES.

    vals = m_G*(x**alpha_G)

    return vals
