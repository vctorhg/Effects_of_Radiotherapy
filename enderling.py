

# IMPORTING LIBRARIES
import numpy as np


# LINEAR QUADRATIC GROWTH FUNCTION
def g_LQ(x, t, D_LQ, alpha_LQ, beta_LQ):
    # COMPUTES GROWTH USING THE LINEAR-QUADRATIC  FUNCTION SEEN
    # IN ENDERLING
    #
    # INPUT:
    #   x : NUMPY ARRAY. CONTAINS PRIMARY TUMOR SIZES AS FLOATS
    #   t : FLOAT. DISCRETE TIME-STEP.
    #   D_LQ : FLOAT. DOSAGE.
    #   alpha_LQ : FLOAT. SINGLE HIT PARAMETER.
    #   beta_LQ : FLOAT. MULTIPLE HIT PARAMETER.
    #
    # OUTPUT:
    #   growth : NUMPY ARRAY. COMPUTED GROWTH VALUES.

    # COMPUTING GROWTH
    growth = -((alpha_LQ * D_LQ) + (2 * beta_LQ * (D_LQ ** 2) * t)) * x

    return growth

# LINEAR QUADRATIC GROWTH DERIVATIVE
def g_x_LQ(t, D_LQ, alpha_LQ, beta_LQ):
    # COMPUTES LINEAR QUADRATIC GROWTH DERIVATIVE
    #
    # INPUT:
    #   t : FLOAT. DISCRETE TIME-STEP.
    #   D_LQ : FLOAT. DOSAGE.
    #   alpha_LQ : FLOAT. SINGLE HIT PARAMETER.
    #   beta_LQ : FLOAT. MULTIPLE HIT PARAMETER.
    #
    # OUTPUT:
    #   growth_der : NUMPY ARRAY. COMPUTED GROWTH DERIVATIVE VALUES.

    # COMPUTING GROWTH DERIVATIVE
    growth_der = -((alpha_LQ*D_LQ) + (2*beta_LQ*(D_LQ**2)*t))

    return growth_der

# ANALYTICAL SOLUTION OF THE LINEAR QUADRATIC GROWTH FUNCTION FOR COMPUTING MAX X.
def x_max_LQ(x_init, t, D_LQ, alpha_LQ, beta_LQ):
    # COMPUTES x_M FOR THE LINEAR QUADRATIC GROWTH FUNCTION.
    #
    # INPUT:
    #   x_init : NUMPY ARRAY. INITIAL VALUES FOR 'x'.
    #   t : FLOAT. DISCRETE TIME-STEP.
    #   D_LQ : FLOAT. DOSAGE.
    #   alpha_LQ : FLOAT. SINGLE HIT PARAMETER.
    #   beta_LQ : FLOAT. MULTIPLE HIT PARAMETER.
    #
    # OUTPUT:
    #   maximum_x : NUMPY ARRAY. COMPUTED MAXIMUM 'x' VALUES FOR THE GIVEN 't'.

    # COMPUTING MAXIMUM PRIMARY TUMOR SIZE
    maximum_x = x_init * np.exp(-(alpha_LQ*D_LQ*t)-(beta_LQ*(D_LQ**2)*(t**2)))

    return maximum_x
