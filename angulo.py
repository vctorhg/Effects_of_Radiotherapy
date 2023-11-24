
# IMPORTING LIBRARIES
import numpy as np
import pandas as pd

import mix
# IMPORTING FILES
import tools
import iwata
import enderling

# ======================================================== GOMPERTZ ===================================================
# INTERMEDIATE TIME-STEP
def intermediateTimeStep_G(x_prev, a_G, b_G, k, x_init):
    # MODIFIED   VERSION  OF  THE  INTERMEDIATE  TIME-STEP  FOR
    # GOMPERTZ TUMOR GROWTH FUNCTION.
    # COMPUTES NEW TIME-LEVEL VALUE FROM GIVEN 't' AND 'k'.
    # COMPUTES NEW 'x' VALUES FROM GIVEN 'x'.
    # COMPUTES NEW DENSITY VALUES FROM GIVEN 'rho' AND FRESH 'x'
    #
    # INPUT:
    #   k : FLOAT. STEP VALUE.
    #   previous_time : FLOAT. DISCRETE TIME-STEP.
    #   x_prev : NUMPY ARRAY. CONTAINS PRIMARY TUMOR SIZES AS FLOATS.
    #   rho_prev : NUMPY ARRAY. DENSITY VALUES.
    #
    # OUTPUT:
    #   x_new : NUMPY ARRAY. COMPUTED TUMOR SIZE.

    # INITIALIZING 'X' VECTOR
    x_new = np.zeros(len(x_prev)+1)
    # COMPUTING 'X' VECTOR ENTRIES
    x_new[0] = x_init
    x_new[1:] = x_prev + ((k/2)*iwata.g_G(x_prev, a_G, b_G))

    return x_new

# FULL TIME-STEP
def fullTimeStep_G(previous_time, x_prev, x_inter, rho_prev, a_G, b_G, k, m_G, alpha_G, x_init):
    # FULL   TIME-STEP  COMPUTATIONS ARE CARRIED OUT AS SEEN IN
    # ANGULO, WITH THE EXCEPTION THAT AN EXTRA TERM IS ADDED TO
    # ACCOUNT FOR THE SEEDING SEEN IN IWATA.
    # COMPUTES  NEW TIME, THEN NEX X-VALUES AND FINALLY DENSITY
    # VALUES.
    #
    # INPUT:
    #   previous_time: FLOAT. DISCRETE TIME-STEP FROM PREVIOUS
    #                  COMPUTATIONS.
    #   x_prev : NUMPY ARRAY. CONTAINS  PRIMARY  TUMOR SIZES AS
    #            FLOATS.
    #   x_inter : NUMPY ARRAY. INTERMEDIATE TIME-STEP x VALUES.
    #   rho_prev : ARRAY. PREVIOUSLY COMPUTED DENSITY VALUES.
    #   a_G : FLOAT. CONSTANT TUMOR GROWTH RATE.
    #   b_G : FLOAT. TUMOR SIZE AT SATURATION LEVEL.
    #   k : FLOAT. STEP VALUE.
    #   m_G, alpha_G : FLOATS. USED IN COLONY RATE COMPUTATION.
    #   x_init : FOLLOWING IWATA, THE INITIAL TUMOR SIZE IS SET
    #            TO '1'.
    #
    # OUTPUT:
    #   t_new : FLOAT. DISCRETE TIME STEP.
    #   x_new : ARRAY. FRESH PRIMARY TUMOR SIZE.
    #   rho_new : ARRAY: FRESH DENSITY VALUES.

    # COMPUTING NEW TIME
    t_new = previous_time + k
    # INITIALIZING 'X' VECTOR
    x_new = np.zeros(len(x_inter))
    # COMPUTING 'X' VECTOR ENTRIES
    x_new[0] = x_init
    x_new[1:] = x_prev + (k*iwata.g_G(x_inter[1:], a_G, b_G))

    # INITIALIZING 'RHO' VECTOR
    rho_new = np.zeros(len(x_inter))
    # COMPUTING 'RHO' VECTOR ENTRIES
    rho_new[1:] = np.multiply(rho_prev, np.exp(-k*iwata.g_x_G(x_inter[1:], a_G, b_G)))
    beta_vals = iwata.beta(x_new, m_G, alpha_G)
    x_new_jmin1, x_new_j = tools.splitShift(x_new,1)
    betaRho = np.multiply(beta_vals[1:], rho_new[1:])
    betaRho_jmin1, betaRho_j = tools.splitShift(betaRho, 1)
    coeff = 2/((2*iwata.g_G(x_new[0], a_G, b_G))-((x_new[1]-x_new[0])*beta_vals[0]))
    first_term = ((x_new[1]-x_new[0])/2)*(beta_vals[1]*rho_new[1])
    sum_term = np.sum(np.multiply(((x_new_j[1:]-x_new_jmin1[1:])/2), (betaRho_jmin1+betaRho_j)))
    additonal_seed = iwata.beta(iwata.x_max_G(x_init, t_new, a_G, b_G), m_G, alpha_G)
    # INCORPORATING IWATA'S SEEDING TERM
    rho_new[0] = coeff * (first_term + sum_term + additonal_seed)
    return t_new, x_new, rho_new

# ITERATION
def iterate_G(number_of_iterations, initial_time, maximum_time, a_G, b_G, k, m_G, alpha_G, x_initial_condition,
              export_filename=None, initialize=True, x_initial=None, rho_initial=None):
    # ITERATES OVER CALLS TO FUNCTIONS  'intermediateTimeStep',
    # AND  'fullTimeStep' AS SEEN IN THE ANGULO PAPERS. IN EACH
    # ITERATION IT COMPUTES NUMERICAL APPROXIMATIONS FOR 'x(t)'
    # /TUMOR SIZE IN CELLS/,  AND 'rho(x,t)' /THE CORRESPONDING
    # DENSITY/. UNLIKE IN THE ANGULO PAPERS, NO DIMENSIONAL RE-
    # DUCTION IS PERFORMED i.e: THE ARRAYS' DIMENSION INCREASES
    # BY  1  EACH ITERATION. THE MORE ITERATIONS, THE FINER THE
    # FINAL MESH ON THE 'x(t)' ARRAY /AND THEREFORE MORE EVALU-
    # ATIONS OF 'rho'/.
    #
    # INPUT:
    #   number_of_iterations : INTEGER. SPECIFIES THE NUMBER OF
    #                          ITERATIONS TO BE PERFORMED.
    #   initial_time ,..., x_initial_condition : ARGUMENTS USED
    #                 BY FUNCTIONS  'intermediateTimeStep', AND
    #                 'fullTimeStep'
    #   export_filename : STRING.  DEFAULTS TO 'None'. FILENAME
    #                     /WITHOUT  FILE-TYPE  APPENDIX/  UNDER
    #                     WHICH  THE  DATAFRAME WILL BE STORED.
    #   initialize : BOOLEAN.  DEFAULTS  TO TRUE.  IF TRUE, 'x'
    #                AND 'rho' VALUES ARE INITIALIZED. IF FALSE
    #                OPTIONAL  ARGUMENTS  'x_initial' AND 'rho_
    #                initial'  ARE  USED  TO  INITIALIZE CORRE-
    #                SPONDING ARRAYS.
    #   x_initial : NUMPY ARRAY.  CONTAINS  TUMOR  SIZES.  ONLY
    #               USED  WHEN  'initialize'  IS  SET TO FALSE.
    #   rho_initial : NUMPY  ARRAY.  CONTAINS   DENSITY VALUES.
    #                 ONLY USED  WHEN  'initialize'  IS  SET TO
    #                 FALSE.
    #
    # HARD-CODED VARIABLES:
    #   export_path : STRING.PATH USED TO SAVE THE DATAFRAME AS
    #                 A '.csv'  FILE. EXPORTS FILE TO  'Output'
    #                 DIRECTORY.
    #
    # OUTPUT:
    #   x_prev : NUMPY ARRAY. CONTAINS COMPUTED TUMOR SIZE VALS.
    #   rho_prev : NUMPY ARRAY.  CONTAINS COMPUTED DENSITY VALS.
    #   export_filename.csv : CSV FILE. STORED
    #                         'numerical_approximation_DF'.

    if initialize:
        # INITIALIZING NUMPY ARRAYS
        t = np.linspace(initial_time, maximum_time, number_of_iterations+1)
        x_initial = iwata.x_max_G(x_initial_condition, t, a_G, b_G)
        rho_initial = iwata.rho_at_x_t0(x_initial, t)

    # INITIALIZING VARIABLES
    t_prev = initial_time
    x_prev = x_initial
    rho_prev = rho_initial

    # ITERATING
    for iteration in range(number_of_iterations):
        # COMPUTING INTERMEDIATE TIME-STEP
        x_inter = intermediateTimeStep_G(x_prev, a_G, b_G, k, x_initial_condition)
        # COMPUTING FULL TIME-STEP
        t_prev, x_prev, rho_prev = fullTimeStep_G(t_prev, x_prev, x_inter, rho_prev, a_G, b_G, k, m_G,
                                                         alpha_G, x_initial_condition)

    # EXPORTING
    if export_filename:
        numerical_approximation_DF = pd.DataFrame({'x': x_prev, 'rho': rho_prev})
        export_path = '/Users/victor/Documents/TUM/Thesis/Output/' + export_filename + '.csv'
        numerical_approximation_DF.to_csv(export_path, index=False)

    return x_prev, rho_prev


# =================================================== GOMPERTZ & LIN-QUAD =============================================
# INTERMEDIATE TIME-STEP
def intermediateTimeStep_GLQ(t_prev, x_prev, a_G, b_G, k, D_LQ, alpha_LQ, beta_LQ, x_init):
    # COMPUTES   MODIFIED  ANGULO  INTERMEDIATE TIME-STEP FOR A
    # GROWTH  RATE  EQUAL TO GOMPERTZ-LIN_QUAD. THIS IS USED TO
    # MODEL CONTINUOUS EPSILON-DOSE RADIOTHERAPY.
    #
    # INPUT:
    #   t_prev : FLOAT. PREVIOUSLY COMPUTED TIME.
    #   x_prev : NUMPY ARRAY.CONTAINS PREVIOUSLY COMPUTED TUMOR
    #            SIZES.
    #   a_G, b_G : FLOATS. PARAMETERS USED FOR GOMPERTZ GROWTH.
    #   k : INTEGER. TIME STEP FOR ANGULO-ALGORITHM ITERATIONS.
    #   D_LQ : FLOAT.  DOSAGE USED IN RADIOTHERAPY. FOR THE GLQ
    #          MODEL  /EPSILON-DOSAGE/ THE VALUES ARE TYPICALLY
    #          SMALL.
    #   alpha_LQ, beta_LQ : FLOATS. PARAMETERS USED IN LIN_QUAD
    #   x_init : INTEGER. INITIAL TUMOR SIZE.
    #
    # OUTPUT:
    #   x_new : NUMPY ARRAY. INTERMEDIATE TIME-STEP TUMOR SIZES

    # INITIALIZING 'X' VECTOR.
    x_new = np.zeros(len(x_prev)+1)
    # COMPUTING 'X' VECTOR ENTRIES
    x_new[0] = x_init
    x_new[1:] = x_prev + ((k / 2) * mix.g_GLQ(t_prev, x_prev, a_G, b_G, D_LQ, alpha_LQ, beta_LQ))

    return x_new

# FULL TIME-STEP
def fullTimeStep_GLQ(t_prev, x_prev, x_inter, rho_prev, a_G, b_G, k, m_G, D_LQ, alpha_G, alpha_LQ, beta_LQ, x_init):
    # COMPUTES FULL TIME-STEP VECTOR VALUES ACCORDING TO ANGULO
    # NUMERICAL ALGORITHM.
    #
    # INPUT:
    #   t_prev : FLOAT. PREVIOUSLY COMPUTED TIME.
    #   x_prev : NUMPY ARRAY.CONTAINS PREVIOUSLY COMPUTED TUMOR
    #            SIZES.
    #   x_inter : NUMPY ARRAY. INTERMEDIATE TIME-STEP x VALUES.
    #   rho_prev : ARRAY.  PREVIOUSLY  COMPUTED DENSITY VALUES.
    #   a_G, b_G, m_G : FLOATS.  PARAMETERS  USED FOR GOMPERTZ.
    #   k : INTEGER. TIME STEP FOR ANGULO-ALGORITHM ITERATIONS.
    #   D_LQ : FLOAT.  DOSAGE USED IN RADIOTHERAPY. FOR THE GLQ
    #          MODEL  /EPSILON-DOSAGE/ THE VALUES ARE TYPICALLY
    #          SMALL.
    #   alpha_LQ, beta_LQ : FLOATS. PARAMETERS USED IN LIN_QUAD
    #   x_init : INTEGER. INITIAL TUMOR SIZE.
    #
    # OUTPUT:
    #   t_new : FLOAT. DISCRETE TIME STEP.
    #   x_new : ARRAY. FRESH PRIMARY TUMOR SIZE.
    #   rho_new : ARRAY: FRESH DENSITY VALUES.

    # COMPUTING NEW TIME
    t_new = t_prev + k
    # INITIALIZING 'X' VECTOR
    x_new = np.zeros(len(x_inter))
    # COMPUTING 'X' VECTOR ENTRIES
    x_new[0] = x_init
    x_new[1:] = x_prev + (k * mix.g_GLQ(t_prev, x_inter[1:], a_G, b_G, D_LQ, alpha_LQ, beta_LQ))

    # INITIALIZING 'RHO' VECTOR
    rho_new = np.zeros(len(x_inter))
    # COMPUTING 'RHO' VECTOR ENTRIES
    rho_new[1:] = np.multiply(rho_prev, np.exp(-k * mix.g_x_GLQ(t_new, x_inter[1:], a_G, b_G, D_LQ, alpha_LQ, beta_LQ)))
    beta_vals = iwata.beta(x_new, m_G, alpha_G)
    x_new_jmin1, x_new_j = tools.splitShift(x_new, 1)
    betaRho = np.multiply(beta_vals[1:], rho_new[1:])
    betaRho_jmin1, betaRho_j = tools.splitShift(betaRho, 1)
    coeff = 2/((2*mix.g_GLQ(t_new, x_new[0], a_G, b_G, D_LQ, alpha_LQ, beta_LQ))-((x_new[1]-x_new[0])*beta_vals[0]))
    first_term = ((x_new[1]-x_new[0])/2)*(beta_vals[1]*rho_new[1])
    sum_term = np.sum(np.multiply(((x_new_j[1:]-x_new_jmin1[1:])/2), (betaRho_jmin1+betaRho_j)))
    additonal_seed = iwata.beta(mix.x_max_GLQ(x_init, t_new, a_G, b_G, D_LQ, alpha_LQ, beta_LQ), m_G, alpha_G)
    # ACCOUNTING FOR ADDITIONAL SEEDING TERM SEEN IN IWATA
    rho_new[0] = coeff * (first_term + sum_term + additonal_seed)

    return t_new, x_new, rho_new

# ITERATION
def iterate_GLQ(number_of_iterations,initial_time, maximum_time, a_G, b_G, k, m_G, D_LQ, alpha_G, alpha_LQ, beta_LQ,
                x_initial_condition, export_filename=None, initialize=True, x_initial=None, rho_initial=None):
    # ITERATES OVER CALLS TO FUNCTIONS  'intermediateTimeStep',
    # AND  'fullTimeStep' AS SEEN IN THE ANGULO PAPERS. IN EACH
    # ITERATION IT COMPUTES NUMERICAL APPROXIMATIONS FOR 'x(t)'
    # /TUMOR SIZE IN CELLS/,  AND 'rho(x,t)' /THE CORRESPONDING
    # DENSITY/. UNLIKE IN THE ANGULO PAPERS, NO DIMENSIONAL RE-
    # DUCTION IS PERFORMED i.e: THE ARRAYS' DIMENSION INCREASES
    # BY  1  EACH ITERATION. THE MORE ITERATIONS, THE FINER THE
    # FINAL MESH ON THE 'x(t)' ARRAY /AND THEREFORE MORE EVALU-
    # ATIONS OF 'rho'/.
    #
    # INPUT:
    #   number_of_iterations : INTEGER. SPECIFIES THE NUMBER OF
    #                          ITERATIONS   TO   BE  PERFORMED.
    #   initial_time ,..., x_initial_condition : ARGUMENTS USED
    #                 BY FUNCTIONS  'intermediateTimeStep', AND
    #                 'fullTimeStep'
    #   export_filename : STRING.  DEFAULTS TO 'None'. FILENAME
    #                     /WITHOUT  FILE-TYPE  APPENDIX/  UNDER
    #                     WHICH  THE  DATAFRAME WILL BE STORED.
    #   initialize : BOOLEAN. DEFUALTS TO TRUE. IF SET TO FALSE
    #                USES     PARAMETERS     'x_initial'    AND
    #                'rho_initial'.
    #
    # HARD-CODED VARIABLES:
    #   export_path : STRING.PATH USED TO SAVE THE DATAFRAME AS
    #                 A '.csv'  FILE. EXPORTS FILE TO  'Output'
    #                 DIRECTORY.
    #
    # OUTPUT:
    #   x_prev : NUMPY ARRAY. CONTAINS COMPUTED TUMOR SIZE VALS.
    #   rho_prev : NUMPY ARRAY.  CONTAINS COMPUTED DENSITY VALS.
    #   export_filename.csv : CSV FILE. STORED
    #                         'numerical_approximation_DF'.

    if initialize:
        # INITIALIZING NUMPY ARRAYS
        t = np.linspace(initial_time, maximum_time, number_of_iterations+1)
        x_initial = mix.x_max_GLQ(x_initial_condition, t, a_G, b_G, D_LQ, alpha_LQ, beta_LQ)
        rho_initial = iwata.rho_at_x_t0(x_initial, t)

    # INITIALIZING VARIABLES
    t_prev = initial_time
    x_prev = x_initial
    rho_prev = rho_initial

    # ITERATING
    for iteration in range(number_of_iterations):
        # COMPUTING INTERMEDIATE TIME-STEP
        x_inter = intermediateTimeStep_GLQ(t_prev, x_prev, a_G, b_G, k, D_LQ, alpha_LQ, beta_LQ, x_initial_condition)
        # COMPUTING FULL TIME-STEP
        t_prev, x_prev, rho_prev = fullTimeStep_GLQ(t_prev, x_prev, x_inter, rho_prev, a_G, b_G, k, m_G,
                                                    D_LQ, alpha_G, alpha_LQ, beta_LQ, x_initial_condition)

    # EXPORTING
    if export_filename:
        numerical_approximation_DF = pd.DataFrame({'x': x_prev, 'rho': rho_prev})
        export_path = '/Users/victor/Documents/TUM/Thesis/Output/' + export_filename + '.csv'
        numerical_approximation_DF.to_csv(export_path, index=False)

    return x_prev, rho_prev


# ======================================================== LIN-QUAD ===================================================
# INTERMEDIATE TIME STEP
def intermediateTimeStep_LQ(t_prev, x_prev, k, D_LQ, alpha_LQ, beta_LQ, x_init):
    # COMPUTES  ANGULO  INTERMEDIATE TIME-STEP FOR THE LIN-QUAD
    # GROWTH  FUNCTION.  THIS  MODEL  IS  USED DURING CLINICAL-
    # STANDARD  THERAPY  SESSIONS  TO ACCOUNT FOR PRIMARY-TUMOR
    # SIZE DECREASE.
    #
    # INPUT:
    #   t_prev : FLOAT. PREVIOUSLY COMPUTED TIME.
    #   x_prev : NUMPY ARRAY.CONTAINS PREVIOUSLY COMPUTED TUMOR
    #            SIZES.
    #   k : INTEGER. TIME STEP FOR ANGULO-ALGORITHM ITERATIONS.
    #   D_LQ : FLOAT.  DOSAGE USED IN RADIOTHERAPY. FOR THE GLQ
    #          MODEL  /EPSILON-DOSAGE/ THE VALUES ARE TYPICALLY
    #          SMALL.
    #   alpha_LQ, beta_LQ : FLOATS. PARAMETERS USED IN LIN_QUAD
    #   x_init : INTEGER. INITIAL TUMOR SIZE.
    #
    # OUTPUT:
    #   x_new : NUMPY ARRAY. COMPUTED PRIMARY TUMOR SIZE.

    # INITIALIZING 'X' VECTOR
    x_new = np.zeros(len(x_prev)+1)
    # COMPUTING 'X' VECTOR ENTRIES.
    x_new[0] = x_init
    x_new[1:] = x_prev + ((k / 2) * enderling.g_LQ(x_prev, t_prev, D_LQ, alpha_LQ, beta_LQ))

    return x_new

# FULL TIME-STEP
def fullTimeStep_LQ(t_prev, x_prev, x_inter, rho_prev, k, m_G, D_LQ, alpha_G, alpha_LQ, beta_LQ, x_init):
    # COMPUTES    ANGULO   FULL   TIME-STEP   FOR  THE LIN-QUAD
    # GROWTH  FUNCTION.  THIS  MODEL  IS  USED DURING CLINICAL-
    # STANDARD  THERAPY  SESSIONS  TO ACCOUNT FOR PRIMARY-TUMOR
    # SIZE DECREASE.
    #
    # INPUT:
    #   t_prev : FLOAT. PREVIOUSLY COMPUTED TIME.
    #   x_prev : NUMPY ARRAY.CONTAINS PREVIOUSLY COMPUTED TUMOR
    #            SIZES.
    #   a_G, b_G : FLOATS. PARAMETERS USED FOR GOMPERTZ GROWTH.
    #   k : INTEGER. TIME STEP FOR ANGULO-ALGORITHM ITERATIONS.
    #   D_LQ : FLOAT.  DOSAGE USED IN RADIOTHERAPY. FOR THE GLQ
    #          MODEL  /EPSILON-DOSAGE/ THE VALUES ARE TYPICALLY
    #          SMALL.
    #   alpha_LQ, beta_LQ : FLOATS. PARAMETERS USED IN LIN_QUAD
    #   x_init : INTEGER. INITIAL TUMOR SIZE.
    #
    # OUTPUT:
    #   t_new : FLOAT. DISCRETE TIME STEP.
    #   x_new : ARRAY. FRESH PRIMARY TUMOR SIZE.
    #   rho_new : ARRAY: FRESH DENSITY VALUES.

    # COMPUTING NEW TIME
    t_new = t_prev + k
    # INITIALIZING 'X' VECTOR
    x_new = np.zeros(len(x_inter))
    # COMPUTING 'X' VECTOR ENTRIES
    x_new[0] = x_init
    x_new[1:] = x_prev + (k * enderling.g_LQ(x_inter[1:], t_prev, D_LQ, alpha_LQ, beta_LQ))

    # INITIALIZING 'RHO' VECTOR
    rho_new = np.zeros(len(x_inter))
    # COMPUTING 'RHO' VECTOR ENTRIES
    rho_new[1:] = np.multiply(rho_prev, np.exp(-k * enderling.g_x_LQ(t_new, D_LQ, alpha_LQ, beta_LQ)))
    beta_vals = iwata.beta(x_new, m_G, alpha_G)
    x_new_jmin1, x_new_j = tools.splitShift(x_new, 1)
    betaRho = np.multiply(beta_vals[1:], rho_new[1:])
    betaRho_jmin1, betaRho_j = tools.splitShift(betaRho, 1)
    coeff = 2/((2*enderling.g_LQ(x_new[0], t_new, D_LQ, alpha_LQ, beta_LQ))-((x_new[1]-x_new[0])*beta_vals[0]))
    first_term = ((x_new[1]-x_new[0])/2)*(beta_vals[1]*rho_new[1])
    sum_term = np.sum(np.multiply(((x_new_j[1:]-x_new_jmin1[1:])/2), (betaRho_jmin1+betaRho_j)))
    additonal_seed = iwata.beta(enderling.x_max_LQ(x_init, t_new, D_LQ, alpha_LQ, beta_LQ), m_G, alpha_G)
    # ACCOUNTING FOR ADDITIONAL SEEDING TERM SEEN IN IWATA.
    rho_new[0] = coeff * (first_term + sum_term + additonal_seed)

    return t_new, x_new, rho_new

# ITERATION
def iterate_LQ(number_of_iterations, initial_time, maximum_time, k, m_G, D_LQ, alpha_G, alpha_LQ, beta_LQ,
               x_initial_condition, export_filename=None):
    # ITERATES OVER CALLS TO FUNCTIONS  'intermediateTimeStep',
    # AND  'fullTimeStep' AS SEEN IN THE ANGULO PAPERS. IN EACH
    # ITERATION IT COMPUTES NUMERICAL APPROXIMATIONS FOR 'x(t)'
    # /TUMOR SIZE IN CELLS/,  AND 'rho(x,t)' /THE CORRESPONDING
    # DENSITY/. UNLIKE IN THE ANGULO PAPERS, NO DIMENSIONAL RE-
    # DUCTION IS PERFORMED i.e: THE ARRAYS' DIMENSION INCREASES
    # BY  1  EACH ITERATION. THE MORE ITERATIONS, THE FINER THE
    # FINAL MESH ON THE 'x(t)' ARRAY /AND THEREFORE MORE EVALU-
    # ATIONS OF 'rho'/.
    #
    # INPUT:
    #   number_of_iterations : INTEGER. SPECIFIES THE NUMBER OF
    #                          ITERATIONS TO BE PERFORMED.
    #   initial_time ,..., x_initial_condition : ARGUMENTS USED
    #                 BY FUNCTIONS  'intermediateTimeStep', AND
    #                 'fullTimeStep'
    #   export_filename : STRING.  DEFAULTS TO 'None'. FILENAME
    #                     /WITHOUT  FILE-TYPE  APPENDIX/  UNDER
    #                     WHICH  THE  DATAFRAME WILL BE STORED.
    #   initialize : BOOLEAN.  DEFAULTS  TO TRUE.  IF TRUE, 'x'
    #                AND 'rho' VALUES ARE INITIALIZED. IF FALSE
    #                OPTIONAL  ARGUMENTS  'x_initial' AND 'rho_
    #                initial'  ARE  USED  TO  INITIALIZE CORRE-
    #                SPONDING ARRAYS.
    #   x_initial : NUMPY ARRAY.  CONTAINS  TUMOR  SIZES.  ONLY
    #               USED  WHEN  'initialize'  IS  SET TO FALSE.
    #   rho_initial : NUMPY  ARRAY.  CONTAINS   DENSITY VALUES.
    #                 ONLY USED  WHEN  'initialize'  IS  SET TO
    #                 FALSE.
    #
    # HARD-CODED VARIABLES:
    #   export_path : STRING.PATH USED TO SAVE THE DATAFRAME AS
    #                 A '.csv'  FILE. EXPORTS FILE TO  'Output'
    #                 DIRECTORY.
    #
    # OUTPUT:
    #   x_prev : NUMPY ARRAY. CONTAINS COMPUTED TUMOR SIZE VALS.
    #   rho_prev : NUMPY ARRAY.  CONTAINS COMPUTED DENSITY VALS.
    #   export_filename.csv : CSV FILE. STORED
    #                         'numerical_approximation_DF'.

    # INITIALIZING NUMPY ARRAYS
    t = np.linspace(initial_time, maximum_time, number_of_iterations+1)
    x_initial = enderling.x_max_LQ(x_initial_condition, t, D_LQ, alpha_LQ, beta_LQ)
    rho_initial = iwata.rho_at_x_t0(x_initial, t)

    # INITIALIZING VARIABLES
    t_prev = initial_time
    x_prev = x_initial
    rho_prev = rho_initial

    # ITERATING
    for iteration in range(number_of_iterations):
        # COMPUTING INTERMEDIATE TIME-STEP
        x_inter = intermediateTimeStep_LQ(t_prev, x_prev, k, D_LQ, alpha_LQ, beta_LQ, x_initial_condition)
        # COMPUTING FULL TIME-STEP
        t_prev, x_prev, rho_prev = fullTimeStep_LQ(t_prev, x_prev, x_inter, rho_prev, k, m_G, D_LQ, alpha_G,
                                                    alpha_LQ, beta_LQ, x_initial_condition)

    # EXPORTING
    if export_filename:
        numerical_approximation_DF = pd.DataFrame({'x': x_prev, 'rho': rho_prev})
        export_path = '/Users/victor/Documents/TUM/Thesis/Output/' + export_filename + '.csv'
        numerical_approximation_DF.to_csv(export_path, index=False)

    return x_prev, rho_prev


# ====================================================== RADIOTHERAPY =================================================
# TUMOR DENSITY /RADIOTHERAPY/
def iterate_G_with_Radiotherapy(therapy_type, time_at_therapy_start, therapy_days, rest_days, therapy_sessions,
                                initial_time, total_time, a_G, b_G, k, m_G, D_LQ, alpha_G, alpha_LQ, beta_LQ,
                                x_initial_condition, export_filename_list=[None, None, None]):
    # COMPUTES TUMOR DENSITY WHEN THERAPY IS ADMINISTERED.  THE
    # THERAPY  CAN  BE ADMINISTERED OVER A CONTINUOUS PERIOD OF
    # TIME  USING  MICRO-DOSAGE, OR INTENSELY OVER A SHORT, AND
    # DEFINED PERIOD OF TIME.PERFORMS TRUNCATION OF TUMOR SIZES
    # AND CORRESPONDING DENSITIES WHEN THEIR SIZE IS BELOW '1'.
    # THIS  IS  DONE TO ENSURE THAT ERADICATED METASTASIS DON'T
    # SEED NEW TUMORS.
    #
    # INPUT:
    #   therapy_type : STRING. DETERMINES THERAPY TYPE TO MODEL
    #                  POSSIBLES  VALUES ARE:  'continuous' AND
    #                  'interval'.
    #   time_at_therapy_start : INTEGER. NUMBER OF DAY IN WHICH
    #                           RADIOTHERAPY  WAS  RECEIVED FOR
    #                           THE FIRST TIME.
    #   therapy_days : INTEGER.NUMBER  OF UNINTERRUPTED DAYS IN
    #                  WHICH THE PATIENT RECEIVES RADIOTHERAPY.
    #   rest_days : INTEGER.  NUMBER  OF  UNINTERRUPTED DAYS IN
    #               WHICH  THE PATIENT RESTS FROM RADIOTHERAPY.
    #   therapy_sessions : INTEGER. NUMBER OF TIMES THE PATIENT
    #                      WILL   RECEIVE   RADIOTHERAPY.  EACH
    #                      SESSION HAS  A TOTAL LENGTH EQUAL TO
    #                      'therapy_days' + 'rest_days'.
    #   initial_time : INTEGER. DAY ON WHICH SIMULATION BEGINS.
    #   total_time : INTEGER. TOTAL NUMBER OF DAYS SPANNED BY
    #                THE SIMULATION.
    #   a_G,...,x_initial_condition : PARAMETERS USED IN TUMOR-
    #                                 COMPUTATIONS.
    #   export_filename_list : LIST.  ENTRIES USED AS FOLLOWS :
    #                         -[0] : TUMOR GROWTH UNTIL THERAPY
    #                         -[1] : TUMOR GROWTH  THERAPY REST
    #                         -[2] : TUMOR GROWTH AFTER THERAPY
    #
    # OUTPUT:
    #   x_G : NUMPY ARRAY. COMPUTED TUMOR SIZES.
    #   rho_G : NUMPY ARRAY. COMPUTED DENSITY VALUES.

    # COMPUTING 'x' AND 'rho' VALUES USING GOMPERTZ GROWTH RATE UNTIL 'time_at_therapy_start'
    iterations = round(time_at_therapy_start/k)
    x_G, rho_G = iterate_G(number_of_iterations=iterations, initial_time=initial_time,
                           maximum_time=time_at_therapy_start, a_G=a_G, b_G=b_G, k=k, m_G=m_G, alpha_G=alpha_G,
                           x_initial_condition=x_initial_condition, export_filename=export_filename_list[0])

    # CLINICALLY ESTABLISHED RADIOTHERAPY PROCEDURE
    if therapy_type == 'interval':
        # COMPUTING RADIOTHERAPY EFFECTS
        for session in range(therapy_sessions):
            # COMPUTING TUMOR DECREASE DUE TO RADIOTHERAPY
            radiotherapy_effect = enderling.x_max_LQ(x_init=x_G, t=therapy_days, D_LQ=D_LQ,
                                                     alpha_LQ=alpha_LQ, beta_LQ=beta_LQ)
            # DELETING METASTASIS BELOW THRESHOLD
            treated_x_G, treated_rho_G = tools.deleteMetastasis(x_array=radiotherapy_effect, rho_array=rho_G,
                                                                threshold=1)
            # COMPUTING TUMOR INCREASE DURING RADIOTHERAPY REST
            iterations = round(rest_days/k)
            t_start = time_at_therapy_start + (therapy_days*(session+1)) + (rest_days*session)
            x_G, rho_G = iterate_G(number_of_iterations=iterations, initial_time=t_start,
                                   maximum_time=t_start + rest_days, a_G=a_G, b_G=b_G, k=k, m_G=m_G, alpha_G=alpha_G,
                                   x_initial_condition=x_initial_condition, export_filename=export_filename_list[1],
                                   initialize=False, x_initial=treated_x_G, rho_initial=treated_rho_G)

        # COMPUTING TUMOR INCREASE FOR THE REMAINING DAYS AFTER THERAPY
        time_elapsed = time_at_therapy_start + ((therapy_days+rest_days)*therapy_sessions)
        remaining_time = total_time - time_elapsed
        if remaining_time >= 1:
            iterations = round(remaining_time / k)
            x_G, rho_G = iterate_G(number_of_iterations=iterations, initial_time=time_elapsed,
                                   maximum_time=remaining_time, a_G=a_G, b_G=b_G, k=k, m_G=m_G, alpha_G=alpha_G,
                                   x_initial_condition=x_initial_condition, export_filename=export_filename_list[2],
                                   initialize=False, x_initial=x_G, rho_initial=rho_G)

    # CONTINUOUS RADIOTHERAPY MICRODOSAGE
    elif therapy_type == 'continuous':
        remaining_time = total_time - time_at_therapy_start
        iterations = round(remaining_time/k)
        x_G, rho_G = iterate_GLQ(number_of_iterations=iterations, initial_time=time_at_therapy_start,
                                 maximum_time=remaining_time, a_G=a_G, b_G=b_G, k=k, m_G=m_G, D_LQ=D_LQ,
                                 alpha_G=alpha_G, alpha_LQ=alpha_LQ, beta_LQ=beta_LQ,
                                 x_initial_condition=x_initial_condition, export_filename=export_filename_list[1],
                                 initialize=False, x_initial=x_G, rho_initial=rho_G)

    return x_G, rho_G
