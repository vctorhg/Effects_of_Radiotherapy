

# IMPORTING LIBRARIES
import numpy as np
import pandas as pd

# IMPORTING FILES
import iwata
import enderling


# PRIMARY TUMOR SIZE /NO RADIOTHERAPY/
def primaryTumor(initial_time, total_time, a_G, b_G, x_initial_condition, export_filename=None):
    # MODELS PRIMARY TUMOR GROWTH WHEN NO TREATMENT IS PROVIDED
    # FOR  THIS IT USES THE GOMPERTZ GROWTH FUNCTION AS SEEN IN
    # IWATA.
    #
    # INPUT:
    #   initial_time : INTEGER.  PRIMARY  TUMOR  GROWTH  DAY 0.
    #   total_time : INTEGER.  PRIMARY  TUMOR GROWTH FINAL DAY.
    #   a_G :  FLOAT.  PRIMARY  TUMOR  GROWTH  RATE   CONSTANT.
    #   b_G : FLOAT.  PRIMARY  TUMOR  SIZE AT SATURATED GROWTH.
    #   x_initial_condition : PRIMARY  TUMOR  INITIAL  SIZE  IN
    #                         CELLS. FOLLOWING IWATA, IT IS SET
    #                         TO 1.
    #   export_filename : STRING. NAME UNDER WHICH COMPUTATIONS
    #                     WILL BE STORED IN A CSV FILE.
    #
    # OUTPUT:
    #   time : NUMPY ARRAY.  CONTAINS  INTEGERS /DAYS/; ONE FOR
    #          EACH DAY BETWEEN /AND INCLUDING/  'initial_time'
    #          AND 'total_time'.
    #   untreated_tumor_size : NUMPY ARRAY.   CONTAINS  FLOATS.
    #                          THAT  CORRESPOND  TO THE PRIMARY
    #                          TUMOR  SIZE  AT EACH OF THE DAYS
    #                          IN THE 'time' ARRAY.

    # COMPUTING TIME AND CORRESPONDING PRIMARY TUMOR SIZE
    time = np.array(range(initial_time, total_time+1))
    untreated_tumor_size = iwata.x_max_G(x_init=x_initial_condition, t=time, a_G=a_G, b_G=b_G)

    # EXPORTING COMPUTATIONS
    if export_filename:
        max_tumor_DF = pd.DataFrame({'t': time, 'x(t) UNTREATED': untreated_tumor_size})
        export_path = '/Users/victor/Documents/TUM/Thesis/Output/' + export_filename + '.csv'
        max_tumor_DF.to_csv(export_path, index=False)

    return time, untreated_tumor_size

# PRIMARY TUMOR SIZE /RADIOTHERAPY/
def primaryTumorWithRadiotherapy(therapy_type, time_at_therapy_start, therapy_days, rest_days, therapy_sessions,
                                 initial_time, total_time, a_G, b_G, D_LQ, alpha_LQ, beta_LQ,
                                 x_initial_condition, export_filename=None):
    # COMPUTES  PRIMARY TUMOR SIZE WHEN PATIENT RECEIVES RADIO-
    # THERAPY.  PARAMETER 'therapy_type' DETERMINES WHETHER THE
    # THERAPY  IS  ADMINISTERED  IN  HIGH DOSES OVER A SPECIFIC
    # PERIOD OF TIME OR IN DAILY MICRO DOSES UNTIL 'total_time'
    #
    # INPUT:
    #   therapy_type : STRING. POSSIBLE VALUES ARE:
    #                  - 'interval' : WHEN  THERAPY  IS APPLIED
    #                                 FOR  'therapy_days'  DAYS
    #                                 FOLLOWED  BY  A  BREAK OF
    #                                 'rest_days' DAYS.STARTING
    #                                 ON 'time_at_therapy_start'
    #                                 THE  CYCLE  CONSISTING OF
    #                                 THERAPY+REST  IS REPEATED
    #                                 'therapy_sessions' TIMES.
    #                  - 'continuous' : WHEN THERAPY IS APPLIED
    #                                   FROM 'time_at_therapy_-
    #                                   start' TO 'total_time'.
    #   time_at_therapy_start : INTEGER. DAY ON WHICH THE FIRST
    #                           DOSE  OF RADIOTHERAPY IS ADMIN-
    #                           ISTERED.
    #   therapy_days :  INTEGER.  UNINTERRUPTED  THERAPY  DAYS.
    #   rest_days: INTEGER.  UNINTERRUPTED  REST DAYS FOLLOWING
    #              'therapy_days'.
    #   therapy_sessions : NUMBER OF TIMES THE 'therapy_days' +
    #                      'rest_days' CYCLE IS REPEATED.
    #   initial_time : INTEGER.  DAY  ON  WHICH  THE SIMULATION
    #                  BEGINS
    #   total_time : TOTAL DURATION OF THE SIMULATION.
    #   a_g,...,x_initial_condition : PARAMETERS.
    #   export_filename : STRING. FILENAME UNDER WHICH COMPUTED
    #                     DATA WILL BE STORED AS A '.csv' FILE.

    # INITIALIZING LISTS
    therapy_time_intervals = []
    therapy_tumor_sizes = []
    untreated_time_intervals = []
    untreated_tumor_sizes = []

    # COMPUTING TUMOR SIZES BEFORE THERAPY BEGINS
    time_until_therapy = np.array(range(initial_time, time_at_therapy_start+1))
    tumor_until_therapy = iwata.x_max_G(x_init=x_initial_condition, t=time_until_therapy, a_G=a_G, b_G=b_G)
    x_initial_condition = tumor_until_therapy[-1]
    # STORING COMPUTATIONS ACCORDINGLY
    untreated_time_intervals.extend(time_until_therapy)
    untreated_tumor_sizes.extend(tumor_until_therapy)


    # COMPUTING TUMOR SIZES FOR INTERVAL THERAPY
    if therapy_type == 'interval':
        # SETTING TIME INTERVALS TO BE USED THROUGHOUT ITERATIONS
        therapy_time = np.array(range(initial_time, therapy_days+1))
        rest_time = np.array(range(initial_time, rest_days+1))
        # COMPUTING TUMOR SIZES FOR EACH THERAPY SESSION
        for session in range(therapy_sessions):
            treated_tumor_size = enderling.x_max_LQ(x_init=x_initial_condition, t=therapy_time, D_LQ=D_LQ,
                                                    alpha_LQ=alpha_LQ, beta_LQ=beta_LQ)
            x_initial_condition = treated_tumor_size[-1]
            rest_tumor_size = iwata.x_max_G(x_init=x_initial_condition, t=rest_time, a_G=a_G, b_G=b_G)
            x_initial_condition = rest_tumor_size[-1]
            # STORING COMPUTATIONS ACCORDINGLY
            therapy_time_intervals.extend(therapy_time+time_at_therapy_start+((therapy_days+rest_days)*session))
            therapy_tumor_sizes.extend(treated_tumor_size)
            untreated_time_intervals.extend(rest_time+time_at_therapy_start+
                                            (therapy_days*(session+1))+(rest_days*session))
            untreated_tumor_sizes.extend(rest_tumor_size)

        # COMPUTING TUMOR SIZES AFTER THERAPY
        remaining_time = np.array(range(initial_time,
                                        total_time-
                                        (time_at_therapy_start+((therapy_days+rest_days)*therapy_sessions))+1))
        after_therapy_tumor_size = iwata.x_max_G(x_init=x_initial_condition, t=remaining_time, a_G=a_G, b_G=b_G)
        # STORING COMPUTATIONS ACCORDINGLY
        untreated_time_intervals.extend(remaining_time+
                                        time_at_therapy_start+((therapy_days+rest_days)*therapy_sessions))
        untreated_tumor_sizes.extend(after_therapy_tumor_size)

    # COMPUTING TUMOR SIZES FOR CONTINUOUS THERAPY
    elif therapy_type == 'continuous':
        remaining_time = np.array(range(initial_time, total_time-time_at_therapy_start+1))
        treated_tumor_size = enderling.x_max_LQ(x_init=x_initial_condition, t=remaining_time,
                                                D_LQ=D_LQ, alpha_LQ=alpha_LQ, beta_LQ=beta_LQ)
        # STORING COMPUTATIONS ACCORDINGLY
        therapy_tumor_sizes.extend(treated_tumor_size)
        therapy_time_intervals.extend(remaining_time+time_at_therapy_start)

    # EXPORTING
    if export_filename:
        max_tumor_DF = pd.DataFrame({'UNTREATED TIME': untreated_time_intervals,
                                     'x(t) UNTREATED': untreated_tumor_sizes,
                                     'THERAPY TIME': therapy_time_intervals,
                                     'x(t) THERAPY': therapy_tumor_sizes})
        export_path = '/Users/victor/Documents/TUM/Thesis/Output/' + export_filename + '.csv'
        max_tumor_DF.to_csv(export_path, index=False)


    return untreated_time_intervals, untreated_tumor_sizes, therapy_time_intervals, therapy_tumor_sizes