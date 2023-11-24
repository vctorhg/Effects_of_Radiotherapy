

# IMPORTING LIBRARIES
import numpy as np
import pandas as pd
from scipy import integrate

# IMPORTING FILES
import angulo

# SPLITTING ARRAY
def splitShift(x, shift):
    # TAKES NUMPY ARRAY OR LIST AS ARGUMENT,PERFORMS TWO SPLITS
    #   - FROM: 0 TO: ARRAY LENGTH - SHIFT
    #   - FROM: SHIFT TO: ARRAY LENGTH
    #
    # INPUT:
    #   x : NUMPY ARRAY OR LIST.
    #   shift : INTEGER.SHIFT TO BE USED WHN SPLITTING ARAY 'x'
    #
    # OUTPUT:
    #   start_at_zero : NUMPY ARRAY. FIRST SPLIT.
    #   start_at_shift : NUMPY ARRAY. SECOND SPLIT.

    # FIRST SPLIT
    start_at_zero = x[0:len(x)-shift]
    # SECOND SPLIT
    start_at_shift = x[shift:]

    return start_at_zero, start_at_shift

# FUNCTION TO COMPUTE THE TIME IT TAKES TO REACH A GIVEN PRIMARY TUMOR SIZE
def timeToGivenSize_G(primary_tumor_size, a_G, b_G):
    # BASED  ON  PIRMIN'S  ARTICLE:  'HOW MATHEMATICAL MODELING
    # COULD  CONTRIBUTE  TO  THE  QUANTIFICATION  OF METASTATIC
    # TUMOR BURDEN UNDER THERAPY: INSIGHTS IN IMMUNOTHERAPEUTIC
    # TREATMENT  OF  NON-SMALL CELL LUNG CANCER', USES EQ(3) OF
    # THE PAPER TO COMPUTE THE TIME IT WOULD TAKE FOR A PRIMARY
    # TUMOR  TO  REACH  THE  GIVEN  SIZE /'primary_tumor_size'/
    # CONSIDERING  THE  ANALYTICAL  SOLUTION  TO  THE  GOMPERTZ
    # GROWTH RATE AND AN INITIAL CONDITION OF x_0 = 1.
    #
    # INPUT:
    #   primary_tumor_size : INTEGER.  SIZE EXPRESSED IN NUMBER
    #                        OF   CELLS.   SCIENTIFIC  NOTATION
    #                        ACCEPTED.
    #   a_parameter_G : FLOAT. GROWTH RATE CONSTANT.
    #   b_parameter_G : FLOAT. TUMOR SIZE AT SATURATED LEVEL.
    #
    # OUTPUT:
    #   t : INTEGER. NUMBER OF DAYS IT WOULD TAKE TO REACH  THE
    #       GIVEN PRIMARY TUMOR SIZE AND GROWTH PARAMETERS.

    # COMPUTING TIME
    t = -(1/a_G)*np.log(1-(np.log(primary_tumor_size)/np.log(b_G)))

    return round(t)

# SETS DATA IN A SEABORN FRIDNLY FORMAT
def setDataInPlottingFormat(x_data, y_data, y_data_task, label_data=None):
    # SETS RECEIVED DATA IN A 'seaborn.object' FRIENDLY FORMAT.
    # NAMELY, EXTRACTS  ELEMENTS  FROM  'x_data'  AND 'y_data',
    # /'x_data'  AND  'y_data'  ARE  LISTS  WHOSE  ELEMENTS ARE
    # LISTS THEMSELVES/  AND,  CONSIDERING   'y_data_category',
    # CONCATENATES  THEM  ACCORDINGLY.  IT THEN SETS THE CORRE-
    # SPONDING  LABEL  IN  A SEPARATE LIST,  AND USES ALL THREE
    # LISTS  TO  CREATE  A  DATAFRAME THAT THE 'seaborn'
    # LIBRARY CAN HANDLE.
    #
    # INPUT:
    #   x_data : LIST. ELEMENTS ARE LISTS THAT CONTAIN COMPUTED
    #            'x'  VALUES  ACCORDING TO THE ANGULO ALGORITHM
    #            EACH LIST MAY HAVE A DIFFERENT LENGTH FROM THE
    #            OTHER.
    #   y_data : LIST.  ELEMENTS  ARE  LISTS  THAT  CONTAIN THE
    #            COMP.  'rho'  VALUES  ACCORDING  TO THE ANGULO
    #            ALGORITHM.  EACH  LIST  MAY  HAVE  A DIFFERENT
    #            LENGTH  FROM  THE  OTHER,  BUT  MUST  HAVE THE
    #            SAME  LENGTH  AS ITS CORRESPONDING  'x'  VALUE
    #            LIST.
    #   y_data_task : STRING. SPECIFIES  ACTION TO BE PERFORMED
    #                 ON THE DATA CONTAINED IN 'y_list'.POSSIBLE
    #                 VALUES ARE 'normalize' AND 'append'.
    #   label_data : LIST.DEFAULTS TO 'None'.ELSE, ELEMENTS ARE
    #                STRINGS  /POSSIBLY  LATEX  FORMAT/. NUMBER
    #                OF cELEMENTS  MUST  BE  EQUAL  TO  THAT OF
    #                'x_data' AND 'rho_data'.
    #
    # VARIABLES:
    #   x : LIST. CONTAINS CONCATENATED 'x' VALUES.
    #   y : LIST. CONTAINS CONCATENATED 'rho' VALUES.
    #   curve : LIST. CONTAINS CORRESPONDING CURVE TAGS.
    #
    # OUTPUT:
    #   plot_data_DF : PANDAS DATAFRAME. CONSISTS OF  'x',  'y'
    #                  AND 'curve' AS COLUMNS.

    # INITIALIZING LISTS THAT WILL TURN INTO DATAFRAME COLUMNS
    x = []
    y = []
    curve = []
    for idx, x_values in enumerate(x_data):
        # APPENDING NUMERICAL VALUES
        x.extend(x_values)
        if y_data_task == 'normalize':
            # APPENDED RHO VALUES MUST BE NORMALIZED.
            y.extend(y_data[idx]/sum(y_data[idx]))
        elif y_data_task == 'append':
            y.extend(y_data[idx])

        # LABELING VALUES ACCORDINGLY
        if label_data:
            # SETTING CUSTOM LABELS
            curve.extend([label_data[idx]]*len(x_data[idx]))
        else:
            # SETTING DEFAULT LABELS
            curve.extend(['y_' + str(idx)] * len(x_data[idx]))

    # CREATING DATAFRAME
    plot_data_DF = pd.DataFrame({'x': x, 'y': y, 'curve': curve})

    return plot_data_DF


def findLargestMetastasis(x_data, rho_data, n_largest=None, export_filename=None):
    # FLIPS  RHO ARRAY BEFORE COMPUTING THE CUMULATIVE SUM, AND
    # FLOORS THE CUMULATIVE SUM ARRAY ELEMENTS. SPLITS ARRAY IN
    # TWO  PARTS WITH 1-SHIFT TO PERFORM FASTER COMPARISON, AND
    # LOOKS FOR NON-ZERO RESULTS.  ITERATES OVER FOUND NON-ZERO
    # INDEXES TO EXTRACT CORRESPONDING TUMOR SIZE AND NUMBER OF
    # TUMORS.
    # - IF  A  SPECIFIC 'n_largest' NUMBER OF TUMORS WAS GIVEN,
    #   IT TRUNCATES THE 'tumor_size' and 'number_of_tumors'
    #   LISTS AT ENTRY NUMBER 'n_largest'
    # - ELSE, LEAVES LISTS AS THEY ARE.
    # RETURNS LISTS 'tumor_size' AND 'number_of_tumors'.
    #
    # INPUT:
    #   x_data : LIST. COMPUTED 'x' VALUES.
    #   rho_data : LIST. COMPUTED 'rho' VALUES.
    #   n_largest : DEFAULT 'None'. INTEGER. USED TO RETURN THE
    #               FIRST  ENTRIES  OF THE CORRESPONDING LISTS.
    #   export_filename : STRING.  DEFAULTS TO 'None'. FILENAME
    #                     /WITHOUT  FILE-TYPE  APPENDIX/  UNDER
    #                     WHICH  THE  DATAFRAME WILL BE STORED.
    #
    # HARD-CODED VARIABLES:
    #   export_path : STRING.PATH USED TO SAVE THE DATAFRAME AS
    #                 A '.csv'  FILE. EXPORTS FILE TO  'Output'
    #                 DIRECTORY.
    #
    # OUTPUT:
    #   tumor_size : LIST. CONTAINS TUMOR SIZES /NO. OF CELLS/.
    #   number_of_tumors : LIST. CONTAINS INTEGERS.
    #   export_filename.csv : CSV FILE.
    #                         STORED 'largest_metatasis_DF'.

    flipped_rho = np.flip(rho_data)
    cumulative_rho = np.floor(np.cumsum(flipped_rho))
    # SPLITTING COMPUTED CUMULATIVE SUM TO OPTIMIZE ENTRY COMPARISON.
    cum_rho_jmin1, cum_rho_j = splitShift(cumulative_rho, 1)
    tumor_count_difference = cum_rho_j - cum_rho_jmin1
    # NON-ZERO ENTRIES IN DIFFERENCE CORRESPOND TO INTEGER JUMPS BETWEEN ENTRIES IN RHO.
    nonzero_tumor_count = list(np.nonzero(tumor_count_difference)[0])
    # FLIPPING 'x_data' SO FOUND INDEXES CORRESPOND.
    flipped_x = list(np.flip(x_data))
    # CREATING OUTPUT LISTS
    tumor_size = [flipped_x[idx+1] for idx in nonzero_tumor_count]
    number_of_tumors = [cumulative_rho[idx+1] for idx in nonzero_tumor_count]

    # TRUNCATING IF NECESSARY
    if n_largest:
        tumor_size = tumor_size[0:n_largest]
        number_of_tumors = number_of_tumors[0:n_largest]

    # EXPORTING
    if export_filename:
        largest_metastasis_DF = pd.DataFrame({'TUMOR SIZE': tumor_size, 'NUMBER OF TUMORS': number_of_tumors})
        export_path = '/Users/victor/Documents/TUM/Thesis/Output/' + export_filename + '.csv'
        largest_metastasis_DF.to_csv(export_path, index=False)

    # ----- TEMPORARY DATAFRAME :: DELETE THIS FOR FINAL VERSION.
    #raw_metastasis_df = pd.DataFrame({'TUMOR SIZE': np.flip(x_data), 'RHO': np.flip(rho_data), 'CUMSUM': np.cumsum(np.flip(rho_data))})
    #print(raw_metastasis_df.to_markdown())
    #print(largest_metastasis_df.to_markdown())
    # ----- TEMPORARY DATAFRAME

    return tumor_size, number_of_tumors

# ITERATING UNTIL CURVE (C) CHANGES CONCAVITY
def findTimeWhenConcavityCahnges(initial_time_c, initial_maximum_time, a_G_c, b_G_c, k_days_c,
                                 m_G_c, alpha_G_c, x_initial_condition_c,
                                 threshold, step_increase, max_cycles=1000, export_filename=None):
    # INITIALIZES  LISTS  USED TO PRODUCE THE OUTPUT DATAFRAME.
    # INITIALIZES  VARIABLES USED DURING THE ITERATIVE PROCESS.
    # PERFORMS  THE  ITERATIVE  PROCESS  AS LONG AS THE MAXIMUM
    # NUMERICAL   DERIVATIVE   IS  BELOW  THE  GIVEN  THRESHOLD
    # /ARGUMENT 'threshold'/  AND  THE MAXIMUM NUMBER OF CYCLES
    # HAS  NOT  BEEN  SURPASSED.  DURING  EACH  CYCLE, FUNCTION
    # 'iterate_G'  IS CALLED TO COMPUTE 'x' AND 'rho' ACCORDING
    # TO THE SPECIFIED PARAMETERS, THEN BUILT-IN NUMPY FUNCTION
    # 'gradient' IS USED TO COMPUTE THE DERIVATIVE OF THE 'rho'
    # VALUES.  MINIMA  AND MAXIMA ARE FOUND AND STORED IN THEIR
    # RESPECTIVE  LISTS ALONG WITH THE CORRESPONDING MAX. TIME.
    #   + IF THE MAXIMUM NUMBER OF CYCLES HAS BEEN SURPASSED, A
    #     MESSAGE, NOTIFYING OF THE EVENT, IS PRINTED.
    #   + IF A 'filename' HAS BEEN PROVIDED, THE PRODUCED DATA-
    #     FRAME IS EXPORTED.
    #
    # INPUT:
    #   initial_time_c,...,x_initial_condition_c:  FUNCTION PA-
    #                    RAMETERS USEB BY FUNCTION 'iterate()'.
    #   threshold : FLOAT.  DERIVATIVE VALUE TO BE ACHIEVED AND
    #               SURPASSED.
    #   step_increase : INTEGER.  NUMBER OF DAYS TO BE ADDED TO
    #                   'iterative_t_max'  FOR  THE NEXT CYCLE.
    #   max_cycles : INTEGER. DEFAULTS TO 1000.  MAXIMUM NUMBER
    #                OF CYCLES TO BE PERFORMED.
    #   export_filename : STRING.  DEFAULTS TO 'None'. FILENAME
    #                     /WITHOUT  FILE-TYPE  APPENDIX/  UNDER
    #                     WHICH  THE  DATAFRAME WILL BE STORED.
    #
    # HARD-CODED VARIABLES:
    #   export_path : STRING.PATH USED TO SAVE THE DATAFRAME AS
    #                 A '.csv'  FILE. EXPORTS FILE TO  'Output'
    #                 DIRECTORY.
    #
    # OUTPUT:
    #   derivative_DF : PANDAS DATAFRAME.  CONTAINS THE MINIMUM
    #                   AND MAXIMUM NUMERICAL DERIVATIVE VALUES
    #                   AS  WELL AS THEIR CORRESPONDING MAXIMUM
    #                   TIME.
    #   export_filename.csv : CSV FILE. STORED 'derivative_DF'.

    # INITIALIZING OUTPUT LISTS
    max_time_list = []
    max_num_derivative_list = []
    min_num_derivative_list = []

    # INITIALIZING CYCLE VARIABLES
    cycles_performed = 0
    max_numerical_derivative = -1
    iterative_t_max = initial_maximum_time - step_increase

    while max_numerical_derivative <= threshold and cycles_performed <= max_cycles:
        # MAXIMUM TIME FOR COMPUTING 'rho' MUST BE ADJUSTED FOR EACH RUN OF THE CYCLE.
        iterative_t_max += step_increase
        cycles_performed += 1
        iterations_c = round(iterative_t_max / k_days_c)
        # COMPUTING 'x' AND 'rho'
        x_c, rho_c = angulo.iterate_G(number_of_iterations=iterations_c, initial_time=initial_time_c,
                                      maximum_time=iterative_t_max, a_G=a_G_c, b_G=b_G_c,
                                      k=k_days_c, m_G=m_G_c, alpha_G=alpha_G_c,
                                      x_initial_condition=x_initial_condition_c)
        # COMPUTING THE NUMERICAL DERIVATIVE
        gradient_array = np.gradient(rho_c, x_c)
        # FINDING AND STORING DERIVATIVE MAXIMA AND MINIMA
        max_numerical_derivative = max(gradient_array)
        min_numerical_derivative = min(gradient_array)
        max_num_derivative_list.append(max_numerical_derivative)
        min_num_derivative_list.append(min_numerical_derivative)
        # APPENDING THE CORRESPONDING MAXIMUM TIME
        max_time_list.append(iterative_t_max)

    # LETTING USER KNOW THAT THE MAXIMUM NUMBER OF CYCLES WAS REACHED
    if cycles_performed >= max_cycles:
        print('REACHED MAXIMUM NUMBER OF CYCLES.')

    # CREATING RETURN DATAFRAME
    derivative_DF = pd.DataFrame({'MIN DERIVATIVE': min_num_derivative_list,
                                  'MAX DERIVATIVE': max_num_derivative_list,
                                  'TIME': max_time_list})

    # ----- TEMPORARY DATAFRAME PRINT
    #print(derivative_DF.to_markdown())
    # ----- TEMPORARY DATAFRAME PRINT

    # EXPORTING DATAFRAME
    if export_filename:
        export_path = '/Users/victor/Documents/TUM/Thesis/Output/' + export_filename + '.csv'
        derivative_DF.to_csv(export_path, index=False)

    return derivative_DF

# DELETING METASTASIS WITH SIZE BELOW '1'
def deleteMetastasis(x_array, rho_array, threshold):
    # METASTASIS WHOSE SIZE DROPS BELOW THE GIVEN THRESHOLD MAY
    # BE DISCARDED.  IF  THE THRESHOLD IS '1', RADIOTHERAPY HAS
    # SUCCESSFULLY  ATTACKED SMALL METASTASIS, AND THEY MUST BE
    # TAKEN  OUT  OF THE 'x_array' TO PREVENT THEM FROM SEEDING
    # IN FUTURE RUNS /IT IS IMPOSSIBLE FOR A DEAD METASTASIS TO
    # SEED/.
    #
    # INPUT:
    #   x_array : NUMPY ARRAY.  CONTAINS TUMOR SIZE VALUES COM-
    #             PUTED USING THE ANGULO ALGORITHM.
    #   rho_array : NUMPY  ARRAY  OR  LIST. CONTAINS THE CORRE-
    #               SPONDING DENSITY VALUES.
    #   threshold : INTEGER /POSSIBLY FLOAT/.  MINIMUM SIZE FOR
    #               A METASTASIS TO BE ALIVE.
    #
    # OUTPUT:
    #   truncated_x : NUMPY ARRAY.  'x_array' TRUNCATED SO THAT
    #                 ITS FIRST ENTRY IS LARGER OR EQUAL TO THE
    #                 GIVEN THRESHOLD.
    #   truncated_rho : NUMPY  ARRAY OR LIST. TRUNCATED DENSITY
    #                   ARRAY  OF SAME LENGTH AS 'truncated_x'.

    # FINDING INDEX TO TRUNCATE
    truncation_idx = np.argmax(x_array >= threshold)
    # TRUNCATING
    truncated_x = x_array[truncation_idx:]
    truncated_rho = rho_array[truncation_idx:]

    return truncated_x, truncated_rho
