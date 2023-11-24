

# IMPORTING LIBRARIES
import pandas as pd
import seaborn.objects as so

# IMPORTING FILES
import tools


def plotAllCurves(x_list, y_list, y_task, plot_type, x_axis_label, y_axis_label, plot_title,
                  curve_label_list=None, export_plot_name=None):
    # CALLS FUNCTION  'setDataInPlottingFormat'  TO EXTRACT AND
    # FORMAT INFORMATION IN 'x_list'  AND  'y_list'.  THEN USES
    # 'seaborn' LIBRARY TO PRODUCE THE NECESSARY PLOT ACCORDING
    # TO ARGUMENT  'y_category'.ALL OTHER ARGUMENTS SPECIFY THE
    # PLOT'S LABELS.
    #
    # INPUT:
    #   x_list : LIST. CONTAINS LISTS.EACH HOLDING COMPUTED 'x'
    #            VALUES.
    #   y_list : LIST. CONTAINS LISTS.EACH HOLDING THE CORRESP.
    #            'y' VALUES TO THOSE IN 'x_list'.
    #   y_task : STRING. SPECIFIES  ACTION  TO  BE PERFORMED ON
    #            THE DATA CONTAINED IN 'y_list'.POSSIBLE VALUES
    #            ARE 'normalize' AND 'append'.
    #   plot_type : STRING. INDICATES WHICH PLOT CONFIGURATION
    #               TO USE. POSSIBLE VALUES ARE:
    #               - 'density'   FOR  IWATA'S  DENSITY  PLOTS.
    #               - 'number_of_tumors' FOR IWATAS ACCUMULATED
    #                                    NUMBER OF TUMORS.
    #               - 'max_tumor_size'  FOR MAXIMUM TUMOR SIZE.
    #   x_axis_label : STRING. LABEL FOR THE 'x' AXIS.
    #   y_axis_label : STRING. LABEL FOR THE 'y' AXIS.
    #   plot_title : STRING. TITLE FOR THE PLOT.
    #   curve_label_list : LIST. CONTAINS STRINGS. EACH BEING A
    #                      LABEL  FOR  THE CURVES PRODUCED WHEN
    #                      USING  THE  INFORMATION  IN 'x_list'
    #                      AND 'y_list',THUS LENGTH MUST BE THE
    #                      EQUAL  AS THAT OF THE AFOREMENTIONED
    #                      LISTS.
    #   export_plot_name : STRING. FILENAME  /WITHOUT FILE-TYPE
    #                      APPENDIX/  UNDER WHICH THE PLOT WILL
    #                      BE STORED.
    #
    # HARD-CODED VARIABLES:
    #   export_path : STRING. PATH  USED  TO  SAVE  THE PLOT AS
    #                 A '.png'  FILE. EXPORTS FILE TO  'Images'
    #                 DIRECTORY.
    #
    # OUTPUT:
    #   complete_plot : SEABORN OBJECT.PLOT SHOWING ALL CURVES.
    #   export_plot_name.png : PNG FILE.STORED 'complete_plot'.

    # SETTING EXPORT PATH
    if not export_plot_name:
        export_plot_name = plot_type
    export_path = '/Users/victor/Documents/TUM/Thesis/Images/' + export_plot_name + '.png'

    # SETTING NUMERICAL VALUES IN A 'seaborn.objects' FRIENDLY FORMAT:
    plot_data_DF = tools.setDataInPlottingFormat(x_data=x_list, y_data=y_list,
                                                 y_data_task=y_task, label_data=curve_label_list)
    # PLOTTING
    if plot_type == 'no_treatment_density':
        complete_plot = (
            so.Plot(data=plot_data_DF, x='x', y='y', color='curve')
                .add(so.Lines())
                .scale(y='log')
                .limit(x=(1e-1, 8e10), y=(1e-18, 1e-8))
                .label(title=plot_title, x=x_axis_label, y=y_axis_label, color='')
                .layout(size=(16, 9))
                .save(export_path, bbox_inches='tight')
        )

    elif plot_type == 'continuous_therapy_density':
        complete_plot = (
            so.Plot(data=plot_data_DF, x='x', y='y', color='curve')
            .add(so.Lines())
            .scale(y='log')
            .limit(x=(1e-1, 8e10), y=(1e-18, 1e-8))
            .label(title=plot_title, x=x_axis_label, y=y_axis_label, color='')
            .layout(size=(16, 9))
            .save(export_path, bbox_inches='tight')
        )

    if plot_type == 'interval_therapy_density':
        complete_plot = (
            so.Plot(data=plot_data_DF, x='x', y='y', color='curve')
                .add(so.Lines())
                .scale(y='log')
                .limit(x=(1e-1, 8e10), y=(1e-20, 1e-8))
                .label(title=plot_title, x=x_axis_label, y=y_axis_label, color='')
                .layout(size=(16, 9))
                .save(export_path, bbox_inches='tight')
        )

    elif plot_type == 'number_of_tumors':
        complete_plot = (
            so.Plot(data=plot_data_DF, x='x', y='y', color='curve')
                .add(so.Lines())
                .scale(x='log', y='log')
                .label(title=plot_title, x=x_axis_label, y=y_axis_label, color='')
                .layout(size=(16, 9))
                .save(export_path, bbox_inches='tight')
        )

    elif plot_type == 'max_tumor_size':
        complete_plot = (
            so.Plot(data=plot_data_DF, x='x', y='y', color='curve')
                .add(so.Lines())
                .label(title=plot_title, x=x_axis_label, y=y_axis_label, color='')
                .layout(size=(16, 9))
                .save(export_path, bbox_inches='tight')
        )

    complete_plot.show()


def plotConcavityChange(derivative_dataframe, export_plot_name=None):
    # MELTS  COLUMNS  IN  THE  DATAFRAME  RETURNED  BY FUNCTION
    # 'findTimeWhenConcavityChanges' SO THAT THE 'seaborn' LIB.
    # PRODUCES THE CORRESPONDING PLOTS.IT THEN CALLS THE AFORE-
    # MENTIONED LIBRARY TO PRODUCE AND STORE THE PLOT.
    #
    # INPUT:
    #  derivative_dataframe : PANDAS DATAFRAME.CONSISTS OF COLS
    #                         'TIME',   'MIN DERIVATIVE',   AND
    #                         'MAX DERIVATIVE'.DATAFRAME IS THE
    #                         OUTPUT FROM FUNCION
    #                         'findTimeWhenConcavityChanges'.
    #   export_plot_name : STRING. FILENAME  /WITHOUT FILE-TYPE
    #                      APPENDIX/  UNDER WHICH THE PLOT WILL
    #                      BE STORED.
    #
    # HARD-CODED VARIABLES:
    #   export_path : STRING. PATH  USED  TO  SAVE  THE PLOT AS
    #                 A '.png'  FILE. EXPORTS FILE TO  'Images'
    #                 DIRECTORY.
    #
    # OUTPUT:
    #   concavity_plot : SEABORN OBJECT.PLOT SHOWING CURVES FOR
    #                    MAXIMUM AND MINIMUM DERIVATIVE VALUES.
    #   export_plot_name.png : PNG FILE.STORED 'concavity_plot'

    # SETTING EXPORT PATH
    if not export_plot_name:
        export_plot_name = 'Concavity'
    export_path = '/Users/victor/Documents/TUM/Thesis/Images/' + export_plot_name + '.png'

    # MELTING DATAFRAME
    plot_data_DF = pd.melt(frame=derivative_dataframe, id_vars=['TIME'],
                           value_vars=['MIN DERIVATIVE', 'MAX DERIVATIVE'],
                           var_name='CURVE', value_name='DERIVATIVE VALUE', ignore_index=True)
    # PLOTTING
    concavity_plot = (
        so.Plot(data=plot_data_DF, x='TIME', y='DERIVATIVE VALUE', color='CURVE')
           .add(so.Lines())
           .label(title='DERIVATIVE VALUE EVOLUTION THROUGH TIME',
                  x='TIME', y='NUMERICAL VALUE', color='')
           .limit(y=(-5, 5))
           .layout(size=(16, 9))
           .save(export_path, bbox_inches='tight')
    )

    concavity_plot.show()




