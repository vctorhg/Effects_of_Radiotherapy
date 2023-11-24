
# ========================== NUMERICAL APPROXIMATION FOR MODELLING CANCEROUS NEOPLASM GROWTH ==========================
# TODO: ADD INTRODUCTION
# =====================================================================================================================


# IMPORTING LIBRARIES
import numpy as np

# IMPORTING FILES
import tools
import iwata
import angulo
import plotter
import tumorSize


# ===================================================== VARIABLES =====================================================
# DECLARING INITIAL CONDITIONS AS SEEN IN IWATA
rho_0 = 0
x_0 = 1

# DECLARING CONSTANTS AS SEEN IN IWATA
a_parameter_G = 0.00286
a2_parameter_G = 0.0143
b_parameter_G = 7.3e10
b2_parameter_G = 3.65e10
m_parameter_G = 5.3e-8
alpha_parameter_G = 2/3
alpha2_parameter_G = 0.4
alpha3_parameter_G = 0.8

# DECLARING CONSTANTS AS SEEN IN ENDERLING
D1_parameter_LQ = 0.01
D2_parameter_LQ = 2
alpha_parameter_LQ = 0.1
beta_parameter_LQ = 0.01

# DECLARING AUXILIARY VARS
t_0 = 0
t_therapy = 5
t_rest = 2
sessions = 7
t_diagnosis = tools.timeToGivenSize_G(primary_tumor_size=1.55e9, a_G=a_parameter_G, b_G=b_parameter_G)
t_therapy_start = t_diagnosis + 639
t_max = t_therapy_start + ((t_therapy + t_rest)*sessions) + 60
k_days = 1
iterations = round(t_max/k_days)
# =====================================================================================================================


# ================================================ COMPUTING DENSITIES ================================================
# COMPUTING DENSITY CURVES (a)-(e) /NO RADIOTHERAPY/
# CURVE (a)
x_a, rho_a = angulo.iterate_G(number_of_iterations=iterations, initial_time=t_0, maximum_time=t_max,
                              a_G=a_parameter_G, b_G=b_parameter_G, k=k_days, m_G=m_parameter_G,
                              alpha_G=alpha_parameter_G, x_initial_condition=x_0,
                              export_filename='angulo_G_a_'+str(t_max))

# CURVE (b)
x_b, rho_b = angulo.iterate_G(number_of_iterations=iterations, initial_time=t_0, maximum_time=t_max,
                              a_G=a2_parameter_G, b_G=b_parameter_G, k=k_days, m_G=m_parameter_G,
                              alpha_G=alpha_parameter_G, x_initial_condition=x_0,
                              export_filename='angulo_G_b_'+str(t_max))

# CURVE (c)
x_c, rho_c = angulo.iterate_G(number_of_iterations=iterations, initial_time=t_0, maximum_time=t_max,
                              a_G=a_parameter_G, b_G=b_parameter_G, k=k_days, m_G=m_parameter_G,
                              alpha_G=alpha2_parameter_G, x_initial_condition=x_0,
                              export_filename='angulo_G_c_'+str(t_max))

# CURVE (d)
x_d, rho_d = angulo.iterate_G(number_of_iterations=iterations, initial_time=t_0, maximum_time=t_max,
                              a_G=a_parameter_G, b_G=b_parameter_G, k=k_days, m_G=m_parameter_G,
                              alpha_G=alpha3_parameter_G, x_initial_condition=x_0,
                              export_filename='angulo_G_d_'+str(t_max))

# CURVE (e)
x_e, rho_e = angulo.iterate_G(number_of_iterations=iterations, initial_time=t_0, maximum_time=t_max,
                              a_G=a_parameter_G, b_G=b2_parameter_G, k=k_days, m_G=m_parameter_G,
                              alpha_G=alpha_parameter_G, x_initial_condition=x_0,
                              export_filename='angulo_G_e_'+str(t_max))

# PLOTTING ALL CURVES TOGETHER
plotter.plotAllCurves(x_list=[x_a, x_b, x_c, x_d, x_e],
                      y_list=[rho_a, rho_b, rho_c, rho_d, rho_e],
                      y_task='normalize',
                      plot_type='no_treatment_density',
                      x_axis_label='TUMOR SIZE',
                      y_axis_label='DENSITY',
                      plot_title='TUMOR DENSITY ACCORDING TO SIZE FOR ' + str(t_max) + ' DAYS',
                      curve_label_list=['(a)',
                                        '(b): ' + r'$a_{2}$',
                                        '(c): ' + r'$\alpha_{2}$',
                                        '(d): ' + r'$\alpha_{3}$',
                                        '(e): ' + r'$b_{2}$'],
                      export_plot_name='All_G_Density_'+str(t_max))


# COMPUTING DENSITY CURVES (a)-(e) /CONTINUOUS RADIOTHERAPY/
# CURVE (a)
x_a_continuous, rho_a_continuous = angulo.iterate_G_with_Radiotherapy(therapy_type='continuous',
                                                                      time_at_therapy_start=t_therapy_start,
                                                                      therapy_days=t_therapy, rest_days=t_rest,
                                                                      therapy_sessions=sessions, initial_time=t_0,
                                                                      total_time=t_max, a_G=a_parameter_G,
                                                                      b_G=b_parameter_G, k=k_days, m_G=m_parameter_G,
                                                                      D_LQ=D1_parameter_LQ, alpha_G=alpha_parameter_G,
                                                                      alpha_LQ=alpha_parameter_LQ,
                                                                      beta_LQ=beta_parameter_LQ,
                                                                      x_initial_condition=x_0,
                                                                      export_filename_list=[None, None])
# CURVE (b)
x_b_continuous, rho_b_continuous = angulo.iterate_G_with_Radiotherapy(therapy_type='continuous',
                                                                      time_at_therapy_start=t_therapy_start,
                                                                      therapy_days=t_therapy, rest_days=t_rest,
                                                                      therapy_sessions=sessions, initial_time=t_0,
                                                                      total_time=t_max, a_G=a2_parameter_G,
                                                                      b_G=b_parameter_G, k=k_days, m_G=m_parameter_G,
                                                                      D_LQ=D1_parameter_LQ, alpha_G=alpha_parameter_G,
                                                                      alpha_LQ=alpha_parameter_LQ,
                                                                      beta_LQ=beta_parameter_LQ,
                                                                      x_initial_condition=x_0,
                                                                      export_filename_list=[None, None])
# CURVE (c)
x_c_continuous, rho_c_continuous = angulo.iterate_G_with_Radiotherapy(therapy_type='continuous',
                                                                      time_at_therapy_start=t_therapy_start,
                                                                      therapy_days=t_therapy, rest_days=t_rest,
                                                                      therapy_sessions=sessions, initial_time=t_0,
                                                                      total_time=t_max, a_G=a_parameter_G,
                                                                      b_G=b_parameter_G, k=k_days, m_G=m_parameter_G,
                                                                      D_LQ=D1_parameter_LQ,  alpha_G=alpha2_parameter_G,
                                                                      alpha_LQ=alpha_parameter_LQ,
                                                                      beta_LQ=beta_parameter_LQ,
                                                                      x_initial_condition=x_0,
                                                                      export_filename_list=[None, None])
# CURVE (d)
x_d_continuous, rho_d_continuous = angulo.iterate_G_with_Radiotherapy(therapy_type='continuous',
                                                                      time_at_therapy_start=t_therapy_start,
                                                                      therapy_days=t_therapy, rest_days=t_rest,
                                                                      therapy_sessions=sessions, initial_time=t_0,
                                                                      total_time=t_max, a_G=a_parameter_G,
                                                                      b_G=b_parameter_G, k=k_days, m_G=m_parameter_G,
                                                                      D_LQ=D1_parameter_LQ, alpha_G=alpha3_parameter_G,
                                                                      alpha_LQ=alpha_parameter_LQ,
                                                                      beta_LQ=beta_parameter_LQ,
                                                                      x_initial_condition=x_0,
                                                                      export_filename_list=[None, None])
# CURVE (e)
x_e_continuous, rho_e_continuous = angulo.iterate_G_with_Radiotherapy(therapy_type='continuous',
                                                                      time_at_therapy_start=t_therapy_start,
                                                                      therapy_days=t_therapy, rest_days=t_rest,
                                                                      therapy_sessions=sessions, initial_time=t_0,
                                                                      total_time=t_max, a_G=a_parameter_G,
                                                                      b_G=b2_parameter_G, k=k_days, m_G=m_parameter_G,
                                                                      D_LQ=D1_parameter_LQ, alpha_G=alpha_parameter_G,
                                                                      alpha_LQ=alpha_parameter_LQ,
                                                                      beta_LQ=beta_parameter_LQ,
                                                                      x_initial_condition=x_0,
                                                                      export_filename_list=[None, None])


# PLOTTING ALL CURVES TOGETHER
plotter.plotAllCurves(x_list=[x_a_continuous, x_b_continuous, x_c_continuous, x_d_continuous, x_e_continuous],
                      y_list=[rho_a_continuous, rho_b_continuous, rho_c_continuous, rho_d_continuous, rho_e_continuous],
                      y_task='normalize',
                      plot_type='continuous_therapy_density',
                      x_axis_label='TUMOR SIZE',
                      y_axis_label='DENSITY',
                      plot_title='TUMOR DENSITY ACCORDING TO SIZE FOR ' + str(t_max) +
                                 ' DAYS AND ' + r'$D_{T} = $' + str(D1_parameter_LQ*(t_max-t_therapy_start)) +
                                 ' STARTING AT: '+str(t_therapy_start),
                      curve_label_list=['(a)',
                                        '(b): ' + r'$a_{2}$',
                                        '(c): ' + r'$\alpha_{2}$',
                                        '(d): ' + r'$\alpha_{3}$',
                                        '(e): ' + r'$b_{2}$'],
                      export_plot_name='All_Cont_Density_'+str(t_max)+'_DT_'+
                                       str(D1_parameter_LQ*(t_max-t_therapy_start)) +
                                       '_D_'+str(D1_parameter_LQ).split('.')[1])


# COMPUTING DENSITY CURVES (a)-(e) /RADIOTHERAPY IN A SPECIFIC TIME INTERVAL/
# CURVE (a)
x_a_interval, rho_a_interval = angulo.iterate_G_with_Radiotherapy(therapy_type='interval',
                                                                  time_at_therapy_start=t_therapy_start,
                                                                  therapy_days=t_therapy, rest_days=t_rest,
                                                                  therapy_sessions=sessions, initial_time=t_0,
                                                                  total_time=t_max, a_G=a_parameter_G,
                                                                  b_G=b_parameter_G, k=k_days, m_G=m_parameter_G,
                                                                  D_LQ=D2_parameter_LQ,alpha_G=alpha_parameter_G,
                                                                  alpha_LQ=alpha_parameter_LQ,
                                                                  beta_LQ=beta_parameter_LQ, x_initial_condition=x_0,
                                                                  export_filename_list=[None, None, None])
# CURVE (b)
x_b_interval, rho_b_interval = angulo.iterate_G_with_Radiotherapy(therapy_type='interval',
                                                                  time_at_therapy_start=t_therapy_start,
                                                                  therapy_days=t_therapy, rest_days=t_rest,
                                                                  therapy_sessions=sessions, initial_time=t_0,
                                                                  total_time=t_max, a_G=a2_parameter_G,
                                                                  b_G=b_parameter_G, k=k_days, m_G=m_parameter_G,
                                                                  D_LQ=D2_parameter_LQ, alpha_G=alpha_parameter_G,
                                                                  alpha_LQ=alpha_parameter_LQ,
                                                                  beta_LQ=beta_parameter_LQ, x_initial_condition=x_0,
                                                                  export_filename_list=[None, None, None])
# CURVE (c)
x_c_interval, rho_c_interval = angulo.iterate_G_with_Radiotherapy(therapy_type='interval',
                                                                  time_at_therapy_start=t_therapy_start,
                                                                  therapy_days=t_therapy, rest_days=t_rest,
                                                                  therapy_sessions=sessions, initial_time=t_0,
                                                                  total_time=t_max, a_G=a_parameter_G,
                                                                  b_G=b_parameter_G, k=k_days, m_G=m_parameter_G,
                                                                  D_LQ=D2_parameter_LQ, alpha_G=alpha2_parameter_G,
                                                                  alpha_LQ=alpha_parameter_LQ,
                                                                  beta_LQ=beta_parameter_LQ, x_initial_condition=x_0,
                                                                  export_filename_list=[None, None, None])
# CURVE (d)
x_d_interval, rho_d_interval = angulo.iterate_G_with_Radiotherapy(therapy_type='interval',
                                                                  time_at_therapy_start=t_therapy_start,
                                                                  therapy_days=t_therapy, rest_days=t_rest,
                                                                  therapy_sessions=sessions, initial_time=t_0,
                                                                  total_time=t_max, a_G=a_parameter_G,
                                                                  b_G=b_parameter_G, k=k_days, m_G=m_parameter_G,
                                                                  D_LQ=D2_parameter_LQ, alpha_G=alpha3_parameter_G,
                                                                  alpha_LQ=alpha_parameter_LQ,
                                                                  beta_LQ=beta_parameter_LQ, x_initial_condition=x_0,
                                                                  export_filename_list=[None, None, None])
# CURVE (e)
x_e_interval, rho_e_interval = angulo.iterate_G_with_Radiotherapy(therapy_type='interval',
                                                                  time_at_therapy_start=t_therapy_start,
                                                                  therapy_days=t_therapy, rest_days=t_rest,
                                                                  therapy_sessions=sessions, initial_time=t_0,
                                                                  total_time=t_max, a_G=a_parameter_G,
                                                                  b_G=b2_parameter_G, k=k_days, m_G=m_parameter_G,
                                                                  D_LQ=D2_parameter_LQ, alpha_G=alpha_parameter_G,
                                                                  alpha_LQ=alpha_parameter_LQ,
                                                                  beta_LQ=beta_parameter_LQ, x_initial_condition=x_0,
                                                                  export_filename_list=[None, None, None])
# PLOTTING ALL CURVES TOGETHER
plotter.plotAllCurves(x_list=[x_a_interval, x_b_interval, x_c_interval, x_d_interval, x_e_interval],
                      y_list=[rho_a_interval, rho_b_interval, rho_c_interval, rho_d_interval, rho_e_interval],
                      y_task='normalize',
                      plot_type='interval_therapy_density',
                      x_axis_label='TUMOR SIZE',
                      y_axis_label='DENSITY',
                      plot_title='TUMOR DENSITY ACCORDING TO SIZE FOR ' + str(t_max) +
                                 ' DAYS AND '+r'$D_{T}=$'+str(D2_parameter_LQ*sessions*t_therapy)+' FROM ' +
                                 str(t_therapy_start) + ' TO ' + str(t_therapy_start+((t_therapy+t_rest)*sessions)),
                      curve_label_list=['(a)',
                                        '(b): ' + r'$a_{2}$',
                                        '(c): ' + r'$\alpha_{2}$',
                                        '(d): ' + r'$\alpha_{3}$',
                                        '(e): ' + r'$b_{2}$'],
                      export_plot_name='All_Interval_Density_'+str(t_max)+'_T_'+str(t_therapy_start) +
                                       '_DT_'+str(D2_parameter_LQ*sessions*t_therapy)+'_S_'+str(sessions) +
                                       '_D_'+str(D2_parameter_LQ))


# ============================================= COMPUTING LARGEST TUMORS ==============================================
# FINDING LARGEST TUMORS /NO RADIOTHERAPY/
# CURVE (a)
tumor_size_a, tumor_numbers_a = tools.findLargestMetastasis(x_data=x_a, rho_data=rho_a, export_filename='largest_a')
# CURVE (b)
tumor_size_b, tumor_numbers_b = tools.findLargestMetastasis(x_data=x_b, rho_data=rho_b, export_filename='largest_b')
# CURVE (c)
tumor_size_c, tumor_numbers_c = tools.findLargestMetastasis(x_data=x_c, rho_data=rho_c, export_filename='largest_c')
# CURVE (d)
tumor_size_d, tumor_numbers_d = tools.findLargestMetastasis(x_data=x_d, rho_data=rho_d, export_filename='largest_d')
# CURVE (e)
tumor_size_e, tumor_numbers_e = tools.findLargestMetastasis(x_data=x_e, rho_data=rho_e, export_filename='largest_e')

# PLOTTING
plotter.plotAllCurves(x_list=[tumor_size_a, tumor_size_b, tumor_size_c, tumor_size_d, tumor_size_e],
                      y_list=[tumor_numbers_a, tumor_numbers_b, tumor_numbers_c, tumor_numbers_d, tumor_numbers_e],
                      y_task='append',
                      plot_type='number_of_tumors',
                      x_axis_label='TUMOR SIZE',
                      y_axis_label='NUMBER OF TUMORS',
                      plot_title='NUMBER OF TUMORS ACCORDING TO SIZE',
                      curve_label_list=['(a)',
                                        '(b): ' + r'$a_{2}$',
                                        '(c): ' + r'$\alpha_{2}$',
                                        '(d): ' + r'$\alpha_{3}$',
                                        '(e): ' + r'$b_{2}$'],
                      export_plot_name='Largest_Metastasis_'+str(t_max))


# FINDING LARGEST TUMORS /CONTINUOUS RADIOTHERAPY/
# CURVE (a)
tumor_size_a_continuous, \
    tumor_numbers_a_continuous = tools.findLargestMetastasis(x_data=x_a_continuous, rho_data=rho_a_continuous,
                                                             export_filename='largest_a_continuous')
# CURVE (b)
tumor_size_b_continuous, \
    tumor_numbers_b_continuous = tools.findLargestMetastasis(x_data=x_b_continuous, rho_data=rho_b_continuous,
                                                             export_filename='largest_b_continuous')
# CURVE (c)
tumor_size_c_continuous,\
    tumor_numbers_c_continuous = tools.findLargestMetastasis(x_data=x_c_continuous, rho_data=rho_c_continuous,
                                                             export_filename='largest_c_continuous')
# CURVE (d)
tumor_size_d_continuous, \
    tumor_numbers_d_continuous = tools.findLargestMetastasis(x_data=x_d_continuous, rho_data=rho_d_continuous,
                                                             export_filename='largest_d_continuous')
# CURVE (e)
tumor_size_e_continuous, \
    tumor_numbers_e_continuous = tools.findLargestMetastasis(x_data=x_e_continuous, rho_data=rho_e_continuous,
                                                             export_filename='largest_e_continuous')

# PLOTTING
plotter.plotAllCurves(x_list=[tumor_size_a_continuous, tumor_size_b_continuous, tumor_size_c_continuous,
                              tumor_size_d_continuous, tumor_size_e_continuous],
                      y_list=[tumor_numbers_a_continuous, tumor_numbers_b_continuous, tumor_numbers_c_continuous,
                              tumor_numbers_d_continuous, tumor_numbers_e_continuous],
                      y_task='append',
                      plot_type='number_of_tumors',
                      x_axis_label='TUMOR SIZE',
                      y_axis_label='NUMBER OF TUMORS',
                      plot_title='NUMBER OF TUMORS ACCORDING TO SIZE AFTER CONTINUOUS RADIOTHERAPY',
                      curve_label_list=['(a)',
                                        '(b): ' + r'$a_{2}$',
                                        '(c): ' + r'$\alpha_{2}$',
                                        '(d): ' + r'$\alpha_{3}$',
                                        '(e): ' + r'$b_{2}$'],
                      export_plot_name='Largest_Metastasis_Continuous_'+str(t_max))


# FINDING LARGEST TUMORS /RADIOTHERAPY IN A SPECIFIC TIME INTERVAL/
# CURVE (a)
tumor_size_a_interval, \
    tumor_numbers_a_interval = tools.findLargestMetastasis(x_data=x_a_interval, rho_data=rho_a_interval,
                                                           export_filename='largest_a_interval')
# CURVE (b)
tumor_size_b_interval, \
    tumor_numbers_b_interval = tools.findLargestMetastasis(x_data=x_b_interval, rho_data=rho_b_interval,
                                                           export_filename='largest_b_interval')
# CURVE (c)
tumor_size_c_interval,\
    tumor_numbers_c_interval = tools.findLargestMetastasis(x_data=x_c_interval, rho_data=rho_c_interval,
                                                           export_filename='largest_c_interval')
# CURVE (d)
tumor_size_d_interval, \
    tumor_numbers_d_interval = tools.findLargestMetastasis(x_data=x_d_interval, rho_data=rho_d_interval,
                                                           export_filename='largest_d_interval')
# CURVE (e)
tumor_size_e_interval, \
    tumor_numbers_e_interval = tools.findLargestMetastasis(x_data=x_e_interval, rho_data=rho_e_interval,
                                                           export_filename='largest_e_interval')

# PLOTTING
plotter.plotAllCurves(x_list=[tumor_size_a_interval, tumor_size_b_interval, tumor_size_c_interval,
                              tumor_size_d_interval, tumor_size_e_interval],
                      y_list=[tumor_numbers_a_interval, tumor_numbers_b_interval, tumor_numbers_c_interval,
                              tumor_numbers_d_interval, tumor_numbers_e_interval],
                      y_task='append',
                      plot_type='number_of_tumors',
                      x_axis_label='TUMOR SIZE',
                      y_axis_label='NUMBER OF TUMORS',
                      plot_title='NUMBER OF TUMORS ACCORDING TO SIZE AFTER INTERVAL RADIOTHERAPY',
                      curve_label_list=['(a)',
                                        '(b): ' + r'$a_{2}$',
                                        '(c): ' + r'$\alpha_{2}$',
                                        '(d): ' + r'$\alpha_{3}$',
                                        '(e): ' + r'$b_{2}$'],
                      export_plot_name='Largest_Metastasis_Interval_'+str(t_max))


# ========================================== COMPUTING TUMORS' MAXIMUM SIZE ===========================================
# MAXIMUM TUMOR SIZE /NO RADIOTHERAPY/
no_treatment_time, no_treatment_tumor = tumorSize.primaryTumor(initial_time=t_0, total_time=t_max, a_G=a_parameter_G,
                                                               b_G=b_parameter_G, x_initial_condition=x_0,
                                                               export_filename=None)

# MAXIMUM TUMOR SIZE /CONTINUOUS RADIOTHERAPY/
untreated_cont_time, untreated_cont_tumor, therapy_cont_time, therapy_cont_tumor = \
    tumorSize.primaryTumorWithRadiotherapy(therapy_type='continuous', time_at_therapy_start=t_therapy_start,
                                           therapy_days=t_therapy, rest_days=t_rest, therapy_sessions=sessions,
                                           initial_time=t_0, total_time=t_max, a_G=a_parameter_G, b_G=b_parameter_G,
                                           D_LQ=D1_parameter_LQ, alpha_LQ=alpha_parameter_LQ, beta_LQ=beta_parameter_LQ,
                                           x_initial_condition=x_0, export_filename=None)

# MAXIMUM TUMOR SIZE /RADIOTHERAPY IN A SPECIFIC TIME INTERVAL/
untreated_interval_time, untreated_interval_tumor, therapy_interval_time, therapy_interval_tumor = \
    tumorSize.primaryTumorWithRadiotherapy(therapy_type='interval', time_at_therapy_start=t_therapy_start,
                                           therapy_days=t_therapy, rest_days=t_rest, therapy_sessions=sessions,
                                           initial_time=t_0, total_time=t_max, a_G=a_parameter_G, b_G=b_parameter_G,
                                           D_LQ=D2_parameter_LQ, alpha_LQ=alpha_parameter_LQ, beta_LQ=beta_parameter_LQ,
                                           x_initial_condition=x_0, export_filename=None)

# PLOTTING
plotter.plotAllCurves(x_list=[no_treatment_time, untreated_interval_time, therapy_interval_time,
                              untreated_cont_time, therapy_cont_time],
                      y_list=[no_treatment_tumor, untreated_interval_tumor, therapy_interval_tumor,
                              untreated_cont_tumor, therapy_cont_tumor],
                      y_task='append', plot_type='max_tumor_size', x_axis_label='TIME /DAYS/',
                      y_axis_label='TUMOR SIZE /CELLS/', plot_title='MAXIMUM TUMOR SIZES FOR DIFFERENT THERAPIES',
                      curve_label_list=['NO THERAPY', 'UNTREATED INTERVAL', 'RADIOTHERAPY INTERVAL',
                                        'PRIOR TO TREATMENT', 'CONTINUOUS RADIOTHERAPY'],
                      export_plot_name='Max_Primary_Tumor_'+str(t_max))


# ================================= COMPUTING MAXIMUM TUMORS FOR DIFFERENT PARAMETERS =================================
# COMPUTING MAXIMUM TUMOR SIZES FOR DIFFERENT VALUES OF PARAMETERS 'a' AND 'b'.
max_tumor_a = iwata.x_max_G(x_init=x_0, t=np.array(range(t_0, t_max+1, k_days)), a_G=a_parameter_G, b_G=b_parameter_G)
max_tumor_b = iwata.x_max_G(x_init=x_0, t=np.array(range(t_0, t_max+1, k_days)), a_G=a2_parameter_G, b_G=b_parameter_G)
max_tumor_e = iwata.x_max_G(x_init=x_0, t=np.array(range(t_0, t_max+1, k_days)), a_G=a_parameter_G, b_G=b2_parameter_G)
max_tumor_f = iwata.x_max_G(x_init=x_0, t=np.array(range(t_0, t_max+1, k_days)), a_G=a2_parameter_G, b_G=b2_parameter_G)
max_tumor = [max_tumor_a, max_tumor_b, max_tumor_e, max_tumor_f]
t = [list(range(t_0, t_max+1, k_days)), list(range(t_0, t_max+1, k_days)),
     list(range(t_0, t_max+1, k_days)), list(range(t_0, t_max+1, k_days))]

# PLOTTING
plotter.plotAllCurves(x_list=t, y_list=max_tumor, y_task='append', plot_type='max_tumor_size',
                      x_axis_label='TIME (DAYS)', y_axis_label='MAX TUMOR SIZE (CELLS)',
                      plot_title='MAXIMUM TUMOR SIZE OVER TIME', curve_label_list=[r'$a_{1},$' + ' ' + r'$b_{1}$',
                                                                                   r'$a_{2},$' + ' ' + r'$b_{1}$',
                                                                                   r'$a_{1},$' + ' ' + r'$b_{2}$',
                                                                                   r'$a_{2},$' + ' ' + r'$b_{2}$'],
                      export_plot_name='Max_Tumor_Size_'+str(t_max))


# LOOKING FOR THE TIME /IN DAYS/ WHEN CONCAVITY CHANGES.
numerical_derivative_DF = tools.findTimeWhenConcavityCahnges(initial_time_c=t_0, initial_maximum_time=365,
                                                             a_G_c=a_parameter_G, b_G_c=b_parameter_G,
                                                             k_days_c=k_days, m_G_c=m_parameter_G,
                                                             alpha_G_c=alpha2_parameter_G, x_initial_condition_c=x_0,
                                                             threshold=5, step_increase=120,
                                                             max_cycles=10000, export_filename='derivative_DF')

# PLOTTING
plotter.plotConcavityChange(derivative_dataframe=numerical_derivative_DF, export_plot_name='Derivative_Plot')

