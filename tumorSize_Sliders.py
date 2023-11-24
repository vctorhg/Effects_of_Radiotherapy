
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider


total_days = 50    # WE TAKE 50 BASED ON FIG.3 FROM THE LIN_QUAD PAPER IN THESIS>SOURCES.
hourly_mesh = (total_days * 24) + 1    # SETTING HOURLY MESH.
# DEFINING THE TIME VECTOR.
total_t = np.linspace(0, total_days, hourly_mesh)
non_therapy_days = 2
therapy_days = 5

# DEFINING ALL NECESSARY PARAMETERS.
x_0 = 10**9
a = 0.00286
b = 7.3 *(10**10)
alpha = .01
#beta = .02
beta = alpha/10
D = 1.8

rho_0 = .001


# DEFINING FUNCTIONS TO BE PLOTTED:
def therapy(x_0, alpha, beta, D, t):
    # COMPUTES TUMOR SIZE UNDER THERAPY FOR THE GIVEN TIME INTERVAL.
    #
    # INPUT:
    #   x_0 : INT. INITIAL TUMOR SIZE.
    #   alpha : FLOAT.
    #   beta : FLOAT.
    #   D : FLOAT. DOSAGE (?)
    #   t : ARRAY. CONTAINS TIME VALUES IN THE RANGE [0, therapy_days] WITH HOURLY MESH.
    #
    # OUTPUT:
    #   ARRAY. CONTAINS TUMOR SIZE VALUES THAT FOLLOW THE LIN-QUAD GROWTH RATE FROM ENDERLING.

    return x_0 * np.exp(-((alpha*D*t) + (beta*(D**2)*(t**2))))
def noTherapy(x_0, a, b, t):
    # COMPUTES TUMOR SIZE WHEN NO THERAPY IS APPLIED FOR THE GIVEN TIME INTERVAL.
    #
    # INPUT:
    #   x_0 : INT. INITIAL TUMOR SIZE.
    #   a : FLOAT. GROWTH RATE CONSTANT /FROM IWATA/.
    #   b : FLOAT. TUMOR SIZE AT SATURATED LEVEL /FROM IWATA/.
    #   t : ARRAY. CONTAINS TIME VALUES IN THE RANGE [0, non_therapy_days] WITH HOURLY MESH.
    #
    # OUTPUT:
    #   ARRAY. CONTAINS TUMOR SIZE VALUES THAT FOLLOW THE GOMPERTZ GROWTH RATE FROM IWATA.
    return (x_0**(np.exp(-a*t)))*(b**(1-np.exp(-a*t)))
def timePeriod (days):
    # DROPS LAST ENTRY OF TIME ARRAY IN ORDER OT HAVE TIME INTERVALS OF THE FORM: [t_0, t_n)
    #
    # INPUT:
    #   days : INT. NUMBER OF DAYS THE TIME ARRAY WILL SPAN.
    #
    # OUTPUT:
    #   ARRAY. CONTAINS TIME VALUES IN THE RANGE [0, days) WITH HOURLY MESH.
    return np.linspace(0, days, (days*24)+1)[:days*24]
def smartAppend(tumor_size, tumor_vals, total_t_length):
    # MAKES SURE THAT 'tumor_size' LENGTH DOES NOT EXCEED THAT OF 'total_t' AFTER
    # APPENDING 'tumor_vals' TO THE 'tumor_size' ARRAY.
    #
    # INPUT:
    #   tumor_size : ARRAY. CONTAINS ACCUMULATED TUMOR SIZE VALUES.
    #   tumor_vals : ARRAY. CONTAINS RECENTLY COMPUTED TUMOR SIZE VALUES.
    #   total_t_length : INT. LENGTH OF THE TOTAL TIME ARRAY.
    #
    # OUTPUT:
    #   tumor_size : ARRAY. CONSISTS OF 'tumor_size' CONCATENATED WITH AS MANY VALUES OF
    #                'tumor_vals' AS POSSIBLE SUCH THAT THE NEW LENGTH DOES NOT EXCEED
    #                'total_t_length'.
    #   from_idx :  INT. STARTING INDEX OF THE APPENDED VALUES.
    #   to_idx: INT. ENDING INDEX OF THE APPENDED VALUES.

    from_idx = len(tumor_size)
    tumor_size = np.append(tumor_size, tumor_vals)
    if len(tumor_size) > total_t_length:
        to_idx = total_t_length
        tumor_size = tumor_size[0:to_idx]
    else:
        to_idx = len(tumor_size)
    return tumor_size, from_idx, to_idx
def tumorSize(x_0, a, b, alpha, beta, D, total_days, therapy_days, non_therapy_days):
    # COMPUTES THE TUMOR SIZE OVER 'total days' CONSIDERING THE NECESSARY SWITCH
    # IN GROWTH RATES ACCORDING TO 'therapy_days', AND 'non_therapy_days'.
    #
    # INPUT:
    #   x_0 : INT. INITIAL CONDITION
    #   a :
    #   b :
    #   alpha :
    #   beta :
    #   D :
    #   total_days :
    #   therapy_days :
    #   non_therapy_days :
    #
    # OUTPUT:
    #   tumor_size : ARRAY. CONTAINS TUMOR SIZE VALUES AT EACH TIME POINT IN THE
    #                RANGE [0, total, days] WITH AN HOURLY GRID.

    therapy_indicator = False
    tumor_size = np.array([])
    total_t = np.linspace(0, total_days, (total_days*24)+1)
    total_t_length = len(total_t)


    while len(tumor_size) < (total_days*24) +1:
        if therapy_indicator == False:
            t = timePeriod(non_therapy_days)
            try:
                tumor_vals = noTherapy(tumor_size[-1], a, b, t)
            except:
                tumor_vals = noTherapy(x_0, a, b, t)
            tumor_size, from_idx, to_idx = smartAppend(tumor_size, tumor_vals, total_t_length)
            therapy_indicator = True

        else:
            t = timePeriod(therapy_days)
            tumor_vals = therapy(tumor_size[-1], alpha, beta, D, t)
            tumor_size, from_idx, to_idx = smartAppend(tumor_size, tumor_vals, total_t_length)
            therapy_indicator = False

    return tumor_size

tumor_size = tumorSize(x_0, a, b, alpha, beta, D, total_days, therapy_days, non_therapy_days)

# SETTING MAIN FIGURE CHARACTERISTICS
fig, ax = plt.subplots()
tumor_plot, = ax.plot(total_t, tumorSize(x_0, a, b, alpha, beta, D, total_days, therapy_days, non_therapy_days))
ax.set_xlabel('Time [d]')
ax.set_ylabel('Tumor Size')
ax.set_title('Tumor Size Evolution over Time')

# SETTING MARGIN SO PLOT DOESN'T OVERLAP WITH SLIDER
fig.subplots_adjust(bottom=.4)

# DECLARING SLIDERS
aSlider = fig.add_axes([0.15, 0.1, 0.70, 0.01])
a_slider = Slider(ax=aSlider, label='Growth Rate',
                valmin=0, valmax=0.5, valinit=a)

bSlider = fig.add_axes([0.15, 0.15, 0.70, 0.01])
b_slider = Slider(ax=bSlider, label='Satuareted\nTumor Size',
                valmin=3.65*(10**10), valmax=1*(10**12), valinit=b)

alphaSlider = fig.add_axes([0.15, 0.2, 0.70, 0.01])
alpha_slider = Slider(ax=alphaSlider, label='alpha',
                valmin=0.025, valmax=0.036, valinit=alpha)

betaSlider = fig.add_axes([0.15, 0.25, 0.70, 0.01])
beta_slider = Slider(ax=betaSlider, label='beta',
                valmin=0.0025, valmax=0.0036, valinit=beta)

DSlider = fig.add_axes([0.15, 0.3, 0.70, 0.01])
D_slider = Slider(ax=DSlider, label='D',
                valmin=0, valmax=1, valinit=D)

# DECLARING UPDATE FUNCTION
def update(val):
    tumor_plot.set_ydata(tumorSize(x_0, a_slider.val, b_slider.val,
                                    alpha_slider.val, beta_slider.val,
                                    D_slider.val, total_days, therapy_days, non_therapy_days))
    fig.canvas.draw_idle()

def updateBeta(val):
    beta_slider.set_val(alpha_slider.val/10)
    tumor_plot.set_ydata(tumorSize(x_0, a_slider.val, b_slider.val,
                                    alpha_slider.val, beta_slider.val,
                                    D_slider.val, total_days, therapy_days, non_therapy_days))
    fig.canvas.draw_idle()

def updateAlpha(val):
    alpha_slider.set_val(beta_slider.val*10)
    tumor_plot.set_ydata(tumorSize(x_0, a_slider.val, b_slider.val,
                                    alpha_slider.val, beta_slider.val,
                                    D_slider.val, total_days, therapy_days, non_therapy_days))
    fig.canvas.draw_idle()
#def linkSlider(val):
#    beta.set_val(alpha_slider.val/10)
#    update(val)

# CALLING UPDATE WHENEVER A SLIDER VALUE IS CHANGED
a_slider.on_changed(update)
b_slider.on_changed(update)
alpha_slider.on_changed(updateBeta)
beta_slider.on_changed(updateAlpha)
D_slider.on_changed(update)

plt.show()