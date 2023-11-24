

import pandas as pd
import seaborn.objects as so


# FILE CREATED ON 25.09.2023 WITH THE INTENTION OF PRODUCING PLOTS FOR THE
# DEATH DATA TAKEN FROM THE GERMAN FEDERAL STATISTICS OFFICE DATABASE.

# IMPORTING DATA FROM CSV FILE
female_death_df = pd.read_csv('/Users/victor/Documents/TUM/Thesis/Data/Deaths_W_Germany_2022.txt', sep=';')

death_plot = \
    (
    so.Plot(data=female_death_df, x='Age', y='Deaths')
        .add(so.Bar(color='steelblue'))
        .label(title='WOMENS\' AGE AT DEATH',
               x='AGE (YEARS)', y='REGISTERED NUMBER OF DEATHS')
    )

death_plot.show()


