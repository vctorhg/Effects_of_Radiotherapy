

import pandas as pd
import seaborn.objects as so


# FILE CREATED ON 25.09.2023 WITH THE INTENTION OF PRODUCING PLOTS FOR THE
# BIRTH DATA TAKEN FROM THE GERMAN FEDERAL STATISTICS OFFICE DATABASE.

birth_df = pd.read_csv('/Users/victor/Documents/TUM/Thesis/Data/Births_Germany_2022_Formatted.txt', sep=';')

birth_plot = \
    (
    so.Plot(data=birth_df, x='Age', y='Births')
        .add(so.Bar(color='steelblue'))
        .label(title='MOTHER\'S AGE WHEN GIVING BIRTH TO FIRST CHILD',
               x='AGE (YEARS)', y='REGISTERED NUMBER OF BIRTHS')
    )

birth_plot.show()