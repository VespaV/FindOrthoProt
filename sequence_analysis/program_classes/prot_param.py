
from Bio import SeqIO
from FindOrthoProt.generate_file_name import generate_id_file
from collections import Counter
import pandas as pd
import numpy as np
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import plotly.express as px
import plotly.graph_objects as go
import json
import plotly.io as pio
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from sequence_analysis.constants import *

class ProteinPropertiesAnalyzer:
    def __init__(self, homologous_fasta, predicted_fasta):
        self.homologous_fasta = homologous_fasta
        self.predicted_fasta = predicted_fasta
        self.amino_acids_to_keep = ["M", "P", "L", "D", "C", "S", "N", "V", "G", "F", "E", "I", "Y", "W", "A", "R", "T", "K", "Q", "H"]
        self.predicted_params = self.prot_param(self.predicted_fasta)
        self.homologous_params = self.prot_param(self.homologous_fasta)

    def count_all_seq_amino_acids(self):
        df_homologous = self.count_amino_acids(self.homologous_fasta)
        df_predicted = self.count_amino_acids(self.predicted_fasta)
        joined_df = pd.concat([df_predicted, df_homologous])
        joined_df = joined_df[[col for col in joined_df.columns if col != 'Total'] + ['Total']]
        print(joined_df)

        return joined_df

    def count_amino_acids(self, input_file):
        sequences = list(SeqIO.parse(input_file, "fasta"))
        aa_composition = {}
        sequence_names = []
        total_sum = []

        for seq_record in sequences:
            sequence_name = seq_record.id
            sequence_names.append(sequence_name)
            sequence = str(seq_record.seq)
            total_sum.append(len(sequence))
            aa_counts = Counter(sequence)
            aa_composition[sequence_name] = aa_counts

        df_col = pd.DataFrame(aa_composition).T.fillna(0)
        df_col["Total"] = total_sum

        print(df_col)

        return df_col

    def amino_gistogram(self):
        df_joined = self.count_all_seq_amino_acids()
        columns_to_keep = self.amino_acids_to_keep + ["Total"]
        df = df_joined[columns_to_keep]
        df = (df.div(df["Total"].values, axis=0) * 100).round(2)

        colors = px.colors.qualitative.Alphabet[:len(df.columns) - 1]
        print(len(df.columns))

        sequence_names = df.index

        fig = go.Figure()

        for i, amino_acid in enumerate(df.columns[:-1]):
            amino_acid_counts = df[amino_acid]
            cumulative_counts = amino_acid_counts.groupby(level=0).cumsum()

            if i == 0:
                fig.add_trace(
                    go.Bar(x=sequence_names, y=cumulative_counts, name=amino_acid, marker=dict(color=colors[i])))
            else:
                fig.add_trace(go.Bar(x=sequence_names, y=cumulative_counts, name=amino_acid,
                                     marker=dict(color=colors[i], line=dict(color=colors[i]))))

        fig.update_layout(
            plot_bgcolor='white',
            autosize=False,
            width=1200,
            height=800,
            margin=dict(l=200, r=200, b=150, t=100, pad=5),
            xaxis=dict(
                title=dict(text='Sequence Names', standoff=30, font=dict(color='black')),
                tickmode='array',
                tickvals=list(range(len(sequence_names))) + [len(sequence_names)],
                tickfont=dict(size=16, color='black'),
            ),
            yaxis=dict(tickfont=dict(size=16, color='black')),

            xaxis_title=dict(text='Sequence Names', font=dict(size=20, color='black')),
            yaxis_title=dict(text='Percentage', font=dict(size=20, color='black')),
            title=dict(text='Accumulated Histogram for Amino Acid Counts', font=dict(size=25, color='black'), x=0.5, y=0.9),
            barmode='stack',
            showlegend=True,
            legend=dict(
                x=1.05,
                y=0.5,
                font=dict(size=12),
                traceorder='reversed',
                orientation='v',
            ),
        )

        picture_path = generate_id_file('amino_count', 'png', AMINO_COUNT_PICTURES)
        pio.write_image(fig, picture_path)
        print(picture_path)
        df_joined = df_joined[~df_joined.index.duplicated(keep='first')]

        json_data = df_joined.to_json(orient='index')

        return "AminoGist", picture_path, json_data

    def prot_param(self, file):
        params = {}
        molecular_weight = {}
        aromaticity = {}
        instability_index = {}
        gravy = {}
        iso_point = {}
        secondary_structure_fraction = {}

        for record in SeqIO.parse(file, 'fasta'):
            filtered_seq = ''.join(aa for aa in record.seq if aa in self.amino_acids_to_keep)
            sequence = ProteinAnalysis(filtered_seq)
            molecular_weight[record.id] = "%0.2f" % sequence.molecular_weight()
            aromaticity[record.id] = "%0.2f" % sequence.aromaticity()
            instability_index[record.id] = "%0.2f" % sequence.instability_index()
            gravy[record.id] = "%0.2f" % sequence.gravy()
            iso_point[record.id] = sequence.isoelectric_point()
            secondary_structure_fraction[record.id] = sequence.secondary_structure_fraction()

        params['molecular_weight'] = molecular_weight
        params['aromaticity'] = aromaticity
        params['instability_index'] = instability_index
        params['gravy'] = gravy
        params['iso_point'] = iso_point
        params['secondary_structure_fraction'] = secondary_structure_fraction

        return params

    def combine_params(self):
        # Combine parameters from predicted and homologous sequences
        all_params = {}
        for param_name in self.predicted_params.keys():
            all_params[param_name] = {}
            for seq_id, value in self.predicted_params[param_name].items():
                # Extracting the first element of the tuple and converting to float
                all_params[param_name][seq_id] = round(float(value[0] if isinstance(value, tuple) else value), 2)

            # Handle homologous values
            homologous_values = [float(val[0]) if isinstance(val, tuple) else float(val) for val in
                                 self.homologous_params[param_name].values()]
            homologous_mean = round(np.mean(homologous_values), 2)
            homologous_std = round(np.std(homologous_values), 2)
            all_params[param_name]['Mean_homologous'] = homologous_mean
            all_params[param_name]['Std_homologous'] = homologous_std

        # Create a DataFrame
        df = pd.DataFrame(all_params)
        df = df.iloc[:, :-1]
        print('Общий df', df)

        json_data = df.to_json(orient='index')

        return json_data

    def bar_chart(self, where_safe, param, title, y_label):
        print('start bar_chart')

        predict_dict_values = self.predicted_params[param]
        homologues_dict_values = self.homologous_params[param]
        values_dict1 = [float(value) for value in predict_dict_values.values()]
        values_dict2 = [float(value) for value in homologues_dict_values.values()]

        # Построение столбчатой диаграммы
        fig, ax = plt.subplots(figsize=((len(values_dict1) + len(values_dict2)) * 0.4, 6))
        # plt.tight_layout()
        # Столбцы для значений из predict_dict_values
        ax.bar(np.arange(len(values_dict1)), values_dict1, label='Predicted sequences')

        # Единый столбец для среднего значения из homologues_dict_values
        mean_dict2 = np.mean(values_dict2)
        ax.bar(len(values_dict1), mean_dict2, label='Mean homologues', color='#f56042')

        # Уса стандартного отклонения для homologues_dict_values
        std_dict2 = np.std(values_dict2)
        ax.errorbar(len(values_dict1), mean_dict2, yerr=std_dict2, fmt='o', color='black', capsize=5, linewidth=1,
                    markersize=3)

        ax.axhline(mean_dict2 + std_dict2, linestyle='dashed', color='black', linewidth=1, alpha=0.3,
                   label='Upper bound')
        ax.axhline(mean_dict2 - std_dict2, linestyle='dashed', color='black', linewidth=1, alpha=0.3,
                   label='Lower bound')

        # Настройка меток оси x
        ax.set_xticks(list(range(len(values_dict1) + 1)))
        ax.set_xticklabels(list(predict_dict_values.keys()) + ['Mean homologues'], fontdict={'family': 'Arial'})
        ax.tick_params(axis='y', labelsize=6)

        ax.spines['bottom'].set_color('lightgray')
        ax.spines['top'].set_color('lightgray')
        ax.spines['right'].set_color('lightgray')
        ax.spines['left'].set_color('lightgray')

        ax.set_axisbelow(True)
        ax.yaxis.grid(color='lightgray', linestyle='-', linewidth=0.4, zorder=-2)

        # Убрать горизонтальные линии ticks на оси y
        ax.yaxis.set_ticks_position('none')

        # Убрать вертикальные линии ticks на оси x
        ax.xaxis.set_ticks_position('none')

        # Добавление легенды
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), frameon=False, prop={'size': 8, 'weight': 'normal'})

        # Настройка заголовка и меток осей
        ax.set_title(title, fontdict={'size': 15, 'weight': 'normal'})
        ax.set_xlabel('Sequence', fontdict={'fontsize': 10, 'weight': 'normal'}, labelpad=10)
        ax.set_ylabel(y_label, fontdict={'fontsize': 10, 'weight': 'normal'}, labelpad=10)

        # Настройка меток оси x с поворотом
        ax.set_xticks(list(range(len(values_dict1) + 1)))
        ax.set_xticklabels(list(predict_dict_values.keys()) + ['Mean homologues'],
                           rotation=45, ha='right', fontsize=6, weight='normal')  # Поворот меток оси x на 45 градусов вправо

        fig.subplots_adjust(left=0.18, right=0.70, bottom=0.18, top=0.85)

        picture_path = generate_id_file('amino_count', 'png', where_safe)
        canvas = FigureCanvas(fig)
        canvas.print_png(picture_path)

        return picture_path

    def all_params(self):
        results = {}

        molecular_weight = self.molecular_weight()
        aromaticity = self.aromaticity()
        instability_index = self.instability_index()
        gravy = self.gravy()
        iso_point = self.iso_point()

        results["molecular_weight"] = molecular_weight
        results["aromaticity"] = aromaticity
        results["instability_index"] = instability_index
        results["gravy"] = gravy
        results["iso_point"] = iso_point

        combine_params = self.combine_params()

        return "AllParams", json.dumps(results), combine_params

    def molecular_weight(self):
        return self.bar_chart(
            MOLECULAR_WEIGHT_PICTURES,
            "molecular_weight", "Molecular weight", "MW values, g/mol")

    def aromaticity(self):
        return self.bar_chart(AROMACITY_PICTURES, "aromaticity", "Lobry aromaticity values",
                              "Relative frequency Phe+Trp+Tyr")

    def instability_index(self):
        return self.bar_chart(INSTABILITY_INDEX_PICTURES, "instability_index", "Guruprasad instability index",
                              "Index value")

    def gravy(self):
        return self.bar_chart(GRAVY_PICTURES, "gravy", "GRAVY according to Kyte and Doolittle", "Index value")

    def iso_point(self):
        return self.bar_chart(ISO_POINT_PICTURES, "iso_point", "Isoelectric point (pI)", "pH value")

    def build_df_secondary_structure_fraction(self):
        predict_dict_values = self.prot_param(self.predicted_fasta)['secondary_structure_fraction']
        homologues_dict_values = self.prot_param(self.homologous_fasta)['secondary_structure_fraction']

        df_pred = pd.DataFrame.from_dict(predict_dict_values, orient='index', columns=['helix', 'turn', 'sheet'])
        df_hom = pd.DataFrame.from_dict(homologues_dict_values, orient='index', columns=['helix', 'turn', 'sheet'])
        print(df_pred)
        print(df_hom)

        joined_df = pd.concat([df_pred, df_hom])
        print('joined_df', joined_df)

        return df_pred, df_hom, joined_df

    def secondary_structure_gistogram(self):
        print('secondary_structure_gistogram start')
        df_combined = self.build_df_secondary_structure_fraction()[2]

        fig = px.bar(df_combined, barmode='stack', labels={'index': 'Sequence', 'value': 'Fractions'},
                     title='Fractions of secondary structures',
                     category_orders={'index': list(df_combined.index)})

        fig.update_layout(
            margin=dict(l=200, r=200, b=150, t=100, pad=5),
            legend=dict(x=1.02, y=0.5, font=dict(size=20, color='black'),),
            legend_title=' ',
            font=dict(family='Arial'),
            width=1200,
            height=800,
            plot_bgcolor='white',
            yaxis=dict(showgrid=True, gridcolor='lightgray', gridwidth=2, tickfont=dict(size=14, color='black')),
            showlegend=True,
            xaxis=dict(
                title=dict(text='Sequence Names', standoff=30, font=dict(size=20, color='black')),
                tickmode='array',
                tickfont=dict(size=16, color='black')
            ),
            yaxis_title=dict(text='Cumulative Count', font=dict(size=20, color='black')),
            title=dict(text='Accumulated Histogram for Amino Acid Counts', font=dict(size=25, color='black'), x=0.5, y=0.9),
        )

        picture_path = generate_id_file('amino_count', 'png', SECONDARY_STUCTURE_PICTURES)
        fig.write_image(picture_path)

        df_combined = df_combined[~df_combined.index.duplicated(keep='first')]

        json_data = df_combined.to_json(orient='index')

        return "SecStruct", picture_path, json_data

