import subprocess
import pandas as pd
from FindOrthoProt.generate_file_name import generate_id_file
from sequence_analysis.constants import MEME_OUTPUT, FIMO_OUTPUT
import json
import os
import glob

class SearchMotifs:
    def __init__(self, homologous_file, predicted_file, num_motif):
        self.homologous_file = homologous_file
        self.predicted_file = predicted_file
        self.num_motif = num_motif
        self.output_meme = (generate_id_file('meme', '', MEME_OUTPUT)).rstrip('.')
        self.output_fimo = generate_id_file('fimo', '', FIMO_OUTPUT).rstrip('.')

    def find_motifs_with_meme(self):
        print('self.output_meme', self.output_meme)
        print('self.output_fimo',  self.output_fimo)
        meme_command = f" meme '{self.homologous_file}' -o '{self.output_meme}' -nmotifs {self.num_motif}"
        subprocess.run(meme_command, shell=True)

        motifs_of_homologous = {'id': [], 'sequence': [], 'e-value': [], 'width': []}

        with open(f"{self.output_meme}/meme.txt", 'r') as file:
            lines = file.readlines()
            motif_lines = [line.strip() for line in lines if line.startswith('MOTIF')]

        for line in motif_lines:
            result_list = line.split()
            motifs_of_homologous['id'].append(result_list[2])
            motifs_of_homologous['sequence'].append(result_list[1])
            motifs_of_homologous['e-value'].append(result_list[14])
            motifs_of_homologous['width'].append(result_list[5])

        df_motifs_of_homologous = pd.DataFrame(motifs_of_homologous)
        print(df_motifs_of_homologous)
        meme_html_file = self.output_meme + '/' +'meme.html'
        print('meme finished')
        print('meme_html_file', meme_html_file)

        # Преобразуйте DataFrame в список словарей
        motifs_list = df_motifs_of_homologous.to_dict(orient='records')
        # Сериализуйте список словарей в формат JSON
        motifs_json = json.dumps(motifs_list)

        all_files = glob.glob(os.path.join(self.output_meme, '*'))
        # Удаляем все файлы, кроме определенного
        for file in all_files:
            if not any(file.endswith(ext) for ext in ['.txt', '.html']):
                os.remove(file)

        return 'Meme', motifs_json, meme_html_file

    def fimo_run(self):
        print('Starting fimo run')
        meme_file = self.output_meme + '/' + 'meme.txt'
        print('meme_file', meme_file)
        fimo_command = f"fimo --oc '{self.output_fimo}' --verbosity 1 '{meme_file}' '{self.predicted_file}'"
        subprocess.run(fimo_command, shell=True)

        fimo_result = self.output_fimo + '/' + 'fimo.tsv'

        # Загрузите данные из файла fimo.tsv в DataFrame
        df_fimo = pd.read_csv(fimo_result, sep='\t')

        # Очистка и обработка данных
        df_fimo = df_fimo.iloc[:-3]
        df_fimo = df_fimo.drop(df_fimo.columns[[0, 3, 4, 5]], axis=1)
        df_fimo = df_fimo[['sequence_name'] + [col for col in df_fimo.columns if col != 'sequence_name']]
        df_fimo = df_fimo.sort_values(by=['sequence_name', 'motif_alt_id'])
        df_fimo = df_fimo.reset_index(drop=True)
        df_fimo['p-value'] = df_fimo['p-value'].apply(lambda x: '{:.2e}'.format(x))
        df_fimo['q-value'] = df_fimo['q-value'].apply(lambda x: '{:.2e}'.format(x))
        print(df_fimo)

        fimo_html_file = self.output_fimo + '/' + 'fimo.html'

        # Преобразуйте DataFrame в список словарей
        df_fimo_list = df_fimo.to_dict(orient='records')

        # Сериализуйте список словарей в формат JSON
        df_fimo_json = json.dumps(df_fimo_list)

        all_files = glob.glob(os.path.join(self.output_fimo, '*'))
        # Удаляем все файлы, кроме определенного
        for file in all_files:
            if not any(file.endswith(ext) for ext in ['.tsv', '.html']):
                os.remove(file)

        return 'Fimo', df_fimo_json, fimo_html_file