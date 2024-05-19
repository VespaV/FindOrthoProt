import subprocess
import pandas as pd
from pymemesuite.fimo import FIMO
import pymemesuite
from pymemesuite.common import MotifFile
from pymemesuite.common import Sequence
import Bio.SeqIO

# Получение списка всех атрибутов в pymemesuite
all_attributes = dir(pymemesuite)

class SearchMotifs:
    def __init__(self, homologous_file, predicted_file):
        self.homologous_file = homologous_file
        self.predicted_file = predicted_file
        self.output_directory = '/Users/maria/projects/FindOrthoProt/sequence_analysis/output_files/motif_search/new_meme'
        self.output_fimo = '/Users/maria/projects/FindOrthoProt/sequence_analysis/output_files/motif_search/new_fimo'

    def find_motifs_with_meme(self):
        meme_command = f"meme '{self.homologous_file}' -o {self.output_directory} -nmotifs 9"
        subprocess.run(meme_command, shell=True)

        motifs_of_homologous = {'id': [], 'sequence': [], 'e-value': [], 'width': []}

        with open(f"{self.output_directory}/meme.txt", 'r') as file:
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

        return df_motifs_of_homologous

    def fimo_run(self):
        motifs_file_st = f"{self.output_directory}/meme.txt"
        with MotifFile(motifs_file_st) as motif_file:
            motif = motif_file.read()

        print(motif.name.decode())
        print(motif.consensus)

        for row in motif.frequencies:
            print(" ".join(f'{freq:.2f}' for freq in row))
            
        
        sequences = [
            Sequence(str(record.seq), name=record.id.encode())
            for record in Bio.SeqIO.parse(self.predicted_file, "fasta")
        ]

        fimo = FIMO(both_strands=False)
        pattern = fimo.score_motif(motif, sequences, motif_file.background)

        for m in pattern.matched_elements:
            print(
                m.source.accession.decode(),
                m.start,
                m.stop,
                m.strand,
                m.score,
                m.pvalue,
                m.qvalue
            )

