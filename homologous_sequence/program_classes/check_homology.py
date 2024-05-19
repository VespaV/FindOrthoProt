from Bio import Align
from Bio import SeqIO
import pandas as pd
import subprocess
from Bio import AlignIO
from Bio.Align import AlignInfo
from concurrent.futures import ThreadPoolExecutor
import tempfile
import os
import itertools
from homologous_sequence.constants import *
from Bio.Seq import Seq
from FindOrthoProt.generate_file_name import generate_id_file
from scipy.stats import shapiro, ranksums, iqr

class Homology:
    def __init__(self, sequences_type, search_seq_file, sequences_file, num_blast_results = 5, 
                amount=10, algorithm = 'hmmer3', algorithm_blast = 'psiblast'):
        self.sequences_type = sequences_type
        self.search_seq_file = search_seq_file
        self.sequences_file = sequences_file
        self.score_output_file = generate_id_file('score_output_file', 'txt', SCORE_OUTPUT_FOLDER)
        self.amount = amount
        self.output_blast_file = generate_id_file(name='BLAST_out', resolution='xlsx', folder=BLAST_XLSX_FOLDER)
        self.num_blast_results = str(num_blast_results)
        self.algorithm = algorithm
        self.output_alignment_file = generate_id_file('output_alignment', 'aln', ALIGNMENT_FOLDER)

        self.algorithm_blast = algorithm_blast
        self.new_fasta = generate_id_file('output_fasta', 'txt', FASTA_OUTPUT_FOLDER)
        self.result_blast = generate_id_file('output_blast_cvs', 'txt', BLAST_CVS_FOLDER)

    def run_muscle(self):
        input_file = self.search_seq_file
        output_file = self.output_alignment_file
        command = f"muscle -in '{input_file}' -out '{output_file}'"
        subprocess.run(command, shell=True)
        
        print ('output_align', output_file)

        return output_file
        
    def _consensus(self):
        alignment = AlignIO.read(open(self.output_alignment_file), "fasta")
        summary_align = AlignInfo.SummaryInfo(alignment)
        consensus_seq = summary_align.dumb_consensus()
        print('Consensus seq:', consensus_seq)
    
        return consensus_seq
    
    def search_homologous(self):
        answer = None
        try:
            if  self.algorithm == 'hmmer3':
                answer = self.HMMER3()
            elif self.algorithm == "pair_aligment":
                answer = self.main_pair_aligment()
        except Exception as e:
            print(f"Произошла ошибка: {e}")
        
        print ('answer', answer)
        return answer, self.new_fasta
         
    def _pairwise_alignment(self, seq1, reference_sequence):
        aligner = Align.PairwiseAligner()
        aligner.mode = 'global'

        if self.sequences_type == 'Protein':
            aligner.open_gap_score = -3
            aligner.extend_gap_score = -1
            aligner.substitution_matrix = Align.substitution_matrices.load("BLOSUM62")

        if self.sequences_type == 'DNA':
            aligner.open_gap_score = -4
            aligner.extend_gap_score = -2
            aligner.substitution_matrix = Align.substitution_matrices.load("NUC.4.4")

        alignments = aligner.align(seq1.seq, reference_sequence)
        alignment = alignments[0]
        print(seq1.name, 'Ok')

        return alignment.score

    def _process_chunk(self, chunk, reference_sequence):
        results = {}
        for seq_record in chunk:
            alignment_result = self._pairwise_alignment(seq_record, reference_sequence)
            results[seq_record.id] = alignment_result
            print("results", results)
        return results

    def main_pair_aligment(self):
        input_file = self.sequences_file
        reference_sequence = self._consensus()
        output_cvs = self.score_output_file

        # Разделение на части по 50 последовательностей
        chunk_size = 50
        chunks = []
        with open(input_file, "r") as handle:
            records_iter = SeqIO.parse(handle, "fasta")
            while True:
                chunk_records = list(itertools.islice(records_iter, chunk_size))
                if not chunk_records:
                    break  # Выход из цикла, если достигнут конец файла
                chunks.append(chunk_records)

        # Обработка остатка (если он есть)
        remaining_records = list(records_iter)
        if remaining_records:
            chunks.append(remaining_records)

        # Параллельное выполнение операций с использованием ThreadPoolExecutor
        with ThreadPoolExecutor() as executor:
            results = list(executor.map(lambda chunk: self._process_chunk(chunk, reference_sequence), chunks))

        result_list = []


        for result_chunk in results:
            for seq_id, alignment_result in result_chunk.items():
                result_list.append((seq_id, alignment_result))  # Используем score из alignment_result

        print("result_list", result_list)
        df = pd.DataFrame(result_list, columns=["ID", "Score"])
        print('df', df)
        df.sort_values(by="Score", ascending=False, inplace=True)
        df.reset_index(drop=True, inplace=True)

        print("result_list", result_list)

        stat, p_value_shapiro = shapiro(df["Score"])
        print(f"Shapiro-Wilk Test Statistic: {stat}, P-Value: {p_value_shapiro}")

        # Проверка нормальности
        if p_value_shapiro > 0.05:
            print("Данные вероятно распределены нормально.")
            df["Z-value"] = self._calculate_z_value(df["Score"])


        else:
            print("Данные вероятно не распределены нормально.")
            df = df[df["Score"] >= 0]
            median_value = df["Score"].median()
            iqr_value = iqr(df["Score"])
            df['median'] = median_value
            df['IQR'] = iqr_value

            df['upper_threshold'] = median_value + 2 * iqr_value


        top_sequences = df.head(int(self.amount))["ID"].tolist()
        print('top_sequences', top_sequences)

        df.to_csv(output_cvs, index=False, sep='\t')
        print(df)

        # with open(self.new_fasta, "w") as output_handle:
        #     with open(self.sequences_file, "r") as input_handle:
        #         for record in SeqIO.parse(input_handle, "fasta"):
        #             if record.id in top_sequences:
        #                 SeqIO.write(record, output_handle, "fasta")

        # Перезапись файла self.new_fasta с учетом порядка из top_sequences
        with open(self.new_fasta, "w") as output_handle:
            # Перебор идентификаторов последовательностей из списка top_sequences
            for seq_id in top_sequences:
                # Перебор записей из входного файла последовательностей
                with open(self.sequences_file, "r") as input_handle:
                    for record in SeqIO.parse(input_handle, "fasta"):
                        # Проверка, совпадает ли идентификатор записи с текущим идентификатором из top_sequences
                        if record.id == seq_id:
                            # Если да, записываем текущую последовательность в файл
                            SeqIO.write(record, output_handle, "fasta")
                            break  # Выходим из цикла по записям, так как уже записали текущую последовательность

        return output_cvs
        
    def _calculate_z_value(self, data):
            mean = data.mean()
            std_dev = data.std()
            z_values = (data - mean) / std_dev
            return z_values
        
            
    def aln_to_stockholm(self):
        alignment = AlignIO.read(self.output_alignment_file, 'fasta')
        alignment_file_sto = self.output_alignment_file[:-3] + "sto" 
        print(alignment_file_sto)
        AlignIO.write([alignment], alignment_file_sto, 'stockholm')
        
        return alignment_file_sto
    

    def HMMER3(self):
            output_align = self.aln_to_stockholm()
            result_file_name = self.score_output_file 

            with tempfile.NamedTemporaryFile(delete=False, dir=TEMP_FILES_DIR, 
                                             suffix='.hmm') as profile_file:
                profile_file_name = profile_file.name              
                fasta_sequences_file = self.sequences_file

                try:
                    build_command = f'hmmbuild "{profile_file_name}" "{output_align}"'
                    subprocess.run(build_command, shell=True, check=True)
                except subprocess.CalledProcessError as e:
                    print(f"Error running command: {e}")
                    print(f"Command output: {e.output}")
                    raise
            
                search_command = f"hmmsearch --tblout '{result_file_name}' '{profile_file_name}' '{fasta_sequences_file}'"
                subprocess.run(search_command, shell=True, check=True)
                os.remove(profile_file_name)
            
                df = pd.read_csv(result_file_name, sep='\s+', skiprows=3, header=None)
                df = df.iloc[:, [0, 4, 5]]
                df = df.iloc[:-10]
                df.columns = ['ID', 'E-value', 'Score']
                df.to_csv(self.score_output_file, sep='\t', index=False)
                
                top_sequences = df.head(int(self.amount))
                records_to_write = []
                    
                for ind, row in top_sequences.iterrows():
                    with open(self.sequences_file, "r") as input_handle:
                            for record in SeqIO.parse(input_handle, "fasta"):
                                    if record.id == row['ID']:
                                         records_to_write.append(record)

                with open(self.new_fasta, "w") as output_handle:
                    SeqIO.write(records_to_write, output_handle, "fasta")


            return self.score_output_file
    
    
    def dna_to_protein(self, file_path):
        sequences = SeqIO.parse(file_path, 'fasta')
        new_records = []

        for record in sequences:
            # Переводим нуклеотидную последовательность в белковую
            dna_sequence = Seq(str(record.seq))
            protein_sequence = dna_sequence.translate()

            # Заменяем нуклеотидную последовательность белковой
            record.seq = protein_sequence
            new_records.append(record)

        # Записываем измененные записи обратно в файл
        with open(file_path, 'w') as out_handle:
            SeqIO.write(new_records, out_handle, 'fasta')
            
                    
            
    def run_BLAST(self):

        if self.sequences_type == 'DNA':
            self.dna_to_protein(self.new_fasta)

        if self.algorithm_blast == 'psiblast':
            self.psiblast_run()
        elif self.algorithm_blast == 'blastp':
            self.blastp_run()
            

        try: 
            df = pd.read_csv(self.result_blast, delimiter='\t', header=None, names=['ID', 'Sequence', 'Description', 'E-Value', 'Score'])
            df['Description'] = df['Description'].str.rsplit('=', n=1).str[-1]
            json_data = df.to_json(orient='records')
            df.to_excel(self.output_blast_file, index=False, header=True)

            return self.output_blast_file, json_data
        except: return "Error"

        
    
    def psiblast_run(self):
        command = f'psiblast -query "{self.new_fasta}" -db {DB_PROTEIN_FOR_BLAST}swissprot -outfmt "6 qseqid sseqid salltitles evalue bitscore" -out {self.result_blast} -max_target_seqs {self.num_blast_results}'
        subprocess.run(command, shell=True)
        print('PsiBLAST finished')
          

    
    def blastp_run(self):
        command = f'blastp -query "{self.new_fasta}" -db {DB_PROTEIN_FOR_BLAST}swissprot -outfmt "6 qseqid sseqid salltitles evalue bitscore" -out {self.result_blast} -max_target_seqs {self.num_blast_results}'
        subprocess.run(command, shell=True)
        print('BLASTp finished')
     

        

      
    
    
    
