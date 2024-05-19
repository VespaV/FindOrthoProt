from django.http import JsonResponse
from Bio import SeqIO
import os

def IUPAC_define(record):
    seq_alphabet = record.seq

    dna_chars= set("ACGTRYKMSWN")
    protein_chars = set("ACDEFGHIKLMNPQRSTVWYBXZJUO*")

    if set(seq_alphabet.upper()) <= dna_chars:
        return 'DNA'
    elif set(seq_alphabet.upper()) <= protein_chars:
        return 'Protein'
    else: 
        return 'type_error'


def check_format_view(request):     
    if request.method == 'POST':
        file1 = list(request.FILES.values())[0]
        file2 = list(request.FILES.values())[1]
        answer_type = None 
          
        try:
            file1_path = os.path.join("/tmp", file1.name)
            file2_path = os.path.join("/tmp",file2.name)

            with open(file1_path, 'wb') as homologous_tempfile:
                homologous_tempfile.write(file1.read())

            with open(file2_path, 'wb') as annotated_tempfile:
                annotated_tempfile.write(file2.read())
            
            if file1_path==file2_path:
                answer_format = 'Format error'

            file1_record = next(SeqIO.parse(file1_path, "fasta"))
            file2_record = next(SeqIO.parse(file2_path, "fasta"))

            type1 = IUPAC_define(file1_record)
            type2 = IUPAC_define(file2_record)
            
            if type1 == 'type_error' or type2 == 'type_error':
                answer_type = 'type_error'
            elif type1==type2:
                answer_type = type2
            else: 
                answer_type = 'type_error'
            
            answer_format = 'Correct fasta'   
            
        except:
            answer_format = 'Format error'
        finally:
            # Удаляем временные файлы
            os.remove(file1_path)
            os.remove(file2_path)
            
        response_data = {'format': answer_format, 'type': answer_type}
        return JsonResponse(response_data)
    
    else:
        return JsonResponse({'error': 'Метод не разрешен'})       
    
            