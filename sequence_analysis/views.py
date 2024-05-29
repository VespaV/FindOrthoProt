from django.shortcuts import render
from Bio import SeqIO
from sequence_analysis.constants import SELECTED_PREDICTED_PROTEINS
from homologous_sequence.constants import FASTA_OUTPUT_FOLDER, HOMOLOGOUS_FOLDER
import os
import uuid
from FindOrthoProt.generate_file_name import generate_id_file
import json
from django.http import HttpResponse

def analysis_seq_view(request):
    homologous_seq = None
    select_predicted_seq = None

    if request.method == 'GET':
        data = request.GET.dict()
        print("-------- data ------------")
        print(data)
        print("-------- data ------------")


        homologous_seq = data['homologous_seq'] if 'homologous_seq' in data else None
        selected_ids = data['selectedIDs'].split(',') if 'selectedIDs' in data else None
        predicted_seq = data['predicted_seq'] if 'predicted_seq' in data else None

        if predicted_seq is not None:
            select_predicted_seq = filter_fasta_by_ids(FASTA_OUTPUT_FOLDER+predicted_seq, selected_ids)
        else: select_predicted_seq = None


    return render(request, 'form_analysis.html', {'homologous_seq': homologous_seq, 'select_predicted_seq': select_predicted_seq})


def analysis_success_view(request):
    homologous_seq = None
    select_predicted_seq = None
    if request.method == 'POST':
        data = request.POST.dict()
        print("-------- data ------------")
        print(data)
        print("-------- data ------------")
        
        homologous_sequences_fasta_file = request.FILES.get('homologous_seq', False)
        if not homologous_sequences_fasta_file: 
            homologous_sequences_filename = request.POST.get('homologous_seq')
        else:
            homologous_sequences_filename = generate_id_file(homologous_sequences_fasta_file, 'txt' ,HOMOLOGOUS_FOLDER)
            with open(homologous_sequences_filename, 'wb+') as destination:
                for chunk in homologous_sequences_fasta_file.chunks():
                    destination.write(chunk)
        
        select_predicted_seq = request.FILES.get('select_predicted_seq', False)

        if not select_predicted_seq: 
            select_predicted_seq_filename = request.POST.get('select_predicted_seq')
        else:
            select_predicted_seq_filename = generate_id_file(select_predicted_seq, 'txt' ,SELECTED_PREDICTED_PROTEINS)
            with open(select_predicted_seq_filename, 'wb+') as destination:
                for chunk in select_predicted_seq.chunks():
                    destination.write(chunk)
        
        chem_phys_properties = request.POST.get('chemPhysProperties', 'off')
        search_motifs = request.POST.get('searchMotifs', 'off')
        num_motifs = request.POST.get('numMotifs', None)
        build_phylo_tree = request.POST.get('buildPhyloTree', 'off')
        alignment_algorithm = request.POST.get('alignmentAlgorithm', None)
        phylo_algorithm = request.POST.get('phyloAlgoritm', None)
        
        # Convert string representations to actual boolean values
   
        form_data = {
            'homologous_sequences_fasta_file': homologous_sequences_filename,
            'select_predicted_seq': select_predicted_seq_filename,
            'chem_phys_properties': chem_phys_properties,
            'search_motifs': search_motifs,
            'num_motifs': num_motifs,
            'build_phylo_tree': build_phylo_tree,
            'alignment_algorithm': alignment_algorithm,
            'phylo_algorithm': phylo_algorithm,
        }
        
        print('view form_data', form_data)

        return render(request, 'success_analysis.html',  context=form_data)
    return render(request,
                  'form_analysis.html', {'homologous_seq': homologous_seq, 'select_predicted_seq': select_predicted_seq})

def filter_fasta_by_ids(input_fasta, target_ids):
    records_to_write = []
    output_fasta = generate_id_file('autopredict', 'txt', SELECTED_PREDICTED_PROTEINS)

    with open(input_fasta, 'r') as input_file:
        for record in SeqIO.parse(input_file, 'fasta'):
            if record.id in target_ids:
                records_to_write.append(record)

    with open(output_fasta, 'w') as output_file:
        SeqIO.write(records_to_write, output_file, 'fasta')

    return output_fasta

def meme_view(request):
    # Получаем значение параметра html_link из запроса
    html_link = request.GET.get('html_link')

    # Проверяем, существует ли файл по указанному пути
    if html_link and os.path.exists(html_link):
        # Открываем файл и передаем его содержимое в контекст представления
        with open(html_link, 'r') as file:
            html_content = file.read()

        # Возвращаем содержимое файла в виде HTTP-ответа
        return HttpResponse(html_content, content_type='text/html')
    else:
        # Если файл не найден, возвращаем ошибку или редирект
        return HttpResponse('File not found', status=404)



