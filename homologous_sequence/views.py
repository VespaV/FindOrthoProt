from django.shortcuts import render
from homologous_sequence.sequenceForm import SequenceForm
from homologous_sequence.constants import HOMOLOGOUS_FOLDER, ANNOTATED_FOLDER
from FindOrthoProt.generate_file_name import generate_id_file


def homologous_sequence_view(request):
  if request.method == 'POST':

    form = SequenceForm(request.POST, request.FILES)
    if form.is_valid():
      print("Cleaned Data:", form.cleaned_data)

      sequences_type = form.cleaned_data['sequences_type']
      if sequences_type[0] == 'Format error' or sequences_type[1] == 'type_error':
           return render(request, 'homologous_sequence.html', {'form': form})
          
      homologous_sequences_fasta_file= form.cleaned_data['homologous_sequences_fasta_file']
      homologous_sequences_filename = generate_id_file(homologous_sequences_fasta_file, 'fasta', HOMOLOGOUS_FOLDER)

      with open(homologous_sequences_filename, 'wb+') as destination:
                for chunk in homologous_sequences_fasta_file.chunks():
                    destination.write(chunk)

      annotated_sequences_fasta_file = form.cleaned_data['annotated_sequences_fasta_file']
      annotated_sequences_filename = generate_id_file(annotated_sequences_fasta_file, 'fasta', ANNOTATED_FOLDER)
      
      with open(annotated_sequences_filename, 'wb+') as destination:
                for chunk in annotated_sequences_fasta_file.chunks():
                    destination.write(chunk)
                    
                    
     
      selected_sequences_count = form.cleaned_data['selected_sequences_count']
      blastp_results_count = form.cleaned_data['blastp_results_count']
      algorithm = form.cleaned_data['algorithm']
      algorithm_blast = form.cleaned_data['algorithm_blast']
      
      return render(request, 'success_template.html',
                    {'sequences_type':sequences_type,
                    'homologous_sequences_fasta_file': homologous_sequences_filename, 
                    'annotated_sequences_fasta_file': annotated_sequences_filename, 
                    'selected_sequences_count' : selected_sequences_count, 
                    'blastp_results_count': blastp_results_count,
                    'algorithm': algorithm, 
                    'algorithm_blast': algorithm_blast})
  else:
    form = SequenceForm()

  return render(request,
          'homologous_sequence.html',
          {'form': form})


