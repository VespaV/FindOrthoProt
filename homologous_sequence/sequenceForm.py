from django import forms
from django.core.validators import MinValueValidator

class SequenceForm(forms.Form):
   sequences_type = forms.CharField(required=True, initial="new_type", widget=forms.HiddenInput())
   homologous_sequences_fasta_file = forms.FileField(required=True, label="1. Файл с гомологичными последовательностями", help_text="Аминокислотные или нуклеотидные последовательности в формате FASTA")
   annotated_sequences_fasta_file = forms.FileField(required=True, label="2. Файл  с аннатированными последовательностями", help_text="Аминокислотные или нуклеотидные последовательности в формате FASTA")
   selected_sequences_count = forms.IntegerField(required=True, validators=[MinValueValidator(0)], min_value=1, max_value=99, label="3. Количество последовательностей для анализа BLAST", initial=5)
   blastp_results_count = forms.IntegerField(required=True,  min_value=1, max_value=99, label="4. Количесво результатов BLAST", initial=5)
   algorithm = forms.ChoiceField(choices=(("hmmer3", "С помощью построения HMM3-профиля"),("pair_aligment", "Нидлмана — Вунша (попарное выравнивание с консенсусом)")), required=True, label="5. Алгоритм поиска гомологичных последовательностей")
   algorithm_blast = forms.ChoiceField(choices=(("psiblast", "PSI-BLAST"),("blastp", "BLASTp")), required=True, label="6. Алгоритм BLAST")