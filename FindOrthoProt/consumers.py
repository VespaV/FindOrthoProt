import json
import threading

import asyncio
from channels.generic.websocket import AsyncWebsocketConsumer
from homologous_sequence.program_classes.check_homology import Homology
from sequence_analysis.program_classes.prot_param import ProteinPropertiesAnalyzer
from sequence_analysis.program_classes.phylo_tree import PhyloTree

class MyConsumer(AsyncWebsocketConsumer):
    async def connect(self):
        await self.accept()
        self.data_lock = threading.Lock()
        self.received_data = None

    async def receive(self, text_data):
        data = json.loads(text_data)
        if data['type'] == 'send_form_data':
            await self.handle_form_data(data['data'])
        elif data['type'] == 'form_seq_analysis':
            await self.handle_form_seq_analysis(data['data'])

    async def handle_form_data(self, form_data):
        print('Form data received on the server:', form_data)
        print(type(form_data))

        with self.data_lock:
            data_for_program = Homology(sequences_type = form_data['sequences_type'],
                                        search_seq_file = form_data['homologous_sequences_fasta_file'],
                                        sequences_file = form_data['annotated_sequences_fasta_file'],
                                        amount = form_data['selected_sequences_count'],
                                        num_blast_results = form_data['blastp_results_count'],
                                        algorithm = form_data['algorithm'],
                                        algorithm_blast = form_data['algorithm_blast']
                                        )

            message1 = await self.my_function1(data_for_program)
            await self.send_message_to_browser(message1)

            message2 = await self.my_function2(data_for_program)
            await self.send_message_to_browser(message2)

            message3 = await self.my_function3(data_for_program)
            await self.send_message_to_browser(message3)


    async def my_function1(self, data):
        answer = data.run_muscle()
        await asyncio.sleep(2)
        return answer 


    async def my_function2(self, data):
        answer = data.search_homologous()
        await asyncio.sleep(3)
        return answer

    async def my_function3(self, data):
        answer = data.run_BLAST()
        await asyncio.sleep(4)
        return answer
       
    async def handle_form_seq_analysis(self, form_data):
        print('Form data received on the server 2:', form_data)
        if form_data['chem_phys_properties'] == 'True':
            print('chem_phys_properties = True')
            prot_param = ProteinPropertiesAnalyzer(homologous_fasta = form_data['homologous_sequences_fasta_file'],
                                                   predicted_fasta = form_data['select_predicted_seq'],
                                                   )
            message1 = await self.prot_param_amino_count(prot_param)
            await self.send_message_to_browser(message1)

            message2 = await self.prot_param_phys_chem(prot_param)
            await self.send_message_to_browser(message2)

            message3 = await self.secondary_structure(prot_param)
            await self.send_message_to_browser(message3)

        if form_data['build_phylo_tree'] == 'True':
            print('phylo_algorithm = True')
            phylo_object = PhyloTree(homologous_fasta = form_data['homologous_sequences_fasta_file'],
                                    predicted_fasta = form_data['select_predicted_seq'],
                                    algoritm_multipl_align = form_data['alignment_algorithm'],
                                    algoritm_phylo = form_data['phylo_algorithm'])

            message1 = await self.align_for_tree(phylo_object)
            await self.send_message_to_browser(message1)

            message2 = await self.buid_tree(phylo_object)
            await self.send_message_to_browser(message2)

            massage3 = await self.visual_tree(phylo_object)
            await self.send_message_to_browser(massage3)

    async def prot_param_amino_count(self, prot_object):
        answer = prot_object.amino_gistogram()
        await asyncio.sleep(2)
        return answer

    async def prot_param_phys_chem(self, prot_object):
        answer = prot_object.all_params()
        await asyncio.sleep(3)
        return answer

    async def secondary_structure(self, prot_object):
        answer = prot_object.secondary_structure_gistogram()
        await asyncio.sleep(4)
        return answer

    async def align_for_tree(self, phylo_object):
        answer = phylo_object.make_align()
        await asyncio.sleep(1)
        return answer

    async def buid_tree(self, phylo_object):
        answer = phylo_object.build_tree()
        await asyncio.sleep(2)
        return answer

    async def visual_tree(self, phylo_object):
        answer = phylo_object.visual_tree()
        await asyncio.sleep(3)
        return answer



    async def send_message_to_browser(self, message):
            await self.send(text_data=json.dumps({'message': message}))
        
    
