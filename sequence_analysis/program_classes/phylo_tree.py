import subprocess
import json
from Bio import Phylo
from Bio import SeqIO
from Bio import AlignIO
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio.Phylo.TreeConstruction import DistanceCalculator
from sequence_analysis.constants import PHYLO_ALIGNMENTS, PHYLO_JOINED_FASTA, PHYLO_TREE_NEWIK, PHYLO_TREE_PICTURE
from FindOrthoProt.generate_file_name import generate_id_file
import os
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')

class PhyloTree:
    def __init__(self, homologous_fasta, predicted_fasta, algoritm_multipl_align, algoritm_phylo):
        self.homologous_fasta = homologous_fasta
        self.predicted_fasta = predicted_fasta
        self.joined_fasta = generate_id_file('joined', 'fasta', PHYLO_JOINED_FASTA)
        self.algoritm_multipl_align = algoritm_multipl_align
        self.algoritm_phylo = algoritm_phylo
        self.output_align_file = generate_id_file('align', 'fasta', PHYLO_ALIGNMENTS)
        self.output_newick_file = generate_id_file('tree', 'nwk', PHYLO_TREE_NEWIK)
        self.output_image_tree = generate_id_file('tree', 'png', PHYLO_TREE_PICTURE)

    def join_fasta(self):
        with open(self.homologous_fasta, 'r') as homologous_file, open(self.predicted_fasta, 'r') as predicted_file, open(self.joined_fasta, 'w') as joined_file:
            for line in homologous_file:
                joined_file.write(line)
            
            joined_file.write('\n')

            for record in SeqIO.parse(predicted_file, "fasta"):
                    record.id += '_predicted'
                    SeqIO.write(record, joined_file, "fasta")
                 
    def make_align(self):
        self.join_fasta()
        if self.algoritm_multipl_align == "muscle":
            return 'multiAlign', self._muscle_align_run()
        elif self.algoritm_multipl_align == "clustal_o":
            return 'multiAlign', self.clustal_o_run()
        elif self.algoritm_multipl_align == "mafft":
            return 'multiAlign', self.mafft_run()
        else: 
            return 'Wrong name of aligner'

    def build_tree(self):
        if self.algoritm_phylo == 'FastTree':
            newik = self.fast_tree_build()
            return 'forTree', newik
        elif self.algoritm_phylo == 'RAxML':
            newik = self.raxml_build()
            return 'forTree', newik
        elif self.algoritm_phylo == 'MaximumParsimony':
            newik = self.maximum_parsimony_build()
            return 'forTree', newik
        else: 
            return 'Wrong name of phylology algorithm'
        
    def _muscle_align_run(self):
        command = ["muscle3", "-in", self.joined_fasta, "-out", self.output_align_file]
        subprocess.run(command)
        print(f"Multiple alignment has been performed and saved to {self.output_align_file}")
        return self.output_align_file
        
    def clustal_o_run(self):
        command = ["clustalo", "-i", self.joined_fasta, "--out", self.output_align_file]
        subprocess.run(command)
        print(f"Clustal Omega alignment has been performed and saved to {self.output_align_file}")
        return self.output_align_file
        
    def mafft_run(self):
        command = f"mafft {self.joined_fasta} > {self.output_align_file}"
        subprocess.run(command, shell=True, executable='/bin/bash', check=True, text=True)
        print(f"MAFFT alignment has been performed and saved to {self.output_align_file}")
        return self.output_align_file
    
    def fast_tree_build(self):
        command = f"FastTree {self.output_align_file} > {self.output_newick_file} "
        subprocess.run(command, shell=True, executable='/bin/bash', check=True, text=True)
        print(f"FastTree has been performed and saved to {self.output_newick_file}")
        return self.output_newick_file
        
    def raxml_build(self):
        command = f"raxml-ng --msa '{self.output_align_file}' --model LG+G8+F --prefix '{self.output_newick_file}'"
        subprocess.run(command, shell=True, executable='/bin/bash', check=True, text=True)
        print(f"Raxml has been performed and saved to {self.output_newick_file}")
        
        file_path1 = self.output_newick_file + '.raxml.log'
        file_path2 = self.output_newick_file + '.raxml.rba'
        file_path3 = self.output_newick_file + '.raxml.bestModel'
        file_path4 = self.output_newick_file + '.raxml.startTree'
        file_path5 = self.output_newick_file + '.raxml.bestTreeCollapsed'
        file_path6 = self.output_newick_file + '.raxml.mlTrees'
        file_path7 = self.output_newick_file + '.raxml.reduced.phy'

        if os.path.exists(file_path1):
            os.remove(file_path1)
        if os.path.exists(file_path2):
            os.remove(file_path2)
        if os.path.exists(file_path3):
            os.remove(file_path3)
        if os.path.exists(file_path4):
            os.remove(file_path4)
        if os.path.exists(file_path5):
            os.remove(file_path5)
        if os.path.exists(file_path6):
            os.remove(file_path6)
        if os.path.exists(file_path7):
            os.remove(file_path7)
        
        path_to_newick = self.output_newick_file + '.raxml.bestTree'
        
        os.rename(path_to_newick, self.output_newick_file)
        print(f'Файл успешно переименован в {self.output_newick_file}.')
        return self.output_newick_file
        
    def maximum_parsimony_build(self):
        alignment = AlignIO.read(self.output_align_file, "fasta")
        calculator = DistanceCalculator('identity')
        dist_matrix = calculator.get_distance(alignment)
        
        constructor = DistanceTreeConstructor()
        starting_tree = constructor.nj(dist_matrix)
        
        scorer = Phylo.TreeConstruction.ParsimonyScorer()
        searcher = Phylo.TreeConstruction.NNITreeSearcher(scorer)
        constructor = Phylo.TreeConstruction.ParsimonyTreeConstructor(searcher, starting_tree)
        parsimony_tree = constructor.build_tree(alignment)
        
        Phylo.write(parsimony_tree, self.output_newick_file, "newick")
        print('File was safed as', self.output_newick_file)

        return self.output_newick_file

    def visual_tree(self):
        tree = Phylo.read(self.output_newick_file, 'newick')
        num_leaves = len(tree.get_terminals())

        fig_width = 20 + 0.5 * num_leaves
        fig_height = 8 + 0.3 * num_leaves

        fig, ax = plt.subplots(figsize=(fig_width, fig_height))

        ax.set_xlabel('X-Label', fontsize=22)  # Установка подписи оси X с размером шрифта 18
        ax.set_ylabel('Y-Label', fontsize=22)

        tree.ladderize()

        # Устанавливаем размер шрифта для листьев
        font = {'weight': 'bold',
                'size': 16}
        plt.rc('font', **font)
        plt.tight_layout()

        Phylo.draw(tree, axes=ax)

        plt.savefig(self.output_image_tree)
        print('Tree pic was saved as', self.output_image_tree)
        plt.close()

        return "Tree visualisation", self.output_image_tree
    


