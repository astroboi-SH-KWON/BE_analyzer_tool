from Bio import SeqIO


from astroboi_bio_tools.ToolUtil import ToolUtils
class Utils(ToolUtils):
    def read_FASTQ_join_by_biopython(self, path):
        temp = list(SeqIO.parse(path, 'fastq'))
        return [str(temp[k].seq) for k in range(len(temp))]