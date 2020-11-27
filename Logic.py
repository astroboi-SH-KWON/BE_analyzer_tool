
from astroboi_bio_tools.ToolLogic import ToolLogics
class Logics(ToolLogics):
    def cnt_mismatch(self, fastq_read, rule_seq):
        cnt = 0
        for i in range(len(rule_seq)):
            if not self.checkSeqByChar(fastq_read[i], rule_seq[i]):
                cnt += 1
        return cnt
