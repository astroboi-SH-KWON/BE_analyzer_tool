import time
import os
import platform

import Util
import Logic
############### st env ################
WORK_DIR = os.getcwd() + "/"
PROJECT_NAME = WORK_DIR.split("/")[-2]
SYSTEM_NM = platform.system()

if SYSTEM_NM == 'Linux':
    # REAL
    pass
else:
    # DEV
    WORK_DIR = "D:/000_WORK/JangHyeWon/20201116_BE_analyzer_web_tool/WORK_DIR/"

IN = 'input/'
OU = 'output/'
FASTQ = 'FASTQ/'

TRGT_SEQ_PATH = "pegRNA target DNA seq.txt"
REF_SEQ_PATH = "Reference seq_FAH.txt"

LEN_WIN = 20
LEN_INDI = 15
LEN_BTWN = 10
LEN_PAM = 3
PAM = 'NGG'
############### en env ################


def main():
    util = Util.Utils()
    logic = Logic.Logics()

    ref_list = util.read_tsv_ignore_N_line(WORK_DIR + IN + REF_SEQ_PATH, 0)
    trgt_list = util.read_tsv_ignore_N_line(WORK_DIR + IN + TRGT_SEQ_PATH, 0)

    ref_seq = ref_list[0][0]
    trgt_seq = trgt_list[0][0]
    st_idx_trgt_ref = ref_seq.index(trgt_seq)
    en_idx_trgt_ref = st_idx_trgt_ref + len(trgt_seq)
    indi1_seq = ref_seq[st_idx_trgt_ref - LEN_BTWN - LEN_INDI: st_idx_trgt_ref - LEN_BTWN]
    indi2_seq = ref_seq[en_idx_trgt_ref + LEN_PAM + LEN_BTWN: en_idx_trgt_ref + LEN_PAM + LEN_BTWN + LEN_INDI]

    fastq_read_list = util.read_FASTQ_join_by_biopython(WORK_DIR + FASTQ + 'result.join.fq')

    filtered_seq_list = []
    for fq_read in fastq_read_list:
        if trgt_seq in fq_read:
            st_idx_trgt_fastq = fq_read.index(trgt_seq)
            en_idx_trgt_fastq = st_idx_trgt_fastq + len(trgt_seq)

            # check PAM
            if not logic.match(0, fq_read[en_idx_trgt_fastq: en_idx_trgt_fastq + LEN_PAM], PAM):
                continue

            indi1_fastq = fq_read[st_idx_trgt_fastq - LEN_BTWN - LEN_INDI: st_idx_trgt_fastq - LEN_BTWN]
            indi2_fastq = fq_read[
                          en_idx_trgt_fastq + LEN_PAM + LEN_BTWN: en_idx_trgt_fastq + LEN_PAM + LEN_BTWN + LEN_INDI]
            if len(indi1_fastq) != LEN_INDI or len(indi2_fastq) != LEN_INDI:
                print(fq_read)
                continue
            if logic.cnt_mismatch(indi1_fastq, indi1_seq) <= 1 and logic.cnt_mismatch(indi2_fastq, indi2_seq) <= 1:
                filtered_seq_list.append(fq_read[st_idx_trgt_fastq - LEN_BTWN - LEN_INDI - LEN_WIN: en_idx_trgt_fastq + LEN_PAM + LEN_BTWN + LEN_INDI + LEN_WIN])

    align_result_list = []
    ref_seq = ref_seq[st_idx_trgt_ref - LEN_BTWN - LEN_INDI - LEN_WIN: en_idx_trgt_ref + LEN_PAM + LEN_BTWN + LEN_INDI + LEN_WIN]
    for filtered_seq in filtered_seq_list:
        seq_a, align, seq_b, _ = logic.get_pairwise2_globalds_result(ref_seq, filtered_seq)
        align_result_list.append([seq_a, align, seq_b])

    header = ['ref_seq', 'align', 'fastq_read']
    util.make_excel(WORK_DIR + OU + 'aligned_result', header, align_result_list)


if __name__ == '__main__':
    start_time = time.perf_counter()
    print("start [ " + PROJECT_NAME + " ]>>>>>>>>>>>>>>>>>>")
    main()
    print("::::::::::: %.2f seconds ::::::::::::::" % (time.perf_counter() - start_time))