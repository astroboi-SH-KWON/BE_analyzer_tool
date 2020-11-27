# BE_analyzer_tool

input :
    result1 = joined fastq files by ea-utils
    ref_seq = ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    trgt_seq = ********************


1. ea-utils 사용하여 Fastq join
    설치파일) https://anaconda.org/bioconda/ea-utils
    result1

2. Reference sequence에서 cleavage point를 찾고나서, 양 말단에 있는 15nt indicator sequence를 찾기
     위 indicator 포함하는 것만 모아서 valid sequence로 count. (1nt mismatch까지 허용)
    
    1. find indicators from ref_seq
     indi1_seq, indi2_seq

    2. find trgt_seq in fq_read ( = FASTQ_read)
        fq_read = ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
        ::::::::::::::::::::::::::::::::::::::********************::::::::::::::::::::::::::::::::::::::::::::::::
        
    3. find PAM (must be 'NGG')     
        ::::::::::::::::::::::::::::::::::::::********************PAM:::::::::::::::::::::::::::::::::::::::::::::
        
    4. find both indicators (15 bp) from fq_read ( cnt_mismatch <= 1 )
        :::::::::::::<<<indicator>>>::::::::::********************PAM::::::::::<<<indicator>>>::::::::::::::::::::
        
    5. get seq from fq_read with both 20 extra bp from each indicators     
        :::::::::20bp<<<indicator>>>::::::::::********************PAM::::::::::<<<indicator>>>20bp::::::::::::::::
        
        result2 = 20bp<<<indicator>>>::::::::::********************PAM::::::::::<<<indicator>>>20bp
    

3. EMBOSS needle 이용하여 reference sequence에 align  Insertion, deletion, substitution 으로 분류

    1. align ref_seq, result2
        ref_seq = 20bp<<<indicator>>>||||||||||********************PAM||||||||||<<<indicator>>>20bp
        result2 = 20bp<<<indicator>>>::::::::::********************PAM::::::::::<<<indicator>>>20bp
    
