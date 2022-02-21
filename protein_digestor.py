#Tryptic peptide generator
#Logic credit: Santosh Kumar Behera, beheras40@gmail.com
#Code credit: Sandeep Kasaragod, sandeep.kolya@gmail.com

    
def Trypsin(sequence, missed_clevage, pep_min_len, pep_max_len):
    if 'K' in sequence or 'R' in sequence:
        get_dup_k = [i for i in range(len(sequence)) if sequence.startswith('K', i)]
        get_dup_r = [j for j in range(len(sequence)) if sequence.startswith('R', j)]
        merge_list = sorted(get_dup_k + get_dup_r)
        merge_list_fltrd = [i for i in merge_list if i+1 < len(sequence) and sequence[i + 1] !='P'] #look for KP or RP position
        merge_list_fltrd.append(len(sequence))
        initialize = 0
        for iter_lst in range(len(merge_list_fltrd) - int(missed_clevage)):
            peptide = (sequence[initialize: int(merge_list_fltrd[iter_lst + missed_clevage]) + 1])
            if len(peptide) >= int(pep_min_len) and len(peptide) <= int(pep_max_len):
                yield peptide
            initialize = merge_list_fltrd[iter_lst] + 1

def Lysc(sequence, missed_clevage, pep_min_len, pep_max_len):
    if 'K' in sequence:
        get_dup_k = [i for i in range(len(sequence)) if sequence.startswith('K', i)]
        merge_list_fltrd = get_dup_k
        merge_list_fltrd.append(len(sequence))
        initialize = 0
        for iter_lst in range(len(merge_list_fltrd) - int(missed_clevage)):
            peptide = (sequence[initialize: int(merge_list_fltrd[iter_lst + missed_clevage]) + 1])
            if len(peptide) >= int(pep_min_len) and len(peptide) <= int(pep_max_len):
                yield peptide
            initialize = merge_list_fltrd[iter_lst] + 1
        
def Chymotrypsin(sequence, missed_clevage, pep_min_len, pep_max_len):
    if 'F' in sequence or 'W' in sequence or 'Y' in sequence:
        get_dup_f = [i for i in range(len(sequence)) if sequence.startswith('F', i)]
        get_dup_w = [j for j in range(len(sequence)) if sequence.startswith('W', j)]
        get_dup_y = [j for j in range(len(sequence)) if sequence.startswith('Y', j)]
        merge_list = sorted(get_dup_f + get_dup_w + get_dup_y)
        merge_list_fltrd = [i for i in merge_list if i+1 < len(sequence) and sequence[i + 1] !='P'] #look for KP or RP position
        merge_list_fltrd.append(len(sequence))
        initialize = 0
        for iter_lst in range(len(merge_list_fltrd) - int(missed_clevage)):
            peptide = (sequence[initialize: int(merge_list_fltrd[iter_lst + missed_clevage]) + 1])
            if len(peptide) >= int(pep_min_len) and len(peptide) <= int(pep_max_len):
                yield peptide
            initialize = merge_list_fltrd[iter_lst] + 1
        
##        
##a = tryptic_peptide_trypsin("MNWTVDIPIDQLPPLPPLPGDLRARLDAALAKPAAQQPSWPTDQATAMRTVLESVPPVTVPSEIVRLQEQLAQVARGEAFLLQGGDCAETFTENTEPHIRGNIRTLLQMAVVLTYGASMPVVKVARIAGQYAKPRSADIDALGLKSYRGDMINGFAPNAAVREHDPSRLVRAYANASAAMNLVRALTSSGLASLHLVHDWNREFVRTSPAGARYEALATEIDRGLRFMSACGVADRNLQTAEIYASHEALVLDYERAMLRLSDSDAGEPQLYDLSAHTVWIGERTRQLDGAHIAFAEVIANPIGVKIGPTITPELAVEYVERLDPHNKPGRLTLVSRMGNSKVRDLLPPIVEKVQATGHQVIWQCDPMHGNTHESSNGYKTRHFDRIVDEVQGFFEVHRALGTHPGGIHVEITGDNVTECLGGAQDISDMDLVGRYETACDPRLNTQQSLELAFLVAEMLRD", 1, 7, 50)
##for i in a:
##    print(i)

        
    
    
