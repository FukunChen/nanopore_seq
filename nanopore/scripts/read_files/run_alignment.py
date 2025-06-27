import parasail
import numpy as np
sub_matrix = parasail.Matrix(r"D:\nanopore\scripts\alignment\didu_matrix_withIU.txt")

def run_alignment(query, ref, match_score, mismatch_score):
    result = parasail.sg_qx_trace_striped_16(query, ref, 10, 5, sub_matrix)
    tb = result.traceback
    
    q_aln = tb.query
    r_aln = tb.ref
    comp = tb.comp
    current_score = result.score

    print("\nAlignment Result:")
    print("Query : ", q_aln)
    print("Comp  : ", comp)
    print("Ref   : ", r_aln)
    print("Score    : ", current_score)

    # 统计符号数量
    comp_array = np.array(list(comp))
    
    for symbol in ['|', ':', '.', ' ']:
        count = np.sum(comp_array == symbol)
        print(f"Symbol '{symbol}': {count} times")
        

# 示例序列
query = "ACGTGCAGTCCCATTTACGACTGAGCTAGCGATTAGCCACA"
ref   = "AGGTACTACNNGUTNNAGC"

run_alignment(query, ref, match_score=4, mismatch_score=-5)

