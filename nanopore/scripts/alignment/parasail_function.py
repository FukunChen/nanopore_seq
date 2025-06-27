import os

try:
    import parasail
    # oligo repeat finding alignment constants
    GAP_OPEN_PENALTY = 8
    GAP_EXTEND_PENALTY = 4
    LOCAL_ALIGN_FUNCTION = parasail.sg_qx_trace_striped_16
    MATRIX = parasail.Matrix(r"D:\nanopore\scripts\alignment\didu_matrix_withIU.txt")
except ImportError:
    print('No alignment library')