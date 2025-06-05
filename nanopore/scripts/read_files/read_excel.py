import pandas as pd

def write_fasta_from_excel(excel_file, output_file, sheet_name=0):
    df = pd.read_excel(excel_file, sheet_name=sheet_name)
    df = df.dropna(subset=['sample_id', 'sequence'])

    with open(output_file, 'w') as fasta:
        for _, row in df.iterrows():
            sample_id = str(row['sample_id']).strip()
            sequence = str(row['sequence']).strip()
            sequence = str(row['sequence']).strip().replace(' ', '')  
            fasta.write(f'>{sample_id}\n{sequence}\n')

if __name__ == "__main__":
    input_excel = rf'D:\nanopore\idealDNA.xlsx' 
    output_fasta = rf'D:\nanopore\data\fasta\idealDNA.fasta'  

    write_fasta_from_excel(input_excel, output_fasta)
    print("FASTA has generated")