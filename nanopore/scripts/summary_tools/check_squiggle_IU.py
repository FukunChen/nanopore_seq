import os
import json
from collections import defaultdict

# path of json file
json_folder = r"D:\nanopore\data\aligned_reference"  

# get all the bases
u_i_matches = defaultdict(list)  # {"U": [base1, base2, ...], "I": [...]}

for filename in os.listdir(json_folder):
    if filename.endswith(".json"):
        filepath = os.path.join(json_folder, filename)

        with open(filepath, "r") as f:
            try:
                data = json.load(f)

                if isinstance(data, list):
                    records = data
                elif isinstance(data, dict):
                    records = [data]
                else:
                    print(f"{filename} is wrong")
                    continue

                for entry in records:
                    q_aln = entry["query_aln"]
                    r_aln = entry["ref_aln"]

                    for q_base, r_base in zip(q_aln, r_aln):
                        if r_base in ("U", "I") and q_base != "-":
                            u_i_matches[r_base].append(q_base)

            except Exception as e:
                print(f"{filename} error is: {e}")

        for q_base, r_base in zip(q_aln, r_aln):
            if r_base in ("U", "I") and q_base != "-":
                u_i_matches[r_base].append(q_base)

for base_type in ['U', 'I']:
    print(f"\n=== {base_type} Matches (total num: {len(u_i_matches[base_type])}) ===")
    for b in sorted(set(u_i_matches[base_type])):
        count = u_i_matches[base_type].count(b)
        print(f"{b}: {count}")