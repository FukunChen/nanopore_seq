input_bed = r"D:\nanopore\data\bed\didumix.bed"
output_txt = "uuid_list.txt"

with open(input_bed, "r") as bed, open(output_txt, "w") as out:
    for line in bed:
        if line.strip():  # 跳过空行
            uuid = line.split()[0]  # 提取第 1 列
            out.write(uuid + "\n")

print(f"UUID 已提取并保存到 {output_txt}")
