import re
import csv
import os

def parse_perf_logs(input_files):
    # Data containers
    gapbs_data = []
    omp_data = []
    omp2_data = [] # New container for fused version
    seq_data = []

    # Updated Regex: (GAPBS|OMP|OMP2|SEQ) to catch the new label
    header_re = re.compile(r"(GAPBS|OMP2|OMP|SEQ)\s+(?:bind=(\w+)\s+place=(\w+)\s+)?threads=(\d+)\s+dataset=(\w+)")

    for file_path in input_files:
        if not os.path.exists(file_path):
            print(f"Warning: {file_path} not found. Skipping.")
            continue

        with open(file_path, 'r') as f:
            content = f.read()

        # Split file by our custom headers (including OMP2)
        blocks = re.split(r"(GAPBS|OMP2|OMP|SEQ)", content)
        
        for i in range(1, len(blocks), 2):
            label = blocks[i]
            block_content = blocks[i+1]
            
            header_match = header_re.search(label + block_content)
            if not header_match:
                continue
                
            label, bind, place, threads, dataset = header_match.groups()

            def extract_val(pattern):
                m = re.search(pattern, block_content)
                return m.group(1).replace(',', '') if m else ""

            row = {
                "dataset": dataset,
                "place": place if place else "",
                "bind": bind if bind else "",
                "threads": threads,
                "execution time": extract_val(r"([\d.]+)\s+seconds time elapsed"),
                "cache misses": extract_val(r"([\d,]+)\s+cache-misses"),
                "cache refs": extract_val(r"([\d,]+)\s+cache-references"),
                "instructions": extract_val(r"([\d,]+)\s+instructions"),
                "cycles": extract_val(r"([\d,]+)\s+cycles")
            }

            # Route to the correct list
            if label == "GAPBS":
                gapbs_data.append(row)
            elif label == "OMP2":
                omp2_data.append(row)
            elif label == "OMP":
                omp_data.append(row)
            elif label == "SEQ":
                seq_data.append(row)

    return gapbs_data, omp_data, omp2_data, seq_data

def save_to_csv(data, filename):
    if not data:
        print(f"No data found for {filename}. Skipping.")
        return
    
    os.makedirs(os.path.dirname(filename), exist_ok=True)
    headers = ["dataset", "place", "bind", "threads", "execution time", "cache misses", "cache refs", "instructions", "cycles"]
    
    with open(filename, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=headers)
        writer.writeheader()
        writer.writerows(data)
    print(f"Successfully created: {filename}")

# --- Execution ---
input_logs = ["gapbs_stats.txt", "custom_stats.txt"]
gapbs, omp, omp2, seq = parse_perf_logs(input_logs)

save_to_csv(gapbs, "results/gapbs.csv")
save_to_csv(omp, "results/omp.csv")
save_to_csv(omp2, "results/omp2.csv")
save_to_csv(seq, "results/seq.csv")