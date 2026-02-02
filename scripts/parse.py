import re
import csv
import os

def parse_perf_logs(input_files):
    # Data containers
    gapbs_data = []
    omp_data = []
    seq_data = []

    # Regex to match headers: supports GAPBS, OMP, and SEQ
    # Group 1: Label, Group 2: bind, Group 3: place, Group 4: threads, Group 5: dataset
    header_re = re.compile(r"(GAPBS|OMP|SEQ)\s+(?:bind=(\w+)\s+place=(\w+)\s+)?threads=(\d+)\s+dataset=(\w+)")

    for file_path in input_files:
        if not os.path.exists(file_path):
            print(f"Warning: {file_path} not found. Skipping.")
            continue

        with open(file_path, 'r') as f:
            content = f.read()

        # Split file by our custom headers
        blocks = re.split(r"(GAPBS|OMP|SEQ)", content)
        
        # re.split with capturing groups keeps the delimiter. 
        # Blocks will look like: ['', 'GAPBS', ' rest of block...', 'OMP', ' rest...']
        for i in range(1, len(blocks), 2):
            label = blocks[i]
            block_content = blocks[i+1]
            
            # Extract header info
            header_match = header_re.search(label + block_content)
            if not header_match:
                continue
                
            label, bind, place, threads, dataset = header_match.groups()

            # Helper to extract perf values and remove commas
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
            elif label == "OMP":
                omp_data.append(row)
            elif label == "SEQ":
                seq_data.append(row)

    return gapbs_data, omp_data, seq_data

def save_to_csv(data, filename):
    if not data:
        print(f"No data found for {filename}. Skipping.")
        return
    
    headers = ["dataset", "place", "bind", "threads", "execution time", "cache misses", "cache refs", "instructions", "cycles"]
    
    with open(filename, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=headers)
        writer.writeheader()
        writer.writerows(data)
    print(f"Successfully created: {filename}")


input_logs = ["gapbs_stats.txt", "custom_stats.txt"]
gapbs, omp, seq = parse_perf_logs(input_logs)

save_to_csv(gapbs, "gapbs.csv")
save_to_csv(omp, "omp.csv")
save_to_csv(seq, "seq.csv")
