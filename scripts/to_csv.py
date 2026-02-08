import re
import csv
import statistics
import os

def parse_stats(input_files):
    # Prepare data containers
    data = {
        'seq': [],
        'omp': [],
        'gapbs': [],
        'omp2': []
    }

    for filename in input_files:
        if not os.path.exists(filename):
            print(f"Warning: {filename} not found. Skipping.")
            continue
            
        with open(filename, 'r') as f:
            content = f.read()
            # Split by sections (looking for the header tags like SEQ, OMP, etc.)
            sections = re.split(r'\n(?=SEQ|OMP|GAPBS|OMP2)', content)
            
            for section in sections:
                if not section.strip(): continue
                
                # Extract Header Info
                header_match = re.search(r'(SEQ|OMP2|OMP|GAPBS)\s+(?:bind=(\S+)\s+place=(\S+)\s+)?threads=(\d+)\s+dataset=(\S+)', section)
                if not header_match: continue
                
                type_tag = header_match.group(1).lower()
                bind = header_match.group(2)
                place = header_match.group(3)
                threads = header_match.group(4)
                dataset = header_match.group(5)

                # Extract Performance Counters (Handling commas and spaces)
                cache_ref = re.search(r'([\d,]+)\s+cache-references', section)
                cache_miss = re.search(r'([\d,]+)\s+cache-misses', section)
                cycles = re.search(r'([\d,]+)\s+cycles', section)
                instructions = re.search(r'([\d,]+)\s+instructions', section)

                def clean_num(match):
                    return match.group(1).replace(',', '') if match else ""

                # Extract Timing
                avg_time = ""
                med_time = ""
                
                if type_tag == 'gapbs':
                    # Calculate median from Trial Times
                    trials = re.findall(r'Trial Time:\s+([\d.]+)', section)
                    if trials:
                        trial_floats = [float(t) for t in trials]
                        med_time = statistics.median(trial_floats)
                    
                    avg_match = re.search(r'Average Time:\s+([\d.]+)', section)
                    avg_time = avg_match.group(1) if avg_match else ""
                else:
                    # OMP, OMP2, and SEQ usually have this line (SEQ might be missing it in some logs)
                    time_match = re.search(r'Average time:\s+([\d.]+)\s+-\s+Median time:\s+([\d.]+)', section)
                    if time_match:
                        avg_time = time_match.group(1)
                        med_time = time_match.group(2)

                row = {
                    'dataset': dataset,
                    'threads': threads,
                    'bind': bind,
                    'place': place,
                    'cache ref': clean_num(cache_ref),
                    'cache miss': clean_num(cache_miss),
                    'cycles': clean_num(cycles),
                    'instructions': clean_num(instructions),
                    'average time': avg_time,
                    'median time': med_time
                }
                
                data[type_tag].append(row)

    return data

def save_csv(filename, data, fields):
    os.makedirs(os.path.dirname(filename), exist_ok=True)
    with open(filename, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for row in data:
            # Filter row to only include relevant fields for this specific CSV
            filtered_row = {k: row[k] for k in fields if k in row}
            writer.writerow(filtered_row)

if __name__ == "__main__":
    log_files = ['results/custom_stats.txt', 'results/gapbs_stats.txt', 'results/omp2_stats.txt']
    parsed_results = parse_stats(log_files)

    # Field configurations
    seq_fields = ['dataset', 'cache ref', 'cache miss', 'cycles', 'instructions', 'average time', 'median time']
    parallel_fields = ['dataset', 'threads', 'bind', 'place', 'cache ref', 'cache miss', 'cycles', 'instructions', 'average time', 'median time']

    # Write outputs
    save_csv('results/seq.csv', parsed_results['seq'], seq_fields)
    save_csv('results/omp.csv', parsed_results['omp'], parallel_fields)
    save_csv('results/gapbs.csv', parsed_results['gapbs'], parallel_fields)
    save_csv('results/omp2.csv', parsed_results['omp2'], parallel_fields)

    print("CSV files generated successfully in the 'result/' directory.")