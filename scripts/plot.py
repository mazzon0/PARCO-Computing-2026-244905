import re
import pandas as pd
import matplotlib.pyplot as plt

OUTPATH = 'results/gapbs_plots.png'

def parse_perf_log(log_text):
    data = []
    # Pattern to match the headers and extract config
    # Matches: GAPBS bind=spread place=threads threads=1 dataset=stanford
    config_pattern = r"GAPBS bind=(\w+) place=(\w+) threads=(\d+) dataset=(\w+)"
    
    # Split the log into blocks based on the GAPBS header
    blocks = re.split(r"(GAPBS bind=)", log_text)
    
    for i in range(1, len(blocks), 2):
        header_content = blocks[i] + blocks[i+1]
        
        # Extract metadata
        meta = re.search(config_pattern, header_content)
        if not meta: continue
        
        bind, place, threads, dataset = meta.groups()
        
        # Extract metrics using regex
        cache_misses = re.search(r"([\d,]+)\s+cache-misses", header_content)
        cache_refs = re.search(r"([\d,]+)\s+cache-references", header_content)
        cycles = re.search(r"([\d,]+)\s+cycles", header_content)
        instructions = re.search(r"([\d,]+)\s+instructions", header_content)
        elapsed = re.search(r"([\d.]+)\s+seconds time elapsed", header_content)
        
        # Clean numeric data (remove commas)
        def clean(val): return float(val.group(1).replace(',', '')) if val else None

        c_miss = clean(cache_misses)
        c_ref = clean(cache_refs)
        inst = clean(instructions)
        cyc = clean(cycles)
        time = clean(elapsed)
        
        data.append({
            'bind': bind,
            'place': place,
            'threads': int(threads),
            'dataset': dataset,
            'time': time,
            'miss_rate': (c_miss / c_ref) * 100 if c_miss and c_ref else 0,
            'ipc': inst / cyc if inst and cyc else 0
        })
        
    return pd.DataFrame(data)

# --- Main Script ---
log_data = open('results/stats.txt').read()

df = parse_perf_log(log_data)

df = df.sort_values('threads')
for (ds, bnd), group in df.groupby(['dataset', 'bind']):
    t1 = group[group['threads'] == 1]['time'].values[0]
    df.loc[group.index, 'speedup'] = t1 / group['time']

fig, axes = plt.subplots(2, 2, figsize=(14, 10))
fig.suptitle('GAPBS Performance Analysis', fontsize=16)

# Strong Scaling
for label, group in df.groupby(['dataset', 'bind']):
    axes[0,0].plot(group['threads'], group['time'], marker='o', label=f"{label[0]} ({label[1]})")
axes[0,0].set_title('Strong Scaling (Lower is Better)')
axes[0,0].set_ylabel('Seconds')
axes[0,0].legend()

# Speedup
for label, group in df.groupby(['dataset', 'bind']):
    axes[0,1].plot(group['threads'], group['speedup'], marker='s', label=f"{label[0]} ({label[1]})")
axes[0,1].plot(df['threads'].unique(), df['threads'].unique(), color='black', linestyle='--', label='Ideal')
axes[0,1].set_title('Strong Scaling Speedup')
axes[0,1].set_ylabel('Speedup Factor')
axes[0,1].legend()

# Cache Miss Rate
for label, group in df.groupby(['dataset', 'bind']):
    axes[1,0].plot(group['threads'], group['miss_rate'], marker='^', label=f"{label[0]} ({label[1]})")
axes[1,0].set_title('Cache Miss Rate')
axes[1,0].set_ylabel('Percentage (%)')
axes[1,0].legend()

# IPC
for label, group in df.groupby(['dataset', 'bind']):
    axes[1,1].plot(group['threads'], group['ipc'], marker='x', label=f"{label[0]} ({label[1]})")
axes[1,1].set_title('Instructions Per Cycle (IPC)')
axes[1,1].set_ylabel('IPC Value')
axes[1,1].legend()

plt.tight_layout(rect=[0, 0.03, 1, 0.95])
plt.savefig(OUTDIR, dpi=300)
plt.close()