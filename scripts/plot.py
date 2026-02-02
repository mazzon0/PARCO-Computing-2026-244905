import pandas as pd
import matplotlib.pyplot as plt
import os

# Configuration
GAPBS_CSV = 'results/gapbs.csv'
OMP_CSV = 'results/omp.csv'
OMP2_CSV = 'results/omp2.csv'
SEQ_CSV = 'results/seq.csv'
OUTPATH = 'results/comparison_plots.png'

def load_and_prepare(file_path, label):
    if not os.path.exists(file_path):
        print(f"Warning: {file_path} not found.")
        return pd.DataFrame()
    
    df = pd.read_csv(file_path)
    df['impl'] = label
    # Calculate derived metrics
    df['miss_rate'] = (df['cache misses'] / df['cache refs']) * 100
    df['ipc'] = df['instructions'] / df['cycles']
    return df

# 1. Load Data
df_gapbs = load_and_prepare(GAPBS_CSV, 'GAPBS')
df_omp = load_and_prepare(OMP_CSV, 'OMP-Standard')
df_omp2 = load_and_prepare(OMP2_CSV, 'OMP-Fused')
df_seq = load_and_prepare(SEQ_CSV, 'Sequential')

# Combine them (excluding Sequential for the main trend lines, we'll use it as T1)
df_parallel = pd.concat([df_gapbs, df_omp, df_omp2], ignore_index=True)
df_parallel = df_parallel.sort_values('threads')

# 2. Calculate Speedup relative to the Sequential version (the "Gold Standard")
for (ds), group in df_parallel.groupby(['dataset']):
    # Get sequential time for this specific dataset
    seq_match = df_seq[df_seq['dataset'] == ds]
    if not seq_match.empty:
        t_seq = seq_match['execution time'].values[0]
        df_parallel.loc[group.index, 'speedup'] = t_seq / group['execution time']
    else:
        # Fallback: Speedup relative to own 1-thread if SEQ is missing
        for impl, sub_group in group.groupby('impl'):
            t1 = sub_group[sub_group['threads'] == 1]['execution time'].min()
            df_parallel.loc[sub_group.index, 'speedup'] = t1 / sub_group['execution time']

# 3. Plotting
fig, axes = plt.subplots(2, 2, figsize=(18, 14))
fig.suptitle('PageRank Benchmark: GAPBS vs Custom OMP vs Fused OMP', fontsize=22)

metrics = [
    ('execution time', 'Execution Time (s)', 'Time Comparison (Lower is Better)'),
    ('speedup', 'Speedup Factor', 'Strong Scaling (Higher is Better)'),
    ('miss_rate', 'Cache Miss Rate (%)', 'Memory Pressure (Lower is Better)'),
    ('ipc', 'IPC (Insn/Cycle)', 'Instruction Throughput (Higher is Better)')
]

# Distinct styles for different implementations
styles = {
    'GAPBS': {'ls': '-', 'marker': 'o', 'alpha': 0.8},
    'OMP-Standard': {'ls': '--', 'marker': 's', 'alpha': 0.6},
    'OMP-Fused': {'ls': '-.', 'marker': '^', 'alpha': 0.9}
}

for i, (col, ylabel, title) in enumerate(metrics):
    ax = axes[i//2, i%2]
    
    for (ds, impl, bind), group in df_parallel.groupby(['dataset', 'impl', 'bind']):
        # Only plot 'spread' or 'close' to avoid cluttering, or plot all with labels
        style = styles.get(impl, {'ls': '-', 'marker': 'x', 'alpha': 0.5})
        label = f"{impl} | {ds} ({bind})"
        
        ax.plot(group['threads'], group[col], label=label, 
                linestyle=style['ls'], marker=style['marker'], alpha=style['alpha'])

    if col == 'speedup':
        ax.plot(df_parallel['threads'].unique(), df_parallel['threads'].unique(), 
                color='black', linestyle=':', label='Ideal Linear Speedup')
    
    ax.set_title(title, fontsize=15, fontweight='bold')
    ax.set_ylabel(ylabel)
    ax.set_xlabel('Thread Count')
    ax.set_xscale('log', base=2) # Logical for 1, 2, 4, 8, 16, 32, 64
    ax.set_xticks(df_parallel['threads'].unique())
    ax.get_xaxis().set_major_formatter(plt.ScalarFormatter())
    ax.grid(True, which="both", ls="-", alpha=0.2)
    ax.legend(prop={'size': 7}, ncol=2)

plt.tight_layout(rect=[0, 0.03, 1, 0.95])
plt.savefig(OUTPATH, dpi=300)
print(f"Refined comparison plots saved to {OUTPATH}")