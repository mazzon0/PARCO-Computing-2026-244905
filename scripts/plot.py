import pandas as pd
import matplotlib.pyplot as plt
import os

# Configuration
GAPBS_CSV = 'gapbs.csv'
OMP_CSV = 'omp.csv'
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
df_omp = load_and_prepare(OMP_CSV, 'Custom-OMP')

# Combine them
df = pd.concat([df_gapbs, df_omp], ignore_index=True)
df = df.sort_values('threads')

# 2. Calculate Speedup
# We calculate speedup relative to the 1-thread time of the SAME implementation and dataset
for (impl, ds), group in df.groupby(['impl', 'dataset']):
    t1_row = group[group['threads'] == 1]
    if not t1_row.empty:
        t1 = t1_row['execution time'].values[0]
        df.loc[group.index, 'speedup'] = t1 / group['execution time']

# 3. Plotting
fig, axes = plt.subplots(2, 2, figsize=(16, 12))
fig.suptitle('Performance Comparison: GAPBS vs Custom OMP', fontsize=20)

# We want to compare datasets, so we'll use different colors/markers for Dataset+Impl
metrics = [
    ('execution time', 'Execution Time (s)', 'Time Comparison'),
    ('speedup', 'Speedup Factor', 'Strong Scaling Speedup'),
    ('miss_rate', 'Cache Miss Rate (%)', 'Cache Efficiency'),
    ('ipc', 'IPC (Insn/Cycle)', 'Instruction Throughput')
]

# Helper to plot each subplot
for i, (col, ylabel, title) in enumerate(metrics):
    ax = axes[i//2, i%2]
    
    # Group by Dataset and Implementation
    # We'll filter for a specific affinity (e.g., spread-threads) to keep the graph clean,
    # or you can plot everything. Here we plot all unique pairs.
    for (ds, impl, bind), group in df.groupby(['dataset', 'impl', 'bind']):
        label = f"{impl} - {ds} ({bind})"
        # Use solid lines for GAPBS, dashed for OMP
        linestyle = '-' if impl == 'GAPBS' else '--'
        ax.plot(group['threads'], group[col], marker='o', label=label, linestyle=linestyle)

    if col == 'speedup':
        ax.plot(df['threads'].unique(), df['threads'].unique(), color='black', linestyle=':', label='Ideal')

    ax.set_title(title, fontsize=14)
    ax.set_ylabel(ylabel)
    ax.set_xlabel('Threads')
    ax.grid(True, alpha=0.3)
    ax.legend(prop={'size': 8})

plt.tight_layout(rect=[0, 0.03, 1, 0.95])
os.makedirs('results', exist_ok=True)
plt.savefig(OUTPATH, dpi=300)
print(f"Comparison plot saved to {OUTPATH}")