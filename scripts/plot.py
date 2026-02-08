import pandas as pd
import matplotlib.pyplot as plt
import os

# Create results directory if it doesn't exist
os.makedirs('results', exist_ok=True)

def load_and_prepare_data():
    files = {
        'seq': 'results/seq.csv',
        'omp': 'results/omp.csv',
        'omp2': 'results/omp2.csv',
        'gapbs': 'results/gapbs.csv'
    }
    
    dfs = {}
    for name, path in files.items():
        if os.path.exists(path):
            df = pd.read_csv(path)
            
            if 'dataset' in df.columns:
                # Standardize to lowercase and remove extensions
                df['dataset'] = (df['dataset'].astype(str)
                                 .str.replace(r'\.(csr|el)$', '', regex=True)
                                 .str.strip().str.lower())

            if 'threads' in df.columns:
                df['threads'] = pd.to_numeric(df['threads'], errors='coerce').fillna(1).astype(int)

            cols_to_fix = ['cache ref', 'cache miss', 'cycles', 'instructions', 'average time', 'median time']
            for col in cols_to_fix:
                if col in df.columns:
                    df[col] = pd.to_numeric(df[col], errors='coerce')
            
            if 'instructions' in df.columns and 'cycles' in df.columns:
                df['ipc'] = df['instructions'] / df['cycles']
            if 'cache miss' in df.columns and 'cache ref' in df.columns:
                df['miss_rate'] = (df['cache miss'] / df['cache ref']) * 100
            
            dfs[name] = df
            
    return dfs

def plot_per_method(dfs):
    methods = [m for m in ['omp', 'omp2', 'gapbs'] if m in dfs]
    
    # Define the same metrics used in the dataset plots
    metrics = [
        ('median time', 'Execution Time (s)'),
        ('speedup', 'Speedup'),
        ('instructions', 'Instructions'),
        ('cycles', 'Clock Cycles'),
        ('cache ref', 'Cache References'),
        ('cache miss', 'Cache Misses'),
        ('ipc', 'Instructions Per Clock (IPC)'),
        ('miss_rate', 'Cache Miss Rate (%)')
    ]

    for method in methods:
        df = dfs[method]
        datasets = sorted(df['dataset'].unique())
        
        # Create a 4x2 grid, similar to plot_per_dataset
        fig, axes = plt.subplots(4, 2, figsize=(15, 22))
        fig.suptitle(f'Performance Overview: {method.upper()}', fontsize=18, fontweight='bold')
        axes = axes.flatten()
        
        plot_happened = False

        for i, (metric, ylabel) in enumerate(metrics):
            ax = axes[i]
            
            for ds in datasets:
                subset = df[df['dataset'] == ds].sort_values('threads')
                if subset.empty: continue
                
                if metric == 'speedup':
                    # Calculate self-relative speedup for each dataset
                    t1_row = subset[subset['threads'] == 1]
                    if not t1_row.empty:
                        base_time = t1_row['median time'].values[0]
                        if base_time > 0:
                            y = base_time / subset['median time']
                            ax.plot(subset['threads'], y, marker='o', label=ds)
                            plot_happened = True
                else:
                    if metric in subset.columns:
                        ax.plot(subset['threads'], subset[metric], marker='o', label=ds)
                        plot_happened = True

            # Add ideal line specifically to the speedup subplot
            if metric == 'speedup':
                max_t = df['threads'].max()
                ax.plot([1, max_t], [1, max_t], color='red', linestyle='--', alpha=0.5, label='Ideal')

            ax.set_title(ylabel, fontweight='bold')
            ax.set_xlabel('Threads')
            ax.grid(True, linestyle='--', alpha=0.6)
            
            # Add legend to every subplot
            handles, labels = ax.get_legend_handles_labels()
            if handles:
                ax.legend(fontsize='x-small', loc='best')

        if plot_happened:
            plt.tight_layout(rect=[0, 0.03, 1, 0.96])
            out_path = f'results/{method}.png'
            plt.savefig(out_path)
            print(f"Saved Method Multi-plot: {out_path}")
        
        plt.close()

def plot_per_dataset(dfs):
    # Ensure seq_baseline is calculated within the scope
    seq_baseline = {}
    if 'seq' in dfs:
        seq_baseline = dfs['seq'].groupby('dataset')['median time'].median().to_dict()

    datasets = set()
    for df in dfs.values():
        if 'dataset' in df.columns:
            datasets.update(df['dataset'].unique())

    methods_to_compare = {'omp': 'omp', 'omp2': 'omp2', 'gapbs': 'gapbs'}
    metrics = [
        ('median time', 'Execution Time (s)'),
        ('speedup', 'Speedup'),
        ('instructions', 'Instructions'),
        ('cycles', 'Clock Cycles'),
        ('cache ref', 'Cache References'),
        ('cache miss', 'Cache Misses'),
        ('ipc', 'Instructions Per Clock (IPC)'),
        ('miss_rate', 'Cache Miss Rate (%)')
    ]

    for ds in datasets:
        print(f"Processing dataset: {ds}...")
        fig, axes = plt.subplots(4, 2, figsize=(15, 20))
        fig.suptitle(f'Metrics for Dataset: {ds}', fontsize=16)
        axes = axes.flatten()
        plot_has_any_data = False

        for i, (metric, ylabel) in enumerate(metrics):
            ax = axes[i]
            for method_key, label in methods_to_compare.items():
                if method_key not in dfs: continue
                
                df_m = dfs[method_key]
                mask = (df_m['dataset'] == ds)
                
                # Check for env filters
                if 'bind' in df_m.columns and 'spread' in df_m['bind'].values:
                    mask &= (df_m['bind'] == 'spread')
                
                subset = df_m[mask].sort_values('threads')
                if subset.empty: continue

                if metric == 'speedup':
                    t1_row = subset[subset['threads'] == 1]
                    t1_time = None
                    label_suffix = ""
                    
                    if not t1_row.empty:
                        t1_time = t1_row['median time'].values[0]
                        label_suffix = "(Self)"
                    elif ds in seq_baseline:
                        t1_time = seq_baseline[ds]
                        label_suffix = "(vs Seq)"

                    if t1_time:
                        y = t1_time / subset['median time']
                        ax.plot(subset['threads'], y, marker='o', label=f"{label} {label_suffix}")
                        plot_has_any_data = True
                else:
                    ax.plot(subset['threads'], subset[metric], marker='o', label=label)
                    plot_has_any_data = True
            
            ax.set_title(ylabel)
            ax.set_xlabel('Threads')
            ax.grid(True, linestyle='--', alpha=0.6)
            ax.legend()

        if plot_has_any_data:
            plt.tight_layout(rect=[0, 0.03, 1, 0.95])
            plt.savefig(f'results/{ds}.png')
            print(f"Saved Dataset Plot: results/{ds}.png")
        plt.close()

if __name__ == "__main__":
    data_frames = load_and_prepare_data()
    if data_frames:
        plot_per_dataset(data_frames)
        plot_per_method(data_frames)