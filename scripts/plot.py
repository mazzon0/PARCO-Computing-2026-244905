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
            
            # Clean dataset names (remove .csr, .el, or whitespace)
            if 'dataset' in df.columns:
                df['dataset'] = df['dataset'].astype(str).str.replace(r'\.(csr|el)$', '', regex=True).str.strip()

            # Ensure numeric types
            cols_to_fix = ['cache ref', 'cache miss', 'cycles', 'instructions', 'average time', 'median time', 'threads']
            for col in cols_to_fix:
                if col in df.columns:
                    df[col] = pd.to_numeric(df[col], errors='coerce')
            
            # Derived metrics
            if 'instructions' in df.columns and 'cycles' in df.columns:
                df['ipc'] = df['instructions'] / df['cycles']
            if 'cache miss' in df.columns and 'cache ref' in df.columns:
                df['miss_rate'] = (df['cache miss'] / df['cache ref']) * 100
            
            dfs[name] = df
            
    return dfs

def plot_per_dataset(dfs):
    if 'seq' not in dfs:
        print("Missing seq.csv. Speedup cannot be calculated.")
        return

    # Create baseline dictionary: { 'dataset_name': median_time_float }
    seq_df = dfs['seq']
    seq_baseline = seq_df.set_index('dataset')['median time'].to_dict()
    
    datasets = set()
    for df in dfs.values():
        datasets.update(df['dataset'].unique())

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

    methods_to_compare = {
        'omp': 'omp (spread threads)',
        'omp2': 'omp2',
        'gapbs': 'gapbs'
    }

    for ds in datasets:
        fig, axes = plt.subplots(4, 2, figsize=(15, 20))
        fig.suptitle(f'Metrics for Dataset: {ds}', fontsize=16)
        axes = axes.flatten()

        for i, (metric, ylabel) in enumerate(metrics):
            ax = axes[i]
            for method_key, label in methods_to_compare.items():
                if method_key not in dfs: continue
                
                df_method = dfs[method_key]
                subset = df_method[
                    (df_method['dataset'] == ds) & 
                    (df_method['bind'] == 'spread') & 
                    (df_method['place'] == 'threads')
                ].sort_values('threads')
                
                if subset.empty: continue
                
                x = subset['threads']
                if metric == 'speedup':
                    base_time = seq_baseline.get(ds)
                    # Check if base_time exists and is a valid number
                    if base_time and pd.notnull(base_time) and base_time > 0:
                        y = base_time / subset['median time']
                    else:
                        continue # Skip if no sequential baseline for this dataset
                else:
                    y = subset[metric]
                
                ax.plot(x, y, marker='o', label=label)
            
            ax.set_title(ylabel)
            ax.set_xlabel('Threads')
            ax.set_ylabel(ylabel)
            ax.legend()
            ax.grid(True, linestyle='--', alpha=0.7)

        plt.tight_layout(rect=[0, 0.03, 1, 0.95])
        plt.savefig(f'results/{ds}.png')
        print(f'Saved results/{ds}.png')
        plt.close()

# Note: plot_per_method should be updated with the same dataset cleaning and speedup logic
def plot_per_method(dfs):
    seq_baseline = dfs['seq'].set_index('dataset')['median time'].to_dict() if 'seq' in dfs else {}

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

    for method_name, df_method in dfs.items():
        datasets = df_method['dataset'].unique()
        fig, axes = plt.subplots(4, 2, figsize=(15, 20))
        fig.suptitle(f'Metrics for Method: {method_name}', fontsize=16)
        axes = axes.flatten()

        for i, (metric, ylabel) in enumerate(metrics):
            ax = axes[i]
            for ds in datasets:
                if method_name == 'seq':
                    subset = df_method[df_method['dataset'] == ds]
                    x = [1]
                else:
                    subset = df_method[
                        (df_method['dataset'] == ds) & 
                        (df_method['bind'] == 'spread') & 
                        (df_method['place'] == 'threads')
                    ].sort_values('threads')
                    x = subset['threads']
                
                if subset.empty: continue

                if metric == 'speedup':
                    base_time = seq_baseline.get(ds)
                    if method_name == 'seq':
                        y = [1.0] * len(subset)
                    elif base_time and pd.notnull(base_time):
                        y = base_time / subset['median time']
                    else:
                        continue
                else:
                    y = subset[metric]
                
                ax.plot(x, y, marker='s', label=ds)
            
            ax.set_title(ylabel)
            ax.set_xlabel('Threads' if method_name != 'seq' else 'N/A')
            ax.set_ylabel(ylabel)
            ax.legend()
            ax.grid(True, linestyle='--', alpha=0.7)

        plt.tight_layout(rect=[0, 0.03, 1, 0.95])
        plt.savefig(f'results/{method_name}.png')
        print(f'Saved results/{method_name}.png')
        plt.close()

if __name__ == "__main__":
    data_frames = load_and_prepare_data()
    if data_frames:
        plot_per_dataset(data_frames)
        plot_per_method(data_frames)