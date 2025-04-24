import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde

# —— 学术图表参数 —— #
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.size': 10,
    'axes.labelsize': 12,
    'axes.titlesize': 14,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'legend.fontsize': 8,
    'axes.linewidth': 1.2,
    'grid.linewidth': 0.8,
    'grid.alpha': 0.25,
    'lines.linewidth': 1.8,
    'xtick.direction': 'in',
    'ytick.direction': 'in',
    'xtick.top': True,
    'ytick.right': True
})

def read_chain_lengths(filepath: str) -> np.ndarray:
    """Read chain-length list from plain-text file."""
    with open(filepath, 'r') as f:
        data = f.read().split()
    return np.array([int(x) for x in data])

# 转化率采样点
conv_points = ['10', '15', '20', '25', '30', '35', '40', '45', '50']
colors = plt.cm.viridis(np.linspace(0, 1, len(conv_points)))

# 顶层文件夹（指数级别）
bases = ['1e5', '1e7']
# 需要遍历的子目录
subfolders = ['2', '4', '8', '20']

for sub in subfolders:
    # 创建一个 2×1 的共享 X 轴图
    fig, axes = plt.subplots(2, 1, figsize=(6, 8), sharex=True)
    
    for ax, base in zip(axes, bases):
        folder_path = os.path.join(base, sub)  # e.g. "1e5/500"
        
        for conv, color in zip(conv_points, colors):
            fname = os.path.join(folder_path, f'all_chain_conv_{conv}.out')
            if not os.path.exists(fname):
                print(f'[Warning] File not found: {fname}')
                continue
            
            lengths = read_chain_lengths(fname)
            if lengths.size == 0:
                print(f'[Warning] Empty data in {fname}')
                continue
            
            Mn = np.mean(lengths)
            Mw = np.sum(lengths**2) / np.sum(lengths)
            PDI = Mw / Mn

            kde = gaussian_kde(lengths, weights=lengths)
            x = np.linspace(lengths.min(), lengths.max(), 500)
            y = kde(x)
            y /= np.trapz(y, x)  # 归一化至面积 1

            ax.plot(x, y, color=color,
                    label=fr'{conv}%: $M_n$={Mn:.1f}, PDI={PDI:.2f}')
        
        ax.set_ylabel('Weight Fraction')
        ax.legend(loc='upper right', frameon=True)
        ax.set_title(f'{base} – {sub} times')

    axes[-1].set_xlabel('Degree of Polymerization')
    plt.tight_layout(pad=2.0)

    # 保存图像
    out_name = f'MWD_1e5_vs_1e7_{sub} times.png'
    plt.savefig(out_name, dpi=300)
    print(f'[Saved] {out_name}')
    plt.close(fig)
