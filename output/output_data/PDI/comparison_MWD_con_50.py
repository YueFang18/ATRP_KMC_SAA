import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
from numpy import trapz # Explicitly import trapz for numerical integration

def read_chain_lengths(filepath):
    """
    Reads chain length data from the given file and returns a numpy array.
    Assumes numbers separated by whitespace. Includes error handling.
    """
    if not os.path.exists(filepath):
        print(f"Error: File not found at {filepath}")
        return None
    try:
        with open(filepath, 'r') as file:
            data = file.read().split()
            lengths = [int(length) for length in data if length.strip().isdigit()]
        if not lengths:
             print(f"Warning: No numeric data found in {filepath}")
             return np.array([])
        return np.array(lengths)
    except ValueError as e:
        print(f"Error converting data to integers in {filepath}: {e}")
        return None
    except Exception as e:
        print(f"An unexpected error occurred while reading {filepath}: {e}")
        return None

def process_and_get_mwd_data(filepath, label_prefix):
    """
    Reads, filters chain lengths, calculates stats, performs KDE,
    calculates weight fraction density, and returns plot data.

    Args:
        filepath (str): Full path to the data file.
        label_prefix (str): Prefix for the legend label (e.g., '1e5/2').

    Returns:
        dict: Contains 'x_vals', 'y_vals' (weight fraction density), 'label',
              or None if processing fails.
    """
    print(f"Processing: {filepath}...")
    chain_lengths = read_chain_lengths(filepath)

    if chain_lengths is None:
        print(f"  Failed to read or process file: {filepath}")
        return None

    if len(chain_lengths) == 0:
        print("  File is empty or contains no valid numeric chain lengths.")
        return None

    total_chains = len(chain_lengths)
    filtered_chain_lengths = chain_lengths[chain_lengths >= 3]
    removed_chains = total_chains - len(filtered_chain_lengths)

    if len(filtered_chain_lengths) == 0:
        print(f"  Original chains: {total_chains}")
        print("  No valid chain lengths (>=3) after filtering.")
        return None

    sum_lengths = np.sum(filtered_chain_lengths)
    if sum_lengths == 0:
        print(f"  Sum of filtered chain lengths is zero. Cannot calculate Mw.")
        return None

    Mn = np.mean(filtered_chain_lengths)
    Mw = np.sum(filtered_chain_lengths.astype(np.float64)**2) / sum_lengths
    PDI = Mw / Mn

    print(f"  File: {os.path.basename(filepath)} (from {label_prefix})")
    print(f"    Original chains: {total_chains}")
    print(f"    Filtered chains (>=3): {len(filtered_chain_lengths)}")
    print(f"    Removed chains (<3): {removed_chains}")
    print(f"    Mn: {Mn:.2f}")
    print(f"    Mw: {Mw:.2f}")
    print(f"    PDI: {PDI:.2f}\n")

    try:
        # Perform KDE for number density
        kde = gaussian_kde(filtered_chain_lengths)
        min_val = max(0, np.min(filtered_chain_lengths) - 0.1 * Mn)
        max_val = np.max(filtered_chain_lengths) + 0.1 * Mn
        x_vals = np.linspace(min_val, max_val, 500)
        kde_vals = kde(x_vals) # This is proportional to number density f_N(DP)

        # --- Calculate Weight Fraction Density ---
        # 1. Unnormalized weight density = DP * f_N(DP)
        # Use np.maximum to avoid potential issues if x_vals starts exactly at 0
        # Although kde_vals should be ~0 there anyway.
        w_unnormalized = np.maximum(x_vals, 1e-9) * kde_vals

        # 2. Calculate normalization factor (area under unnormalized curve)
        norm_factor = trapz(w_unnormalized, x_vals)

        # 3. Normalize to get weight fraction density
        if norm_factor <= 1e-9: # Check for very small or zero area
            print(f"  Warning: Normalization factor for weight fraction is near zero ({norm_factor}). Check KDE output. Setting weight fraction to zero.")
            w_vals_normalized = np.zeros_like(w_unnormalized)
        else:
            w_vals_normalized = w_unnormalized / norm_factor
        # --- End of Weight Fraction Calculation ---

    except Exception as e:
        print(f"  Error during KDE or Weight Fraction calculation for {filepath}: {e}")
        return None

    plot_label = f'{label_prefix} (Mn={Mn:.1f}, PDI={PDI:.2f})'

    return {
        "x_vals": x_vals,
        "y_vals": w_vals_normalized, # Return normalized weight fraction density
        "label": plot_label,
        "raw_data": filtered_chain_lengths, # Return the filtered data for histogram
        "stats": {
            "Mn": Mn,
            "Mw": Mw,
            "PDI": PDI
        }
    }

def calculate_weight_fraction_histogram(data, bins=25):
    """
    手动计算加权的直方图，将数目分数转换为重量分数
    
    Args:
        data (np.array): 聚合度数据
        bins (int): 直方图的bin数量
    
    Returns:
        tuple: (bin_edges, weight_fraction_values, bin_centers)
    """
    # 计算直方图的bin和计数
    counts, bin_edges = np.histogram(data, bins=bins)
    
    # 计算bin的中心位置
    bin_centers = 0.5 * (bin_edges[1:] + bin_edges[:-1])
    
    # 将每个bin的计数乘以该bin的中心DP值，得到重量分布
    weight_counts = counts * bin_centers
    
    # 计算总重量用于归一化
    total_weight = np.sum(weight_counts)
    
    if total_weight <= 0:
        print("Warning: Total weight is zero or negative. Cannot normalize.")
        return bin_edges, np.zeros_like(counts), bin_centers
    
    # 归一化为重量分数
    bin_widths = np.diff(bin_edges)
    weight_fraction = weight_counts / (total_weight * bin_widths)  # 除以bin宽度使得总面积为1
    
    return bin_edges, weight_fraction, bin_centers

def create_plot_for_folder(folder_name):
    """
    为指定文件夹创建MWD图表
    
    Args:
        folder_name (str): 子文件夹名称 (500, 1000, 1500, 2500)
    
    Returns:
        bool: 是否成功创建图表
    """
    datasets_info = [
        {'base_dir': f'1e5/{folder_name}', 'label_prefix': '1e5', 'style': 'histogram', 'color': 'lightgreen'},
        {'base_dir': f'1e7/{folder_name}', 'label_prefix': '1e7', 'style': 'line', 'color': 'green'},
    ]
    common_filename = 'all_chain_conv_50.out'

    plt.figure(figsize=(6, 4.2), dpi=200)
    plot_successful_count = 0
    
    # Process all datasets first to determine common x-axis range
    processed_datasets = []
    
    for info in datasets_info:
        filepath = os.path.join(info['base_dir'], common_filename)
        processed_data = process_and_get_mwd_data(filepath, info['label_prefix'])
        
        if processed_data:
            processed_data['style'] = info['style']
            processed_data['color'] = info['color']
            processed_datasets.append(processed_data)
    
    # Now plot each dataset with appropriate style
    for data in processed_datasets:
        if data['style'] == 'line':
            # Plot data as solid line (using KDE weight fraction)
            plt.plot(
                data["x_vals"],
                data["y_vals"],
                linewidth=2,
                label=data["label"],
                color=data['color']
            )
        elif data['style'] == 'histogram':
            # 计算加权的直方图（weight fraction）
            num_bins = 25  # 调整数量
            bin_edges, weight_frac_values, bin_centers = calculate_weight_fraction_histogram(
                data['raw_data'], bins=num_bins
            )
            
            # 使用plt.bar绘制加权直方图
            plt.bar(
                bin_centers,
                weight_frac_values,
                width=np.diff(bin_edges),  # 使用完整宽度，柱状图之间没有空隙
                alpha=0.6,
                label=data["label"],
                color=data['color']
            )
        
        plot_successful_count += 1

    if plot_successful_count > 0:
        plt.xlabel('Degree of Polymerization', fontsize=12)
        plt.ylabel('Weight Fraction', fontsize=12)
        plt.grid(True, linestyle='--', alpha=0.5)
        
        # 获取当前Y轴范围
        y_min, y_max = plt.ylim()
        # 增加Y轴上限，留出空间放置图例
        plt.ylim(bottom=0, top=y_max*1.2)  # 上限增加20%
        
        # 将图例放在右上角
        plt.legend(fontsize=10, loc='upper right')
        
        # 在标题中加入子文件夹信息
        plt.title(f'Molecular Weight Distribution - {folder_name}', fontsize=14)

        output_filename = f'Combined_Conv50_MWD_WeightFraction-{folder_name}.png'
        plt.savefig(output_filename, dpi=300, bbox_inches='tight')
        print(f"\nCombined plot for {folder_name} saved as {output_filename}")
        plt.close()  # 关闭图表避免内存泄漏
        return True
    else:
        print(f"\nNo data was successfully processed and plotted for {folder_name}.")
        plt.close()  # 关闭图表避免内存泄漏
        return False

def main():
    # 要处理的所有子文件夹
    folders = ['2', '4', '8', '20']
    successful_plots = 0
    
    for folder in folders:
        print(f"\n{'='*50}")
        print(f"Processing folder: {folder}")
        print(f"{'='*50}")
        
        if create_plot_for_folder(folder):
            successful_plots += 1
    
    print(f"\nSummary: Successfully created {successful_plots} out of {len(folders)} plots.")

if __name__ == '__main__':
    main()
