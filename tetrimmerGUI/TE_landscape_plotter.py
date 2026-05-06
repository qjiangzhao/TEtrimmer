import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches


def generate_unified_landscape_r_style(rm_output_file, max_div=40):
    # 1. Define column names
    cols = ['sw_score', 'perc_div', 'perc_del', 'perc_ins', 'query', 'q_start', 'q_end',
            'q_left', 'strand', 'repeat', 'class_family', 'r_start', 'r_end', 'r_left', 'id']

    # 2. Load the file
    df = pd.read_csv(
        rm_output_file,
        sep='\s+',
        skiprows=3,
        names=cols,
        engine='python',
        usecols=range(15)
    )

    # 3. Clean and Extract Family
    df['Length'] = df['q_end'] - df['q_start']
    df['Main_Class'] = df['class_family'].apply(lambda x: str(x).split('/')[-1])

    # 4. Binning
    df['div_bin'] = df['perc_div'].round(0).astype(int)

    # Aggregate data: Sum length by Bin and Family
    plot_data = df.groupby(['div_bin', 'Main_Class'])['Length'].sum().reset_index()
    plot_data = plot_data[plot_data['div_bin'] <= max_div]

    # 5. Pivot for Stacked Plotting
    # This transforms the data so each Main_Class is a column (required for stacked bar)
    pivot_df = plot_data.pivot(index='div_bin', columns='Main_Class', values='Length').fillna(0)

    # 6. Setup Plot Aesthetics (Mimicking theme_bw and ggsci)
    plt.style.use('seaborn-v0_8-white')
    fig, ax = plt.subplots(figsize=(12, 7))

    # AAAS/ggsci Palette replication
    # Colors: Navy, Red, Green, Orange, Purple, Light Blue, Yellow, Magenta
    aaas_colors = ["#3b4992", "#ee0000", "#008b45", "#631879", "#008280", "#bb0021", "#5f559b", "#a20056"]

    # Plotting
    pivot_df.plot(
        kind='bar',
        stacked=True,
        ax=ax,
        width=0.9,
        color=aaas_colors,
        edgecolor='white',
        linewidth=0.3
    )

    # --- THEME CUSTOMIZATION (R theme() equivalent) ---
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_color('black')
    ax.spines['bottom'].set_color('black')

    # Scientific notation off (options(scipen=999))
    ax.ticklabel_format(style='plain', axis='y')

    # Font sizing to match your R axis.title and axis.text
    plt.xlabel("Divergence (%)", fontsize=15, fontweight='bold')
    plt.ylabel("Sequence length (bp)", fontsize=15, fontweight='bold')
    plt.xticks(fontsize=10, fontweight='bold')
    plt.yticks(fontsize=10, fontweight='bold')

    # Legend Customization
    plt.legend(
        title='TE Family',
        bbox_to_anchor=(1.05, 1),
        loc='upper left',
        fontsize=10,
        title_fontsize=12,
        frameon=False
    )

    plt.xlim(-0.5, max_div + 0.5)
    plt.tight_layout()

    output_path = rm_output_file.replace(".out", "_unified_R_style.pdf")
    plt.savefig(output_path, dpi=300)
    print(f"Unified plot saved to: {output_path}")
    plt.show()


if __name__ == "__main__":
    file_path = "/Users/panstrugamacbook/Documents/TE_Trimmer/Boundary_model_development/Bhordei_RepeatMasker/RM_Bhordei_bgh_dh14_v4.fa.rm.out"
    generate_unified_landscape_r_style(file_path)