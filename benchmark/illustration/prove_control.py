import matplotlib.pyplot as plt
import pandas as pd


def plot_equation():
    """
    Plot y = 1 + (5/x)(z-1)
    where x ranges from 0.001 to 0.1
    and z ranges from 1 to 5
    """
    import numpy as np
    
    # Create x values (from 0.0001 to 0.1)
    x = np.arange(5, 3000, 5)
    
    # Create different z values
    z_values = [1.5, 2, 2.5,  3]
    
    # Create figure
    plt.figure(figsize=(10, 6))
    
    # Plot for each z value
    for z in z_values:
        y = 1 + (5 / x) * (z - 1)
        plt.plot(x, y, marker='', linewidth=2, label=f'modification effect ratio = {z}')
    
    # Print values at specific x points
    print("\nValues at specific x points:")
    print("=" * 50)
    x_points = [50, 100, 150]
    for x_val in x_points:
        print(f"\nWhen x = {x_val}:")
        for z in z_values:
            y_val = 1 + (5 / x_val) * (z - 1)
            print(f"  modification effect ratio = {z}: y = {y_val:.4f}")
    print("=" * 50)
    
    plt.xlabel('Metagenome length (Mbp)', fontsize=12)
    plt.ylabel('Control IPD deviations', fontsize=12)
    plt.title('', fontsize=14)
    plt.legend(fontsize=10)
    plt.grid(True, alpha=0.3)
    
    # Add more x-axis ticks for denser labeling
    x_ticks = np.arange(0, 3100, 100)
    plt.xticks(x_ticks, rotation=45)
    
    plt.tight_layout()
    
    # Save the plot
    plt.savefig('../../tmp/figures/framework/equation_plot.png', dpi=300)
    plt.savefig('../../tmp/figures/framework/equation_plot.pdf')
    print("Plot saved as equation_plot.png and equation_plot.pdf")
    plt.show()


plot_equation()