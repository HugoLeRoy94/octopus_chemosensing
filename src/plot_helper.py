import matplotlib.pyplot as plt
import numpy as np

def plot_grid_heatmap(data, cmap='jet', ax=None,figsize=(12,10)):
    """Core function to plot a matrix as a grid heatmap."""
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    
    # Render the matrix
    im = ax.imshow(data, aspect='auto', cmap=cmap, interpolation='none')
    
    # Draw the grid lines at the edges of the cells
    rows, cols = data.shape
    for x in range(cols + 1):
        ax.axvline(x - 0.5, color='k', linewidth=0.5)
    for y in range(rows + 1):
        ax.axhline(y - 0.5, color='k', linewidth=0.5)
    plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    return im, ax

def apply_fancy_labels(im, ax, x_labels, y_labels, title, y_title):
    """Applies specific ticks, labels, and colorbars to an existing heatmap."""
    # Set Ticks
    ax.set_xticks(np.arange(len(x_labels)))
    ax.set_xticklabels(x_labels, rotation=90, fontsize=8)
    
    ax.set_yticks(np.arange(len(y_labels)))
    ax.set_yticklabels(y_labels, fontsize=10)
    
    # Decorations
    
    ax.set_title(title, fontsize=14)
    ax.set_xlabel('Sorted by mean ligand response', fontsize=12)
    ax.set_ylabel(y_title, fontsize=12)
    plt.tight_layout()