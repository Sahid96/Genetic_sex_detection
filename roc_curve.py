#!/usr/bin/env python3

# stylish_roc_curve.py

import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, auc
from scipy.ndimage import gaussian_filter1d
import joblib  # for loading the saved model, if needed

def plot_stylish_roc(model, x_test, y_test, output_file="Roc_curve.svg"):
    # Get the predicted probabilities for the positive class
    y_prob = model.predict_proba(x_test)[:, 1]

    # Compute ROC curve
    fpr, tpr, thresholds = roc_curve(y_test, y_prob)

    # Compute AUC (Area Under Curve)
    roc_auc = auc(fpr, tpr)

    # Apply Gaussian filter for smooth curve
    tpr_smooth = gaussian_filter1d(tpr, sigma=2)

    # Create figure with dark background
    fig, ax = plt.subplots(figsize=(7, 5), dpi=300, facecolor="black")
    ax.set_facecolor("black")

    # Define colors
    curve_color = "#00CFFF"      # Neon Blue for ROC Curve
    fill_color = "#00CFFF40"     # Transparent Blue for Shaded Area
    diagonal_color = "#FF6F61"   # Soft Red for the diagonal line
    text_color = "white"         # Ensure text is visible on dark background

    # Plot the smoothed ROC curve
    ax.plot(fpr, tpr_smooth, color=curve_color, lw=3, label=f'ROC Curve (AUC = {roc_auc:.3f})')

    # Fill under the curve for a modern effect
    ax.fill_between(fpr, 0, tpr_smooth, color=fill_color, alpha=0.6)

    # Plot the random chance diagonal line
    ax.plot([0, 1], [0, 1], linestyle="dashed", lw=2, color=diagonal_color, label="Chance Level (AUC = 0.5)")

    # Customize labels and titles
    ax.set_xlabel("False Positive Rate (FPR)", fontsize=12, color=text_color)
    ax.set_ylabel("True Positive Rate (TPR)", fontsize=12, color=text_color)
    ax.set_title("Stylish Receiver Operating Characteristic (ROC) Curve", fontsize=14, color=text_color)

    # Customize ticks
    ax.tick_params(axis='both', colors=text_color)

    # Add the AUC value as text inside the plot
    ax.text(0.6, 0.2, f"AUC = {roc_auc:.3f}", fontsize=14, color=text_color,
            bbox=dict(facecolor="black", edgecolor="white", boxstyle="round,pad=0.3"))

    # Add a legend with a visible white border
    legend = ax.legend(loc="lower right", fontsize=10, facecolor="black", edgecolor="white", framealpha=0.9)
    for text in legend.get_texts():
        text.set_color(text_color)

    # Save the stylish ROC curve as an SVG file
    fig.savefig(output_file, format="svg", bbox_inches="tight", transparent=True)

    # Show the plot
    plt.show()
    print(f"ROC curve saved as '{output_file}'")


def main():
    # Replace these with actual data/model loading as needed
    # Example:
    # model = joblib.load("your_model.pkl")
    # x_test = np.load("x_test.npy")
    # y_test = np.load("y_test.npy")

    print("This script requires a trained model and test data.")
    print("Please load your model and data inside the main() function.")


if __name__ == "__main__":
    main()
