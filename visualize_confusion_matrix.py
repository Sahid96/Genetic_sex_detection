#!/usr/bin/env python3

# visualize_confusion_matrix.py

import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.metrics import confusion_matrix
import numpy as np

# Example usage assumes model, x_test, and y_test are already defined
def plot_confusion_matrix(model, x_test, y_test, output_file="confusion_matrix.svg"):
    # Generate predictions
    y_pred = model.predict(x_test)

    # Compute confusion matrix
    cm = confusion_matrix(y_test, y_pred)

    # Define labels for the matrix
    group_names = ["True Negative (TN)", "False Positive (FP)",
                   "False Negative (FN)", "True Positive (TP)"]
    group_counts = [f"{value}" for value in cm.flatten()]
    group_percentages = [f"{100 * value / sum(cm.flatten()):.1f}%" for value in cm.flatten()]

    labels = [f"{name}\n{count}\n{percent}" for name, count, percent in zip(group_names, group_counts, group_percentages)]
    labels = np.array(labels).reshape(2, 2)

    # Create a figure and axis for proper saving
    fig, ax = plt.subplots(figsize=(6, 4), dpi=300)

    sns.heatmap(cm, annot=labels, fmt='', cmap='Blues', xticklabels=['Pred 0', 'Pred 1'], 
                yticklabels=['Actual 0', 'Actual 1'], cbar=False, linewidths=1, linecolor='black', ax=ax)

    ax.set_xlabel("Predicted Label")
    ax.set_ylabel("True Label")
    ax.set_title("Confusion Matrix with TP, TN, FP, FN")

    # Save the figure properly as SVG
    fig.savefig(output_file, format="svg", bbox_inches="tight", transparent=True)
    plt.show()

    print(f"Confusion matrix saved as '{output_file}'")


def main():
    print("This script requires a trained model and test data.")
    print("Please define and pass `model`, `x_test`, and `y_test` to the `plot_confusion_matrix` function.")


if __name__ == "__main__":
    main()
