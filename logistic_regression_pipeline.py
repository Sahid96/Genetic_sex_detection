#!/usr/bin/env python3

# logistic_regression_pipeline.py

# Import the necessary libraries
import numpy as np
import pandas as pd
from sklearn.metrics import confusion_matrix
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression

# Import the CSV file
path = r"E:\logistic_regression.csv"
df = pd.read_csv(path)

# Display the first 10 rows of the dataset
print("First 10 rows of the dataset:")
print(df.head(10))

# Check for null values in the dataset
df.info()

# Define the independent (X) and dependent (y) variables
x = df[["Y/X_Ratio"]]
y = df["Sex"]

# Display first 10 values of X and y
print("First 10 rows of X:")
print(x.head(10))
print("First 10 rows of y:")
print(y.head(10))

# Split data into training and testing sets (80% train, 20% test)
x_train, x_test, y_train, y_test = train_test_split(x, y, test_size=0.2, random_state=42)

# Check the train and test data
print("Training data sample:")
print(x_train.head())
print("Testing data sample:")
print(x_test.head())

# Train the logistic regression model
model = LogisticRegression()
model.fit(x_train, y_train)

# Print the test data
print("Test Data (X_test):")
print(x_test)
print("Test Labels (y_test):")
print(y_test)

# Predict probability of a single value (optional)
sample_value = [[0.011104]]
probability = model.predict_proba(sample_value)
print(f"Prediction probability for {sample_value}: {probability}")

# Check model accuracy
accuracy = model.score(x_test, y_test)
print(f"Model Accuracy: {accuracy}")

# Generate predictions on test data
y_pred = model.predict(x_test)

# Calculate and print confusion matrix
cm = confusion_matrix(y_test, y_pred)
print("Confusion Matrix:")
print(cm)

# Save training data to CSV
train_data = pd.DataFrame(x_train)
train_data['target'] = y_train
train_data.to_csv("train_data.csv", index=False)
print("Training data saved as train_data.csv")

# Save testing data to CSV
test_data = pd.DataFrame(x_test)
test_data['target'] = y_test
test_data.to_csv("test_data.csv", index=False)
print("Testing data saved as test_data.csv")
