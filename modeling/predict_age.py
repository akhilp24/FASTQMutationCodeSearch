import joblib
import pandas as pd
import sys
from features import extract_features_from_fastq

if len(sys.argv) < 2:
    print("Usage: python predict_age.py <fastq_file>")
    sys.exit(1)

fastq_file = sys.argv[1]
model = joblib.load('age_predictor.joblib')
features = extract_features_from_fastq(fastq_file)
X_new = pd.DataFrame([features])

# Select only the features the model was trained on
X_new = X_new[model.feature_names_in_]

predicted_age = model.predict(X_new)
print("Predicted Age:", predicted_age[0])
