import pandas as pd
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_absolute_error, r2_score
import joblib
import matplotlib.pyplot as plt

df = pd.read_csv('/Users/akhilpeddikuppa/FieldLab/GreiderCodeSearch/analysis/telomere_analysis.csv')

df = df[df['FileName'].str.startswith('JH')]

df = df.dropna(subset=['Age'])
feature_cols = [
    'g_strand_mutations_G>A_g3_per_1k',
    'g_strand_mutations_G>C_g1_per_1k',
    'g_strand_mutations_G>C_g3_per_1k',
    'g_strand_mutations_G>T_g3_per_1k',
    'g_strand_mutations_T>C_t1_per_1k',
    'g_strand_mutations_T>G_t1_per_1k',
    'total_mutations_over_total_g_strand_3xrepeats_per_1k',
]
X = df[feature_cols]
y = df['Age'].astype(float)

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)
model = GradientBoostingRegressor()
model.fit(X_train, y_train)
y_pred = model.predict(X_test)
print("MAE:", mean_absolute_error(y_test, y_pred))
print("R^2:", r2_score(y_test, y_pred))

# Plot predicted vs actual ages
plt.figure(figsize=(6,6))
plt.scatter(y_test, y_pred, alpha=0.7)
plt.plot([y_test.min(), y_test.max()], [y_test.min(), y_test.max()], 'r--')
plt.xlabel('Actual Age')
plt.ylabel('Predicted Age')
plt.title('Predicted vs Actual Age')
plt.tight_layout()
plt.savefig('predicted_vs_actual.png')

joblib.dump(model, 'age_predictor.joblib')