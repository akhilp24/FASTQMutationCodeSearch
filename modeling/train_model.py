import pandas as pd
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_absolute_error, r2_score
import joblib
import matplotlib.pyplot as plt
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import GridSearchCV

csv_path = '/Users/akhilpeddikuppa/FieldLab/GreiderCodeSearch/analysis/telomere_analysis.csv'

df = pd.read_csv(csv_path)

print("Reading CSV at path: ", csv_path)

# Filter to only rows where Age is NOT empty (not NaN)
print("Filtering out rows where Age is empty")
df = df[df['Age'].notna()]



# Use all per_1k features
feature_cols = [
    'g_strand_mutations_T>C_t1_per_1k',
   
    'Telomere_Length',
]


print("Feature selection", feature_cols)


X = df[feature_cols]
y = df['Age'].astype(int)

# Hyperparameter grid
param_grid = {
    'n_estimators': [100, 200, 500],
    'max_depth': [None, 10, 20, 30],
    'min_samples_split': [2, 5, 10]
}

rf = RandomForestRegressor(random_state=42)
grid_search = GridSearchCV(rf, param_grid, cv=5, scoring='neg_mean_absolute_error')
grid_search.fit(X, y)

print("Best params:", grid_search.best_params_)
print("Best MAE:", -grid_search.best_score_)

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.8, random_state=42)
print("Training model...")
model = RandomForestRegressor(random_state=42)

# 5-fold cross-validation, using negative mean absolute error as the scoring metric
cv_scores = cross_val_score(model, X, y, cv=5, scoring='neg_mean_absolute_error')

# Convert negative MAE to positive
mae_scores = -cv_scores

print("Cross-validated MAE scores:", mae_scores)
print("Mean CV MAE:", mae_scores.mean())
print("Std CV MAE:", mae_scores.std())

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

importances = grid_search.best_estimator_.feature_importances_
feat_names = X.columns
plt.figure(figsize=(8,6))
plt.barh(feat_names, importances)
plt.xlabel('Feature Importance')
plt.title('Random Forest Feature Importances')
plt.tight_layout()
plt.savefig('feature_importances.png')