import numpy as np
import polars as pl
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_squared_error
from sklearn.preprocessing import StandardScaler
import joblib
import matplotlib.pyplot as plt
import constants

file_path = f"{constants.DATA_FOLDER}/chembl_gpcr.parquet"
activities = pl.read_parquet(file_path)
activities = activities.drop_nans()

features = activities["fingerprint"].to_list()
label = activities.select(["standard_value"]).to_numpy().ravel()
X_train, X_test, y_train, y_test = train_test_split(features, label, test_size=0.2, random_state=42)

#print(X_train)
print(y_train)

# # Scale features
# scaler = StandardScaler()
# X_train = scaler.fit_transform(X_train)
# X_test = scaler.transform(X_test)

model = RandomForestRegressor(n_estimators=100, random_state=42)
model.fit(X_train, y_train)
joblib.dump(model, f"{constants.DATA_FOLDER}/standard_value_predictor.pkl")
print("training done")


with open(f"{constants.DATA_FOLDER}/standard_value_predictor.pkl", "rb") as file:
    model = joblib.load(file)

y_pred = model.predict(X_test)

plt.figure(figsize=(8, 5))
plt.scatter(y_test, y_pred, alpha=0.5)
plt.plot([min(y_test), max(y_test)], [min(y_test), max(y_test)], color='red', linestyle='dashed')  # y = x line
plt.xlabel("Actual Values")
plt.ylabel("Predicted Values")
plt.title("Actual vs. Predicted Values")
plt.savefig(f"{constants.DATA_FOLDER}/RF_bioact_fp_regres.png")


rmse = np.sqrt(mean_squared_error(y_test, y_pred))
print(f"RMSE: {rmse}")
print(y_test)
print(y_pred)

# Save model