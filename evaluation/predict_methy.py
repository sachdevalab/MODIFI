from sklearn.datasets import make_regression
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.ensemble import RandomForestRegressor
from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import SVR
from sklearn.model_selection import train_test_split
import pandas as pd
import numpy as np
import os
from imblearn.under_sampling import RandomUnderSampler

def handle_features(loaded_feature):
    field = loaded_feature[0]
    ## field[0] is like '103123230211', transform it to a data array with each element has a int in field[0]
    base_feature = []
    for i in range(len(field[0])):
        base_feature.append(int(field[0][i]))
    ipd_fearure = loaded_feature[1:]
    ## transform the ipd_feature to a list of float
    ipd_feature = [float(i) for i in ipd_fearure]
    return base_feature + ipd_feature

def load_data(data_file):
    ## load the csv data without header, and sep  by tab
    df = pd.read_csv(data_file, header=None, sep='\t')
    ### randaomly subsample the data
    df = df.sample(frac=0.01)

    raw_X = df[1].tolist()
    ## each feature is separated by space
    X = []
    for i in raw_X:
        loaded_feature = i.split(' ')
        X.append(handle_features(loaded_feature))
    y = df[0].tolist()
    print (len(X), len(y))
    return X, y

def model(X, y, seed):
    # Assuming X and y are defined and properly shaped
    X = np.array(X).reshape(-1, 1)
    y = np.array(y)
    # balance the data by subsampling

    rus = RandomUnderSampler(random_state=42)
    X, y = rus.fit_resample(X, y )

    X_train, X_test, y_train, y_test = train_test_split(
        X, y, random_state=seed, test_size=0.1)

    print (len(X_train), len(X_test), len(y_train), len(y_test))
    clf = RandomForestClassifier(
                n_estimators=100,
                max_depth=12,
                random_state=seed,
        )
    clf.fit(X_train, y_train)
    print (clf.score(X_test, y_test))
    ## claculate the prediction accuracy and AUC
    y_pred = clf.predict(X_test)
    from sklearn.metrics import accuracy_score
    print (accuracy_score(y_test, y_pred))




data_file = 'borg/train_data.csv'
X, y = load_data(data_file)
print (X[:2], y[:2])
model(X, y, 1)