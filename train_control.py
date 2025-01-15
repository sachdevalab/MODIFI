from sklearn.datasets import make_regression
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.ensemble import RandomForestRegressor
from sklearn.svm import SVR
from sklearn.model_selection import train_test_split
import pandas as pd
import numpy as np
import os


def model(X, y, seed):
    # Assuming X and y are defined and properly shaped
    # X = np.array(X).reshape(-1, 1)
    X = np.array(X)
    y = np.array(y)
    print (len(X), len(y))
    for i in range(10):
        print (X[i], y[i])

    X_train, X_test, y_train, y_test = train_test_split(
        X, y, random_state=seed, test_size=0.1)

    reg = HistGradientBoostingRegressor(
        learning_rate=0.35,
        max_iter=100000,
        max_depth=12,
        random_state=seed,
    )

    # reg = RandomForestRegressor(
    #         n_estimators=100,
    #         max_depth=12,
    #         random_state=seed,
    # )

    # reg = SVR(    ### very slow
    #         kernel='rbf',
    #         C=1.0,
    #         epsilon=0.1
    #     )

    reg.fit(X_train, y_train)
    # reg.predict(X_test)
    result = reg.score(X_test, y_test)
    print (result)
    print (reg.predict(X_test[1:5]))
    print (y_test[1:5])
    for i in range(10):
        print (X_test[i], reg.predict(X_test[i].reshape(1, -1)), y_test[i])

def load_data(input_file, downsample_rate=0.1):
    df = pd.read_csv(input_file)
    ## shuffle the data
    df = df.sample(frac=1)
    ## downsample the data
    df = df.sample(frac=downsample_rate)
    ## the file has 8 columns, the final column is the target, and other columns are features
    X = df.iloc[:, :-2]
    y = df.iloc[:, -1]
    S = df.iloc[:, -2]
    ## extract the last two column into a new df
    df2 = df.iloc[:, -2:]
    
    return X, y, S, df2

def random_split(df, train_frac, validation_frac):
    # Shuffle the entire DataFrame
    df = df.sample(frac=1, random_state=123).reset_index(drop=True)

    # Calculate split indices
    train_end = int(len(df) * train_frac)
    validation_end = train_end + int(len(df) * validation_frac)

    # Split the DataFrame
    train_df = df[:train_end]
    validation_df = df[train_end:validation_end]
    test_df = df[validation_end:]

    return train_df, validation_df, test_df

def get_bet(df, bert_dir):
    ## downsample the data
    df = df.sample(frac=0.001)

    ## if the folder does not exist, create it
    if not os.path.exists(bert_dir):
        os.makedirs(bert_dir)
    train_df, validation_df, test_df = random_split(df, 0.7, 0.1)
    # Test size is implied to be 0.2 as the remainder

    train_df.to_csv(bert_dir+ "train.csv", index=None, sep=',')
    validation_df.to_csv(bert_dir + "dev.csv", index=None, sep=',')
    test_df.to_csv(bert_dir + "test.csv", index=None,   sep=',')



if __name__ == "__main__":
    seed = 1
    folder = "../data/"

    input_file = folder + "/control_train_data.csv"
    bert_dir = folder + "/bert/"

    
    X, y, S , df2= load_data(input_file, 1)
    print (X.shape, y.shape)
    get_bet(df2, bert_dir)
    
    # model(X, y, seed)