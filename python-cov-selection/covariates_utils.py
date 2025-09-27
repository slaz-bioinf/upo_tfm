# -*- coding: utf-8 -*-
"""
Created on Wed Jun 25 20:20:20 2025

TFM: Endometrial transcriptomic signatures of female infertility uncovered 
through RNA-seq and covariate-adjusted differential expression analysis

Author: Sofia Lazaro


Functions for processing metadata and performing feature selection

"""


__version__ = "1.0"

######## Metadata processing functions ########

import unicodedata

def remove_accents(text):
    """Remove accents from a string"""
    # Check if the input is of type str
    if isinstance(text, str): 
        # Split letters and accents as individual characters
        text_normalized = unicodedata.normalize('NFD', text) 
        #Generate the empty str
        result = '' 
        # Iterate through the str (character by character)
        for char in text_normalized: 
            # Evaluate if the character is NOT an accent
            if unicodedata.category(char) != 'Mn': 
                # If NOT, add it to the str result
                result += char
        return result
    else:
        # If the input is not a str, it returns the original value "text"
        return text 
    
def replace_decimal(text):
    """Replace "," by "." for numeric values that are str type"""
    # Check if the input is of type str
    if isinstance(text, str):
        # Replace decimals symbol
        correct_num = text.replace(",", ".")
    else:
        correct_num = text
    return correct_num


######## sklearn feature selection functions ########

import matplotlib.pyplot as plt
import numpy as np
from sklearn.feature_selection import SelectKBest, mutual_info_classif
import pandas as pd


def univariate_ft_sel(X_encoded, y_encoded, encoded_feature_names, k=30, method=mutual_info_classif):
    """Perform univariate feature selection using SelectKBest.
    X_encoded: array with encoded values of covariates
    y_encoded: array with encoded values of response variable
    encoded feature names: list of feature names
    k: top k features that SelectKBest will return (default = 30)
    method: scoring function (default = mutual_info_classif), options: f_classif, chi2, mutual_info_classif
    
    The function returns a data.frame with selected covariates and prints and
    elbow plot of the selected features and their respective scores"""
    
    # Create the feature selector
    selector = SelectKBest(score_func=lambda X, y: method(X, y, random_state=42), k=k)
    
    # Fit and transform the selector with the input data
    X_new = selector.fit_transform(X_encoded, y_encoded)

    # Get the selected feature indices
    selected_indices = selector.get_support(indices=True)
    
    # Get the feature scores
    scores = selector.scores_
    # Store selected features and scores
    selected_features = []
    selected_scores = []
    for i in selected_indices:
        selected_features.append(encoded_feature_names[i])
        selected_scores.append(scores[i])


    # Create the data frame with the selected features and sccores
    ft_univ = pd.DataFrame({"Feature": selected_features, "Score": selected_scores})
    
    # Sort features by score (descending)
    ft_univ_sorted=ft_univ.sort_values(ascending=False, by="Score").reset_index()

    # Split encoded feature names into covariate and value (for one-hot encoded variables)
    ft_univ_sorted[['Covariate', 'Value']] = ft_univ_sorted['Feature'].str.rsplit('_', n=1, expand=True)
    
    # Display top 10 selected features
    pd.set_option('display.max_columns', None)
    print("Top features selected by SelectKBest:")
    print(ft_univ_sorted.head(10))

    # Plot feature scores in descending order (elbow method visualization)
    scores = ft_univ_sorted["Score"].values
    plt.plot(range(len(scores)), sorted(scores, reverse=True), marker="o")
    plt.xlabel("Covariates number")
    plt.ylabel("Score")
    plt.title("Elbow method for top covariates selection: SelectKBest")
    plt.show()

    # Return final data.frame
    return(ft_univ_sorted)


from sklearn.feature_selection import RFECV
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold

def rfecv_ft_sel(X_encoded, y_encoded, encoded_feature_names, min_feat=1, cv=10):
    
    """ Perform feature selection using Recursive Feature Elimination with Cross-Validation (RFECV).
    X_encoded: array with encoded values of covariates
    y_encoded: array with encoded values of response variable
    encoded feature names: list of feature names
    min_feat: Minimum number of features to select (default=1)
    cv:  Number of cross-validation folds (StratifiedKFold),(default=10)
    The function returns a data.frame with selected covariates and prints a
    plot with the number of features and the accuracy of the model"""
    
    # Convert feature names to NumPy array for indexing
    feat_names = np.array(encoded_feature_names)
    
    # Minimum number of features to include
    min_features_to_select = min_feat  
    
    # Classification algorithm (Logistic Regression).
    # Automatically adapts to binary or multinomial depending on the target.
    clf = LogisticRegression(random_state=42) 
    
    # StratifiedKFold maintains the same class distribution in each fold
    cv = StratifiedKFold(cv)
    
    # Define RFECV object
    rfecv = RFECV(
        estimator=clf,
        step=1,
        cv=cv,
        scoring="accuracy",
        min_features_to_select=min_features_to_select
    )

    # Fit model and run recursive feature elimination
    rfecv.fit(X_encoded, y_encoded)

    print(f"Optimal number of features: {rfecv.n_features_}")


    # Extract relevant cross-validation results
    data = {
        key: value
        for key, value in rfecv.cv_results_.items()
        if key in ["n_features", "mean_test_score", "std_test_score"]
    }
    cv_results = pd.DataFrame(data)
    
    # Plot accuracy vs number of features
    plt.figure()
    plt.xlabel("Number of features selected")
    plt.ylabel("Mean test accuracy")
    plt.errorbar(
        x=cv_results["n_features"],
        y=cv_results["mean_test_score"],
        yerr=cv_results["std_test_score"],
    )
    plt.title("Recursive Feature Elimination \nwith correlated features")
    plt.show()

    # Get names of selected features
    selected_features = feat_names[rfecv.support_]
    print("Selected features after RFECV:", selected_features)


    # Create data.frame with covariates and importance scores
    ft_rfcev = pd.DataFrame({"Covariate": selected_features, "RFECV importance": rfecv.estimator_.coef_[0]})
    
    # Sort features by importance score (descending)
    ft_rfcev = ft_rfcev.sort_values(by="RFECV importance", ascending=False).reset_index(drop=True)


    return(ft_rfcev)

def trees_ft_importance(X_encoded, y_encoded, encoded_feature_names, algorithm, n_estimators = 10000):
    
    """ Create tree-based models and obtain feature importance 
    X_encoded: array with encoded values of covariates
    y_encoded: array with encoded values of response variable
    encoded feature names: list of feature names
    algorithm: Algorithm to use: "ExtraTreesClassifier" or "RandomForestClassifier"
    n_estimators: Number of trees in the ensemble (default=10000)
    The function returns a data.frame with selected covariates and prints and
    elbow plot of the selected features and their respective scores"""
    
    # Select model depending on the algorithm argument
    if algorithm == "ExtraTreesClassifier":
        from sklearn.ensemble import ExtraTreesClassifier
        model = ExtraTreesClassifier(n_estimators=n_estimators, random_state=42)
    elif algorithm == "RandomForestClassifier":
        from sklearn.ensemble import RandomForestClassifier
        model = RandomForestClassifier(n_estimators=n_estimators, random_state=42)

    else:
        print("Only ExtraTreesClassifier or RandomForestClassifier are supported")
        return None
    
    # Fit model and print feature importances
    model.fit(X_encoded, y_encoded)
    print(model.feature_importances_)
    
    # Create DataFrame with features and their importance scores
    ft_imp = pd.DataFrame({"Feature": encoded_feature_names, "Ft importance": model.feature_importances_})

    # Sort by importance (descending)
    ft_imp = ft_imp.sort_values(by="Ft importance", ascending=False).reset_index(drop=True)
    
    # Split encoded feature names into covariate and value (for one-hot encoded variables)
    ft_imp[['Covariate', 'Value']] = ft_imp['Feature'].str.rsplit('_', n=1, expand=True)

    # Display top 10 features
    pd.set_option('display.max_columns', None)
    print(f"Top features selected from {algorithm} model:")
    print(ft_imp.head(10))

    # Plot elbow represtantion with importances and number of features
    importances = ft_imp["Ft importance"].values
    plt.plot(range(len(importances)), sorted(importances, reverse=True), marker="o")
    plt.xlabel("Number of covariates")
    plt.ylabel("Feature importance")
    plt.title(f"Elbow method for top covariates selection: {algorithm}")
    plt.plot(range(len(importances)), sorted(importances, reverse=True), marker="o")
    plt.show()
    
    return(ft_imp)
    

from sklearn.model_selection import GridSearchCV    

def lasso_ft_sel(X_encoded, y_encoded, encoded_feature_names, class_type = "binomial"):
    """ Perform feature selection using LASSO regularization (Logistic Regression with L1 penalty).
    Hyperparameter C is optimized with cross-validated grid search. 
    X_encoded: array with encoded values of covariates
    y_encoded: array with encoded values of response variable
    encoded feature names: list of feature names
    class_typee: Type of classification of response variable: "binomial" or "multinomial
    
    The function returns a data.frame with selected covariates."""
    
    # Set classification mode for Logistic Regression
    if class_type == "binomial":
        multi_class = "ovr"
    elif class_type == "multinomial":
        multi_class = "multinomial"
    
    # Define search range for regularization parameter C 
    param_grid = {'C': np.logspace(-3, 3, 15)}  # from 0.001 to 1000

    # Create the Logistic Regression model with LASSO penalty
    model = LogisticRegression(
        penalty='l1', #LASSO
        solver='saga',
        multi_class = multi_class,  # 'ovr' para binaria, 'multinomial' para multinomial
        max_iter=10000,
        random_state=42
    )

    # Grid search with cross-validation to find optimal C
    grid_search = GridSearchCV(
        model,
        param_grid,
        cv=5,
        scoring='accuracy',
        n_jobs=-1,
        verbose=1
    )

    # Fit model
    grid_search.fit(X_encoded, y_encoded)

    # Extract the best model after grid search
    best_model = grid_search.best_estimator_

    # Coefficients vector (1D)
    coef_vector = best_model.coef_[0]  # vector 1D con coeficientes

    # Create DataFrame with feature names and coefficients
    ft_LASSO = pd.DataFrame({
        'Covariate': encoded_feature_names,
        'Coefficients': coef_vector
    })

    # Select features with non-zero coefficients (selected by LASSO)
    ft_LASSO_sel = ft_LASSO[ft_LASSO['Coefficients'] != 0].sort_values(by='Coefficients', ascending=False)

    # Display selected features
    print("Selected features by LASSO: ")
    print(ft_LASSO_sel.reset_index(drop=True))


    return(ft_LASSO_sel)

from sklearn.inspection import permutation_importance

def permutation_ft_sel(X_encoded, y_encoded, encoded_feature_names):
    """ Perform feature selection using Permutation Importance with Logistic Regression.
    X_encoded: array with encoded values of covariates
    y_encoded: array with encoded values of response variable
    encoded feature names: list of feature names 
    The function returns a data.frame with covariates and their mean permutation importance,
        sorted by importance in descending order."""
    
    # Base model: Logistic Regression with L2 regularization
    model = LogisticRegression(penalty='l2', solver='liblinear', max_iter=1000, random_state=42)
    
    # Fit model to data
    model.fit(X_encoded, y_encoded)
    
    # Estimate feature importance via permutations
    # Shuffle each feature column and see how much accuracy drops
    result = permutation_importance(model, X_encoded, y_encoded, n_repeats=10, random_state=42)
    
    # Create data.frame with features and their mean importance scores
    ft_perm = pd.DataFrame({
    "Covariate": encoded_feature_names,
    "Importance": result.importances_mean,
    })
    
    # Sort by importance (descending)
    ft_perm = ft_perm.sort_values(by="Importance", ascending=False).reset_index(drop=True)
    
    # Print top 10 features
    print("Top covariates with Permutation Importance (LogisticRegression):")
    print(ft_perm.head(10))


    # Elbow plot using importance value vs. number of features
    importances = ft_perm["Importance"].values
    plt.plot(range(len(importances)), sorted(importances, reverse=True), marker="o")
    plt.xlabel("Number of covariates")
    plt.ylabel("Importance")
    plt.title("Elbow method for top covariates selection: Permutation importance")
    plt.plot(range(len(importances)), sorted(importances, reverse=True), marker="o")
    plt.show()    

    return ft_perm