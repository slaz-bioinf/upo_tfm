# -*- coding: utf-8 -*-
"""
Created on Thu May  1 11:12:54 2025

TFM: Endometrial transcriptomic signatures of female infertility uncovered 
through RNA-seq and covariate-adjusted differential expression analysis

Author: Sofia Lazaro

Covariates selection pipeline:
This pipeline processes metadata and selects relevant covariates using feature 
selection algorithms for RNA-seq analysis

"""
import pandas as pd
import numpy as np
from sklearn.preprocessing import LabelEncoder
from sklearn.preprocessing import OneHotEncoder
from sklearn.preprocessing import OrdinalEncoder
from sklearn.preprocessing import StandardScaler
from sklearn.impute import KNNImputer
from covariates_utils import remove_accents, replace_decimal

######### DATA PREPROCESSING #########
# %% ----Import and pre-processing of metadata----#

# Define metadata file
metadata = "metadata_29-07.csv"

# Read and store metadata
raw_metadata = pd.read_csv(metadata, sep=";", index_col="Name",encoding = "ANSI") #ANSI is the encoding used by Excel2010

# Process accents in values, columns and index
raw_metadata = raw_metadata.applymap(remove_accents)
raw_metadata.columns = raw_metadata.columns.map(remove_accents)
raw_metadata.index = raw_metadata.index.map(remove_accents)

# Replace "," with "." as decimal separators in values
raw_metadata = raw_metadata.map(replace_decimal)

# Check imported data and data types
for col, dtype in raw_metadata.dtypes.items():
    print(f"{col}: {dtype}")
    
# Define factor covariates
factor_cols = ["Menstrual Cycle Phase"]

# Transform 'object' type columns to 'numeric' 
for col in raw_metadata.columns:
    raw_metadata[col] = pd.to_numeric(raw_metadata[col], errors="ignore")
    
# Transform 'object' type columns to 'category'    
for col in raw_metadata.columns:
    if raw_metadata[col].dtype == "object":
        raw_metadata[col] = raw_metadata[col].astype("category")
    
# Force transformation of factor_cols
raw_metadata[factor_cols] = raw_metadata[factor_cols].astype('category')

#Check types after correction         
for col, dtype in raw_metadata.dtypes.items():
    print(f"{col}: {dtype}")
    
# Ordinal variables
# Define dictionary with levels for each covariate
ordinal_covariates = {
    "Opinion vida": ["MUY INSATISFECHA", "INSATISFECHA", "UN POCO DE AMBOS", "ALGO SATISFECHA", "SATISFECHA", "MUY SATISFECHA"],
    "Estado salud actual": ["NO MUY BUENO", "BUENO", "MUY BUENO"],
    "Desarrollo trabajo": ["Mayormente sentado","Mayormente sentado y caminando","Mayormente caminando","Mayormente caminando y trabajo fisico pesado","Trabajo fisico pesado"],
    "Frecuencia alcohol": ["NUNCA", "UNA O MENOS VECES AL MES", "2-4 VECES AL MES", "2-3 VECES A LA SEMANA", "4 O MAS VECES A LA SEMANA"],
    "Nº consumiciones": ["1 O 2", "3 O 4", "5 O 6"],
    "Frecuencia >6 bebidas": ["NUNCA", "MENOS DE UNA VEZ AL MES", "MENSUALMENTE"],
    "Frecuencia sexo ultimo mes": ["NO HE MANTENIDO", "1 VEZ POR MES", "2-3 VECES AL MES", "1-2 VECES POR SEMANA", "1 VEZ POR SEMANA", "3-4 VECES POR SEMANA"],
    "Condicion fisica general": ["MALA", "ACEPTABLE", "BUENA", "MUY BUENA"],
    "Categoria SUENO": ["Muy malo", "Bastante malo", "Bastante bueno", "Muy bueno"],
    "Categoria ESTRES": ["Bajo", "Medio", "Alto"],
    "Categoria ANSIEDAD": ["Baja", "Media"],
    "Cucharadas AO": ["<4", "AL MENOS 4"],
    "Verduras/hortalizas": ["<2", "AL MENOS 2"],
    "Fruta": ["<3", "AL MENOS 3"],
    "Carne": ["<1", "AL MENOS 1"],
    "Mantequillas": ["<1", "AL MENOS 1"],
    "Bebidas carbonatadas y azucaradas": ["<1", "AL MENOS 1"],
    "Vino": ["NADA", "<7"],
    "Legumbres": ["<3", "AL MENOS 3"],
    "Pescado/marisco": ["<3", "AL MENOS 3"],
    "Reposteria comercial": ["<2", "AL MENOS 2"],
    "Frutos secos": ["<3", "AL MENOS 3"],
    "Sofritos": ["<2", "AL MENOS 2"],
    "Omega-3": ["NO", "OCASIONALMENTE", "A DIARIO"],
    "Acido folico": ["NO", "OCASIONALMENTE", "A DIARIO"],
    "Vitaminas/minerales": ["NO", "OCASIONALMENTE", "A DIARIO"],
    "Hierro": ["NO", "OCASIONALMENTE", "A DIARIO"],
    "Comida rapida": ["NUNCA", "POCAS VECES AL MES", "1-3 VECES AL MES", "1 VEZ POR SEMANA"],
}

#  Apply order in the ordinal variables and check the result
for col, categories in ordinal_covariates.items():
    if col in raw_metadata.columns:
        print(raw_metadata[col].cat.categories)
        raw_metadata[col] = pd.Categorical(raw_metadata[col], categories=categories, ordered=True)
        print(raw_metadata[col].cat.categories)

    
# Create a copy of the metadata data.frame
processed_data = raw_metadata.copy()

# %% ---- Missing values ----#

# Remove covariates with more than a defined number of missing values
threshold_na = 0.30 # proportion of NA

proportion_NA_col = []
colsNA_threshold = []

# Evaluate each column for proportion of NAs
for col in processed_data.columns:
    NA = round(processed_data[col].isnull().mean(), 2)
    # proportion_NA_col.append((col, NA))
    if NA >= threshold_na:
        colsNA_threshold.append((col,NA)) # Store columns with >threshold NA
        processed_data.drop(labels = col, axis=1, inplace = True) # Delete those columns

# Print the deleted columns
print("The following covariates have more", f"{threshold_na:.2%} of missing values and have been deleted from the dataset:")
for col, NA in colsNA_threshold:
    print(f"{col}: {NA:.2%} missing values")
    
# %% ----Encoding----#

# Select only covariates (drop group columns)
X = processed_data.drop(columns=["Group2","Group4"])

# Make a copy of the data frame
X_transformed = X.copy()

# Create lists for each type of column to be processed accordingly
ordinal_cols = []
binary_cols = []
multi_cols = []
numeric_cols = []

# Iterate through each column and encode according to its type
for col in X_transformed.columns:
    # If covariate is ordinal (defined in ordinal_covariates dictionary)
    if col in ordinal_covariates.keys() and col in X_transformed.columns:
        ordinal_cols.append(col)
        encoder_ord = OrdinalEncoder(categories=[ordinal_covariates[col]],handle_unknown='use_encoded_value', unknown_value=np.nan)
        X_transformed[[col]] = encoder_ord.fit_transform(X_transformed[[col]])    
    
    # If categorical but not ordinal
    elif X_transformed[col].dtype == "category" and col not in ordinal_cols: 
        levels = X_transformed[col].nunique() # Count unique values
        
        if levels == 2: # If it is binary
            binary_cols.append(col) # include in binary list
            X_transformed[col] = LabelEncoder().fit_transform(X_transformed[col]) # Transform with label encoder
        
        # Multi-class categorical variable   
        else:
            multi_cols.append(col)
            # Replace missing values with the mode
            X_transformed[col] = X_transformed[col].fillna(X_transformed[col].mode()[0])
    # If numeric
    elif X_transformed[col].dtype in ["int64", "float64"]:
        numeric_cols.append(col)      
    


# Create the transformer for categorical
# sparse_out = False produce a dense array (all options in each attribute are included). Otherwise it delete '0'.
encoder = OneHotEncoder(sparse_output=False) 

# Fit and transform the array for multi-class columns (nominal cols)
X_encoded_multi = encoder.fit_transform(X_transformed[multi_cols])

# Store feature names of multi_cols
multi_feature_names = encoder.get_feature_names_out(multi_cols)

# Transform to data.frame
X_encoded_multi_df = pd.DataFrame(X_encoded_multi, columns = multi_feature_names, index=X_transformed.index)

# Create the scaler for num cols
scaler = StandardScaler()
X_num_scaled = scaler.fit_transform(X_transformed[numeric_cols])
numeric_feature_names =  scaler.get_feature_names_out(numeric_cols)

X_num_scaled_df = pd.DataFrame(X_num_scaled, columns=numeric_feature_names, index=X_transformed.index)

# Concat with the others df
X_final = pd.concat([X_transformed[binary_cols + ordinal_cols], X_encoded_multi_df, X_num_scaled_df], axis=1)

# %% ----Assessment of missing values with KNN imputation----#
imputer = KNNImputer(n_neighbors=5)
X_imputed = imputer.fit_transform(X_final)
X_imputed_df = pd.DataFrame(X_imputed, columns=X_final.columns, index=X_final.index)

# Round binary and ordinal values to restore integer representation
X_imputed_df[binary_cols + ordinal_cols] = X_imputed_df[binary_cols + ordinal_cols].round().astype(int)

# Inverse scale numeric columns to get original scale
X_imputed_df[numeric_cols] = scaler.inverse_transform(X_imputed_df[numeric_cols])

# %% ----Array for feature selection methods----#
# Convert to numpy array
X_encoded = X_imputed_df.to_numpy(dtype=float)

# Store feature names (binary + ordinal + one-hot + numeric)
bin_ord_cols = binary_cols + ordinal_cols
encoded_feature_names = bin_ord_cols + list(multi_feature_names) + list(numeric_feature_names)

# Print encoded matrix
pd.set_option('display.max_columns', None)
print(X_encoded)    
    


# %%----Feature selection----#

# %% Group 2
# Obtain target column, y
y2 = processed_data["Group2"]

#Encode target column into numerical labels (required by sklearn algorithms)
if y2.dtype == 'category':
    y_encoder=LabelEncoder() #For binary class/target column LabelEncoder() is used 
    y_encoded = y_encoder.fit_transform(y2) 

# Apply feature selection algorithms

###### Univariate feature selection SelectKBest #######
from covariates_utils import univariate_ft_sel
ft_univ2 = univariate_ft_sel(X_encoded, y_encoded, encoded_feature_names, k=69)


###### Perform RCEFV wrapped selection ####### Recursive feature elimination with cross-validation
from covariates_utils import rfecv_ft_sel
ft_rfcev2 = rfecv_ft_sel(X_encoded, y_encoded, encoded_feature_names)


###### Feature importance #######
from covariates_utils import trees_ft_importance
ft_imp_RF2 = trees_ft_importance(X_encoded, y_encoded, encoded_feature_names, algorithm="RandomForestClassifier")


###### LASSO #######
from covariates_utils import lasso_ft_sel
ft_lasso2 = lasso_ft_sel(X_encoded, y_encoded, encoded_feature_names, class_type="binomial")

###### Permutation importance #######
from covariates_utils import permutation_ft_sel
ft_perm2=permutation_ft_sel(X_encoded, y_encoded, encoded_feature_names)


####### Feature selection summary #########
# From the top-k feature lists, select features according to the elbow plot
methods_selection2 = [ft_univ2[0:18],ft_rfcev2,ft_imp_RF2[0:4],ft_lasso2,ft_perm2[0:7]]

# Create a dictionary of covariates and the number of methods that selected them
cov_count2 = dict()
for col in X_transformed.columns:
    cov_count2[col]=0
    for df in methods_selection2:
        if col in df["Covariate"].values:
            cov_count2[col]=cov_count2[col]+1
            
#Print the covariates selected in >= 3 methods
for key,value in cov_count2.items():
    if value >= 3:
        print(f"{key} was selected in {value} out of 5 feature selection methods")    
    

   

# %% Group 4
# Obtain target column, y
y4 = processed_data["Group4"]

# Encode target column into numerical labels (required by sklearn algorithms)
if y4.dtype == 'category':
    y_encoder=LabelEncoder()
    y4_encoded = y_encoder.fit_transform(y4) #For class/target column LabelEncoder() is used 

# Apply feature selection algorithms    

###### Perform univariate feature selection: SelectKBest #######
from covariates_utils import univariate_ft_sel
ft_univ4 = univariate_ft_sel(X_encoded, y4_encoded, encoded_feature_names)


###### Perform RCEFV wrapped selection ####### Recursive feature elimination with cross-validation
from covariates_utils import rfecv_ft_sel
ft_rfcev4 = rfecv_ft_sel(X_encoded, y4_encoded, encoded_feature_names)


###### Feature importance #######
from covariates_utils import trees_ft_importance
ft_imp_RF4 = trees_ft_importance(X_encoded, y4_encoded, encoded_feature_names, algorithm="RandomForestClassifier")


###### LASSO #######
from covariates_utils import lasso_ft_sel
ft_lasso4 = lasso_ft_sel(X_encoded, y4_encoded, encoded_feature_names, class_type="multinomial")


###### Permutation importance #######
from covariates_utils import permutation_ft_sel
ft_perm4=permutation_ft_sel(X_encoded, y4_encoded, encoded_feature_names)


####### Feature selection summary #########
# From the top-k feature lists, select features according to the elbow plot
methods_selection4 = [ft_univ4[0:16],ft_rfcev4,ft_imp_RF4[0:15],ft_lasso4,ft_perm4[0:10]]

# Create a dictionary of covariates and the number of methods that selected them
cov_count4 = dict()
for col in X_transformed.columns:
    cov_count4[col]=0
    for df in methods_selection4:
        if col in df["Covariate"].values:
            cov_count4[col]=cov_count4[col]+1
            
#Print the covariates selected in >= 3 methods
for key,value in cov_count4.items():
    if value >= 3:
        print(f"{key} was selected in {value} out of 5 feature selection methods")     




# %% ----Export final metadata----#

# Targets and covariates for RNA seq limma-voom analysis
# Check index to join the dataframes
print(X_imputed_df.index.equals(y2.index))
print(X_imputed_df.index.equals(y4.index))
print(y2.index.equals(y4.index))

final_metadata = pd.concat([y2, y4, X_imputed_df], axis=1)

# Rename columns
final_metadata.index.name = "Filename"
final_metadata = final_metadata.rename(columns={"Group2":"DiseaseClass"})
final_metadata = final_metadata.rename(columns={"Group4":"Type"})

# Change samples name to be the same as in RNAseq rawdata
final_metadata.index = final_metadata.index.str.replace('G0', 'G', regex=False)
final_metadata["Name"]=final_metadata.index
final_metadata["Name_Group"]=final_metadata["Name"].astype(str) + "_" + final_metadata["Type"].astype(str)

# Export the final table with covariates
final_metadata.to_csv("targets_cov.csv", sep = ",", index=True, encoding = "utf-8")
