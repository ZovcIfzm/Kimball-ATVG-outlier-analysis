import argparse
import pandas as pd
import numpy as np
import sklearn
from sklearn.cluster import AgglomerativeClustering
from sklearn.cluster import DBSCAN
from sklearn.decomposition import PCA
from sklearn.preprocessing import normalize
import matplotlib.pyplot as plt
from scipy.stats import multivariate_normal, norm
from scipy.spatial.distance import mahalanobis
import seaborn as sns

# Constants
columns = ['Geneid', 'DIO1HighReads', 'DIO1LowReads',
        'DIO2HighReads', 'DIO2LowReads',
        'ND1HighReads', 'ND1LowReads',
        'ND2HighReads', 'ND2LowReads']


# python3 clustering_outliers.py --counts_file count_output.txt --show_plots 0

parser = argparse.ArgumentParser()
parser.add_argument('--counts_file', type=str, default='', help='Path to counts_file')
parser.add_argument('--show_plots', type=int, default=0, help='show plots or not')
# parser.add_argument('--output_file', type=str, default="output_cell.csv", help="path for output file")
args = parser.parse_args()

raw_counts_df = pd.read_csv(args.counts_file, header=1, sep='\t')

# Drop unimportant information
raw_counts_df = raw_counts_df.drop(["Chr", "Start", "End", "Strand", "Length"], axis=1)
raw_counts_df.columns = columns

# Drop genes with less than 10 counts for all cells
raw_counts_df = raw_counts_df.loc[(raw_counts_df[columns[1:]] > 10).all(axis=1)]

#log fold change (dataset with colums of label, feature_1, feature_2)
log_change_df = pd.DataFrame({"Gene_ID": raw_counts_df["Geneid"], \
    "LowLFC": np.log(raw_counts_df["DIO1LowReads"]+raw_counts_df["DIO2LowReads"]+1)-np.log(raw_counts_df["ND1LowReads"]+raw_counts_df["ND2LowReads"]+1), \
    "HighLFC": np.log(raw_counts_df["DIO1HighReads"]+raw_counts_df["DIO2HighReads"]+1)-np.log(raw_counts_df["ND1HighReads"]+raw_counts_df["ND2HighReads"]+1), \
        "Counts": np.log(raw_counts_df["DIO1LowReads"]+raw_counts_df["DIO2LowReads"]+raw_counts_df["ND1LowReads"]+raw_counts_df["ND2LowReads"] \
            + raw_counts_df["DIO1HighReads"]+raw_counts_df["DIO2HighReads"] + raw_counts_df["ND1HighReads"]+raw_counts_df["ND2HighReads"])})


if args.show_plots == 1:
    # Show scatter plot of log fold change
    log_change_df.plot.scatter(x="LowLFC", y="HighLFC", alpha=0.5)
    plt.title("HighLFC vs LowLFC of DIO and ND")
    plt.show()

    # show density plot of log fold change
    sns.kdeplot(log_change_df.LowLFC, log_change_df.HighLFC, bw=.15)
    plt.title("HighLFC vs LowLFC density plot")
    plt.show()

# Data preprocessing complete
### Find outliers
# Convert data to an array of x,y values
data = np.concatenate((np.vstack(log_change_df["LowLFC"]), np.vstack(log_change_df["HighLFC"])), axis=1)

# Find mean and covariance
mean = np.matrix(np.mean(data, axis=0))
cov = np.matrix(np.cov(data, rowvar=0))


# Define mahalnobus distance function
def mahalnobus_dist(value, means, covariance):
    diff = np.matrix(value-means)
    return np.sqrt(diff*np.linalg.inv(covariance)*diff.T)

# Find outliers, is outlier if mahalnobus distance is above 3 (extreme value analysis)
outliers = np.asarray([i for i in range(len(data)) if (mahalnobus_dist(data[i], mean, cov) > 3)])
print("Number of outliers found using extreme value analysis:", outliers.shape)

# Find outliers using DBSCAN, epsilon changed until nubmer of outliers meets that of extreme value analysis
clustering = DBSCAN(eps=0.04657036).fit(data)
print("Number of outliers found using DBScan, eps=0.04657036:", np.count_nonzero(clustering.labels_ == -1))

# Compare outliers
DBSCAN_indices = list(np.nonzero(clustering.labels_ == -1)[0])
Extreme_indices = list(outliers)

jaccard_index = len(list(set(DBSCAN_indices) & set(Extreme_indices)))/len(list(set(DBSCAN_indices) | set(Extreme_indices)))
print("Outliers Jaccard Index:", jaccard_index)

if args.show_plots == 1:
    # Create outlier plots
    DBSCAN_colors = []
    Extreme_colors = []
    for i in range(len(log_change_df)):
        indices = set(DBSCAN_indices)
        if i in indices:
            DBSCAN_colors.append(1)
        else:
            DBSCAN_colors.append(0)

    for i in range(len(log_change_df)):
        indices = set(Extreme_indices)
        if i in indices:
            Extreme_colors.append(1)
        else:
            Extreme_colors.append(0)

            
    DBSCAN_colors = np.asarray(DBSCAN_colors)
    Extreme_colors = np.asarray(Extreme_colors)
    print(np.count_nonzero(DBSCAN_colors == 1))

    colored_plot = log_change_df
    colored_plot["Outlier"] = DBSCAN_colors
    colored_plot.plot.scatter(x="LowLFC", y="HighLFC", c="Outlier", colormap="viridis", alpha=0.5)
    plt.title("Outlier plot according to DBSCAN")
    plt.show()

    colored_plot["Outlier"] = Extreme_colors
    colored_plot.plot.scatter(x="LowLFC", y="HighLFC", c="Outlier", colormap="viridis", alpha=0.5)
    plt.title("Outlier plot according to Extreme value analysis")
    plt.show()

# Find which components explain most variance
pca = PCA(n_components=2)
pca.fit(data)
print("Variance explained by each component:", pca.explained_variance_ratio_)

'''
pca_counts = PCA(n_components=8)
pca_counts.fit(raw_counts_df)
print("Variance explained by components:", pca_counts.explained_variance_ratio_)
'''

# Find outliers using extreme value analysis, using only the first component
single_data = np.asarray(log_change_df["LowLFC"]).flatten()
single_mean = np.mean(single_data)
single_std = np.std(single_data)

outliers = np.asarray([val for val in single_data if (val-single_mean)/single_std > 3])
print("Number of outliers using only first component, std > 3:", outliers.shape)

#clustering = AgglomerativeClustering().fit(data)
#log_change_df["Counts"] = clustering.labels_

'''

log_change_df.plot.scatter(x="LowLFC", y="HighLFC", alpha=0.5)
plt.title("HighLFC vs LowLFC of DIO and ND")
plt.show()

'''