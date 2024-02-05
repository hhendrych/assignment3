"""Script created with Python 3, started on 30th of January 2024
@author: Jan Hendrych
S-number: S5687217
Institution: University of Groningen
Email: j.hendrych@student.rug.nl

This script aims to compare isotopic profiles (d13C and d15N) gathered from different locations"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import colors, ticker
from scipy.stats import shapiro, f_oneway

# Read data from two CSV files
df1 = pd.read_csv('isotopes.tsv', sep="\t")
df2 = pd.read_csv('locations.tsv', sep="\t")

columns_df1 = ["stable_isotope_assay_is", "marecon_sample_id"]
columns_df2 = ["stable_isotope_assay_is", "d13C", "d15N"]

merged_df = pd.merge(df1[columns_df1], df2[columns_df2], on="stable_isotope_assay_is", how='inner')
print(merged_df)
merged_df.to_csv("merged_df.tsv", sep="\t", index=False)

# Read data from the merged file
df = pd.read_csv('merged_df.tsv', sep="\t")

# Select specific columns
columns_df = ["stable_isotope_assay_is", "marecon_sample_id", "d13C", "d15N"]

# Extract the values for the x and y axes
x = df["d13C"]
y = df["d15N"]

# Extract the first two characters of the second column for color
color_values = df["marecon_sample_id"].str[:2]

# Map string values to colors using a colormap rainbow
colormap = plt.cm.get_cmap('rainbow', len(color_values.unique()))
normalized_values = colors.Normalize(vmin=0, vmax=len(color_values.unique()))
colors_mapped = colormap(normalized_values(range(len(color_values.unique()))))

# Map the color to the corresponding colors
color_mapping = dict(zip(color_values.unique(), colors_mapped))

# Use the mapped colors in scatter plot
colors_for_scatter = color_values.map(color_mapping)

# Scatter plot with original data, colored based on the first two characters
scatter = plt.scatter(x, y, c=colors_for_scatter)

# Legend
legend_labels = []
for key in color_values.unique():
    legend_labels.append(plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=color_mapping[key],
                                    markersize=8, label=key))

plt.legend(handles=legend_labels, loc="lower left", fontsize="x-small")

# Titles and labels
plt.title("Fin Whales - Isotopic values")
plt.xlabel("δ13C (‰)")
plt.ylabel("δ15N (‰)")

# Plot
plt.savefig("scatter_all_locations.png", dpi=600)

plt.show()

# Numeric values for output file
with open("output_means.txt", "w") as output:
    locations = ["GL", "IT", "GM", "PT", "EI", "SP", "NE", "IL"]

    output.write("Location\tSample size\t\tStatistics\t\t\tMean d13C \t\t\t\t\tMean d15N \n")

    for location in locations:
        # Filter rows where "location" is present in the "marecon_sample_id" column
        filtered_df = df[df["marecon_sample_id"].str.contains(location)]

        # Calculate the mean of the "d13C" and "d15N" values in the filtered_df
        mean_d13C = filtered_df["d13C"].mean()
        std_dev_d13C = filtered_df["d13C"].std()
        mean_d15N = filtered_df["d15N"].mean()
        std_dev_d15N = filtered_df["d15N"].std()

        # Results with 9 decimal places
        output.write(f"{location}\t\t\t{len(filtered_df)}\t\t\t\tMean+-SD Range\t"
                     f"\t{mean_d13C:.9f}+-{std_dev_d13C:.9f}\t{mean_d15N:.9f}+-{std_dev_d15N:.9f}\n")

grouped_locations = {"Western Atlantic": ["GM", "GL"], "Central Atlantic": ["IL"],
                     "Eastern Atlantic": ["NE", "PT", "EI"], "Mediterranean Sea": ["IT"], "Spain": ["SP"]}

# Read the merged_df.tsv file into a DataFrame
df_mean = pd.read_csv('merged_df.tsv', sep='\t')

with open("output_means.txt", "a") as output:
    for group_name, location_list in grouped_locations.items():
        # Filter rows where 'marecon_sample_id' contains any location in the current group
        filtered_df_mean = df_mean[df_mean['marecon_sample_id'].str.contains('|'.join(location_list))]

        # Extract the relevant columns for the dataset
        selected_columns = ['d13C', 'd15N']
        filtered_dataset = filtered_df_mean[selected_columns]

        # Calculate the mean for each column
        mean_values = filtered_dataset.mean()

        # Write results
        output.write(f"\nGroup: {group_name}\n")
        output.write("Sample size\tStatistics\tMean d13C \t\t\t\tMean d15N \n")

        # Calculate the mean values for each group
        mean_d13C_group = df_mean[df_mean['marecon_sample_id'].str.contains('|'.join(location_list))]['d13C'].mean()
        mean_d15N_group = df_mean[df_mean['marecon_sample_id'].str.contains('|'.join(location_list))]['d15N'].mean()

        # Write mean values for the group
        output.write(f"{len(filtered_df_mean)}\t\t\tMean\t\t{mean_d13C_group}\t\t{mean_d15N_group}\n")

# Plot with grouped locations
colormap = plt.cm.get_cmap('coolwarm', len(grouped_locations))
normalized_values = colors.Normalize(vmin=0, vmax=len(grouped_locations))
colors_mapped = colormap(normalized_values(range(len(grouped_locations))))

# Map the group values to the corresponding colors in the colormap
group_mapping = dict(zip(grouped_locations.keys(), colors_mapped))

# Scatter plots for each group
for group, locations in grouped_locations.items():
    group_df = df[df["marecon_sample_id"].str[:2].isin(locations)]

    # Extract the values for the x and y axes
    x = group_df["d13C"]
    y = group_df["d15N"]

    # Use the mapped color for the scatter plot
    scatter = plt.scatter(x, y, label=group, c=group_mapping[group])

    # Axis ticks, nine decimal places
    plt.gca().xaxis.set_major_formatter(ticker.StrMethodFormatter("{x:.1f}"))
    plt.gca().yaxis.set_major_formatter(ticker.StrMethodFormatter("{x:.1f}"))

# Plot titles and labels
plt.title("Fin Whales - Isotopic values - grouped")
plt.xlabel("δ13C (‰)")
plt.ylabel("δ15N (‰)")

# Legend
plt.legend(fontsize="small")

# Plot
plt.savefig("scatter_grouped_locations.png", dpi=600)

plt.show()

# Shapiro-Wilk test
def shapiro_test(data, column_name):
    with open("output_means.txt", "a") as output:

        stat, p_value = shapiro(data)
        output.write(f"\n\nShapiro-Wilk Test for {column_name}:\nStatistic: {stat}\nP-value: {p_value}")
        if p_value < 0.05:
            output.write(f"\nThe data does not follow a normal distribution.")
        else:
            output.write(f"\nThe data follows a normal distribution.")

# Combine data from all groups
combined_df = df[df["marecon_sample_id"].str[:2].isin([location for locations in grouped_locations.values()
                                                       for location in locations])]

# Extract the values for the x and y axes
x_combined = combined_df["d13C"]
y_combined = combined_df["d15N"]

# Perform Shapiro-Wilk test for both columns on the combined data
print(shapiro_test(x_combined, "d13C"))
print(shapiro_test(y_combined, "d15N"))

# One-way ANOVA
# Read data from two CSV files
merged_df.to_csv("merged_df.tsv", sep="\t", index=False)

# Read data from the merged file
df = pd.read_csv('merged_df.tsv', sep="\t")

# Drop rows with missing values
df = df.dropna(subset=['d13C', 'd15N'])

# Define groups based on marecon_sample_id
grouped_locations = {
    'Western Atlantic': ['GM', 'GL'],
    'Central Atlantic': ['IL'],
    'Eastern Atlantic': ['NE', 'PT', 'EI'],
    'Mediterranean Sea': ['IT', 'SP']
}

df['group'] = df['marecon_sample_id'].str[:2]

# One-way ANOVA
f_statistic, p_value = f_oneway(df['d13C'][df['group'].isin(grouped_locations['Western Atlantic'])],
                                df['d13C'][df['group'].isin(grouped_locations['Central Atlantic'])],
                                df['d13C'][df['group'].isin(grouped_locations['Eastern Atlantic'])],
                                df['d13C'][df['group'].isin(grouped_locations['Mediterranean Sea'])],
                                df['d15N'][df['group'].isin(grouped_locations['Western Atlantic'])],
                                df['d15N'][df['group'].isin(grouped_locations['Central Atlantic'])],
                                df['d15N'][df['group'].isin(grouped_locations['Eastern Atlantic'])],
                                df['d15N'][df['group'].isin(grouped_locations['Mediterranean Sea'])])

# Print results
print(f"One-way ANOVA F-statistic: {f_statistic}")
print(f"P-value: {p_value}")

# Compare p-value
if p_value < 0.05:
    print("There is a significant difference in means between at least two groups.")
else:
    print("There is no significant difference in means between groups.")

# Randomization
observed_diff_within_between = 0.0

for group, locations in grouped_locations.items():
    # Calculate within-group variability (e.g., variance)
    within_group_variability = np.var(df[['d13C', 'd15N']][[x in locations for x in df['marecon_sample_id'].str[:2]]],
                                      axis=0).sum()

    # Calculate between-group variability (e.g., variance)
    between_group_variability = np.var(df[['d13C', 'd15N']][[x not in locations
                                                             for x in df['marecon_sample_id'].str[:2]]], axis=0).sum()

    # Accumulate the observed difference
    observed_diff_within_between += within_group_variability - between_group_variability

# Number of permutations for the randomization test
num_permutations = 1000

# Initialize an array to store permuted differences
permuted_diffs = np.zeros(num_permutations)

# Perform the randomization test
for i in range(num_permutations):
    # Randomly permute the group labels
    permuted_group_labels = np.random.permutation(df['marecon_sample_id'].str[:2])

    # Calculate within-group and between-group variability for permuted labels
    permuted_diff = 0.0
    for group, locations in grouped_locations.items():
        # Calculate within-group variability for permuted labels
        permuted_within_group_variability = np.var(df[['d13C', 'd15N']]
                                                   [[x in locations for x in permuted_group_labels]], axis=0).sum()

        # Calculate between-group variability for permuted labels
        permuted_between_group_variability = np.var(df[['d13C', 'd15N']]
                                                    [[x not in locations for x in permuted_group_labels]], axis=0).sum()

        # Accumulate the permuted difference
        permuted_diff += permuted_within_group_variability - permuted_between_group_variability

    # Store the permuted difference
    permuted_diffs[i] = permuted_diff

# Calculate the p-value
p_value_randomization = np.mean(permuted_diffs >= observed_diff_within_between)

# Print results
print(f"Observed Difference (Within-Group vs. Between-Group Variability): {observed_diff_within_between}")
print(f"P-value (Randomization Test): {p_value_randomization}")

# Compare p-value
if p_value_randomization < 0.05:
    print("\n\nThe isotopic variation within groups is significantly different from between groups.")
else:
    print("There is no significant difference in isotopic variation between and within groups.")

print("\n\nAll done!")
