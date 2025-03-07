import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Read the data from the timing results file into a pandas DataFrame
data = []
with open("timing_results.txt", "r") as file:
    for line in file.readlines():
        # Parse each line to extract N, Threads, and Elapsed time
        parts = line.strip().split(", ")
        N = parts[0].split(" = ")[1]
        threads = int(parts[1].split(" = ")[1])
        elapsed_time = parts[3].split(" = ")[1]

        # Convert elapsed time to seconds for easier comparison
        minutes, seconds = elapsed_time.split(":")
        elapsed_seconds = int(minutes) * 60 + float(seconds)

        data.append([N, threads, elapsed_seconds])

# Create a DataFrame from the parsed data
df = pd.DataFrame(data, columns=["N", "Threads", "Elapsed Time (s)"])

# Set up the seaborn style
sns.set(style="whitegrid")

# Create the plot
plt.figure(figsize=(10, 6))
sns.lineplot(data=df, x="Threads", y="Elapsed Time (s)", hue="N", marker="o", palette="tab10")

# Set plot labels and title
plt.title("Elapsed Time vs. Threads for Different Ellipse Sizes", fontsize=16)
plt.xlabel("Number of Threads", fontsize=12)
plt.ylabel("Elapsed Time (seconds)", fontsize=12)

# Show the legend and plot
plt.legend(title="Ellipse N", bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()

# Show the plot
plt.show()