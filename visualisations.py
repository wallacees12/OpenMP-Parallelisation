import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

sns.set_theme(style="whitegrid")

# --- Data for Normal Test (steps=200) for N=3000 ---
# Thread counts (from 1 to 14)
threads = np.arange(1, 15)

# Openmp Normal Test times in seconds
openmp_normal = np.array([10.25, 5.24, 3.53, 2.69, 2.16, 1.81, 1.58, 1.43, 1.39, 1.34, 1.33, 1.39, 1.35, 1.36])
# Pthreads Normal Test times in seconds
pthreads_normal = np.array([7.10, 5.34, 4.02, 3.16, 2.62, 2.22, 1.94, 1.75, 1.58, 1.46, 1.40, 1.36, 1.41, 1.33])
# Vanilla Normal Test (only one value, so we assume it’s independent of thread count)
vanilla_normal = 9.73

# # --- Plot 1: Wall Clock Time for Normal Test (N=3000) ---
# plt.figure(figsize=(8,6))
# plt.plot(threads, openmp_normal, marker='o', label='Openmp (Normal)')
# plt.plot(threads, pthreads_normal, marker='s', label='Pthreads (Normal)')
# # Plot Vanilla as a horizontal dashed line
# plt.hlines(vanilla_normal, threads[0], threads[-1], colors='red', linestyles='dashed', label='Vanilla (Normal)')
# plt.xlabel("Number of Threads")
# plt.ylabel("Wall Clock Time (s)")
# plt.title("Wall Clock Time for N=3000 (Normal Test, steps=200)")
# plt.legend()
# plt.show()

# # --- Compute Speedups for Normal Test ---
# # Speedup is defined relative to the 1-thread time.
# openmp_speedup = openmp_normal[0] / openmp_normal
# pthreads_speedup = pthreads_normal[0] / pthreads_normal

# # --- Plot 2: Speedup for Normal Test (N=3000) ---
# plt.figure(figsize=(8,6))
# plt.plot(threads, openmp_speedup, marker='o', label='Openmp Speedup')
# plt.plot(threads, pthreads_speedup, marker='s', label='Pthreads Speedup')
# # Plot the ideal (linear) speedup: ideal speedup equals the number of threads.
# plt.plot(threads, threads, 'k--', label='Ideal Speedup')
# plt.xlabel("Number of Threads")
# plt.ylabel("Speedup (Relative to 1 Thread)")
# plt.title("Speedup for N=3000 (Normal Test, steps=200)")
# plt.legend()
# plt.show()

# # --- Data for Final Test (steps=400) for N=3000 ---
# openmp_final = np.array([41.21, 20.95, 14.16, 10.70, 8.64, 7.26, 6.33, 5.69, 5.54, 5.57, 5.43, 5.39, 5.44, 5.42])
# pthreads_final = np.array([28.32, 21.44, 15.94, 12.64, 10.48, 8.90, 7.77, 6.94, 6.32, 5.85, 5.75, 5.31, 5.27, 5.10])
# vanilla_final = 39.08

# # --- Plot 3: Wall Clock Time for Final Test (N=3000, steps=400) ---
# plt.figure(figsize=(8,6))
# plt.plot(threads, openmp_final, marker='o', label='Openmp (Final)')
# plt.plot(threads, pthreads_final, marker='s', label='Pthreads (Final)')
# plt.hlines(vanilla_final, threads[0], threads[-1], colors='red', linestyles='dashed', label='Vanilla (Final)')
# plt.xlabel("Number of Threads")
# plt.ylabel("Wall Clock Time (s)")
# plt.title("Wall Clock Time for N=3000 (Final Test, steps=400)")
# plt.legend()
# plt.show()

# openmp_speedup = openmp_final[0] / openmp_final
# pthreads_speedup = pthreads_final[0] / pthreads_final

# # --- Plot 4: Speedup for Normal Test (N=3000) steps 400 ---
# plt.figure(figsize=(8,6))
# plt.plot(threads, openmp_speedup, marker='o', label='Openmp Speedup')
# plt.plot(threads, pthreads_speedup, marker='s', label='Pthreads Speedup')
# # Plot the ideal (linear) speedup: ideal speedup equals the number of threads.
# plt.plot(threads, threads, 'k--', label='Ideal Speedup')
# plt.xlabel("Number of Threads")
# plt.ylabel("Speedup (Relative to 1 Thread)")
# plt.title("Speedup for N=3000 (Normal Test, steps=400)")
# plt.legend()
# plt.show()

# # ================================
# # N=1000 Normal Test (steps=200)
# # ================================
# openmp_normal_1000 = np.array([2.30, 1.18, 0.81, 0.62, 0.51, 0.43, 0.38, 0.36, 0.36, 0.36, 0.43, 0.38, 0.37, 0.37])
# pthreads_normal_1000 = np.array([1.57, 1.20, 0.90, 0.71, 0.60, 0.51, 0.45, 0.40, 0.36, 0.34, 0.31, 0.29, 0.28, 0.27])
# vanilla_normal_1000 = np.full(len(threads), 2.17)

# openmp_speedup_1000 = openmp_normal_1000[0] / openmp_normal_1000
# pthreads_speedup_1000 = pthreads_normal_1000[0] / pthreads_normal_1000
# vanilla_speedup_1000 = np.full(len(threads), 1.0)

# # Plot for N=1000 Normal Test
# fig, axs = plt.subplots(1, 2, figsize=(14, 6))

# # Wall Clock Time plot
# axs[0].plot(threads, openmp_normal_1000, marker='o', label='Openmp')
# axs[0].plot(threads, pthreads_normal_1000, marker='s', label='Pthreads')
# axs[0].hlines(vanilla_normal_1000[0], threads[0], threads[-1], colors='r', linestyles='dashed', label='Vanilla')
# axs[0].set_xlabel("Number of Threads")
# axs[0].set_ylabel("Time (s)")
# axs[0].set_title("N=1000 Normal Test (steps=200): Wall Clock Time")
# axs[0].legend()

# # Speedup plot
# axs[1].plot(threads, openmp_speedup_1000, marker='o', label='Openmp')
# axs[1].plot(threads, pthreads_speedup_1000, marker='s', label='Pthreads')
# axs[1].plot(threads, threads, 'k--', label='Ideal')
# axs[1].hlines(1, threads[0], threads[-1], colors='r', linestyles='dashed', label='Vanilla')
# axs[1].set_xlabel("Number of Threads")
# axs[1].set_ylabel("Speedup")
# axs[1].set_title("N=1000 Normal Test (steps=200): Speedup")
# axs[1].legend()

# plt.tight_layout()
# plt.show()


y = np.array(["Vanilla", "Pthreads", "OpenMP"])
x = np.array([39.08, 5.10, 5.42])

# Set Seaborn style
sns.set(style="whitegrid")

# Create the bar plot
plt.figure(figsize=(10, 6))
ax = sns.barplot(x=x, y=y, palette="viridis", edgecolor="black")

# Add title and labels
plt.title("Wall Clock Timing for N = 3000, Steps =", fontsize=16, pad=20)
plt.xlabel("Time (seconds)", fontsize=14)
plt.ylabel("Implementation", fontsize=14)

# Annotate the bars with their values
for i, value in enumerate(x):
    ax.text(value, i, f"{value:.2f}s", va="center", ha="left", fontsize=12, color="black")

# Customize the grid
ax.grid(True, linestyle="--", alpha=0.7)

# Remove spines
sns.despine(left=True, bottom=True)

# Show the plot
plt.tight_layout()
plt.show()
