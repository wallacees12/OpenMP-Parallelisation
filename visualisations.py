import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np 

num_threads = np.array([i for i in range(1,13)])
optimal_speedup = np.array([i for i in range(1,13)])
wall_clock_3000 = np.array([14.00,10.69,8.00,6.34,5.31,4.47,3.96,3.50,3.23,2.93,2.63,2.44])
speedup_3000 = wall_clock_3000[0]/wall_clock_3000

sns.set_theme(style = 'whitegrid')

plt.figure(figsize=(10,6))
sns.lineplot(x=num_threads, y=speedup_3000, marker = 's', label = 'N = 3000', color = 'orange')
sns.lineplot(x=num_threads, y=optimal_speedup, linestyle ='--',color='black', label="Ideal Speedup (Linear)")

# Labels and legend
plt.xlabel("Number of Threads")
plt.ylabel("Speedup (Relative to 1 Thread)")
plt.title("Parallel Speedup for N-body Simulation")
plt.legend()
plt.show()