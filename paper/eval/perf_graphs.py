rust_sasa_points = [8.071]
freesasa_points = [54.914]
biopython_points = [368.025]

import matplotlib.pyplot as plt

# Data for plotting
libraries = ["rust-sasa", "freesasa", "biopython"]
times = [8.071, 54.914, 368.025]
errors = [0.361, 0.455, 51.156]

# Create the plot
plt.figure(figsize=(10, 6))
plt.rcParams.update({"font.size": 16})
plt.bar(libraries, times, yerr=errors, color=["#33BBEE", "#EE7733", "#EE3377"], capsize=5)
plt.ylabel("Time (seconds)")
plt.title("Performance Comparison of SASA Libraries on E. coli proteome")
# plt.yscale("log")  # Using log scale due to large differences in values

# Add value labels on bars
for i, v in enumerate(times):
    plt.text(i, v + errors[i], f"{v}s", ha="center", va="bottom")

# Add buffer space above bars to prevent text from hitting figure top
max_height = max(times[i] + errors[i] for i in range(len(times)))
plt.ylim(0, max_height * 1.10)

plt.tight_layout()
plt.savefig("performance_comparison.pdf")
