import pandas as pd
import matplotlib.pyplot as plt

# Load results
df = pd.read_csv("output/results.csv")

# Keep only rows where objective improved compared to all previous runs
best_so_far = -float("inf")
improving_rows = []

for idx, row in df.iterrows():
    if row["objective"] > best_so_far:
        improving_rows.append({
            "run": idx + 1,
            "objective": row["objective"]
        })
        best_so_far = row["objective"]

improving_df = pd.DataFrame(improving_rows)

# Plot
plt.figure(figsize=(10, 6))
plt.plot(
    improving_df["run"],
    improving_df["objective"],
    marker="o"
)

# Annotate points
for _, row in improving_df.iterrows():
    plt.text(
        row["run"],
        row["objective"],
        f"{row['objective']:.3f}",
        fontsize=9
    )

plt.xlabel("Run Number")
plt.ylabel("Best Objective Score")
plt.title("Objective Score Improvement Over Runs")
plt.grid(True)

plt.tight_layout()
plt.savefig("output/objective_improvement.png", dpi=300)
plt.show()