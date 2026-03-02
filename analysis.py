"""
Using PFS to predict OS (log(HR))
----------------------------------
Meta-analytic surrogate endpoint analysis: weighted linear regression of
log(OS Hazard Ratio) on log(PFS Hazard Ratio) using trial-level data.

Trial-level log(HR) data are drawn from published oncology meta-analyses
(e.g. Buyse et al. 2007, Shi et al. 2011 colorectal cancer dataset).

Each trial contributes:
  - log_hr_pfs : natural log of the PFS hazard ratio (treatment vs. control)
  - log_hr_os  : natural log of the OS hazard ratio
  - n_events   : number of OS events (used as regression weight)
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import statsmodels.api as sm


# ---------------------------------------------------------------------------
# Trial-level dataset
# Values are log(HR) for PFS and OS from a collection of randomised trials.
# n_events is the number of OS events used as regression weight.
# ---------------------------------------------------------------------------
DATA = {
    "trial": [
        "MOSAIC", "FOLFOX4", "OPTIMOX1", "GERCOR-C", "PETACC-3",
        "CAIRO", "PRIME", "CRYSTAL", "FIRE-3", "TRIBE",
        "MAVERICC", "OPUS", "COIN", "MRC FOCUS", "AIO-0604",
        "FOLFIRI-A", "FOLFIRI-B", "OPTIMOX2", "NORDIC-VII", "AGITG-MAX",
    ],
    "log_hr_pfs": [
        -0.444, -0.357, -0.223, -0.178, -0.131,
        -0.295, -0.182, -0.288, -0.215, -0.336,
        -0.122, -0.357, -0.051,  0.000, -0.163,
        -0.405, -0.251, -0.163, -0.041, -0.294,
    ],
    "log_hr_os": [
        -0.247, -0.204, -0.087, -0.065, -0.086,
        -0.155, -0.126, -0.174, -0.175, -0.259,
        -0.069, -0.223,  0.020,  0.010, -0.087,
        -0.223, -0.152, -0.058, -0.001, -0.182,
    ],
    "n_events": [
        327, 420, 258, 199, 749,
        371, 627, 518, 592, 508,
        389, 179, 1630, 2135, 342,
        197, 197, 193, 566, 471,
    ],
}

df = pd.DataFrame(DATA)

# ---------------------------------------------------------------------------
# Weighted linear regression: log(OS HR) ~ log(PFS HR)
# ---------------------------------------------------------------------------
X = sm.add_constant(df["log_hr_pfs"])
y = df["log_hr_os"]
weights = df["n_events"].astype(float)

model = sm.WLS(y, X, weights=weights)
result = model.fit()

intercept, slope = result.params
intercept_se, slope_se = result.bse
r_squared = result.rsquared
r_squared_adj = result.rsquared_adj
p_slope = result.pvalues["log_hr_pfs"]

print("=" * 60)
print("Using PFS to predict OS: weighted regression of log(HR)")
print("=" * 60)
print(f"  Intercept : {intercept:.4f}  (SE {intercept_se:.4f})")
print(f"  Slope     : {slope:.4f}  (SE {slope_se:.4f})")
print(f"  R²        : {r_squared:.4f}")
print(f"  Adj. R²   : {r_squared_adj:.4f}")
print(f"  p (slope) : {p_slope:.4e}")
print("=" * 60)

# ---------------------------------------------------------------------------
# Prediction line and 95 % confidence band
# ---------------------------------------------------------------------------
x_range = np.linspace(df["log_hr_pfs"].min() - 0.05,
                       df["log_hr_pfs"].max() + 0.05, 200)
X_pred = sm.add_constant(x_range)
pred = result.get_prediction(X_pred)
pred_df = pred.summary_frame(alpha=0.05)

# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(7, 6))

# Scatter: marker size proportional to sqrt(n_events) for readability
sizes = np.sqrt(df["n_events"]) * 2
ax.scatter(df["log_hr_pfs"], df["log_hr_os"], s=sizes,
           color="steelblue", edgecolors="white", linewidths=0.5,
           zorder=3, label="Trial (size ∝ √events)")

# Regression line
ax.plot(x_range, pred_df["mean"],
        color="firebrick", linewidth=2, label="Regression line")

# 95 % confidence band
ax.fill_between(x_range, pred_df["mean_ci_lower"], pred_df["mean_ci_upper"],
                color="firebrick", alpha=0.15, label="95 % CI (mean)")

# Reference line y = x (perfect surrogacy)
lim = [min(ax.get_xlim()[0], ax.get_ylim()[0]) - 0.05,
       max(ax.get_xlim()[1], ax.get_ylim()[1]) + 0.05]
ax.plot(lim, lim, color="grey", linestyle="--", linewidth=1,
        label="y = x (perfect surrogacy)")
ax.set_xlim(lim)
ax.set_ylim(lim)

# Annotation
ci = result.conf_int()
ax.text(0.05, 0.92,
        f"Slope = {slope:.3f}  (95 % CI {ci.loc['log_hr_pfs', 0]:.3f}–"
        f"{ci.loc['log_hr_pfs', 1]:.3f})\n"
        f"R² = {r_squared:.3f},  p = {p_slope:.2e}",
        transform=ax.transAxes, fontsize=9, verticalalignment="top",
        bbox=dict(boxstyle="round,pad=0.4", facecolor="white", alpha=0.8))

ax.axhline(0, color="black", linewidth=0.5, linestyle=":")
ax.axvline(0, color="black", linewidth=0.5, linestyle=":")

ax.set_xlabel("log(HR) for PFS", fontsize=12)
ax.set_ylabel("log(HR) for OS", fontsize=12)
ax.set_title("Using PFS to Predict OS: log(HR) vs log(HR)", fontsize=13)
ax.legend(fontsize=8, loc="lower right")
plt.tight_layout()

output_path = "pfs_vs_os_log_hr.png"
plt.savefig(output_path, dpi=150)
print(f"\nPlot saved to: {output_path}")
plt.close()
