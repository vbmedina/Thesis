from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

# Paths and thresholds
input   = Path("./pp.csv")
out     = Path("./pp_with_labels.csv")
hist    = Path("./activity_pIC50_hist.png")
main_treshold        = 6.0   # primary “active” threshold  (pIC50 ≥ 6.0)
upperclass_threshold = 7.5   # secondary “upper-class” threshold (pIC50 ≥ 7.5)

def main() -> None:
    # Load data 
    if not input.exists():
        raise FileNotFoundError(f"Input file not found: {input.resolve()}")
    df = pd.read_csv(input)

    if "pIC50" not in df.columns:
        raise KeyError("Column 'pIC50' is missing from the input CSV.")

    # Plot histogram
    fig, ax = plt.subplots(figsize=(6, 4))
    df["pIC50"].hist(bins=60, ax=ax, edgecolor="black")
    ax.axvline(main_treshold,        ls="--", lw=1.5, label=f"main ≥ {main_treshold}")
    ax.axvline(upperclass_threshold, ls=":",  lw=1.5, label=f"upper ≥ {upperclass_threshold}")
    ax.set_xlabel("pIC50")
    ax.set_ylabel("Count")
    ax.set_title("pIC50 distribution with activity thresholds")
    ax.legend()
    fig.tight_layout()
    fig.savefig(hist, dpi=300)
    plt.close(fig)

    # Boolean label columns 
    df["active_6"] = df["pIC50"] >= main_treshold       
    df["elite_75"] = df["pIC50"] >= upperclass_threshold 

    # Save 
    df.to_csv(out, index=False)

    # Print summary 
    print(f"Wrote labelled data - {out.resolve()}")
    print(f"Wrote histogram - {hist.resolve()}\n")

    print("Label balance summary:")
    print(f"  active_6  True : {df['active_6'].sum():7d} "
          f"({df['active_6'].mean():.1%} of rows)")
    print(f"  elite_75  True : {df['elite_75'].sum():7d} "
          f"({df['elite_75'].mean():.1%} of rows)")

if __name__ == "__main__":
    main()
