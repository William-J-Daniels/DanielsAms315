import pandas as pd
import matplotlib.pyplot as plt
import os


listing = os.listdir("data/Project1_PartB_data/")
listing.sort()

num_plots = len(listing)
fig, axs = plt.subplots(num_plots, layout="tight",
                        figsize=(6.4, 4.8*num_plots*0.75))

for file, ax in zip(listing, axs):
    data = pd.read_csv("data/Project1_PartB_data/"+file)
    
    ax.set_title(file)
    ax.set_xlabel("x")
    ax.set_ylabel("y")

    ax.scatter(data["x"], data["y"],
               marker=".", color="black")
    
fig.savefig("figures/all_part_bs.png")