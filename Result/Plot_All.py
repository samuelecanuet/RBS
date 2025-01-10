from ROOT import *
from root2mpl import *

filename = "RBS_All_Results.root"
file = TFile(filename, "READ")
canvas = file.Get("All")

#fucntion to convert counter 0 to 7 to 4x2 grid
def convert_counter_to_grid(counter):
    return counter//2, counter%2


#do a list of tuple (min, max) for each window  
windows = [
    (800, 1100),
    (2100, 2800),
    (200, 1200),
    (1900, 2700),
]

counter=-1
fig, ax = plt.subplots(4, 2)
#adding space between subplots
fig.tight_layout(pad=2.0)
for i in range(0, canvas.GetListOfPrimitives().GetSize()):
    pad = canvas.GetListOfPrimitives().At(i)
    counter+=1
    for j in range(0, pad.GetListOfPrimitives().GetSize()):
        sub_obj = pad.GetListOfPrimitives().At(j)
        if sub_obj.InheritsFrom("TH1"): 
            if counter%2==0:
                ylabel="Counts/keV"
            else:
                ylabel=None
            if counter == 7 or counter == 6:
                xlabel="Energy [keV]"
            else:
                xlabel=None 
            DisplayTH1D(sub_obj, ax=ax[convert_counter_to_grid(counter)], lw=0.8, xlim=windows[counter//2], titlesize=10, xlabel=xlabel, ylabel=ylabel, labelsize=15, ticksize=10)
        if sub_obj.InheritsFrom("TLegend"):
            DisplayTLegend(sub_obj, ax=ax[convert_counter_to_grid(counter)], loc="upper right")
        if sub_obj.InheritsFrom("TLatex"):
            DisplayTLatex(sub_obj, ax=ax[convert_counter_to_grid(counter)], color="red", fontsize=14, x=0.015, y=0.75)

# plt.show()
plt.savefig("RBS_All_Results.png", dpi=300)