import pandas as pd
import matplotlib.pyplot as plt

# Load the data
data = pd.read_csv("bounding_boxes.csv")

# Plot each net's bounding box
for index, row in data.iterrows():
    plt.plot([row['Xmin'], row['Xmax'], row['Xmax'], row['Xmin'], row['Xmin']],
             [row['Ymin'], row['Ymin'], row['Ymax'], row['Ymax'], row['Ymin']], marker='o')
plt.title('Net Bounding Boxes')
plt.show()
