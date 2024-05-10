#' Wyatt Petryshen
#' Yale University
#' May 10th, 2024
#'
#' Python program to interactively define a habitat matrix

from tkinter import *
import numpy as np
from tkinter import ttk

# Gui functions

# Confiugre Window
root = Tk()
root.title("Create Landscape")
root.configure(background = "white")
root.minsize(500,500)
root.maxsize(1500,1000)
root.geometry("1250x750")

# Internal functions
# Create numpy grid
def create_grid(nrow,ncol):
    print("Generate")

# Assign 1 value to grid cell if selected
def assign_grid_values(nrow,ncol):
    print("Habitatable")

# Save function
def save_grid(nrow,ncol):
    print("Landscape Saved")

# Create Input Frame
input_frame = Frame(root, width = 200, height = 200, bg = "lightgrey")
input_frame.grid(row=0, column=0, padx=10, pady=5,sticky=E)

# Assign user inputs
input_section = Label(input_frame, text = "Initial Landscape Parameters:", bg = "white")
input_section.grid(row = 0, column = 0, columnspan=4, padx=5, pady=5)

# ncol, nrow, save path
Label(input_frame, text="Number of Rows:", bg="white").grid(row=2, column=0, padx=5, pady=5)
nrow = Entry(input_frame, bd=3)
nrow.grid(row=2, column=1, padx=5, pady=5)

Label(input_frame, text="Number of Cols:", bg="white").grid(row=3, column=0, padx=5, pady=5)
ncol = Entry(input_frame, bd=3)
ncol.grid(row=3, column=1, padx=5, pady=5)

Label(input_frame, text="Save Path:", bg="white").grid(row=4, column=0, padx=5, pady=5)
savepath = Entry(input_frame, bd=3)
savepath.grid(row=4, column=1, padx=5, pady=5)

# Add generate grid button
generate_land = Button(input_frame, text = "Generate Landscape")
generate_land.grid(column=1,row = 5)

# Add save button
save = Button(input_frame, text = "Save Landscape")
save.grid(column=1, row = 6)

# Create Landscape Grid Frame
input_frame = Frame(root, width = 900, height = 750, bg = "teal")
input_frame.grid(row=0, column=10, padx=10, pady=5)


# Initalize Window
root.mainloop()
