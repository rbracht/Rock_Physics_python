from Tkinter import *

def sel():
	test = Tk()
	selection = "You selected the option " + str( var.get() )
	label.config(text = selection)
	button = Button( test, text="This is a test exit", command=test.destroy).pack()
	test.mainloop
	
root = Tk()

var = IntVar()

selection = "Nothing selected yet"

# def sel():
#    selection = "You selected the option " + str( var.get() )
#    label.config(text = selection)

for i in range( 4 ):
	R1 = Radiobutton(root, text="Option " + str( i+1 ), variable=var, value=i+1, command=sel )
	R1.pack( anchor = W )

# R1 = Radiobutton(root, text="Option 1", variable=var, value=1,
#                   command=sel)
# R1.pack( anchor = CENTER )
# 
# R2 = Radiobutton(root, text="Option 2", variable=var, value=2,
#                   command=sel)
# R2.pack( anchor = W )
# 
# R3 = Radiobutton(root, text="Option 3", variable=var, value=3,
#                   command=sel)
# R3.pack( anchor = W)

label = Label(root)
label.config(text = selection)
label.pack()

root.mainloop()