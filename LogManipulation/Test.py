from Tkinter import *
from time import *

def blink(rectangle, canvas):
    for i in range(3):
        canvas.itemconfigure(rectangle, fill = "red")
        sleep(1)
        canvas.itemconfigure(rectangle, fill = "white")
        sleep(1)


root = Tk()
fr = Frame(root)
fr.pack()
canv = Canvas(fr, height = 100, width = 100)
canv.pack()
rect = canv.create_rectangle(25, 25, 75, 75, fill = "white")
blink(rect, canv)
canv.show()
root.mainloop()
