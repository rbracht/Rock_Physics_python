#! /usr/bin/python

from Tkinter import *
import tkMessageBox
import tkFileDialog
import sys
from MyGuiMod import *


def displayText():
    """ Display the Entry text value. """
    branch1 = Tk()

    listbox = Listbox(branch1)
    listbox.pack(fill=BOTH, expand=1)
	
    global FileEntryWidget, skipLineWidget

    if skipLineWidget.get().strip() == "":
        tkMessageBox.showerror("Second Entry Widget", "Enter number of lines to skip")
    else:
        tkMessageBox.showinfo("/path/filename Entry Widget", skipLineWidget.get().strip())

    if FileEntryWidget.get().strip() == "":
        tkMessageBox.showerror("First Entry Widget", "Enter a /path/filename")
    else:
        f = open(FileEntryWidget.get().strip(), "r")
        count = 0
        for line in f:
	        count = count + 1
	        if count > int(skipLineWidget.get().strip()):
		        listbox.insert(END,str(line))
		
def findFileTest():
    """
    Trying out the GUI file finder in TK returns path and filename
    """
    branch2 = Tk()
    branch2.withdraw()
    branch3 = Tk()
    listbox2 = Listbox(branch3)

    pathFile = tkFileDialog.askopenfile(parent=branch2,title='Choose a file')
    if pathFile != None:
        myStdErrorMessage("File amd path", str(pathFile))
        for line in pathFile:
            print line
            listbox2.insert(END,str(line))
        return pathFile
    else:
        myStdErrorMessage("File read error", "No file was detected")

def ExitNow():
    sys.exit()
		
if __name__ == "__main__":
    root = Tk()
 
    root.title("/path/filename Entry Widget")
    root["padx"] = 20
    root["pady"] = 20 

    # Create a text frame to hold the text Label and the Entry widget

    textFrame1 = Frame(root)
    FileEntryWidget = myStdLableAndBox(textFrame1,"Enter /path/filename to be read:",100)
    textFrame1.pack()

    textFrame2 = Frame(root)
    skipLineWidget = myStdLableAndBox(textFrame2,"Enter number of lines to skip:    ",5)
    textFrame2.pack(side=LEFT,anchor=NW)
    textFrame2.pack()

    textFrame3 = Frame(root)
    myStdLable(textFrame3,"----------------------------------------------------------------------")
    textFrame3.pack(anchor=S)

    textFrame4 = Frame(root)
    button = Button(textFrame4, text="Submit", command=displayText)
    button.pack(anchor=S)
    textFrame4.pack()

    textFrame5 = Frame(root)
    button = Button(textFrame5, text="Find File", command=findFileTest)
    button.pack(anchor=S)
    textFrame5.pack()

    textFrame6 = Frame(root)
    button = Button(textFrame6, text="Exit", command=ExitNow)
    button.pack(anchor=S)
    textFrame6.pack()

    root.mainloop()