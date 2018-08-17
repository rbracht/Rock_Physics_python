# My Pyton GUI Modules
# 

from Tkinter import *
from datetime import datetime
import shelve
import tkMessageBox
import tkFileDialog

def findFileGui():
	"""
	Opens up a dialog box for finding files and returns the path
	to the detected file.
	"""
	branch1 = Tk()
	branch1.withdraw()
	
	pathFile = tkFileDialog.askopenfile( parent = branch1, title = 'Choose a file' )
	
	if pathFile != None:
		return pathFile
	else:
		myStdErrorMessage("File read error", "No file was detected")
		
	branch1.destroy

def myStdLable( myTxtFrame, myTxt, mySide=LEFT ):
    """Create a Label within user specified textFrame"""
    myLabel = Label(myTxtFrame)
    myLabel["text"] = myTxt
    myLabel.pack(side=mySide)
    return myLabel

def myStdEntryBox( myTxtFrame, myWidth=50, mySide=LEFT):
    """Create an text entry box within user specified textFrame"""
    myEntryBox = Entry(myTxtFrame)
    myEntryBox["width"] = myWidth
    myEntryBox.pack(side=mySide)

def myStdLableAndBox( myTxtFrame, myTxt="", myWidth=50, mySideTxt=LEFT, mySideBox=LEFT):
    """Create an entry box with text description with user specified textFrame"""
    myText = myStdLable( myTxtFrame, myTxt, mySideTxt )
    myBox = myStdEntryBox( myTxtFrame, myWidth, mySideBox)
    return myBox

def myStdErrorMessage( myError, myMessage ):
    """Create a message box and display the error messasge"""
    tkMessageBox.showerror( str(myError), str(myMessage) ) 

def exitNow():
	"""
	
	Function to exit a program using sys.exit
	
	"""
	
	import sys
	
	branch = Tk()
	branch.withdraw()
	branch.wm_title("Exit Dialogue")
	
	branch1 = Toplevel(branch)
	TerminationLabel = Label( branch1, text="Are you sure you want to exit.", 
						font=("Comic Sans MS",16), fg="Blue" )
	TerminationLabel.grid( columnspan=2, sticky=E+W, padx=20, pady=10 )
	
	R1 = Button( branch1, text="Yes", command=sys.exit, font=("Comic Sans MS",14), fg = "red")
	R1.grid( row=1, column=0, sticky=E+W, padx=20 )
	
	R2 = Button( branch1, text="No", command=branch1.destroy, font=("Comic Sans MS",14), fg = "green" )
	R2.grid( row=1, column=1, sticky=E+W, padx=20 )
	
	L1 = Label(branch1, text="        ").grid( row=2, columnspan=2, padx=20, pady=10 )
	
	branch1.mainloop


def saveSession():
	"""
	
	The function attempts to save all the variables within a particular
	run so it can be restored at a future date. A directory path can be 
	entered below for consistancy of where the files will be stored.
	
	"""

	jobSaveDir = "/Users/browntsunami/.loggerjob/"								# defining session saving directory path
	
	print "+--------------------------------------------------+"
	print "|                                                  |"
	print "|        Entering session saving function.         |"
	print "|                                                  |"
	print "+--------------------------------------------------+"
	
	def shelveSession(fn):
		print "Attempting to save to: "+str(fn)
		my_shelf = shelve.open(fn,'n') # 'n' for new

		for key in dir():
		    try:
		        my_shelf[key] = globals()[key]
		    except TypeError:
		        #
		        # __builtins__, my_shelf, and imported modules can not be shelved.
		        #
		        print('ERROR shelving: {0}'.format(key))
		my_shelf.close()
		print "Session (" + str(fn) + ") saved."
		branch1.destroy
	

	[date,time] = str(datetime.now()).split()
	date = date.replace("-","_")
	time = time.replace(":","_")
	
	fileName = str(jobSaveDir) + "Session_" + str(date) + "_" + str(time)
	print "Default session path/name: "+str(fileName)

	branch = Tk()
	branch.withdraw()
	branch.wm_title("Save Session Dialogue")
	
	branch1 = Toplevel(branch)
	
	rowCounter = 0
	SaveSessionLabel = Label( branch1, text="Enter your session name: ", 
						font=("Comic Sans MS",16), fg="Blue" )
	SaveSessionLabel.grid( row = rowCounter, column = 0, sticky=E+W, padx=20, pady=10 )
	
	SessionNameEntry = Entry( branch1)							# Entry box for a user determined session name 
	SessionNameEntry.delete(0, END)
	SessionNameEntry.insert(0, fileName)						# Set default session filename
	SessionNameEntry.grid ( row=rowCounter, column=1, sticky=E+W+N+S, pady=10, padx=20 )
	rowCounter += 1
	
	if SessionNameEntry.get() != fileName:
		fileName = str(jobSaveDir) + str(SessionNameEntry.get() )
	
	saveButton = Button( branch1, text="Save", command=lambda: shelveSession(fileName) )
	saveButton.grid( row=rowCounter, columnspan=2, padx=20, pady=20 )
	
	ExitButton = Button( branch1, text="  Exit Process	  ", command=branch1.destroy, font=("Comic Sans MS",14), fg="Red" )
	ExitButton.grid(row=rowCounter, column=1, padx=20  )
	rowCounter += 1
	
	branch1.mainloop
	
def restoreJob():
	import os
	"""
	
	This function attempts to restore variables from a previously 
	saved session. This assumes you have the same directory path defined
	below as the one used to save the session.
	
	"""
	
	print "+-------------------------------------------------------+"
	print "|                                                       |"
	print "|        Entering session restoration function.         |"
	print "|                                                       |"
	print "+-------------------------------------------------------+"
	
	
	jobSaveDir = "/Users/browntsunami/.loggerjob/"								# defining session saving directory path
	
	def restoreSession( myfile ):
		my_shelf = shelve.open( myfile )
		for key in my_shelf:
		    globals()[key]=my_shelf[key]
		my_shelf.close()
		print "session (" + str(myfile) + ") restored"
		branch1.destroy
		
	
	rowCounter = 0
	
	branch = Tk()
	branch.withdraw()
	branch.wm_title("Restore Session Dialogue")
	
	branch1 = Toplevel(branch)
	
	
		
	if len(os.listdir( str(jobSaveDir) ) ) == 0:
		noFilesLabel = Label( branch1, "There appears to be no saved sessions in your directory" )
	else:
		chooseFileLabel = Label( branch1, "Choose on of the following session files to restore.")
		chooseFileLabel.grid( row=rowCounter, column=0, font=("Comic Sans MS",16), fg="Blue" )
		rowCounter += 1
		
		i=0
		for fn in os.listdir( str(jobSaveDir) ):
			r = radioButton( branch1, text=str(fn), font=("Comic Sans MS",14), command=lambda: restoreSession( str(jobSaveDir) + str(fn)) )	
			r.grid( row=rowCounter, column=0, padx=20, pady=5 )
	
	branch1.mainloop()	