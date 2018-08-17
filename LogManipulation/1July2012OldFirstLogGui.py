#!/usr/bin/python

from Tkinter import *
from MyModule import *
from FirstLogGuiMods import *

# My main Log manipulation GUI

if __name__ == "__main__":
	"""
	I will attempt to create a simple GUI that allows the user to find
	a LAS file read in a log from that file and run a process on that log
	"""
	global logsFile
	logsFile = BunchOfLogs()
	
	root = Tk()
	root.wm_title("Main Log Menue")
	
	# Defining the main widgets

	spaceLabel1 = Label( root, text="")
	
	findFileLabel = Label( root, text="Find and read in a LAS file: ")
	findFileButton = Button( root, text="Find File", command=lambda: findReadLogLasFile(logsFile) )
	
	processLogLobel = Label( root, text="Run a predefined process on a log: ")
	processLogButton = Button( root, text="Process Log")#, command=procLog(LogsFile) )
	
	outputLogLabel = Label( root, text="Output log to ASCII file: ")
	outputLogButton = Button( root, text="Output ASCII")#, command=outputLog(LogsFile) )
	
	displayAvailableLogsLabel = Label( root, text="Display Available logs information: ")
	displayAvailiableLogsButton = Button( root, text="Display Logs", command=lambda: DisplayContentsOfBunchOfLogs(logsFile) )	

	plotAvailableLogLabel = Label( root, text="Plot Available logs: ")
	plotAvailiableLogButton = Button( root, text="Plot Log", command=lambda: PlotALogWithinBunchOfLogs(logsFile) )	
	
	saveSessionLabel = Label ( root, text="Save variables of current session: ")
	saveSessionButton = Button( root, text="Save Session", command=saveSession )

	restoreSessionLabel = Label ( root, text="Restore variables of saved session: ")
	restoreSessionButton = Button( root, text="Restore Session", command=restoreJob )
	
	exitProcessLabel = Label( root, text="Exit this program: " )
	exitProcessButton = Button( root, text="Exit", command=exitNow )
	
	spaceLabel2 = Label( root, text="")	
	
	
	# Defining display geometry for widgets
	
	rowCounter = 0
	
	spaceLabel1.grid( sticky=E+W, row=0, columnspan=2 )
	rowCounter += 1
	
	findFileLabel.grid( sticky=W, column=0, row=rowCounter, padx=20, pady=5 )
	findFileButton.grid( sticky=W+E, column=1, row=rowCounter, padx = 20, pady=5 )
	rowCounter += 1
	
	processLogLobel.grid( sticky=W, column=0, row=rowCounter, padx=20, pady=5 )
	processLogButton.grid( sticky=W+E, column=1, row=rowCounter, padx = 20, pady=5 )
	rowCounter += 1
		
	outputLogLabel.grid( sticky=W, column=0, row=rowCounter, padx = 20, pady=5 )
	outputLogButton.grid( sticky=W+E, column=1, row=rowCounter, padx = 20, pady=5 )
	rowCounter += 1	
	
	displayAvailableLogsLabel.grid( sticky=W, column=0, row=rowCounter, padx = 20, pady=5 )
	displayAvailiableLogsButton.grid( sticky=W+E, column=1, row=rowCounter, padx = 20, pady=5 )
	rowCounter += 1

	plotAvailableLogLabel.grid( sticky=W, column=0, row=rowCounter, padx = 20, pady=5 )
	plotAvailiableLogButton.grid( sticky=W+E, column=1, row=rowCounter, padx = 20, pady=5 )
	rowCounter += 1
	
	saveSessionLabel.grid( sticky=W, column=0, row=rowCounter, padx = 20, pady=5 )
	saveSessionButton.grid( sticky=W+E, column=1, row=rowCounter, padx = 20, pady=5 )
	rowCounter += 1
	
	restoreSessionLabel.grid( sticky=W, column=0, row=rowCounter, padx = 20, pady=5 )
	restoreSessionButton.grid( sticky=W+E, column=1, row=rowCounter, padx = 20, pady=5 )
	rowCounter += 1	
	
	exitProcessLabel.grid( sticky=W, column=0, row=rowCounter, padx = 20, pady=5 )
	exitProcessButton.grid( sticky=W+E, column=1, row=rowCounter, padx = 20, pady=5 )
	rowCounter += 1	
	
	spaceLabel2.grid( sticky=E+W, row=rowCounter, columnspan=2 )
	rowCounter += 1	
	
	mainloop()
