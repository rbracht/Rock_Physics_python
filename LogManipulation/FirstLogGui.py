#!/usr/bin/python

from Tkinter import *
from MyModule import *
from FirstLogGuiMods import *
import processLog as pl

# My main Log manipulation GUI

if __name__ == "__main__":
	"""
	I will attempt to create a simple GUI that allows the user to find
	a LAS file read in a log from that file and run a process on that log
	"""
	global logsFile
	
	logsFile = BunchOfLogs()										# This will be the container for all logs in this session
	
	root = Tk()
	root.wm_title("Main Log Menue")
	root.withdraw
	
	# Defining the main widgets

	spaceLabel1 = Label( root, text="")
	
	findFileLabel = Label( root, text="Find and read in an unwrapped LAS file: ", font=("Comic Sans MS",16), fg="Blue"  )
	findFileButton = Button( root, text="Find File", 
							command = lambda: findReadLogLasFile( logsFile ), 
							font=("Comic Sans MS",14), 
							fg="Black"  )
	
	processLogLobel = Label( root, text="Run a predefined process on a log: ", font=("Comic Sans MS",16), fg="Blue" )
	processLogButton = Button( root, text="Process Log", 
							command = lambda: pl.processLog( root, logsFile ), 
							font=("Comic Sans MS",14), 
							fg="Black"  )
	
	outputLogLabel = Label( root, text="Output log to ASCII file: ", font=("Comic Sans MS",16), fg="Blue" )
	outputLogButton = Button( root, 
							text="Output ASCII")#, command=outputLog(LogsFile) )
	
	displayAvailableLogsLabel = Label( root, text="List currently available logs and information: ", font=("Comic Sans MS",16), fg="Blue" )
	displayAvailiableLogsButton = Button( root, 
							text="List Logs", 
							command = lambda: DisplayContentsOfBunchOfLogs( logsFile ), 
							font=("Comic Sans MS",14), 
							fg="Black"  )	

	plotAvailableLogLabel = Label( root, text="Plot Available logs: ", font=("Comic Sans MS",16), fg="Blue" )
	plotAvailiableLogButton = Button( root, text="Plot Log", 
							command = lambda: cp.crossPlotOfLogs( root, logsFile ), 
							font=("Comic Sans MS",14), 
							fg="Black"  )	
	
	saveSessionLabel = Label ( root, text="Save variables of current session: ", font=("Comic Sans MS",16), fg="Blue" )
	saveSessionButton = Button( root, text="Save Session", 
							command=saveSession, 
							font=("Comic Sans MS",14), 
							fg="Black"  )

	restoreSessionLabel = Label ( root, text="Restore variables of saved session: ", font=("Comic Sans MS",16), fg="Blue" )
	restoreSessionButton = Button( root, text="Restore Session", 
							command=restoreJob, 
							font=("Comic Sans MS",14), 
							fg="Black"  )
	
	exitProcessLabel = Label( root, text="Exit this program: ", font=("Comic Sans MS",16), fg="Blue" )
	exitProcessButton = Button( root, 
							text="Exit", 
							command=exitNow, 
							font=("Comic Sans MS",14), 
							fg="Black"  )
	
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
