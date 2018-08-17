# Temporary storage for cross-plotting functions

import Tkinter as tk
import MyModule as mm
import SimPy.SimPlot as sp

def crossplotLogs( aBunchOfLogs = mm.BunchOfLogs() ):
	
	branch = tk.Tk()
	branch.wm_title( "Cross plotting logs" )
	branch.withdraw
	
	class simPlotXY:
		def __init__( self, X = [], Y = []):
			self.Z = zip( self.X, self.Y )
			self.plt = SimPlot()
			plt.plotLine( self.Z )
			plt.mainloop()
	
	class set_X_VariablesIn:
		def __init__( self, aBunchOfLogs = mm.BunchOfLogs, logNum = -999 ):
			mainCrossplotFrame.xIdVariable.set( str( aBunchOfLogs.MyLogLabel[LogNum] ) )
			mainCrossplotFrame.xTypeVariable.set( str( aBunchOfLogs.LogType[LogNum] ) )
			mainCrossplotFrame.xLog = aBunchOfLogs.TheLog[LogNum].Value
			
	class set_Y_VariablesIn:
		def __init__( self, aBunchOfLogs = mm.BunchOfLogs, logNum = -999 ):
			mainCrossplotFrame.yIdVariable.set( str( aBunchOfLogs.MyLogLabel[LogNum] ) )
			mainCrossplotFrame.yTypeVariable.set( str( aBunchOfLogs.LogType[LogNum] ) )
			mainCrossplotFrame.yLog = aBunchOfLogs.TheLog[LogNum].Value
	
	class getXLog:
		def __init__( self, myParent, aBunchOfLogs = mm.BunchOfLogs() ):
			
			self.getXLogFrame = tk.Frame( myParent )
			
			if aBunchOfLogs.NumLogs == 0:
				
				self.rowCounter = 0
				
				self.noLogsLabel = tk.Label( self.getXLogFrame, text = " There appears to be no logs for this session ", font=("Comic Sans MS",16), fg="Blue" )
				self.noLogsLabel.grid( row = self.rowCounter, column = 0)
				
				self.rowCounter += 1
				
				self.exitButton = tk.Button( self.getXLogFrame,
										text = "Exit", 
										font = ("Comic Sans MS",14), 
										fg="Red", 
										command = myParent.destroy )
				self.exitButton.grid( row = self.rowCounter, columnspan = 3, stick = E+W )
				
			else:
				
				self.rowCounter = 0
				
				self.pickLogLabel = tk.Label( self.getXLogFrame, text = "Pick a log from the list below.", font=("Comic Sans MS",16), fg="Blue" )
				self.pickLogLabel.grid( row = self.rowCounter, columnspan = 3 )
				
				self.rowCounter += 1
				
				self.idLabel = tk.Label( self.getXLogFrame, text = "Log ID", font=("Comic Sans MS",14), fg="Green" )
				self.idLabel.grid( row = self.rowCounter, column = 0 )
				self.typeLabel = tk.Label( self.getXLogFrame, text = "Log Type", font=("Comic Sans MS",14), fg="Green" )
				self.typeLabel.grid( row = self.rowCounter, column = 1 )
				self.unitLabel = tk.Label( self.getXLogFrame, text = "Log Unit", font=("Comic Sans MS",14), fg="Green" )
				self.unitLabel.grid( row = self.rowCounter, column = 2 )
				
				self.rowCounter += 1
				
				for self.i in range( aBunchOfLogs.NumLogs):
					
					self.LogIdLabel = tk.Label( self.getXLogFrame, 
										text = str( aBunchOfLogs.MyLogLabel[self.i] ), 
										font=("Comic Sans MS",14), 
										fg="Black" )
					self.LogTypeRadioButton = tk.Radiobutton( self.getXLogFrame, 
										text = str( aBunchOfLogs.LogType[self.i] ), 
										font=("Comic Sans MS",14), 
										fg="Black", 
										command = lambda: set_X_VariablesIn( aBunchOfLogs, self.i ) )
					self.LogUnitLabel = tk.Label( self.getXLogFrame, 
										text = str( aBunchOfLogs.TheLog[self.i].ValueUnit ), 
										font=("Comic Sans MS",14), 
										fg="Black" )
					
					self.rowCounter += 1
					
				self.exitButton = tk.Button( self.getXLogFrame,
										text = "Exit", 
										font = ("Comic Sans MS",14), 
										fg="Red", 
										command = myParent.destroy )
				self.exitButton.grid( row = self.rowCounter, columnspan = 3, stick = E+W )
				
		class getYLog:
			def __init__( self, myParent, aBunchOfLogs = mm.BunchOfLogs() ):

				self.getYLogFrame = tk.Frame( myParent )

				if aBunchOfLogs.NumLogs == 0:

					self.rowCounter = 0

					self.noLogsLabel = tk.Label( self.getYLogFrame, text = " There appears to be no logs for this session ", font=("Comic Sans MS",16), fg="Blue" )
					self.noLogsLabel.grid( row = self.rowCounter, column = 0)

					self.rowCounter += 1

					self.exitButton = tk.Button( self.getYLogFrame,
											text = "Exit", 
											font = ("Comic Sans MS",14), 
											fg="Red", 
											command = myParent.destroy )
					self.exitButton.grid( row = self.rowCounter, columnspan = 3, stick = E+W )

				else:

					self.rowCounter = 0

					self.pickLogLabel = tk.Label( self.getYLogFrame, text = "Pick a log from the list below.", font=("Comic Sans MS",16), fg="Blue" )
					self.pickLogLabel.grid( row = self.rowCounter, columnspan = 3 )

					self.rowCounter += 1

					self.idLabel = tk.Label( self.getYLogFrame, text = "Log ID", font=("Comic Sans MS",14), fg="Green" )
					self.idLabel.grid( row = self.rowCounter, column = 0 )
					self.typeLabel = tk.Label( self.getYLogFrame, text = "Log Type", font=("Comic Sans MS",14), fg="Green" )
					self.typeLabel.grid( row = self.rowCounter, column = 1 )
					self.unitLabel = tk.Label( self.getYLogFrame, text = "Log Unit", font=("Comic Sans MS",14), fg="Green" )
					self.unitLabel.grid( row = self.rowCounter, column = 2 )
					
					self.rowCounter += 1

					for self.i in range( aBunchOfLogs.NumLogs):

						self.LogIdLabel = tk.Label( self.getYLogFrame, text = str( aBunchOfLogs.MyLogLabel[self.i] ), font=("Comic Sans MS",14), fg="Black" )
						self.LogTypeRadioButton = tk.Radiobutton( self.getYLogFrame, 
												text = str( aBunchOfLogs.LogType[self.i] ), 
												font=("Comic Sans MS",14), 
												fg="Black",
												command = lambda: set_Y_VariablesIn( aBunchOfLogs, self.i ) )
						self.LogUnitLabel = tk.Label( self.getYLogFrame, text = str( aBunchOfLogs.TheLog[self.i].ValueUnit ), font=("Comic Sans MS",14), fg="Black" )

						self.rowCounter += 1

					self.exitButton = tk.Button( self.getYLogFrame,
											text = "Exit", 
											font = ("Comic Sans MS",14), 
											fg="Red", 
											command = myParent.destroy )
					self.exitButton.grid( row = self.rowCounter, columnspan = 3, stick = E+W )
	
	class mainCrossplotFrame:
		def __init__( self, myParent, aBunchOfLogs ):
			
			"""
			
			The main definitions of a cross-plotting object
			
			"""
			# setting up the main variables for this class of object
			mainCrossplotFrame.xIdVariable = tk.StrVar()
			mainCrossplotFrame.xIdVariable.set( "None" )
			mainCrossplotFrame.xTypeVariable = tk.StrVar()
			mainCrossplotFrame.xTypeVariable.set( "None" )
			mainCrossplotFrame.xLog = []
			
			mainCrossplotFrame.yIdVariable = tk.StrVar()
			mainCrossplotFrame.yIdVariable.set( "None" )
			mainCrossplotFrame.yTypeVariable = tk.StrVar()
			mainCrossplotFrame.yTypeVariable.set( "None" )
			mainCrossplotFrame.yLog = []
			
			self.xplotFrame = tk.Frame( myParent )
			self.xplotFrame.pack()
			
			self.rowCounter = 0
			
			self.mainXplotLabel = tk.Label( self.xplotFrame, 
									text = "Log (x,y) plotting function", 
									font=("Comic Sans MS",16), 
									fg="Blue" )
			self.mainXplotLabel.grid( row = self.rowCounter, columnspan = 3, sticky = E+W )
			
			self.rowCounter += 1
			
			self.clickLable = tk.Label( self.xplotFrame,
									text = "Click button below",
									font = ("Comic Sans MS",14), 
									fg = "Green" )
			self.clickLable.grid( row = self.rowCounter, column = 0 )
			self.idLable = tk.Label( row = self.rowCounter,
									text = "Log ID",
									font = ("Comic Sans MS",14), 
									fg = "Green" )
			self.idLable.grid( row = self.rowCounter, column = 1 )
			self.typeLable = tk.Label( row = self.rowCounter,
									text = "Log ID",
									font = ("Comic Sans MS",14), 
									fg = "Green" )
			self.typeLable.grid( row = self.rowCounter, column = 2 )
			
			self.rowCounter += 1
			
			self.xLogChooseButton = tk.Button( self.xplotFrame, 
									text = "Choose X Log", 
									font = ("Comic Sans MS",14), 
									fg = "Black", 
									command = lambda: getXLog( myParent, aBunchOfLogs ) )
			self.xLogChooseButton.grid( row = self.rowCounter, column = 0 )
			self.xLogIdLable = tk.Label(	self.xplotFrame, 
									textvariable = self.xIdVariable, 
									font = ("Comic Sans MS",14), 
									fg="Black" )
			self.xLogIdLable.grid( row = self.rowCounter, column = 1 )
			self.xLogTypeLable = tk.Label(	self.xplotFrame, 
									textvariable = self.xIdVariable, 
									font = ("Comic Sans MS",14), 
									fg="Black" )
			self.xLogTypeLable.grid( row = self.rowCounter, column = 2 )
			
			self.rowCounter += 1
			
			self.yLogChooseButton = tk.Button( self.xplotFrame, 
									text = "Choose Y Log", 
									font = ("Comic Sans MS",14), 
									fg = "Black", 
									command = lambda: getYLog( myParent, aBunchOfLogs ) )
			self.yLogChooseButton.grid( row = self.rowCounter, column = 0 )
			self.yLogIdLable = tk.Label( self.xplotFrame, 
									textvariable = self.yIdVariable, 
									font = ("Comic Sans MS",14), 
									fg="Black" )
			self.yLogIdLable.grid( row = self.rowCounter, column = 1 )
			self.yLogTypeLable = tk.Label( self.xplotFrame, 
									textvariable = self.yIdVariable, 
									font = ("Comic Sans MS",14), 
									fg="Black" )
			self.yLogTypeLable.grid( row = self.rowCounter, column = 2 )
			
			self.rowCounter += 1
			
			self.submitPlotButton = tk.Button( self.xplotFrame,
									text = "Submit",
									font = ("Comic Sans MS",14), 
									fg="Black",
									command = lambda: simPlotXY( self.xLog, self.yLog ) )
			self.submitPlotButton.grid( row = self.rowCounter, columnspan = 3, sticky = E+W )
			
			self.rowCounter += 1
			
			self.exitButton = tk.Button( self.xplotFrame, 
									text = "Exit", 
									font = ("Comic Sans MS",14), 
									fg="Red", 
									command = myParent.destroy )
			self.exitButton.grid( row = self.rowCounter, columnspan = 3, stick = E+W )
			
			self.rowCounter +=1
	
	