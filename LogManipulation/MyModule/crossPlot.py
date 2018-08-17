# 2nd class trial for cross-plotting logs

import Tkinter as tk
import MyModule as mm
# import SimPy.SimPlot as sp
import matplotlib
matplotlib.use('TkAgg')
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
# implement the default mpl key bindings
from matplotlib.backend_bases import key_press_handler
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
import numpy as np

class crossPlotOfLogs:
	"""
	
	A cross plot of two user defined logs.
	A graph of ( x=Log1, y=Log2 ) pairs, assuming the pairing is appropriate.
	One thing to keep in mind is the relavance of the sequential pairing, since the
	two logs may not even have the same depth or correct stratigraphic relationship.
	This problem is generally not an issue if you are comparing different logs from the same
	wellbore referenced to the same depths.
	
	"""
	
	def __init__( self, parent = tk.Tk(), aBunchOfLogs = mm.BunchOfLogs() ):
		"""
		
		Sets up the main control panel for the cross plot
		
		"""
		
		print "+---------------------------------------------+"
		print "|                                             |"
		print "|     Start of log cross-ploting function     |"
		print "|                                             |"
		print "+---------------------------------------------+"
		
		self.parent = parent
		self.localBunchOfLogs = aBunchOfLogs
		
		print "\n Currently you have (" + str(self.localBunchOfLogs.numLogs) + ") logs available.\n"
		
		self.bigBrother = tk.Toplevel( self.parent )
		self.bigBrother.wm_title( "Cross plotting logs function")
		
		if self.localBunchOfLogs.numLogs < 2:
			
			print " You need at least 2 logs to make a cross-plot"
			
			self.rowCounter = 0 #===================================================================

			self.mainXplotLabel = tk.Label( self.bigBrother, 
									text = "Need 2 logs to cross-plot. You have (" + str(self.localBunchOfLogs.numLogs) + ") Please clik Exit button.", 
									font = ("Comic Sans MS",18), 
									fg = "Blue" )
			self.mainXplotLabel.grid( row = self.rowCounter, columnspan = 3, sticky = tk.E+tk.W, padx = 20, pady = 10 )
			
			self.rowCounter += 1 #---------------------------------------------------------
			
			self.exitButton = tk.Button( self.bigBrother, 
									text = "Exit", 
									font = ("Comic Sans MS",16), 
									fg = "Red", 
									command = self.bigBrother.destroy )
			self.exitButton.grid( row = self.rowCounter, column = 1, sticky = tk.E+tk.W, padx = 20, pady = 10 )

			self.rowCounter +=1 #---------------------------------------------------------
		else:
			
			self.xIdVariable = tk.StringVar( self.bigBrother )
			self.xIdVariable.set( "None" )
			self.xTypeVariable = tk.StringVar( self.bigBrother )
			self.xTypeVariable.set( "None")
			self.xLogNumber = tk.IntVar( self.bigBrother )
			self.xLogNumber.set( 0 )
			self.xLog = []

			self.yIdVariable = tk.StringVar( self.bigBrother )
			self.yIdVariable.set( "None" )
			self.yTypeVariable = tk.StringVar( self.bigBrother )
			self.yTypeVariable.set( "None" )
			self.yLogNumber = tk.IntVar( self.bigBrother )
			self.yLogNumber.set( 1 )
			self.yLog = []
		
			self.rowCounter = 0 #===================================================================
		
			self.mainXplotLabel = tk.Label( self.bigBrother, 
								text = "Log (x,y) plotting function", 
								font = ("Comic Sans MS",18), 
								fg = "Blue" )
			self.mainXplotLabel.grid( row = self.rowCounter, columnspan = 3, sticky = tk.E+tk.W, padx = 20, pady = 10 )
			
			self.rowCounter += 1 #---------------------------------------------------------
			
			self.clickLabel = tk.Label( self.bigBrother,
								text = "Click button below",
								font = ("Comic Sans MS",16), 
								fg = "Green" )
			self.clickLabel.grid( row = self.rowCounter, column = 0, padx = 10, pady = 5 )
		
			self.idLabel = tk.Label( self.bigBrother,
								text = "Log ID",
								font = ("Comic Sans MS",16), 
								fg = "Green" )
			self.idLabel.grid( row = self.rowCounter, column = 1, padx = 10, pady = 5 )
			
			self.typeLabel = tk.Label( self.bigBrother,
								text = "Log Type",
								font = ("Comic Sans MS",16), 
								fg = "Green" )
			self.typeLabel.grid( row = self.rowCounter, column = 2, padx = 10, pady = 5 )
		
			self.rowCounter += 1 #---------------------------------------------------------
		
			self.xLogChooseButton = tk.Button( self.bigBrother, 
								text = "Choose X Log", 
								font = ("Comic Sans MS",14), 
								fg = "Black", 
								command = self.getXLog )
			self.xLogChooseButton.grid( row = self.rowCounter, column = 0 )
		
			print " Current chosen x-log ID  : "+str( self.xIdVariable.get() )
			self.xLogIdLabel = tk.Label( self.bigBrother, 
								textvariable = self.xIdVariable,
								font = ("Comic Sans MS",14), 
								fg = "Black" )
			self.xLogIdLabel.grid( row = self.rowCounter, column = 1 )
		
			print " Current chosen x-log Type: "+str( self.xTypeVariable.get() )+"\n"
			self.xlogTypeLabel = tk.Label(	self.bigBrother, 
								textvariable = self.xTypeVariable, 
								font = ("Comic Sans MS",14), 
								fg = "Black" )
			self.xlogTypeLabel.grid( row = self.rowCounter, column = 2 )
		
			self.rowCounter += 1 #---------------------------------------------------------
		
			self.yLogChooseButton = tk.Button( self.bigBrother, 
								text = "Choose Y Log", 
								font = ("Comic Sans MS",14), 
								fg = "Black", 
								command = self.getYLog )
			self.yLogChooseButton.grid( row = self.rowCounter, column = 0 )
		
			print " Current chosen y-log ID  : "+str( self.yIdVariable.get() )
			self.yLogIdLabel = tk.Label( self.bigBrother, 
								textvariable = self.yIdVariable,
								font = ("Comic Sans MS",14), 
								fg = "Black" )
			self.yLogIdLabel.grid( row = self.rowCounter, column = 1 )
		
			print " Current chosen y-log Type: "+str( self.yTypeVariable.get() )+"\n"
			self.ylogTypeLabel = tk.Label( self.bigBrother, 
								textvariable = self.yTypeVariable,
								font = ("Comic Sans MS",14), 
								fg = "Black" )
			self.ylogTypeLabel.grid( row = self.rowCounter, column = 2 )
		
			self.rowCounter += 1 #---------------------------------------------------------
			
			self.plotType = tk.StringVar( self.bigBrother )
			self.plotType.set( "line" )
			self.plotTypeLabel = tk.Label( self.bigBrother,
								text = "Choose a plot type below",
								font = ("Comic Sans MS",16), 
								fg = "Green")
			self.plotTypeLabel.grid( row = self.rowCounter, columnspan = 3, sticky = tk.E+tk.W, padx = 20, pady = 10 )
			
			self.rowCounter += 1 #---------------------------------------------------------
			
			self.lineButton = tk.Radiobutton( self.bigBrother,
								text = "Line graph",
								variable = self.plotType,
								value = "line",
								font = ("Comic Sans MS",14), 
								fg = "Black" )
			self.lineButton.grid( row = self.rowCounter, column = 0, sticky = tk.W+tk.E, padx = 20 )
			
			self.pointButton = tk.Radiobutton( self.bigBrother,
								text = "Scatter plot",
								variable = self.plotType,
								value = "point",
								font = ("Comic Sans MS",14), 
								fg = "Black" )
			self.pointButton.grid( row = self.rowCounter, column = 1, sticky = tk.W+tk.E, padx = 20  )
			
			self.stepButton = tk.Radiobutton( self.bigBrother,
								text = "Step graph",
								variable = self.plotType,
								value = "step",
								font = ("Comic Sans MS",14), 
								fg = "Black" )
			self.stepButton.grid( row = self.rowCounter, column = 2, sticky = tk.W+tk.E, padx = 20  )
								
			self.rowCounter += 1 #---------------------------------------------------------
			
			self.axisScaleLabel = tk.Label( self.bigBrother,
								text = "Use sliders below to set relative scale of x and y axis (pixels)",
								font = ("Comic Sans MS",16), 
								fg = "Green")
			self.axisScaleLabel.grid( row = self.rowCounter, columnspan = 3, sticky = tk.E+tk.W, padx = 20, pady = 10 )
			
			self.rowCounter += 1 #---------------------------------------------------------
			
			self.ySlideLabel = tk.Label( self.bigBrother, 
								text = "Y-axis slider",
								font = ("Comic Sans MS",14), 
								fg = "Black" )
			self.ySlideLabel.grid( row = self.rowCounter, column = 0, sticky = tk.E + tk.W )
			
			self.xSlideLabel = tk.Label( self.bigBrother, 
								text = "X-axis slider",
								font = ("Comic Sans MS",14), 
								fg = "Black" )
			self.xSlideLabel.grid( row = self.rowCounter, column = 1, sticky = tk.E + tk.W )
			
			self.rowCounter += 1 #---------------------------------------------------------
			
			self.minNumberOfPixels = 100
			self.maxNumberOfPixels = 1000
			self.xScale = tk.IntVar( self.bigBrother )
			self.xScale.set(500)
			self.yScale = tk.IntVar( self.bigBrother )
			self.yScale.set(500)
			
			self.yScaleSlider = tk.Scale( self.bigBrother,
			 					variable = self.yScale,
								from_ = self.minNumberOfPixels, 
								to = self.maxNumberOfPixels )
			self.yScaleSlider.set(500)
			self.yScaleSlider.grid( row = self.rowCounter, column = 0, sticky = tk.E+tk.W, padx=20 )
			
			self.xScaleSlider = tk.Scale( self.bigBrother,
			 					variable = self.xScale,
								from_ = self.minNumberOfPixels, 
								to = self.maxNumberOfPixels,
								orient = tk.HORIZONTAL )
			self.xScaleSlider.set(500)
			self.xScaleSlider.grid( row = self.rowCounter, column = 1, sticky = tk.E+tk.W+tk.S, padx=20 )
			
			self.rowCounter += 1 #---------------------------------------------------------
			
			self.submitPlotButton = tk.Button( self.bigBrother,
								text = "Submit",
								font = ("Comic Sans MS",16), 
								fg = "Black",
								command = self.pyPlotXY )
			self.submitPlotButton.grid( row = self.rowCounter, column = 1, sticky = tk.E+tk.W, padx = 20, pady = 5 )
		
			self.rowCounter += 1 #---------------------------------------------------------
		
			self.exitButton = tk.Button( self.bigBrother, 
								text = "Exit", 
								font = ("Comic Sans MS",16), 
								fg = "Red", 
								command = self.bigBrother.destroy )
			self.exitButton.grid( row = self.rowCounter, column = 1, sticky = tk.E+tk.W, padx = 20, pady = 10 )
		
			self.rowCounter +=1 #---------------------------------------------------------
		
		self.bigBrother.mainloop()
		
	def getXLog( self ):
		"""
		Function to allow the selection of a log from a bunch of logs to be used
		as the x value for cross plotting purposes
		"""
		
		self.littleBrotherX = tk.Toplevel( self.parent )
		self.littleBrotherX.wm_title( "Choose log for values along the x-axis" )
		
		if self.localBunchOfLogs.numLogs == 0:
			
			self.rowCounter = 0 #===================================================================
			
			self.noLogsLabel = tk.Label( self.littleBrotherX, 
									text = " There appears to be no logs for this session ", 
									font = ("Comic Sans MS",16), 
									fg = "Blue" )
			self.noLogsLabel.grid( row = self.rowCounter, column = 0, padx = 20, pady = 10 )
			
			self.rowCounter += 1 #---------------------------------------------------------
			
			self.exitButton = tk.Button( self.littleBrotherX,
									text = "Exit", 
									font = ("Comic Sans MS",14), 
									fg = "Red", 
									command = self.littleBrotherX.destroy )
			self.exitButton.grid( row = self.rowCounter, columnspan = 0 )
			
			self.rowCounter += 1 #---------------------------------------------------------
			
		else:
			
			self.rowCounter = 0 #===================================================================
			
			self.pickLogLabel = tk.Label( self.littleBrotherX, 
									text = "Pick a log for x-values from the list below.", 
									font = ("Comic Sans MS",16), 
									fg = "Blue" )
			self.pickLogLabel.grid( row = self.rowCounter, columnspan = 3, padx = 20, pady = 10 )
			
			self.rowCounter += 1 #---------------------------------------------------------
			
			self.idLabel = tk.Label( self.littleBrotherX, 
									text = "Log ID", 
									font = ("Comic Sans MS",14), 
									fg = "Green" )
			self.idLabel.grid( row = self.rowCounter, column = 0, padx = 10, pady = 5 )
			
			self.typeLabel = tk.Label( self.littleBrotherX, 
									text = "Log Type", 
									font = ("Comic Sans MS",14), 
									fg = "Green" )
			self.typeLabel.grid( row = self.rowCounter, column = 1, padx = 10, pady = 5 )
			
			self.unitLabel = tk.Label( self.littleBrotherX, 
									text = "Log Unit", 
									font = ("Comic Sans MS",14), 
									fg = "Green" )
			self.unitLabel.grid( row = self.rowCounter, column = 2, padx = 10, pady = 5 )
			
			self.rowCounter += 1 #---------------------------------------------------------
			
			for i in range( self.localBunchOfLogs.numLogs):
				
				self.LogIdLabel = tk.Label( self.littleBrotherX, 
									text = str( self.localBunchOfLogs.myLogLabel[i] ), 
									font = ("Comic Sans MS",14), 
									fg = "Black" )
				self.LogIdLabel.grid( row = self.rowCounter, column = 0 )
				
				self.logTypeRadioButton = tk.Radiobutton( self.littleBrotherX, 
									text = str( self.localBunchOfLogs.logType[i] ),
									variable = self.xLogNumber,
									value = int(i), 
									font = ("Comic Sans MS",14), 
									fg = "Black", 
									command = self.set_X_Variables )
				self.logTypeRadioButton.grid( row = self.rowCounter, column = 1, sticky = tk.W, padx = 10 )
				
				self.LogUnitLabel = tk.Label( self.littleBrotherX, 
									text = str( self.localBunchOfLogs.theLog[i].ValueUnit ), 
									font = ("Comic Sans MS",14), 
									fg = "Black" )
				self.LogUnitLabel.grid( row = self.rowCounter, column = 2, padx = 10 )
				
				self.rowCounter += 1 #---------------------------------------------------------
				
			self.exitButton = tk.Button( self.littleBrotherX,
									text = "Exit", 
									font = ("Comic Sans MS",14), 
									fg = "Red", 
									command = self.littleBrotherX.destroy )
			self.exitButton.grid( row = self.rowCounter, column = 1, sticky = tk.E+tk.W, padx = 20, pady = 10 )
			
	def getYLog( self ):
		"""
		Function to allow the selection of a log from a bunch of logs to be used
		as the y value for cross plotting purposes
		"""

		self.littleBrotherY = tk.Toplevel( self.parent )
		self.littleBrotherY.wm_title( "Choose log of y-axis" )
		
		if self.localBunchOfLogs.numLogs == 0:

			self.rowCounter = 0 #===================================================================

			self.noLogsLabel = tk.Label( self.littleBrotherY, 
									text = " There appears to be no logs for this session ", 
									font = ("Comic Sans MS",16), 
									fg = "Blue" )
			self.noLogsLabel.grid( row = self.rowCounter, column = 0, padx = 20, pady = 10)

			self.rowCounter += 1 #---------------------------------------------------------

			self.exitButton = tk.Button( self.littleBrotherY,
									text = "Exit", 
									font = ("Comic Sans MS",14), 
									fg = "Red", 
									command = self.littleBrother.destroy )
			self.exitButton.grid( row = self.rowCounter, columnspan = 3, sticky = tk.E+tk.W, padx = 20, pady = 10 )

			self.rowCounter += 1 #---------------------------------------------------------

		else:

			self.rowCounter = 0 #===================================================================

			self.pickLogLabel = tk.Label( self.littleBrotherY, 
									text = "Pick a log from the list below.", 
									font = ("Comic Sans MS",16), 
									fg = "Blue" )
			self.pickLogLabel.grid( row = self.rowCounter, columnspan = 3, padx = 20, pady = 10 )

			self.rowCounter += 1 #---------------------------------------------------------

			self.idLabel = tk.Label( self.littleBrotherY, 
									text = "Log ID", 
									font = ("Comic Sans MS",14), 
									fg = "Green" )
			self.idLabel.grid( row = self.rowCounter, column = 0, padx = 10, pady = 5 )

			self.typeLabel = tk.Label( self.littleBrotherY, 
									text = "Log Type", 
									font = ("Comic Sans MS",14), 
									fg = "Green" )
			self.typeLabel.grid( row = self.rowCounter, column = 1, padx = 10, pady = 5 )
			self.unitLabel = tk.Label( self.littleBrotherY, 
									text = "Log Unit", 
									font = ("Comic Sans MS",14), 
									fg = "Green" )
			self.unitLabel.grid( row = self.rowCounter, column = 2, padx = 10, pady = 5 )

			self.rowCounter += 1 #---------------------------------------------------------

			for i in range( self.localBunchOfLogs.numLogs):

				self.LogIdLabel = tk.Label( self.littleBrotherY, 
									text = str( self.localBunchOfLogs.myLogLabel[i] ),
									font = ("Comic Sans MS",14), 
									fg = "Black" )
				self.LogIdLabel.grid( row = self.rowCounter, column = 0, padx = 10 )

				self.logTypeRadioButton = tk.Radiobutton( self.littleBrotherY, 
									text = str( self.localBunchOfLogs.logType[i] ),
									variable = self.yLogNumber,
									value = int(i), 
									font = ("Comic Sans MS",14), 
									fg = "Black", 
									command = self.set_Y_Variables )
				self.logTypeRadioButton.grid( row = self.rowCounter, column = 1, sticky = tk.W, padx = 10 )

				self.LogUnitLabel = tk.Label( self.littleBrotherY, 
									text = str( self.localBunchOfLogs.theLog[i].ValueUnit ), 
									font = ("Comic Sans MS",14), 
									fg = "Black" )
				self.LogUnitLabel.grid( row = self.rowCounter, column = 2, padx = 10 )

				self.rowCounter += 1 #---------------------------------------------------------

			self.exitButton = tk.Button( self.littleBrotherY,
									text = "Exit", 
									font = ("Comic Sans MS",14), 
									fg = "Red", 
									command = self.littleBrotherY.destroy )
			self.exitButton.grid( row = self.rowCounter, columnspan = 3, sticky = tk.E+tk.W )
			
		# self.littleBrotherY.mainloop()
			
	def set_X_Variables( self ):
		"""
		Set class variables associated with a cross plot for x-axis.
		"""
		print "\n You have choosen log number (" + str(self.xLogNumber.get()) + ") for x values\n"
		print "      x-Log ID  : " + str( self.localBunchOfLogs.myLogLabel[self.xLogNumber.get()] )
		print "      x-Log Type: " + str( self.localBunchOfLogs.logType[self.xLogNumber.get()] )
		print "      x-Log Unit: " + str( self.localBunchOfLogs.theLog[self.xLogNumber.get()].ValueUnit ) + "\n"
		
		self.xIdVariable.set( str( self.localBunchOfLogs.myLogLabel[self.xLogNumber.get()] ) )
		self.xTypeVariable.set( str( self.localBunchOfLogs.logType[self.xLogNumber.get()] ) )
		self.xLog = self.localBunchOfLogs.theLog[self.xLogNumber.get()].Value
		
	def set_Y_Variables( self ):
		"""
		Set class variables associated with a cross plot for y-axis.
		"""
		print "\n You have choosen log number (" + str(self.yLogNumber.get()) + ") for x values\n"
		print "      y-Log ID  : " + str( self.localBunchOfLogs.myLogLabel[self.yLogNumber.get()] )
		print "      y-Log Type: " + str( self.localBunchOfLogs.logType[self.yLogNumber.get()] )
		print "      y-Log Unit: " + str( self.localBunchOfLogs.theLog[self.yLogNumber.get()].ValueUnit ) + "\n"
		
		self.yIdVariable.set( str( self.localBunchOfLogs.myLogLabel[self.yLogNumber.get()] ) )
		self.yTypeVariable.set( str( self.localBunchOfLogs.logType[self.yLogNumber.get()] ) )
		self.yLog = self.localBunchOfLogs.theLog[self.yLogNumber.get()].Value
		
	# def simPlotXY( self ):
	# 		"""
	# 		Simple (x.y) plot using basic plotting package SimPy.SimPlot
	# 		"""
	# 		print"\n *** Entering plotting function with choosen logs. ***\n"
	# 		print "     X-axis set to : " + str(self.xScale) + " pixels"
	# 		print "     Y-axis set to : " + str(self.yScale) + " pixels"
	# 		
	# 		self.zip_Logs = zip( self.xLog, self.yLog )
	# 		
	# 		if self.plotType.get() == "line":
	# 			
	# 			plt=sp.SimPlot()
	# 			plt.plotLine( self.zip_Logs, 
	# 								color = "blue",
	# 								windowsize = ( self.xScale.get(), self.yScale.get() ) )
	# 			plt.mainloop()
	# 			
	# 		elif self.plotType.get() == "point":
	# 			
	# 			plt=sp.SimPlot()
	# 			plt.plotScatter( self.zip_Logs, 
	# 								size = 1, 
	# 								marker = "circle", 
	# 								color = "blue",
	# 								windowsize = ( self.xScale.get(), self.yScale.get() ) )
	# 			plt.mainloop()
	# 			
	# 		elif self.plotType.get() == "step":
	# 			
	# 			plt=sp.SimPlot()
	# 			plt.plotStep( self.zip_Logs, 
	# 								color = "blue",
	# 								windowsize = ( self.xScale.get(), self.yScale.get() ) )
	# 			plt.mainloop()
			
	def on_key_event( self, event ):
		print('you pressed %s'%event.key)
		key_press_handler(event, self.canvas, self.toolbar)
	
	def _quit( self ):
		self.sister.quit()
		self.sister.destroy()
	
	def pyPlotXY( self ):
		"""
		This ploting routine is ment to eventually replace simPlotXY
		"""
		self.sister = tk.Tk()
		self.sister.wm_title("Log Cross-Plotting.")
		
		self.dpi = 70
		self.xLengthInches = self.xScale.get()/self.dpi
		self.yLengthInches = self.yScale.get()/self.dpi
		

		print"\n *** Entering plotting function with choosen logs. ***\n"
		print "     X-axis set to : " + str(self.xLengthInches) + " inches"
		print "     Y-axis set to : " + str(self.yLengthInches) + " inches"
		
		self.crossPlot = Figure( figsize=(self.xLengthInches,self.yLengthInches), 
								dpi=self.dpi )
								
		self.graph = self.crossPlot.add_subplot(1,1,1)
		
		if self.plotType.get() == "line":
			self.graph.plot( self.xLog, self.yLog )
		elif self.plotType.get() == "point":
			self.graph.plot( self.xLog, self.yLog, 'ro' )
		elif self.plotType.get() == "step":
			self.graph.step( self.xLog, self.yLog )
			
		self.graph.set_xlabel( str( self.localBunchOfLogs.theLog[self.xLogNumber.get()].ValueUnit ) )	
		self.graph.set_ylabel( str( self.localBunchOfLogs.theLog[self.yLogNumber.get()].ValueUnit ) )
		self.graph.set_title( str( self.localBunchOfLogs.logType[self.xLogNumber.get()] ) + 
												" vs. " + 
												str( self.localBunchOfLogs.logType[self.yLogNumber.get()] ) )
										
		self.canvas = FigureCanvasTkAgg( self.crossPlot, master = self.sister )
		self.canvas.show()
		self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
		
		self.toolbar = NavigationToolbar2TkAgg( self.canvas, self.sister )
		self.toolbar.update()
		self.canvas._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=1)
		
		
		self.canvas.mpl_connect('key_press_event', self.on_key_event)
		
		
			
		self.exitButton = tk.Button( master = self.sister, 
							text = "Quit", 
							command = self._quit )
		self.exitButton.pack( side = tk.BOTTOM )
		
		self.sister.mainloop()
