import Tkinter
import tkMessageBox

top = Tkinter.Tk()

def helloCallBack():
   tkMessageBox.showinfo( "Hello Python", "Hello World you know not what is to be coming next")

B = Tkinter.Button(top, width=7, padx=20, pady=20, bg ="black", activebackground="green", activeforeground="yellow", fg="white", text ="Hello", command = helloCallBack)

B.pack()
top.mainloop()