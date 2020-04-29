import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
try :
    import subprocess
except:
    print("install subprocess")
    exit(0)
try :
    import Tkinter as Tk
except ImportError :
    import tkinter as Tk

import os

cwd = os.getcwd()
# https://www.blog.pythonlibrary.org/2012/07/26/tkinter-how-to-show-hide-a-window/

class showPlotClass(Tk.Toplevel):
    def __init__(self, original):
        """Constructor"""
        self.original_frame = original
        self.data = original.data
        self.Dict = original.Dict
        Tk.Toplevel.__init__(self)
        self.geometry("800x600")
        self.title("Show Plot")
        

        self.showPlotMenu()
       
        btn = Tk.Button(self, text="Close", command=self.onClose)
        btn.pack()          

    def onClose(self):
        self.destroy()
        self.original_frame.show()

        #plotting menu
    def showPlotMenu(self):
 
        total = 0
        for col in (self.data).columns:
            self.Dict[total] = col
            if (total != 0):
                Tk.Button(self,text = col,command=(lambda total = total:(self.showPlot(total)))).pack()
            total += 1      
    
    def showPlot(self, choice_show):
        plt.close()
        # figure = plt.Figure()
        # ax[choice_show] = plt.gca()
        print(choice_show)
        print(self.Dict.get(choice_show))
        # (self.data).plot(kind = 'line', x = 't', y = ((self.Dict).get(choice_show)), ax = ax[choice_show])
        (self.data).plot(kind = 'line', x = 't', y = ((self.Dict).get(choice_show)))
        plt.show()

class addFrame(Tk.Toplevel):
    def __init__(self, original, frameName):
        self.original_frame = original
        Tk.Toplevel.__init__(self)
        self.geometry("400x300")
        self.title(frameName) 
        self.Row = 0

    def onClose(self):
        self.destroy()
        self.original_frame.show()

    def getRow(self):
        self.Row += 1
        return (self.Row - 1)    


class addComponentClass(Tk.Toplevel):
    def __init__(self, original, inputFile):
        self.original_frame = original
        self.data = original.data
        self.Dict = original.Dict
        self.Row = 0
        Tk.Toplevel.__init__(self)
        self.geometry("800x600")
        self.title("addComponentClass")        
        self.file = open(original.inputFile, "w+")

        btnres =Tk.Button(self, text="Add Resistance", command=self.addResistance)
        btnres.grid(row = self.getRow()) 

        btncap =Tk.Button(self, text="Add Capacitance", command=self.addCapacitance)
        btncap.grid(row = self.getRow())

        btnind =Tk.Button(self, text="Add Inductance", command=self.addInductance)
        btnind.grid(row = self.getRow())                           

        btnvolt =Tk.Button(self, text="Add Constant Voltage Source", command=self.addVoltageSource)
        btnvolt.grid(row = self.getRow())                                   

        btncurr =Tk.Button(self, text="Add Constant Current Source", command=self.addCurrentSource)
        btncurr.grid(row = self.getRow())         

        # btnvoltAc =Tk.Button(self, text="Add AC Voltage Source", command=self.addACVoltageSource)
        # btnvoltAc.grid(row = self.getRow())        

        btncurrAc =Tk.Button(self, text="Add AC Current Source", command=self.addACCurrentSource)
        btncurrAc.grid(row = self.getRow())                 

        btnvoltAc =Tk.Button(self, text="Add AC Voltage Source", command=self.addACVoltageSource)
        btnvoltAc.grid(row = self.getRow())           
                                 
       
        btn = Tk.Button(self, text="Save", command=self.onClose)
        btn.grid(row = self.getRow())          

    def hide(self):
        self.withdraw()
    def show(self):
        """"""
        self.update()
        self.deiconify()    

    def onClose(self):
        self.destroy()
        self.original_frame.show()

    def getRow(self):
        self.Row += 1
        return (self.Row - 1)

    #Add resistance Component
    def addResistance(self):
        self.hide()
        otherFrame = addFrame(self, "addResistance")
        resValueEntry = Tk.Entry(otherFrame)
        
        resNode1Entry = Tk.Entry(otherFrame)
        
        resNode2Entry = Tk.Entry(otherFrame)

        rowStart = otherFrame.getRow()
        resValueLabel=Tk.Label(otherFrame,text = "Nilai resistor(Ohm):")        
        resNode1Label=Tk.Label(otherFrame,text = "Nilai Node 1:")
        resNode2Label=Tk.Label(otherFrame,text = "Nilai Node 2:")

        resValueLabel.grid(row=rowStart, column = 1)
        resNode1Label.grid(row=rowStart+1, column = 1)        
        resNode2Label.grid(row=rowStart+2, column = 1)

        resNode1Entry.grid(row=rowStart+1, column = 2)        
        resValueEntry.grid(row=rowStart, column = 2)
        resNode2Entry.grid(row=rowStart+2, column = 2)    

        submit = Tk.Button(otherFrame, text = "Submit")
        cancel = Tk.Button(otherFrame, text = "Cancel")
        submit['command'] =lambda binst =submit, binst2 = cancel:self.putResistance(1,resValueEntry,
                resNode1Entry,resNode2Entry, resValueLabel, resNode1Label, resNode2Label, binst, binst2, otherFrame)
        cancel['command'] = lambda binst = cancel,  binst2 = submit:self.putResistance(0,resValueEntry,
                resNode1Entry,resNode2Entry, resValueLabel, resNode1Label, resNode2Label, binst, binst2, otherFrame)
        cancel.grid(row = rowStart+4, column = 4)
        submit.grid(row = rowStart+4, column = 5)

    def putResistance(self, state ,resValueEntry, resNode1Entry, resNode2Entry, resValueLabel, resNode1Label, resNode2Label,binst, binst2, otherFrame):
        if state == 1:
            self.file.write("R")
            self.file.write("\n")
            self.file.write(str(resValueEntry.get()))
            self.file.write("\n")       
            self.file.write(str(resNode1Entry.get()))
            self.file.write("\n")       
            self.file.write(str(resNode2Entry.get()))
            self.file.write("\n")
   
        otherFrame.onClose()
        

    #Add Capacitor Component
    def addCapacitance(self):
        self.hide()
        otherFrame = addFrame(self, "addResistance")

        capValueEntry = Tk.Entry(otherFrame)
        
        capNode1Entry = Tk.Entry(otherFrame)
        
        capNode2Entry = Tk.Entry(otherFrame)

        capVoltEntry = Tk.Entry(otherFrame)


        rowStart = otherFrame.getRow()
        capValueLabel=Tk.Label(otherFrame,text = "Nilai capacitor:")        
        capNode1Label=Tk.Label(otherFrame,text = "Nilai Node 1:")
        capNode2Label=Tk.Label(otherFrame,text = "Nilai Node 2:")
        capVoltLabel=Tk.Label(otherFrame,text = "Tegangan awal kapasitor:")

        capValueLabel.grid(row=rowStart, column = 1)
        capNode1Label.grid(row=rowStart+1, column = 1)        
        capNode2Label.grid(row=rowStart+2, column = 1)
        capVoltLabel.grid(row = rowStart+3, column = 1)

        capNode1Entry.grid(row=rowStart+1, column = 2)        
        capValueEntry.grid(row=rowStart, column = 2)
        capNode2Entry.grid(row=rowStart+2, column = 2)    
        capVoltEntry.grid(row=rowStart+3, column = 2)

        submit = Tk.Button(otherFrame, text = "Submit")
        cancel = Tk.Button(otherFrame, text = "Cancel")
        submit['command'] =lambda binst =submit, binst2 = cancel:self.putCapacitance(1,capValueEntry,
                capNode1Entry,capNode2Entry, capValueLabel, capNode1Label, capNode2Label,capVoltEntry ,binst, binst2, otherFrame)
        cancel['command'] = lambda binst = cancel,  binst2 = submit:self.putCapacitance(0,capValueEntry,
                capNode1Entry,capNode2Entry, capValueLabel, capNode1Label, capNode2Label,capVoltEntry ,binst, binst2, otherFrame)
        cancel.grid(row = rowStart+4, column = 4)
        submit.grid(row = rowStart+4, column = 5)

    def putCapacitance(self, state ,capValueEntry, capNode1Entry, capNode2Entry, capValueLabel, capNode1Label, capNode2Label,capVoltEntry,binst, binst2, otherFrame):
        if state == 1:
            self.file.write("C")
            self.file.write("\n")
            self.file.write(str(capValueEntry.get()))
            self.file.write("\n")       
            self.file.write(str(capNode1Entry.get()))
            self.file.write("\n")       
            self.file.write(str(capNode2Entry.get()))
            self.file.write("\n")
            self.file.write(str(capVoltEntry.get()))
            self.file.write("\n")

        capValueEntry.destroy()
        capNode1Entry.destroy()     
        capNode2Entry.destroy()  
        capValueLabel.destroy()
        capNode1Label.destroy()     
        capNode2Label.destroy()          

        binst.destroy()   
        binst2.destroy() 
        otherFrame.onClose()

    def addInductance(self):
        self.hide()
        otherFrame = addFrame(self, "Add Inductance")

        indValueEntry = Tk.Entry(otherFrame)
        indNode1Entry = Tk.Entry(otherFrame)
        indNode2Entry = Tk.Entry(otherFrame)
        indCurrEntry  = Tk.Entry(otherFrame)

        rowStart = otherFrame.getRow()
        indValueLabel=Tk.Label(otherFrame,text = "Nilai induktor(H):")        
        indNode1Label=Tk.Label(otherFrame,text = "Nilai Node 1:")
        indNode2Label=Tk.Label(otherFrame,text = "Nilai Node 2:")
        indCurrLabel=Tk.Label(otherFrame,text = "Nilai Arus Awal(A):")

        indValueLabel.grid(row=rowStart, column = 1)
        indNode1Label.grid(row=rowStart+1, column = 1)        
        indNode2Label.grid(row=rowStart+2, column = 1)
        indCurrLabel.grid(row=rowStart+3, column = 1)

        indNode1Entry.grid(row=rowStart+1, column = 2)        
        indValueEntry.grid(row=rowStart, column = 2)
        indNode2Entry.grid(row=rowStart+2, column = 2)    
        indCurrEntry.grid(row=rowStart+3, column = 2)    

        submit = Tk.Button(otherFrame, text = "Submit")
        cancel = Tk.Button(otherFrame, text = "Cancel")
        submit['command'] =lambda binst =submit, binst2 = cancel:self.putInductance(1,indValueEntry,
                indNode1Entry,indNode2Entry, indValueLabel, indNode1Label, indNode2Label,indCurrEntry ,binst, binst2, otherFrame)
        cancel['command'] = lambda binst = cancel,  binst2 = submit:self.putInductance(0,indValueEntry,
                indNode1Entry,indNode2Entry, indValueLabel, indNode1Label, indNode2Label,indCurrEntry ,binst, binst2, otherFrame)
        cancel.grid(row = rowStart+4, column = 4)
        submit.grid(row = rowStart+4, column = 5)


    def putInductance(self, state ,indValueEntry, indNode1Entry, indNode2Entry, indValueLabel, indNode1Label, indNode2Label,indCurrEntry,binst, binst2, otherFrame):
        if state == 1:
            self.file.write("L")
            self.file.write("\n")
            self.file.write(str(indValueEntry.get()))
            self.file.write("\n")       
            self.file.write(str(indNode1Entry.get()))
            self.file.write("\n")       
            self.file.write(str(indNode2Entry.get()))
            self.file.write("\n")
            self.file.write(str(indCurrEntry.get()))
            self.file.write("\n")

        indValueEntry.destroy()
        indNode1Entry.destroy()     
        indNode2Entry.destroy()  
        indValueLabel.destroy()
        indNode1Label.destroy()     
        indNode2Label.destroy()          

        binst.destroy()   
        binst2.destroy()  
        otherFrame.onClose()

    def addVoltageSource(self):
        self.hide()
        otherFrame = addFrame(self, "Add Constant Voltage Source")
        voltValueEntry = Tk.Entry(otherFrame)
        
        voltNode1Entry = Tk.Entry(otherFrame)
        
        voltNode2Entry = Tk.Entry(otherFrame)

        rowStart = otherFrame.getRow()
        voltValueLabel=Tk.Label(otherFrame,text = "Nilai Sumber Tegangan:")        
        voltNode1Label=Tk.Label(otherFrame,text = "Nilai Node 1:")
        voltNode2Label=Tk.Label(otherFrame,text = "Nilai Node 2:")

        voltValueLabel.grid(row=rowStart, column = 1)
        voltNode1Label.grid(row=rowStart+1, column = 1)        
        voltNode2Label.grid(row=rowStart+2, column = 1)

        voltNode1Entry.grid(row=rowStart+1, column = 2)        
        voltValueEntry.grid(row=rowStart, column = 2)
        voltNode2Entry.grid(row=rowStart+2, column = 2)    

        submit = Tk.Button(otherFrame, text = "Submit")
        cancel = Tk.Button(otherFrame, text = "Cancel")
        submit['command'] =lambda binst =submit, binst2 = cancel:self.putVoltageSource(1,voltValueEntry,
                voltNode1Entry,voltNode2Entry, voltValueLabel, voltNode1Label, voltNode2Label, binst, binst2, otherFrame)
        cancel['command'] = lambda binst = cancel,  binst2 = submit:self.putVoltageSource(0,voltValueEntry,
                voltNode1Entry,voltNode2Entry, voltValueLabel, voltNode1Label, voltNode2Label, binst, binst2, otherFrame)
        cancel.grid(row = rowStart+4, column = 4)
        submit.grid(row = rowStart+4, column = 5)

    def putVoltageSource(self, state ,voltValueEntry, voltNode1Entry, voltNode2Entry, voltValueLabel, voltNode1Label, voltNode2Label,binst, binst2, otherFrame):
        if state == 1:
            self.file.write("V")
            self.file.write("\n")
            self.file.write(str(voltValueEntry.get()))
            self.file.write("\n")       
            self.file.write(str(voltNode1Entry.get()))
            self.file.write("\n")       
            self.file.write(str(voltNode2Entry.get()))
            self.file.write("\n")

        otherFrame.onClose()

    def addCurrentSource(self):
        self.hide()
        otherFrame = addFrame(self, "Add Current Source")
        currValueEntry = Tk.Entry(otherFrame)
        
        currNode1Entry = Tk.Entry(otherFrame)
        
        currNode2Entry = Tk.Entry(otherFrame)

        rowStart = otherFrame.getRow()
        currValueLabel=Tk.Label(otherFrame,text = "Nilai Sumber Arus:")        
        currNode1Label=Tk.Label(otherFrame,text = "Nilai Node 1:")
        currNode2Label=Tk.Label(otherFrame,text = "Nilai Node 2:")

        currValueLabel.grid(row=rowStart, column = 1)
        currNode1Label.grid(row=rowStart+1, column = 1)        
        currNode2Label.grid(row=rowStart+2, column = 1)

        currNode1Entry.grid(row=rowStart+1, column = 2)        
        currValueEntry.grid(row=rowStart, column = 2)
        currNode2Entry.grid(row=rowStart+2, column = 2)    

        submit = Tk.Button(otherFrame, text = "Submit")
        cancel = Tk.Button(otherFrame, text = "Cancel")
        submit['command'] =lambda binst =submit, binst2 = cancel:self.putCurrentSource(1,currValueEntry,
                currNode1Entry,currNode2Entry, currValueLabel, currNode1Label, currNode2Label, binst, binst2, otherFrame)
        cancel['command'] = lambda binst = cancel,  binst2 = submit:self.putCurrentSource(0,currValueEntry,
                currNode1Entry,currNode2Entry, currValueLabel, currNode1Label, currNode2Label, binst, binst2, otherFrame)
        cancel.grid(row = rowStart+4, column = 4)
        submit.grid(row = rowStart+4, column = 5)

    def putCurrentSource(self, state ,currValueEntry, currNode1Entry, currNode2Entry, currValueLabel, currNode1Label, currNode2Label,binst, binst2, otherFrame):
        if state == 1:
            self.file.write("I")
            self.file.write("\n")
            self.file.write(str(currValueEntry.get()))
            self.file.write("\n")       
            self.file.write(str(currNode1Entry.get()))
            self.file.write("\n")       
            self.file.write(str(currNode2Entry.get()))
            self.file.write("\n")

        otherFrame.onClose()

    def addACCurrentSource(self):
        self.hide()
        otherFrame = addFrame(self, "Add AC Current Source")
        currValueEntry = Tk.Entry(otherFrame)
        currFreqEntry = Tk.Entry(otherFrame)        
        currShiftEntry = Tk.Entry(otherFrame)
        currNode1Entry = Tk.Entry(otherFrame)        
        currNode2Entry = Tk.Entry(otherFrame)

        rowStart = otherFrame.getRow()
        currValueLabel=Tk.Label(otherFrame,text = "Nilai Amplitude Sumber Arus(A):")        
        currFreqLabel=Tk.Label(otherFrame,text = "Nilai Frekuensi (Hz):")        
        currShiftLabel=Tk.Label(otherFrame,text = "Nilai Shifting (s):")        
        currNode1Label=Tk.Label(otherFrame,text = "Nilai Node 1:")
        currNode2Label=Tk.Label(otherFrame,text = "Nilai Node 2:")
        currTypeLabel = Tk.Label(otherFrame,text = "Tipe Sumber:")

        currValueLabel.grid(row=rowStart, column = 1)
        currFreqLabel.grid(row = rowStart + 1, column = 1)
        currNode1Label.grid(row=rowStart+2, column = 1)        
        currNode2Label.grid(row=rowStart+3, column = 1)
        currShiftLabel.grid(row=rowStart+4, column = 1)

        currValueEntry.grid(row=rowStart, column = 2)
        currFreqEntry.grid(row = rowStart + 1, column = 2)
        currNode1Entry.grid(row=rowStart+2, column = 2)                
        currNode2Entry.grid(row=rowStart+3, column = 2)    
        currShiftEntry.grid(row=rowStart+4, column = 2)

        typeAC = ["Sinus", "Kotak"]
        variable = Tk.StringVar(otherFrame)
        variable.set(typeAC[0])

        typeOptions = apply(Tk.OptionMenu, (otherFrame, variable)+tuple(typeAC))
        currTypeLabel.grid(row = rowStart+5, column = 1)
        typeOptions.grid(row = rowStart+5, column = 2)

        submit = Tk.Button(otherFrame, text = "Submit")
        cancel = Tk.Button(otherFrame, text = "Cancel")
        submit['command'] =lambda binst =submit, binst2 = cancel:self.putCurrentACSource(1,currValueEntry,
                currNode1Entry,currNode2Entry,currFreqEntry,currShiftEntry,variable, binst, binst2, otherFrame)
        cancel['command'] = lambda binst = cancel,  binst2 = submit:self.putCurrentACSource(0,currValueEntry,
                currNode1Entry,currNode2Entry,currFreqEntry,currShiftEntry,variable, binst, binst2, otherFrame)
        cancel.grid(row = rowStart+6, column = 4)
        submit.grid(row = rowStart+6, column = 5)

    def putCurrentACSource(self, state ,currValueEntry, currNode1Entry, currNode2Entry,currFreqEntry,currShiftEntry,variable,binst, binst2, otherFrame):
        if state == 1:
            self.file.write("J")
            self.file.write("\n")
            self.file.write(str(currValueEntry.get())) #A
            self.file.write("\n")       
            self.file.write(str(currNode1Entry.get()))#1
            self.file.write("\n")       
            self.file.write(str(currNode2Entry.get()))#2
            self.file.write("\n")
            if variable.get() == "Sinus":
                self.file.write("0")
            else :
                self.file.write("1")
            self.file.write("\n")
            self.file.write(str(currFreqEntry.get()))
            self.file.write("\n")            
            self.file.write(str(currShiftEntry.get()))
            self.file.write("\n")            

        otherFrame.onClose()        

    def addACVoltageSource(self):
        self.hide()
        otherFrame = addFrame(self, "Add AC Voltage Source")
        voltValueEntry = Tk.Entry(otherFrame)
        voltFreqEntry = Tk.Entry(otherFrame)        
        voltShiftEntry = Tk.Entry(otherFrame)
        voltNode1Entry = Tk.Entry(otherFrame)        
        voltNode2Entry = Tk.Entry(otherFrame)

        rowStart = otherFrame.getRow()
        voltValueLabel=Tk.Label(otherFrame,text = "Nilai Amplitude Sumber Tegangan(V):")        
        voltFreqLabel=Tk.Label(otherFrame,text = "Nilai Frekuensi (Hz):")        
        voltShiftLabel=Tk.Label(otherFrame,text = "Nilai Shifting (s):")        
        voltNode1Label=Tk.Label(otherFrame,text = "Nilai Node 1:")
        voltNode2Label=Tk.Label(otherFrame,text = "Nilai Node 2:")
        voltTypeLabel = Tk.Label(otherFrame,text = "Tipe Sumber:")

        voltValueLabel.grid(row=rowStart, column = 1)
        voltFreqLabel.grid(row = rowStart + 1, column = 1)
        voltNode1Label.grid(row=rowStart+2, column = 1)        
        voltNode2Label.grid(row=rowStart+3, column = 1)
        voltShiftLabel.grid(row=rowStart+4, column = 1)

        voltValueEntry.grid(row=rowStart, column = 2)
        voltFreqEntry.grid(row = rowStart + 1, column = 2)
        voltNode1Entry.grid(row=rowStart+2, column = 2)                
        voltNode2Entry.grid(row=rowStart+3, column = 2)    
        voltShiftEntry.grid(row=rowStart+4, column = 2)

        typeAC = ["Sinus", "Kotak"]
        variable = Tk.StringVar(otherFrame)
        variable.set(typeAC[0])

        typeOptions = apply(Tk.OptionMenu, (otherFrame, variable)+tuple(typeAC))
        voltTypeLabel.grid(row = rowStart+5, column = 1)
        typeOptions.grid(row = rowStart+5, column = 2)

        submit = Tk.Button(otherFrame, text = "Submit")
        cancel = Tk.Button(otherFrame, text = "Cancel")
        submit['command'] =lambda binst =submit, binst2 = cancel:self.putACtVoltageSource(1,voltValueEntry,
                voltNode1Entry,voltNode2Entry,voltFreqEntry,voltShiftEntry,variable, binst, binst2, otherFrame)
        cancel['command'] = lambda binst = cancel,  binst2 = submit:self.putACVoltageSource(0,voltValueEntry,
                voltNode1Entry,voltNode2Entry,voltFreqEntry,voltShiftEntry,variable, binst, binst2, otherFrame)
        cancel.grid(row = rowStart+6, column = 4)
        submit.grid(row = rowStart+6, column = 5)

    def putACVoltageSource(self, state ,voltValueEntry, voltNode1Entry, voltNode2Entry,voltFreqEntry,voltShiftEntry,variable,binst, binst2, otherFrame):
        if state == 1:
            self.file.write("W")
            self.file.write("\n")
            self.file.write(str(voltValueEntry.get())) #A
            self.file.write("\n")       
            self.file.write(str(voltNode1Entry.get()))#1
            self.file.write("\n")       
            self.file.write(str(voltNode2Entry.get()))#2
            self.file.write("\n")
            if variable.get() == "Sinus":
                self.file.write("0")
            else :
                self.file.write("1")
            self.file.write("\n")
            self.file.write(str(voltFreqEntry.get()))
            self.file.write("\n")            
            self.file.write(str(voltShiftEntry.get()))
            self.file.write("\n")            

        otherFrame.onClose()        


class simulateCircuit(Tk.Toplevel):
    def __init__(self, original):
        """Constructor"""
        self.original_frame = original
        self.file = open("timefile.txt", "w+")
        self.nameFileGrou = "groundfile.txt"
        
        self.Row = 0
        Tk.Toplevel.__init__(self)
        self.geometry("800x600")
        self.title("Simulate")

        rowStart = self.getRow()

        timestartEntry = Tk.Entry(self)        
        timendEntry = Tk.Entry(self)

        timestartLabel=Tk.Label(self,text = "Time Start :")        
        timeendLabel=Tk.Label(self,text = "Time End :")

        timestartEntry.grid(row = rowStart, column = 2)
        timestartLabel.grid(row=rowStart, column = 1)

        timendEntry.grid(row=rowStart+1, column = 2)
        timeendLabel.grid(row=rowStart+1, column = 1)
        self.getRow()

        btngrou =Tk.Button(self, text="Add Ground", command=self.addGround)
        btngrou.grid(row = self.getRow())                

        btnsim = Tk.Button(self, text="Simulate", command=lambda :self.simulateNow(timestartEntry, timendEntry))
        btnsim.grid(row = self.getRow(), column = 0)        

        btn = Tk.Button(self, text="Close", command=self.onClose)
        btn.grid(row = self.getRow(), column = 0) 
    def hide(self):
        self.withdraw()

    def show(self):
        """"""
        self.update()
        self.deiconify()           

    def onClose(self):
        self.destroy()
        self.original_frame.show()    

    def getRow(self):
        self.Row += 1
        return (self.Row - 1)

    def addGround(self):
        self.hide()
        otherFrame = addFrame(self, "Add Ground")
        grouValueEntry = Tk.Entry(otherFrame)

        rowStart = otherFrame.getRow()
        grouValueLabel=Tk.Label(otherFrame,text = "Nilai Node Ground:")        

        grouValueLabel.grid(row=rowStart, column = 1)

        grouValueEntry.grid(row=rowStart, column = 2)

        submit = Tk.Button(otherFrame, text = "Submit")
        cancel = Tk.Button(otherFrame, text = "Cancel")
        submit['command'] =lambda binst =submit, binst2 = cancel:self.putGround(1,grouValueEntry,otherFrame)
        cancel['command'] = lambda binst = cancel,  binst2 = submit:self.putGround(0,grouValueEntry,otherFrame)
        cancel.grid(row = rowStart+4, column = 4)
        submit.grid(row = rowStart+4, column = 5)

    def putGround(self, state ,grouValueEntry, otherFrame):
        if state == 1:
            self.fileGrou = open("groundfile.txt", "w+")
            self.fileGrou.write(str(grouValueEntry.get()))
            self.file.write("\n")
            self.fileGrou.close()

        otherFrame.onClose() 

    def simulateNow(self, timestartEntry, timendEntry):
        self.file.write(str(timestartEntry.get()))
        self.file.write("\n")
        self.file.write(str(timendEntry.get()))
        self.file.write("\n")

        os.chdir(cwd)

        process = subprocess.Popen(self.original_frame.simulateexecute, stdout=subprocess.PIPE)
        print process
        process.wait()

        self.onClose()

    
class MyApp(object):
    
    def __init__(self, parent, fileName):
        self.inputFile = "infile.txt"
        self.simulateexecute = "program.exe"
        self.root = parent
        self.root.title("Main Menu")
        self.frame = Tk.Frame(parent)
        self.frame.pack()        

        self.data = pd.read_csv(fileName)
        self.Dict = {}        
        
                
        btnShowPlot = Tk.Button(self.frame, text="Show Plot", command=self.showPlotInMenu)
        btnShowPlot.pack()
        btnAddComponent = Tk.Button(self.frame, text="Add Component", command=self.insertComponentInMenu)
        btnAddComponent.pack()
        btnSimulate = Tk.Button(self.frame, text="Simulate", command=self.simulate)
        btnSimulate.pack()

        btnExit = Tk.Button(self.frame, text = "Exit", command = self.root.destroy)
        btnExit.pack()
        
    def hide(self):
        self.root.withdraw()        
        
    def onCloseOtherFrame(self, otherFrame):
        otherFrame.destroy()
        self.show()
        
    def show(self):
        """"""
        self.root.update()
        self.root.deiconify()
    #showplot
    def showPlotInMenu(self):
        self.hide()
        subFrame = showPlotClass(self)   

    # insert component
    def insertComponentInMenu(self)   :
        self.hide()
        subFrame = addComponentClass(self, self.inputFile)

    def simulate(self):
        self.hide()
        subFrame = simulateCircuit(self)


     
    
#----------------------------------------------------------------------
if __name__ == "__main__":
    root = Tk.Tk()
    root.geometry("800x600")
    app = MyApp(root, "outfile.csv")
    root.mainloop()