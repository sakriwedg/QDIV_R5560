# -*- coding: utf-8 -*-
"""
List of PSD detectors which can be selected in QDIV-gui

Created on Sat Jul 11 07:59:24 2026

@author: marchalj
"""
# Detector class
class detector:
    """
    Defines detector parameters

    Attributes
    ----------
    name : string
        detector name
    nb_of_tubes : int
        number of PSD tubes in the detector
    tube_length : float
        PSD tube length in mm
    tube_spacing : float
        spacing between PSD tubes in mm

    Methods
    -------
    showParameters()
        print detector parameters.
    """
    
    def __init__(self,name,nb_of_tubes,tube_length,tube_spacing):
  
        self.name=name
        self.nb_of_tubes=nb_of_tubes
        self.tube_length=tube_length
        self.tube_spacing=tube_spacing
        self.SN=[]
        self.orientation='Horizontal'
        self.flipud=False
        self.fliplr=False
        self.aspect_ratio=self.nb_of_tubes*self.tube_spacing/self.tube_length
        
    def showParameters(self):
        print("###########")
        print("Name: "+self.name)
        print("Number of tubes: "+str(self.nb_of_tubes))
        print("Tube length in mm: "+str(self.tube_length))
        print("Tube spacing in mm:"+str(self.tube_spacing))

# Build a list of detectors
list=[]
# CSPEC module
CSPEC_module=detector('CSPEC',32,3500.0,25.0)
CSPEC_module.SN=1
list.append(CSPEC_module)
# CSPEC detector
CSPEC_det=detector('CSPEC detector',32*12,3500.0,25.0)
list.append(CSPEC_det)

# Function showing the list of detectors available
def show():
    print("List of detectors available in the GUI")
    [detector.showParameters() for n,detector in enumerate(list)]

    

