# Michael Jermyn - final

import sys
import vtk
from vtk.qt4.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor
from PyQt4.QtCore import *
from PyQt4.QtGui import *

# TOP LEFT
class VTK_Widget1(QWidget):
    
    def __init__(self, parent=None):

        super(VTK_Widget1, self).__init__(parent)
        self.source_is_connected = False
        self.source2_is_connected = False
        
        self.axis = 0; # 0 is z, 1 is y, 2 is x
        
        # vtk to point data
        self.c2p = vtk.vtkCellDataToPointData()
        self.opacityTransferFunction = vtk.vtkPiecewiseFunction()
        self.colorTransferFunction = vtk.vtkColorTransferFunction()

        # create a volume property for describing how the data will look
        self.volumeProperty = vtk.vtkVolumeProperty()
        self.volumeProperty.SetColor(self.colorTransferFunction)
        self.volumeProperty.SetScalarOpacity(self.opacityTransferFunction)
        self.volumeProperty.ShadeOn()
        self.volumeProperty.SetInterpolationTypeToLinear()

        # create a ray cast mapper
        self.compositeFunction = vtk.vtkUnstructuredGridBunykRayCastFunction()
        self.volumeMapper = vtk.vtkUnstructuredGridVolumeRayCastMapper()
        self.volumeMapper.SetRayCastFunction(self.compositeFunction)
        self.volumeMapper.SetInputConnection(self.c2p.GetOutputPort())
        
        # create a volume
        self.volume = vtk.vtkVolume()
        self.volume.SetMapper(self.volumeMapper)
        self.volume.SetProperty(self.volumeProperty)
        self.volume.VisibilityOff()
        
        # cutters
        self.cutPlane = vtk.vtkPlane()
        self.cutPlane.SetNormal(0, 0, 1)
        
        self.cutter = vtk.vtkCutter()
        self.cutter.SetCutFunction(self.cutPlane)
        self.cutter.SetValue(0,0)

        self.cutterMapper=vtk.vtkPolyDataMapper()
        self.cutterMapper.SetInputConnection(self.cutter.GetOutputPort())

        self.cutterActor=vtk.vtkActor()
        self.cutterActor.SetMapper(self.cutterMapper)
        self.cutterActor.VisibilityOff()
        
        # create the VTK widget for rendering
        self.vtkw=QVTKRenderWindowInteractor(self)
        self.ren = vtk.vtkRenderer()
        self.vtkw.GetRenderWindow().AddRenderer(self.ren)
        self.ren.AddVolume(self.volume)
        self.ren.AddActor(self.cutterActor)
        
        # we want now to have a slider for setting the cut plane
        self.cutPlaneSlider = QSlider(Qt.Vertical)
        self.cutPlaneSlider.setValue(50)
        self.cutPlaneSlider.setRange(0,100)
        self.cutPlaneSlider.setTickPosition(QSlider.NoTicks) 
        self.connect(self.cutPlaneSlider,SIGNAL("valueChanged(int)"),self.AdjustCutPlane)
        
        # layout manager
        self.layout = QHBoxLayout()
        self.layout.addWidget(self.vtkw)
        self.layout.addSpacing(10)
        self.layout.addWidget(self.cutPlaneSlider)
        self.setLayout(self.layout)
        
        # initialize the interactor
        self.vtkw.Initialize()
        self.vtkw.Start()
        
        
    def SetSource(self,source):   

        self.source = source
        self.c2p.SetInput(self.source)
        self.volume.VisibilityOn()
        
        # the volume will be made completely transparent for values below 5%, 
        # somewhat transparent up to 65% of the scalar range
        # and to opaque for values between 65% and 100%
        range = source.GetScalarRange()
        zero_pc  = range[0]
        fifty_pc  = range[0]+(range[1]-range[0])*0.50
        hundred_pc  = range[1]
             
        self.opacityTransferFunction.AddPoint(zero_pc, 0.01)
        self.opacityTransferFunction.AddPoint(fifty_pc, 0.01)
        self.opacityTransferFunction.AddPoint(fifty_pc+1e-6, 0.2) # anything > 65% is transparent    
        
        self.colorTransferFunction.AddRGBPoint(zero_pc, 0.0, 0.0, 1.0)
        self.colorTransferFunction.AddRGBPoint(fifty_pc, 1.0, 0.5, 0.0)
        self.colorTransferFunction.AddRGBPoint(hundred_pc, 1.0, 0.0, 0.0)
         
        self.ren.ResetCamera() 
        self.vtkw.GetRenderWindow().Render()
        self.source_is_connected = True
        
    def SetSource2(self,source):   

        self.source2 = source.GetOutput()

        self.cutter.SetInput(self.source2)
        
        center = self.source2.GetCenter()
        self.cutPlane.SetOrigin(center)
        self.cutPlaneSlider.setValue(50)
        self.cutterActor.VisibilityOn()
        self.cutterMapper.SetScalarRange(self.source2.GetScalarRange())
        
        self.ren.ResetCamera() 
        self.vtkw.GetRenderWindow().Render()
        self.source2_is_connected = True
        
    def AdjustCutPlane(self):
        
        if self.source2_is_connected:
        
            slider_pos = self.cutPlaneSlider.value()
            center = self.source2.GetCenter()
            bounds = self.source2.GetBounds() 
            
            if self.axis == 0:
                cut_z_pos = bounds[4]+(bounds[5]-bounds[4])*(slider_pos/100.0)
                self.cutPlane.SetOrigin(center[0],center[1],cut_z_pos)
                self.cutPlane.SetNormal(0, 0, 1)
            elif self.axis == 1:
                cut_y_pos = bounds[2]+(bounds[3]-bounds[2])*(slider_pos/100.0)
                self.cutPlane.SetOrigin(center[0],cut_y_pos,center[2]) 
                self.cutPlane.SetNormal(0, 1, 0)
            elif self.axis == 2:
                cut_x_pos = bounds[0]+(bounds[1]-bounds[0])*(slider_pos/100.0)
                self.cutPlane.SetOrigin(cut_x_pos,center[1],center[2]) 
                self.cutPlane.SetNormal(1, 0, 0)
                
            self.vtkw.GetRenderWindow().Render() 
            
    def ToggleCutPlane(self):
        
        if self.sender().checkState() == Qt.Checked:
            if self.source2_is_connected:
                self.cutterActor.VisibilityOn()
                self.vtkw.GetRenderWindow().Render()
        else:
            self.cutterActor.VisibilityOff()
            self.vtkw.GetRenderWindow().Render()
            
    def SetAxis(self):
        
        self.axis = self.sender().currentIndex()
        self.AdjustCutPlane()
        
    def SetProperty(self,property):
        
        if self.source_is_connected:
            self.source.GetPointData().SetActiveScalars(property)
            self.vtkw.GetRenderWindow().Render()
        
        
# TOP RIGHT
class VTK_Widget2(QWidget):
    def __init__(self, parent=None):
        
        super(VTK_Widget2, self).__init__(parent)

        self.source_is_connected = False
             
        self.cutPlane = vtk.vtkPlane()
        self.cutPlane.SetNormal(0, 0, 1) # x-y plane in this view
        
        self.cutter = vtk.vtkCutter()
        self.cutter.SetCutFunction(self.cutPlane)
        self.cutter.SetValue(0,0)

        self.cutterMapper=vtk.vtkPolyDataMapper()
        self.cutterMapper.SetInputConnection(self.cutter.GetOutputPort())
        
        # colorbar
        self.colorbar = vtk.vtkScalarBarActor()
        self.colorbar.SetLookupTable(self.cutterMapper.GetLookupTable())
        self.colorbar.SetNumberOfLabels(4)
        self.colorbar.SetMaximumWidthInPixels(50)
        
        # lookup table
        self.lookupTable = vtk.vtkLookupTable()
        self.lookupTable.Build()
        k = 40.0
        self.lookupTable.SetNumberOfTableValues(k*3)
        for i in range(0,int(k-1)):
            self.lookupTable.SetTableValue(i,(i/k,0.0,0.0,1.0))
            self.lookupTable.SetTableValue(i+int(k),(1.0,i/k,0.0,1.0))
            self.lookupTable.SetTableValue(i+2*int(k),(1.0,1.0,i/k,1.0))
        self.lookupTable.SetTableValue(119,(1.0,1.0,1.0,1.0))
        self.cutterMapper.SetLookupTable(self.lookupTable)
        self.colorbar.SetLookupTable(self.lookupTable)

        self.cutterActor=vtk.vtkActor()
        self.cutterActor.SetMapper(self.cutterMapper)
                
        self.vtkw=QVTKRenderWindowInteractor(self)
      
        self.ren = vtk.vtkRenderer()
        self.vtkw.GetRenderWindow().AddRenderer(self.ren)
        self.ren.AddActor(self.cutterActor)
        self.ren.AddActor2D(self.colorbar)
        
        self.cutterActor.VisibilityOff()
           
        # we want now to have a slider for setting the cut plane
        self.cutPlaneSlider = QSlider(Qt.Vertical)
        self.cutPlaneSlider.setValue(50) # the range is set 0 to 100, and the initial in middle
        self.cutPlaneSlider.setRange(0,100)
        self.cutPlaneSlider.setTickPosition(QSlider.NoTicks) 
        # connect now the slider "valueChanged(int)" signal to the method 
        # AdjustCutPlane of this class (if we need we can connect also to methods of other classes)
        self.connect(self.cutPlaneSlider,SIGNAL("valueChanged(int)"),self.AdjustCutPlane)
        
        self.layout = QHBoxLayout()
        self.layout.addWidget(self.vtkw)
        self.layout.addSpacing(10)
        self.layout.addWidget(self.cutPlaneSlider)
        self.setLayout(self.layout)
        
        self.vtkw.Initialize()
        self.vtkw.Start()
        self.vtkw.Disable()
        
    def SetSource(self,source):
        
        self.source=source
        self.cutter.SetInput(self.source)

        center = self.source.GetCenter()
        self.cutPlane.SetOrigin(center)
        self.cutPlaneSlider.setValue(50)
        self.cutterActor.VisibilityOn()
        self.cutterMapper.SetScalarRange(source.GetScalarRange())
        
        self.ren.GetActiveCamera().ParallelProjectionOn()
        self.ren.GetActiveCamera().SetViewUp(0,1,0)
        self.ren.GetActiveCamera().SetFocalPoint(center)
        
        x=center[0]
        y=center[1]
        z=center[2]+1
        self.ren.GetActiveCamera().SetPosition(x,y,z)
        self.ren.ResetCamera()
        cam_pos = self.ren.GetActiveCamera().GetPosition()
        bounds = self.source.GetBounds() 
        clip_near = cam_pos[2]-bounds[5]
        clip_far = cam_pos[2]-bounds[4]
        self.ren.GetActiveCamera().SetClippingRange(clip_near,clip_far)
          
        self.vtkw.GetRenderWindow().Render()
        self.source_is_connected = True
        
    def AdjustCutPlane(self):
        
        if self.source_is_connected:
        
            slider_pos = self.sender().value() 
            center = self.source.GetCenter()
            bounds = self.source.GetBounds() 
            
            cut_z_pos = bounds[4]+(bounds[5]-bounds[4])*(slider_pos/100.0) 
            self.cutPlane.SetOrigin(center[0],center[1],cut_z_pos)

            self.vtkw.GetRenderWindow().Render() 
        
# BOTTOM LEFT        
class  VTK_Widget3(QWidget):
    def __init__(self, parent=None):

        super(VTK_Widget3, self).__init__(parent)
        
        self.source_is_connected = False 
        
        self.cutPlane = vtk.vtkPlane()
        self.cutPlane.SetNormal(1, 0, 0) # y-z plane in this window
        
        self.cutter = vtk.vtkCutter()
        self.cutter.SetCutFunction(self.cutPlane)
        self.cutter.SetValue(0,0)

        self.cutterMapper=vtk.vtkPolyDataMapper()
        self.cutterMapper.SetInputConnection(self.cutter.GetOutputPort())

        # colorbar
        self.colorbar = vtk.vtkScalarBarActor()
        self.colorbar.SetLookupTable(self.cutterMapper.GetLookupTable())
        self.colorbar.SetNumberOfLabels(4)
        self.colorbar.SetMaximumWidthInPixels(50)
        
        # lookup table
        self.lookupTable = vtk.vtkLookupTable()
        self.lookupTable.Build()
        k = 40.0
        self.lookupTable.SetNumberOfTableValues(k*3)
        for i in range(0,int(k-1)):
            self.lookupTable.SetTableValue(i,(i/k,0.0,0.0,1.0))
            self.lookupTable.SetTableValue(i+int(k),(1.0,i/k,0.0,1.0))
            self.lookupTable.SetTableValue(i+2*int(k),(1.0,1.0,i/k,1.0))
        self.lookupTable.SetTableValue(119,(1.0,1.0,1.0,1.0))
        self.cutterMapper.SetLookupTable(self.lookupTable)
        self.colorbar.SetLookupTable(self.lookupTable)

        self.cutterActor=vtk.vtkActor()
        self.cutterActor.SetMapper(self.cutterMapper)
                
        self.vtkw=QVTKRenderWindowInteractor(self)
      
        self.ren = vtk.vtkRenderer()
        self.vtkw.GetRenderWindow().AddRenderer(self.ren)
        self.ren.AddActor(self.cutterActor)
        self.ren.AddActor2D(self.colorbar)
        
        self.cutterActor.VisibilityOff()
   
        self.cutPlaneSlider = QSlider(Qt.Vertical)
        self.cutPlaneSlider.setValue(50)
        self.cutPlaneSlider.setRange(0,100)
        self.cutPlaneSlider.setTickInterval(5)
        self.cutPlaneSlider.setTickPosition(QSlider.NoTicks)
        
        self.connect(self.cutPlaneSlider,SIGNAL("valueChanged(int)"),self.AdjustCutPlane)
        
        self.layout = QHBoxLayout()
        self.layout.addWidget(self.vtkw)
        self.layout.addSpacing(10)
        self.layout.addWidget(self.cutPlaneSlider)
        self.setLayout(self.layout)
        
        self.vtkw.Initialize()
        self.vtkw.Start()
        self.vtkw.Disable() 
        
    def SetSource(self,source):
        
        self.source=source
        
        self.cutter.SetInput(self.source)
        
        center = self.source.GetCenter()
        self.cutPlane.SetOrigin(center)
        self.cutPlaneSlider.setValue(50)
        self.cutterActor.VisibilityOn()
        
        self.cutterMapper.SetScalarRange(source.GetScalarRange())
        
        self.ren.GetActiveCamera().ParallelProjectionOn()
        self.ren.GetActiveCamera().SetViewUp(0,0,1)
        self.ren.GetActiveCamera().SetFocalPoint(center)

        x=center[0]+1
        y=center[1]
        z=center[2]
        self.ren.GetActiveCamera().SetPosition(x,y,z)
        self.ren.ResetCamera()
          
        cam_pos = self.ren.GetActiveCamera().GetPosition()
        bounds = self.source.GetBounds() 
        clip_near = cam_pos[0]-bounds[1]
        clip_far = cam_pos[0]-bounds[0]
        self.ren.GetActiveCamera().SetClippingRange(clip_near,clip_far)
        
        self.vtkw.GetRenderWindow().Render()
        self.source_is_connected = True
        
    def AdjustCutPlane(self):
        
        if self.source_is_connected:
        
            slider_pos = self.sender().value()
            
            center = self.source.GetCenter()
            bounds = self.source.GetBounds() 
            
            cut_x_pos = bounds[0]+(bounds[1]-bounds[0])*(slider_pos/100.0) 
                    
            self.cutPlane.SetOrigin(cut_x_pos,center[1],center[2])    
            self.vtkw.GetRenderWindow().Render() 
        
# BOTTOM RIGHT   
class VTK_Widget4(QWidget):
    def __init__(self, parent=None):

        super(VTK_Widget4, self).__init__(parent)
        
        self.source_is_connected = False
  
        self.cutPlane = vtk.vtkPlane()
        self.cutPlane.SetNormal(0, 1, 0) # x-z plane in this window
        
        self.cutter = vtk.vtkCutter()
        self.cutter.SetCutFunction(self.cutPlane)
        self.cutter.SetValue(0,0)

        self.cutterMapper=vtk.vtkPolyDataMapper()
        self.cutterMapper.SetInputConnection(self.cutter.GetOutputPort())

        # colorbar
        self.colorbar = vtk.vtkScalarBarActor()
        self.colorbar.SetLookupTable(self.cutterMapper.GetLookupTable())
        self.colorbar.SetNumberOfLabels(4)
        self.colorbar.SetMaximumWidthInPixels(50)
        
        # lookup table
        self.lookupTable = vtk.vtkLookupTable()
        self.lookupTable.Build()
        k = 40.0
        self.lookupTable.SetNumberOfTableValues(k*3)
        for i in range(0,int(k-1)):
            self.lookupTable.SetTableValue(i,(i/k,0.0,0.0,1.0))
            self.lookupTable.SetTableValue(i+int(k),(1.0,i/k,0.0,1.0))
            self.lookupTable.SetTableValue(i+2*int(k),(1.0,1.0,i/k,1.0))
        self.lookupTable.SetTableValue(119,(1.0,1.0,1.0,1.0))
        self.cutterMapper.SetLookupTable(self.lookupTable)
        self.colorbar.SetLookupTable(self.lookupTable)

        self.cutterActor=vtk.vtkActor()
        self.cutterActor.SetMapper(self.cutterMapper)
                
        self.vtkw=QVTKRenderWindowInteractor(self)
      
        self.ren = vtk.vtkRenderer()
        self.vtkw.GetRenderWindow().AddRenderer(self.ren)
        self.ren.AddActor(self.cutterActor)
        self.ren.AddActor2D(self.colorbar)
        #self.ren.SetBackground(1.0, 1.0, 1.0)
        
        self.cutterActor.VisibilityOff()
        
        self.cutPlaneSlider = QSlider(Qt.Vertical)
        self.cutPlaneSlider.setValue(50)
        self.cutPlaneSlider.setRange(0,100)
        self.cutPlaneSlider.setTickInterval(5)
        self.cutPlaneSlider.setTickPosition(QSlider.NoTicks)
   
        self.connect(self.cutPlaneSlider,SIGNAL("valueChanged(int)"),self.AdjustCutPlane)
        
        self.layout = QHBoxLayout()
        self.layout.addWidget(self.vtkw)
        self.layout.addSpacing(10)
        self.layout.addWidget(self.cutPlaneSlider)
        self.setLayout(self.layout)
        
        self.vtkw.Initialize()
        self.vtkw.Start()
        self.vtkw.Disable() 
        
        
    def SetSource(self,source):
        
        self.source=source
        
        self.cutter.SetInput(self.source)
        
        center = self.source.GetCenter()
        self.cutPlane.SetOrigin(center)
        self.cutPlaneSlider.setValue(50)
        self.cutterActor.VisibilityOn()
        
        self.cutterMapper.SetScalarRange(source.GetScalarRange())
        
        self.ren.GetActiveCamera().ParallelProjectionOn()
        self.ren.GetActiveCamera().SetViewUp(0,0,1)
        self.ren.GetActiveCamera().SetFocalPoint(center)
     
        x=center[0]
        y=center[1]+1 
        z=center[2] 
        self.ren.GetActiveCamera().SetPosition(x,y,z)
        self.ren.ResetCamera()
        
        cam_pos = self.ren.GetActiveCamera().GetPosition()
        bounds = self.source.GetBounds() 
        clip_near = cam_pos[1]-bounds[3]
        clip_far = cam_pos[1]-bounds[2]
        self.ren.GetActiveCamera().SetClippingRange(clip_near,clip_far)
        
        self.vtkw.GetRenderWindow().Render()
        self.source_is_connected = True
        
    def AdjustCutPlane(self):
        
        if self.source_is_connected:
        
            slider_pos = self.sender().value()
            center = self.source.GetCenter()
            bounds = self.source.GetBounds() 
            cut_y_pos = bounds[2]+(bounds[3]-bounds[2])*(slider_pos/100.0) 
            self.cutPlane.SetOrigin(center[0],cut_y_pos,center[2])

            self.vtkw.GetRenderWindow().Render() 
        
# MAIN WINDOW                          
class MainVizWindow(QMainWindow):    
    def __init__(self, parent=None):
         
         QMainWindow.__init__(self, parent)
           
         self.setWindowTitle(self.tr("Nirfast"))
         
         # splitters are used for generating the four views
         self.VSplitter = QSplitter(Qt.Vertical)
         self.HSplitterTop = QSplitter(Qt.Horizontal)
         self.HSplitterBottom = QSplitter(Qt.Horizontal)
         
         # one instance of each of the VTK_Widget classes
         self.vtk_widget_1 = VTK_Widget1(self)
         self.vtk_widget_2 = VTK_Widget2(self)
         self.vtk_widget_3 = VTK_Widget3(self)
         self.vtk_widget_4 = VTK_Widget4(self)
         
         # the VTK widgets are added to the splitters
         self.VSplitter.addWidget(self.HSplitterTop)
         self.VSplitter.addWidget(self.HSplitterBottom)
         self.HSplitterTop.addWidget(self.vtk_widget_1)
         self.HSplitterTop.addWidget(self.vtk_widget_2)
         self.HSplitterBottom.addWidget(self.vtk_widget_3)
         self.HSplitterBottom.addWidget(self.vtk_widget_4)
         
         # the top splitter (vertical) is set as central widget
         self.setCentralWidget(self.VSplitter)
         
         # we embed a reader in the main window, which will fan out the data to all VTK views
         self.reader = vtk.vtkUnstructuredGridReader()
         self.reader2 = vtk.vtkDICOMImageReader()
         self.reader.SetFileName('')
         self.reader2.SetDirectoryName('')
         
         # we declare a file open action
         self.fileOpenAction = QAction("&Open Solution",self)
         self.fileOpenAction.setShortcut("Ctrl+O")
         self.fileOpenAction.setToolTip("Opens a VTK volume file")
         self.fileOpenAction.setStatusTip("Opens a VTK volume file")
         
         self.fileOpenAction2 = QAction("&Open DICOM",self)
         self.fileOpenAction2.setShortcut("Ctrl+D")
         self.fileOpenAction2.setToolTip("Opens a set of DICOMs")
         self.fileOpenAction2.setStatusTip("Opens a set of DICOMs")
     
         self.connect(self.fileOpenAction, SIGNAL("triggered()"),self.fileOpen)
         self.connect(self.fileOpenAction2, SIGNAL("triggered()"),self.fileOpen2)
         
         self.fileMenu = self.menuBar().addMenu("&File")
         self.fileMenu.addAction(self.fileOpenAction)   
         self.fileMenu.addAction(self.fileOpenAction2)
         
         # property label
         self.label_property = QLabel("Property: ")
         
         # property dropdown
         self.dropdown_property = QComboBox()
         
         # spacing label
         self.label_spacing = QLabel("      ")
         
         # dicom slice checkbox
         self.cbox_cutplane=QCheckBox("DICOM Slice")
         self.cbox_cutplane.setCheckState(Qt.Unchecked) 
         
         # dicom slice axis dropdown
         self.dropdown_axis = QComboBox()
         self.dropdown_axis.addItem('z axis')
         self.dropdown_axis.addItem('y axis')
         self.dropdown_axis.addItem('x axis')
         
         # toolbar
         self.viewToolbar = self.addToolBar("View")
         self.viewToolbar.setObjectName("ViewToolbar")
         self.viewToolbar.addWidget(self.label_property)
         self.viewToolbar.addWidget(self.dropdown_property)
         self.viewToolbar.addWidget(self.label_spacing)
         self.viewToolbar.addWidget(self.cbox_cutplane)
         self.viewToolbar.addWidget(self.dropdown_axis)
		 
         self.connect(self.cbox_cutplane, SIGNAL("stateChanged(int)"), self.vtk_widget_1.ToggleCutPlane) 
         self.connect(self.dropdown_axis, SIGNAL("currentIndexChanged(int)"), self.vtk_widget_1.SetAxis)
         self.connect(self.dropdown_property, SIGNAL("currentIndexChanged(int)"), self.SetProperty)

            
    def setSource(self,source):
        
        self.reader.SetFileName(source)
        self.reader.ReadAllScalarsOn()
        self.reader.Update()    
        pointdata = self.reader.GetOutput().GetPointData()
        for i in range(pointdata.GetNumberOfArrays()):      
            self.dropdown_property.addItem(pointdata.GetArrayName(i))
        mainwindow.vtk_widget_1.SetSource(self.reader.GetOutput())
        mainwindow.vtk_widget_2.SetSource(self.reader.GetOutput())
        mainwindow.vtk_widget_3.SetSource(self.reader.GetOutput())
        mainwindow.vtk_widget_4.SetSource(self.reader.GetOutput())
           
    # define a method for handling the file open action     
    def fileOpen(self):
        
        dir ="."
        format = "*.vtk"
        fname = unicode(QFileDialog.getOpenFileName(self,"Open VTK File",dir,format))
                        
        if (len(fname)>0):         
            self.setSource(fname)
    
    def fileOpen2(self):
        
        dir ="."
        fname = unicode(QFileDialog.getExistingDirectory(self,"Select DICOM Directory",dir))
                        
        if (len(fname)>0):         
            self.reader2.SetDirectoryName(fname)
            self.reader2.Update()          
            mainwindow.vtk_widget_1.SetSource2(self.reader2)
            self.cbox_cutplane.setCheckState(Qt.Checked) 
            
    def SetProperty(self):
        
        property = str(self.sender().currentText())
        pointdata = self.reader.GetOutput().GetPointData()
        pointdata.SetActiveScalars(property)
        self.vtk_widget_1.hide()
        self.vtk_widget_2.hide()
        self.vtk_widget_3.hide()
        self.vtk_widget_4.hide()
        self.vtk_widget_1 = VTK_Widget1(self)
        self.vtk_widget_2 = VTK_Widget2(self)
        self.vtk_widget_3 = VTK_Widget3(self)
        self.vtk_widget_4 = VTK_Widget4(self)
        self.HSplitterTop.addWidget(self.vtk_widget_1)
        self.HSplitterTop.addWidget(self.vtk_widget_2)
        self.HSplitterBottom.addWidget(self.vtk_widget_3)
        self.HSplitterBottom.addWidget(self.vtk_widget_4)
        mainwindow.vtk_widget_1.SetSource(self.reader.GetOutput())
        mainwindow.vtk_widget_2.SetSource(self.reader.GetOutput())
        mainwindow.vtk_widget_3.SetSource(self.reader.GetOutput())
        mainwindow.vtk_widget_4.SetSource(self.reader.GetOutput())
        
            
# START APPLICATION    
app = QApplication(sys.argv)
mainwindow = MainVizWindow()
mainwindow.show()
if sys.argv.__len__() > 1:
    source = sys.argv[1]
    mainwindow.setSource(source)
sys.exit(app.exec_())
