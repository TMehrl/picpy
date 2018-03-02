#!/usr/bin/env python3
from __future__ import print_function
import h5py
from vtk import *
from ROOT import PData, PDataHiP, PGlobals
import sys, argparse
import numpy as np


# Command argument line parser
parser = argparse.ArgumentParser(description='3D renderer using VTK.')
parser.add_argument('sim', nargs='?', default='none', help='simulation name')
parser.add_argument('-b', action='store_true', help='run in batch mode without graphics')
parser.add_argument('-t', type=int, help='simulation time step')
#parser.add_argument('-i', type=int, dest='tstart', help='time start')
#parser.add_argument('-f', type=int, dest='tend', help='time end')
#parser.add_argument('-s', type=int, dest='tdelta', default=1,help='time delta (step between dumps)')
parser.add_argument('-z', type=float, dest='zoom', default=1,help='zoom')
parser.add_argument('-ar', type=float, dest='aratio', default=1,help='Aspect ratio trans/long')
parser.add_argument('-azi', type=float, dest='azimuth', default=0,help='azimuth')
parser.add_argument('-ele', type=float, dest='elevation', default=0,help='elevation')
parser.add_argument('--white', action='store_true', default=0,help='white background')
parser.add_argument('--surf', action='store_true', default=0,help='draw surfaces')
parser.add_argument('--log', action='store_true', default=0,help='log scale')
parser.add_argument('--axes', action='store_true', default=0,help='show axes')
parser.add_argument('--cbar', action='store_true', default=0,help='show color maps')
parser.add_argument('--transx', action='store_true', default=0,help='transpose transverse dimensions')
parser.add_argument('--lden', action='store_true', default=0,help='use local density')
parser.add_argument('--laser', action='store_true', default=0,help='show laser fields')
parser.add_argument('--txbck', action='store_true', default=0,help='textured background')
parser.add_argument('--test', action='store_true', default=0,help='run a direct example')
parser.add_argument('--nowin', action='store_true', default=0, help='no windows output (to run in batch)')

try:
    args = parser.parse_args()
except:
    parser.print_help()
    sys.exit(0)

if ('none' in args.sim) & (args.test == 0) :
    parser.print_help()
    sys.exit(0)

# End of command line setup

# Get data files
hfl = []
if args.test :
    hfl.append(h5py.File('data/charge-plasma-000026.h5','r'))
    hfl.append(h5py.File('data/charge-beam-driver-000026.h5','r'))
    hfl.append(h5py.File('data/charge-He-electrons-000026.h5','r'))
else :
    pData = PData(args.sim)
    if pData.isHiPACE() :
        del(pData)
        pData = PDataHiP(args.sim)

    pData.LoadFileNames(args.t)
    for i in range(0,pData.NSpecies()) :
        hfl.append(h5py.File(pData.GetChargeFileName(i).c_str(),'r'))
        print('File %i = %s' % (i,pData.GetChargeFileName(i).c_str()))

    if args.laser :
        hfl.append(h5py.File(pData.GetEfieldFileName(1).c_str(),'r'))

nfiles = len(hfl)
# -----------


# Get data, build multi-component vtkVolume,
# Set colors and opacities for the components...
data = []
npdata = []

# Multi component (single) volume
volumeprop = vtk.vtkVolumeProperty()
#volumeprop.SetIndependentComponents(ncomp)
volumeprop.IndependentComponentsOn()
volumeprop.SetInterpolationTypeToLinear()

# This block reads the colormaps file (if any)

fcolname = './%s/%s.col' % (pData.GetPath(),args.sim)
opafile = []
colfile = []
minfile = []
maxfile = []
stypefile = []

try:
    fcolmap = open(fcolname,'r')
except IOError:
    print('\nNo color map file found')
else:
    print('\nReading color maps file ',fcolname)

    for line in fcolmap: # Loop over the lines
        pos = line.find('Species kind : ')
        if pos>=0 :
            stypefile.append(line[pos+len('Species kind : '):].replace('\n',''))
            print('Type : ',stypefile[len(stypefile)-1])

        pos = line.find('# MinMax')
        if pos>=0 :
            line = fcolmap.next()
            minfile.append(float(line))
            line = fcolmap.next()
            maxfile.append(float(line))
            print('Min = %.2f  Max = %.2f' % (minfile[len(minfile)-1],maxfile[len(minfile)-1]))

        pos = line.find('# Opacities')
        if pos>=0 :
            chunk = []
            line = fcolmap.next()
            while line.find('#') < 0 :
                chunk.append([float(s) for s in line.split()])
                line = fcolmap.next()

            opafile.append(chunk)
            print('Opacities : ', opafile[len(opafile)-1])

        pos = line.find('# Colors')
        if pos>=0 :
            chunk = []
            line = fcolmap.next()
            while (line.find('#') < 0) & (len(line.rstrip())>0) :
                chunk.append([float(s) for s in line.split()])
                line = fcolmap.next()

            colfile.append(chunk)
            print('Colors : ', colfile[len(colfile)-1])

    fcolmap.close()

# end of reading color maps

# Loop over the files
opacity = []
color = []
Min = []
Max = []
stype = []
baseden = 1
localden = 1
for i, hf in enumerate(hfl):

    name = hf.attrs['NAME']
    data.append(hf.get(name[0]))
    data[i] = np.absolute(data[i])

    if pData.isHiPACE() :
        xmin = hf.attrs['XMIN']
        xmax = hf.attrs['XMAX']

        axisz = [xmin[0],xmax[0]]
        axisx = [xmin[1],xmax[1]]
        axisy = [xmin[2],xmax[2]]

        data[i] = np.transpose(data[i],(2, 1, 0))

    else:
        axisz = hf.get('AXIS/AXIS1')
        axisx = hf.get('AXIS/AXIS2')
        axisy = hf.get('AXIS/AXIS3')

    dz = (axisz[1]-axisz[0])/data[i].shape[2]
    dx = (axisx[1]-axisx[0])/data[i].shape[0]
    dy = (axisy[1]-axisy[0])/data[i].shape[1]

    uaxzmin = axisz[0]
    uaxzmax = axisz[1]
    uaxymin = axisy[0]
    uaxymax = axisy[1]
    uaxxmin = axisx[0]
    uaxxmax = axisx[1]

    udz = dz
    udx = dx
    udy = dy

    if args.aratio != 1 :
        udx = args.aratio * dx
        udy = args.aratio * dy
        uaxxmin = axisx[0] * args.aratio
        uaxxmax = uaxxmin + udx * data[i].shape[0]
        uaxymin = axisy[0] * args.aratio
        uaxymax = uaxymin + udy * data[i].shape[1]

    print('\nFilename : ',hf.filename)
    print('Axis z range: [%.2f,%.2f]  Nbins = %i  dz = %.4f' % (axisz[0],axisz[1],data[i].shape[2],dz) )
    print('Axis x range: [%.2f,%.2f]  Nbins = %i  dx = %.4f' % (axisx[0],axisx[1],data[i].shape[0],dx) )
    print('Axis y range: [%.2f,%.2f]  Nbins = %i  dy = %.4f' % (axisy[0],axisy[1],data[i].shape[1],dy) )

    maxvalue = np.amax(data[i])
    if maxvalue < 1E-5 : continue

    stype.append('default')

    npdata.append(np.array(data[i]))
    j = len(npdata) - 1

    if args.transx :
        npdata[i] = np.transpose(npdata[i],(1, 0, 2))

    Max.append(np.amax(npdata[j]))
    Min.append(0.01*Max[j])

    if args.test == 0 :
        if pData.GetDenMax(j) > 0 : Max[j] = pData.GetDenMax(j)
        if pData.GetDenMin(j) > 0 : Min[j] = pData.GetDenMin(j)
        else : Min[j] = 0.01 * Max[j]

    if len(opafile) == nfiles :
        Min[j] = minfile[j]
        Max[j] = maxfile[j]

    if Min[j] > Max[j] :
        Min[j] = 0.1 * Max[j]

    print('Maximum = %.2f' % maxvalue)
    print('Min = %.3e  Max = %.3e ' % (Min[j],Max[j]))

    if 'plasma' in hf.filename :
        # Min[j] = baseden
        if args.lden :
            localden = npdata[j][npdata[j].shape[0]-10][npdata[j].shape[1]-10][npdata[j].shape[2]-2]
            print('Local density = %f' % (localden))
            Min[j] *= localden/baseden
            Max[j] *= localden/baseden
    elif 'e2' in hf.filename:
        Min[j] = -Max[j]
    else :
        if Max[j] > maxvalue :
            Max[j] = maxvalue
            print('New max = %.3e' % Max[j])

    def lscale(x) :
#        x = x - Min[j] + 1
        x = np.where(x <= 1E-6,1E-6, x)
        x = np.log10(x)
        return x

    if args.log & ('e2' not in hf.filename) :
        npdata[j] = lscale(npdata[j])
        Min[j] = lscale(Min[j])
        Max[j] = lscale(Max[j])
        maxvalue = lscale(maxvalue)
        localden = lscale(localden)
        # stop = lscale(stop)

    # Intermediate points in the density scale
    def stop(x):
        x = Min[j] + (x/100) * (Max[j] - Min[j])
        return x

    # average value for non-zero data points
    # denavg = np.average(npdata[j],weights=npdata[j].astype(bool))
    # print('Average = %.2f' % denavg)

    # Opacity and color scales
    opacity.append(vtk.vtkPiecewiseFunction())
    color.append(vtk.vtkColorTransferFunction())

    if len(opafile) == nfiles :
        for k, opa in enumerate(opafile[j]) :
            opacity[j].AddPoint(stop(opa[0]),opa[1])
        for k, col in enumerate(colfile[j]) :
            color[j].AddRGBPoint(stop(col[0]),col[1],col[2],col[3])

        if stypefile[j] == 'plasma' :
            opacity[j].AddPoint( 0.98 * localden ,0.0)
            opacity[j].AddPoint( 1.02 * localden ,0.0)

    else :

        opalist = []
        colorlist = []

        if 'plasma' in hf.filename :
            stype[j] = 'plasma'

            opalist.append( [0,0.0] )
            opalist.append( [1,0.8] )
            opalist.append( [10,0.2] )
            opalist.append( [50,0.2] )
            opalist.append( [100,0.1] )

            colorlist.append( [0,0.7,0.7,0.7] )
            colorlist.append( [100,1.0,1.0,1.0] )

        elif 'beam' in hf.filename :
            stype[j] = 'beam'

            opalist.append( [0,0.0] )
            opalist.append( [100,0.9] )

            colorlist.append( [0,0.220, 0.039, 0.235] )
            colorlist.append( [20,0.390, 0.050, 0.330] )
            colorlist.append( [40,0.700, 0.200, 0.300] )
            colorlist.append( [100,1.00, 1.00, 0.20] )

        elif ('He-electrons' in hf.filename) or ('He-' in hf.filename) or ('N-electrons' in hf.filename):
            stype[j] = 'witness'

            opalist.append( [0,0.0] )
            opalist.append( [1,0.4] )
            opalist.append( [5,0.6] )
            opalist.append( [10,0.8] )
            opalist.append( [100,0.98] )

            colorlist.append( [0,0.627, 0.125, 0.235] )
            colorlist.append( [10,0.700, 0.200, 0.300] )
            colorlist.append( [100,1.00, 1.00, 0.20] )

        elif "e2" in hf.filename:
            stype[j] = 'laser'

            opalist.append( [0,0.6] )
            opalist.append( [35,0.2] )
            opalist.append( [40,0.0] )
            opalist.append( [60,0.0] )
            opalist.append( [65,0.2] )
            opalist.append( [100,0.6] )

            colorlist.append( [0,0.427, 0.588, 0.811] )
            colorlist.append( [100,0.972, 0.647, 0.145] )


        for k, opa in enumerate(opalist) :
            opacity[j].AddPoint(stop(opa[0]),opa[1])
        for k, col in enumerate(colorlist) :
            color[j].AddRGBPoint(stop(col[0]),col[1],col[2],col[3])

        fcolmap = open(fcolname,'a')
        #np.savetxt(fcolmap,stype[j])
        fcolmap.write('Species kind : %s\n' % str(stype[j]))
        # Go back to original values of MinMax
        if args.log & (args.lden & (stype[j] == 'plasma')) :
            localden = np.power(10,localden)
            Min[j] = np.power(10,Min[j]) * (baseden/localden)
            Max[j] = np.power(10,Max[j]) * (baseden/localden)
        elif args.lden & (stype[j] == 'plasma') :
            Min[j] = Min[j] * (baseden/localden)
            Max[j] = Max[j] * (baseden/localden)

        np.savetxt(fcolmap,(Min[j],Max[j]),fmt='%.3e',header='MinMax')
        np.savetxt(fcolmap,opalist,fmt='%7.3f',header='Opacities')
        np.savetxt(fcolmap,colorlist,fmt='%7.3f',header='Colors')
        fcolmap.write('\n')

    volumeprop.SetColor(j,color[j])
    volumeprop.SetScalarOpacity(j,opacity[j])
    volumeprop.ShadeOff(j)
    #volumeprop.ShadeOn(j)


fcolmap.close()

# Add data components as a 4th dimension
# npdatamulti = np.stack((npdata[:]),axis=3)
# Alternative way compatible with earlier versions of numpy
npdatamulti = np.concatenate([aux[...,np.newaxis] for aux in npdata], axis=3)
print('\nShape of the multi-component array: ', npdatamulti.shape,' Type: ',npdatamulti.dtype)

# For VTK to be able to use the data, it must be stored as a VTKimage.
# vtkImageImport imports raw data and stores it in the image.
dataImport = vtk.vtkImageImport()
dataImport.SetImportVoidPointer(npdatamulti)
dataImport.SetDataScalarTypeToFloat()
# Number of scalar components
dataImport.SetNumberOfScalarComponents(len(npdata))
# The following two functions describe how the data is stored
# and the dimensions of the array it is stored in.
dataImport.SetDataExtent(0, npdatamulti.shape[2]-1, 0, npdatamulti.shape[0]-1, 0, npdatamulti.shape[1]-1)
dataImport.SetWholeExtent(0, npdatamulti.shape[2]-1, 0, npdatamulti.shape[0]-1, 0, npdatamulti.shape[1]-1)
#dataImport.SetDataSpacing(dz,dy,dx)
#dataImport.SetDataOrigin(0.0,(axisy[1]+axisy[0])/(2.0),(axisx[1]+axisx[0])/(2.0))
dataImport.SetDataSpacing(udz,udx,udy)
dataImport.SetDataOrigin(0.0,uaxxmin,uaxymin)
dataImport.Update()

mapper = vtk.vtkGPUVolumeRayCastMapper()
mapper.SetAutoAdjustSampleDistances(1)
#mapper.SetSampleDistance(0.1)
#mapper.SetBlendModeToMaximumIntensity();

# Add data to the mapper
mapper.SetInputConnection(dataImport.GetOutputPort())

# The class vtkVolume is used to pair the previously declared volume
# as well as the properties to be used when rendering that volume.
volume = vtk.vtkVolume()
volume.SetMapper(mapper)
volume.SetProperty(volumeprop)


planeClip = vtk.vtkPlane()
planeClip.SetOrigin((axisz[0]+axisz[1])/2.0-axisz[0],0.0,0.0)
planeClip.SetNormal(0.0, 0.0, -1.0)
#mapper.AddClippingPlane(planeClip)

light = vtk.vtkLight()
light.SetColor(1.0, 0.0, 0.0)
light.SwitchOn()
light.SetIntensity(1)
# renderer.AddLight(light)

renderer = vtk.vtkRenderer()
# Set background
if args.white :
    renderer.SetBackground(0.9,0.9,0.9) # almost white
else :
    renderer.SetBackground(0,0,0) # black
# renderer.SetBackground(0.09,0.10,0.12)
# renderer.SetBackground(.1,.2,.3) # Background dark blue
# renderer.SetBackground(0.1,0.1,0.1)

# Other colors
nc = vtk.vtkNamedColors()
#renderer.SetBackground(nc.GetColor3d('MidnightBlue'))
#renderer.SetBackground(nc.GetColor3d('DeepSkyBlue'))
#renderer.SetBackground(0.09,0.10,0.12)

if args.txbck :
    renderer.TexturedBackgroundOn()

# Add the volume to the renderer ...
renderer.AddVolume(volume)


if args.axes :
    axes = vtk.vtkAxesActor()
    #    axes.SetTotalLength(0.2*(axisz[1]-axisz[0]),0.2*(axisy[1]-axisy[0]),0.2*(axisx[1]-axisx[0]))
    axes.SetTotalLength(6.28/3.0,6.28/3.0,6.28/3.0)
    axes.SetShaftTypeToLine()
    axes.SetNormalizedShaftLength(1, 1, 1)
    axes.SetNormalizedTipLength(0.1, 0.1, 0.1)
    propA = vtkTextProperty()
    propA.SetFontFamilyToArial()
    propA.ItalicOff()
    propA.BoldOff()
    propA.SetFontSize(24)
    if args.white :
        propA.SetColor(0.0,0.0,0.0)

    # the actual text of the axis label can be changed:
    axes.SetXAxisLabelText("z");
    axes.SetYAxisLabelText("y");
    axes.SetZAxisLabelText("x");

    axcolor = [0.9,0.9,0.9]
    if args.white :
        axcolor = [0.1,0.1,0.1]

    axes.GetXAxisShaftProperty().SetColor(axcolor);
    axes.GetXAxisTipProperty().SetColor(axcolor);
    axes.GetZAxisShaftProperty().SetColor(axcolor);
    axes.GetZAxisTipProperty().SetColor(axcolor);
    axes.GetYAxisShaftProperty().SetColor(axcolor);
    axes.GetYAxisTipProperty().SetColor(axcolor);

    axes.GetXAxisCaptionActor2D().GetTextActor().SetTextScaleModeToNone()
    axes.GetXAxisCaptionActor2D().SetCaptionTextProperty(propA)
    axes.GetZAxisCaptionActor2D().GetTextActor().SetTextScaleModeToNone()
    axes.GetZAxisCaptionActor2D().SetCaptionTextProperty(propA)
    axes.GetYAxisCaptionActor2D().GetTextActor().SetTextScaleModeToNone()
    axes.GetYAxisCaptionActor2D().SetCaptionTextProperty(propA)

    # The axes are positioned with a user transform
    # transform = vtk.vtkTransform()
    # transform.Translate(0.0,(axisy[1]-axisy[0])/(2.0),(axisx[1]-axisx[0])/(2.0))
    # axes.SetUserTransform(transform)

    renderer.AddActor(axes)



if args.cbar :
    scalarBar = []
    nbar = len(npdata)
    step = 1.0/nbar
    for i in range(0,len(npdata)):
        # Adding the scalar bar color palette
        scalarBar.append(vtkScalarBarActor())
        if 'e2' in hfl[i].filename :
            scalarBar[i].SetTitle("Ex")
        else :
            if args.log :
                scalarBar[i].SetTitle("log(n/np)")
            else :
                scalarBar[i].SetTitle("n/np")

        scalarBar[i].SetLookupTable(color[i]);
        scalarBar[i].SetOrientationToHorizontal();
        scalarBar[i].GetPositionCoordinate().SetCoordinateSystemToNormalizedViewport()
        scalarBar[i].SetPosition(i*step,0.02)
        scalarBar[i].SetPosition2(step,0.08)
        #print('%5.2f %5.2f' % (i*step,(i+1)*step-0.01) )
        propT = vtkTextProperty()
        propT.SetFontFamilyToArial()
        propT.ItalicOff()
        propT.BoldOn()
        propT.SetFontSize(28)
        propL = vtkTextProperty()
        propL.SetFontFamilyToArial()
        propL.ItalicOff()
        propL.BoldOff()
        propL.SetFontSize(22)
        if args.white :
            propT.SetColor(0.0,0.0,0.0)
            propL.SetColor(0.0,0.0,0.0)

        scalarBar[i].UseOpacityOff()
        #scalarBar[i].AnnotationTextScalingOff()
        scalarBar[i].SetTitleTextProperty(propT);
        scalarBar[i].SetLabelTextProperty(propL);
        scalarBar[i].SetLabelFormat("%5.2f")

        renderer.AddActor(scalarBar[i])


# Set window
window = vtk.vtkRenderWindow()
# ... and set window size.
window.SetSize(1280, 800)
if args.nowin :
    window.SetOffScreenRendering(1)

# add renderer
window.AddRenderer(renderer)

# Mouse control
interactor = vtk.vtkRenderWindowInteractor()
#style = vtkInteractorStyleTrackballCamera();
#interactor.SetInteractorStyle(style);
if args.nowin == 0 :
    interactor.SetRenderWindow(window)

# Camera
renderer.ResetCamera()
camera = renderer.GetActiveCamera()
#camera.SetPosition((axisz[1]-axisz[0])/2.0,74.0,0.0);
camera.SetFocalPoint((axisz[1]-axisz[0])/2.0,0.0,0.0);
camera.Zoom(args.zoom)
#camera.Roll(45)
camera.Elevation(args.elevation)
camera.Azimuth(args.azimuth)
#camera.ParallelProjectionOn()

#print('\nCamera position:    ', camera.GetPosition())
#print('\nCamera focal point: ', camera.GetFocalPoint())

window.Render()

if args.nowin == 0 :
    interactor.Initialize()
    interactor.Start()

# Output file
if args.test :
    foutname = './snapshot3d.png'
else :
    foutname = './%s/Plots/Snapshots3D/Snapshot3D-%s_%i.png' % (pData.GetPath(),args.sim,args.t)
PGlobals.mkdir(foutname)

# Write an EPS file.
# exp = vtk.vtkGL2PSExporter()  # Not working with openGL2 yet
#exp.SetRenderWindow(window)
#exp.SetFileFormatToPDF()
#exp.CompressOn()
#exp.SetSortToBSP()
#exp.SetFilePrefix("screenshot")
#exp.DrawBackgroundOn()
#exp.Write()

# Write to PNG file
w2if = vtk.vtkWindowToImageFilter()
w2if.SetInput(window)
w2if.Update();

writer = vtk.vtkPNGWriter()
writer.SetFileName(foutname)
writer.SetInputConnection(w2if.GetOutputPort())
writer.Write()

print('\n%s has been created.' % (foutname) )

# END
