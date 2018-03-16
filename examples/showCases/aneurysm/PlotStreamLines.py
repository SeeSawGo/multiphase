#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# create a new 'XML Image Data Reader'
yz_0vti = XMLImageDataReader(FileName=['/home/ss/palabos-v2.0r0/examples/showCases/aneurysm/xy_0.vti'])
yz_0vti.PointArrayStatus = ['p', 'u', 'voxel']

# set active source
SetActiveSource(yz_0vti)

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [1209, 692]

# show data in view
yz_0vtiDisplay = Show(yz_0vti, renderView1)

# trace defaults for the display properties.
yz_0vtiDisplay.Representation = 'Outline'
yz_0vtiDisplay.ColorArrayName = [None, '']
yz_0vtiDisplay.OSPRayScaleArray = 'p'
yz_0vtiDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
yz_0vtiDisplay.SelectOrientationVectors = 'None'
yz_0vtiDisplay.ScaleFactor = 0.10800000000000001
yz_0vtiDisplay.SelectScaleArray = 'None'
yz_0vtiDisplay.GlyphType = 'Arrow'
yz_0vtiDisplay.GlyphTableIndexArray = 'None'
yz_0vtiDisplay.DataAxesGrid = 'GridAxesRepresentation'
yz_0vtiDisplay.PolarAxes = 'PolarAxesRepresentation'
yz_0vtiDisplay.ScalarOpacityUnitDistance = 0.04804010758694445
yz_0vtiDisplay.Slice = 19
yz_0vtiDisplay.GaussianRadius = 0.054000000000000006
yz_0vtiDisplay.SetScaleArray = ['POINTS', 'p']
yz_0vtiDisplay.ScaleTransferFunction = 'PiecewiseFunction'
yz_0vtiDisplay.OpacityArray = ['POINTS', 'p']
yz_0vtiDisplay.OpacityTransferFunction = 'PiecewiseFunction'

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
yz_0vtiDisplay.ScaleTransferFunction.Points = [-0.02143377810716629, 0.0, 0.5, 0.0, 0.031839434057474136, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
yz_0vtiDisplay.OpacityTransferFunction.Points = [-0.02143377810716629, 0.0, 0.5, 0.0, 0.031839434057474136, 1.0, 0.5, 0.0]

# reset view to fit data
renderView1.ResetCamera()

# show data in view
yz_0vtiDisplay = Show(yz_0vti, renderView1)

# reset view to fit data
renderView1.ResetCamera()

# update the view to ensure updated data information
renderView1.Update()

# set scalar coloring
ColorBy(yz_0vtiDisplay, ('POINTS', 'u', 'X'))

# rescale color and/or opacity maps used to include current data range
yz_0vtiDisplay.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
yz_0vtiDisplay.SetScalarBarVisibility(renderView1, True)

# get color transfer function/color map for 'u'
uLUT = GetColorTransferFunction('u')

# set scalar coloring
ColorBy(yz_0vtiDisplay, ('POINTS', 'u', 'Magnitude'))

# rescale color and/or opacity maps used to exactly fit the current data range
yz_0vtiDisplay.RescaleTransferFunctionToDataRange(False, False)

# Update a scalar bar component title.
UpdateScalarBarsComponentTitle(uLUT, yz_0vtiDisplay)

# change representation type
yz_0vtiDisplay.SetRepresentationType('Surface')

# create a new 'Slice'
slice1 = Slice(Input=yz_0vti)
slice1.SliceType = 'Plane'
slice1.SliceOffsetValues = [0.0]

# init the 'Plane' selected for 'SliceType'
slice1.SliceType.Origin = [0.015, -0.0050000000000000044, 0.05]

# show data in view
slice1Display = Show(slice1, renderView1)

# trace defaults for the display properties.
slice1Display.Representation = 'Surface'
slice1Display.ColorArrayName = ['POINTS', 'u']
slice1Display.LookupTable = uLUT
slice1Display.OSPRayScaleArray = 'p'
slice1Display.OSPRayScaleFunction = 'PiecewiseFunction'
slice1Display.SelectOrientationVectors = 'None'
slice1Display.ScaleFactor = 0.10800000429153443
slice1Display.SelectScaleArray = 'None'
slice1Display.GlyphType = 'Arrow'
slice1Display.GlyphTableIndexArray = 'None'
slice1Display.DataAxesGrid = 'GridAxesRepresentation'
slice1Display.PolarAxes = 'PolarAxesRepresentation'
slice1Display.GaussianRadius = 0.054000002145767216
slice1Display.SetScaleArray = ['POINTS', 'p']
slice1Display.ScaleTransferFunction = 'PiecewiseFunction'
slice1Display.OpacityArray = ['POINTS', 'p']
slice1Display.OpacityTransferFunction = 'PiecewiseFunction'

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
slice1Display.ScaleTransferFunction.Points = [-0.02143377810716629, 0.0, 0.5, 0.0, 0.028186984360218048, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
slice1Display.OpacityTransferFunction.Points = [-0.02143377810716629, 0.0, 0.5, 0.0, 0.028186984360218048, 1.0, 0.5, 0.0]

# hide data in view
Hide(yz_0vti, renderView1)

# show color bar/color legend
slice1Display.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()

# create a new 'Surface Vectors'
surfaceVectors1 = SurfaceVectors(Input=slice1)
surfaceVectors1.SelectInputVectors = ['POINTS', 'u']

# show data in view
surfaceVectors1Display = Show(surfaceVectors1, renderView1)

# trace defaults for the display properties.
surfaceVectors1Display.Representation = 'Surface'
surfaceVectors1Display.ColorArrayName = ['POINTS', 'u']
surfaceVectors1Display.LookupTable = uLUT
surfaceVectors1Display.OSPRayScaleArray = 'p'
surfaceVectors1Display.OSPRayScaleFunction = 'PiecewiseFunction'
surfaceVectors1Display.SelectOrientationVectors = 'u'
surfaceVectors1Display.ScaleFactor = 0.10800000429153443
surfaceVectors1Display.SelectScaleArray = 'None'
surfaceVectors1Display.GlyphType = 'Arrow'
surfaceVectors1Display.GlyphTableIndexArray = 'None'
surfaceVectors1Display.DataAxesGrid = 'GridAxesRepresentation'
surfaceVectors1Display.PolarAxes = 'PolarAxesRepresentation'
surfaceVectors1Display.GaussianRadius = 0.054000002145767216
surfaceVectors1Display.SetScaleArray = ['POINTS', 'p']
surfaceVectors1Display.ScaleTransferFunction = 'PiecewiseFunction'
surfaceVectors1Display.OpacityArray = ['POINTS', 'p']
surfaceVectors1Display.OpacityTransferFunction = 'PiecewiseFunction'

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
surfaceVectors1Display.ScaleTransferFunction.Points = [-0.02143377810716629, 0.0, 0.5, 0.0, 0.028186984360218048, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
surfaceVectors1Display.OpacityTransferFunction.Points = [-0.02143377810716629, 0.0, 0.5, 0.0, 0.028186984360218048, 1.0, 0.5, 0.0]

# hide data in view
Hide(slice1, renderView1)

# show color bar/color legend
surfaceVectors1Display.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()

# create a new 'Mask Points'
maskPoints1 = MaskPoints(Input=surfaceVectors1)
maskPoints1.OnRatio = 1

# show data in view
maskPoints1Display = Show(maskPoints1, renderView1)

# trace defaults for the display properties.
maskPoints1Display.Representation = 'Surface'
maskPoints1Display.ColorArrayName = ['POINTS', 'u']
maskPoints1Display.LookupTable = uLUT
maskPoints1Display.OSPRayScaleArray = 'p'
maskPoints1Display.OSPRayScaleFunction = 'PiecewiseFunction'
maskPoints1Display.SelectOrientationVectors = 'u'
maskPoints1Display.ScaleFactor = 0.10800000429153443
maskPoints1Display.SelectScaleArray = 'None'
maskPoints1Display.GlyphType = 'Arrow'
maskPoints1Display.GlyphTableIndexArray = 'None'
maskPoints1Display.DataAxesGrid = 'GridAxesRepresentation'
maskPoints1Display.PolarAxes = 'PolarAxesRepresentation'
maskPoints1Display.GaussianRadius = 0.054000002145767216
maskPoints1Display.SetScaleArray = ['POINTS', 'p']
maskPoints1Display.ScaleTransferFunction = 'PiecewiseFunction'
maskPoints1Display.OpacityArray = ['POINTS', 'p']
maskPoints1Display.OpacityTransferFunction = 'PiecewiseFunction'

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
maskPoints1Display.ScaleTransferFunction.Points = [-0.02143377810716629, 0.0, 0.5, 0.0, 0.028186984360218048, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
maskPoints1Display.OpacityTransferFunction.Points = [-0.02143377810716629, 0.0, 0.5, 0.0, 0.028186984360218048, 1.0, 0.5, 0.0]

# hide data in view
Hide(surfaceVectors1, renderView1)

# show color bar/color legend
maskPoints1Display.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()

# create a new 'Stream Tracer With Custom Source'
streamTracerWithCustomSource1 = StreamTracerWithCustomSource(Input=surfaceVectors1,
    SeedSource=maskPoints1)
streamTracerWithCustomSource1.Vectors = ['POINTS', 'u']
streamTracerWithCustomSource1.MaximumStreamlineLength = 1.0800000429153442

# show data in view
streamTracerWithCustomSource1Display = Show(streamTracerWithCustomSource1, renderView1)

# trace defaults for the display properties.
streamTracerWithCustomSource1Display.Representation = 'Surface'
streamTracerWithCustomSource1Display.ColorArrayName = ['POINTS', 'u']
streamTracerWithCustomSource1Display.LookupTable = uLUT
streamTracerWithCustomSource1Display.OSPRayScaleArray = 'AngularVelocity'
streamTracerWithCustomSource1Display.OSPRayScaleFunction = 'PiecewiseFunction'
streamTracerWithCustomSource1Display.SelectOrientationVectors = 'Normals'
streamTracerWithCustomSource1Display.ScaleFactor = 0.05379930194467306
streamTracerWithCustomSource1Display.SelectScaleArray = 'AngularVelocity'
streamTracerWithCustomSource1Display.GlyphType = 'Arrow'
streamTracerWithCustomSource1Display.GlyphTableIndexArray = 'AngularVelocity'
streamTracerWithCustomSource1Display.DataAxesGrid = 'GridAxesRepresentation'
streamTracerWithCustomSource1Display.PolarAxes = 'PolarAxesRepresentation'
streamTracerWithCustomSource1Display.GaussianRadius = 0.02689965097233653
streamTracerWithCustomSource1Display.SetScaleArray = ['POINTS', 'AngularVelocity']
streamTracerWithCustomSource1Display.ScaleTransferFunction = 'PiecewiseFunction'
streamTracerWithCustomSource1Display.OpacityArray = ['POINTS', 'AngularVelocity']
streamTracerWithCustomSource1Display.OpacityTransferFunction = 'PiecewiseFunction'

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
streamTracerWithCustomSource1Display.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.1757813367477812e-38, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
streamTracerWithCustomSource1Display.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.1757813367477812e-38, 1.0, 0.5, 0.0]

# hide data in view
Hide(surfaceVectors1, renderView1)

# hide data in view
Hide(maskPoints1, renderView1)

# show color bar/color legend
streamTracerWithCustomSource1Display.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()

#### saving camera placements for all active views

# current camera placement for renderView1
renderView1.CameraPosition = [2.3941473941645692, 0.007364569314465025, 0.2475000098347664]
renderView1.CameraFocalPoint = [0.014999999664723873, -0.004999995231628418, 0.2475000098347664]
renderView1.CameraViewUp = [0.0, 0.0, 1.0]
renderView1.CameraParallelScale = 0.6157769725177792

#### uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).