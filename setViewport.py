# -*- coding:UTF-8 -*-

from abaqus import *
import odbAccess
import visualization
from abaqusConstants import *

#---------------------------------------------------------------------------

#将当前视口中的 ODB 文件赋予变量 vp

vp = session.viewports[session.currentViewportName]

# 将背景色改为白色 white

#session.graphicsOptions.setValues(backgroundColor='#FFFFFF',backgroundStyle=SOLID)

# 将默认图注替换为用户自定义图注 

font="-*-Times New Roman-medium-r-*-*-*-140-*-*-*-*-iso8859-1"
vp.viewportAnnotationOptions.setValues(legendFont=font, titleFont=font,stateFont=font)

#云图连续显示

vp.odbDisplay.contourOptions.setValues(contourStyle=CONTINUOUS)

#不显示网格

vp.odbDisplay.commonOptions.setValues(visibleEdges=FEATURE,
    deformationScaling=UNIFORM, uniformScaleFactor=1)
    
vp.view.fitView()
