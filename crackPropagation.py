# -*- coding:UTF-8 -*-

from abaqus import *
from part import *
from material import *
from section import *
from assembly import *
from step import *
from interaction import *
from load import *
from mesh import *
from job import *
from sketch import *
from visualization import *
from odbAccess import *
from connectorBehavior import *
from caeModules import *
from jobMessage import *
from abaqusConstants import *
import regionToolset, displayGroupMdbToolset as dgm
import mesh
import testUtils
testUtils.setBackwardCompatibility()

#---------------------------------------------------------------------------

def getModel(modelName):
    #获得模型
    
    model=mdb.models[modelName]
    return model
    
#---------------------------------------------------------------------------

def createViewport(modelName):
    #创建新视口显示模型和分析结果
    
    viewport= session.Viewport(name=modelName)
    viewport.makeCurrent()
    viewport.maximize()
    return viewport

#---------------------------------------------------------------------------

def getPart(model):
    #获得部件
    
    partName=model.parts.keys()[-1]
    part=model.parts[partName]
    return part

#---------------------------------------------------------------------------

def getAssembly(model, viewport):
    #获得装配体
    
    assembly=model.rootAssembly
    viewport.setValues(displayedObject=assembly)
    assembly.DatumCsysByDefault(CARTESIAN)
    instanceName=assembly.instances.keys()[-1]
    partInstance=assembly.instances[instanceName]
    return assembly, partInstance

#---------------------------------------------------------------------------

def getFractureArgs(modelName, stepNum, crackNum):
    #从输出数据库中提取应力强度因子、转角、J-积分
    
    odb=openOdb(path=modelName+str(stepNum)+'.odb')   # !!! #加上路径
    
    outputName='CRACK'+str(crackNum)
    crackName='CRACK'+str(crackNum)
    regName='CRACKRG'+str(crackNum)
    
    lastStep=odb.steps.values()[-1]
    
    crackArgsOutputs=lastStep.historyRegions['ElementSet  PIBATCH'].historyOutputs
    
    K1=(crackArgsOutputs['K1 at '+outputName+'_'+crackName+'_'+regName+'-1__Contour_2'].data[0][1]+\
           crackArgsOutputs['K1 at '+outputName+'_'+crackName+'_'+regName+'-1__Contour_3'].data[0][1]+\
           crackArgsOutputs['K1 at '+outputName+'_'+crackName+'_'+regName+'-1__Contour_4'].data[0][1])/3
       
    K2=(crackArgsOutputs['K2 at '+outputName+'_'+crackName+'_'+regName+'-1__Contour_2'].data[0][1]+\
           crackArgsOutputs['K2 at '+outputName+'_'+crackName+'_'+regName+'-1__Contour_3'].data[0][1]+\
           crackArgsOutputs['K2 at '+outputName+'_'+crackName+'_'+regName+'-1__Contour_4'].data[0][1])/3
    
    cpd=(crackArgsOutputs['Cpd at '+outputName+'_'+crackName+'_'+regName+'-1__Contour_2'].data[0][1]+\
            crackArgsOutputs['Cpd at '+outputName+'_'+crackName+'_'+regName+'-1__Contour_3'].data[0][1]+\
            crackArgsOutputs['Cpd at '+outputName+'_'+crackName+'_'+regName+'-1__Contour_4'].data[0][1])/3
       
    jIntegral=(crackArgsOutputs['JKs at '+outputName+'_'+crackName+'_'+regName+'-1__Contour_2'].data[0][1]+\
           crackArgsOutputs['JKs at '+outputName+'_'+crackName+'_'+regName+'-1__Contour_3'].data[0][1]+\
           crackArgsOutputs['JKs at '+outputName+'_'+crackName+'_'+regName+'-1__Contour_4'].data[0][1])/3
           
    odb.close()
    
    return K1, cpd, K2, jIntegral

'''
    print '---------------------------------------------------------------------------'
    print "The crack's fracture arguments:\n"
    print 'SIF K1: %.3f' % K1
    print 'SIF K2: %.3f' % K2
    print 'MTS DIRECTION (DEG): %.3f' % cpd
    print 'J-Integral: %.4f' % jIntegral
'''

#---------------------------------------------------------------------------

def getOldCrackTipCoord(assembly, crackNum):
    #获得裂尖坐标
    
    setName=assembly.sets.keys()
    tipLocation=assembly.sets['crackPt'+str(crackNum)].vertices[-1].pointOn[0]
    return tipLocation

#---------------------------------------------------------------------------

def calcNewCrackTipCoord(cpd, growth):
    #计算新的裂尖位置
    
    cpd=pi/180*cpd    #度转化为弧度
    
    newX=growth*cos(cpd)
    newY=growth*sin(cpd)

    return newX, newY

#---------------------------------------------------------------------------

def calcRadiusPoint(radius, cpd):
    #计算裂尖圆形分割圆周上的一点的坐标
    
    cpd=pi/180*cpd
    
    radiusX=radius*cos(cpd)
    radiusY=radius*sin(cpd)
    
#    print "radiusX=%f,radiusY=%f" %(radiusX, radiusY)
    return radiusX, radiusY

#---------------------------------------------------------------------------

def deleteCrackTipPartition(assembly, crackNum):
    #删除裂尖的圆周分割
    
    for id in range(crackNum, 0, -1):
        assembly.deleteFeatures(featureNames=('innerPartition'+str(id), 
                                                                'outterPartition'+str(id) ) )

#---------------------------------------------------------------------------

def deleteSets(assembly, crackNum):
    #删除集合
    
    del assembly.sets['crackPt'+str(crackNum)]
    del assembly.sets['crackRg'+str(crackNum)]
    del assembly.sets['crackRgOut'+str(crackNum)]
    del assembly.sets['seamIn'+str(crackNum)]
    del assembly.sets['seamOut'+str(crackNum)]
    del assembly.sets['innerPartition'+str(crackNum)+'Down']
    del assembly.sets['innerPartition'+str(crackNum)+'Up']
    del assembly.sets['outterPartition'+str(crackNum)+'Down']
    del assembly.sets['outterPartition'+str(crackNum)+'Up']

#---------------------------------------------------------------------------

def createEdgeSetOnPart(part, setName, coord):
    
    edges1=part.edges
    s1 = part.edges.findAt(coord,)
    side1Edges1 = edges1[s1.index:(s1.index+1)]
    part.Set(edges=side1Edges1, name=setName)

#---------------------------------------------------------------------------

def  createEdgeSet(part, assembly, partInstance, setName, coord):
    #创建线集合
    
    assembly.regenerate()
    
    edges1=partInstance.edges
    s1 = partInstance.edges.findAt(coord,)
    side1Edges1 = edges1[s1.index:(s1.index+1)]
    assembly.Set(edges=side1Edges1, name=setName)

#---------------------------------------------------------------------------

def createPtSet(assembly, partInstance, coord, setName):
    #为裂尖创建集合
    
    assembly.regenerate()
    
    verts1 = partInstance.vertices
    vxx = partInstance.vertices.findAt(coord,)
    tip = verts1[vxx.index:(vxx.index+1)]
    assembly.Set(vertices=tip, name=setName)

#---------------------------------------------------------------------------

def crackRgSet(assembly, partInstance, coord, setName):
    #为区域创建集合
    
    assembly.regenerate()
    
    face = partInstance.faces.findAt((coord, ), )
    assembly.Set(faces=face, name=setName)

#---------------------------------------------------------------------------

def modelCentroid(part):
    #获得模型的形心
    
    face=part.faces.findAt((0.0, 0.0, 0.0), )
    centroid=face.getCentroid()
    centroid=centroid[0][0]
    
    return centroid

#---------------------------------------------------------------------------

def createSeamPartiton(part, model, assembly, partInstance, centroid, oldX, oldY, oldZ, reNewX, reNewY, newX, newY, newZ, seamName):
    #在给定的增长量下创建seam分割
    
    if oldX>centroid:
        orient = LEFT
    elif oldX<centroid:
        orient = RIGHT
    
    face=part.faces.findAt((oldX, oldY, oldZ), )
    t=part.MakeSketchTransform(sketchPlane=face, sketchPlaneSide=SIDE1,
                               sketchOrientation=orient, origin=(oldX, oldY, oldZ))
    mySketch=model.Sketch(name='plateProfile', sheetSize=200.0, transform=t)
    mySketch.setPrimaryObject(option=SUPERIMPOSE)
    part.projectReferencesOntoSketch(sketch=mySketch, filter=COPLANAR_EDGES)
    mySketch.Line(point1=(0.0, 0.0),point2=(reNewX, reNewY) )
    pickedFaces=part.faces.findAt((oldX, oldY, oldZ), )
    part.PartitionFaceBySketch(faces=pickedFaces, sketch=mySketch)
    part.features.changeKey(fromName='Partition face-1', toName=seamName)    #裂纹每扩展一次，seam增加1    
    part.regenerate()
    
    createEdgeSetOnPart(part=part, setName=seamName, coord=(newX, newY, newZ))    #给创建的seam创建集合
    
    mySketch.unsetPrimaryObject()
    del model.sketches['plateProfile']

#---------------------------------------------------------------------------

def createCrackTipPartition(model, part, assembly, partInstance, radius, cpd, partitionName, centroid, oldX, oldY, oldZ, newX, newY, newZ):
    #创建裂尖圆形分割
    # !!!    如果newX>centroid=LEFT；如果newX<centroid,参数sketchOrientation=RIGHT
    # centroid为图形的形心横坐标
    
    if newX>centroid:
        orient = LEFT
    elif newX<centroid:
        orient = RIGHT
    
    radiusX, radiusY=calcRadiusPoint(radius, cpd)
    
    faces1=partInstance.faces.findAt((oldX, oldY, oldZ), )
    t=assembly.MakeSketchTransform(sketchPlane=faces1, sketchPlaneSide=SIDE1,
                                   sketchOrientation=orient, origin=(newX, newY, newZ))
    mySketch=model.Sketch(name='plateProfile', sheetSize=200.0, transform=t)
    assembly.projectReferencesOntoSketch(sketch=mySketch,filter=ALL_EDGES)
    mySketch.CircleByCenterPerimeter(center=(0.0, 0.0),point1=(radiusX, radiusY))    
    pickedFaces=partInstance.faces.findAt((newX, newY, newZ), )
    assembly.PartitionFaceBySketch(faces=pickedFaces, sketch=mySketch)
    assembly.features.changeKey(fromName='Partition face-1', toName=partitionName)
    assembly.regenerate()
    
    #为分割创建集合
    X1, Y1= calcRadiusPoint(radius, cpd+90)
    createEdgeSet(part=part, assembly=assembly, partInstance=partInstance, setName=partitionName+'Up', 
              coord=(newX+X1, newY+Y1, 0.0))
    X2, Y2= calcRadiusPoint(radius, cpd-90)
    createEdgeSet(part=part, assembly=assembly, partInstance=partInstance, setName=partitionName+'Down', 
              coord=(newX+X2, newY+Y2, 0.0))
    
    del model.sketches['plateProfile']    

#---------------------------------------------------------------------------

def assignSeam(assembly, partInstance, crackNum, increaseNum):
    #指定裂纹
    
    for num in range(1, increaseNum+1):
        assembly.engineeringFeatures.assignSeam(regions=partInstance.sets['seam'+str(crackNum)+'-'+str(num)])
    
    assembly.engineeringFeatures.assignSeam(regions=assembly.sets['seamOut'+str(crackNum)])
    assembly.engineeringFeatures.assignSeam(regions=assembly.sets['seamIn'+str(crackNum)])

#---------------------------------------------------------------------------

def createContourIntegral(assembly, partInstance, crackName, oldX, oldY, oldZ, reNewX, reNewY, newX, newY, newZ, crackNum):
    #设置裂纹属性
    
    verts1 = partInstance.vertices
    vxx=partInstance.vertices.findAt((newX, newY, 0.0),)    
    vyy=partInstance.vertices.findAt((oldX+reNewX*2, oldY+reNewY*2, 0.0),)
    crackFront=assembly.sets['crackRg'+str(crackNum)]
    crackTip=assembly.sets['crackPt'+str(crackNum)]
    assembly.engineeringFeatures.ContourIntegral(name=crackName,
         symmetric=OFF, crackFront=crackFront, crackTip=crackTip,    
         extensionDirectionMethod=Q_VECTORS, qVectors=((vxx,vyy),),
         midNodePosition=0.25, collapsedElementAtTip=DUPLICATE_NODES) 

#---------------------------------------------------------------------------

def seedEdge(assembly, partInstance, initCrackLength, circularNum, radialNum, increaseNum, crackNum, flag=0):
    #设置边界上的种子
    
    assembly.seedEdgeByNumber(edges=assembly.sets['seamIn'+str(crackNum)].edges, number=1, constraint=FIXED)
    assembly.seedEdgeByNumber(edges=assembly.sets['seamOut'+str(crackNum)].edges, number=radialNum, constraint=FIXED)
    assembly.seedEdgeByNumber(edges=assembly.sets['innerPartition'+str(crackNum)+'Down'].edges, number=circularNum, constraint=FIXED)
    assembly.seedEdgeByNumber(edges=assembly.sets['innerPartition'+str(crackNum)+'Up'].edges, number=circularNum, constraint=FIXED)
    assembly.seedEdgeByNumber(edges=assembly.sets['outterPartition'+str(crackNum)+'Down'].edges, number=circularNum, constraint=FIXED)
    assembly.seedEdgeByNumber(edges=assembly.sets['outterPartition'+str(crackNum)+'Up'].edges, number=circularNum, constraint=FIXED)
    
    if flag == 0:
        assembly.seedEdgeByNumber(edges=partInstance.sets['seam'+str(crackNum)+'-'+str(1)].edges, number=int(initCrackLength*2), constraint=FINER)
        for num in range(2, increaseNum+1):
            assembly.seedEdgeByNumber(edges=partInstance.sets['seam'+str(crackNum)+'-'+str(num)].edges, number=2, constraint=FINER)
    elif flag == 1:
        assembly.seedEdgeByNumber(edges=assembly.sets['seam'+str(crackNum)].edges, number=int(initCrackLength*2), constraint=FINER)

#---------------------------------------------------------------------------

def assignMeshControl(assembly, partInstance, crackNum):
    #设置网格划分控制参数
    
    assembly.setMeshControls(regions=assembly.sets['wholePart'].faces, elemShape=TRI,
        technique=FREE, algorithm=MEDIAL_AXIS)
    
    assembly.setMeshControls(regions=assembly.sets['crackRgOut'+str(crackNum)].faces, elemShape=QUAD,
        technique=SWEEP, algorithm=MEDIAL_AXIS)

    assembly.setMeshControls(regions=assembly.sets['crackRg'+str(crackNum)].faces, elemShape=TRI,
        technique=FREE)

#---------------------------------------------------------------------------

def assignElementType(assembly, partInstance, crackNum):
    #设置单元类型
    
    elemType1 = mesh.ElemType(elemCode=CPS8, elemLibrary=STANDARD)
    elemType2 = mesh.ElemType(elemCode=CPS6, elemLibrary=STANDARD)
    
    assembly.setElementType(regions=assembly.sets['wholePart'],elemTypes=(elemType1,))
    assembly.setElementType(regions=assembly.sets['crackRgOut'+str(crackNum)],elemTypes=(elemType1,))
    assembly.setElementType(regions=assembly.sets['crackRg'+str(crackNum)],elemTypes=(elemType2,))

#---------------------------------------------------------------------------

def meshInstance(assembly, partInstance):
    #划分网格
    
    partInstances =(partInstance, )
    assembly.generateMesh(regions=partInstances)

#---------------------------------------------------------------------------

def createHistoryOutput(model, outputName, crackName, crackNum, flag=0):
    #为裂纹设置输出
    
    if flag == 0:
        del model.historyOutputRequests['crack'+str(crackNum)]

    stepName=model.steps.values()[-1].name
    model.HistoryOutputRequest(name=outputName, createStepName=stepName,      
                        contourIntegral=crackName, contourType=K_FACTORS, 
                        kFactorDirection=MTS, numberOfContours=4)

#---------------------------------------------------------------------------

def createJob(modelName, increaseNum):
    #创建job
    
    job = mdb.Job(name=modelName+str(increaseNum), model=modelName,
        description='Simple Crack Analysis', numCpus=4, numDomains=4)
    return job

#---------------------------------------------------------------------------

def submitJob(job, modelName, increaseNum):
    #提交分析
    
    job.submit()
    job.waitForCompletion()
    mdb.saveAs(pathName=modelName+str(increaseNum))

#---------------------------------------------------------------------------

def getInitData():
    #获得初始数据：包括CAE文件所在路径，裂纹名称，输出数据库名称，
    #裂纹初始长度，裂纹初始角度，裂纹增长量，径向和周向网格划分数，
    #以及材料的断裂韧度
    
    field1=(('CAE File Path:', ''), ('Radial Mesh Number:', ''), ('Circular Mesh Number:', ''),
             ('the Amount of Growth:', ''), ('KIC:', ''))
    field2=(('Crack Initial Length:', ''), ('Crack Initial Angle:', ''), ('Propagation Direction:', ''))
    cracks=[]
    
    caePath, radialNum, circularNum, growth, KIC=getInputs(fields=field1, label='Init Args')
    crackNum=1
    while True:
        crackName='crack'+str(crackNum)
        initCrackLength, initCpd, direct=getInputs(fields=field2, label='Crack '+str(crackNum)+' Args')
        crackName={}
        crackName['initCrackLength']=float(initCrackLength)
        crackName['initCpd']=float(initCpd)
        crackName['name']='crack'+str(crackNum)
        crackName['direct']=direct
        cracks.append(crackName)
        reply=getWarningReply(message='Is there any cracks?', buttons=(YES, NO))
        if reply == YES:
            crackNum+=1
        elif reply == NO:
            break
    
    caeFilePath=caePath.replace('\\', '/')
    growth=float(growth)
    radialNum=int(radialNum)
    circularNum=int(circularNum)
    KIC=float(KIC)

    
    return caeFilePath, growth, radialNum, circularNum, KIC, crackNum, cracks
    
#---------------------------------------------------------------------------

def calcCompSIF(K1, K2, theta):
    #计算复合应力强度因子
    
    theta=pi/180*theta
    
    sTheta=sin(theta)
    cTheta=cos(theta/2)
    
    Ke=cTheta*(K1*(cTheta**2)-1.5*K2*sTheta)
    
    return Ke

#---------------------------------------------------------------------------

def storeCrackData(path, cracks, increaseNum, num):
    #将每条裂纹的关键参数存入文件
    
    file=open(path+'Data.txt', 'a+')
    try:
        if num+1 == 1:
            file.write('************************************\n')
            file.write('************************************\n')
            file.write('Step '+str(increaseNum)+'\n')
            file.write(' \n')
        
        file.write('Crack '+str(num+1)+':\n')
        file.write(' \n')
        file.write('            K1                                K2                              J-Integral\n')
        file.write('       '+str(cracks[num]['K1'])+'                     '+str(cracks[num]['K2'])+'                    '+str(cracks[num]['jIntegral'])+'\n')
        file.write(' \n')
        file.write('            Ke\n')
        file.write('       '+str(cracks[num]['Ke'])+'\n')
        file.write(' \n')
        file.write('    MTS DIRECTION (DEG)               Absolute Angle (DEG)\n')
        file.write('     '+str(cracks[num]['relativeCpd'])+'                      '+str(cracks[num]['cpd'])+'\n')
        file.write(' \n')
        file.write('    Crack Tip Location:             X                               Y\n')
        file.write('                                 '+str(cracks[num]['oldX'])+'                       '+str(cracks[num]['oldY'])+'\n')
        file.write(' \n')
        file.write('------------------------------------\n')            
        file.write(' \n')
    finally:
        file.close()

#---------------------------------------------------------------------------

def initModel(caeFilePath, modelName, crackNum, cracks, growth, circularNum, radialNum, increaseNum):
    '''模型初始化函数，此函数用于简化建模过程，在运行程序前只需指定裂纹尖端即可,
      裂纹尖端的分割等工作均在此模块中完成。
      建模时，seam的名称需为seam1-1形式，第一个1表示裂纹序号，第二个1表示是
      扩展的第一步，并在 par中创建同名的集合，在assembly中为裂尖创建集合，名称
      应为crackPt1的形式，数字为裂纹序号
    '''
    
    myModel=getModel(modelName)
    
    myViewport=createViewport(modelName)   
    myPart=getPart(model=myModel)
    myAssembly, myPartInstance=getAssembly(model=myModel, viewport=myViewport)
    centroid=modelCentroid(part=myPart) 
    
    for num in range(crackNum):
        cracks[num]['oldX'], cracks[num]['oldY'], cracks[num]['oldZ']=getOldCrackTipCoord(assembly=myAssembly, crackNum=num+1)
            
        createCrackTipPartition(model=myModel, part=myPart, assembly=myAssembly, partInstance=myPartInstance, centroid=centroid, 
                                radius=growth, cpd=cracks[num]['initCpd'], partitionName='outterPartition'+str(num+1), 
                                oldX=cracks[num]['oldX'], oldY=cracks[num]['oldY'], oldZ=cracks[num]['oldZ'], 
                                newX=cracks[num]['oldX'], newY=cracks[num]['oldY'], newZ=cracks[num]['oldZ'])    #外圆分割
        createCrackTipPartition(model=myModel, part=myPart, assembly=myAssembly, partInstance=myPartInstance, centroid=centroid, 
                                radius=growth*0.2, cpd=cracks[num]['initCpd'], partitionName='innerPartition'+str(num+1), 
                                oldX=cracks[num]['oldX'], oldY=cracks[num]['oldY'], oldZ=cracks[num]['oldZ'], 
                                newX=cracks[num]['oldX'], newY=cracks[num]['oldY'], newZ=cracks[num]['oldZ'])    #内圆分割
            
        createEdgeSet(part=myPart, assembly=myAssembly, partInstance=myPartInstance, setName='seamIn'+str(num+1), 
                             coord=(cracks[num]['oldX'], cracks[num]['oldY'], 0.0))
        createEdgeSet(part=myPart, assembly=myAssembly, partInstance=myPartInstance, setName='seamOut'+str(num+1), 
                             coord=(cracks[num]['oldX']-growth/2*cos(cracks[num]['initCpd']*pi/180), 
                                        cracks[num]['oldY']-growth/2*sin(cracks[num]['initCpd']*pi/180), 0.0  ))
        crackRgSet(assembly=myAssembly, partInstance=myPartInstance, setName='crackRg'+str(num+1), 
                             coord=(cracks[num]['oldX'], cracks[num]['oldY'], 0.0) )
        crackRgSet(assembly=myAssembly, partInstance=myPartInstance, setName='crackRgOut'+str(num+1), 
                             coord=(cracks[num]['oldX']-growth/2*cos(cracks[num]['initCpd']*pi/180),
                                        cracks[num]['oldY']-growth/2*sin(cracks[num]['initCpd']*pi/180), 0.0) )
        crackRgSet(assembly=myAssembly, partInstance=myPartInstance, coord=(0.0, 0.0, 0.0), setName='wholePart')
        createEdgeSet(part=myPart, assembly=myAssembly, partInstance=myPartInstance, setName='seam'+str(num+1), 
                             coord=(cracks[num]['oldX']-growth*1.5*cos(cracks[num]['initCpd']*pi/180), 
                                        cracks[num]['oldY']-growth*1.5*sin(cracks[num]['initCpd']*pi/180), 0.0  ))    #此集合仅用于模型初始化中
        
        myAssembly.engineeringFeatures.assignSeam(regions=myPartInstance.sets['seam'+str(num+1)+'-'+str(1)])
        createContourIntegral(assembly=myAssembly, partInstance=myPartInstance, crackName=cracks[num]['name'], 
                             oldX=cracks[num]['oldX'], oldY=cracks[num]['oldY'], oldZ=cracks[num]['oldZ'], 
                             reNewX=growth/2*cos(cracks[num]['initCpd']*pi/180), 
                             reNewY=growth/2*sin(cracks[num]['initCpd']*pi/180), 
                             newX=cracks[num]['oldX'], newY=cracks[num]['oldY'], newZ=cracks[num]['oldZ'], 
                             crackNum=num+1)
        createHistoryOutput(model=myModel, outputName='crack'+str(num+1), 
                             crackName=cracks[num]['name'],  crackNum=num+1, flag=1)     #outputName用crackName加上序号代替
        
        seedEdge(assembly=myAssembly, partInstance=myPartInstance, initCrackLength=cracks[num]['initCrackLength']-growth,
                         circularNum=circularNum, radialNum=radialNum, increaseNum=increaseNum, crackNum=num+1, flag=1)
        assignMeshControl(assembly=myAssembly, partInstance=myPartInstance, crackNum=num+1)
        assignElementType(assembly=myAssembly, partInstance=myPartInstance, crackNum=num+1)
   
    meshInstance(assembly=myAssembly, partInstance=myPartInstance)
      
    myJob=createJob(modelName=modelName, increaseNum=1)
    submitJob(job=myJob, modelName=modelName, increaseNum=1)

#---------------------------------------------------------------------------
    
# ！！！ 创建完seam后建立集合，创建完裂尖分割后创建集合in和out，in和out 用于划分网格，集合seam在下一步分析中应用
# ！！！ 每一个循环开始之后，要删除裂尖分割（2个），集合（9个）
# ！！！ 最后应加上裂缝的初始角度
# ！！！ CAE建模时，part中建立集合seam1
# ！！！ 获得K1，cpd之后要关闭数据库  
# ！！！ 每条裂纹的角度都是从X轴正向开始
# ！！！ 每条裂纹都以crack+数字序号命名，每一个crack对应一段seam
















