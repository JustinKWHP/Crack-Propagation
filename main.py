# -*- coding:UTF-8 -*-

#导入python模块
import sys, pprint, inspect, os
#获取脚本文件所在的路径，将模块crackPropagation加入python的搜索路径
scriptName = inspect.getfile(inspect.currentframe())
scriptPath = os.path.abspath(os.path.dirname(scriptName))
scriptPathChged = scriptPath.replace('\\', '/')
sys.path.append(scriptPathChged)

#导入ABAQUS的python模块
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

#导入自定义模块
from crackPropagation import *

#---------------------------------------------------------------------------

if __name__ == '__main__' :
    
    #-----------------------------------------------------------------------
    #获得初始参数
    
    increaseNum=1    #裂纹每增长一次，值增加1
    countTerminate=0    #记录终止扩展的裂纹条数
    caeFilePath, growth, radialNum, circularNum, KIC, crackNum, cracks=getInitData()
    
    #-----------------------------------------------------------------------
    #计算第一个模型
    openMdb(pathName=caeFilePath+str(1))    #打开CAE数据库
    modelName=mdb.models.keys()[-1]    #获得模型名称
    initModel(caeFilePath=caeFilePath, modelName=modelName, crackNum=crackNum, cracks=cracks, 
                 growth=growth, circularNum=circularNum, radialNum=radialNum, increaseNum=increaseNum)
    
    #-----------------------------------------------------------------------
    #裂纹每扩展一次，就执行一次循环
    
    flag=True
    while flag:

        openMdb(pathName=caeFilePath+str(increaseNum))    #打开CAE数据库
        modelName=mdb.models.keys()[-1]    #获得模型名称
    
        myModel=getModel(modelName)
        myViewport=createViewport(modelName)   
        myPart=getPart(model=myModel)
        myAssembly, myPartInstance=getAssembly(model=myModel, viewport=myViewport)
        centroid=modelCentroid(part=myPart)    
    
        print "************************************"
        print "************************************"
        print "Step %d:"  %  increaseNum
        print ' '
        
        #-----------------------------------------------------------------------
        #删除裂尖的圆形分割，后生成的分割先删除
        
        deleteCrackTipPartition(assembly=myAssembly, crackNum=crackNum)
        
        #-----------------------------------------------------------------------
        #对模型中的每一个裂纹，执行相应的操作
        
        for num in range(crackNum):
            
            K1=0.0
            K2=0.0
            Ke=0.0
            relativeCpd=0.0
            K1, relativeCpd, K2, jIntegral=getFractureArgs(modelName=modelName, stepNum=increaseNum, crackNum=num+1)    #获取裂纹的参数值
            cracks[num]['K1']=K1
            cracks[num]['relativeCpd']=relativeCpd
            cracks[num]['K2']=K2
            cracks[num]['jIntegral']=jIntegral
            cracks[num]['Ke']=calcCompSIF(K1=cracks[num]['K1'], K2=cracks[num]['K2'],
                                          theta=cracks[num]['relativeCpd'])
    
            #计算裂纹的扩展角度
            if increaseNum==1:
                if cracks[num]['direct']=='+':
                    cracks[num]['cpd']=cracks[num]['initCpd']+cracks[num]['relativeCpd'] 
                elif cracks[num]['direct']=='-':
                    cracks[num]['cpd']=cracks[num]['initCpd']-cracks[num]['relativeCpd']
            else:
                if cracks[num]['direct']=='+':
                    cracks[num]['cpd']=cracks[num]['cpd']+cracks[num]['relativeCpd']
                elif cracks[num]['direct']=='-':
                    cracks[num]['cpd']=cracks[num]['cpd']-cracks[num]['relativeCpd']

            print "Crack %d:"  %  (num+1)
            print "Fracture Args: K1=%.3f, MTS DIRECTION (DEG)=%.3f"  %  (cracks[num]['K1'], cracks[num]['relativeCpd'])
            print "Absolute Angle (DEG):%.3f"  %  cracks[num]['cpd']

            #计算裂尖坐标   
            oldX = oldY = oldZ = newX =newY =0.0
            oldX, oldY, oldZ=getOldCrackTipCoord(assembly=myAssembly, crackNum=num+1)
            newX, newY=calcNewCrackTipCoord(cpd=cracks[num]['cpd'], growth=growth)
            cracks[num]['oldX']=oldX
            cracks[num]['oldY']=oldY
            cracks[num]['oldZ']=oldZ
            cracks[num]['reNewX']=newX
            cracks[num]['reNewY']=newY
            cracks[num]['newX']=oldX+newX
            cracks[num]['newY']=oldY+newY  
            print "oldX=%f, oldY=%f, oldZ=%f" % (cracks[num]['oldX'], cracks[num]['oldY'], cracks[num]['oldZ'])

            #若材料的应力强度因子小于断裂韧度，裂纹不扩展，程序结束
            if cracks[num]['Ke'] > KIC:
                countTerminate=countTerminate+1
                print "Crack %d stop propagation at step %d"  %  (num+1, increaseNum)
                if countTerminate == crackNum: 
                    flag=False
                    print ' '
                    print 'Propagation terminated!'
                    break
                print ' '
                continue
            else:
                print "------------------------------------"
                print "Now begin step %d for crack %d:"  %  (increaseNum+1, num+1)
             
  
            print "RelativeNewX=%f, RelativeNewY=%f" % (cracks[num]['reNewX'], cracks[num]['reNewY'])
            print "newX=%f, newY=%f" % (cracks[num]['newX'], cracks[num]['newY'])
            print "------------------------------------"
            print ' '
            
            #将每条裂纹的关键参数写入文件
            storeCrackData(path=caeFilePath, cracks=cracks, increaseNum=increaseNum, num=num)
    
            deleteSets(assembly=myAssembly, crackNum=num+1)
    
            createSeamPartiton(part=myPart, model=myModel, assembly=myAssembly, partInstance=myPartInstance, centroid=centroid, 
                               oldX=cracks[num]['oldX'], oldY=cracks[num]['oldY'], oldZ=cracks[num]['oldZ'], 
                               reNewX=cracks[num]['reNewX'], reNewY=cracks[num]['reNewY'], 
                               newX=cracks[num]['newX'], newY=cracks[num]['newY'], newZ=cracks[num]['oldZ'], 
                               seamName='seam'+str(num+1)+'-'+str(increaseNum+1))    # 裂缝seam 
        
            createCrackTipPartition(model=myModel, part=myPart, assembly=myAssembly, partInstance=myPartInstance, centroid=centroid, 
                                  radius=growth, cpd=cracks[num]['cpd'], partitionName='outterPartition'+str(num+1), 
                                  oldX=cracks[num]['oldX'], oldY=cracks[num]['oldY'], oldZ=cracks[num]['oldZ'], 
                                  newX=cracks[num]['newX'], newY=cracks[num]['newY'], newZ=cracks[num]['oldZ'])    #外圆分割
            createCrackTipPartition(model=myModel, part=myPart, assembly=myAssembly, partInstance=myPartInstance, centroid=centroid, 
                                  radius=growth*0.2, cpd=cracks[num]['cpd'], partitionName='innerPartition'+str(num+1), 
                                  oldX=cracks[num]['oldX'], oldY=cracks[num]['oldY'], oldZ=cracks[num]['oldZ'], 
                                  newX=cracks[num]['newX'], newY=cracks[num]['newY'], newZ=cracks[num]['oldZ'])    #内圆分割
                                  
            createEdgeSet(part=myPart, assembly=myAssembly, partInstance=myPartInstance, setName='seamIn'+str(num+1), 
                              coord=(cracks[num]['newX'], cracks[num]['newY'], 0.0))
            createEdgeSet(part=myPart, assembly=myAssembly, partInstance=myPartInstance, setName='seamOut'+str(num+1), 
                              coord=(cracks[num]['oldX']+cracks[num]['reNewX']/2, cracks[num]['oldY']+cracks[num]['reNewY']/2, 0.0  ))
            createPtSet(assembly=myAssembly, partInstance=myPartInstance, setName='crackPt'+str(num+1), 
                              coord=(cracks[num]['newX'], cracks[num]['newY'], 0.0) )    
            crackRgSet(assembly=myAssembly, partInstance=myPartInstance, setName='crackRg'+str(num+1), 
                              coord=(cracks[num]['newX'], cracks[num]['newY'], 0.0) )
            crackRgSet(assembly=myAssembly, partInstance=myPartInstance, setName='crackRgOut'+str(num+1), 
                              coord=(cracks[num]['oldX']+cracks[num]['reNewX']*1.5, cracks[num]['oldY']+cracks[num]['reNewY']*1.5, 0.0) )
            crackRgSet(assembly=myAssembly, partInstance=myPartInstance, coord=(0.0, 0.0, 0.0), setName='wholePart')
    
            assignSeam(assembly=myAssembly, partInstance=myPartInstance, crackNum=num+1, increaseNum=increaseNum+1)
            createContourIntegral(assembly=myAssembly, partInstance=myPartInstance, crackName=cracks[num]['name'], 
                              oldX=cracks[num]['oldX'], oldY=cracks[num]['oldY'], oldZ=cracks[num]['oldZ'], 
                              reNewX=cracks[num]['reNewX'], reNewY=cracks[num]['reNewY'], 
                              newX=cracks[num]['newX'], newY=cracks[num]['newY'], newZ=cracks[num]['oldZ'], 
                              crackNum=num+1)
            createHistoryOutput(model=myModel, outputName='crack'+str(num+1), 
                              crackName=cracks[num]['name'],  crackNum=num+1)     #outputName用crackName加上序号代替

            seedEdge(assembly=myAssembly, partInstance=myPartInstance, initCrackLength=cracks[num]['initCrackLength'],
                         circularNum=circularNum, radialNum=radialNum, increaseNum=increaseNum, crackNum=num+1)
            assignMeshControl(assembly=myAssembly, partInstance=myPartInstance, crackNum=num+1)
            assignElementType(assembly=myAssembly, partInstance=myPartInstance, crackNum=num+1)
        
        #对所有裂纹的操作结束，划分网格后提交计算
        #-----------------------------------------------------------------------
    
        if flag == False: break
        
        meshInstance(assembly=myAssembly, partInstance=myPartInstance)
      
        myJob=createJob(modelName=modelName, increaseNum=increaseNum+1)
        submitJob(job=myJob, modelName=modelName, increaseNum=increaseNum+1)
        
        increaseNum+=1
    
    #-----------------------------------------------------------------------

#-----------------------------------------------------------------------
