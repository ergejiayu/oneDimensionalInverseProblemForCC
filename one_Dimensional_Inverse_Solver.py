def long_function_oneDimensionalInverseSolver(var_XNumber,var_X,var_temperatureWater,var_rouS,var_rouL,var_specificHeatS,var_specificHeatL,var_TconductivityS,var_TconductivityL,var_liqTemp,var_SodTemp,var_m,var_castingTemp,var_dis,var_measuredTemp,var_disMeasure,var_numberofMeasuredPoint,var_SCZ_num,var_VcastOriginal,var_deltTime,var_iter_max,var_eps,var_h_initial)
    from one_Dinensional_Diff_Solver import long_function_oneDensonaldiff
    import numpy.linalg as nlg
    import numpy as np
    import math
    h_end=[] # 要求的换热系数
    no_MeasuredPoint=var_SCZ_num-var_numberofMeasuredPoint # 为测量点的个数
    h_new=[0]*var_SCZ_num # 求换热系数时的灵敏度系数
    dh=[0]*var_numberofMeasuredPoint # 换热系数增量的初值
    tl=[0]*10;  # 铸坯凝固时间的初值
    for i in range(10):
        tl[i]=var_dis[i]/var_VcastOriginal
    time_Mold = int((tl[1]-tl[0])/var_deltTime)
    time_SCZ = int((tl[9]-tl[1])/var_deltTime)
    Time_all = time_Mold+time_SCZ
    TemSur = [0]*(Time_all) # 表面温度的初值
    measuredStepPoint=[0]*var_numberofMeasuredPoint # 网格划分以后测量温度在网格中的位置
    calTemp=[0]*var_numberofMeasuredPoint # 计算温度初始化
    calTemp_sensive=[0]*var_numberofMeasuredPoint # 需计算温度的灵敏度系数
    errorTemp=[0]*var_numberofMeasuredPoint # 表面温度计算值与测量值之间的偏差
    initialTemp=[1]*var_XNumber # 整个铸坯温度场的初始温度
    J=[0]*var_iter_max # 迭代过程中铸坯表面温度的测量值与计算值之间的平方差的值
    eye=[([0]*var_numberofMeasuredPoint) for i in range(var_numberofMeasuredPoint)] # 单位矩阵的初始化
    delth=1 # 换热系数的增量
    dudh_1=[([0]*var_numberofMeasuredPoint) for i in range(var_numberofMeasuredPoint)] # 当前换热系数得到的表面温度组成的矩阵初始化
    dudh_2=[([0]*var_numberofMeasuredPoint) for i in range(var_numberofMeasuredPoint)] # 换热系数加入增量以后得到的表面温度组成的矩阵初始化
    dudh=[([0]*var_numberofMeasuredPoint) for i in range(var_numberofMeasuredPoint)] # 灵敏度系数矩阵初始化
    dudh_transport=[([0]*var_numberofMeasuredPoint) for i in range(var_numberofMeasuredPoint)] # 灵敏度系数矩阵转置初始化
    h_delt_all=[([0]*var_SCZ_num) for i in range(var_SCZ_num)] # 换热系数加入增量以后得到的矩阵
    for i in range(var_XNumber):
        initialTemp[i]=var_castingTemp
    for i in range(var_numberofMeasuredPoint):
        measuredStepPoint[i]=int((var_disMeasure[i]/var_VcastOriginal)/var_deltTime)
    MiddleTemp=initialTemp # 当前时刻铸坯界面的温度分布
    NextTemp=initialTemp # 下一时刻铸坯界面的温度分布
    for stepTime in range(1,Time_all):
        if stepTime <= time_Mold:
            tTime = var_deltTime*stepTime;
            h=1000*(0.07128*math.exp(-tTime)+2.328*math.exp(-tTime/9.5)+0.698)
            NextTemp=long_function_oneDensonaldiff(h,var_deltTime,MiddleTemp,var_XNumber,var_X,var_temperatureWater,var_rouS,var_rouL,var_specificHeatS,var_specificHeatL,var_TconductivityS,var_TconductivityL,var_liqTemp,var_SodTemp,var_m)
        else:
            disNow = var_dis[1]+stepTime*var_VcastOriginal*var_deltTime
            if var_dis[1]<=disNow<=var_dis[2]:
                h=var_h_initial[0]          
            if var_dis[2]<disNow<=var_dis[3]:
                h=var_h_initial[1]
            if var_dis[3]<disNow<=var_dis[4]:
                h=var_h_initial[2]
            if var_dis[4]<disNow<=var_dis[5]:
                h=var_h_initial[3]
            if var_dis[5]<disNow<=var_dis[6]:
                h=var_h_initial[4]
            if var_dis[6]<disNow<=var_dis[7]:
                h=var_h_initial[5]
            if var_dis[7]<disNow<=var_dis[8]:
                h=var_h_initial[6]
            if var_dis[8]<disNow:
                h=var_h_initial[7]
            NextTemp=long_function_oneDensonaldiff(h,var_deltTime,MiddleTemp,var_XNumber,var_X,var_temperatureWater,var_rouS,var_rouL,var_specificHeatS,var_specificHeatL,var_TconductivityS,var_TconductivityL,var_liqTemp,var_SodTemp,var_m)
        MiddleTemp=NextTemp
        TemSur[stepTime-1]=MiddleTemp[0]
        for j in range(var_numberofMeasuredPoint):
            calTemp[j]=TemSur[measuredStepPoint[j]]
    for i in range(var_numberofMeasuredPoint):
        for j in range(var_numberofMeasuredPoint):
            if i==j:
                eye[i][j]=1.0
    for iter_num in range(var_iter_max):
        ####  model start ####
        for i in range(var_numberofMeasuredPoint):
            J[iter_num]=J[iter_num]+(calTemp[i]-var_measuredTemp[i])*(calTemp[i]-var_measuredTemp[i])
        #### model end #### where calTemp is calculated by equation
        #print(J[iter_num])
        if J[iter_num]<var_eps:
            break
        J[iter_num]=0
        #### optimation start ############
        for i in range(var_numberofMeasuredPoint):
            for j in range(var_numberofMeasuredPoint):
                dudh_1[i][j]=calTemp[j]
        for i in range(var_SCZ_num):
            for j in range(var_SCZ_num):
                if i==j:
                    h_delt_all[i][j]=var_h_initial[j]+delth
                else:
                    h_delt_all[i][j]=var_h_initial[j]
        for i in range(no_MeasuredPoint):
            h_delt_all[i][i]=h_delt_all[i][i]-delth
        for i in range(no_MeasuredPoint,var_SCZ_num):     
            for j in range(var_SCZ_num):
                h_new[j]=h_delt_all[i][j]
            MiddleTemp=initialTemp
            NextTemp=initialTemp
            for stepTime in range(1,Time_all):
                if stepTime <= time_Mold:
                    tTime = var_deltTime*stepTime;
                    h=1000*(0.07128*math.exp(-tTime)+2.328*math.exp(-tTime/9.5)+0.698)
                    NextTemp=long_function_oneDensonaldiff(h,var_deltTime,MiddleTemp,var_XNumber,var_X,var_temperatureWater,var_rouS,var_rouL,var_specificHeatS,var_specificHeatL,var_TconductivityS,var_TconductivityL,var_liqTemp,var_SodTemp,var_m)
                else:
                    disNow = var_dis[1]+stepTime*var_VcastOriginal*var_deltTime
                    if var_dis[1]<=disNow<=var_dis[2]:
                        h=h_new[0]          
                    if var_dis[2]<disNow<=var_dis[3]:
                        h=h_new[1]
                    if var_dis[3]<disNow<=var_dis[4]:
                        h=h_new[2]
                    if var_dis[4]<disNow<=var_dis[5]:
                        h=h_new[3]
                    if var_dis[5]<disNow<=var_dis[6]:
                        h=h_new[4]
                    if var_dis[6]<disNow<=var_dis[7]:
                        h=h_new[5]
                    if var_dis[7]<disNow<=var_dis[8]:
                        h=h_new[6]
                    if var_dis[8]<disNow:
                        h=h_new[7]
                    NextTemp=long_function_oneDensonaldiff(h,var_deltTime,MiddleTemp,var_XNumber,var_X,var_temperatureWater,var_rouS,var_rouL,var_specificHeatS,var_specificHeatL,var_TconductivityS,var_TconductivityL,var_liqTemp,var_SodTemp,var_m)
                MiddleTemp=NextTemp
                TemSur[stepTime-1]=MiddleTemp[0]
            for j in range(var_numberofMeasuredPoint):
                calTemp_sensive[j]=TemSur[measuredStepPoint[j]]
            for k in range(var_numberofMeasuredPoint):
                dudh_2[i-no_MeasuredPoint][k]=calTemp_sensive[k]
        for i in range(var_numberofMeasuredPoint):
            for j in range(var_numberofMeasuredPoint):
                dudh[i][j]=(dudh_1[i][j]-dudh_2[i][j])/delth
        for i in range(var_numberofMeasuredPoint):
            errorTemp[i]=(calTemp[i]-var_measuredTemp[i])
        dudh_transport=nlg.inv(eye+np.transpose(dudh)*dudh)*dudh
        for i in range(var_numberofMeasuredPoint):
            dh[i]=0
            for j in range(var_numberofMeasuredPoint):
                dh[i]=dh[i]+dudh_transport[i][j]*errorTemp[j]
        for r in range(no_MeasuredPoint,var_SCZ_num):
            var_h_initial[r]=var_h_initial[r]+dh[r-no_MeasuredPoint]
        MiddleTemp=initialTemp
        NextTemp=initialTemp
        for stepTime in range(1,Time_all):
            if stepTime <= time_Mold:
                tTime = var_deltTime*stepTime;
                h=1000*(0.07128*math.exp(-tTime)+2.328*math.exp(-tTime/9.5)+0.698)
                NextTemp=long_function_oneDensonaldiff(h,var_deltTime,MiddleTemp,var_XNumber,var_X,var_temperatureWater,var_rouS,var_rouL,var_specificHeatS,var_specificHeatL,var_TconductivityS,var_TconductivityL,var_liqTemp,var_SodTemp,var_m)
            else:
                disNow = var_dis[1]+stepTime*var_VcastOriginal*var_deltTime
                if var_dis[1]<=disNow<=var_dis[2]:
                    h=var_h_initial[0]          
                if var_dis[2]<disNow<=var_dis[3]:
                    h=var_h_initial[1]
                if var_dis[3]<disNow<=var_dis[4]:
                    h=var_h_initial[2]
                if var_dis[4]<disNow<=var_dis[5]:
                    h=var_h_initial[3]
                if var_dis[5]<disNow<=var_dis[6]:
                    h=var_h_initial[4]
                if var_dis[6]<disNow<=var_dis[7]:
                    h=var_h_initial[5]
                if var_dis[7]<disNow<=var_dis[8]:
                    h=var_h_initial[6]
                if var_dis[8]<disNow:
                    h=var_h_initial[7]
                NextTemp=long_function_oneDensonaldiff(h,var_deltTime,MiddleTemp,var_XNumber,var_X,var_temperatureWater,var_rouS,var_rouL,var_specificHeatS,var_specificHeatL,var_TconductivityS,var_TconductivityL,var_liqTemp,var_SodTemp,var_m)
            MiddleTemp=NextTemp
            TemSur[stepTime-1]=MiddleTemp[0]
            for j in range(var_numberofMeasuredPoint):
                calTemp[j]=TemSur[measuredStepPoint[j]]
        #### optimation end ############
    h_end = var_h_initial
    return h_end
   
