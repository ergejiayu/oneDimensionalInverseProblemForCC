def long_function_oneDensonaldiff(h,var_deltTime,MiddleTemp,var_XNumber,var_X,var_temperatureWater,var_rouS,var_rouL,var_specificHeatS,var_specificHeatL,var_TconductivityS,var_TconductivityL,var_liqTemp,var_SodTemp,var_m):
    NextTemp=[0]*var_XNumber
    deltX=var_X/var_XNumber # 铸坯在厚度方向的空间间隔
    Boltzman = 0.000000056684 # 玻尔兹曼常数
    emissivity=0.8 # 辐射系数
    #### equation strat #####
    for i in range(var_XNumber):
        if MiddleTemp[i]>=var_liqTemp:
            rou=var_rouL
            specificHeat=var_specificHeatL
            Tconductivity=var_TconductivityL
        if  MiddleTemp[i]<=var_SodTemp:
            rou=var_rouS
            specificHeat=var_specificHeatS
            Tconductivity=var_TconductivityS
        if  (var_SodTemp<MiddleTemp[i]) & (var_liqTemp>MiddleTemp[i]):
            rou=(var_rouS-var_rouL)*(var_liqTemp-MiddleTemp[i])/(var_liqTemp-var_SodTemp)+var_rouL
            Tconductivity=(var_TconductivityS)*(var_liqTemp-MiddleTemp[i])/(var_liqTemp-var_SodTemp)+var_m*(1-(var_liqTemp-MiddleTemp[i])/(var_liqTemp-var_SodTemp))*var_TconductivityL
            specificHeat=(var_specificHeatS-var_specificHeatL)*(var_liqTemp-MiddleTemp[i])/(var_liqTemp-var_SodTemp)+var_specificHeatL+274950/(var_liqTemp-var_SodTemp)
        a1=Tconductivity/(rou*specificHeat)
        a= (Tconductivity*var_deltTime)/(rou*specificHeat*deltX*deltX)
        if i==0:
            NextTemp[i]=(1-2*var_deltTime*h/(rou*specificHeat*deltX)-2*var_deltTime*a1/(deltX*deltX))*MiddleTemp[i]+2*var_deltTime*a1*MiddleTemp[i+1]/(deltX*deltX)+2*var_deltTime*h*var_temperatureWater/(rou*specificHeat*deltX)
        if i==var_XNumber-1:
            NextTemp[i]=(1-2*var_deltTime*h/(rou*specificHeat*deltX)-2*var_deltTime*a1/(deltX*deltX))*MiddleTemp[i]+2*var_deltTime*a1*MiddleTemp[i-1]/(deltX*deltX)+2*var_deltTime*h*var_temperatureWater/(rou*specificHeat*deltX)
        if 0<i<var_XNumber-1:
            NextTemp[i]=a*MiddleTemp[i+1]+a*MiddleTemp[i-1]+(1-2*a)*MiddleTemp[i]
    #### equation end #####
    return NextTemp
        
