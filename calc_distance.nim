import json
import os,strutils,tables,times
import Tables
import math
import algorithm

let start = cputime()

type Vincenty = ref object
    lat1 : float
    lon1 : float
    lat2 : float
    lon2 : float
    ellipsoid : int
    ELLIPSOID_GRS80 : int
    ELLIPSOID_WGS84 : int
    GEODETIC_DATUM :  Table[int, seq[float]]
    ITERATION_LIMIT : int

proc vincenty_inverse(calculator:Vincenty): float = 

    # 差異が無ければ0.0を返す
    if calculator.lat1 == calculator.lat2 and calculator.lon1 == calculator.lon2:
        return 0.0
    # 計算時に必要な長軸半径(a)と扁平率(ƒ)を定数から取得し、短軸半径(b)を算出する
    # 楕円体が未指定の場合はGRS80の値を用いる
    var a = calculator.GEODETIC_DATUM[calculator.ELLIPSOID_GRS80][0]
    var ƒ = calculator.GEODETIC_DATUM[calculator.ELLIPSOID_GRS80][1]
    var b = (1 - ƒ) * a

    var φ1 = degToRad(calculator.lat1)
    var φ2 = degToRad(calculator.lat2)
    var λ1 = degToRad(calculator.lon1)
    var λ2 = degToRad(calculator.lon2)

    # 更成緯度(補助球上の緯度)
    var U1 = arctan((1 - ƒ) * tan(φ1))
    var U2 = arctan((1 - ƒ) * tan(φ2))

    var sinU1 = sin(U1)
    var sinU2 = sin(U2)
    var cosU1 = cos(U1)
    var cosU2 = cos(U2)

    # 2点間の経度差
    var L = λ2 - λ1

    # λをLで初期化
    var λ = L

    var sinλ:float
    var cosλ:float
    var sinσ:float
    var cosσ:float
    var σ:float
    var sinα:float
    var cos2α:float
    var cos2σm:float
    var C:float
    var λʹ:float
    # 以下の計算をλが収束するまで反復する
    # 地点によっては収束しないことがあり得るため、反復回数に上限を設ける
    for count in 0..calculator.ITERATION_LIMIT:
        sinλ = sin(λ)
        cosλ = cos(λ)
        sinσ = sqrt((cosU2 * sinλ) ^ 2 + (cosU1 * sinU2 - sinU1 * cosU2 * cosλ) ^ 2)
        cosσ = sinU1 * sinU2 + cosU1 * cosU2 * cosλ
        σ = arctan2(sinσ, cosσ)
        sinα = cosU1 * cosU2 * sinλ / sinσ
        cos2α = 1 - sinα ^ 2
        cos2σm = cosσ - 2 * sinU1 * sinU2 / cos2α
        C = ƒ / 16 * cos2α * (4 + ƒ * (4 - 3 * cos2α))
        λʹ = λ
        λ = L + (1 - C) * ƒ * sinα * (σ + C * sinσ * (cos2σm + C * cosσ * (-1 + 2 * cos2σm ^ 2)))

        # 偏差が.000000000001以下ならbreak
        if abs(λ - λʹ) <= 1e-12:
            break

    # λが所望の精度まで収束したら以下の計算を行う
    var u2 = cos2α * (a ^ 2 - b ^ 2) / (b ^ 2)
    var A = 1 + u2 / 16384 * (4096 + u2 * (-768 + u2 * (320 - 175 * u2)))
    var B = u2 / 1024 * (256 + u2 * (-128 + u2 * (74 - 47 * u2)))
    var Δσ = B * sinσ * (cos2σm + B / 4 * (cosσ * (-1 + 2 * cos2σm ^ 2) - B / 6 * cos2σm * (-3 + 4 * sinσ ^ 2) * (-3 + 4 * cos2σm ^ 2)))

    # 2点間の楕円体上の距離
    var s = b * A * (σ - Δσ)

    return s


var calculator : Vincenty  = new Vincenty

calculator.lat1 = 0.0
calculator.lon1 = 0.0
calculator.lat2 = 0.0
calculator.lon2 = 0.0

calculator.ellipsoid = 1
calculator.ELLIPSOID_GRS80 = 1
calculator.ELLIPSOID_WGS84 = 2
calculator.GEODETIC_DATUM[calculator.ELLIPSOID_GRS80] = @[6378137.0, 1 / 298.257222101]
calculator.GEODETIC_DATUM[calculator.ELLIPSOID_WGS84] = @[6378137.0, 1 / 298.257223563]

calculator.ITERATION_LIMIT = 1000

type TempType = tuple
    distance: float
    address: string

type CoordinateType = tuple
    address: string
    distance: float

let jsonObj = parseFile("tokyoCoordinates.json")

for key1, value1 in jsonObj:
    let pre_start = cputime()
    var temp_dict: seq[CoordinateType]
    
    for key2, value2 in jsonObj:

        calculator.lat1 = value1[0].getFloat()
        calculator.lon1 = value1[1].getFloat()
        calculator.lat2 = value2[0].getFloat()
        calculator.lon2 = value2[1].getFloat()

        if key1 == key2:
            discard
        elif round(calculator.vincenty_inverse(), 3) == 0.0:
            discard
        else:
            temp_dict.add((key2,round(calculator.vincenty_inverse(), 3)))

    var ranking : seq[TempType]
    for i,kv in temp_dict:
        ranking.add((kv[1],kv[0]))
    ranking.sort

    var newkey = key1 & " との距離が近いランキング TOP10"
    
    var coordinate : CoordinateType
    
    block addTop10:
        var f : File = open("calc_distance.txt" ,fmAppend)
        for i,kv in ranking:
            if i < 10:
                coordinate = (kv[1],kv[0])
                f.writeLine(newkey," : ",coordinate," m")
            else:
                defer: f.close()
                break addTop10
                
            
    echo key1, " . ", round(cputime() - pre_start, 3), " sec"
echo "Finished. ", round(cputime() - start, 3), " sec"