#Tungsten Ionization Function
def ionization(temperature):
    Z_W = 1
    temperature = temperature * 1000
    if temperature < 16.37:
        Z_W = 1
    if 16.37 < temperature < 26.0:
        Z_W = 2
    if 26.0 < temperature < 38.2:
        Z_W = 3
    if 38.2 < temperature < 51.6:
        Z_W = 4
    if 51.6 < temperature < 64.77:
        Z_W = 5
    if 64.77 < temperature < 122.01:
        Z_W = 6
    if 122.01 < temperature < 141.2:
        Z_W = 7
    if 141.2 < temperature < 160.2:
        Z_W = 8
    if 160.2 < temperature < 179.0:
        Z_W = 9
    if 179.0 < temperature < 208.9:
        Z_W = 10
    if 208.9 < temperature < 231.6:
        Z_W = 11
    if 231.6 < temperature < 258.3:
        Z_W = 12
    if 258.3 < temperature < 290.7:
        Z_W = 13
    if 290.7 < temperature < 325.3:
        Z_W = 14
    if 325.3 < temperature < 361.9:
        Z_W = 15
    if 361.9 < temperature < 387.9:
        Z_W = 16
    if 387.9 < temperature < 420.7:
        Z_W = 17
    if 420.7 < temperature < 462.1:
        Z_W = 18
    if 462.1 < temperature < 502.6:
        Z_W = 19
    if 502.6 < temperature < 543.4:
        Z_W = 20
    if 543.4 < temperature < 594.5:
        Z_W = 21
    if 594.5 < temperature < 640.6:
        Z_W = 22
    if 640.6 < temperature < 685.6:
        Z_W = 23
    if 685.6 < temperature < 734.1:
        Z_W = 24
    if 734.1 < temperature < 784.4:
        Z_W = 25
    if 784.4 < temperature < 833.4:
        Z_W = 26
    if 833.4 < temperature < 881.4:
        Z_W = 27
    if 881.4 < temperature < 1132.2:
        Z_W = 28
    if 1132.2 < temperature < 1180.0:
        Z_W = 29
    if 1180.0 < temperature < 1230.4:
        Z_W = 30
    if 1230.4 < temperature < 1283.4:
        Z_W = 31
    if 1283.4 < temperature < 1335.1:
        Z_W = 32
    if 1335.1 < temperature < 1386.8:
        Z_W = 33
    if 1386.8 < temperature < 1459.9:
        Z_W = 34
    if 1459.9 < temperature < 1512.4:
        Z_W = 35
    if 1512.4 < temperature < 1569.1:
        Z_W = 36
    if 1569.1 < temperature < 1621.7:
        Z_W = 37
    if 1621.7 < temperature < 1829.8:
        Z_W = 38
    if 1829.8 < temperature < 1882.9:
        Z_W = 39
    if 1882.9 < temperature < 1940.6:
        Z_W = 40
    if 1940.6 < temperature < 1994.8:
        Z_W = 41
    if 1994.8 < temperature < 2149.1:
        Z_W = 42
    if 2149.1 < temperature < 2210.0:
        Z_W = 43
    if 2210.0 < temperature < 2354.5:
        Z_W = 44
    if 2354.5 < temperature < 2414.1:
        Z_W = 45
    if 2414.1 < temperature < 4057:
        Z_W = 46
    if 4057 < temperature < 4180:
        Z_W = 47
    if 4180 < temperature < 4309:
        Z_W = 48
    if 4309 < temperature < 4446:
        Z_W = 49
    if 4446 < temperature < 4578:
        Z_W = 50
    if 4578 < temperature < 4709:
        Z_W = 51
    if 4709 < temperature < 4927:
        Z_W = 52
    if 4927 < temperature < 5063:
        Z_W = 53
    if 5063 < temperature < 5209:
        Z_W = 54
    if 5209 < temperature < 5348:
        Z_W = 55
    if 5348 < temperature < 5719:
        Z_W = 56
    if 5719 < temperature < 5845:
        Z_W = 57
    if 5845 < temperature < 5970:
        Z_W = 58
    if 5970 < temperature < 6093:
        Z_W = 59
    if 6093 < temperature < 6596:
        Z_W = 60
    if 6596 < temperature < 6735:
        Z_W = 61
    if 6735 < temperature < 7000:
        Z_W = 62
    if 7000 < temperature < 7130:
        Z_W = 63
    if 7130 < temperature < 15566:
        Z_W = 64
    if 15566 < temperature < 15896:
        Z_W = 65
    if 15896 < temperature < 16252:
        Z_W = 66
    if 16252 < temperature < 16588:
        Z_W = 67
    if 16588 < temperature < 18476:
        Z_W = 68
    if 18476 < temperature < 18872:
        Z_W = 69
    if 18872 < temperature < 19362:
        Z_W = 70
    if 19362 < temperature < 19686.74:
        Z_W = 71
    if 19686.74 < temperature < 79181.94:
        Z_W = 72
    if 79181.94 < temperature < 80755.6:
        Z_W = 73
    return Z_W
