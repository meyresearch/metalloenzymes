%Chk=vim2_large_mk.chk
%Mem=12000MB
%NProcShared=8
# B3LYP/6-31G* Integral=(Grid=UltraFine) Opt Pop(MK,ReadRadii)
IOp(6/33=2,6/42=6)
 
CLR
 
2  2
H       0  -27.450   14.262    0.186
C      -1  -27.705   13.657   -0.684
H       0  -27.324   14.220   -1.536
H       0  -27.225   12.679   -0.647
C      -1  -29.205   13.484   -0.900
O      -1  -29.597   13.080   -1.999
N      -1  -30.056   13.758    0.090
H       0  -29.847   13.771    0.924
C      -1  -31.461   14.066   -0.204
H       0  -31.506   14.613   -1.004
C      -1  -32.284   12.804   -0.554
H       0  -33.140   13.085   -0.913
H       0  -31.824   12.323   -1.260
C      -1  -32.535   11.852    0.586
N      -1  -31.609   10.914    0.988
H       0  -30.827   10.811    0.646
C      -1  -32.107   10.198    1.980
H       0  -31.668    9.512    2.428
N      -1  -33.335   10.622    2.228
C      -1  -33.626   11.661    1.370
H       0  -34.422   12.141    1.335
C      -1  -32.046   14.863    0.956
O      -1  -31.347   15.201    1.914
N      -1  -33.333   15.193    0.850
H       0  -33.885   14.823    0.304
C      -1  -33.909   16.245    1.687
H       0  -33.093   16.948    1.854
H       0  -34.742   16.734    1.182
C      -1  -34.383   15.762    3.049
O      -1  -34.830   16.594    3.852
N      -1  -34.311   14.461    3.324
H       0  -33.960   13.876    2.800
C      -1  -34.837   13.919    4.568
H       0  -35.699   14.337    4.722
C      -1  -35.034   12.400    4.467
H       0  -34.238   12.003    4.080
H       0  -35.128   12.034    5.360
C      -1  -36.223   12.003    3.650
N      -1  -36.387   10.736    3.136
C      -1  -37.525   10.682    2.467
H       0  -37.859    9.937    2.022
N      -1  -38.108   11.864    2.536
H       0  -38.862   12.071    2.177
C      -1  -37.316   12.708    3.273
H       0  -37.496   13.598    3.476
C      -1  -33.921   14.246    5.749
O      -1  -32.756   14.622    5.599
N      -1  -34.483   14.082    6.944
H       0  -35.253   13.716    7.060
C      -1  -33.837   14.513    8.183
H       0  -33.538   15.556    8.078
H       0  -34.556   14.402    8.995
C      -1  -32.544   13.748    8.477
O      -1  -31.664   14.272    9.170
N      -1  -32.398   12.503    7.991
H       0  -32.979   12.063    7.534
C      -1  -31.132   11.827    8.263
H       0  -30.875   12.065    9.168
C      -1  -31.296   10.294    8.197
H       0  -30.460    9.877    8.459
H       0  -31.964   10.021    8.845
C      -1  -31.700    9.771    6.818
O      -1  -32.186   10.545    5.963
O      -1  -31.563    8.537    6.617
C      -1  -30.006   12.304    7.353
O      -1  -28.862   11.859    7.510
N      -1  -30.291   13.228    6.438
H       0  -31.083   13.536    6.305
C      -1  -29.279   13.807    5.572
H       0  -28.351   13.301    5.840
H       0  -29.536   13.647    4.525
H       0  -29.134   14.877    5.718
H       0  -35.306    7.105   -6.658
C      -1  -35.247    6.120   -6.196
H       0  -35.800    5.396   -6.795
H       0  -34.207    5.798   -6.154
C      -1  -35.884    6.155   -4.815
O      -1  -37.114    6.184   -4.677
N      -1  -35.047    6.146   -3.785
H       0  -34.196    6.038   -3.849
C      -1  -35.542    6.326   -2.429
H       0  -36.168    5.620   -2.204
C      -1  -34.353    6.265   -1.473
H       0  -34.049    5.346   -1.409
H       0  -33.622    6.780   -1.849
C      -1  -34.644    6.773   -0.100
N      -1  -35.184    5.976    0.889
H       0  -35.395    5.147    0.801
C      -1  -35.326    6.689    1.993
H       0  -35.669    6.380    2.800
N      -1  -34.896    7.918    1.754
C      -1  -34.466    7.995    0.452
H       0  -34.118    8.746    0.029
C      -1  -36.302    7.641   -2.302
O      -1  -37.300    7.727   -1.576
N      -1  -35.842    8.670   -3.010
H       0  -35.100    8.637   -3.444
C      -1  -36.519    9.951   -3.137
H       0  -37.590    9.832   -2.972
H       0  -36.177   10.657   -2.380
H       0  -36.340   10.319   -4.147
H       0  -26.198    6.538   -0.147
C      -1  -26.991    7.011    0.431
H       0  -26.615    7.290    1.415
H       0  -27.308    7.930   -0.061
C      -1  -28.178    6.070    0.580
O      -1  -28.076    4.846    0.438
N      -1  -29.336    6.667    0.836
H       0  -29.468    7.517    0.819
C      -1  -30.502    5.860    1.176
H       0  -30.223    5.130    1.750
C      -1  -31.498    6.751    1.919
H       0  -32.293    6.250    2.161
H       0  -31.787    7.487    1.357
S      -1  -30.645    7.372    3.409
C      -1  -31.128    5.176   -0.032
O      -1  -32.062    4.372    0.134
N      -1  -30.623    5.453   -1.233
H       0  -30.035    6.064   -1.375
C      -1  -31.034    4.729   -2.426
H       0  -32.121    4.646   -2.402
H       0  -30.718    5.282   -3.311
H       0  -30.646    3.711   -2.457
H       0  -27.003    4.074    3.039
C      -1  -27.315    4.842    3.747
H       0  -26.628    5.686    3.684
H       0  -28.300    5.211    3.462
C      -1  -27.365    4.349    5.177
O      -1  -27.155    5.125    6.111
N      -1  -27.638    3.058    5.375
H       0  -27.722    2.486    4.738
C      -1  -27.820    2.496    6.708
H       0  -27.372    3.057    7.360
C      -1  -29.305    2.448    7.092
H       0  -29.764    1.820    6.512
H       0  -29.389    2.112    7.998
C      -1  -29.974    3.781    6.999
N      -1  -30.022    4.653    8.061
H       0  -29.690    4.507    8.841
C      -1  -30.660    5.754    7.694
H       0  -30.825    6.492    8.235
N      -1  -31.022    5.626    6.429
C      -1  -30.605    4.400    5.970
H       0  -30.732    4.059    5.114
C      -1  -27.228    1.104    6.717
O      -1  -27.548    0.293    5.840
N      -1  -26.362    0.833    7.695
H       0  -26.065    1.409    8.261
C      -1  -25.852   -0.515    7.861
H       0  -25.755   -0.716    8.928
H       0  -26.593   -1.213    7.471
H       0  -24.896   -0.729    7.384
Zn     -1  -34.834    9.422    3.139
Zn     -1  -31.821    7.098    5.287
C      -1  -34.220    3.393    8.083
C      -1  -33.874    4.740    8.102
C      -1  -34.884    5.693    8.127
C      -1  -34.646    7.058    7.483
C      -1  -36.192    5.293    8.463
C      -1  -37.204    6.225    8.768
C      -1  -38.480    5.782    9.094
C      -1  -38.759    4.419    9.121
C      -1  -37.763    3.500    8.821
C      -1  -36.469    3.947    8.490
F      -1  -36.936    7.562    8.743
F      -1  -40.013    3.986    9.441
O      -1  -33.376    2.569    7.963
O      -1  -34.636    8.347    5.043
O      -1  -33.470    6.176    5.147
O      -1  -35.934    6.255    5.179
O      -1  -35.508    3.020    8.200
P      -1  -34.673    6.954    5.625
H       0  -35.360    7.690    7.785
H       0  -33.756    7.407    7.777
H       0  -32.913    5.018    8.098
H       0  -39.199    6.443    9.309
H       0  -37.964    2.521    8.840
H       0  -36.745    6.820    5.366
 
Zn 1.395
 
Zn 1.395
 
 
