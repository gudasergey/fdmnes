! FDMNES subroutines
! Package with different tables

!***********************************************************************

! Table with edge energies and shift to make the correspondance with the Fermi energy

subroutine esdata(Eseuil,icheck,jseuil,nbseuil,nseuil,numat,mpirank) 

  use declarations
  implicit none

  integer, parameter:: nm1 = 18
  integer, parameter:: nm4 = 30
  integer, parameter:: nn1 = 36
  integer, parameter:: nn4 = 48
  integer, parameter:: nn6 = 58
  integer, parameter:: no1 = 54
  integer, parameter:: no4 = 80
  integer, parameter:: np1 = 86
  integer, parameter:: np2 = 87

  integer:: icheck, ipr, jseuil, mpirank, nbseuil, nseuil, numat, Zm
  
  real(kind=db):: Shift, Shift_edge
  real(kind=db), dimension(nbseuil):: Eseuil
  real(kind=db), dimension(Z_Mendeleiev_max):: ek1, el1, el2, el3
  real(kind=db), dimension(nm1:Z_Mendeleiev_max):: em1, em2, em3
  real(kind=db), dimension(nm4:Z_Mendeleiev_max):: em4, em5
  real(kind=db), dimension(nn1:Z_Mendeleiev_max):: en1, en2, en3
  real(kind=db), dimension(nn4:Z_Mendeleiev_max):: en4, en5
  real(kind=db), dimension(nn6:Z_Mendeleiev_max):: en6, en7
  real(kind=db), dimension(no1:Z_Mendeleiev_max):: eo1, eo2, eo3
  real(kind=db), dimension(no4:Z_Mendeleiev_max):: eo4, eo5
  real(kind=db), dimension(np1:Z_Mendeleiev_max):: ep1
  real(kind=db), dimension(np2:Z_Mendeleiev_max):: ep2, ep3

! Compile par Gwyn Williams,
! http://xray.uu.se/hypertext/EBindEnergies.html
! Pu K: Canada et al, Nucl. Instr. Meth. (1973)
! Actinides : http://www.ruppweb.org/Xray/elements.html

  data ek1/     13.6,    24.6,    54.7,   111.5,   188.0,   284.2,   409.9,   543.1,   696.7,   870.2, &
              1070.8,  1303.0,  1559.0,  1839.0,  2145.5,  2472.0,  2822.4,  3205.9,  3608.4,  4038.5, &
              4492.0,  4966.0,  5465.0,  5989.0,  6539.0,  7112.0,  7709.0,  8333.0,  8979.0,  9659.0, &
             10367.0, 11103.0, 11867.0, 12658.0, 13474.0, 14326.0, 15200.0, 16105.0, 17038.0, 17998.0, &
             18986.0, 20000.0, 21044.0, 22117.0, 23220.0, 24350.0, 25514.0, 26711.0, 27940.0, 29200.0, &
             30491.0, 31814.0, 33169.0, 34561.0, 35985.0, 37441.0, 38925.0, 40443.0, 41991.0, 43569.0, &
             45184.0, 46834.0, 48519.0, 50239.0, 51996.0, 53789.0, 55618.0, 57486.0, 59390.0, 61332.0, &
             63314.0, 65351.0, 67416.0, 69525.0, 71676.0, 73871.0, 76111.0, 78395.0, 80725.0, 83102.0, &
             85530.0, 88005.0, 90526.0, 93105.0, 95730.0, 98404.0,101137.0,103922.0,106755.0,109651.0, &
            112601.0,115606.0,118678.0,121795.0,125027.0,128200.0,131590.0,135960.0,139490.0,143090.0, &
            146780.0,150540.0,154380.0/

  data el1/      0.0,     0.0,     0.0,     0.0,     0.0,     0.0,    37.3,    41.6,     0.0,    48.5, &
                63.5,    88.6,   117.8,   149.7,   189.0,   230.9,   270.0,   326.3,   378.6,   438.4, &
               498.0,   560.9,   626.7,   696.0,   769.1,   844.6,   925.1,  1008.6,  1096.7,  1196.2, &
              1299.0,  1414.6,  1527.0,  1652.0,  1782.0,  1921.0,  2065.0,  2216.0,  2373.0,  2532.0, &
              2698.0,  2866.0,  3043.0,  3224.0,  3412.0,  3604.0,  3806.0,  4018.0,  4238.0,  4465.0, &
              4698.0,  4939.0,  5188.0,  5453.0,  5714.0,  5989.0,  6266.0,  6548.0,  6835.0,  7126.0, &
              7428.0,  7737.0,  8052.0,  8376.0,  8708.0,  9046.0,  9394.0,  9751.0, 10116.0, 10486.0, &
             10870.0, 11271.0, 11682.0, 12100.0, 12527.0, 12968.0, 13419.0, 13880.0, 14353.0, 14839.0, &
             15347.0, 15861.0, 16388.0, 16939.0, 17493.0, 18049.0, 18639.0, 19237.0, 19840.0, 20472.0, &
             21104.6, 21757.4, 22426.8, 23097.2, 23772.9, 24460.0, 25275.0, 26110.0, 26900.0, 27700.0, &
             28530.0, 29380.0, 30240.0/

  data el2/      0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,    21.7, &
                30.4,    49.6,    72.9,    99.8,   136.0,   163.6,   202.0,   250.6,   297.3,   349.7, &
               403.6,   460.2,   519.8,   583.8,   649.9,   719.9,   793.2,   870.0,   952.3,  1044.9, &
              1143.2,  1248.1,  1359.1,  1474.3,  1596.0,  1730.9,  1864.0,  2007.0,  2156.0,  2307.0, &
              2465.0,  2625.0,  2793.0,  2967.0,  3146.0,  3330.0,  3524.0,  3727.0,  3938.0,  4156.0, &
              4380.0,  4612.0,  4852.0,  5107.0,  5359.0,  5624.0,  5891.0,  6164.0,  6440.0,  6722.0, &
              7013.0,  7312.0,  7617.0,  7930.0,  8252.0,  8581.0,  8918.0,  9264.0,  9617.0,  9978.0, &
             10349.0, 10739.0, 11136.0, 11544.0, 11959.0, 12385.0, 12824.0, 13273.0, 13734.0, 14209.0, &
             14698.0, 15200.0, 15711.0, 16244.0, 16785.0, 17337.0, 17907.0, 18484.0, 19083.0, 19693.0, &
             20313.7, 20947.6, 21600.5, 22266.2, 22944.0, 23779.0, 24385.0, 25250.0, 26020.0, 26810.0, &
             27610.0, 28440.0, 29280.0/

  data el3/      0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,    21.6, &
                30.5,    49.2,    72.5,    99.2,   135.0,   162.5,   200.0,   248.4,   294.6,   346.2, &
               398.7,   453.8,   512.1,   574.1,   638.7,   706.8,   778.1,   852.7,   932.7,  1021.8, &
              1116.4,  1217.0,  1323.6,  1433.9,  1550.0,  1678.4,  1804.0,  1940.0,  2080.0,  2223.0, &
              2371.0,  2520.0,  2677.0,  2838.0,  3004.0,  3173.0,  3351.0,  3538.0,  3730.0,  3929.0, &
              4132.0,  4341.0,  4557.0,  4786.0,  5012.0,  5247.0,  5483.0,  5723.0,  5964.0,  6208.0, &
              6459.0,  6716.0,  6977.0,  7243.0,  7510.0,  7790.0,  8071.0,  8358.0,  8648.0,  8944.0, &
              9244.0,  9561.0,  9881.0, 10207.0, 10535.0, 10871.0, 11215.0, 11564.0, 11919.0, 12284.0, &
             12658.0, 13035.0, 13419.0, 13814.0, 14214.0, 14619.0, 15031.0, 15444.0, 15871.0, 16300.0, &
             16733.1, 17166.3, 17610.0, 18056.8, 18504.1, 18930.0, 19452.0, 19930.0, 20410.0, 20900.0, &
             21390.0, 21880.0, 22360.0/

  data em1/                                                                    29.3,    34.8,    44.3, &
                51.1,    58.7,    66.3,    74.1,    82.3,    91.3,   101.0,   110.8,   122.5,   139.8, &
               159.5,   180.1,   204.7,   229.6,   257.0,   292.8,   326.7,   358.7,   392.0,   430.3, &
               466.6,   506.3,   544.0,   586.1,   628.1,   671.6,   719.0,   772.0,   827.2,   884.7, &
               940.0,  1006.0,  1072.0,  1148.7,  1211.0,  1293.0,  1362.0,  1436.0,  1511.0,  1575.0, &
                 0.0,  1723.0,  1800.0,  1881.0,  1968.0,  2047.0,  2128.0,  2206.0,  2307.0,  2398.0, &
              2491.0,  2601.0,  2708.0,  2820.0,  2932.0,  3049.0,  3174.0,  3296.0,  3425.0,  3562.0, &
              3704.0,  3851.0,  3999.0,  4149.0,  4317.0,  4482.0,  4652.0,  4822.0,  5002.0,  5182.0, &
              5367.0,  5548.0,  5723.0,  5933.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0, &
                 0.0,     0.0,     0.0/

  data em2/                                                                    15.9,    18.3,    25.4, &
                28.3,    32.6,    37.2,    42.2,    47.2,    52.7,    58.9,    68.0,    77.3,    91.4, &
               103.5,   124.9,   146.2,   166.5,   189.0,   222.2,   248.7,   280.3,   310.6,   343.5, &
               376.1,   411.6,   447.6,   483.3,   521.3,   559.9,   603.8,   652.6,   703.2,   756.5, &
               812.7,   870.8,   931.0,  1002.1,  1071.0,  1137.0,  1209.0,  1274.0,  1337.0,  1403.0, &
              1471.4,  1541.0,  1614.0,  1688.0,  1768.0,  1842.0,  1923.0,  2006.0,  2090.0,  2173.0, &
              2264.0,  2365.0,  2469.0,  2575.0,  2682.0,  2792.0,  2909.0,  3027.0,  3148.0,  3279.0, &
              3416.0,  3554.0,  3696.0,  3854.0,  4008.0,  4159.0,  4327.0,  4490.0,  4656.0,  4830.0, &
              5001.0,  5182.0,  5366.0,  5541.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0, &
                 0.0,     0.0,     0.0/

  data em3/                                                                    15.7,    18.3,    25.4, &
                28.3,    32.6,    37.2,    42.2,    47.2,    52.7,    59.9,    66.2,    75.1,    88.6, &
               100.0,   120.8,   141.2,   160.7,   182.0,   214.4,   239.1,   270.0,   298.8,   329.8, &
               360.6,   394.0,   417.7,   461.5,   496.5,   532.3,   573.0,   618.4,   665.3,   714.6, &
               766.4,   820.8,   875.0,   940.6,  1003.0,  1063.0,  1128.0,  1187.0,  1242.0,  1297.0, &
              1357.0,  1419.8,  1481.0,  1544.0,  1611.0,  1676.0,  1741.0,  1812.0,  1885.0,  1950.0, &
              2024.0,  2107.0,  2194.0,  2281.0,  2367.0,  2457.0,  2551.0,  2645.0,  2743.0,  2847.0, &
              2957.0,  3066.0,  3177.0,  3302.0,  3426.0,  3538.0,  3663.0,  3792.0,  3909.0,  4046.0, &
              4174.0,  4304.0,  4435.0,  4557.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0, &
                 0.0,     0.0,     0.0/

  data em4/                                                                                      10.2, &
                18.7,    29.8,    41.7,    55.5,    70.0,    95.0,   113.0,   136.0,   157.7,   181.1, &
               205.0,   231.1,   257.6,   284.2,   311.9,   340.5,   374.0,   411.9,   451.4,   493.2, &
               537.5,   583.4,   630.8,   689.0,   740.5,   795.7,   853.0,   902.4,   948.3,  1003.3, &
              1052.0,  1110.9,  1158.6,  1221.9,  1276.9,  1333.0,  1392.0,  1453.0,  1515.0,  1576.0, &
              1639.0,  1716.0,  1793.0,  1949.0,  1949.0,  2031.0,  2116.0,  2202.0,  2291.0,  2385.0, &
              2485.0,  2586.0,  2688.0,  2798.0,  2909.0,  3022.0,  3136.0,  3248.0,  3370.0,  3491.0, &
              3611.0,  3728.0,  3851.0,  3973.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0, &
                 0.0,     0.0,     0.0/

  data em5/                                                                                      10.1, &
                18.7,    29.2,    41.7,    54.6,    69.0,    93.8,   112.0,   134.2,   155.8,   178.8, &
               202.3,   227.9,   253.9,   280.0,   307.2,   335.2,   368.3,   405.2,   443.9,   484.9, &
               528.2,   573.0,   619.3,   676.4,   726.6,   780.5,   836.0,   883.8,   928.8,   980.4, &
              1027.0,  1083.4,  1127.5,  1189.6,  1241.1,  1292.0,  1351.0,  1409.0,  1468.0,  1528.0, &
              1589.0,  1662.0,  1735.0,  1809.0,  1883.0,  1960.0,  2040.0,  2122.0,  2206.0,  2295.0, &
              2389.0,  2484.0,  2580.0,  2683.0,  2787.0,  2892.0,  3000.0,  3105.0,  3219.0,  3332.0, &
              3442.0,  3552.0,  3666.0,  3775.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0, &
                 0.0,     0.0,     0.0/

  data en1/     27.5,    30.5,    38.9,    43.8,    50.6,    56.4,    63.2,    69.5,    75.0,    81.4, &
                87.1,    97.0,   109.8,   122.9,   137.1,   153.2,   169.4,   186.0,   213.2,   232.3, &
               253.5,   274.7,   291.0,   304.5,   319.2,     0.0,   347.2,   360.0,   378.6,   396.0, &
               414.2,   432.4,   449.8,   470.9,   480.5,   506.8,   538.0,   563.4,   594.1,   625.4, &
               658.2,   691.1,   725.4,   762.1,   802.2,   846.2,   891.8,   939.0,   995.0,  1042.0, &
              1097.0,  1153.0,  1208.0,  1269.0,  1330.0, &
              1387.0,  1439.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0, &
                 0.0,     0.0,     0.0/

  data en2/     14.1,    16.3,    21.6,    24.4,    28.5,    32.6,    37.6,    42.3,    46.3,    50.5, &
                55.7,    63.7,    63.9,    73.5,    83.6,    95.6,   103.3,   123.0,   146.7,   172.4, &
               192.0,   205.8,   223.2,   236.3,   243.3,   242.0,   265.6,   284.0,   286.0,   322.4, &
               333.5,   343.5,   366.2,   385.9,   388.7,   412.4,   438.2,   463.4,   490.4,   518.7, &
               549.1,   577.8,   609.1,   642.7,   680.2,   720.5,   761.9,   805.2,   851.0,   886.0, &
               929.0,   980.0,  1058.0,  1080.0,  1168.0, &
              1224.0,  1271.0,     0.0,     0.0,     0.0,    0.0,     0.0,     0.0,     0.0,     0.0, &
                 0.0,     0.0,     0.0/

  data en3/     14.1,    15.3,    20.1,    23.1,    27.1,    30.8,    35.5,    39.9,    43.2,    47.3, &
                50.9,    58.3,    63.9,    73.5,    83.6,    95.6,   103.3,   123.0,   145.5,   161.3, &
               178.6,   196.0,   206.5,   217.6,   224.6,   242.0,   247.4,   257.0,   271.0,   284.1, &
               293.2,   308.2,   320.2,   332.6,   339.7,   359.2,   380.7,   400.9,   423.6,   446.8, &
               470.7,   495.8,   519.4,   546.3,   576.6,   609.5,   643.5,   678.8,   705.0,   740.0, &
               768.0,   810.0,   879.0,   890.0,   966.4, &
              1007.0,  1043.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,   &
                 0.0,     0.0,     0.0/

  data en4/                       11.7,    17.7,    24.9,    33.3,    41.9,    50.6,    69.5,    79.8, &
                92.6,   105.3,   109.0,   115.1,   120.5,   120.0,   129.0,   133.0,     0.0,   150.5, &
               153.6,   160.0,   167.6,   175.5,   191.2,   206.1,   220.0,   237.9,   255.9,   273.9, &
               293.1,   311.9,   331.6,   353.2,   378.2,   405.7,   434.3,   464.0,   500.0,   533.0, &
               567.0,   603.0,   636.0,   675.0,   712.1, &
               743.0,   778.3,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0, &
                 0.0,     0.0,     0.0/

  data en5/                       10.7,    16.9,    23.9,    32.1,    40.4,    48.9,    67.5,    77.5, &
                89.9,   102.5,     0.0,   115.1,   120.5,   120.0,   129.0,   127.7,   142.6,   150.5, &
               153.6,   160.0,   167.6,   175.5,   182.4,   196.3,   211.5,   226.4,   243.5,   260.5, &
               278.5,   296.3,   314.6,   335.1,   358.8,   385.0,   412.2,   440.1,   473.0,   507.0, &
               541.0,   577.0,   603.0,   639.0,   675.2, &
               708.0,   736.2,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0, &
                 0.0,     0.0,     0.0/

  data en6/                        0.1,     2.0,     1.5,     0.0,     5.2,     0.0,     8.6,     7.7, &
                 8.0,     8.6,     0.0,     0.0,     2.5,     8.9,    15.9,    23.5,    33.6,    42.9, &
                53.4,    63.8,    74.5,    87.6,   104.0,   122.2,   141.7,   162.3,   184.0,   210.0, &
               238.0,   268.0,   299.0,   319.0,   342.4, &
               371.0,   388.2,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,   &
                 0.0,     0.0,     0.0/

  data en7/                        0.1,     2.0,     1.5,     0.0,     5.2,     0.0,     8.6,     2.4, &
                 4.3,     5.2,     4.7,     4.6,     1.3,     7.5,    14.2,    21.6,    31.4,    40.5, &
                50.7,    60.8,    71.2,    83.9,    99.9,   117.8,   136.9,   157.0,   184.0,   210.0, &
               238.0,   268.0,   299.0,   319.0,   333.1, &
               360.0,   377.4,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,   &
                 0.0,     0.0,     0.0/

  data eo1/                            23.3,   22.7,  30.3,    34.3,   37.8,   37.4,   37.5, &
                0.0,   37.4,   32.0,   36.0,   45.6,  49.9,    49.3,   50.6,   54.7,   52.0, &
               57.3,   64.2,   69.7,   75.6,   83.0,  84.0,    95.2,  101.7,  107.2,  127.0, &
              136.0,  147.0,  159.3,  177.0,  195.0, 214.0,   234.0,  254.0,  272.0,  290.0, &
              310.0,  321.0,    0.0,    0.0,    0.0,   0.0,     0.0,    0.0,    0.0,    0.0, &
                0.0,    0.0,    0.0/

  data eo2/                            13.4,   14.2,  17.0,    19.3,   19.8,   22.3,   21.1, &
                0.0,   21.3,   22.0,   20.0,   28.7,  26.3,    30.8,   31.4,   31.8,   30.3, &
               33.6,   38.0,   42.2,   45.3,   45.6,  58.0,    63.0,   65.3,   74.2,   83.1, &
               94.6,  106.4,  119.0,  132.0,  148.0, 164.0,   182.0,  200.0,  215.0,  229.0, &
              232.0,  257.0,    0.0,    0.0,    0.0,   0.0,     0.0,    0.0,    0.0,    0.0, &
                0.0,    0.0,    0.0/

  data eo3/                            12.1,   12.1,   14.8,   16.8,   17.0,   22.3,   21.1, &
                0.0,   21.3,   22.0,   20.0,   22.6,   26.3,   24.1,   24.7,   25.0,   24.1, &
               26.7,   29.9,   32.7,   36.8,   34.6,   44.5,   48.0,   51.7,   57.2,   64.5, &
               73.5,   83.3,   92.6,  104.0,  115.0,  127.0,  140.0,  153.0,  167.0,  182.0, &
              232.0,  192.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0, &
                0.0,    0.0,    0.0/

  data eo4/                                                                             9.6, &
               14.7,   20.7,   26.9,   31.0,   40.0,   48.0,   58.0,   68.0,   80.0,   92.5, &
               94.0,  102.8,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0, &
                0.0,    0.0,    0.0/

  data eo5/                                                                             7.8, &
               12.5,   18.1,   23.8,   31.0,   40.0,   48.0,   58.0,   68.0,   80.0,   85.4, &
               94.0,   94.2,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0, &
                0.0,    0.0,    0.0/

  data ep1/                                            26.0,   34.0,   44.0,    0.0,   41.4, &
                0.0,   43.9,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0, &
                0.0,    0.0,    0.0/

  data ep2/                                                    15.0,   19.0,    0.0,   24.5, &
                0.0,   26.8,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0, &
                0.0,    0.0,    0.0/

  data ep3/                                                    15.0,   19.0,    0.0,   16.6, &
                0.0,   16.8,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0, &
                0.0,    0.0,    0.0/

  Eseuil(:) = 0._db

  if( ( numat > Z_Mendeleiev_max .or. numat <= 0 ) .and. mpirank == 0 ) then
    call write_error
    do ipr = 3,9,3
      write(ipr,105) numat, Z_Mendeleiev_max
    end do
    stop
  endif

  select case(nseuil)

    case(1)
      eseuil(1) = ek1( numat )

    case(2)
      select case(jseuil)
        case(1)
          eseuil(1) = el1( numat )
        case(2)
          eseuil(1) = el2( numat )
          if( nbseuil == 2 ) eseuil(2) = el3( numat )
        case(3)
          eseuil(1) = el3( numat )
      end select

    case(3)
      if( ( numat >= nm1 .and. jseuil < 4 ) .or. ( numat >= nm4 .and. jseuil >= 4 ) ) then
        select case(jseuil)
          case(1)
            eseuil(1) = em1( numat )
          case(2)
            eseuil(1) = em2( numat )
            if( nbseuil == 2 ) eseuil(2) = em3( numat )
          case(3)
            eseuil(1) = em3( numat )
          case(4)
            eseuil(1) = em4( numat )
            if( nbseuil == 2 ) eseuil(2) = em5( numat )
          case(5)
            eseuil(1) = em5( numat )
        end select
      endif

    case(4)
      if( ( numat >= nn1 .and. jseuil < 4 ) .or. ( numat >= nn4 .and. jseuil >= 4 .and. jseuil < 6) .or. &
          ( numat >= nn6 .and. jseuil >= 6 ) ) then
        select case(jseuil)
          case(1)
            eseuil(1) = en1( numat )
          case(2)
            eseuil(1) = en2( numat )
            if( nbseuil == 2 ) eseuil(2) = en3( numat )
          case(3)
            eseuil(1) = en3( numat )
          case(4)
            eseuil(1) = en4( numat )
            if( nbseuil == 2 ) eseuil(2) = en5( numat )
          case(5)
            eseuil(1) = en5( numat )
          case(6)
            eseuil(1) = en6( numat )
            if( nbseuil == 2 ) eseuil(2) = en7( numat )
          case(7)
            eseuil(1) = en7( numat )
        end select
      endif

    case(5)
      if( ( numat >= no1 .and. jseuil < 4 ) .or. ( numat >= no4 .and. jseuil >= 4 ) ) then
        select case(jseuil)
          case(1)
            eseuil(1) = eo1( numat )
          case(2)
            eseuil(1) = eo2( numat )
            if( nbseuil == 2 ) eseuil(2) = eo3( numat )
          case(3)
            eseuil(1) = eo3( numat )
          case(4)
            eseuil(1) = eo4( numat )
            if( nbseuil == 2 ) eseuil(2) = eo5( numat )
          case(5)
            eseuil(1) = eo5( numat )
        end select
      endif

    case(6)
      if( ( numat >= np1 .and. jseuil < 2 ) .or. ( numat >= np2 .and. jseuil >= 2 ) ) then
        select case(jseuil)
          case(1)
            eseuil(1) = ep1( numat )
          case(2)
            eseuil(1) = ep2( numat )
            if( nbseuil == 2 ) eseuil(2) = ep3( numat )
          case(3)
            eseuil(1) = ep3( numat )
        end select
      endif

  end select

! Valeurs non encore tabulees, simplement extrapolees a partir de la derniere valeur connue

  if( nseuil /= 0 .and. Eseuil(1) < eps10 .and. numat > 92 ) then
    if( nseuil == 3 ) then
      Zm = 94
      select case(jseuil)
        case(1)
          eseuil(1) = ( numat - Zm ) * ( em1(Zm) - em1(Zm-1) ) + em1(Zm)
        case(2)
          eseuil(1) = ( numat - Zm ) * ( em2(Zm) - em2(Zm-1) ) + em2(Zm)
        case(3)
          eseuil(1) = ( numat - Zm ) * ( em3(Zm) - em3(Zm-1) ) + em3(Zm)
        case(4)
          eseuil(1) = ( numat - Zm ) * ( em4(Zm) - em4(Zm-1) ) + em4(Zm)
        case(5)
          eseuil(1) = ( numat - Zm ) * ( em5(Zm) - em5(Zm-1) ) + em5(Zm)
      end select
    elseif( nseuil == 4 ) then
      Zm = 92
      select case(jseuil)
        case(1)
          eseuil(1) = ( numat - Zm ) * ( en1(Zm) - en1(Zm-1) ) + en1(Zm)
        case(2)
          eseuil(1) = ( numat - Zm ) * ( en2(Zm) - en2(Zm-1) ) + en2(Zm)
        case(3)
          eseuil(1) = ( numat - Zm ) * ( en3(Zm) - en3(Zm-1) ) + en3(Zm)
        case(4)
          eseuil(1) = ( numat - Zm ) * ( en4(Zm) - en4(Zm-1) ) + en4(Zm)
        case(5)
          eseuil(1) = ( numat - Zm ) * ( en5(Zm) - en5(Zm-1) ) + en5(Zm)
        case(6)
          eseuil(1) = ( numat - Zm ) * ( en6(Zm) - en6(Zm-1) ) + en6(Zm)
        case(7)
          eseuil(1) = ( numat - Zm ) * ( en7(Zm) - en7(Zm-1) ) + en7(Zm)
      end select
    elseif( nseuil == 5 ) then
      Zm = 92
      select case(jseuil)
        case(1)
          eseuil(1) = ( numat - Zm ) * ( eo1(Zm) - eo1(Zm-1) ) + eo1(Zm)
        case(2)
          eseuil(1) = ( numat - Zm ) * ( eo2(Zm) - eo2(Zm-1) ) + eo2(Zm)
        case(3)
          eseuil(1) = ( numat - Zm ) * ( eo3(Zm) - eo3(Zm-1) ) + eo3(Zm)
        case(4)
          eseuil(1) = ( numat - Zm ) * ( eo4(Zm) - eo4(Zm-1) ) + eo4(Zm)
        case(5)
          eseuil(1) = ( numat - Zm ) * ( eo5(Zm) - eo5(Zm-1) ) + eo5(Zm)
      end select
    endif
  endif
  
  if( nseuil /= 0 .and. Eseuil(1) < eps10 .and. mpirank == 0 ) then
    call write_error
    do ipr = 3,9,3
      write(ipr,120)
    end do
    stop
  endif

  Eseuil(:) = nint( 100 * Eseuil(:) ) / 100._db

  if( mpirank == 0 ) then
    do ipr = 3,6,3
      if( ipr == 3 .and. icheck == 0 ) cycle
      if( nbseuil == 1 ) then
        write(ipr,150) Eseuil(:)
      else
        write(ipr,160) Eseuil(:)
      endif
    end do
  endif

  Eseuil(:) = Eseuil(:) / rydb

  if( nseuil /= 0 ) then  
 ! Shift of the conventional edge energy in order its value correspond to the Fermi energy or the HOMO.
    Shift = Shift_edge(icheck,numat)
    Eseuil(:) = Eseuil(:) + Shift
  endif

  return
  105 format(//' Z =',i12,' < 0 or > Z_Mendeleiev_max =',i4,' in routine Esdata !'//)
  120 format(//' Threshold non included in the data of the routine Esdata in the file tab_data.f !',// &
               ' It may be that it does not exist !'//)
  150 format(/' E_edge     =',f9.2,' eV')
  160 format(/' E_edge(1)  =',f9.2,' eV,  E_edge(2) =',f9.2,' eV')
end

!***********************************************************************

function Shift_edge(icheck,Z) 

  use declarations
  implicit none

  integer:: icheck, Z

  real(kind=db):: Shift_edge
  real(kind=db), dimension(Z_Mendeleiev_max):: Shift

  data Shift/ -13.6_db, -24.6_db,   0._db,   0._db,   0._db,   0._db,   0._db,   0._db,   0._db,   0._db, &
                0._db,    0._db,    0._db,   0._db,   0._db,   0._db,   0._db,   0._db,   0._db,   0._db, &
                0._db,    0._db,    0._db,   0._db,   0._db,   0._db,   0._db,   0._db,   0._db,   0._db, &
                0._db,    0._db,    0._db,   0._db,   0._db,   0._db,   0._db,   0._db,   0._db,   0._db, &
                0._db,    0._db,    0._db,   0._db,   0._db,   0._db,   0._db,   0._db,   0._db,   0._db, &
                0._db,    0._db,    0._db,   0._db,   0._db,   0._db,   0._db,   0._db,   0._db,   0._db, &
                0._db,    0._db,    0._db,   0._db,   0._db,   0._db,   0._db,   0._db,   0._db,   0._db, &
                0._db,    0._db,    0._db,   0._db,   0._db,   0._db,   0._db,   0._db,   0._db,   0._db, &
                0._db,    0._db,    0._db,   0._db,   0._db,   0._db,   0._db,   0._db,   0._db,   0._db, &
                0._db,    0._db,    0._db,   0._db,   0._db,   0._db,   0._db,   0._db,   0._db,   0._db, &
                0._db,    0._db,    0._db/
                
  Shift_edge = Shift(Z)
  
  if( icheck > 0 .and. abs(Shift_edge) > eps10 ) write(3,110) Shift_edge

  Shift_edge = Shift_edge / Rydb
  
  return
  110 format(' Shift_edge =',f9.2,' eV')
end

!***********************************************************************

function Workf_val(icheck,numat) 

  use declarations
  implicit none

  integer:: icheck, numat

  real(kind=db):: Workf_val
  real(kind=db), dimension(Z_Mendeleiev_max):: Workfct

! Handbook of chemistry and physics, E-82
! Pour les gaz rares et H, F, N, O, Cl, les seuils sont par rapport au
! niveau du vide. On prend donc le travail de sortie a zero.

  data Workfct/ 0.00_db,  0.00_db,  2.28_db,  3.92_db,  4.50_db,  4.81_db,  0.00_db,  0.00_db,  0.00_db,  0.00_db,  &
                2.28_db,  3.68_db,  4.08_db,  4.37_db,  4.00_db,  4.00_db,  0.00_db,  0.00_db,  2.24_db,  2.706_db, &
                4.00_db,  3.95_db,  3.77_db,  4.37_db,  3.76_db,  4.70_db,  3.90_db,  5.01_db,  4.00_db,  3.08_db,  &
                4.00_db,  4.29_db,  4.72_db,  4.62_db,  4.00_db,  0.00_db,  2.09_db,  2.74_db,  4.00_db,  3.73_db,  &
                2.29_db,  4.15_db,  4.00_db,  4.00_db,  4.57_db,  4.97_db,  4.73_db,  4.07_db,  4.00_db,  4.38_db,  &
                4.01_db,  4.04_db,  4.00_db,  0.00_db,  4.70_db,  2.48_db,  4.00_db,  2.84_db,  4.00_db,  4.00_db,  &
                4.00_db,  4.00_db,  4.00_db,  4.00_db,  4.00_db,  4.00_db,  4.00_db,  4.00_db,  4.00_db,  4.00_db,  &
                4.00_db,  4.00_db,  4.05_db,  4.49_db,  5.00_db,  4.00_db,  4.00_db,  4.09_db,  4.82_db,  4.53_db,  &
                3.68_db,  3.97_db,  4.25_db,  4.00_db,  4.00_db,  0.00_db,  4.00_db,  4.00_db,  4.00_db,  3.38_db,  &
                0.00_db,  3.63_db,  4.00_db,  4.00_db,  4.00_db,  0.00_db,  4.00_db,  4.00_db,  4.00_db,  4.00_db,  &
                4.00_db,  4.00_db,  4.00_db/

  Workf_val = Workfct( numat )
  if( icheck > 0 ) write(3,110) Workf_val

  Workf_val = Workf_val / rydb

  return
 110 format(' WorkF     =',f9.2,' eV')
end

!***********************************************************************

! Determination du sous groupe non magnetique.

function igrpt_sg_so(igrpt,igrpt_nomag)

  select case(igrpt)
    case(39,40)
      igrpt_sg_so = 4
    case(41)
      igrpt_sg_so = 1  ! A ameliorer plus tard...
    case(44)
      igrpt_sg_so = 5
    case(45,46)
      igrpt_sg_so = 16
    case(48)
      igrpt_sg_so = 21
    case(53,58,67)
      igrpt_sg_so = 11
    case(60)
      igrpt_sg_so = 10
    case(70)
      igrpt_sg_so = 17
    case(73,78)
      igrpt_sg_so = 22
    case(84)
      igrpt_sg_so = 23
    case default
      igrpt_sg_so = igrpt_nomag
  end select

  return
end

! **********************************************************************

! Sous groupe acceptant des harmoniques complexes.

function igrpt_sg_cmp(igrpt,igrpt_nomag)

! Sous groupe acceptant les harmoniques complexes
  select case(igrpt)
    case(6,7)             ! C2v (mm2), D2 (222)
      igrpt_sg_cmp = 4    ! C2 (2)
    case(8)               ! D2h (mmm)
      igrpt_sg_cmp = 5    ! C2h (2/m)
    case(12,14)           ! C4v (4mm), D4 (422)
      igrpt_sg_cmp = 9    ! C4 (4)
    case(13)              ! D2d (-42m)
      igrpt_sg_cmp = 10   ! S4 (-4)
    case(15)              ! D4h (4/mmm)
      igrpt_sg_cmp = 11   ! C4h (4/m)
    case default
      igrpt_sg_cmp = igrpt_nomag
  end select

  return
end

! **********************************************************************

  integer function iordresym(is)

  implicit none

  integer:: is
  integer, dimension(64):: ior

  data ior/1,18,21,24,25,28,31,42,49,50,51,52,53,54,55,56, 10,11,22,23,40,41,45,48,57,58,59,60,61,62,63,64, &
           2,3,4,5,6,7,8,9,12,13,14,15,16,17,19,20,26,27,29,30,32, 33,34,35,36,37,38,39,43,44,46,47/

  iordresym = ior(is)

  return
end

!*********************************************************************

! Name of symmetry operations

function nomsym(isym)

  character(len=9):: nomsym
  character(len=9), dimension(64):: nsym

  data nsym/ '        E','   C3_111','  -C3_111','  C3_1-11',' -C3_1-11', &
   '  C3_-111',' -C3_-111','  C3_11-1',' -C3_11-1','   C2_110', &
   '  C2_-110','   C2_101','  C2_-101','   C2_011','  C2_0-11', &
   '      C4x','      C4y','      C4z','     -C4x','     -C4y', &
   '     -C4z','      C2x','      C2y','      C2z','        i', &
   '      S4x','      S4y','      S4z','     -S4x','     -S4y', &
   '     -S4z','  iC3_111',' -iC3_111',' iC3_1-11','-iC3_1-11', &
   ' iC3_-111','-iC3_-111',' iC3_11-1','-iC3_11-1','       mx', &
   '       my','       mz','    d_011','    d_101','    d_110', &
   '   d_01-1','   d_10-1','   d_1-10','      C3z','     -C3z', &
   '      C6z','     -C6z','      S3z','     -S3z','      S6z', &
   '     -S6z','    d_030','   C2_030','    d_060','   C2_060', '    d_120','   C2_120','    d_150','   C2_150'/

  nomsym = nsym(isym)

  return
end

!*********************************************************************

function ptgrname_int_nomag(igrpt)

  use declarations
  implicit none

  integer:: igrpt 
  
  character(len=8):: ptgrname_int_nomag
  character(len=8), dimension(ngrptm+ngrpt_compm):: ptgrname_int_t

  data ptgrname_int_t/ '1       ','-1      ','m       ','2       ','2/m     ', &
      'mm      ','222     ','mmm     ','4       ','-4      ', '4/m     ','4mm     ','-42m    ','422     ','4/mmm   ', &
      '3       ','-3      ','3m      ','32      ','-3m     ', '-6      ','6       ','6/m     ','-6m2    ','6mm     ', &
      '622     ','6/mmm   ','23      ','m3      ','-43m    ', '432     ','m3m     ', &
      'mx      ','my      ','2x      ','2y      ','2/mx    ', '2/my    ','2mm     ','m2m     ','32v     ','dd      ', &
      '3mv     '/

  ptgrname_int_nomag = ptgrname_int_t(igrpt)

  return
end

!***********************************************************************

function ptgrname_int(igrpt)

  use declarations
  implicit none

  integer:: igrpt 
  
  character(len=8):: ptgrname_int, ptgrname_int_nomag
  character(len=8), dimension(ngrptm+1:ngrptmagm+ngrptmag_compm):: ptgrname_mag_t

  data ptgrname_mag_t/ "-1'     ","2'      ","m'      ","2/m'    ","2'/m    ", &
      "2'/m'   ","22'2'   ","2m'm'   ","2'm'm   ","m'm'm'  ", "mmm'    ","m'm'm   ","32'     ","3m'     ","-6'     ", &
      "-6m'2'  ","-6'm2'  ","-6'm'2  ","4'      ","-4'     ", "42'     ","4'2     ","4/m'    ","4'/m'   ","4'/m    ", &
      "4m'm'   ","4'mm'   ","-42'm'  ","-4'2m'  ","-4'2'm  ", "4/m'm'm'","4/m'mm  ","4'/mmm  ","4'/m'm'm","4/mm'm' ", &
      "6'      ","-3'     ","-3m'    ","-3'm    ","-3'm'   ", "62'     ","6'2     ","6/m'    ","6'/m'   ","6'/m    ", &
      "6m'm'   ","6'mm'   ","6'/mm'm ","6'/m'm'm","6/mm'm' ", "6/m'mm  ","6/mm'm' ","m'3     ","-4'3m'  ","4'3     ", &
      "m'3m'   ","m'3m    ","m3m'    ", "22'2'   ","2'22'   ","2'mm'   ","m'mm    ","mm'm    ", &
      "mm'm'   ","m'mm'   ","4'm'm   ","6'2d    ","6'm'm   "/

  if( igrpt > ngrptm ) then
    ptgrname_int = ptgrname_mag_t(igrpt)
  else
    ptgrname_int = ptgrname_int_nomag(igrpt)
  endif

  return
end

!***********************************************************************

function ptgrname_sch(igrpt)

  use declarations
  implicit none

  integer:: igrpt 
  
  character(len=8):: ptgrname_sch
  character(len=8), dimension(ngrptm+ngrpt_compm):: ptgrname_sch_t

  data ptgrname_sch_t/'C1      ','Ci      ','Cs      ','C2      ','C2h     ','C2v     ','D2      ','D2h     ','C4      ', &
                      'S4      ','C4h     ','C4v     ','D2d     ','D4      ','D4h     ','C3      ','S6      ','C3v     ', &
                      'D3      ','D3d     ','C3h     ','C6      ','C6h     ','D3h     ','C6v     ','D6      ','D6h     ', &
                      'T       ','Th      ','Td      ','O       ','Oh      ', &
                      'Csx     ','Csy     ','C2x     ','C2y     ','C2hx    ','C2hy    ','C2vx    ','C2vy    ','D3v     ', &
                      'C2d     ','C3d     '/

  ptgrname_sch = ptgrname_sch_t(igrpt)

  return
end

!***********************************************************************

! Table du nombre des operations de symetrie de chaque groupe ponctuel

function numbops(is)

  use declarations
  implicit none

  integer:: numbops, is 

  integer, dimension(ngrptm+ngrpt_compm):: no_tab

  data no_tab/ 1, 2, 2, 2, 4, 4, 4, 8, 4, 4, 8, 8, 8, 8,16, 3, 6, 6, 6,12, 6, 6,12,12,12,12,24,12,24,24, 24,48, &
               2, 2, 2, 2, 4, 4, 4, 4, 6, 4, 6/

  numbops = no_tab(is)

  return
end

!***********************************************************************

! Table des operations de symetrie de chaque groupe ponctuel

function ios(is)

  integer, dimension(332+40):: os_tab

! L'ordre correspond a celui des tables de caractere.
! C1   1    1
! Ci  -1    2
! Cs   m    3
! C2   2    4
! C2h  2/m  5
! C2v  mm   6
! D2   222  7
! D2h  mmm  8
! C4   4    9
! S4  -4    10
! C4h  4/m  11
! C4v  4mm  12
! D2d -42m  13
! D4   422  14
! D4h 4/mmm 15
! C3   3    16
! S6  -3    17
! C3v  3m   18
! D3   32   19
! D3d -3m   20
! C3h -6    21
! C6   6    22
! C6h  6/m  23
! D3h -6m2  24
! C6v  6mm  25
! D6   622  26
! D6h 6/mmm 27
! T    23   28
! Th   m3   29
! Td  -43m  30
! O    432  31
! Oh   m3m  32
! Csx   mx   33
! Csy   my   34
! C2x   2x   35
! C2y   2y   36
! C2hx  2/mx 37
! C2hy  2/my 38
! C2vx  2mm  39
! C2vy  m2m  40
! D3v   32v  41
! C2d   dd   42
! C3d   3mv  43

  data os_tab/ 1, 1, 25, 1, 42, 1, 24, 1, 24, 25, 42, 1, 24, 41, 40, 1, 22, 24, 23, 1, 24, 23, 22, 25, 42, 41, 40, &
   1, 24, 18, 21, 1, 24, 28, 31, 1, 18, 24, 21, 25, 31, 42, 28, 1, 24, 18,21, 40,41, 45,48, 1, 24, 28,31, 22,23, 45,48, &
   1, 24, 18,21, 22,23, 10,11, 1, 18,21, 24, 22,23, 10,11, 25, 28,31, 42, 40,41, 45,48, 1, 49, 50, 1, 49, 50, 25, 56, 55, &
   1, 49,50, 40,57,63, 1, 49,50, 22,60,62, 1, 49,50, 22,60,62, 25, 55,56, 57,40,63, 1, 49, 50, 42, 53, 54, &
   1, 49, 50, 24, 52, 51, 1, 51, 49, 24, 50, 52, 25, 54, 56, 42, 55, 53, 1, 49,50, 58,23,64, 42, 53,54, 57,40,63, &
   1, 24, 49,50, 51,52, 41,59,61, 57,40,63, 1, 49,50, 22,60,62, 24, 51,52, 58,23,64, &
   1, 51,52, 49,50, 24, 22,60,62, 58,23,64, 25, 53,54, 55,56, 42, 57,40,63, 41,59,61, 1, 22,23,24, 2,4,6,8, 3,5,7,9, &
   1, 2,4,6,8, 3,5,7,9, 22,23,24, 25, 33,35,37,39, 32,34,36,38, 40,41,42, 1, 2,3,4,5,6,7,8,9, 22,23,24, 43,44,45,46,47,48, &
                    26,27,28,29,30,31, 1, 2,3,4,5,6,7,8,9, 22,23,24, 10,11,12,13,14,15, 16,17,18,19,20,21, &
   1, 2,3,4,5,6,7,8,9, 10,11,12,13,14,15, 16,17,18,19,20,21, 22,23,24, 25, 32,33,34,35,36,37,38,39, &
      43,44,45,46,47,48, 26,27,28,29,30,31, 40,41,42, 1, 40, 1, 41, 1, 22, 1, 23, 1, 22, 25, 40, 1, 23, 25, 41, &
   1, 22, 41, 42, 1, 23, 40, 42, 1, 49,50, 58,23,64, 1, 24, 45, 48, 1, 49,50, 41,59,61/

  ios = os_tab(is)

  return
end

!***********************************************************************

function Chemical_Symbol(Z)

  implicit none
  
  integer:: Z

  character(len=2):: Chemical_Symbol
  character(len=2), dimension(118):: Symbol

  data Symbol/ 'H ','He','Li','Be','B ','C ','N ','O ','F ','Ne', 'Na','Mg','Al','Si','P ','S ','Cl','Ar','K ','Ca', &
      'Sc','Ti','V ','Cr','Mn','Fe','Co','Ni','Cu','Zn', 'Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y ','Zr', &
      'Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn', 'Sb','Te','I ','Xe','Cs','Ba','La','Ce','Pr','Nd', &
      'Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb', 'Lu','Hf','Ta','W ','Re','Os','Ir','Pt','Au','Hg', &
      'Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th', 'Pa','U ','Np','Pu','Am','Cm','Bk','Cf','Es','Fm','Md','No','Lw', &
      'Rf','Db','Sg','Bh','Hs','Mt','Ds','Rg','Cn','Nh', 'Fl','Mc','Lv','Ts','Og' /

  Chemical_Symbol = Symbol(Z)

  return
end

!***********************************************************************

function Chemical_Symbol_c(Z)

  implicit none
  integer:: Z

  character(len=2):: Chemical_Symbol_c
  character(len=2), dimension(118):: Symbol

  data Symbol/ 'H ','HE','LI','BE','B ','C ','N ','O ','F ','NE', 'NA','MG','AL','SI','P ','S ','CL','AR','K ','CA', &
      'SC','TI','V ','CR','MN','FE','CO','NI','CU','ZN', 'GA','GE','AS','SE','BR','KR','RB','SR','Y ','ZR', &
      'NB','MO','TC','RU','RH','PD','AG','CD','IN','SN', 'SB','TE','I ','XE','CS','BA','LA','CE','PR','ND', &
      'PM','SM','EU','GD','TB','DY','HO','ER','TM','YB', 'LU','HF','TA','W ','RE','OS','IR','PT','AU','HG', &
      'TL','PB','BI','PO','AT','RN','FR','RA','AC','TH', 'PA','U ','NP','PU','AM','CM','BK','CF','ES','FM','MD','NO','LW', &
      'RF','DB','SG','BH','HS','MT','DS','RG','CN','NH', 'FL','MC','LV','TS','OG' /

  Chemical_Symbol_c = Symbol(Z)

  return
end

!***********************************************************************

function Chemical_Name(Z)

  implicit none
  integer:: Z

  character(len=13):: Chemical_Name
  character(len=13), dimension(118):: Name

  data Name/ 'Hydrogen     ','Helium       ','Lithium      ','Beryllium    ','Boron        ', &
             'Carbon       ','Nitrogen     ','Oxygen       ','Fluorine     ','Neon         ', &
             'Sodium       ','Magnesium    ','Aluminium    ','Silicon      ','Phosphorus   ', &
             'Sulfur       ','Chlorine     ','Argon        ','Potassium    ','Calcium      ', &
             'Scandium     ','Titanium     ','Vanadium     ','Chromium     ','Manganese    ', &
             'Iron         ','Cobalt       ','Nickel       ','Copper       ','Zinc         ', &
             'Gallium      ','Germanium    ','Arsenic      ','Selenium     ','Bromine      ', &
             'Krypton      ','Rubidium     ','Strontium    ','Yttrium      ','Zirconium    ', &
             'Niobium      ','Molybdenum   ','Technetium   ','Ruthenium    ','Rhodium      ', &
             'Palladium    ','Silver       ','Cadmium      ','Indium       ','Tin          ', &
             'Antimony     ','Tellurium    ','Iodine       ','Xenon        ','Cesium       ', &
             'Barium       ','Lanthanum    ','Cerium       ','Praeseodymium','Neodymium    ', &
             'Promethium   ','Samarium     ','Europium     ','Gadolinium   ','Terbium      ', &
             'Dysprosium   ','Holmium      ','Erbium       ','Thulium      ','Ytterbium    ', &
             'Lutetium     ','Hafnium      ','Tantalum     ','Tungsten     ','Rhenium      ', &
             'Osmium       ','Iridium      ','Platinum     ','Gold         ','Mercury      ', &
             'Thallium     ','Lead         ','Bismuth      ','Polonium     ','Astatine     ', &
             'Radon        ','Francium     ','Radium       ','Actinium     ','Thorium      ', &
             'Protactinium ','Uranium      ','Neptunium    ','Plutonium    ','Americium    ', &
             'Curium       ','Berkelium    ','Californium  ','Einsteinium  ','Fermium      ', &
             'Mendelevium  ','Nobelium     ','Lawrencium   ','Rutherfordium','Dubdium      ', &
             'Seagorgium   ','Bohrium      ','Hassium      ','Meitnerium   ','Darmstadtium ', &
             'Roentgenium  ','Copernicium  ','Nihonium     ','Flevorium    ','Moscovium    ', &
             'Livermorium  ','Tennessine   ','Oganesson    '/

  Chemical_Name = Name(Z)

  return
end

!*********************************************************************

! Index l of the orbital on which the Hubbard correction applies

function l_hubbard(Z)

  implicit none
  
  integer:: l_hubbard, Z
  integer, dimension(118):: lh

  data lh/ 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, &
           1, 1, 1, 1, 1, 1, 0, 0, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 0, 0, &
           3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 0, 0, &
           3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1/

  l_hubbard = lh(Z)

  return
end

!*********************************************************************

! Table qui contient le nombre d'orbitales de valence qui sont accesibles
! pour toutes les espèces chimiques: (nonrelativistes, on en a besoin pour nlat)

function nvnonrel(Z)

  implicit none
  
  integer:: nvnonrel, Z
  integer, dimension(118):: nvnr

  data nvnr / 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, &
     3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, &
           4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, &  
           4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4 /

  nvnonrel = nvnr(Z)

  return
end

!*********************************************************************

! Table qui contient le nombre d'orbitales de valence qui sont accessibles
! pour toutes les especes chimiques: (relativistes)

function nvrel(Z)

  implicit none
  
  integer nvrel, Z
  integer, dimension(118):: nvr

  data nvr / 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, &
     5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, &
           7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, &
           7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7 /

  nvrel = nvr(Z)

  return
end

!*********************************************************************

! Number of non relativistic atomic orbitals (by default)

function n_orb_base(Z)

  integer:: n_orb_base, Z
  integer, dimension(118):: n

  data n/ 1, 1, 2, 2, 3, 3, 3, 3, 3, 3, 4, 4, 5, 5, 5, 5, 5, 5, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 8, 8, 8, 8, 8, 8, &
    9, 9,10,10,10,10,10,10,10,10,10,10,11,11,11,11,11,11, 12,12,13,14,14,14,14,14,14,14,14,14,14,14,14,14, &
         14,14,14,14,14,14,14,14,14,14,15,15,15,15,15,15, 16,16,17,18,18,18,18,18,18,18,18,18,18,18,18,18, &
         19,19,19,19,19,19,19,19,19,19,20,20,20,20,20,20/

  n_orb_base = n(Z)

  return
end

!*********************************************************************

! Number of relativistic atomic orbitals
! Need for dimentionning

function n_orb_rel(Z)

  integer:: ipr, n_orb_rel, Z
  integer, dimension(118):: n

  data n/ 1, 1, 2, 2, 4, 4, 4, 4, 4, 4, 5, 5, 7, 7, 7, 7, 7, 7, 8, 8,10,10,10,10,10,10,10,10,10,10,12,12,12,12,12,12, &
   13,13,15,15,15,15,15,15,15,15,15,15,17,17,17,17,17,17,18,18,24,24,24,24,24,24,24,24,24,24,24,24,24,24, &
         24,24,24,24,24,24,24,24,24,24,26,26,26,26,26,26,27,27,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31, &
         31,31,31,31,31,31,31,31,31,33,33,33,33,33,33/

  if( Z < 1 .or. Z > 118 ) then
    call write_error
    do ipr = 3,9,3
      write(ipr,'(/A)') ' Z not possible in function n_orb_rel (in tab_data.f90) !'
    end do
    stop
  endif

  n_orb_rel = n(Z)

  return
end

!*********************************************************************

! Index of the last core orbital. Called by dirgen.
! Order of the occupancy orbital: 1s[He]2s2p[Ne]3s3p[Ar]4s3d4p[Kr]5s4d5p[Xe]6s4f5d6p[Rn]7s6d5f
! inversion: 6s2 4f0 5d1, then one fills the 4f, then the 6p

function n_orb_coeur(Z)

  use declarations
  implicit none
  
  integer:: n_orb_coeur, Z
  integer, dimension(Z_Mendeleiev_max):: nc

  data nc/ 0, 0, 1, 1, 1, 1, 1, 1, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, &
     7, 7, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 9, 9, 9, 9, 9, 10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10, &
             12,12,12,12,12,12,12,12,12,12,12,12,13,13,13, 14,14,14,14,14,14,14,15,15,15,15,15,15,15,15,15,15/

  n_orb_coeur = nc(Z)

  return
end

!*********************************************************************

! Index l of the valence orbital, by default partialy occupied

function l_level_val(Z)

  integer:: Z, l_level_val

  select case( Z )
    case(1,2,3,4,11,12,19,20,37,38,55,56,87,88)
      l_level_val = 0
    case(5,6,7,8,9,13,14,15,16,17,18,31,32,33,34,35,36, 49,50,51,52,53,54,81,82,83,84,85,86)
      l_level_val = 1
    case(21,22,23,24,25,26,27,28,29,30,39,40,41,42,43,44,45, 46,47,48,57,72,73,74,75,76,77,78,79,80,89)
      l_level_val = 2
    case default
      l_level_val = 3
  end select

  return
end

!*********************************************************************

! Atomic mass

function Mass_atom(Z)

  use declarations
  implicit none

  integer:: Z

  real(kind=db):: Mass_atom
  real(kind=db), dimension(Z_Mendeleiev_max):: Mass

  data Mass /  1.0079,  4.0026,  6.941,   9.0122,  10.811,  12.0107, 14.0067, 15.9994, 18.9984,  20.1797, &
              22.9898, 24.3050, 26.9815, 28.0855,  30.9738, 32.066,  35.4527, 39.948,  39.0983,  40.078, &
              44.9559, 47.867,  50.9415, 51.9961,  54.9380, 55.845,  58.9332, 58.6934, 63.546,   65.39, &
              69.72,   72.59,   74.922,  78.96,    79.91,   83.8,    85.47,   87.62,   88.91,    91.22, &
              92.91,   95.94,   98.91,  101.07,   102.9,   106.4,   107.87,  112.4,   114.82,   118.69, &
             121.75,  127.6,   126.9,   131.3,    132.91,  137.34,  138.91,  140.12,  140.91,   144.24, &
             145,     150.35,  151.96,  157.25,   158.92,  162.5,   164.93,  167.26,  168.93,   173.04, &
             174.97,  178.49,  180.95,  183.85,   186.2,   190.2,   192.22,  195.09,  196.97,   200.59, &
             204.37,  207.19,  208.98,  210,      210,     222,     223,     226,     227,      232.04, &
             231,     238.03,  237.05,  244,      243,     247,     247,     251,     252,      257,  258,  259,  262/

  Mass_atom = Mass(Z)

  return
end

!*********************************************************************

! Atom radius
! Covalent, single bond, radius
! From Z = 58 up to 70, and Z > 90, metallic radius

function Atom_radius(Z)

  use declarations
  implicit none

  integer:: Z

  real(kind=db):: Atom_radius
  real(kind=db), dimension(Z_Mendeleiev_max):: Ray

  data Ray/ 0.38_db,  0.32_db,  1.34_db,  0.90_db,  0.82_db,  0.77_db,  0.75_db,  0.73_db,  0.71_db,  0.69_db, &
            1.54_db,  1.30_db,  1.18_db,  1.11_db,  1.06_db,  1.02_db,  0.99_db,  0.97_db,  1.96_db,  1.74_db, &
            1.44_db,  1.36_db,  1.25_db,  1.27_db,  1.39_db,  1.25_db,  1.26_db,  1.21_db,  1.38_db,  1.31_db, &
            1.26_db,  1.22_db,  1.19_db,  1.16_db,  1.14_db,  1.10_db,  2.11_db,  1.92_db,  1.62_db,  1.48_db, &
            1.37_db,  1.45_db,  1.56_db,  1.26_db,  1.35_db,  1.31_db,  1.53_db,  1.48_db,  1.44_db,  1.41_db, &
            1.38_db,  1.35_db,  1.33_db,  1.30_db,  2.25_db,  1.98_db,  1.69_db,  1.818_db, 1.824_db, 1.814_db,&
            1.834_db, 1.804_db, 1.804_db, 1.804_db, 1.773_db, 1.781_db, 1.762_db, 1.761_db, 1.759_db, 1.76_db, &
            1.60_db,  1.50_db,  1.38_db,  1.46_db,  1.59_db,  1.28_db,  1.37_db,  1.28_db,  1.44_db,  1.49_db, &
            1.48_db,  1.47_db,  1.46_db,  1.29_db,  1.38_db,  1.33_db,  1.33_db,  1.59_db,  1.40_db,  1.79_db, &
            1.63_db,  1.56_db,  1.55_db,  1.59_db,  1.73_db,  1.74_db,  1.70_db,  1.86_db,  1.86_db,  1.86_db, &
            1.86_db,  1.86_db,  1.86_db/

  Atom_radius = Ray(Z) / bohr

  return
end

!*********************************************************************

! Ionic radius
! Come from Shannon (1976), most common valence in octahedral coordinance.
! For rare earth (Z = 2,10,18,36,54,86) and carbone (6), it is the atomic radius.
! For astatine (85), Polonium (84) value.
! For Fm (100) and Md (101), Es (99) value.
! For Lr (103), No (102) value.

function RayIon(Z)

  use declarations
  implicit none

  integer:: Z

  real(kind=db):: Rayion
  real(kind=db), dimension(Z_Mendeleiev_max):: Ray

  data Ray/ 0.012_db, 0.49_db,  0.76_db,  0.45_db,  0.27_db,  0.91_db,  1.46_db,  1.40_db,  1.33_db,  0.51_db, &
            1.02_db,  0.67_db,  0.48_db,  0.40_db,  0.44_db,  1.84_db,  1.81_db,  0.88_db,  1.38_db,  1.00_db, &
            0.745_db, 0.86_db,  0.79_db,  0.80_db,  0.83_db,  0.78_db,  0.74_db,  0.69_db,  0.73_db,  0.74_db, &
            0.62_db,  0.73_db,  0.58_db,  1.98_db,  1.96_db,  1.03_db,  1.52_db,  1.18_db,  1.02_db,  0.72_db, &
            0.72_db,  0.69_db,  0.64_db,  0.68_db,  0.665_db, 0.86_db,  1.15_db,  0.95_db,  0.80_db,  0.69_db, &
            0.76_db,  2.21_db,  2.20_db,  1.24_db,  1.67_db,  1.35_db,  1.061_db, 1.034_db, 1.013_db, 0.995_db,&
            0.979_db, 0.964_db, 0.947_db, 0.938_db, 0.923_db, 0.912_db, 0.901_db, 0.881_db, 0.869_db, 0.858_db,&
            0.848_db, 0.71_db,  0.64_db,  0.62_db,  0.56_db,  0.63_db,  0.625_db, 0.625_db, 0.85_db,  1.02_db, &
            1.5_db,   1.19_db,  1.03_db,  2.3_db,   2.3_db,   1.34_db,  1.8_db,   1.43_db,  1.119_db, 0.972_db,&
            0.78_db,  0.52_db,  0.75_db,  0.887_db, 0.982_db, 0.97_db,  0.949_db, 0.934_db, 0.925_db, 0.925_db,&
            0.925_db, 1.1_db,   1.1_db/

  RayIon = Ray(Z)

  return
end

!*********************************************************************

! Debye temperature
! source: Ashcroft and Mermin; the Debye temperatures correspond to the free element
! Unknown Debye temperatures are taken as 1

function Debye_temperature(Z)

  use declarations
  implicit none

  integer:: ipr, Z
  
  real(kind=db):: Debye_temperature 
  real(kind=db), dimension(Z_Mendeleiev_max):: TD

  data TD /  110._db,  26._db, 400._db,1000._db,1250._db,1860._db,  79._db,  46._db,   1._db,  63._db, &
             150._db, 318._db, 394._db, 625._db,   1._db,   1._db,   1._db,  85._db, 100._db, 230._db, &
             359._db, 380._db, 390._db, 460._db, 400._db, 420._db, 385._db, 375._db, 315._db, 234._db, &
             240._db, 360._db, 285._db, 150._db,   1._db,  73._db,  56._db, 147._db, 256._db, 250._db, &
             275._db, 380._db,   1._db, 382._db, 350._db, 275._db, 215._db, 120._db, 129._db, 170._db, &
             200._db, 139._db,   1._db,  55._db,  40._db, 110._db, 132._db, 139._db, 152._db, 157._db, &
               1._db, 160._db, 107._db, 176._db, 188._db, 186._db, 191._db, 196._db, 200._db, 118._db, &
             207._db,   1._db, 225._db, 310._db, 416._db, 400._db, 430._db, 230._db, 170._db, 100._db, &
              96._db,  88._db, 120._db,   1._db,   1._db,   1._db,   1._db,   1._db,   1._db, 100._db, &
               1._db, 210._db, 188._db, 150._db,   1._db,   1._db,   1._db,   1._db,   1._db,   1._db, &
               1._db,   1._db,   1._db/

  if( Z <= 0 .or. Z > Z_Mendeleiev_max ) then
    call write_error
    do ipr = 3,9,3
      write(ipr,110) Z, Z_Mendeleiev_max
    end do
    stop
  endif

  Debye_temperature = TD(Z)

  return
  110 format(//' Z =',i12,' < 0 or > Z_Mendeleiev_max =',i4,' in function Debye_temperature !'//)
end

!*********************************************************************

! lmax value for the COOP

function lmax_coop(Z)

  use declarations
  implicit none

  integer:: ipr, lmax_coop, Z
  integer, dimension(Z_Mendeleiev_max):: lmax_c

  data lmax_c/1,                                  1, &  ! He
              1,1,                      2,2,2,2,2,2, &  ! Ne
              2,2,                      2,2,2,2,2,2, &  ! Ar
              2,2, 2,2,2,2,2,2,2,2,2,2, 2,2,2,2,2,2, &  ! Kr 36  
              2,2, 2,2,2,2,2,2,2,2,2,2, 2,2,2,2,2,2, &  ! Xe 54   
              3,3, 3,3,3,3,3,3,3,3,3,3,3,3,3,3, 3,3,3,3,3,3,3,3,3,3, 3,3,3,3,3,3, & ! Rn 86     
              3,3, 3,3,3,3,3,3,3,3,3,3,3,3,3,3, 3/ ! Lw 103     

  if( Z <= 0 .or. Z > Z_Mendeleiev_max ) then
    call write_error
    do ipr = 3,9,3
      write(ipr,110) Z, Z_Mendeleiev_max
    end do
    stop
  endif

  lmax_coop = lmax_c(Z)

  return
  110 format(//' Z =',i12,' < 0 or > Z_Mendeleiev_max =',i4,' in function lmax_coop !'//)
end

 
