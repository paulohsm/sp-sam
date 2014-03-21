DSET ^/Users/santiago/Tupa/WORK_HOME/Modelos/AGCM-1.0/pos/dataout/TQ0062L028.1997060100SEMIPROG/GPOSNMC1997060100%y4%m2%d2%h2C.fct.TQ0062L028
*
OPTIONS TEMPLATE SEQUENTIAL YREV BIG_ENDIAN
*
UNDEF -2.56E+33
*
TITLE PRESSURE HISTORY    CPTEC AGCM REVIS 1.0 2000  T062L28  COLD
*
XDEF   192 LINEAR    0.000   1.8750000000
YDEF    96 LEVELS 
 -88.57217 -86.72253 -84.86197 -82.99894 -81.13498 -79.27056 -77.40589 -75.54106
 -73.67613 -71.81113 -69.94608 -68.08099 -66.21587 -64.35073 -62.48557 -60.62040
 -58.75521 -56.89001 -55.02481 -53.15960 -51.29438 -49.42915 -47.56393 -45.69869
 -43.83346 -41.96822 -40.10298 -38.23774 -36.37249 -34.50724 -32.64199 -30.77674
 -28.91149 -27.04624 -25.18099 -23.31573 -21.45048 -19.58522 -17.71996 -15.85470
 -13.98945 -12.12419 -10.25893  -8.39367  -6.52841  -4.66315  -2.79789  -0.93263
   0.93263   2.79789   4.66315   6.52841   8.39367  10.25893  12.12419  13.98945
  15.85470  17.71996  19.58522  21.45048  23.31573  25.18099  27.04624  28.91149
  30.77674  32.64199  34.50724  36.37249  38.23774  40.10298  41.96822  43.83346
  45.69869  47.56393  49.42915  51.29438  53.15960  55.02481  56.89001  58.75521
  60.62040  62.48557  64.35073  66.21587  68.08099  69.94608  71.81113  73.67613
  75.54106  77.40589  79.27056  81.13498  82.99894  84.86197  86.72253  88.57217
ZDEF    28 LEVELS  1000  975  950  925  900  875  850  825  800  775
                  750  725  700  675  650  600  500  300  250  200
                  150  100   70   50   30   20   10    3
TDEF     1464 LINEAR 01Z01JUN1997 1HR
*
VARS    60
TOPO    0 99 TOPOGRAPHY                              (M               )
LSMK    0 99 LAND SEA MASK                           (NO DIM          )
PSLC    0 99 SURFACE PRESSURE                        (Mb              )
UVEL   28 99 ZONAL WIND (U)                          (M/Sec           )
VVEL   28 99 MERIDIONAL WIND (V)                     (M/Sec           )
TEMP   28 99 ABSOLUTE TEMPERATURE                    (K               )
UMES   28 99 SPECIFIC HUMIDITY                       (No Dim          )
TSFC    0 99 SURFACE TEMPERATURE                     (K               )
USSL    0 99 SOIL WETNESS OF SURFACE                 (No Dim          )
UZRS    0 99 SOIL WETNESS OF ROOT ZONE               (No Dim          )
UZDS    0 99 SOIL WETNESS OF DRAINAGE ZONE           (No Dim          )
T02M    0 99 TEMPERATURE AT 2-M FROM SURFACE         (K               )
Q02M    0 99 SPECIFIC HUMIDITY AT 2-M FROM SURFACE   (No Dim          )
U10M    0 99 ZONAL WIND AT 10-M FROM SURFACE         (M/Sec           )
V10M    0 99 MERID WIND AT 10-M FROM SURFACE         (M/Sec           )
PSMT    0 99 TIME MEAN SURFACE PRESSURE              (Mb              )
UVMT   28 99 TIME MEAN ZONAL WIND (U)                (M/Sec           )
VVMT   28 99 TIME MEAN MERIDIONAL WIND (V)           (M/Sec           )
GHMT   28 99 TIME MEAN GEOPOTENTIAL HEIGHT           (M               )
SPMT    0 99 TIME MEAN SEA LEVEL PRESSURE            (Mb-1000         )
TSMT    0 99 TIME MEAN SURFACE ABSOLUTE TEMPERATURE  (K               )
ATMT   28 99 TIME MEAN ABSOLUTE TEMPERATURE          (K               )
RSMT    0 99 TIME MEAN SURFACE RELATIVE HUMIDITY     (No Dim          )
RHMT   28 99 TIME MEAN RELATIVE HUMIDITY             (No Dim          )
SHMT   28 99 TIME MEAN SPECIFIC HUMIDITY             (No Dim          )
PCMT    0 99 TIME MEAN PRECIP. WATER                 (Kg M**-2        )
STTM    0 99 TIME MEAN SURFACE TEMPERATURE           (K               )
OMTM   28 99 TIME MEAN OMEGA                         (Cb/Sec          )
PREC    0 99 TOTAL PRECIPITATION                     (Kg M**-2 Day**-1)
PRCV    0 99 CONVECTIVE PRECIPITATION                (Kg M**-2 Day**-1)
NEVE    0 99 SNOWFALL                                (Kg M**-2 Day**-1)
RNOF    0 99 RUNOFF                                  (Kg M**-2 Day**-1)
CSSF    0 99 SENSIBLE HEAT FLUX FROM SURFACE         (W M**-2         )
CLSF    0 99 LATENT HEAT FLUX FROM SURFACE           (W M**-2         )
EVAP    0 99 EVAPORATION                             (Kg M**-2 Day**-1)
USST    0 99 SURFACE ZONAL WIND STRESS               (Pa              )
VSST    0 99 SURFACE MERIDIONAL WIND STRESS          (Pa              )
CBNV    0 99 CLOUD COVER                             (No Dim          )
OLIS    0 99 DOWNWARD LONG WAVE AT BOTTOM            (W M**-2         )
OLES    0 99 UPWARD LONG WAVE AT BOTTOM              (W M**-2         )
ROLE    0 99 OUTGOING LONG WAVE AT TOP               (W M**-2         )
ISWF    0 99 INCIDENT SHORT WAVE FLUX                (W M**-2         )
OCIS    0 99 DOWNWARD SHORT WAVE AT GROUND           (W M**-2         )
OCES    0 99 UPWARD SHORT WAVE AT GROUND             (W M**-2         )
ROCE    0 99 UPWARD SHORT WAVE AT TOP                (W M**-2         )
CVLH   28 99 CONVECTIVE LATENT HEATING               (K/Sec           )
CVMS   28 99 CONVECTIVE MOISTURE SOURCE              (1/Sec           )
OLIC    0 99 DOWNWARD LONG WAVE AT BOTTOM (CLEAR)    (W M**-2         )
LWTC    0 99 OUTGOING LONG WAVE AT TOP (CLEAR)       (W M**-2         )
OCIC    0 99 DOWNWARD SHORT WAVE AT GROUND (CLEAR)   (W M**-2         )
SWGC    0 99 UPWARD SHORT WAVE AT GROUND (CLEAR)     (W M**-2         )
SWTC    0 99 UPWARD SHORT WAVE AT TOP (CLEAR)        (W M**-2         )
TGSC    0 99 GROUND/SURFACE COVER TEMPERATURE        (K               )
VDCC   28 99 VERTICAL DIST TOTAL CLOUD COVER         (No Dim          )
WTNV   28 99 CLOUD LIQUID WATER PATH                 (No Dim          )
T2MT    0 99 TIME MEAN TEMP AT 2-M FROM SFC          (K               )
Q2MT    0 99 TIME MEAN SPEC HUMIDITY AT 2-M FROM SFC (No Dim          )
S10M    0 99 SPEED WIND AT 10-M FROM SURFACE         (M/Sec           )
U10T    0 99 TIME MEAN AT 10 METRE U-WIND COMPONENT  (M/Sec           )
V10T    0 99 TIME MEAN AT 10 METRE V-WIND COMPONENT  (M/Sec           )
ENDVARS
