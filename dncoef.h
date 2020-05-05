c====================== include file "dncoef.h" ========================
c
c @(#) SCCS module: dncoef.h  version: 1.10
c     Creation date: 08/09/96
c
c-----------------------------------------------------------------------
c     Enter normalizing temperatures, salinities and coefficients
c     generated by program "eqstat" which fits 3rd order 
c     polynomials to the equation 'of state for each model level.
c-----------------------------------------------------------------------
c
      data to /
     &        13.4979166,13.4926169,13.4845768,13.4727122, 8.4665976,
     &         8.4484312, 5.9367561, 4.4197725, 2.9018104, 2.8627440,
     &         2.8119188, 2.7478969, 2.6693535, 2.5760430, 2.9494306/
c
      data so /
     &        -0.0022500,-0.0022500,-0.0022500,-0.0022500, 0.0001500,
     &         0.0001500,-0.0001000,-0.0001000,-0.0002500,-0.0002500,
     &        -0.0002000,-0.0002000,-0.0002000,-0.0002000,-0.0002000/
c
      data (c( 1,n),n=1,9)/
     &         -.2017851E-03,0.7729900E+00,-.4916357E-05,-.2021094E-02,
     &         0.3168210E+00,0.3606242E-07,0.3773530E-02,0.3601930E-04,
     &         0.1605978E+01/
c
      data (c( 2,n),n=1,9)/
     &         -.2024842E-03,0.7727355E+00,-.4904976E-05,-.2017314E-02,
     &         0.3160354E+00,0.3591818E-07,0.3768401E-02,0.3595768E-04,
     &         0.1598916E+01/
c
      data (c( 3,n),n=1,9)/
     &         -.2035380E-03,0.7723525E+00,-.4887711E-05,-.2011439E-02,
     &         0.3151638E+00,0.3569715E-07,0.3736588E-02,0.3584786E-04,
     &         0.1627915E+01/
c
      data (c( 4,n),n=1,9)/
     &         -.2050843E-03,0.7717891E+00,-.4862535E-05,-.2003172E-02,
     &         0.3129073E+00,0.3537948E-07,0.3740574E-02,0.3572830E-04,
     &         0.1634156E+01/
c
      data (c( 5,n),n=1,9)/
     &         -.1628465E-03,0.7822322E+00,-.5256112E-05,-.2310488E-02,
     &         0.4441855E+00,0.4386074E-07,0.4209105E-02,0.4284397E-04,
     &         -.3533121E+01/
c
      data (c( 6,n),n=1,9)/
     &         -.1664290E-03,0.7810346E+00,-.5197709E-05,-.2290309E-02,
     &         0.4061066E+00,0.4321431E-07,0.4105291E-02,0.4270665E-04,
     &         -.4451626E+01/
c
      data (c( 7,n),n=1,9)/
     &         -.1446374E-03,0.7851275E+00,-.5443910E-05,-.2465505E-02,
     &         0.1508644E+00,0.4865008E-07,-.2298690E-02,0.4636501E-04,
     &         0.2921848E+02/
c
      data (c( 8,n),n=1,9)/
     &         -.1353437E-03,0.7868303E+00,-.5536851E-05,-.2558586E-02,
     &         0.1452740E+00,0.5092440E-07,-.1494541E-02,0.4538149E-04,
     &         -.2139407E+03/
c
      data (c( 9,n),n=1,9)/
     &         -.1276864E-03,0.7879311E+00,-.5628878E-05,-.2711959E-02,
     &         -.1701061E+01,0.5181664E-07,-.1254796E+00,0.5131700E-04,
     &         -.2181232E+04/
c
      data (c(10,n),n=1,9)/
     &         -.1392459E-03,0.7838731E+00,-.5445710E-05,-.2609743E-02,
     &         -.4988049E+01,0.5096849E-07,0.7100123E-01,0.4921809E-04,
     &         0.4546029E+04/
c
      data (c(11,n),n=1,9)/
     &         -.1529192E-03,0.7808815E+00,-.5195704E-05,-.2513633E-02,
     &         -.7580082E+01,0.4536900E-07,0.1053939E+01,0.1641317E-05,
     &         -.1635490E+05/
c
      data (c(12,n),n=1,9)/
     &         -.1684001E-03,0.7768872E+00,-.4939702E-05,-.2213298E-02,
     &         -.1041420E+02,0.4656355E-07,0.1212965E+01,-.6386872E-04,
     &         -.2994899E+05/
c
      data (c(13,n),n=1,9)/
     &         -.1853807E-03,0.7711285E+00,-.4658209E-05,-.2091224E-02,
     &         -.1515848E+02,0.4140049E-07,0.2341507E+01,0.7015846E-04,
     &         -.3064146E+05/
c
      data (c(14,n),n=1,9)/
     &         -.2033667E-03,0.7655853E+00,-.4382178E-05,-.1989417E-02,
     &         -.2318759E+02,0.3566120E-07,0.3340599E+01,-.1688436E-04,
     &         -.3545847E+05/
c
      data (c(15,n),n=1,9)/
     &         -.2256227E-03,0.7588712E+00,-.3992833E-05,-.1987264E-02,
     &         -.1391063E+02,0.3392166E-07,0.3334040E+01,0.1038726E-04,
     &         -.5790932E+05/
c
c  the above coefficients were calculated using program "eqstat"
c  (employing the UNESCO equation of state)
c  and are valid for the following depths and T and S ranges
c     k      depth        tmin      tmax           smin      smax
c     1      15.00e2     -2.000    29.000        28.5000   37.0000
c     2      53.08e2     -2.000    29.000        28.5000   37.0000
c     3     110.61e2     -2.000    29.000        28.5000   37.0000
c     4     195.04e2     -2.000    29.000        28.5000   37.0000
c     5     315.32e2     -2.000    19.000        33.7000   36.6000
c     6     481.69e2     -2.000    19.000        33.7000   36.6000
c     7     705.13e2     -2.000    14.000        34.0000   35.8000
c     8     996.49e2     -2.000    11.000        34.1000   35.7000
c     9    1365.37e2     -1.000     7.000        34.4000   35.1000
c    10    1818.82e2     -1.000     7.000        34.5000   35.0000
c    11    2360.06e2     -1.000     7.000        34.6000   35.0000
c    12    2987.28e2     -1.000     7.000        34.6000   35.0000
c    13    3693.03e2     -1.000     7.000        34.6000   35.0000
c    14    4464.06e2     -1.000     7.000        34.6000   35.0000
c    15    5281.94e2      0.000     7.000        34.6000   35.0000
c
c  the 3rd order polynomial will return density departures [gm/cm**3] as
c  a function of potential temperature [deg C] & salinity [model units]
c
c       k level reference densities (in sigma units):
c   24.6436 24.8139 25.0711 25.4479 28.7717 29.5230 30.7254 32.2579
c   34.0132 36.0659 38.5291 41.3037 44.3820 47.6937 51.0399
