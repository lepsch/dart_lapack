      double             BASE, EMAX, EMIN, EPS, PREC, RMAX, RMIN, RND, SFMIN, T;
      // ..
      // .. External Functions ..
      //- double             DLAMCH;
      // EXTERNAL DLAMCH

      EPS   = dlamch( 'Epsilon' );
      SFMIN = dlamch( 'Safe minimum' );
      BASE  = dlamch( 'Base' );
      PREC  = dlamch( 'Precision' );
      T     = dlamch( 'Number of digits in mantissa' );
      RND   = dlamch( 'Rounding mode' );
      EMIN  = dlamch( 'Minimum exponent' );
      RMIN  = dlamch( 'Underflow threshold' );
      EMAX  = dlamch( 'Largest exponent' );
      RMAX  = dlamch( 'Overflow threshold' );

      WRITE( 6, * )' Epsilon                      = ', EPS;
      WRITE( 6, * )' Safe minimum                 = ', SFMIN;
      WRITE( 6, * )' Base                         = ', BASE;
      WRITE( 6, * )' Precision                    = ', PREC;
      WRITE( 6, * )' Number of digits in mantissa = ', T;
      WRITE( 6, * )' Rounding mode                = ', RND;
      WRITE( 6, * )' Minimum exponent             = ', EMIN;
      WRITE( 6, * )' Underflow threshold          = ', RMIN;
      WRITE( 6, * )' Largest exponent             = ', EMAX;
      WRITE( 6, * )' Overflow threshold           = ', RMAX;
      WRITE( 6, * )' Reciprocal of safe minimum   = ', 1 / SFMIN;

      }
