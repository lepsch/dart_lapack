      int                NMAX, ITS;
      const              NMAX = 1000, ITS = 50000 ;
      // ..
      // .. Local Scalars ..
      int                I, J;
      REAL               ALPHA, AVG, T1, T2, TNOSEC, TOTAL;
      // ..
      // .. Local Arrays ..
      REAL               X( NMAX ), Y( NMAX );
      // ..
      // .. External Functions ..
      REAL               SECOND;
      // EXTERNAL SECOND
      // ..
      // .. External Subroutines ..
      // EXTERNAL MYSUB
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC REAL
      // ..
      // .. Executable Statements ..

// .. Figure TOTAL flops ..
      TOTAL = REAL(NMAX) * REAL(ITS) * 2.0;

      // Initialize X and Y

      for (I = 1; I <= NMAX; I++) { // 10
         X( I ) = REAL( 1 ) / REAL( I );
         Y( I ) = REAL( NMAX-I ) / REAL( NMAX );
      } // 10
      ALPHA = 0.315;

      // Time TOTAL SAXPY operations

      T1 = SECOND( );
      for (J = 1; J <= ITS; J++) { // 30
         for (I = 1; I <= NMAX; I++) { // 20
            Y( I ) = Y( I ) + ALPHA*X( I );
         } // 20
         ALPHA = -ALPHA;
      } // 30
      T2 = SECOND( );
      TNOSEC = T2 - T1;
      WRITE( 6, 9999 )TOTAL, TNOSEC;
      if ( TNOSEC > 0.0 ) {
         WRITE( 6, 9998 )(TOTAL/1.0e6)/TNOSEC;
      } else {
         WRITE( 6, 9994 );
      }

      // Time TOTAL SAXPY operations with SECOND in the outer loop

      T1 = SECOND( );
      for (J = 1; J <= ITS; J++) { // 50
         for (I = 1; I <= NMAX; I++) { // 40
            Y( I ) = Y( I ) + ALPHA*X( I );
         } // 40
         ALPHA = -ALPHA;
         T2 = SECOND( );
      } // 50

      // Compute the time used in milliseconds used by an average call
      // to SECOND.

      WRITE( 6, 9997 )T2 - T1;
      AVG = ( ( T2-T1 ) - TNOSEC ) * 1000.0e+00/REAL( ITS );
      if (AVG > 0.0) WRITE( 6, 9996 )AVG;

      // Compute the equivalent number of floating point operations used
      // by an average call to SECOND.

      if(( AVG > 0.0 ) && ( TNOSEC > 0.0 )) WRITE( 6, 9995 )(AVG/1000) * TOTAL / TNOSEC;

 9999 FORMAT( ' Time for ', G10.3,' SAXPY ops = ', G10.3, ' seconds' );
 9998 FORMAT( ' SAXPY performance rate        = ', G10.3, ' mflops ' );
 9997 FORMAT( ' Including SECOND, time        = ', G10.3, ' seconds' );
 9996 FORMAT( ' Average time for SECOND       = ', G10.3, ' milliseconds' );
 9995 FORMAT( ' Equivalent floating point ops = ', G10.3, ' ops' );
 9994 FORMAT( ' *** Warning:  Time for operations was less or equal', ' than zero => timing in TESTING might be dubious' );
      mysub(NMAX,X,Y);
      }
      void mysub(N,X,Y) {
      int     N;
      REAL X(N), Y(N);
      return;
      }
