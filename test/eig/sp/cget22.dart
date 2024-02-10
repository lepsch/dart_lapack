      void cget22(final int TRANSA, final int TRANSE, final int TRANSW, final int N, final Matrix<double> A, final int LDA, final Matrix<double> E, final int LDE, final int W, final Array<double> _WORK, final Array<double> RWORK, final int RESULT) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             TRANSA, TRANSE, TRANSW;
      int                LDA, LDE, N;
      double               RESULT( 2 ), RWORK( * );
      Complex            A( LDA, * ), E( LDE, * ), W( * ), WORK( * );
      // ..

      double               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      Complex            CZERO, CONE;
      const              CZERO = ( 0.0, 0.0 ), CONE = ( 1.0, 0.0 ) ;
      String             NORMA, NORME;
      int                ITRNSE, ITRNSW, J, JCOL, JOFF, JROW, JVEC;
      double               ANORM, ENORM, ENRMAX, ENRMIN, ERRNRM, TEMP1, ULP, UNFL;
      Complex            WTEMP;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- REAL               CLANGE, SLAMCH;
      // EXTERNAL lsame, CLANGE, SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL CGEMM, CLASET
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, AIMAG, CONJG, MAX, MIN, REAL

      // Initialize RESULT (in case N=0)

      RESULT[1] = ZERO;
      RESULT[2] = ZERO;
      if (N <= 0) return;

      UNFL = SLAMCH( 'Safe minimum' );
      ULP = SLAMCH( 'Precision' );

      ITRNSE = 0;
      ITRNSW = 0;
      NORMA = 'O';
      NORME = 'O';

      if ( lsame( TRANSA, 'T' ) || lsame( TRANSA, 'C' ) ) {
         NORMA = 'I';
      }

      if ( lsame( TRANSE, 'T' ) ) {
         ITRNSE = 1;
         NORME = 'I';
      } else if ( lsame( TRANSE, 'C' ) ) {
         ITRNSE = 2;
         NORME = 'I';
      }

      if ( lsame( TRANSW, 'C' ) ) {
         ITRNSW = 1;
      }

      // Normalization of E:

      ENRMIN = ONE / ULP;
      ENRMAX = ZERO;
      if ( ITRNSE == 0 ) {
         for (JVEC = 1; JVEC <= N; JVEC++) { // 20
            TEMP1 = ZERO;
            for (J = 1; J <= N; J++) { // 10
               TEMP1 = max( TEMP1, ABS( double( E( J, JVEC ) ) )+ ABS( AIMAG( E( J, JVEC ) ) ) );
            } // 10
            ENRMIN = min( ENRMIN, TEMP1 );
            ENRMAX = max( ENRMAX, TEMP1 );
         } // 20
      } else {
         for (JVEC = 1; JVEC <= N; JVEC++) { // 30
            RWORK[JVEC] = ZERO;
         } // 30

         for (J = 1; J <= N; J++) { // 50
            for (JVEC = 1; JVEC <= N; JVEC++) { // 40
               RWORK[JVEC] = max( RWORK( JVEC ), ABS( double( E( JVEC, J ) ) )+ ABS( AIMAG( E( JVEC, J ) ) ) );
            } // 40
         } // 50

         for (JVEC = 1; JVEC <= N; JVEC++) { // 60
            ENRMIN = min( ENRMIN, RWORK( JVEC ) );
            ENRMAX = max( ENRMAX, RWORK( JVEC ) );
         } // 60
      }

      // Norm of A:

      ANORM = max( CLANGE( NORMA, N, N, A, LDA, RWORK ), UNFL );

      // Norm of E:

      ENORM = max( CLANGE( NORME, N, N, E, LDE, RWORK ), ULP );

      // Norm of error:

      // Error =  AE - EW

      claset('Full', N, N, CZERO, CZERO, WORK, N );

      JOFF = 0;
      for (JCOL = 1; JCOL <= N; JCOL++) { // 100
         if ( ITRNSW == 0 ) {
            WTEMP = W( JCOL );
         } else {
            WTEMP = CONJG( W( JCOL ) );
         }

         if ( ITRNSE == 0 ) {
            for (JROW = 1; JROW <= N; JROW++) { // 70
               WORK[JOFF+JROW] = E( JROW, JCOL )*WTEMP;
            } // 70
         } else if ( ITRNSE == 1 ) {
            for (JROW = 1; JROW <= N; JROW++) { // 80
               WORK[JOFF+JROW] = E( JCOL, JROW )*WTEMP;
            } // 80
         } else {
            for (JROW = 1; JROW <= N; JROW++) { // 90
               WORK[JOFF+JROW] = CONJG( E( JCOL, JROW ) )*WTEMP;
            } // 90
         }
         JOFF = JOFF + N;
      } // 100

      cgemm(TRANSA, TRANSE, N, N, N, CONE, A, LDA, E, LDE, -CONE, WORK, N );

      ERRNRM = CLANGE( 'One', N, N, WORK, N, RWORK ) / ENORM;

      // Compute RESULT(1) (avoiding under/overflow)

      if ( ANORM > ERRNRM ) {
         RESULT[1] = ( ERRNRM / ANORM ) / ULP;
      } else {
         if ( ANORM < ONE ) {
            RESULT[1] = ONE / ULP;
         } else {
            RESULT[1] = min( ERRNRM / ANORM, ONE ) / ULP;
         }
      }

      // Compute RESULT(2) : the normalization error in E.

      RESULT[2] = max( ( ENRMAX-ONE ).abs(), ( ENRMIN-ONE ).abs() ) / ( REAL( N )*ULP );

      }
