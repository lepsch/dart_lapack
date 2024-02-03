      SUBROUTINE DGET22( TRANSA, TRANSE, TRANSW, N, A, LDA, E, LDE, WR, WI, WORK, RESULT );

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             TRANSA, TRANSE, TRANSW;
      int                LDA, LDE, N;
      // ..
      // .. Array Arguments ..
      double             A( LDA, * ), E( LDE, * ), RESULT( 2 ), WI( * ), WORK( * ), WR( * );
      // ..

// =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      // ..
      // .. Local Scalars ..
      String             NORMA, NORME;
      int                IECOL, IEROW, INCE, IPAIR, ITRNSE, J, JCOL, JVEC;
      double             ANORM, ENORM, ENRMAX, ENRMIN, ERRNRM, TEMP1, ULP, UNFL;
      // ..
      // .. Local Arrays ..
      double             WMAT( 2, 2 );
      // ..
      // .. External Functions ..
      bool               LSAME;
      double             DLAMCH, DLANGE;
      // EXTERNAL LSAME, DLAMCH, DLANGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL DAXPY, DGEMM, DLASET
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, MAX, MIN
      // ..
      // .. Executable Statements ..

      // Initialize RESULT (in case N=0)

      RESULT( 1 ) = ZERO;
      RESULT( 2 ) = ZERO;
      if (N <= 0) RETURN;

      UNFL = DLAMCH( 'Safe minimum' );
      ULP = DLAMCH( 'Precision' );

      ITRNSE = 0;
      INCE = 1;
      NORMA = 'O';
      NORME = 'O';

      if ( LSAME( TRANSA, 'T' ) || LSAME( TRANSA, 'C' ) ) {
         NORMA = 'I';
      }
      if ( LSAME( TRANSE, 'T' ) || LSAME( TRANSE, 'C' ) ) {
         NORME = 'I';
         ITRNSE = 1;
         INCE = LDE;
      }

      // Check normalization of E

      ENRMIN = ONE / ULP;
      ENRMAX = ZERO;
      if ( ITRNSE == 0 ) {

         // Eigenvectors are column vectors.

         IPAIR = 0;
         for (JVEC = 1; JVEC <= N; JVEC++) { // 30
            TEMP1 = ZERO;
            IF( IPAIR == 0 && JVEC < N && WI( JVEC ) != ZERO ) IPAIR = 1;
            if ( IPAIR == 1 ) {

               // Complex eigenvector

               for (J = 1; J <= N; J++) { // 10
                  TEMP1 = MAX( TEMP1, ABS( E( J, JVEC ) )+ ABS( E( J, JVEC+1 ) ) );
               } // 10
               ENRMIN = MIN( ENRMIN, TEMP1 );
               ENRMAX = MAX( ENRMAX, TEMP1 );
               IPAIR = 2;
            } else if ( IPAIR == 2 ) {
               IPAIR = 0;
            } else {

               // Real eigenvector

               for (J = 1; J <= N; J++) { // 20
                  TEMP1 = MAX( TEMP1, ABS( E( J, JVEC ) ) );
               } // 20
               ENRMIN = MIN( ENRMIN, TEMP1 );
               ENRMAX = MAX( ENRMAX, TEMP1 );
               IPAIR = 0;
            }
         } // 30

      } else {

         // Eigenvectors are row vectors.

         for (JVEC = 1; JVEC <= N; JVEC++) { // 40
            WORK( JVEC ) = ZERO;
         } // 40

         for (J = 1; J <= N; J++) { // 60
            IPAIR = 0;
            for (JVEC = 1; JVEC <= N; JVEC++) { // 50
               IF( IPAIR == 0 && JVEC < N && WI( JVEC ) != ZERO ) IPAIR = 1;
               if ( IPAIR == 1 ) {
                  WORK( JVEC ) = MAX( WORK( JVEC ), ABS( E( J, JVEC ) )+ABS( E( J, JVEC+1 ) ) );
                  WORK( JVEC+1 ) = WORK( JVEC );
               } else if ( IPAIR == 2 ) {
                  IPAIR = 0;
               } else {
                  WORK( JVEC ) = MAX( WORK( JVEC ), ABS( E( J, JVEC ) ) );
                  IPAIR = 0;
               }
            } // 50
         } // 60

         for (JVEC = 1; JVEC <= N; JVEC++) { // 70
            ENRMIN = MIN( ENRMIN, WORK( JVEC ) );
            ENRMAX = MAX( ENRMAX, WORK( JVEC ) );
         } // 70
      }

      // Norm of A:

      ANORM = MAX( DLANGE( NORMA, N, N, A, LDA, WORK ), UNFL );

      // Norm of E:

      ENORM = MAX( DLANGE( NORME, N, N, E, LDE, WORK ), ULP );

      // Norm of error:

      // Error =  AE - EW

      dlaset('Full', N, N, ZERO, ZERO, WORK, N );

      IPAIR = 0;
      IEROW = 1;
      IECOL = 1;

      for (JCOL = 1; JCOL <= N; JCOL++) { // 80
         if ( ITRNSE == 1 ) {
            IEROW = JCOL;
         } else {
            IECOL = JCOL;
         }

         IF( IPAIR == 0 && WI( JCOL ) != ZERO ) IPAIR = 1;

         if ( IPAIR == 1 ) {
            WMAT( 1, 1 ) = WR( JCOL );
            WMAT( 2, 1 ) = -WI( JCOL );
            WMAT( 1, 2 ) = WI( JCOL );
            WMAT( 2, 2 ) = WR( JCOL );
            dgemm(TRANSE, TRANSW, N, 2, 2, ONE, E( IEROW, IECOL ), LDE, WMAT, 2, ZERO, WORK( N*( JCOL-1 )+1 ), N );
            IPAIR = 2;
         } else if ( IPAIR == 2 ) {
            IPAIR = 0;

         } else {

            daxpy(N, WR( JCOL ), E( IEROW, IECOL ), INCE, WORK( N*( JCOL-1 )+1 ), 1 );
            IPAIR = 0;
         }

      } // 80

      dgemm(TRANSA, TRANSE, N, N, N, ONE, A, LDA, E, LDE, -ONE, WORK, N );

      ERRNRM = DLANGE( 'One', N, N, WORK, N, WORK( N*N+1 ) ) / ENORM;

      // Compute RESULT(1) (avoiding under/overflow)

      if ( ANORM > ERRNRM ) {
         RESULT( 1 ) = ( ERRNRM / ANORM ) / ULP;
      } else {
         if ( ANORM < ONE ) {
            RESULT( 1 ) = ONE / ULP;
         } else {
            RESULT( 1 ) = MIN( ERRNRM / ANORM, ONE ) / ULP;
         }
      }

      // Compute RESULT(2) : the normalization error in E.

      RESULT( 2 ) = MAX( ABS( ENRMAX-ONE ), ABS( ENRMIN-ONE ) ) / ( DBLE( N )*ULP );

      return;

      // End of DGET22

      }
