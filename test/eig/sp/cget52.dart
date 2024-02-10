      void cget52(LEFT, N, final Matrix<double> A, final int LDA, final Matrix<double> B, final int LDB, final Matrix<double> E, final int LDE, ALPHA, BETA, WORK, final Array<double> RWORK, RESULT ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      bool               LEFT;
      int                LDA, LDB, LDE, N;
      double               RESULT( 2 ), RWORK( * );
      Complex            A( LDA, * ), ALPHA( * ), B( LDB, * ), BETA( * ), E( LDE, * ), WORK( * );
      // ..

      double               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      Complex            CZERO, CONE;
      const              CZERO = ( 0.0, 0.0 ), CONE = ( 1.0, 0.0 ) ;
      String             NORMAB, TRANS;
      int                J, JVEC;
      double               ABMAX, ALFMAX, ANORM, BETMAX, BNORM, ENORM, ENRMER, ERRNRM, SAFMAX, SAFMIN, SCALE, TEMP1, ULP;
      Complex            ACOEFF, ALPHAI, BCOEFF, BETAI, X;
      // ..
      // .. External Functions ..
      //- REAL               CLANGE, SLAMCH;
      // EXTERNAL CLANGE, SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL CGEMV
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, AIMAG, CONJG, MAX, REAL
      // ..
      // .. Statement Functions ..
      double               ABS1;
      // ..
      // .. Statement Function definitions ..
      ABS1[X] = ( double( X ) ).abs() + ( AIMAG( X ) ).abs();

      RESULT[1] = ZERO;
      RESULT[2] = ZERO;
      if (N <= 0) return;

      SAFMIN = SLAMCH( 'Safe minimum' );
      SAFMAX = ONE / SAFMIN;
      ULP = SLAMCH( 'Epsilon' )*SLAMCH( 'Base' );

      if ( LEFT ) {
         TRANS = 'C';
         NORMAB = 'I';
      } else {
         TRANS = 'N';
         NORMAB = 'O';
      }

      // Norm of A, B, and E:

      ANORM = max( CLANGE( NORMAB, N, N, A, LDA, RWORK ), SAFMIN );
      BNORM = max( CLANGE( NORMAB, N, N, B, LDB, RWORK ), SAFMIN );
      ENORM = max( CLANGE( 'O', N, N, E, LDE, RWORK ), ULP );
      ALFMAX = SAFMAX / max( ONE, BNORM );
      BETMAX = SAFMAX / max( ONE, ANORM );

      // Compute error matrix.
      // Column i = ( b(i) A - a(i) B ) E(i) / max( |a(i) B|, |b(i) A| )

      for (JVEC = 1; JVEC <= N; JVEC++) { // 10
         ALPHAI = ALPHA( JVEC );
         BETAI = BETA( JVEC );
         ABMAX = max( ABS1( ALPHAI ), ABS1( BETAI ) );
         if ( ABS1( ALPHAI ) > ALFMAX || ABS1( BETAI ) > BETMAX || ABMAX < ONE ) {
            SCALE = ONE / max( ABMAX, SAFMIN );
            ALPHAI = SCALE*ALPHAI;
            BETAI = SCALE*BETAI;
         }
         SCALE = ONE / max( ABS1( ALPHAI )*BNORM, ABS1( BETAI )*ANORM, SAFMIN );
         ACOEFF = SCALE*BETAI;
         BCOEFF = SCALE*ALPHAI;
         if ( LEFT ) {
            ACOEFF = CONJG( ACOEFF );
            BCOEFF = CONJG( BCOEFF );
         }
         cgemv(TRANS, N, N, ACOEFF, A, LDA, E( 1, JVEC ), 1, CZERO, WORK( N*( JVEC-1 )+1 ), 1 );
         cgemv(TRANS, N, N, -BCOEFF, B, LDA, E( 1, JVEC ), 1, CONE, WORK( N*( JVEC-1 )+1 ), 1 );
      } // 10

      ERRNRM = CLANGE( 'One', N, N, WORK, N, RWORK ) / ENORM;

      // Compute RESULT(1)

      RESULT[1] = ERRNRM / ULP;

      // Normalization of E:

      ENRMER = ZERO;
      for (JVEC = 1; JVEC <= N; JVEC++) { // 30
         TEMP1 = ZERO;
         for (J = 1; J <= N; J++) { // 20
            TEMP1 = max( TEMP1, ABS1( E( J, JVEC ) ) );
         } // 20
         ENRMER = max( ENRMER, ( TEMP1-ONE ).abs() );
      } // 30

      // Compute RESULT(2) : the normalization error in E.

      RESULT[2] = ENRMER / ( REAL( N )*ULP );

      }
