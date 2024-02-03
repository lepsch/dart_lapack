      SUBROUTINE CGET52( LEFT, N, A, LDA, B, LDB, E, LDE, ALPHA, BETA, WORK, RWORK, RESULT )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      bool               LEFT;
      int                LDA, LDB, LDE, N;
      // ..
      // .. Array Arguments ..
      REAL               RESULT( 2 ), RWORK( * )
      COMPLEX            A( LDA, * ), ALPHA( * ), B( LDB, * ), BETA( * ), E( LDE, * ), WORK( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE
      const              ZERO = 0.0E+0, ONE = 1.0E+0 ;
      COMPLEX            CZERO, CONE
      const              CZERO = ( 0.0E+0, 0.0E+0 ), CONE = ( 1.0E+0, 0.0E+0 ) ;
      // ..
      // .. Local Scalars ..
      String             NORMAB, TRANS;
      int                J, JVEC;
      REAL               ABMAX, ALFMAX, ANORM, BETMAX, BNORM, ENORM, ENRMER, ERRNRM, SAFMAX, SAFMIN, SCALE, TEMP1, ULP
      COMPLEX            ACOEFF, ALPHAI, BCOEFF, BETAI, X
      // ..
      // .. External Functions ..
      REAL               CLANGE, SLAMCH
      // EXTERNAL CLANGE, SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL CGEMV
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, AIMAG, CONJG, MAX, REAL
      // ..
      // .. Statement Functions ..
      REAL               ABS1
      // ..
      // .. Statement Function definitions ..
      ABS1( X ) = ABS( REAL( X ) ) + ABS( AIMAG( X ) )
      // ..
      // .. Executable Statements ..

      RESULT( 1 ) = ZERO
      RESULT( 2 ) = ZERO
      if (N.LE.0) RETURN;

      SAFMIN = SLAMCH( 'Safe minimum' )
      SAFMAX = ONE / SAFMIN
      ULP = SLAMCH( 'Epsilon' )*SLAMCH( 'Base' )

      if ( LEFT ) {
         TRANS = 'C'
         NORMAB = 'I'
      } else {
         TRANS = 'N'
         NORMAB = 'O'
      }

      // Norm of A, B, and E:

      ANORM = MAX( CLANGE( NORMAB, N, N, A, LDA, RWORK ), SAFMIN )
      BNORM = MAX( CLANGE( NORMAB, N, N, B, LDB, RWORK ), SAFMIN )
      ENORM = MAX( CLANGE( 'O', N, N, E, LDE, RWORK ), ULP )
      ALFMAX = SAFMAX / MAX( ONE, BNORM )
      BETMAX = SAFMAX / MAX( ONE, ANORM )

      // Compute error matrix.
      // Column i = ( b(i) A - a(i) B ) E(i) / max( |a(i) B|, |b(i) A| )

      for (JVEC = 1; JVEC <= N; JVEC++) { // 10
         ALPHAI = ALPHA( JVEC )
         BETAI = BETA( JVEC )
         ABMAX = MAX( ABS1( ALPHAI ), ABS1( BETAI ) )
         if ( ABS1( ALPHAI ).GT.ALFMAX || ABS1( BETAI ).GT.BETMAX || ABMAX < ONE ) {
            SCALE = ONE / MAX( ABMAX, SAFMIN )
            ALPHAI = SCALE*ALPHAI
            BETAI = SCALE*BETAI
         }
         SCALE = ONE / MAX( ABS1( ALPHAI )*BNORM, ABS1( BETAI )*ANORM, SAFMIN )
         ACOEFF = SCALE*BETAI
         BCOEFF = SCALE*ALPHAI
         if ( LEFT ) {
            ACOEFF = CONJG( ACOEFF )
            BCOEFF = CONJG( BCOEFF )
         }
         cgemv(TRANS, N, N, ACOEFF, A, LDA, E( 1, JVEC ), 1, CZERO, WORK( N*( JVEC-1 )+1 ), 1 );
         cgemv(TRANS, N, N, -BCOEFF, B, LDA, E( 1, JVEC ), 1, CONE, WORK( N*( JVEC-1 )+1 ), 1 );
      } // 10

      ERRNRM = CLANGE( 'One', N, N, WORK, N, RWORK ) / ENORM

      // Compute RESULT(1)

      RESULT( 1 ) = ERRNRM / ULP

      // Normalization of E:

      ENRMER = ZERO
      for (JVEC = 1; JVEC <= N; JVEC++) { // 30
         TEMP1 = ZERO
         for (J = 1; J <= N; J++) { // 20
            TEMP1 = MAX( TEMP1, ABS1( E( J, JVEC ) ) )
         } // 20
         ENRMER = MAX( ENRMER, ABS( TEMP1-ONE ) )
      } // 30

      // Compute RESULT(2) : the normalization error in E.

      RESULT( 2 ) = ENRMER / ( REAL( N )*ULP )

      RETURN

      // End of CGET52

      }
