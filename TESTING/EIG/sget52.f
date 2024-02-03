      SUBROUTINE SGET52( LEFT, N, A, LDA, B, LDB, E, LDE, ALPHAR, ALPHAI, BETA, WORK, RESULT )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      bool               LEFT;
      int                LDA, LDB, LDE, N;
      // ..
      // .. Array Arguments ..
      REAL               A( LDA, * ), ALPHAI( * ), ALPHAR( * ), B( LDB, * ), BETA( * ), E( LDE, * ), RESULT( 2 ), WORK( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE, TEN
      const              ZERO = 0.0, ONE = 1.0, TEN = 10.0 ;
      // ..
      // .. Local Scalars ..
      bool               ILCPLX;
      String             NORMAB, TRANS;
      int                J, JVEC;
      REAL               ABMAX, ACOEF, ALFMAX, ANORM, BCOEFI, BCOEFR, BETMAX, BNORM, ENORM, ENRMER, ERRNRM, SAFMAX, SAFMIN, SALFI, SALFR, SBETA, SCALE, TEMP1, ULP
      // ..
      // .. External Functions ..
      REAL               SLAMCH, SLANGE
      // EXTERNAL SLAMCH, SLANGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL SGEMV
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, REAL
      // ..
      // .. Executable Statements ..

      RESULT( 1 ) = ZERO
      RESULT( 2 ) = ZERO
      IF( N.LE.0 ) RETURN

      SAFMIN = SLAMCH( 'Safe minimum' )
      SAFMAX = ONE / SAFMIN
      ULP = SLAMCH( 'Epsilon' )*SLAMCH( 'Base' )

      if ( LEFT ) {
         TRANS = 'T'
         NORMAB = 'I'
      } else {
         TRANS = 'N'
         NORMAB = 'O'
      }

      // Norm of A, B, and E:

      ANORM = MAX( SLANGE( NORMAB, N, N, A, LDA, WORK ), SAFMIN )
      BNORM = MAX( SLANGE( NORMAB, N, N, B, LDB, WORK ), SAFMIN )
      ENORM = MAX( SLANGE( 'O', N, N, E, LDE, WORK ), ULP )
      ALFMAX = SAFMAX / MAX( ONE, BNORM )
      BETMAX = SAFMAX / MAX( ONE, ANORM )

      // Compute error matrix.
      // Column i = ( b(i) A - a(i) B ) E(i) / max( |a(i) B|, |b(i) A| )

      ILCPLX = .FALSE.
      for (JVEC = 1; JVEC <= N; JVEC++) { // 10
         if ( ILCPLX ) {

            // 2nd Eigenvalue/-vector of pair -- do nothing

            ILCPLX = .FALSE.
         } else {
            SALFR = ALPHAR( JVEC )
            SALFI = ALPHAI( JVEC )
            SBETA = BETA( JVEC )
            if ( SALFI.EQ.ZERO ) {

               // Real eigenvalue and -vector

               ABMAX = MAX( ABS( SALFR ), ABS( SBETA ) )
               if ( ABS( SALFR ).GT.ALFMAX .OR. ABS( SBETA ).GT. BETMAX .OR. ABMAX.LT.ONE ) {
                  SCALE = ONE / MAX( ABMAX, SAFMIN )
                  SALFR = SCALE*SALFR
                  SBETA = SCALE*SBETA
               }
               SCALE = ONE / MAX( ABS( SALFR )*BNORM, ABS( SBETA )*ANORM, SAFMIN )
               ACOEF = SCALE*SBETA
               BCOEFR = SCALE*SALFR
               sgemv(TRANS, N, N, ACOEF, A, LDA, E( 1, JVEC ), 1, ZERO, WORK( N*( JVEC-1 )+1 ), 1 )                CALL SGEMV( TRANS, N, N, -BCOEFR, B, LDA, E( 1, JVEC ), 1, ONE, WORK( N*( JVEC-1 )+1 ), 1 );
            } else {

               // Complex conjugate pair

               ILCPLX = .TRUE.
               if ( JVEC.EQ.N ) {
                  RESULT( 1 ) = TEN / ULP
                  RETURN
               }
               ABMAX = MAX( ABS( SALFR )+ABS( SALFI ), ABS( SBETA ) )
               if ( ABS( SALFR )+ABS( SALFI ).GT.ALFMAX .OR. ABS( SBETA ).GT.BETMAX .OR. ABMAX.LT.ONE ) {
                  SCALE = ONE / MAX( ABMAX, SAFMIN )
                  SALFR = SCALE*SALFR
                  SALFI = SCALE*SALFI
                  SBETA = SCALE*SBETA
               }
               SCALE = ONE / MAX( ( ABS( SALFR )+ABS( SALFI ) )*BNORM, ABS( SBETA )*ANORM, SAFMIN )
               ACOEF = SCALE*SBETA
               BCOEFR = SCALE*SALFR
               BCOEFI = SCALE*SALFI
               if ( LEFT ) {
                  BCOEFI = -BCOEFI
               }

               sgemv(TRANS, N, N, ACOEF, A, LDA, E( 1, JVEC ), 1, ZERO, WORK( N*( JVEC-1 )+1 ), 1 )                CALL SGEMV( TRANS, N, N, -BCOEFR, B, LDA, E( 1, JVEC ), 1, ONE, WORK( N*( JVEC-1 )+1 ), 1 )                CALL SGEMV( TRANS, N, N, BCOEFI, B, LDA, E( 1, JVEC+1 ), 1, ONE, WORK( N*( JVEC-1 )+1 ), 1 );

               sgemv(TRANS, N, N, ACOEF, A, LDA, E( 1, JVEC+1 ), 1, ZERO, WORK( N*JVEC+1 ), 1 )                CALL SGEMV( TRANS, N, N, -BCOEFI, B, LDA, E( 1, JVEC ), 1, ONE, WORK( N*JVEC+1 ), 1 )                CALL SGEMV( TRANS, N, N, -BCOEFR, B, LDA, E( 1, JVEC+1 ), 1, ONE, WORK( N*JVEC+1 ), 1 );
            }
         }
      } // 10

      ERRNRM = SLANGE( 'One', N, N, WORK, N, WORK( N**2+1 ) ) / ENORM

      // Compute RESULT(1)

      RESULT( 1 ) = ERRNRM / ULP

      // Normalization of E:

      ENRMER = ZERO
      ILCPLX = .FALSE.
      for (JVEC = 1; JVEC <= N; JVEC++) { // 40
         if ( ILCPLX ) {
            ILCPLX = .FALSE.
         } else {
            TEMP1 = ZERO
            if ( ALPHAI( JVEC ).EQ.ZERO ) {
               for (J = 1; J <= N; J++) { // 20
                  TEMP1 = MAX( TEMP1, ABS( E( J, JVEC ) ) )
               } // 20
               ENRMER = MAX( ENRMER, ABS( TEMP1-ONE ) )
            } else {
               ILCPLX = .TRUE.
               for (J = 1; J <= N; J++) { // 30
                  TEMP1 = MAX( TEMP1, ABS( E( J, JVEC ) )+ ABS( E( J, JVEC+1 ) ) )
               } // 30
               ENRMER = MAX( ENRMER, ABS( TEMP1-ONE ) )
            }
         }
      } // 40

      // Compute RESULT(2) : the normalization error in E.

      RESULT( 2 ) = ENRMER / ( REAL( N )*ULP )

      RETURN

      // End of SGET52

      }
