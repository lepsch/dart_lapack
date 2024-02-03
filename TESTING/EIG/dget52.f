      SUBROUTINE DGET52( LEFT, N, A, LDA, B, LDB, E, LDE, ALPHAR, ALPHAI, BETA, WORK, RESULT )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      bool               LEFT;
      int                LDA, LDB, LDE, N;
      // ..
      // .. Array Arguments ..
      double             A( LDA, * ), ALPHAI( * ), ALPHAR( * ), B( LDB, * ), BETA( * ), E( LDE, * ), RESULT( 2 ), WORK( * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, ONE, TEN;
      const              ZERO = 0.0D0, ONE = 1.0D0, TEN = 10.0D0 ;
      // ..
      // .. Local Scalars ..
      bool               ILCPLX;
      String             NORMAB, TRANS;
      int                J, JVEC;
      double             ABMAX, ACOEF, ALFMAX, ANORM, BCOEFI, BCOEFR, BETMAX, BNORM, ENORM, ENRMER, ERRNRM, SAFMAX, SAFMIN, SALFI, SALFR, SBETA, SCALE, TEMP1, ULP;
      // ..
      // .. External Functions ..
      double             DLAMCH, DLANGE;
      // EXTERNAL DLAMCH, DLANGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL DGEMV
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, MAX
      // ..
      // .. Executable Statements ..

      RESULT( 1 ) = ZERO
      RESULT( 2 ) = ZERO
      if (N.LE.0) RETURN;

      SAFMIN = DLAMCH( 'Safe minimum' )
      SAFMAX = ONE / SAFMIN
      ULP = DLAMCH( 'Epsilon' )*DLAMCH( 'Base' )

      if ( LEFT ) {
         TRANS = 'T'
         NORMAB = 'I'
      } else {
         TRANS = 'N'
         NORMAB = 'O'
      }

      // Norm of A, B, and E:

      ANORM = MAX( DLANGE( NORMAB, N, N, A, LDA, WORK ), SAFMIN )
      BNORM = MAX( DLANGE( NORMAB, N, N, B, LDB, WORK ), SAFMIN )
      ENORM = MAX( DLANGE( 'O', N, N, E, LDE, WORK ), ULP )
      ALFMAX = SAFMAX / MAX( ONE, BNORM )
      BETMAX = SAFMAX / MAX( ONE, ANORM )

      // Compute error matrix.
      // Column i = ( b(i) A - a(i) B ) E(i) / max( |a(i) B|, |b(i) A| )

      ILCPLX = false;
      for (JVEC = 1; JVEC <= N; JVEC++) { // 10
         if ( ILCPLX ) {

            // 2nd Eigenvalue/-vector of pair -- do nothing

            ILCPLX = false;
         } else {
            SALFR = ALPHAR( JVEC )
            SALFI = ALPHAI( JVEC )
            SBETA = BETA( JVEC )
            if ( SALFI == ZERO ) {

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
               dgemv(TRANS, N, N, ACOEF, A, LDA, E( 1, JVEC ), 1, ZERO, WORK( N*( JVEC-1 )+1 ), 1 );
               dgemv(TRANS, N, N, -BCOEFR, B, LDA, E( 1, JVEC ), 1, ONE, WORK( N*( JVEC-1 )+1 ), 1 );
            } else {

               // Complex conjugate pair

               ILCPLX = true;
               if ( JVEC == N ) {
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

               dgemv(TRANS, N, N, ACOEF, A, LDA, E( 1, JVEC ), 1, ZERO, WORK( N*( JVEC-1 )+1 ), 1 );
               dgemv(TRANS, N, N, -BCOEFR, B, LDA, E( 1, JVEC ), 1, ONE, WORK( N*( JVEC-1 )+1 ), 1 );
               dgemv(TRANS, N, N, BCOEFI, B, LDA, E( 1, JVEC+1 ), 1, ONE, WORK( N*( JVEC-1 )+1 ), 1 );

               dgemv(TRANS, N, N, ACOEF, A, LDA, E( 1, JVEC+1 ), 1, ZERO, WORK( N*JVEC+1 ), 1 );
               dgemv(TRANS, N, N, -BCOEFI, B, LDA, E( 1, JVEC ), 1, ONE, WORK( N*JVEC+1 ), 1 );
               dgemv(TRANS, N, N, -BCOEFR, B, LDA, E( 1, JVEC+1 ), 1, ONE, WORK( N*JVEC+1 ), 1 );
            }
         }
      } // 10

      ERRNRM = DLANGE( 'One', N, N, WORK, N, WORK( N**2+1 ) ) / ENORM

      // Compute RESULT(1)

      RESULT( 1 ) = ERRNRM / ULP

      // Normalization of E:

      ENRMER = ZERO
      ILCPLX = false;
      for (JVEC = 1; JVEC <= N; JVEC++) { // 40
         if ( ILCPLX ) {
            ILCPLX = false;
         } else {
            TEMP1 = ZERO
            if ( ALPHAI( JVEC ) == ZERO ) {
               for (J = 1; J <= N; J++) { // 20
                  TEMP1 = MAX( TEMP1, ABS( E( J, JVEC ) ) )
               } // 20
               ENRMER = MAX( ENRMER, ABS( TEMP1-ONE ) )
            } else {
               ILCPLX = true;
               for (J = 1; J <= N; J++) { // 30
                  TEMP1 = MAX( TEMP1, ABS( E( J, JVEC ) )+ ABS( E( J, JVEC+1 ) ) )
               } // 30
               ENRMER = MAX( ENRMER, ABS( TEMP1-ONE ) )
            }
         }
      } // 40

      // Compute RESULT(2) : the normalization error in E.

      RESULT( 2 ) = ENRMER / ( DBLE( N )*ULP )

      RETURN

      // End of DGET52

      }
