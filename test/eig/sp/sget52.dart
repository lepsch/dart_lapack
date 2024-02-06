      void sget52(LEFT, N, A, LDA, B, LDB, E, LDE, ALPHAR, ALPHAI, BETA, WORK, RESULT ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      bool               LEFT;
      int                LDA, LDB, LDE, N;
      // ..
      // .. Array Arguments ..
      double               A( LDA, * ), ALPHAI( * ), ALPHAR( * ), B( LDB, * ), BETA( * ), E( LDE, * ), RESULT( 2 ), WORK( * );
      // ..

// =====================================================================

      // .. Parameters ..
      double               ZERO, ONE, TEN;
      const              ZERO = 0.0, ONE = 1.0, TEN = 10.0 ;
      // ..
      // .. Local Scalars ..
      bool               ILCPLX;
      String             NORMAB, TRANS;
      int                J, JVEC;
      double               ABMAX, ACOEF, ALFMAX, ANORM, BCOEFI, BCOEFR, BETMAX, BNORM, ENORM, ENRMER, ERRNRM, SAFMAX, SAFMIN, SALFI, SALFR, SBETA, SCALE, TEMP1, ULP;
      // ..
      // .. External Functions ..
      //- REAL               SLAMCH, SLANGE;
      // EXTERNAL SLAMCH, SLANGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL SGEMV
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, REAL
      // ..
      // .. Executable Statements ..

      RESULT[1] = ZERO;
      RESULT[2] = ZERO;
      if (N <= 0) return;

      SAFMIN = SLAMCH( 'Safe minimum' );
      SAFMAX = ONE / SAFMIN;
      ULP = SLAMCH( 'Epsilon' )*SLAMCH( 'Base' );

      if ( LEFT ) {
         TRANS = 'T';
         NORMAB = 'I';
      } else {
         TRANS = 'N';
         NORMAB = 'O';
      }

      // Norm of A, B, and E:

      ANORM = max( SLANGE( NORMAB, N, N, A, LDA, WORK ), SAFMIN );
      BNORM = max( SLANGE( NORMAB, N, N, B, LDB, WORK ), SAFMIN );
      ENORM = max( SLANGE( 'O', N, N, E, LDE, WORK ), ULP );
      ALFMAX = SAFMAX / max( ONE, BNORM );
      BETMAX = SAFMAX / max( ONE, ANORM );

      // Compute error matrix.
      // Column i = ( b(i) A - a(i) B ) E(i) / max( |a(i) B|, |b(i) A| )

      ILCPLX = false;
      for (JVEC = 1; JVEC <= N; JVEC++) { // 10
         if ( ILCPLX ) {

            // 2nd Eigenvalue/-vector of pair -- do nothing

            ILCPLX = false;
         } else {
            SALFR = ALPHAR( JVEC );
            SALFI = ALPHAI( JVEC );
            SBETA = BETA( JVEC );
            if ( SALFI == ZERO ) {

               // Real eigenvalue and -vector

               ABMAX = max( ( SALFR ).abs(), ( SBETA ).abs() );
               if ( ( SALFR ).abs() > ALFMAX || ( SBETA ).abs() > BETMAX || ABMAX < ONE ) {
                  SCALE = ONE / max( ABMAX, SAFMIN );
                  SALFR = SCALE*SALFR;
                  SBETA = SCALE*SBETA;
               }
               SCALE = ONE / max( ( SALFR ).abs()*BNORM, ( SBETA ).abs()*ANORM, SAFMIN );
               ACOEF = SCALE*SBETA;
               BCOEFR = SCALE*SALFR;
               sgemv(TRANS, N, N, ACOEF, A, LDA, E( 1, JVEC ), 1, ZERO, WORK( N*( JVEC-1 )+1 ), 1 );
               sgemv(TRANS, N, N, -BCOEFR, B, LDA, E( 1, JVEC ), 1, ONE, WORK( N*( JVEC-1 )+1 ), 1 );
            } else {

               // Complex conjugate pair

               ILCPLX = true;
               if ( JVEC == N ) {
                  RESULT[1] = TEN / ULP;
                  return;
               }
               ABMAX = max( ( SALFR ).abs()+( SALFI ).abs(), ( SBETA ).abs() );
               if ( ( SALFR ).abs()+( SALFI ).abs() > ALFMAX || ( SBETA ).abs() > BETMAX || ABMAX < ONE ) {
                  SCALE = ONE / max( ABMAX, SAFMIN );
                  SALFR = SCALE*SALFR;
                  SALFI = SCALE*SALFI;
                  SBETA = SCALE*SBETA;
               }
               SCALE = ONE / max( ( ( SALFR ).abs()+( SALFI ).abs() )*BNORM, ( SBETA ).abs()*ANORM, SAFMIN );
               ACOEF = SCALE*SBETA;
               BCOEFR = SCALE*SALFR;
               BCOEFI = SCALE*SALFI;
               if ( LEFT ) {
                  BCOEFI = -BCOEFI;
               }

               sgemv(TRANS, N, N, ACOEF, A, LDA, E( 1, JVEC ), 1, ZERO, WORK( N*( JVEC-1 )+1 ), 1 );
               sgemv(TRANS, N, N, -BCOEFR, B, LDA, E( 1, JVEC ), 1, ONE, WORK( N*( JVEC-1 )+1 ), 1 );
               sgemv(TRANS, N, N, BCOEFI, B, LDA, E( 1, JVEC+1 ), 1, ONE, WORK( N*( JVEC-1 )+1 ), 1 );

               sgemv(TRANS, N, N, ACOEF, A, LDA, E( 1, JVEC+1 ), 1, ZERO, WORK( N*JVEC+1 ), 1 );
               sgemv(TRANS, N, N, -BCOEFI, B, LDA, E( 1, JVEC ), 1, ONE, WORK( N*JVEC+1 ), 1 );
               sgemv(TRANS, N, N, -BCOEFR, B, LDA, E( 1, JVEC+1 ), 1, ONE, WORK( N*JVEC+1 ), 1 );
            }
         }
      } // 10

      ERRNRM = SLANGE( 'One', N, N, WORK, N, WORK( N**2+1 ) ) / ENORM;

      // Compute RESULT(1)

      RESULT[1] = ERRNRM / ULP;

      // Normalization of E:

      ENRMER = ZERO;
      ILCPLX = false;
      for (JVEC = 1; JVEC <= N; JVEC++) { // 40
         if ( ILCPLX ) {
            ILCPLX = false;
         } else {
            TEMP1 = ZERO;
            if ( ALPHAI( JVEC ) == ZERO ) {
               for (J = 1; J <= N; J++) { // 20
                  TEMP1 = max( TEMP1, ( E( J, JVEC ) ).abs() );
               } // 20
               ENRMER = max( ENRMER, ( TEMP1-ONE ).abs() );
            } else {
               ILCPLX = true;
               for (J = 1; J <= N; J++) { // 30
                  TEMP1 = max( TEMP1, ( E( J, JVEC ) ).abs()+ ( E( J, JVEC+1 ) ).abs() );
               } // 30
               ENRMER = max( ENRMER, ( TEMP1-ONE ).abs() );
            }
         }
      } // 40

      // Compute RESULT(2) : the normalization error in E.

      RESULT[2] = ENRMER / ( REAL( N )*ULP );

      return;
      }
