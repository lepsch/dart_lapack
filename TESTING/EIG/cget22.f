      SUBROUTINE CGET22( TRANSA, TRANSE, TRANSW, N, A, LDA, E, LDE, W, WORK, RWORK, RESULT )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             TRANSA, TRANSE, TRANSW;
      int                LDA, LDE, N;
      // ..
      // .. Array Arguments ..
      REAL               RESULT( 2 ), RWORK( * )
      COMPLEX            A( LDA, * ), E( LDE, * ), W( * ), WORK( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE
      const              ZERO = 0.0E+0, ONE = 1.0E+0 ;
      COMPLEX            CZERO, CONE
      const              CZERO = ( 0.0E+0, 0.0E+0 ), CONE = ( 1.0E+0, 0.0E+0 ) ;
      // ..
      // .. Local Scalars ..
      String             NORMA, NORME;
      int                ITRNSE, ITRNSW, J, JCOL, JOFF, JROW, JVEC;
      REAL               ANORM, ENORM, ENRMAX, ENRMIN, ERRNRM, TEMP1, ULP, UNFL
      COMPLEX            WTEMP
      // ..
      // .. External Functions ..
      bool               LSAME;
      REAL               CLANGE, SLAMCH
      // EXTERNAL LSAME, CLANGE, SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL CGEMM, CLASET
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, AIMAG, CONJG, MAX, MIN, REAL
      // ..
      // .. Executable Statements ..

      // Initialize RESULT (in case N=0)

      RESULT( 1 ) = ZERO
      RESULT( 2 ) = ZERO
      IF( N.LE.0 ) RETURN

      UNFL = SLAMCH( 'Safe minimum' )
      ULP = SLAMCH( 'Precision' )

      ITRNSE = 0
      ITRNSW = 0
      NORMA = 'O'
      NORME = 'O'

      if ( LSAME( TRANSA, 'T' ) .OR. LSAME( TRANSA, 'C' ) ) {
         NORMA = 'I'
      }

      if ( LSAME( TRANSE, 'T' ) ) {
         ITRNSE = 1
         NORME = 'I'
      } else if ( LSAME( TRANSE, 'C' ) ) {
         ITRNSE = 2
         NORME = 'I'
      }

      if ( LSAME( TRANSW, 'C' ) ) {
         ITRNSW = 1
      }

      // Normalization of E:

      ENRMIN = ONE / ULP
      ENRMAX = ZERO
      if ( ITRNSE.EQ.0 ) {
         for (JVEC = 1; JVEC <= N; JVEC++) { // 20
            TEMP1 = ZERO
            for (J = 1; J <= N; J++) { // 10
               TEMP1 = MAX( TEMP1, ABS( REAL( E( J, JVEC ) ) )+ ABS( AIMAG( E( J, JVEC ) ) ) )
   10       CONTINUE
            ENRMIN = MIN( ENRMIN, TEMP1 )
            ENRMAX = MAX( ENRMAX, TEMP1 )
   20    CONTINUE
      } else {
         for (JVEC = 1; JVEC <= N; JVEC++) { // 30
            RWORK( JVEC ) = ZERO
   30    CONTINUE

         for (J = 1; J <= N; J++) { // 50
            for (JVEC = 1; JVEC <= N; JVEC++) { // 40
               RWORK( JVEC ) = MAX( RWORK( JVEC ), ABS( REAL( E( JVEC, J ) ) )+ ABS( AIMAG( E( JVEC, J ) ) ) )
   40       CONTINUE
   50    CONTINUE

         for (JVEC = 1; JVEC <= N; JVEC++) { // 60
            ENRMIN = MIN( ENRMIN, RWORK( JVEC ) )
            ENRMAX = MAX( ENRMAX, RWORK( JVEC ) )
   60    CONTINUE
      }

      // Norm of A:

      ANORM = MAX( CLANGE( NORMA, N, N, A, LDA, RWORK ), UNFL )

      // Norm of E:

      ENORM = MAX( CLANGE( NORME, N, N, E, LDE, RWORK ), ULP )

      // Norm of error:

      // Error =  AE - EW

      claset('Full', N, N, CZERO, CZERO, WORK, N );

      JOFF = 0
      for (JCOL = 1; JCOL <= N; JCOL++) { // 100
         if ( ITRNSW.EQ.0 ) {
            WTEMP = W( JCOL )
         } else {
            WTEMP = CONJG( W( JCOL ) )
         }

         if ( ITRNSE.EQ.0 ) {
            for (JROW = 1; JROW <= N; JROW++) { // 70
               WORK( JOFF+JROW ) = E( JROW, JCOL )*WTEMP
   70       CONTINUE
         } else if ( ITRNSE.EQ.1 ) {
            for (JROW = 1; JROW <= N; JROW++) { // 80
               WORK( JOFF+JROW ) = E( JCOL, JROW )*WTEMP
   80       CONTINUE
         } else {
            for (JROW = 1; JROW <= N; JROW++) { // 90
               WORK( JOFF+JROW ) = CONJG( E( JCOL, JROW ) )*WTEMP
   90       CONTINUE
         }
         JOFF = JOFF + N
  100 CONTINUE

      cgemm(TRANSA, TRANSE, N, N, N, CONE, A, LDA, E, LDE, -CONE, WORK, N );

      ERRNRM = CLANGE( 'One', N, N, WORK, N, RWORK ) / ENORM

      // Compute RESULT(1) (avoiding under/overflow)

      if ( ANORM.GT.ERRNRM ) {
         RESULT( 1 ) = ( ERRNRM / ANORM ) / ULP
      } else {
         if ( ANORM.LT.ONE ) {
            RESULT( 1 ) = ONE / ULP
         } else {
            RESULT( 1 ) = MIN( ERRNRM / ANORM, ONE ) / ULP
         }
      }

      // Compute RESULT(2) : the normalization error in E.

      RESULT( 2 ) = MAX( ABS( ENRMAX-ONE ), ABS( ENRMIN-ONE ) ) / ( REAL( N )*ULP )

      RETURN

      // End of CGET22

      }
