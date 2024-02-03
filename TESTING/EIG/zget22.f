      SUBROUTINE ZGET22( TRANSA, TRANSE, TRANSW, N, A, LDA, E, LDE, W, WORK, RWORK, RESULT )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             TRANSA, TRANSE, TRANSW;
      int                LDA, LDE, N;
      // ..
      // .. Array Arguments ..
      double             RESULT( 2 ), RWORK( * );
      COMPLEX*16         A( LDA, * ), E( LDE, * ), W( * ), WORK( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0D+0, ONE = 1.0D+0 ;
      COMPLEX*16         CZERO, CONE
      const              CZERO = ( 0.0D+0, 0.0D+0 ), CONE = ( 1.0D+0, 0.0D+0 ) ;
      // ..
      // .. Local Scalars ..
      String             NORMA, NORME;
      int                ITRNSE, ITRNSW, J, JCOL, JOFF, JROW, JVEC;
      double             ANORM, ENORM, ENRMAX, ENRMIN, ERRNRM, TEMP1, ULP, UNFL;
      COMPLEX*16         WTEMP
      // ..
      // .. External Functions ..
      bool               LSAME;
      double             DLAMCH, ZLANGE;
      // EXTERNAL LSAME, DLAMCH, ZLANGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZGEMM, ZLASET
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, DCONJG, DIMAG, MAX, MIN
      // ..
      // .. Executable Statements ..

      // Initialize RESULT (in case N=0)

      RESULT( 1 ) = ZERO
      RESULT( 2 ) = ZERO
      if (N.LE.0) RETURN;

      UNFL = DLAMCH( 'Safe minimum' )
      ULP = DLAMCH( 'Precision' )

      ITRNSE = 0
      ITRNSW = 0
      NORMA = 'O'
      NORME = 'O'

      if ( LSAME( TRANSA, 'T' ) || LSAME( TRANSA, 'C' ) ) {
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
      if ( ITRNSE == 0 ) {
         for (JVEC = 1; JVEC <= N; JVEC++) { // 20
            TEMP1 = ZERO
            for (J = 1; J <= N; J++) { // 10
               TEMP1 = MAX( TEMP1, ABS( DBLE( E( J, JVEC ) ) )+ ABS( DIMAG( E( J, JVEC ) ) ) )
            } // 10
            ENRMIN = MIN( ENRMIN, TEMP1 )
            ENRMAX = MAX( ENRMAX, TEMP1 )
         } // 20
      } else {
         for (JVEC = 1; JVEC <= N; JVEC++) { // 30
            RWORK( JVEC ) = ZERO
         } // 30

         for (J = 1; J <= N; J++) { // 50
            for (JVEC = 1; JVEC <= N; JVEC++) { // 40
               RWORK( JVEC ) = MAX( RWORK( JVEC ), ABS( DBLE( E( JVEC, J ) ) )+ ABS( DIMAG( E( JVEC, J ) ) ) )
            } // 40
         } // 50

         for (JVEC = 1; JVEC <= N; JVEC++) { // 60
            ENRMIN = MIN( ENRMIN, RWORK( JVEC ) )
            ENRMAX = MAX( ENRMAX, RWORK( JVEC ) )
         } // 60
      }

      // Norm of A:

      ANORM = MAX( ZLANGE( NORMA, N, N, A, LDA, RWORK ), UNFL )

      // Norm of E:

      ENORM = MAX( ZLANGE( NORME, N, N, E, LDE, RWORK ), ULP )

      // Norm of error:

      // Error =  AE - EW

      zlaset('Full', N, N, CZERO, CZERO, WORK, N );

      JOFF = 0
      for (JCOL = 1; JCOL <= N; JCOL++) { // 100
         if ( ITRNSW == 0 ) {
            WTEMP = W( JCOL )
         } else {
            WTEMP = DCONJG( W( JCOL ) )
         }

         if ( ITRNSE == 0 ) {
            for (JROW = 1; JROW <= N; JROW++) { // 70
               WORK( JOFF+JROW ) = E( JROW, JCOL )*WTEMP
            } // 70
         } else if ( ITRNSE == 1 ) {
            for (JROW = 1; JROW <= N; JROW++) { // 80
               WORK( JOFF+JROW ) = E( JCOL, JROW )*WTEMP
            } // 80
         } else {
            for (JROW = 1; JROW <= N; JROW++) { // 90
               WORK( JOFF+JROW ) = DCONJG( E( JCOL, JROW ) )*WTEMP
            } // 90
         }
         JOFF = JOFF + N
      } // 100

      zgemm(TRANSA, TRANSE, N, N, N, CONE, A, LDA, E, LDE, -CONE, WORK, N );

      ERRNRM = ZLANGE( 'One', N, N, WORK, N, RWORK ) / ENORM

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

      RESULT( 2 ) = MAX( ABS( ENRMAX-ONE ), ABS( ENRMIN-ONE ) ) / ( DBLE( N )*ULP )

      RETURN

      // End of ZGET22

      }
