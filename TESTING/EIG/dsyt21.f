      SUBROUTINE DSYT21( ITYPE, UPLO, N, KBAND, A, LDA, D, E, U, LDU, V, LDV, TAU, WORK, RESULT )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                ITYPE, KBAND, LDA, LDU, LDV, N;
      // ..
      // .. Array Arguments ..
      double             A( LDA, * ), D( * ), E( * ), RESULT( 2 ), TAU( * ), U( LDU, * ), V( LDV, * ), WORK( * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, ONE, TEN;
      const              ZERO = 0.0D0, ONE = 1.0D0, TEN = 10.0D0 ;
      // ..
      // .. Local Scalars ..
      bool               LOWER;
      String             CUPLO;
      int                IINFO, J, JCOL, JR, JROW;
      double             ANORM, ULP, UNFL, VSAVE, WNORM;
      // ..
      // .. External Functions ..
      bool               LSAME;
      double             DLAMCH, DLANGE, DLANSY;
      // EXTERNAL LSAME, DLAMCH, DLANGE, DLANSY
      // ..
      // .. External Subroutines ..
      // EXTERNAL DGEMM, DLACPY, DLARFY, DLASET, DORM2L, DORM2R, DSYR, DSYR2
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE, MAX, MIN
      // ..
      // .. Executable Statements ..

      RESULT( 1 ) = ZERO
      IF( ITYPE.EQ.1 ) RESULT( 2 ) = ZERO       IF( N.LE.0 ) RETURN

      if ( LSAME( UPLO, 'U' ) ) {
         LOWER = .FALSE.
         CUPLO = 'U'
      } else {
         LOWER = .TRUE.
         CUPLO = 'L'
      }

      UNFL = DLAMCH( 'Safe minimum' )
      ULP = DLAMCH( 'Epsilon' )*DLAMCH( 'Base' )

      // Some Error Checks

      if ( ITYPE.LT.1 .OR. ITYPE.GT.3 ) {
         RESULT( 1 ) = TEN / ULP
         RETURN
      }

      // Do Test 1

      // Norm of A:

      if ( ITYPE.EQ.3 ) {
         ANORM = ONE
      } else {
         ANORM = MAX( DLANSY( '1', CUPLO, N, A, LDA, WORK ), UNFL )
      }

      // Compute error matrix:

      if ( ITYPE.EQ.1 ) {

         // ITYPE=1: error = A - U S U**T

         dlaset('Full', N, N, ZERO, ZERO, WORK, N );
         dlacpy(CUPLO, N, N, A, LDA, WORK, N );

         for (J = 1; J <= N; J++) { // 10
            dsyr(CUPLO, N, -D( J ), U( 1, J ), 1, WORK, N );
   10    CONTINUE

         if ( N.GT.1 .AND. KBAND.EQ.1 ) {
            DO 20 J = 1, N - 1
               dsyr2(CUPLO, N, -E( J ), U( 1, J ), 1, U( 1, J+1 ), 1, WORK, N );
   20       CONTINUE
         }
         WNORM = DLANSY( '1', CUPLO, N, WORK, N, WORK( N**2+1 ) )

      } else if ( ITYPE.EQ.2 ) {

         // ITYPE=2: error = V S V**T - A

         dlaset('Full', N, N, ZERO, ZERO, WORK, N );

         if ( LOWER ) {
            WORK( N**2 ) = D( N )
            DO 40 J = N - 1, 1, -1
               if ( KBAND.EQ.1 ) {
                  WORK( ( N+1 )*( J-1 )+2 ) = ( ONE-TAU( J ) )*E( J )
                  DO 30 JR = J + 2, N
                     WORK( ( J-1 )*N+JR ) = -TAU( J )*E( J )*V( JR, J )
   30             CONTINUE
               }

               VSAVE = V( J+1, J )
               V( J+1, J ) = ONE
               dlarfy('L', N-J, V( J+1, J ), 1, TAU( J ), WORK( ( N+1 )*J+1 ), N, WORK( N**2+1 ) );
               V( J+1, J ) = VSAVE
               WORK( ( N+1 )*( J-1 )+1 ) = D( J )
   40       CONTINUE
         } else {
            WORK( 1 ) = D( 1 )
            DO 60 J = 1, N - 1
               if ( KBAND.EQ.1 ) {
                  WORK( ( N+1 )*J ) = ( ONE-TAU( J ) )*E( J )
                  DO 50 JR = 1, J - 1
                     WORK( J*N+JR ) = -TAU( J )*E( J )*V( JR, J+1 )
   50             CONTINUE
               }

               VSAVE = V( J, J+1 )
               V( J, J+1 ) = ONE
               dlarfy('U', J, V( 1, J+1 ), 1, TAU( J ), WORK, N, WORK( N**2+1 ) );
               V( J, J+1 ) = VSAVE
               WORK( ( N+1 )*J+1 ) = D( J+1 )
   60       CONTINUE
         }

         for (JCOL = 1; JCOL <= N; JCOL++) { // 90
            if ( LOWER ) {
               for (JROW = JCOL; JROW <= N; JROW++) { // 70
                  WORK( JROW+N*( JCOL-1 ) ) = WORK( JROW+N*( JCOL-1 ) ) - A( JROW, JCOL )
   70          CONTINUE
            } else {
               for (JROW = 1; JROW <= JCOL; JROW++) { // 80
                  WORK( JROW+N*( JCOL-1 ) ) = WORK( JROW+N*( JCOL-1 ) ) - A( JROW, JCOL )
   80          CONTINUE
            }
   90    CONTINUE
         WNORM = DLANSY( '1', CUPLO, N, WORK, N, WORK( N**2+1 ) )

      } else if ( ITYPE.EQ.3 ) {

         // ITYPE=3: error = U V**T - I

         IF( N.LT.2 ) RETURN
         dlacpy(' ', N, N, U, LDU, WORK, N );
         if ( LOWER ) {
            dorm2r('R', 'T', N, N-1, N-1, V( 2, 1 ), LDV, TAU, WORK( N+1 ), N, WORK( N**2+1 ), IINFO );
         } else {
            dorm2l('R', 'T', N, N-1, N-1, V( 1, 2 ), LDV, TAU, WORK, N, WORK( N**2+1 ), IINFO );
         }
         if ( IINFO.NE.0 ) {
            RESULT( 1 ) = TEN / ULP
            RETURN
         }

         for (J = 1; J <= N; J++) { // 100
            WORK( ( N+1 )*( J-1 )+1 ) = WORK( ( N+1 )*( J-1 )+1 ) - ONE
  100    CONTINUE

         WNORM = DLANGE( '1', N, N, WORK, N, WORK( N**2+1 ) )
      }

      if ( ANORM.GT.WNORM ) {
         RESULT( 1 ) = ( WNORM / ANORM ) / ( N*ULP )
      } else {
         if ( ANORM.LT.ONE ) {
            RESULT( 1 ) = ( MIN( WNORM, N*ANORM ) / ANORM ) / ( N*ULP )
         } else {
            RESULT( 1 ) = MIN( WNORM / ANORM, DBLE( N ) ) / ( N*ULP )
         }
      }

      // Do Test 2

      // Compute  U U**T - I

      if ( ITYPE.EQ.1 ) {
         dgemm('N', 'C', N, N, N, ONE, U, LDU, U, LDU, ZERO, WORK, N );

         for (J = 1; J <= N; J++) { // 110
            WORK( ( N+1 )*( J-1 )+1 ) = WORK( ( N+1 )*( J-1 )+1 ) - ONE
  110    CONTINUE

         RESULT( 2 ) = MIN( DLANGE( '1', N, N, WORK, N, WORK( N**2+1 ) ), DBLE( N ) ) / ( N*ULP )
      }

      RETURN

      // End of DSYT21

      }
