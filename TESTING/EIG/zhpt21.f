      SUBROUTINE ZHPT21( ITYPE, UPLO, N, KBAND, AP, D, E, U, LDU, VP, TAU, WORK, RWORK, RESULT )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                ITYPE, KBAND, LDU, N;
      // ..
      // .. Array Arguments ..
      double             D( * ), E( * ), RESULT( 2 ), RWORK( * );
      COMPLEX*16         AP( * ), TAU( * ), U( LDU, * ), VP( * ), WORK( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, ONE, TEN;
      const              ZERO = 0.0D+0, ONE = 1.0D+0, TEN = 10.0D+0 ;
      double             HALF;
      const              HALF = 1.0D+0 / 2.0D+0 ;
      COMPLEX*16         CZERO, CONE
      const              CZERO = ( 0.0D+0, 0.0D+0 ), CONE = ( 1.0D+0, 0.0D+0 ) ;
      // ..
      // .. Local Scalars ..
      bool               LOWER;
      String             CUPLO;
      int                IINFO, J, JP, JP1, JR, LAP;
      double             ANORM, ULP, UNFL, WNORM;
      COMPLEX*16         TEMP, VSAVE
      // ..
      // .. External Functions ..
      bool               LSAME;
      double             DLAMCH, ZLANGE, ZLANHP;
      COMPLEX*16         ZDOTC
      // EXTERNAL LSAME, DLAMCH, ZLANGE, ZLANHP, ZDOTC
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZAXPY, ZCOPY, ZGEMM, ZHPMV, ZHPR, ZHPR2, ZLACPY, ZLASET, ZUPMTR
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE, DCMPLX, MAX, MIN
      // ..
      // .. Executable Statements ..

      // Constants

      RESULT( 1 ) = ZERO
      IF( ITYPE.EQ.1 ) RESULT( 2 ) = ZERO       IF( N.LE.0 ) RETURN

      LAP = ( N*( N+1 ) ) / 2

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
         ANORM = MAX( ZLANHP( '1', CUPLO, N, AP, RWORK ), UNFL )
      }

      // Compute error matrix:

      if ( ITYPE.EQ.1 ) {

         // ITYPE=1: error = A - U S U**H

         CALL ZLASET( 'Full', N, N, CZERO, CZERO, WORK, N )
         CALL ZCOPY( LAP, AP, 1, WORK, 1 )

         DO 10 J = 1, N
            CALL ZHPR( CUPLO, N, -D( J ), U( 1, J ), 1, WORK )
   10    CONTINUE

         if ( N.GT.1 .AND. KBAND.EQ.1 ) {
            DO 20 J = 2, N - 1
               CALL ZHPR2( CUPLO, N, -DCMPLX( E( J ) ), U( 1, J ), 1, U( 1, J-1 ), 1, WORK )
   20       CONTINUE
         }
         WNORM = ZLANHP( '1', CUPLO, N, WORK, RWORK )

      } else if ( ITYPE.EQ.2 ) {

         // ITYPE=2: error = V S V**H - A

         CALL ZLASET( 'Full', N, N, CZERO, CZERO, WORK, N )

         if ( LOWER ) {
            WORK( LAP ) = D( N )
            DO 40 J = N - 1, 1, -1
               JP = ( ( 2*N-J )*( J-1 ) ) / 2
               JP1 = JP + N - J
               if ( KBAND.EQ.1 ) {
                  WORK( JP+J+1 ) = ( CONE-TAU( J ) )*E( J )
                  DO 30 JR = J + 2, N
                     WORK( JP+JR ) = -TAU( J )*E( J )*VP( JP+JR )
   30             CONTINUE
               }

               if ( TAU( J ).NE.CZERO ) {
                  VSAVE = VP( JP+J+1 )
                  VP( JP+J+1 ) = CONE
                  CALL ZHPMV( 'L', N-J, CONE, WORK( JP1+J+1 ), VP( JP+J+1 ), 1, CZERO, WORK( LAP+1 ), 1 )                   TEMP = -HALF*TAU( J )*ZDOTC( N-J, WORK( LAP+1 ), 1, VP( JP+J+1 ), 1 )                   CALL ZAXPY( N-J, TEMP, VP( JP+J+1 ), 1, WORK( LAP+1 ), 1 )                   CALL ZHPR2( 'L', N-J, -TAU( J ), VP( JP+J+1 ), 1, WORK( LAP+1 ), 1, WORK( JP1+J+1 ) )

                  VP( JP+J+1 ) = VSAVE
               }
               WORK( JP+J ) = D( J )
   40       CONTINUE
         } else {
            WORK( 1 ) = D( 1 )
            DO 60 J = 1, N - 1
               JP = ( J*( J-1 ) ) / 2
               JP1 = JP + J
               if ( KBAND.EQ.1 ) {
                  WORK( JP1+J ) = ( CONE-TAU( J ) )*E( J )
                  DO 50 JR = 1, J - 1
                     WORK( JP1+JR ) = -TAU( J )*E( J )*VP( JP1+JR )
   50             CONTINUE
               }

               if ( TAU( J ).NE.CZERO ) {
                  VSAVE = VP( JP1+J )
                  VP( JP1+J ) = CONE
                  CALL ZHPMV( 'U', J, CONE, WORK, VP( JP1+1 ), 1, CZERO, WORK( LAP+1 ), 1 )                   TEMP = -HALF*TAU( J )*ZDOTC( J, WORK( LAP+1 ), 1, VP( JP1+1 ), 1 )                   CALL ZAXPY( J, TEMP, VP( JP1+1 ), 1, WORK( LAP+1 ), 1 )                   CALL ZHPR2( 'U', J, -TAU( J ), VP( JP1+1 ), 1, WORK( LAP+1 ), 1, WORK )
                  VP( JP1+J ) = VSAVE
               }
               WORK( JP1+J+1 ) = D( J+1 )
   60       CONTINUE
         }

         DO 70 J = 1, LAP
            WORK( J ) = WORK( J ) - AP( J )
   70    CONTINUE
         WNORM = ZLANHP( '1', CUPLO, N, WORK, RWORK )

      } else if ( ITYPE.EQ.3 ) {

         // ITYPE=3: error = U V**H - I

         IF( N.LT.2 ) RETURN
         CALL ZLACPY( ' ', N, N, U, LDU, WORK, N )
         CALL ZUPMTR( 'R', CUPLO, 'C', N, N, VP, TAU, WORK, N, WORK( N**2+1 ), IINFO )
         if ( IINFO.NE.0 ) {
            RESULT( 1 ) = TEN / ULP
            RETURN
         }

         DO 80 J = 1, N
            WORK( ( N+1 )*( J-1 )+1 ) = WORK( ( N+1 )*( J-1 )+1 ) - CONE
   80    CONTINUE

         WNORM = ZLANGE( '1', N, N, WORK, N, RWORK )
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

      // Compute  U U**H - I

      if ( ITYPE.EQ.1 ) {
         CALL ZGEMM( 'N', 'C', N, N, N, CONE, U, LDU, U, LDU, CZERO, WORK, N )

         DO 90 J = 1, N
            WORK( ( N+1 )*( J-1 )+1 ) = WORK( ( N+1 )*( J-1 )+1 ) - CONE
   90    CONTINUE

         RESULT( 2 ) = MIN( ZLANGE( '1', N, N, WORK, N, RWORK ), DBLE( N ) ) / ( N*ULP )
      }

      RETURN

      // End of ZHPT21

      }
