      SUBROUTINE SQRT03( M, N, K, AF, C, CC, Q, LDA, TAU, WORK, LWORK, RWORK, RESULT )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      int                K, LDA, LWORK, M, N;
*     ..
*     .. Array Arguments ..
      REAL               AF( LDA, * ), C( LDA, * ), CC( LDA, * ), Q( LDA, * ), RESULT( * ), RWORK( * ), TAU( * ), WORK( LWORK )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ONE
      PARAMETER          ( ONE = 1.0E0 )
      REAL               ROGUE
      PARAMETER          ( ROGUE = -1.0E+10 )
*     ..
*     .. Local Scalars ..
      String             SIDE, TRANS;
      int                INFO, ISIDE, ITRANS, J, MC, NC;
      REAL               CNORM, EPS, RESID
*     ..
*     .. External Functions ..
      bool               LSAME;
      REAL               SLAMCH, SLANGE
      EXTERNAL           LSAME, SLAMCH, SLANGE
*     ..
*     .. External Subroutines ..
      EXTERNAL           SGEMM, SLACPY, SLARNV, SLASET, SORGQR, SORMQR
*     ..
*     .. Local Arrays ..
      int                ISEED( 4 );
*     ..
*     .. Intrinsic Functions ..
      // INTRINSIC MAX, REAL
*     ..
*     .. Scalars in Common ..
      String             SRNAMT;
*     ..
*     .. Common blocks ..
      COMMON             / SRNAMC / SRNAMT
*     ..
*     .. Data statements ..
      DATA               ISEED / 1988, 1989, 1990, 1991 /
*     ..
*     .. Executable Statements ..
*
      EPS = SLAMCH( 'Epsilon' )
*
*     Copy the first k columns of the factorization to the array Q
*
      CALL SLASET( 'Full', M, M, ROGUE, ROGUE, Q, LDA )
      CALL SLACPY( 'Lower', M-1, K, AF( 2, 1 ), LDA, Q( 2, 1 ), LDA )
*
*     Generate the m-by-m matrix Q
*
      SRNAMT = 'SORGQR'
      CALL SORGQR( M, M, K, Q, LDA, TAU, WORK, LWORK, INFO )
*
      DO 30 ISIDE = 1, 2
         IF( ISIDE.EQ.1 ) THEN
            SIDE = 'L'
            MC = M
            NC = N
         ELSE
            SIDE = 'R'
            MC = N
            NC = M
         END IF
*
*        Generate MC by NC matrix C
*
         DO 10 J = 1, NC
            CALL SLARNV( 2, ISEED, MC, C( 1, J ) )
   10    CONTINUE
         CNORM = SLANGE( '1', MC, NC, C, LDA, RWORK )
         IF( CNORM.EQ.0.0 ) CNORM = ONE
*
         DO 20 ITRANS = 1, 2
            IF( ITRANS.EQ.1 ) THEN
               TRANS = 'N'
            ELSE
               TRANS = 'T'
            END IF
*
*           Copy C
*
            CALL SLACPY( 'Full', MC, NC, C, LDA, CC, LDA )
*
*           Apply Q or Q' to C
*
            SRNAMT = 'SORMQR'
            CALL SORMQR( SIDE, TRANS, MC, NC, K, AF, LDA, TAU, CC, LDA, WORK, LWORK, INFO )
*
*           Form explicit product and subtract
*
            IF( LSAME( SIDE, 'L' ) ) THEN
               CALL SGEMM( TRANS, 'No transpose', MC, NC, MC, -ONE, Q, LDA, C, LDA, ONE, CC, LDA )
            ELSE
               CALL SGEMM( 'No transpose', TRANS, MC, NC, NC, -ONE, C, LDA, Q, LDA, ONE, CC, LDA )
            END IF
*
*           Compute error in the difference
*
            RESID = SLANGE( '1', MC, NC, CC, LDA, RWORK )
            RESULT( ( ISIDE-1 )*2+ITRANS ) = RESID / ( REAL( MAX( 1, M ) )*CNORM*EPS )
*
   20    CONTINUE
   30 CONTINUE
*
      RETURN
*
*     End of SQRT03
*
      END
